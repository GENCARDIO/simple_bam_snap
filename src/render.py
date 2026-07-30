"""Matplotlib renderer.

Replaces the old approach of shelling out to `samtools tview`, capturing its
text table, and re-parsing that text with string splits. Here every read is
drawn from its own parsed CIGAR blocks, so insertions/deletions/soft-clips/
mismatches are geometrically exact instead of guessed from column spacing.
"""
from __future__ import annotations

from typing import List, Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
from matplotlib.ticker import MaxNLocator, FuncFormatter

from src.read_model import AlignedRead
from src.reference import ReferenceWindow

# Colors: base identity follows the standard genome-browser convention
# (A green / C blue / G orange / T red); hues are drawn from the
# colorblind-validated categorical set rather than picked freehand.
BASE_COLORS = {
    "A": "#008300",
    "C": "#2a78d6",
    "G": "#eb6834",
    "T": "#e34948",
    "N": "#898781",
}
NORMAL_FILL = "#b0b0b0"  # concordant / non-discordant reads: plain grey, not strand-colored
INSERTION_COLOR = "#4a3aa7"
DELETION_COLOR = "#0b0b0b"
SKIP_COLOR = "#898781"
SOFTCLIP_COLOR = "#1baf7a"  # fallback fill when per-base clip coloring isn't available
COVERAGE_COLOR = "#a3a3a3"
GRIDLINE = "#e1e0d9"
AXIS_INK = "#898781"
PRIMARY_INK = "#0b0b0b"
SECONDARY_INK = "#52514e"

# Discordant-pair fill colors, IGV-equivalent: same categories/roles IGV's
# "color by insert size and pair orientation" mode uses (red = long insert,
# blue = short insert, a blue family for same-strand pairs, green for everted,
# per-chromosome hue for inter-chromosomal mates), drawn from the same
# colorblind-validated categorical set as everything else here rather than
# IGV's undocumented internal hex values.
PAIR_LARGE_INSERT_COLOR = "#e34948"  # red    - insert size bigger than expected
PAIR_SMALL_INSERT_COLOR = "#2a78d6"  # blue   - insert size smaller than expected (darkest of the blue family)
PAIR_FF_COLOR = "#9ec5f4"            # light blue  - FF orientation
PAIR_RR_COLOR = "#5598e7"            # medium blue - RR orientation (darker than FF, lighter than small_insert)
PAIR_EVERTED_COLOR = "#008300"       # green  - RF orientation
CHROM_PALETTE = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]

PAIR_CATEGORY_LABELS = {
    "large_insert": "Large insert (discordant)",
    "small_insert": "Small insert (discordant)",
    "ff": "FF pair (same-strand)",
    "rr": "RR pair (same-strand)",
    "everted": "Everted pair (RF)",
    "interchrom": "Inter-chromosomal pair",
}

MAPQ_ALPHA_FLOOR = 0.15  # fill alpha at mapq==0; scales up to full alpha at mapq_cap

ROW_HEIGHT_IN = 0.22
ROW_MARGIN = 0.12


def chrom_color(chrom: Optional[str]) -> str:
    """Stable per-chromosome hue (same idea as IGV's karyotype coloring) for
    inter-chromosomal mate pairs. Deterministic within a run, not globally
    fixed across all human chromosome names."""
    if not chrom:
        return AXIS_INK
    h = 0
    for ch in chrom:
        h = (h * 31 + ord(ch)) & 0xFFFFFFFF
    return CHROM_PALETTE[h % len(CHROM_PALETTE)]


def compute_coverage(reads: List[AlignedRead], start: int, end: int) -> List[int]:
    """Per-base depth across [start, end), counting only reference-consuming
    match bases (a deleted base is not "covered", matching `samtools depth`)."""
    depth = [0] * max(0, end - start)
    for read in reads:
        for blk in read.blocks:
            if blk.op not in ("M", "=", "X"):
                continue
            lo = max(blk.ref_pos, start)
            hi = min(blk.ref_pos + blk.length, end)
            for pos in range(lo, hi):
                depth[pos - start] += 1
    return depth


def _nice_tick_positions(start: int, end: int, target: int = 8) -> List[int]:
    locator = MaxNLocator(nbins=target, steps=[1, 2, 5, 10])
    ticks = [t for t in locator.tick_values(start, end) if start <= t <= end]
    return [int(t) for t in ticks]


class AlignmentRenderer:
    def __init__(
        self,
        fig_width: float = 14.0,
        dpi: int = 150,
        show_coverage: bool = True,
        annotate_gap: bool = True,
        max_mismatch_render_span: int = 5000,
        pair_colors: bool = True,
        shade_by_mapq: bool = True,
        mapq_cap: int = 60,
    ):
        self.fig_width = fig_width
        self.dpi = dpi
        self.show_coverage = show_coverage
        self.annotate_gap = annotate_gap
        self.max_mismatch_render_span = max_mismatch_render_span
        self.pair_colors = pair_colors
        self.shade_by_mapq = shade_by_mapq
        self.mapq_cap = mapq_cap

    def _legend_row_count(self) -> int:
        n = 3 + (6 if self.pair_colors else 0) + 4
        ncol = min(n, 7)
        return -(-n // ncol)  # ceil

    def _read_style(self, read: AlignedRead):
        """(fill_color, alpha) for a read's main body: hue encodes pair
        discordance category (when enabled), alpha encodes mapping quality -
        low-MAPQ reads get a lighter/more washed-out fill, same idea as IGV's
        "shade by mapping quality"."""
        if self.pair_colors and read.pair_category == "interchrom":
            color = chrom_color(read.mate_chrom)
        elif self.pair_colors and read.pair_category == "large_insert":
            color = PAIR_LARGE_INSERT_COLOR
        elif self.pair_colors and read.pair_category == "small_insert":
            color = PAIR_SMALL_INSERT_COLOR
        elif self.pair_colors and read.pair_category == "ff":
            color = PAIR_FF_COLOR
        elif self.pair_colors and read.pair_category == "rr":
            color = PAIR_RR_COLOR
        elif self.pair_colors and read.pair_category == "everted":
            color = PAIR_EVERTED_COLOR
        else:
            color = NORMAL_FILL

        alpha = 0.5 if (read.is_secondary or read.is_duplicate) else 0.9
        if self.shade_by_mapq and self.mapq_cap > 0:
            mapq_frac = min(max(read.mapq, 0), self.mapq_cap) / self.mapq_cap
            alpha *= MAPQ_ALPHA_FLOOR + (1 - MAPQ_ALPHA_FLOOR) * mapq_frac
        return color, alpha

    def render(
        self,
        rows: List[List[AlignedRead]],
        chrom: str,
        window_start: int,
        window_end: int,
        reference: Optional[ReferenceWindow],
        out_path: str,
        title: str = "",
        layout: str = "pack",
        dropped_reads: int = 0,
        all_reads_for_coverage: Optional[List[AlignedRead]] = None,
    ) -> None:
        span = window_end - window_start
        n_rows = max(len(rows), 1)
        show_ref_track = bool(reference and reference.available)
        show_letters = show_ref_track and (self.fig_width * self.dpi / max(span, 1)) >= 7
        render_base_detail = span <= self.max_mismatch_render_span

        tracks = ["ruler"]
        ratios = [0.35]
        if show_ref_track:
            tracks.append("reference")
            ratios.append(0.32)
        if self.show_coverage:
            tracks.append("coverage")
            ratios.append(1.4)
        tracks.append("alignments")
        ratios.append(max(n_rows * ROW_HEIGHT_IN, ROW_HEIGHT_IN))

        top_margin_in = 0.5   # region title + subtitle
        bottom_margin_in = 0.35 + 0.24 * self._legend_row_count()  # legend
        fig_height = sum(ratios) + top_margin_in + bottom_margin_in
        fig, axes = plt.subplots(
            nrows=len(tracks),
            ncols=1,
            figsize=(self.fig_width, fig_height),
            dpi=self.dpi,
            gridspec_kw={"height_ratios": ratios, "hspace": 0.15},
            sharex=True,
        )
        if len(tracks) == 1:
            axes = [axes]
        ax_by_track = dict(zip(tracks, axes))

        for ax in axes:
            ax.set_xlim(window_start, window_end)
            for spine in ("top", "right", "left"):
                ax.spines[spine].set_visible(False)
            ax.spines["bottom"].set_visible(False)
            ax.tick_params(left=False, labelleft=False)

        tick_positions = _nice_tick_positions(window_start, window_end)

        # --- ruler -----------------------------------------------------
        ruler = ax_by_track["ruler"]
        ruler.set_ylim(0, 1)
        ruler.set_xticks(tick_positions)
        ruler.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{int(x):,}"))
        ruler.tick_params(labelsize=8, labelbottom=True, bottom=True, length=3, colors=AXIS_INK)
        fig.text(
            0.01, 0.995, f"{chrom}:{window_start + 1:,}-{window_end:,} ({span:,} bp)",
            fontsize=10.5, color=PRIMARY_INK, fontweight="bold", va="top", ha="left",
        )
        subtitle = title
        if dropped_reads:
            subtitle = (subtitle + " -- " if subtitle else "") + (
                f"{dropped_reads} lower-priority read(s) not shown (--max-rows)"
            )
        if subtitle:
            fig.text(0.01, 0.995 - 0.42 / fig_height, subtitle, fontsize=8.5,
                     color=SECONDARY_INK, va="top", ha="left")

        # --- reference -----------------------------------------------------
        if show_ref_track:
            ref_ax = ax_by_track["reference"]
            ref_ax.set_ylim(0, 1)
            for pos in range(window_start, window_end):
                base = reference.base_at(pos) or "N"
                color = BASE_COLORS.get(base, BASE_COLORS["N"])
                ref_ax.add_patch(Rectangle((pos, 0), 1, 1, facecolor=color, alpha=0.25, edgecolor="none"))
                if show_letters:
                    ref_ax.text(pos + 0.5, 0.5, base, ha="center", va="center",
                                fontsize=7, color=color, fontweight="bold")

        # --- coverage --------------------------------------------------
        if self.show_coverage:
            cov_ax = ax_by_track["coverage"]
            cov_reads = all_reads_for_coverage if all_reads_for_coverage is not None else [
                r for row in rows for r in row
            ]
            depth = compute_coverage(cov_reads, window_start, window_end)
            max_depth = max(depth) if depth else 0
            xs = [window_start + i + 0.5 for i in range(len(depth))]
            cov_ax.bar(xs, depth, width=1.0, color=COVERAGE_COLOR, alpha=0.85, linewidth=0)
            cov_ax.set_ylim(0, max(max_depth, 1) * 1.15)
            cov_ax.set_yticks([0, max(max_depth, 1)])
            cov_ax.tick_params(left=True, labelleft=True, labelsize=6, colors=AXIS_INK, length=3)
            cov_ax.spines["left"].set_visible(True)
            cov_ax.spines["left"].set_color(AXIS_INK)
            cov_ax.spines["left"].set_linewidth(0.8)
            cov_ax.text(0.0, 1.05, "coverage", transform=cov_ax.transAxes, fontsize=7,
                        color=SECONDARY_INK, va="bottom")

        # --- alignments --------------------------------------------------
        aln_ax = ax_by_track["alignments"]
        aln_ax.set_ylim(n_rows, 0)
        for tick in tick_positions:
            for track in ("reference", "coverage", "alignments"):
                if track in ax_by_track:
                    ax_by_track[track].axvline(tick, color=GRIDLINE, lw=0.6, zorder=0)

        for row_idx, row in enumerate(rows):
            y0 = row_idx + ROW_MARGIN
            h = 1 - 2 * ROW_MARGIN
            for read in row:
                self._draw_read(aln_ax, read, y0, h, render_base_detail)
                if layout == "expand" and self.annotate_gap:
                    label = read.gap_label()
                    if label:
                        aln_ax.text(
                            1.005, y0 + h / 2, label, transform=aln_ax.get_yaxis_transform(),
                            fontsize=6.5, va="center", ha="left", color=SECONDARY_INK, clip_on=False,
                        )

        if not rows:
            aln_ax.text(
                0.5, 0.5, "No alignments in this region", transform=aln_ax.transAxes,
                ha="center", va="center", fontsize=10, color=SECONDARY_INK,
            )

        # --- legend -----------------------------------
        legend_handles = self._legend_handles()
        ncol = min(len(legend_handles), 7)
        fig.legend(
            handles=legend_handles, loc="lower center", ncol=ncol,
            fontsize=7, frameon=False, bbox_to_anchor=(0.5, 0.0),
        )

        fig.subplots_adjust(left=0.05, right=0.92, top=1 - top_margin_in / fig_height,
                            bottom=bottom_margin_in / fig_height)
        fig.savefig(out_path)
        plt.close(fig)

    def _legend_handles(self) -> list:
        handles = [
            Patch(facecolor=NORMAL_FILL, edgecolor="none", label="Normal / concordant"),
            Patch(facecolor=INSERTION_COLOR, edgecolor="none", label="Insertion"),
            Line2D([0], [0], color=DELETION_COLOR, lw=1.5, label="Deletion"),
        ]
        if self.pair_colors:
            handles += [
                Patch(facecolor=PAIR_LARGE_INSERT_COLOR, edgecolor="none", label=PAIR_CATEGORY_LABELS["large_insert"]),
                Patch(facecolor=PAIR_SMALL_INSERT_COLOR, edgecolor="none", label=PAIR_CATEGORY_LABELS["small_insert"]),
                Patch(facecolor=PAIR_FF_COLOR, edgecolor="none", label=PAIR_CATEGORY_LABELS["ff"]),
                Patch(facecolor=PAIR_RR_COLOR, edgecolor="none", label=PAIR_CATEGORY_LABELS["rr"]),
                Patch(facecolor=PAIR_EVERTED_COLOR, edgecolor="none", label=PAIR_CATEGORY_LABELS["everted"]),
                Patch(facecolor=CHROM_PALETTE[0], edgecolor="none", label="Inter-chromosomal (hue = mate chrom)"),
            ]
        # Shown unconditionally: these colors identify mismatch bases AND
        # soft-clipped bases (attached to the read, colored by identity),
        # neither of which requires a reference FASTA to draw.
        handles += [Patch(facecolor=BASE_COLORS[b], edgecolor="none", label=b) for b in "ACGT"]
        return handles

    def render_multi(
        self,
        panels: List[dict],
        chrom: str,
        window_start: int,
        window_end: int,
        reference: Optional[ReferenceWindow],
        out_path: str,
        suptitle: str = "",
    ) -> None:
        """Stack several BAMs' snapshots in one figure, sharing one genomic
        x-axis - the side-by-side comparison view for "does aligner A produce
        more gapped alignments than aligner B here". Each panel is a dict with
        keys: label, rows, all_reads_for_coverage, dropped_reads (optional).
        """
        span = window_end - window_start
        show_ref_track = bool(reference and reference.available)
        show_letters = show_ref_track and (self.fig_width * self.dpi / max(span, 1)) >= 7
        render_base_detail = span <= self.max_mismatch_render_span

        tracks = ["ruler"]
        ratios = [0.35]
        if show_ref_track:
            tracks.append("reference")
            ratios.append(0.32)

        panel_track_names = []
        for i, panel in enumerate(panels):
            n_rows = max(len(panel["rows"]), 1)
            cov_name, aln_name = f"coverage_{i}", f"alignments_{i}"
            panel_track_names.append((cov_name, aln_name))
            if self.show_coverage:
                tracks.append(cov_name)
                ratios.append(1.2)
            tracks.append(aln_name)
            ratios.append(max(n_rows * ROW_HEIGHT_IN, ROW_HEIGHT_IN))

        top_margin_in = 0.5
        bottom_margin_in = 0.35 + 0.24 * self._legend_row_count()
        fig_height = sum(ratios) + top_margin_in + bottom_margin_in
        fig, axes = plt.subplots(
            nrows=len(tracks), ncols=1, figsize=(self.fig_width, fig_height), dpi=self.dpi,
            gridspec_kw={"height_ratios": ratios, "hspace": 0.2}, sharex=True,
        )
        ax_by_track = dict(zip(tracks, axes))
        for ax in axes:
            ax.set_xlim(window_start, window_end)
            for spine in ("top", "right", "left"):
                ax.spines[spine].set_visible(False)
            ax.spines["bottom"].set_visible(False)
            ax.tick_params(left=False, labelleft=False)

        tick_positions = _nice_tick_positions(window_start, window_end)

        ruler = ax_by_track["ruler"]
        ruler.set_ylim(0, 1)
        ruler.set_xticks(tick_positions)
        ruler.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{int(x):,}"))
        ruler.tick_params(labelsize=8, labelbottom=True, bottom=True, length=3, colors=AXIS_INK)
        fig.text(
            0.01, 0.995, f"{chrom}:{window_start + 1:,}-{window_end:,} ({span:,} bp)",
            fontsize=10.5, color=PRIMARY_INK, fontweight="bold", va="top", ha="left",
        )
        if suptitle:
            fig.text(0.01, 0.995 - 0.42 / fig_height, suptitle, fontsize=8.5,
                     color=SECONDARY_INK, va="top", ha="left")

        if show_ref_track:
            ref_ax = ax_by_track["reference"]
            ref_ax.set_ylim(0, 1)
            for pos in range(window_start, window_end):
                base = reference.base_at(pos) or "N"
                color = BASE_COLORS.get(base, BASE_COLORS["N"])
                ref_ax.add_patch(Rectangle((pos, 0), 1, 1, facecolor=color, alpha=0.25, edgecolor="none"))
                if show_letters:
                    ref_ax.text(pos + 0.5, 0.5, base, ha="center", va="center",
                                fontsize=7, color=color, fontweight="bold")

        for track in tracks:
            if track not in ("ruler",):
                for tick in tick_positions:
                    ax_by_track[track].axvline(tick, color=GRIDLINE, lw=0.6, zorder=0)

        for i, panel in enumerate(panels):
            rows = panel["rows"]
            layout = panel.get("layout", "pack")
            n_rows = max(len(rows), 1)
            cov_name, aln_name = panel_track_names[i]

            if self.show_coverage:
                cov_ax = ax_by_track[cov_name]
                cov_reads = panel.get("all_reads_for_coverage") or [r for row in rows for r in row]
                depth = compute_coverage(cov_reads, window_start, window_end)
                max_depth = max(depth) if depth else 0
                xs = [window_start + j + 0.5 for j in range(len(depth))]
                cov_ax.bar(xs, depth, width=1.0, color=COVERAGE_COLOR, alpha=0.85, linewidth=0)
                cov_ax.set_ylim(0, max(max_depth, 1) * 1.15)
                cov_ax.set_yticks([0, max(max_depth, 1)])
                cov_ax.tick_params(left=True, labelleft=True, labelsize=6, colors=AXIS_INK, length=3)
                cov_ax.spines["left"].set_visible(True)
                cov_ax.spines["left"].set_color(AXIS_INK)
                cov_ax.spines["left"].set_linewidth(0.8)
                cov_ax.text(0.0, 1.05, panel.get("label", f"bam{i+1}"), transform=cov_ax.transAxes,
                            fontsize=9, color=PRIMARY_INK, fontweight="bold", va="bottom")

            aln_ax = ax_by_track[aln_name]
            aln_ax.set_ylim(n_rows, 0)
            for row_idx, row in enumerate(rows):
                y0 = row_idx + ROW_MARGIN
                h = 1 - 2 * ROW_MARGIN
                for read in row:
                    self._draw_read(aln_ax, read, y0, h, render_base_detail)
                    if layout == "expand" and self.annotate_gap:
                        label = read.gap_label()
                        if label:
                            aln_ax.text(1.005, y0 + h / 2, label, transform=aln_ax.get_yaxis_transform(),
                                        fontsize=6.5, va="center", ha="left", color=SECONDARY_INK, clip_on=False)
            if not rows:
                aln_ax.text(0.5, 0.5, "No alignments in this region", transform=aln_ax.transAxes,
                            ha="center", va="center", fontsize=9, color=SECONDARY_INK)
            if not self.show_coverage:
                aln_ax.text(0.0, 1.05, panel.get("label", f"bam{i+1}"), transform=aln_ax.transAxes,
                            fontsize=9, color=PRIMARY_INK, fontweight="bold", va="bottom")

        legend_handles = self._legend_handles()
        ncol = min(len(legend_handles), 7)
        fig.legend(handles=legend_handles, loc="lower center", ncol=ncol,
                   fontsize=7, frameon=False, bbox_to_anchor=(0.5, 0.0))

        fig.subplots_adjust(left=0.05, right=0.92, top=1 - top_margin_in / fig_height,
                            bottom=bottom_margin_in / fig_height)
        fig.savefig(out_path)
        plt.close(fig)

    def _draw_read(self, ax, read: AlignedRead, y0: float, h: float, render_base_detail: bool) -> None:
        base_fill, alpha = self._read_style(read)

        if not read.blocks:
            ax.add_patch(Rectangle((read.ref_start, y0), max(read.ref_end - read.ref_start, 1), h,
                                    facecolor=base_fill, alpha=alpha, edgecolor="white", linewidth=0.3))
            return

        n_blocks = len(read.blocks)
        for i, blk in enumerate(read.blocks):
            if blk.op in ("M", "=", "X"):
                ax.add_patch(Rectangle((blk.ref_pos, y0), blk.length, h, facecolor=base_fill,
                                        alpha=alpha, edgecolor="white", linewidth=0.3, zorder=2))
            elif blk.op in ("D", "N"):
                color = DELETION_COLOR if blk.op == "D" else SKIP_COLOR
                style = "-" if blk.op == "D" else "--"
                ax.plot([blk.ref_pos, blk.ref_pos + blk.length], [y0 + h / 2, y0 + h / 2],
                        color=color, linestyle=style, linewidth=1.3, zorder=3, solid_capstyle="butt")
                if blk.length >= 3:
                    ax.text(blk.ref_pos + blk.length / 2, y0 + h / 2, f"{blk.length}",
                            fontsize=5.5, color=color, ha="center", va="bottom", zorder=4)
            elif blk.op == "I":
                width = 0.4
                ax.add_patch(Rectangle((blk.ref_pos - width / 2, y0), width, h, facecolor=INSERTION_COLOR,
                                        edgecolor="none", zorder=5))
                if blk.length >= 3:
                    ax.text(blk.ref_pos, y0, f"+{blk.length}", fontsize=5.5, color=INSERTION_COLOR,
                            ha="center", va="top", zorder=6)
            elif blk.op == "S":
                is_left = i == 0
                is_right = i == n_blocks - 1
                if not (is_left or is_right):
                    continue
                x0 = blk.ref_pos - blk.length if is_left else blk.ref_pos
                # Soft-clipped bases are real query sequence, just unaligned to the
                # reference - draw them attached to the aligned block, each base
                # colored by its own identity (same convention as mismatches),
                # rather than one flat "clip" color.
                if render_base_detail and read.query_sequence:
                    for offset in range(blk.length):
                        cbase = read.query_sequence[blk.query_pos + offset].upper()
                        ax.add_patch(Rectangle((x0 + offset, y0), 1, h,
                                                facecolor=BASE_COLORS.get(cbase, BASE_COLORS["N"]),
                                                alpha=alpha, edgecolor="none", zorder=2))
                else:
                    ax.add_patch(Rectangle((x0, y0), blk.length, h, facecolor=SOFTCLIP_COLOR,
                                            alpha=alpha, edgecolor="white", linewidth=0.3, zorder=2))
            # 'H' hard clips consume no query bases and are not drawn.

        if render_base_detail:
            for rpos, qbase in read.mismatches:
                ax.add_patch(Rectangle((rpos, y0), 1, h, facecolor=BASE_COLORS.get(qbase, BASE_COLORS["N"]),
                                        edgecolor="none", zorder=7))
