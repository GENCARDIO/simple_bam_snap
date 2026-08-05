"""Matplotlib renderer.

Replaces the old approach of shelling out to `samtools tview`, capturing its
text table, and re-parsing that text with string splits. Here every read is
drawn from its own parsed CIGAR blocks, so insertions/deletions/soft-clips/
mismatches are geometrically exact instead of guessed from column spacing.
"""
from __future__ import annotations

from dataclasses import dataclass
from math import ceil
from typing import Any, Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Polygon, Rectangle
from matplotlib.ticker import MaxNLocator

from src.annotations import (
    BAF_TRACK_FORMATS,
    CNV_TRACK_FORMATS,
    LoadedAnnotationTrack,
)
from src.config import (
    DEFAULT_ALIGNMENT_COLORS,
    DEFAULT_BASE_COLORS,
    DEFAULT_CHROMOSOME_PALETTE,
    DEFAULT_CYTOBAND_COLORS,
    DEFAULT_HAPLOTYPE_COLORS,
    DEFAULT_INSERTION_COLOR,
    DEFAULT_STYLES,
    DEFAULT_VISUAL_COLORS,
    load_config,
)
from src.cytobands import Cytoband
from src.read_model import AlignedRead
from src.reference import ReferenceWindow

# Colors: base identity follows the standard genome-browser convention
# (A green / C blue / G orange / T red); hues are drawn from the
# colorblind-validated categorical set rather than picked freehand.
BASE_COLORS = DEFAULT_BASE_COLORS
NORMAL_FILL = DEFAULT_ALIGNMENT_COLORS["normal"]
INSERTION_COLOR = DEFAULT_INSERTION_COLOR
DELETION_COLOR = DEFAULT_VISUAL_COLORS["deletion"]
SKIP_COLOR = DEFAULT_VISUAL_COLORS["reference_skip"]
SOFTCLIP_COLOR = DEFAULT_VISUAL_COLORS["softclip"]
COVERAGE_COLOR = DEFAULT_VISUAL_COLORS["coverage"]
DEFAULT_COVERAGE_VAF_THRESHOLD = 0.20
GRIDLINE = DEFAULT_VISUAL_COLORS["gridline"]
AXIS_INK = DEFAULT_VISUAL_COLORS["axis"]
PRIMARY_INK = DEFAULT_VISUAL_COLORS["primary_text"]
SECONDARY_INK = DEFAULT_VISUAL_COLORS["secondary_text"]
LEGEND_EDGE = DEFAULT_VISUAL_COLORS["legend_edge"]

# Discordant-pair fill colors, IGV-equivalent: same categories/roles IGV's
# "color by insert size and pair orientation" mode uses (red = long insert,
# blue = short insert, a blue family for same-strand pairs, green for everted,
# per-chromosome hue for inter-chromosomal mates), drawn from the same
# colorblind-validated categorical set as everything else here rather than
# IGV's undocumented internal hex values.
CHROM_PALETTE = DEFAULT_CHROMOSOME_PALETTE
HAPLOTYPE_COLORS = DEFAULT_HAPLOTYPE_COLORS

PAIR_CATEGORY_LABELS = {
    "large_insert": "Large insert",
    "small_insert": "Small insert",
    "ff": "FF (same strand)",
    "rr": "RR (same strand)",
    "everted": "Everted (RF)",
    "interchrom": "Inter-chromosomal",
}

MAPQ_ALPHA_FLOOR = DEFAULT_STYLES["mapq_alpha_floor"]

ROW_HEIGHT_IN = DEFAULT_STYLES["row_height_in"]
SQUISH_ROW_HEIGHT_IN = DEFAULT_STYLES["squish_row_height_in"]
ROW_MARGIN = DEFAULT_STYLES["row_margin"]
SQUISH_ROW_MARGIN = DEFAULT_STYLES["squish_row_margin"]
ANNOTATION_ROW_HEIGHT_IN = DEFAULT_STYLES["annotation_row_height_in"]
CNV_TRACK_HEIGHT_IN = DEFAULT_STYLES["cnv_track_height_in"]
BAF_TRACK_HEIGHT_IN = DEFAULT_STYLES["baf_track_height_in"]
IDEOGRAM_HEIGHT_IN = DEFAULT_STYLES["ideogram_height_in"]
PANEL_HEADER_HEIGHT_IN = DEFAULT_STYLES["panel_header_height_in"]
REFERENCE_HEIGHT_IN = DEFAULT_STYLES["reference_height_in"]
DEFAULT_MAX_REFERENCE_SPAN = 250
IDEOGRAM_COLOR = DEFAULT_VISUAL_COLORS["ideogram"]
IDEOGRAM_WINDOW_COLOR = DEFAULT_VISUAL_COLORS["ideogram_window"]
CYTOBAND_COLORS = DEFAULT_CYTOBAND_COLORS
CNV_GAIN_COLOR = DEFAULT_VISUAL_COLORS["cnv_gain"]
CNV_LOSS_COLOR = DEFAULT_VISUAL_COLORS["cnv_loss"]


@dataclass
class SnvEvidence:
    count: int = 0
    forward: int = 0
    reverse: int = 0
    base_quality_sum: int = 0
    mapq_sum: int = 0


def chrom_color(chrom: Optional[str], palette: Optional[List[str]] = None) -> str:
    """Stable per-chromosome hue (same idea as IGV's karyotype coloring) for
    inter-chromosomal mate pairs. Deterministic within a run, not globally
    fixed across all human chromosome names."""
    palette = palette or CHROM_PALETTE
    if not chrom:
        return AXIS_INK
    h = 0
    for ch in chrom:
        h = (h * 31 + ord(ch)) & 0xFFFFFFFF
    return palette[h % len(palette)]


def haplotype_color(
    haplotype: Optional[str], colors: Optional[Dict[str, str]] = None,
    palette: Optional[List[str]] = None,
) -> str:
    colors = colors or HAPLOTYPE_COLORS
    palette = palette or CHROM_PALETTE
    label = str(haplotype) if haplotype is not None else "untagged"
    if label in colors:
        return colors[label]
    if label.isdigit():
        return palette[(int(label) - 1) % len(palette)]
    value = 0
    for character in label:
        value = (value * 31 + ord(character)) & 0xFFFFFFFF
    return palette[value % len(palette)]


def compute_coverage(reads: List[AlignedRead], start: int, end: int) -> List[int]:
    """Per-base depth across [start, end), counting only reference-consuming
    match bases (a deleted base is not "covered", matching `samtools depth`)."""
    span = max(0, end - start)
    changes = [0] * (span + 1)
    for read in reads:
        for blk in read.blocks:
            if blk.op not in ("M", "=", "X"):
                continue
            lo = max(blk.ref_pos, start)
            hi = min(blk.ref_pos + blk.length, end)
            if lo < hi:
                changes[lo - start] += 1
                changes[hi - start] -= 1
    depth = []
    running = 0
    for change in changes[:-1]:
        running += change
        depth.append(running)
    return depth


def compute_snv_counts(
    reads: List[AlignedRead], start: int, end: int
) -> Dict[int, Dict[str, int]]:
    """Count reference-backed A/C/G/T mismatches at each covered position."""
    counts: Dict[int, Dict[str, int]] = {}
    for read in reads:
        for position, base in read.mismatches:
            if start <= position < end and base in BASE_COLORS and base != "N":
                position_counts = counts.setdefault(position, {})
                position_counts[base] = position_counts.get(base, 0) + 1
    return counts


def compute_snv_evidence(
    reads: List[AlignedRead], start: int, end: int,
    min_baseq: int = 0, min_mapq: int = 0,
) -> Tuple[List[int], Dict[int, Dict[str, SnvEvidence]]]:
    """Return quality-filtered nucleotide depth and alternate-SNV evidence."""
    depth = [0] * max(0, end - start)
    evidence: Dict[int, Dict[str, SnvEvidence]] = {}
    for read in reads:
        if getattr(read, "mapq", 0) < min_mapq:
            continue
        sequence = getattr(read, "query_sequence", "") or ""
        qualities = getattr(read, "query_qualities", []) or []
        for block in read.blocks:
            if block.op not in ("M", "=", "X"):
                continue
            lo = max(block.ref_pos, start)
            hi = min(block.ref_pos + block.length, end)
            for position in range(lo, hi):
                query_index = block.query_pos + position - block.ref_pos
                base = sequence[query_index].upper() if query_index < len(sequence) else "N"
                base_quality = qualities[query_index] if query_index < len(qualities) else 0
                if base in "ACGT" and base_quality >= min_baseq:
                    depth[position - start] += 1

        details = getattr(read, "mismatch_details", None)
        if details is None:
            details = []
            for position, base in getattr(read, "mismatches", []):
                details.append((position, base, 0))
        for position, base, base_quality in details:
            if not start <= position < end or base not in "ACGT" or base_quality < min_baseq:
                continue
            base_evidence = evidence.setdefault(position, {}).setdefault(base, SnvEvidence())
            base_evidence.count += 1
            if getattr(read, "is_reverse", False):
                base_evidence.reverse += 1
            else:
                base_evidence.forward += 1
            base_evidence.base_quality_sum += base_quality
            base_evidence.mapq_sum += getattr(read, "mapq", 0)
    return depth, evidence


def nice_tick_positions(start: int, end: int, target: int = 8) -> List[int]:
    locator = MaxNLocator(nbins=target, steps=[1, 2, 5, 10])
    ticks = []
    for value in locator.tick_values(start, end):
        if start <= value <= end:
            ticks.append(int(value))
    return ticks


def left_margin_fraction(
    fig_width: float, genomic_tracks: List[LoadedAnnotationTrack]
) -> float:
    """Reserve enough physical space for annotation labels outside the axes."""
    longest_label = max((len(track.label) for track in genomic_tracks), default=0)
    margin_in = max(0.70, 0.25 + longest_label * 0.065)
    return min(margin_in / fig_width, 0.22)


def ellipsize(text: str, max_chars: int) -> str:
    """Fit a label to an approximate character budget without spilling."""
    if max_chars <= 0:
        return ""
    if len(text) <= max_chars:
        return text
    if max_chars == 1:
        return "…"
    return text[:max_chars - 1].rstrip() + "…"


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
        alignment_colors: Optional[Dict[str, Optional[str]]] = None,
        display_mode: str = "expand",
        show_ideogram: bool = True,
        max_reference_span: int = DEFAULT_MAX_REFERENCE_SPAN,
        view_as_pairs: bool = False,
        coverage_vaf_threshold: float = DEFAULT_COVERAGE_VAF_THRESHOLD,
        min_baseq: int = 0,
        min_variant_mapq: int = 0,
        show_variant_counts: bool = False,
        haplotype_view: str = "none",
        visual_config: Optional[Dict[str, Any]] = None,
        sort_base_position: Optional[int] = None,
        sort_reference_base: Optional[str] = None,
        show_center_guide: bool = False,
    ):
        theme = visual_config or load_config()
        self.base_colors = dict(theme["base_colors"])
        self.visual_colors = dict(theme["visual_colors"])
        self.haplotype_colors = dict(theme["haplotype_colors"])
        self.cytoband_colors = dict(theme["cytoband_colors"])
        self.chromosome_palette = list(theme["chromosome_palette"])
        self.styles = dict(theme["styles"])
        self.sort_base_position = sort_base_position
        self.sort_reference_base = sort_reference_base
        self.active_sort_base_position = sort_base_position
        self.active_sort_reference_base = sort_reference_base
        self.show_center_guide = show_center_guide
        self.fig_width = fig_width
        self.dpi = dpi
        self.show_coverage = show_coverage
        self.annotate_gap = annotate_gap
        self.max_mismatch_render_span = max_mismatch_render_span
        self.pair_colors = pair_colors
        self.shade_by_mapq = shade_by_mapq
        self.mapq_cap = mapq_cap
        self.alignment_colors = dict(theme["alignment_colors"])
        if alignment_colors:
            self.alignment_colors.update(alignment_colors)
        if display_mode not in ("collapse", "expand", "squish"):
            raise ValueError(
                f"Unknown display mode '{display_mode}'. Choose collapse, expand, or squish."
            )
        self.display_mode = display_mode
        self.show_ideogram = show_ideogram
        self.max_reference_span = max_reference_span
        self.view_as_pairs = view_as_pairs
        if not 0 <= coverage_vaf_threshold <= 1:
            raise ValueError("Coverage VAF threshold must be between 0 and 1.")
        self.coverage_vaf_threshold = coverage_vaf_threshold
        if min_baseq < 0 or min_variant_mapq < 0:
            raise ValueError("Variant base-quality and MAPQ filters cannot be negative.")
        self.min_baseq = min_baseq
        self.min_variant_mapq = min_variant_mapq
        self.show_variant_counts = show_variant_counts
        if haplotype_view not in ("none", "color", "split"):
            raise ValueError("Haplotype view must be none, color, or split.")
        self.haplotype_view = haplotype_view
        self.row_height_in = (
            self.styles["squish_row_height_in"] if display_mode == "squish"
            else self.styles["row_height_in"]
        )
        self.row_margin = (
            self.styles["squish_row_margin"] if display_mode == "squish"
            else self.styles["row_margin"]
        )
        if fig_width >= 9:
            self.legend_height_in = 0.78
        elif fig_width >= 6:
            self.legend_height_in = 1.90 if pair_colors or haplotype_view != "none" else 1.15
        else:
            self.legend_height_in = 2.75 if pair_colors or haplotype_view != "none" else 1.55
        self.legend_bottom_in = 0.04
        self.legend_plot_gap_in = 0.12
        self.legend_tick_clearance_in = 0.30
        self.legend_margin_in = (
            self.legend_bottom_in + self.legend_height_in + self.legend_plot_gap_in
            + self.legend_tick_clearance_in
        )

    def read_style(self, read: AlignedRead):
        """(fill_color, alpha) for a read's main body: hue encodes pair
        discordance category (when enabled), alpha encodes mapping quality -
        low-MAPQ reads get a lighter/more washed-out fill, same idea as IGV's
        "shade by mapping quality"."""
        if self.haplotype_view in ("color", "split"):
            color = haplotype_color(
                getattr(read, "haplotype", None), self.haplotype_colors,
                self.chromosome_palette,
            )
        elif self.pair_colors and read.pair_category == "interchrom":
            color = self.alignment_colors["interchrom"] or chrom_color(
                read.mate_chrom, self.chromosome_palette
            )
        elif self.pair_colors and read.pair_category == "large_insert":
            color = self.alignment_colors["large_insert"]
        elif self.pair_colors and read.pair_category == "small_insert":
            color = self.alignment_colors["small_insert"]
        elif self.pair_colors and read.pair_category == "ff":
            color = self.alignment_colors["ff"]
        elif self.pair_colors and read.pair_category == "rr":
            color = self.alignment_colors["rr"]
        elif self.pair_colors and read.pair_category == "everted":
            color = self.alignment_colors["everted"]
        else:
            color = self.alignment_colors["normal"]

        alpha = (
            self.styles["secondary_alignment_alpha"]
            if (read.is_secondary or read.is_duplicate)
            else self.styles["alignment_alpha"]
        )
        if self.shade_by_mapq and self.mapq_cap > 0:
            mapq_frac = min(max(read.mapq, 0), self.mapq_cap) / self.mapq_cap
            floor = self.styles["mapq_alpha_floor"]
            alpha *= floor + (1 - floor) * mapq_frac
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
        downsampled_reads: int = 0,
        all_reads_for_coverage: Optional[List[AlignedRead]] = None,
        genomic_tracks: Optional[List[LoadedAnnotationTrack]] = None,
        contig_length: Optional[int] = None,
        cytobands: Optional[List[Cytoband]] = None,
    ) -> None:
        span = window_end - window_start
        n_rows = max(len(rows), 1)
        show_ref_track = bool(
            reference and reference.available and
            self.max_reference_span > 0 and span <= self.max_reference_span
        )
        render_base_detail = span <= self.max_mismatch_render_span

        tracks = []
        ratios = []
        if self.show_ideogram and contig_length:
            tracks.append("ideogram")
            ratios.append(self.styles["ideogram_height_in"])
        if show_ref_track:
            tracks.append("reference")
            ratios.append(self.styles["reference_height_in"])
        genomic_tracks = genomic_tracks or []
        for index, annotation in enumerate(genomic_tracks):
            tracks.append(f"annotation_{index}")
            if annotation.kind in CNV_TRACK_FORMATS:
                ratios.append(self.styles["cnv_track_height_in"])
            elif annotation.kind in BAF_TRACK_FORMATS:
                ratios.append(self.styles["baf_track_height_in"])
            else:
                ratios.append(
                    max(len(annotation.rows), 1) * self.styles["annotation_row_height_in"]
                )
        if self.show_coverage:
            tracks.append("coverage")
            ratios.append(1.4)
        tracks.append("alignments")
        ratios.append(max(n_rows * self.row_height_in, self.row_height_in))

        top_margin_in = 0.72  # dedicated region-title and subtitle lanes
        bottom_margin_in = self.legend_margin_in
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

        tick_positions = nice_tick_positions(window_start, window_end)

        fig.text(
            0.01, 1 - 0.06 / fig_height,
            f"{chrom}:{window_start + 1:,}-{window_end:,} ({span:,} bp)",
            fontsize=10.5, color=self.visual_colors["primary_text"], fontweight="bold", va="top", ha="left",
        )
        subtitle = title
        if dropped_reads:
            subtitle = (subtitle + " -- " if subtitle else "") + (
                f"{dropped_reads} lower-priority read(s) not shown (--max_rows)"
            )
        if downsampled_reads:
            subtitle = (subtitle + " -- " if subtitle else "") + (
                f"{downsampled_reads} alignment(s) downsampled"
            )
        if subtitle:
            fig.text(
                0.01, 1 - 0.34 / fig_height,
                ellipsize(subtitle, max(30, int(self.fig_width * 15))),
                fontsize=8.5, color=self.visual_colors["secondary_text"], va="top", ha="left",
            )

        if "ideogram" in ax_by_track:
            self.draw_ideogram(
                ax_by_track["ideogram"], chrom, window_start, window_end, contig_length,
                cytobands,
            )

        # --- reference -----------------------------------------------------
        if show_ref_track:
            self.draw_reference_track(
                ax_by_track["reference"], reference, window_start, window_end,
                available_width_in=self.fig_width,
            )

        for index, annotation in enumerate(genomic_tracks):
            self.draw_annotation_track(
                ax_by_track[f"annotation_{index}"], annotation, window_start, window_end
            )

        # --- coverage --------------------------------------------------
        if self.show_coverage:
            cov_ax = ax_by_track["coverage"]
            cov_reads = all_reads_for_coverage
            if cov_reads is None:
                cov_reads = []
                for row in rows:
                    cov_reads.extend(row)
            self.draw_coverage_track(cov_ax, cov_reads, window_start, window_end)

        # --- alignments --------------------------------------------------
        aln_ax = ax_by_track["alignments"]
        aln_ax.set_ylim(n_rows, 0)
        self.draw_haplotype_lanes(aln_ax, rows)
        aln_ax.set_xticks(tick_positions)
        aln_ax.tick_params(
            bottom=True, labelbottom=True, labelsize=9,
            length=3, colors=self.visual_colors["primary_text"],
        )
        for tick in tick_positions:
            for track in tracks:
                if track == "ideogram":
                    continue
                if track in ax_by_track:
                    ax_by_track[track].axvline(
                        tick, color=self.visual_colors["gridline"],
                        lw=self.styles["grid_line_width"], zorder=0,
                    )

        for row_idx, row in enumerate(rows):
            y0 = row_idx + self.row_margin
            h = 1 - 2 * self.row_margin
            self.draw_alignment_row(
                aln_ax, row, y0, h, render_base_detail, layout
            )

        if not rows:
            aln_ax.text(
                0.5, 0.5, "No alignments in this region", transform=aln_ax.transAxes,
                ha="center", va="center", fontsize=10, color=self.visual_colors["secondary_text"],
            )

        for track, ax in ax_by_track.items():
            if track != "ideogram":
                self.draw_center_guide(ax, window_start, window_end)

        # --- legend -----------------------------------
        plot_left = left_margin_fraction(self.fig_width, genomic_tracks)
        if self.haplotype_view == "split":
            plot_left = max(plot_left, min(1.15 / self.fig_width, 0.25))
        plot_right = 0.92
        fig.subplots_adjust(left=plot_left, right=plot_right,
                            top=1 - top_margin_in / fig_height,
                            bottom=bottom_margin_in / fig_height)
        legends = self.draw_legends(fig, fig_height, plot_left, plot_right)
        self.separate_legend_from_plots(fig, axes, legends)
        fig.savefig(out_path)
        plt.close(fig)

    def draw_coverage_track(
        self, ax, reads: List[AlignedRead], start: int, end: int
    ) -> None:
        """Draw depth with qualifying SNV allele fractions stacked in base colours."""
        depth = compute_coverage(reads, start, end)
        max_depth = max(depth) if depth else 0
        positions = []
        for index in range(len(depth)):
            positions.append(start + index + 0.5)
        ax.bar(
            positions, depth, width=1.0, color=self.visual_colors["coverage"],
            alpha=self.styles["coverage_alpha"], linewidth=0,
        )

        evidence_depth, snv_evidence = compute_snv_evidence(
            reads, start, end,
            min_baseq=self.min_baseq, min_mapq=self.min_variant_mapq,
        )
        labels = []
        for position, base_counts in snv_evidence.items():
            position_depth = evidence_depth[position - start]
            if position_depth <= 0:
                continue
            bottom = 0
            for base in "ACGT":
                allele = base_counts.get(base)
                if allele is None or allele.count / position_depth <= self.coverage_vaf_threshold:
                    continue
                ax.bar(
                    position + 0.5, allele.count, width=1.0, bottom=bottom,
                    color=self.base_colors[base], linewidth=0, zorder=2,
                )
                bottom += allele.count
                labels.append((position, base, allele, position_depth))

        ax.set_ylim(0, max(max_depth, 1) * 1.15)
        ax.set_yticks([0, max(max_depth, 1)])
        ax.tick_params(left=True, labelleft=True, labelsize=6, colors=self.visual_colors["axis"], length=3)
        ax.spines["left"].set_visible(True)
        ax.spines["left"].set_color(self.visual_colors["axis"])
        ax.spines["left"].set_linewidth(0.8)
        ax.text(
            0.005, 0.98, "coverage", transform=ax.transAxes, fontsize=7,
            color=self.visual_colors["secondary_text"], va="top",
        )

        base_spacing_px = ax.get_window_extent().width / max(end - start, 1)
        if self.show_variant_counts and base_spacing_px >= 5.5:
            for position, base, allele, position_depth in labels:
                mean_baseq = allele.base_quality_sum / allele.count
                mean_mapq = allele.mapq_sum / allele.count
                label = (
                    f"{base} {allele.count}/{position_depth} "
                    f"{allele.count / position_depth:.0%} "
                    f"F{allele.forward}/R{allele.reverse} "
                    f"BQ{mean_baseq:.0f} MQ{mean_mapq:.0f}"
                )
                ax.text(
                    position + 0.5, 0.97, label,
                    transform=ax.get_xaxis_transform(), rotation=90,
                    ha="right", va="center", fontsize=4.8,
                    color=self.visual_colors["primary_text"], clip_on=True, zorder=4,
                    bbox={
                        "facecolor": self.visual_colors["label_background"],
                        "edgecolor": "none", "alpha": 0.68, "pad": 0.15,
                    },
                )

    def draw_center_guide(self, ax, start: int, end: int) -> None:
        """Draw the optional IGV-like guide at the exact window midpoint."""
        if not self.show_center_guide:
            return
        ax.axvline(
            start + (end - start) / 2,
            color=self.visual_colors["center_guide"],
            alpha=self.styles["center_guide_alpha"],
            linewidth=self.styles["center_guide_width"],
            linestyle=self.styles["center_guide_line_style"],
            zorder=20,
        )

    def draw_ideogram(
        self,
        ax,
        chrom: str,
        window_start: int,
        window_end: int,
        contig_length: int,
        cytobands: Optional[List[Cytoband]] = None,
    ) -> None:
        """Draw a UCSC-style cytoband ideogram with the current window in red."""
        ax.set_ylim(0, 1)
        ax.set_xlim(window_start, window_end)
        # Use the genomic x-axis itself, not figure-relative coordinates. This
        # locks both chromosome ends to the exact plot boundaries in single,
        # comparison, and mate layouts even when their margins differ.
        chromosome_transform = ax.get_xaxis_transform()
        bar_x = window_start
        bar_y = 0.28
        bar_width = max(window_end - window_start, 1)
        bar_height = 0.44
        chromosome_vertices = None
        neck_position = None
        middle_y = None
        neck_half_height = None
        p_centromeres = []
        q_centromeres = []
        for band in cytobands or []:
            if band.stain != "acen":
                continue
            if band.name.startswith("p"):
                p_centromeres.append(band)
            elif band.name.startswith("q"):
                q_centromeres.append(band)
        if p_centromeres and q_centromeres:
            p_centromere = max(p_centromeres, key=lambda band: band.end)
            q_centromere = min(q_centromeres, key=lambda band: band.start)
            p_shoulder = bar_x + p_centromere.start / contig_length * bar_width
            q_shoulder = bar_x + q_centromere.end / contig_length * bar_width
            neck_position = (
                p_centromere.end + q_centromere.start
            ) / 2 / contig_length * bar_width + bar_x
            middle_y = bar_y + bar_height / 2
            neck_half_height = bar_height * 0.22
            chromosome_vertices = [
                (bar_x, bar_y),
                (p_shoulder, bar_y),
                (neck_position, middle_y - neck_half_height),
                (q_shoulder, bar_y),
                (bar_x + bar_width, bar_y),
                (bar_x + bar_width, bar_y + bar_height),
                (q_shoulder, bar_y + bar_height),
                (neck_position, middle_y + neck_half_height),
                (p_shoulder, bar_y + bar_height),
                (bar_x, bar_y + bar_height),
            ]
            chromosome_clip = Polygon(
                chromosome_vertices, closed=True,
                transform=chromosome_transform, facecolor=self.visual_colors["ideogram"],
                edgecolor=self.visual_colors["axis"], linewidth=0.6, zorder=2,
            )
        else:
            chromosome_clip = Rectangle(
                (bar_x, bar_y), bar_width, bar_height,
                transform=chromosome_transform, facecolor=self.visual_colors["ideogram"],
                edgecolor=self.visual_colors["axis"], linewidth=0.6, zorder=2,
            )
        ax.add_patch(chromosome_clip)

        for band in cytobands or []:
            band_start = min(max(band.start, 0), contig_length)
            band_end = min(max(band.end, band_start), contig_length)
            if band_end <= band_start:
                continue
            x0 = bar_x + band_start / contig_length * bar_width
            x1 = bar_x + band_end / contig_length * bar_width
            color = self.cytoband_colors.get(band.stain, self.visual_colors["ideogram"])
            patch = Rectangle(
                (x0, bar_y), x1 - x0, bar_height,
                transform=chromosome_transform, facecolor=color,
                edgecolor=self.visual_colors["cytoband_edge"], linewidth=0.15, zorder=2.1,
            )
            patch.set_clip_path(chromosome_clip)
            ax.add_patch(patch)

        if chromosome_vertices:
            bridge_width = max(bar_width * 0.004, 0.4)
            centromere_bridge = Rectangle(
                (neck_position - bridge_width / 2, middle_y - neck_half_height),
                bridge_width, neck_half_height * 2,
                transform=chromosome_transform,
                facecolor=self.visual_colors["centromere"], edgecolor="none", zorder=2.15,
            )
            centromere_bridge.set_clip_path(chromosome_clip)
            ax.add_patch(centromere_bridge)

        if cytobands:
            if chromosome_vertices:
                outline = Polygon(
                    chromosome_vertices, closed=True,
                    transform=chromosome_transform, facecolor="none",
                    edgecolor=self.visual_colors["axis"], linewidth=0.6, zorder=2.2,
                )
            else:
                outline = Rectangle(
                    (bar_x, bar_y), bar_width, bar_height,
                    transform=chromosome_transform, facecolor="none",
                    edgecolor=self.visual_colors["axis"], linewidth=0.6, zorder=2.2,
                )
            ax.add_patch(outline)

        clamped_start = min(max(window_start, 0), contig_length)
        clamped_end = min(max(window_end, clamped_start), contig_length)
        relative_start = clamped_start / contig_length
        relative_width = max((clamped_end - clamped_start) / contig_length, 0.0)
        # Base-pair windows on chromosome-scale bars would otherwise disappear.
        marker_width = max(relative_width * bar_width, 0.004)
        marker_center = bar_x + (relative_start + relative_width / 2) * bar_width
        marker_x = min(max(marker_center - marker_width / 2, bar_x), bar_x + bar_width - marker_width)
        ax.add_patch(Rectangle(
            (marker_x, bar_y - 0.08), marker_width, bar_height + 0.16,
            transform=chromosome_transform, facecolor=self.visual_colors["ideogram_window"],
            edgecolor=self.visual_colors["contrast_edge"], linewidth=0.35, zorder=3,
        ))
        ax.text(
            -0.012, 0.5, chrom, transform=ax.transAxes, ha="right", va="center",
            fontsize=7, color=self.visual_colors["primary_text"], fontweight="bold",
            clip_on=False,
        )
        ax.text(
            1.012, 0.5, f"{contig_length / 1_000_000:.1f} Mb",
            transform=ax.transAxes, ha="left", va="center", fontsize=6.5,
            color=self.visual_colors["secondary_text"], clip_on=False,
        )
        ax.text(
            marker_x + marker_width / 2, 0.93, "window",
            transform=chromosome_transform,
            ha="center", va="bottom", fontsize=5.5, color=self.visual_colors["ideogram_window"],
            clip_on=False,
        )

    def draw_annotation_track(
        self,
        ax,
        track: LoadedAnnotationTrack,
        window_start: int,
        window_end: int,
        shared_row_count: Optional[int] = None,
    ) -> None:
        """Draw a UCSC-like BED or transcript annotation track."""
        if track.kind in CNV_TRACK_FORMATS:
            self.draw_cnv_track(ax, track, window_start, window_end)
            return
        if track.kind in BAF_TRACK_FORMATS:
            self.draw_baf_track(ax, track, window_start, window_end)
            return
        row_count = max(shared_row_count or len(track.rows), 1)
        ax.set_ylim(row_count, 0)
        margin_in = min(
            max(0.70, 0.25 + len(track.label) * 0.065), self.fig_width * 0.22
        )
        label_capacity = max(5, int((margin_in - 0.20) / 0.065))
        ax.text(
            -0.012, 0.5, ellipsize(track.label, label_capacity),
            transform=ax.transAxes, ha="right", va="center",
            fontsize=7, color=track.color, fontweight="bold", clip_on=False,
        )
        if not track.rows:
            ax.text(
                0.01, 0.5, "No features", transform=ax.transAxes, ha="left", va="center",
                fontsize=6.5, color=self.visual_colors["axis"],
            )
            return

        for row_index, row in enumerate(track.rows):
            # Reserve the upper quarter of each row for the feature name.
            center = row_index + 0.64
            for item in row:
                line_start = max(item.start, window_start)
                line_end = min(item.end, window_end)
                if line_start >= line_end:
                    continue
                ax.plot(
                    [line_start, line_end], [center, center], color=track.color,
                    linewidth=(
                        self.styles["primary_gene_line_width"] if item.primary_rank is not None
                        else self.styles["gene_line_width"]
                    ),
                    zorder=2, solid_capstyle="butt",
                )
                for block_start, block_end in item.blocks:
                    lo, hi = max(block_start, window_start), min(block_end, window_end)
                    if lo < hi:
                        ax.add_patch(Rectangle(
                            (lo, center - 0.23), hi - lo, 0.46,
                            facecolor=track.color, edgecolor=track.color,
                            linewidth=self.styles["alignment_edge_width"], zorder=3,
                        ))
                for utr_start, utr_end in item.utrs:
                    lo, hi = max(utr_start, window_start), min(utr_end, window_end)
                    if lo < hi:
                        ax.add_patch(Rectangle(
                            (lo, center - 0.12), hi - lo, 0.24,
                            facecolor=track.color, edgecolor=track.color,
                            linewidth=self.styles["alignment_edge_width"], zorder=3,
                        ))

                # Repeated small arrows on introns make transcript direction
                # readable without competing with exon blocks.
                merged_intervals: List[List[int]] = []
                for feature_start, feature_end in sorted(item.blocks + item.utrs):
                    if feature_end <= feature_start:
                        continue
                    if not merged_intervals or feature_start > merged_intervals[-1][1]:
                        merged_intervals.append([feature_start, feature_end])
                    else:
                        merged_intervals[-1][1] = max(
                            merged_intervals[-1][1], feature_end
                        )
                introns = []
                for left, right in zip(merged_intervals, merged_intervals[1:]):
                    if left[1] < right[0]:
                        introns.append((left[1], right[0]))
                if item.strand in ("+", "-") and introns:
                    axes_width_px = max(ax.get_window_extent().width, 1)
                    arrow_spacing = max(
                        (window_end - window_start) / max(axes_width_px / 13.0, 1),
                        1.0,
                    )
                    marker = ">" if item.strand == "+" else "<"
                    for intron_start, intron_end in introns:
                        lo = max(intron_start, window_start)
                        hi = min(intron_end, window_end)
                        if lo >= hi:
                            continue
                        position = lo + arrow_spacing / 2
                        while position < hi:
                            ax.plot(
                                position, center, marker=marker, markersize=2.2,
                                color=track.color, markeredgewidth=0, zorder=4,
                            )
                            position += arrow_spacing
                elif item.strand == "+" and item.end <= window_end:
                    ax.plot(item.end, center, marker=">", markersize=2.5,
                            color=track.color, markeredgewidth=0, zorder=4)
                elif item.strand == "-" and item.start >= window_start:
                    ax.plot(item.start, center, marker="<", markersize=2.5,
                            color=track.color, markeredgewidth=0, zorder=4)

                visible_fraction = (line_end - line_start) / max(window_end - window_start, 1)
                name_capacity = int(visible_fraction * 105)
                item_label = (
                    item.transcript_label
                    if track.display_mode != "collapse" and item.transcript_label
                    else item.name
                )
                if item.primary_rank is not None:
                    item_label += " ★"
                display_name = ellipsize(item_label, name_capacity) if name_capacity >= 4 else ""
                if display_name:
                    ax.text(
                        line_start, row_index + 0.05, display_name, ha="left", va="top",
                        fontsize=5.5, color=track.color, clip_on=True, zorder=5,
                    )

    def draw_cnv_track(
        self, ax, track: LoadedAnnotationTrack, window_start: int, window_end: int
    ) -> None:
        """Draw segmented or binned log2 copy-number values around a zero baseline."""
        largest_value = max(
            (abs(item.value) for item in track.items if item.value is not None),
            default=0.0,
        )
        value_limit = max(0.5, ceil(largest_value * 1.15 * 2) / 2)
        ax.set_ylim(-value_limit, value_limit)
        ax.axhline(0, color=self.visual_colors["axis"], linewidth=0.65, zorder=1)

        use_sign_colors = track.color_by_sign
        for item in track.items:
            if item.value is None:
                continue
            lo, hi = max(item.start, window_start), min(item.end, window_end)
            if lo >= hi:
                continue
            if use_sign_colors:
                color = (
                    self.visual_colors["cnv_gain"] if item.value > 0
                    else self.visual_colors["cnv_loss"] if item.value < 0
                    else track.color
                )
            else:
                color = track.color
            ax.add_patch(Rectangle(
                (lo, min(0, item.value)), hi - lo, abs(item.value),
                facecolor=color, edgecolor="none",
                alpha=self.styles["cnv_fill_alpha"], zorder=2,
            ))
            ax.plot(
                [lo, hi], [item.value, item.value], color=color,
                linewidth=1.25, solid_capstyle="butt", zorder=3,
            )

        margin_in = min(
            max(0.70, 0.25 + len(track.label) * 0.065), self.fig_width * 0.22
        )
        label_capacity = max(5, int((margin_in - 0.20) / 0.065))
        ax.text(
            -0.012, 0.5, ellipsize(track.label, label_capacity),
            transform=ax.transAxes, ha="right", va="center", fontsize=7,
            color=track.color, fontweight="bold", clip_on=False,
        )
        sample_names = list(dict.fromkeys(
            item.sample for item in track.items if item.sample
        ))
        if sample_names:
            sample_label = ", ".join(sample_names[:3])
            if len(sample_names) > 3:
                sample_label += f" +{len(sample_names) - 3}"
            ax.text(
                0.005, 0.97, ellipsize(sample_label, max(12, int(self.fig_width * 8))),
                transform=ax.transAxes, ha="left", va="top", fontsize=6,
                color=self.visual_colors["secondary_text"], clip_on=True,
            )
        elif not track.items:
            ax.text(
                0.01, 0.5, "No CNV data", transform=ax.transAxes,
                ha="left", va="center", fontsize=6.5, color=self.visual_colors["axis"],
            )
        ax.set_yticks([-value_limit, 0, value_limit])
        ax.set_yticklabels([
            f"{-value_limit:g}", "0", f"{value_limit:g}",
        ])
        ax.tick_params(
            right=True, labelright=True, left=False, labelleft=False,
            labelsize=6, colors=self.visual_colors["axis"], length=2,
        )
        ax.spines["right"].set_visible(True)
        ax.spines["right"].set_color(self.visual_colors["axis"])
        ax.spines["right"].set_linewidth(0.6)
        ax.yaxis.set_label_position("right")
        ax.set_ylabel(
            "log2", rotation=90, labelpad=18, fontsize=5.5,
            color=self.visual_colors["secondary_text"], va="center",
        )

    def draw_baf_track(
        self, ax, track: LoadedAnnotationTrack, window_start: int, window_end: int
    ) -> None:
        """Draw heterozygous-SNV B-allele fractions on a zero-to-one scale."""
        ax.set_ylim(0, 1)
        ax.axhline(0.5, color=self.visual_colors["axis"], linewidth=0.65, linestyle="--", zorder=1)
        positions = []
        values = []
        for item in track.items:
            if item.value is not None:
                positions.append(item.start + 0.5)
                values.append(item.value)
        if positions:
            ax.scatter(
                positions, values, s=11, color=track.color,
                edgecolors=self.visual_colors["contrast_edge"], linewidths=0.25,
                alpha=self.styles["baf_alpha"], zorder=3,
            )

        margin_in = min(
            max(0.70, 0.25 + len(track.label) * 0.065), self.fig_width * 0.22
        )
        label_capacity = max(5, int((margin_in - 0.20) / 0.065))
        ax.text(
            -0.012, 0.5, ellipsize(track.label, label_capacity),
            transform=ax.transAxes, ha="right", va="center", fontsize=7,
            color=track.color, fontweight="bold", clip_on=False,
        )
        sample_names = list(dict.fromkeys(
            item.sample for item in track.items if item.sample
        ))
        if sample_names:
            ax.text(
                0.005, 0.97, ellipsize(", ".join(sample_names), max(12, int(self.fig_width * 8))),
                transform=ax.transAxes, ha="left", va="top", fontsize=6,
                color=self.visual_colors["secondary_text"], clip_on=True,
            )
        elif not track.items:
            ax.text(
                0.01, 0.5, "No heterozygous SNPs with AD/AF",
                transform=ax.transAxes, ha="left", va="center",
                fontsize=6.5, color=self.visual_colors["axis"],
            )
        ax.set_yticks([0, 0.5, 1])
        ax.tick_params(
            right=True, labelright=True, left=False, labelleft=False,
            labelsize=6, colors=self.visual_colors["axis"], length=2,
        )
        ax.spines["right"].set_visible(True)
        ax.spines["right"].set_color(self.visual_colors["axis"])
        ax.spines["right"].set_linewidth(0.6)
        ax.yaxis.set_label_position("right")
        ax.set_ylabel(
            "BAF", rotation=90, labelpad=18, fontsize=5.5,
            color=self.visual_colors["secondary_text"], va="center",
        )

    def draw_reference_track(
        self,
        ax,
        reference: ReferenceWindow,
        window_start: int,
        window_end: int,
        available_width_in: float,
    ) -> None:
        """Draw one lightly coloured cell per FASTA base, with letters when legible."""
        ax.set_ylim(0, 1)
        span = max(window_end - window_start, 1)
        show_letters = available_width_in * 72 / span >= 6.2
        for pos in range(window_start, window_end):
            base = reference.base_at(pos) or "N"
            color = self.base_colors.get(base, self.base_colors["N"])
            ax.add_patch(Rectangle(
                (pos, 0.08), 1, 0.84, facecolor=color,
                alpha=self.styles["reference_base_alpha"],
                edgecolor=self.visual_colors["contrast_edge"], linewidth=0.25, zorder=2,
            ))
            if show_letters:
                ax.text(
                    pos + 0.5, 0.5, base, ha="center", va="center",
                    fontsize=7, color=color, fontweight="bold", zorder=3,
                    clip_on=True,
                )
        ax.text(
            0.005, 0.98, "reference", transform=ax.transAxes,
            fontsize=7, color=self.visual_colors["secondary_text"], va="top", zorder=4,
        )

    def draw_legends(
        self, fig, fig_height: float, plot_left: float = 0.05,
        plot_right: float = 0.95,
    ) -> list:
        """Draw one plot-aligned legend rail divided into topic compartments."""
        alignment_handles = [
            Patch(facecolor=self.alignment_colors["normal"], edgecolor="none", label="Normal / concordant"),
            Patch(facecolor=self.alignment_colors["small_insert"], edgecolor="none", label="Insertion"),
            Line2D([0], [0], color=self.visual_colors["deletion"], lw=1.5, label="Deletion"),
        ]
        if self.view_as_pairs:
            alignment_handles.append(
                Line2D([0], [0], color=self.visual_colors["secondary_text"], lw=1.0, label="Mate link")
            )
        pair_handles = []
        if self.pair_colors and self.haplotype_view == "none":
            pair_handles = [
                Patch(facecolor=self.alignment_colors["large_insert"], edgecolor="none", label=PAIR_CATEGORY_LABELS["large_insert"]),
                Patch(facecolor=self.alignment_colors["small_insert"], edgecolor="none", label=PAIR_CATEGORY_LABELS["small_insert"]),
                Patch(facecolor=self.alignment_colors["ff"], edgecolor="none", label=PAIR_CATEGORY_LABELS["ff"]),
                Patch(facecolor=self.alignment_colors["rr"], edgecolor="none", label=PAIR_CATEGORY_LABELS["rr"]),
                Patch(facecolor=self.alignment_colors["everted"], edgecolor="none", label=PAIR_CATEGORY_LABELS["everted"]),
                Patch(
                    facecolor=self.alignment_colors["interchrom"] or self.chromosome_palette[0],
                    edgecolor="none",
                    label=(PAIR_CATEGORY_LABELS["interchrom"] if self.alignment_colors["interchrom"] else
                           "Inter-chromosomal (mate hue)"),
                ),
            ]
        haplotype_handles = []
        if self.haplotype_view in ("color", "split"):
            haplotype_handles = [
                Patch(facecolor=haplotype_color("1", self.haplotype_colors, self.chromosome_palette), edgecolor="none", label="HP 1"),
                Patch(facecolor=haplotype_color("2", self.haplotype_colors, self.chromosome_palette), edgecolor="none", label="HP 2"),
                Patch(facecolor=haplotype_color("3", self.haplotype_colors, self.chromosome_palette), edgecolor="none", label="Other HP"),
                Patch(facecolor=haplotype_color(None, self.haplotype_colors, self.chromosome_palette), edgecolor="none", label="Untagged"),
            ]
        middle_handles = haplotype_handles or pair_handles
        middle_title = "Haplotype" if haplotype_handles else "Pair evidence"
        base_handles = []
        for base in "ACGT":
            base_handles.append(Patch(
                facecolor=self.base_colors[base], edgecolor="none", label=base
            ))

        legend_ax = fig.add_axes([
            plot_left, self.legend_bottom_in / fig_height, plot_right - plot_left,
            self.legend_height_in / fig_height,
        ])
        legend_ax.set_xlim(0, 1)
        legend_ax.set_ylim(0, 1)
        legend_ax.set_axis_off()

        if self.fig_width >= 9:
            if middle_handles:
                groups = [
                    ("Alignment events", alignment_handles, 2, 0.00, 0.25, 0.00, 1.00),
                    (middle_title, middle_handles, 2 if haplotype_handles else 3, 0.25, 0.80, 0.00, 1.00),
                    ("Base identity", base_handles, 2, 0.80, 1.00, 0.00, 1.00),
                ]
            else:
                groups = [
                    ("Alignment events", alignment_handles, 2, 0.00, 0.58, 0.00, 1.00),
                    ("Base identity", base_handles, 2, 0.58, 1.00, 0.00, 1.00),
                ]
        elif middle_handles:
            compact_columns = 1 if self.fig_width < 6 else 2
            groups = [
                ("Alignment events", alignment_handles, compact_columns, 0.00, 1.00, 0.68, 1.00),
                (middle_title, middle_handles, compact_columns, 0.00, 1.00, 0.25, 0.68),
                ("Base identity", base_handles, 4, 0.00, 1.00, 0.00, 0.25),
            ]
        else:
            compact_columns = 1 if self.fig_width < 6 else 2
            groups = [
                ("Alignment events", alignment_handles, compact_columns, 0.00, 1.00, 0.50, 1.00),
                ("Base identity", base_handles, 4, 0.00, 1.00, 0.00, 0.50),
            ]

        for index, group in enumerate(groups):
            title, handles, columns, x0, x1, y0, y1 = group
            if index % 2:
                legend_ax.add_patch(Rectangle(
                    (x0, y0), x1 - x0, y1 - y0,
                    transform=legend_ax.transAxes,
                    facecolor=self.visual_colors["legend_background"],
                    edgecolor="none", zorder=0,
                ))
            legend = legend_ax.legend(
                handles=handles, title=title, loc="center", ncol=columns,
                fontsize=7, title_fontsize=7.5, frameon=False,
                bbox_to_anchor=((x0 + x1) / 2, (y0 + y1) / 2),
                bbox_transform=legend_ax.transAxes,
                borderaxespad=0, borderpad=0.25,
                columnspacing=1.2, handletextpad=0.45,
                labelspacing=0.35,
            )
            legend.get_title().set_fontweight("bold")
            legend.get_title().set_color(self.visual_colors["primary_text"])
            legend_ax.add_artist(legend)

        for left_group, right_group in zip(groups, groups[1:]):
            if self.fig_width >= 9:
                divider = left_group[4]
                legend_ax.plot(
                    [divider, divider], [0, 1], transform=legend_ax.transAxes,
                    color=self.visual_colors["legend_edge"], linewidth=0.65, zorder=3,
                )
            else:
                divider = left_group[5]
                legend_ax.plot(
                    [0, 1], [divider, divider], transform=legend_ax.transAxes,
                    color=self.visual_colors["legend_edge"], linewidth=0.65, zorder=3,
                )

        legend_ax.add_patch(Rectangle(
            (0, 0), 1, 1, transform=legend_ax.transAxes,
            facecolor="none", edgecolor=self.visual_colors["legend_edge"],
            linewidth=0.75, clip_on=False, zorder=4,
        ))
        return list(legend_ax.artists)

    def separate_legend_from_plots(self, fig, plot_axes, legends: list) -> None:
        """Enforce a physical gap between plot content and the rendered legend."""
        if not legends:
            return
        fig.canvas.draw()
        canvas_renderer = fig.canvas.get_renderer()
        legend_ax = legends[0].axes
        legend_top = max(
            legend_ax.get_window_extent(canvas_renderer).y1,
            max(legend.get_window_extent(canvas_renderer).y1 for legend in legends),
        )
        plot_bottom = min(
            ax.get_tightbbox(canvas_renderer).y0 for ax in plot_axes if ax is not legend_ax
        )
        required_bottom = legend_top + self.legend_plot_gap_in * fig.dpi
        if plot_bottom >= required_bottom:
            return
        shortfall_fraction = (required_bottom - plot_bottom) / fig.bbox.height
        new_bottom = fig.subplotpars.bottom + shortfall_fraction
        maximum_bottom = fig.subplotpars.top - 0.05
        fig.subplots_adjust(bottom=min(new_bottom, maximum_bottom))
        fig.canvas.draw()

    def render_multi(
        self,
        panels: List[dict],
        chrom: str,
        window_start: int,
        window_end: int,
        reference: Optional[ReferenceWindow],
        out_path: str,
        suptitle: str = "",
        genomic_tracks: Optional[List[LoadedAnnotationTrack]] = None,
        contig_length: Optional[int] = None,
        cytobands: Optional[List[Cytoband]] = None,
    ) -> None:
        """Stack several BAMs' snapshots in one figure, sharing one genomic
        x-axis - the comparison view for "does aligner A produce
        more gapped alignments than aligner B here". Each panel is a dict with
        keys: label, rows, all_reads_for_coverage, dropped_reads and
        downsampled_reads (optional).
        """
        span = window_end - window_start
        show_ref_track = bool(
            reference and reference.available and
            self.max_reference_span > 0 and span <= self.max_reference_span
        )
        render_base_detail = span <= self.max_mismatch_render_span

        tracks = []
        ratios = []
        if self.show_ideogram and contig_length:
            tracks.append("ideogram")
            ratios.append(self.styles["ideogram_height_in"])
        if show_ref_track:
            tracks.append("reference")
            ratios.append(self.styles["reference_height_in"])
        genomic_tracks = genomic_tracks or []
        for index, annotation in enumerate(genomic_tracks):
            tracks.append(f"annotation_{index}")
            if annotation.kind in CNV_TRACK_FORMATS:
                ratios.append(self.styles["cnv_track_height_in"])
            elif annotation.kind in BAF_TRACK_FORMATS:
                ratios.append(self.styles["baf_track_height_in"])
            else:
                ratios.append(
                    max(len(annotation.rows), 1) * self.styles["annotation_row_height_in"]
                )

        panel_track_names = []
        for i, panel in enumerate(panels):
            n_rows = max(len(panel["rows"]), 1)
            header_name = f"panel_header_{i}"
            cov_name, aln_name = f"coverage_{i}", f"alignments_{i}"
            panel_track_names.append((header_name, cov_name, aln_name))
            tracks.append(header_name)
            ratios.append(self.styles["panel_header_height_in"])
            if self.show_coverage:
                tracks.append(cov_name)
                ratios.append(1.2)
            tracks.append(aln_name)
            ratios.append(max(n_rows * self.row_height_in, self.row_height_in))

        top_margin_in = 0.72
        bottom_margin_in = self.legend_margin_in
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

        tick_positions = nice_tick_positions(window_start, window_end)

        fig.text(
            0.01, 1 - 0.06 / fig_height,
            f"{chrom}:{window_start + 1:,}-{window_end:,} ({span:,} bp)",
            fontsize=10.5, color=self.visual_colors["primary_text"], fontweight="bold", va="top", ha="left",
        )
        if suptitle:
            fig.text(
                0.01, 1 - 0.34 / fig_height,
                ellipsize(suptitle, max(30, int(self.fig_width * 15))),
                fontsize=8.5, color=self.visual_colors["secondary_text"], va="top", ha="left",
            )

        if "ideogram" in ax_by_track:
            self.draw_ideogram(
                ax_by_track["ideogram"], chrom, window_start, window_end, contig_length,
                cytobands,
            )

        if show_ref_track:
            self.draw_reference_track(
                ax_by_track["reference"], reference, window_start, window_end,
                available_width_in=self.fig_width,
            )

        for index, annotation in enumerate(genomic_tracks):
            self.draw_annotation_track(
                ax_by_track[f"annotation_{index}"], annotation, window_start, window_end
            )

        for track in tracks:
            if track != "ideogram" and not track.startswith("panel_header_"):
                for tick in tick_positions:
                    ax_by_track[track].axvline(
                        tick, color=self.visual_colors["gridline"],
                        lw=self.styles["grid_line_width"], zorder=0,
                    )

        for i, panel in enumerate(panels):
            rows = panel["rows"]
            layout = panel.get("layout", "pack")
            n_rows = max(len(rows), 1)
            header_name, cov_name, aln_name = panel_track_names[i]

            panel_label = panel.get("label", f"bam{i+1}")
            if panel.get("downsampled_reads"):
                panel_label += f"; {panel['downsampled_reads']} downsampled"
            if panel.get("dropped_reads"):
                panel_label += f"; {panel['dropped_reads']} omitted by --max_rows"
            header_ax = ax_by_track[header_name]
            header_ax.set_ylim(0, 1)
            header_ax.text(
                0.0, 0.45,
                ellipsize(panel_label, max(20, int(self.fig_width * 13))),
                transform=header_ax.transAxes,
                fontsize=9, color=self.visual_colors["primary_text"], fontweight="bold",
                ha="left", va="center", clip_on=True,
            )

            if self.show_coverage:
                cov_ax = ax_by_track[cov_name]
                cov_reads = panel.get("all_reads_for_coverage")
                if not cov_reads:
                    cov_reads = []
                    for row in rows:
                        cov_reads.extend(row)
                self.draw_coverage_track(cov_ax, cov_reads, window_start, window_end)

            aln_ax = ax_by_track[aln_name]
            aln_ax.set_ylim(n_rows, 0)
            self.draw_haplotype_lanes(aln_ax, rows)
            for row_idx, row in enumerate(rows):
                y0 = row_idx + self.row_margin
                h = 1 - 2 * self.row_margin
                self.draw_alignment_row(
                    aln_ax, row, y0, h, render_base_detail, layout
                )
            if not rows:
                aln_ax.text(0.5, 0.5, "No alignments in this region", transform=aln_ax.transAxes,
                            ha="center", va="center", fontsize=9,
                            color=self.visual_colors["secondary_text"])
        for track, ax in ax_by_track.items():
            if track != "ideogram" and not track.startswith("panel_header_"):
                self.draw_center_guide(ax, window_start, window_end)
        bottom_aln_ax = ax_by_track[panel_track_names[-1][2]]
        bottom_aln_ax.set_xticks(tick_positions)
        bottom_aln_ax.tick_params(
            bottom=True, labelbottom=True, labelsize=9,
            length=3, colors=self.visual_colors["primary_text"],
        )
        plot_left = left_margin_fraction(self.fig_width, genomic_tracks)
        if self.haplotype_view == "split":
            plot_left = max(plot_left, min(1.15 / self.fig_width, 0.25))
        plot_right = 0.92
        fig.subplots_adjust(left=plot_left, right=plot_right,
                            top=1 - top_margin_in / fig_height,
                            bottom=bottom_margin_in / fig_height)
        legends = self.draw_legends(fig, fig_height, plot_left, plot_right)
        self.separate_legend_from_plots(fig, axes, legends)
        fig.savefig(out_path)
        plt.close(fig)

    def render_loci(self, panels: List[dict], out_path: str, suptitle: str = "") -> None:
        """Render two independently scaled genomic loci as adjacent panels.

        Unlike :meth:`render_multi`, which stacks BAMs over one shared locus,
        each panel here supplies its own chromosome, bounds, and reference.
        This is the IGV-like mate view used to inspect both sides of a
        discordant or split-read event.
        """
        if len(panels) != 2:
            raise ValueError("Mate view requires exactly two locus panels.")

        max_rows = max(max(len(panel["rows"]), 1) for panel in panels)
        show_ref_track = any(
            bool(
                panel.get("reference") and panel["reference"].available and
                self.max_reference_span > 0 and
                panel["end"] - panel["start"] <= self.max_reference_span
            )
            for panel in panels
        )
        tracks = ["panel_header"]
        ratios = [0.44]
        show_ideogram = self.show_ideogram and any(panel.get("contig_length") for panel in panels)
        if show_ideogram:
            tracks.append("ideogram")
            ratios.append(self.styles["ideogram_height_in"])
        if show_ref_track:
            tracks.append("reference")
            ratios.append(self.styles["reference_height_in"])
        annotation_count = max(len(panel.get("genomic_tracks", [])) for panel in panels)
        annotation_row_counts = []
        for index in range(annotation_count):
            shared_rows = max(
                max(len(panel["genomic_tracks"][index].rows), 1)
                for panel in panels if index < len(panel.get("genomic_tracks", []))
            )
            annotation_row_counts.append(shared_rows)
            tracks.append(f"annotation_{index}")
            is_cnv_track = any(
                panel["genomic_tracks"][index].kind in CNV_TRACK_FORMATS
                for panel in panels if index < len(panel.get("genomic_tracks", []))
            )
            is_baf_track = any(
                panel["genomic_tracks"][index].kind in BAF_TRACK_FORMATS
                for panel in panels if index < len(panel.get("genomic_tracks", []))
            )
            if is_cnv_track:
                ratios.append(self.styles["cnv_track_height_in"])
            elif is_baf_track:
                ratios.append(self.styles["baf_track_height_in"])
            else:
                ratios.append(shared_rows * self.styles["annotation_row_height_in"])
        if self.show_coverage:
            tracks.append("coverage")
            ratios.append(1.4)
        tracks.append("alignments")
        ratios.append(max(max_rows * self.row_height_in, self.row_height_in))

        top_margin_in = 0.5
        bottom_margin_in = self.legend_margin_in
        fig_height = sum(ratios) + top_margin_in + bottom_margin_in
        fig, axes = plt.subplots(
            nrows=len(tracks), ncols=2, squeeze=False,
            figsize=(self.fig_width, fig_height), dpi=self.dpi,
            gridspec_kw={
                "height_ratios": ratios, "hspace": 0.15,
                "wspace": 0.20 if self.haplotype_view == "split" else 0.12,
            },
        )

        if suptitle:
            fig.text(0.01, 0.995, suptitle, fontsize=10.5,
                     color=self.visual_colors["primary_text"],
                     fontweight="bold", va="top", ha="left")

        for panel_idx, panel in enumerate(panels):
            chrom = panel["chrom"]
            start = panel["start"]
            end = panel["end"]
            span = end - start
            rows = panel["rows"]
            layout = panel.get("layout", "pack")
            reference = panel.get("reference")
            self.active_sort_base_position = panel.get(
                "sort_base_position", self.sort_base_position
            )
            self.active_sort_reference_base = panel.get(
                "sort_reference_base", self.sort_reference_base
            )
            render_base_detail = span <= self.max_mismatch_render_span
            axes_by_track = {track: axes[i][panel_idx] for i, track in enumerate(tracks)}
            ticks = nice_tick_positions(start, end, target=4)

            for ax in axes_by_track.values():
                ax.set_xlim(start, end)
                for spine in ("top", "right", "left"):
                    ax.spines[spine].set_visible(False)
                ax.spines["bottom"].set_visible(False)
                ax.tick_params(
                    left=False, labelleft=False, bottom=False, top=False,
                    labelbottom=False, labeltop=False,
                )

            label = panel.get("label", "Primary" if panel_idx == 0 else "Mate")
            dropped = panel.get("dropped_reads", 0)
            dropped_label = f"; {dropped} read(s) omitted" if dropped else ""
            downsampled = panel.get("downsampled_reads", 0)
            downsampled_label = f"; {downsampled} downsampled" if downsampled else ""
            header_ax = axes_by_track["panel_header"]
            header_ax.set_ylim(0, 1)
            header_ax.text(
                0.0, 0.72,
                ellipsize(label, max(16, int(self.fig_width * 6.5))),
                transform=header_ax.transAxes,
                ha="left", va="center", fontsize=8.2,
                color=self.visual_colors["primary_text"],
                fontweight="bold", clip_on=True,
            )
            header_ax.text(
                0.0, 0.16,
                f"{chrom}:{start + 1:,}-{end:,} ({span:,} bp)"
                f"{downsampled_label}{dropped_label}",
                transform=header_ax.transAxes, ha="left", va="center",
                fontsize=7, color=self.visual_colors["secondary_text"], clip_on=True,
            )

            for track in tracks:
                if track not in ("panel_header", "ideogram"):
                    for tick in ticks:
                        axes_by_track[track].axvline(
                            tick, color=self.visual_colors["gridline"],
                            lw=self.styles["grid_line_width"], zorder=0,
                        )

            if show_ref_track:
                ref_ax = axes_by_track["reference"]
                ref_ax.set_ylim(0, 1)
                if (
                    reference and reference.available and
                    self.max_reference_span > 0 and span <= self.max_reference_span
                ):
                    self.draw_reference_track(
                        ref_ax, reference, start, end,
                        available_width_in=self.fig_width / 2,
                    )

            if show_ideogram and panel.get("contig_length"):
                self.draw_ideogram(
                    axes_by_track["ideogram"], chrom, start, end, panel["contig_length"],
                    panel.get("cytobands"),
                )

            panel_annotations = panel.get("genomic_tracks", [])
            for index, annotation in enumerate(panel_annotations):
                self.draw_annotation_track(
                    axes_by_track[f"annotation_{index}"], annotation, start, end,
                    shared_row_count=annotation_row_counts[index],
                )

            if self.show_coverage:
                cov_ax = axes_by_track["coverage"]
                cov_reads = panel.get("all_reads_for_coverage")
                if cov_reads is None:
                    cov_reads = []
                    for row in rows:
                        cov_reads.extend(row)
                self.draw_coverage_track(cov_ax, cov_reads, start, end)

            aln_ax = axes_by_track["alignments"]
            aln_ax.set_ylim(max_rows, 0)
            self.draw_haplotype_lanes(aln_ax, rows)
            aln_ax.set_xticks(ticks)
            aln_ax.tick_params(
                bottom=True, labelbottom=True, labelsize=8,
                length=3, colors=self.visual_colors["primary_text"],
            )
            for row_idx, row in enumerate(rows):
                y0 = row_idx + self.row_margin
                h = 1 - 2 * self.row_margin
                self.draw_alignment_row(
                    aln_ax, row, y0, h, render_base_detail, layout
                )
            if not rows:
                aln_ax.text(0.5, 0.5, "No alignments in this region",
                            transform=aln_ax.transAxes, ha="center", va="center",
                            fontsize=9, color=self.visual_colors["secondary_text"])

            for track, ax in axes_by_track.items():
                if track not in ("panel_header", "ideogram"):
                    self.draw_center_guide(ax, start, end)

        all_genomic_tracks = []
        for panel in panels:
            all_genomic_tracks.extend(panel.get("genomic_tracks", []))
        plot_left = left_margin_fraction(self.fig_width, all_genomic_tracks)
        if self.haplotype_view == "split":
            plot_left = max(plot_left, min(1.15 / self.fig_width, 0.25))
        plot_right = 0.95
        fig.subplots_adjust(left=plot_left, right=plot_right,
                            top=1 - top_margin_in / fig_height,
                            bottom=bottom_margin_in / fig_height)
        legends = self.draw_legends(fig, fig_height, plot_left, plot_right)
        plot_axes = []
        for row in axes:
            plot_axes.extend(row)
        self.separate_legend_from_plots(fig, plot_axes, legends)
        fig.savefig(out_path)
        plt.close(fig)

    def draw_haplotype_lanes(self, ax, rows: List[List[AlignedRead]]) -> None:
        """Shade and label contiguous HP lanes created by split layout mode."""
        if self.haplotype_view != "split" or not rows:
            return
        lanes = []
        lane_start = 0
        lane_haplotype = getattr(rows[0][0], "haplotype", None) if rows[0] else None
        for row_index, row in enumerate(rows[1:], start=1):
            row_haplotype = getattr(row[0], "haplotype", None) if row else None
            if row_haplotype != lane_haplotype:
                lanes.append((lane_start, row_index, lane_haplotype))
                lane_start = row_index
                lane_haplotype = row_haplotype
        lanes.append((lane_start, len(rows), lane_haplotype))

        for index, (start, end, haplotype) in enumerate(lanes):
            color = haplotype_color(
                haplotype, self.haplotype_colors, self.chromosome_palette
            )
            if index % 2 == 0:
                ax.axhspan(
                    start, end, facecolor=color,
                    alpha=self.styles["haplotype_lane_alpha"], zorder=0.2,
                )
            if start:
                ax.axhline(
                    start, color=self.visual_colors["legend_edge"],
                    linewidth=0.7, zorder=1,
                )
            phase_sets = sorted({
                str(read.phase_set)
                for row in rows[start:end] for read in row
                if getattr(read, "phase_set", None) is not None
            })
            label = f"HP {haplotype}" if haplotype is not None else "Untagged"
            if len(phase_sets) == 1:
                label += f" · PS {phase_sets[0]}"
            elif len(phase_sets) > 1:
                label += f" · {len(phase_sets)} PS"
            ax.text(
                -0.012, (start + end) / 2, label,
                transform=ax.get_yaxis_transform(), ha="right", va="center",
                fontsize=6.5, color=color, fontweight="bold", clip_on=False,
            )

    def draw_alignment_row(
        self, ax, row: List[AlignedRead], y0: float, h: float,
        render_base_detail: bool, layout: str,
    ) -> None:
        """Draw one row, including visible mate links and one gap annotation."""
        if self.view_as_pairs:
            pair_members = {}
            for read in row:
                if (
                    read.is_paired and not read.is_secondary
                    and not read.is_supplementary and not read.mate_is_unmapped
                    and read.mate_chrom == read.reference_name
                ):
                    pair_members.setdefault(
                        (read.query_name, read.reference_name), []
                    ).append(read)
            for members in pair_members.values():
                ordered_members = sorted(members, key=lambda read: read.ref_start)
                for index in range(0, len(ordered_members) - 1, 2):
                    left, right = ordered_members[index:index + 2]
                    if left.ref_end < right.ref_start:
                        color, alpha = self.read_style(left)
                        ax.plot(
                            [left.ref_end, right.ref_start],
                            [y0 + h / 2, y0 + h / 2],
                            color=color, alpha=max(alpha, 0.55),
                            linewidth=self.styles["pair_link_width"],
                            zorder=1, solid_capstyle="butt",
                        )

        for read in row:
            self.draw_read(ax, read, y0, h, render_base_detail)

        if layout == "expand" and self.annotate_gap:
            labels = list(dict.fromkeys(
                label for read in row for label in [read.gap_label()] if label
            ))
            if labels:
                ax.text(
                    1.005, y0 + h / 2, " / ".join(labels),
                    transform=ax.get_yaxis_transform(), fontsize=6.5,
                    va="center", ha="left",
                    color=self.visual_colors["secondary_text"], clip_on=False,
                )

    def draw_read(self, ax, read: AlignedRead, y0: float, h: float, render_base_detail: bool) -> None:
        base_fill, alpha = self.read_style(read)

        if not read.blocks:
            ax.add_patch(Rectangle((read.ref_start, y0), max(read.ref_end - read.ref_start, 1), h,
                                    facecolor=base_fill, alpha=alpha,
                                    edgecolor=self.visual_colors["contrast_edge"],
                                    linewidth=self.styles["alignment_edge_width"]))
            return

        n_blocks = len(read.blocks)
        for i, blk in enumerate(read.blocks):
            if blk.op in ("M", "=", "X"):
                ax.add_patch(Rectangle((blk.ref_pos, y0), blk.length, h, facecolor=base_fill,
                                        alpha=alpha, edgecolor=self.visual_colors["contrast_edge"],
                                        linewidth=self.styles["alignment_edge_width"], zorder=2))
            elif blk.op in ("D", "N"):
                color = (
                    self.visual_colors["deletion"] if blk.op == "D"
                    else self.visual_colors["reference_skip"]
                )
                style = "-" if blk.op == "D" else "--"
                ax.plot([blk.ref_pos, blk.ref_pos + blk.length], [y0 + h / 2, y0 + h / 2],
                        color=color, linestyle=style, linewidth=1.3, zorder=3, solid_capstyle="butt")
                if blk.length >= 3:
                    ax.text(blk.ref_pos + blk.length / 2, y0 + h / 2, f"{blk.length}",
                            fontsize=5.5, color=color, ha="center", va="bottom", zorder=4)
            elif blk.op == "I":
                width = 0.4
                insertion_color = self.alignment_colors["small_insert"]
                ax.add_patch(Rectangle((blk.ref_pos - width / 2, y0), width, h, facecolor=insertion_color,
                                        edgecolor="none", zorder=5))
                if blk.length >= 3:
                    ax.text(blk.ref_pos, y0, f"+{blk.length}", fontsize=5.5, color=insertion_color,
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
                                                facecolor=self.base_colors.get(cbase, self.base_colors["N"]),
                                                alpha=alpha, edgecolor="none", zorder=2))
                else:
                    ax.add_patch(Rectangle(
                        (x0, y0), blk.length, h,
                        facecolor=self.visual_colors["softclip"], alpha=alpha,
                        edgecolor=self.visual_colors["contrast_edge"],
                        linewidth=self.styles["alignment_edge_width"],
                        zorder=2,
                    ))
            # 'H' hard clips consume no query bases and are not drawn.

        if render_base_detail:
            for rpos, qbase in read.mismatches:
                ax.add_patch(Rectangle((rpos, y0), 1, h, facecolor=self.base_colors.get(qbase, self.base_colors["N"]),
                                        edgecolor="none", zorder=7))

        sort_position = self.active_sort_base_position
        if sort_position is not None:
            observed = read.base_at(sort_position)
            reference = self.active_sort_reference_base
            if observed in ("A", "C", "G", "T") and observed != reference:
                ax.add_patch(Rectangle(
                    (sort_position, y0), 1, h,
                    facecolor=self.base_colors[observed], edgecolor="none", zorder=8,
                ))
