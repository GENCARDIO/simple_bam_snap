"""Validated YAML configuration for behaviour, colours, and plot styling."""
from __future__ import annotations

from copy import deepcopy
from pathlib import Path
from typing import Any, Dict, Optional

import yaml
from matplotlib.colors import is_color_like


DEFAULT_INSERTION_COLOR = "#4a3aa7"
DEFAULT_ALIGNMENT_COLORS: Dict[str, Optional[str]] = {
    "normal": "#b0b0b0",
    "large_insert": "#e34948",
    "small_insert": DEFAULT_INSERTION_COLOR,
    "ff": "#009696",
    "rr": "#1432c8",
    "everted": "#008300",
    "interchrom": None,
}
DEFAULT_BASE_COLORS = {
    # IGV desktop's default alignment/mismatch nucleotide colours.
    "A": "#00ff00", "C": "#0000ff", "G": "#d17105",
    "T": "#ff0000", "N": "#898781",
}
DEFAULT_TRACK_COLORS = {
    "bed": "#000000",
    "gene": "#17217a",
    "vcf": "#7a1f5c",
    "cnv": "#555555",
    "baf": "#7a1f5c",
    "peak": "#7b3294",
    "signal": "#2c7fb8",
}
DEFAULT_VISUAL_COLORS = {
    "deletion": "#0b0b0b",
    "reference_skip": "#898781",
    "softclip": "#1baf7a",
    "coverage": "#a3a3a3",
    "gridline": "#e1e0d9",
    "axis": "#898781",
    "primary_text": "#0b0b0b",
    "secondary_text": "#52514e",
    "legend_edge": "#cfcec8",
    "legend_background": "#f5f5f2",
    "contrast_edge": "#ffffff",
    "label_background": "#ffffff",
    "cytoband_edge": "#9a9a9a",
    "ideogram": "#c4c4c4",
    "ideogram_window": "#d62728",
    "centromere": "#b84b4b",
    "center_guide": "#52514e",
    "cnv_gain": "#d73027",
    "cnv_loss": "#2878b5",
    "sashimi_combined": "#24557f",
    "sashimi_plus": "#24557f",
    "sashimi_minus": "#8b3f73",
}
DEFAULT_HAPLOTYPE_COLORS = {
    "1": "#3b6fb6", "2": "#e6862d", "untagged": "#9b9b96",
}
DEFAULT_CYTOBAND_COLORS = {
    "gneg": "#ffffff", "gpos25": "#c8c8c8", "gpos50": "#969696",
    "gpos75": "#646464", "gpos100": "#1f1f1f", "gvar": "#dddddd",
    "stalk": "#b8d8d8", "acen": "#d95f5f",
}
DEFAULT_CHROMOSOME_PALETTE = [
    "#2a78d6", "#eb6834", "#1baf7a", "#eda100",
    "#e87ba4", "#008300", "#4a3aa7", "#e34948",
]
# IGV's fixed chromosome colours. In "color alignments by pair orientation and
# insert size" mode, inter-chromosomal reads use the colour of the mate's
# chromosome. Unknown contigs still receive a deterministic fallback hue.
DEFAULT_CHROMOSOME_COLORS = {
    "chr1": "#5050ff", "chr2": "#ce3d32", "chr3": "#749b58",
    "chr4": "#f0e685", "chr5": "#466983", "chr6": "#ba6338",
    "chr7": "#5db1dd", "chr8": "#802268", "chr9": "#6bd76b",
    "chr10": "#d595a7", "chr11": "#924822", "chr12": "#837b8d",
    "chr13": "#c75127", "chr14": "#d58f5c", "chr15": "#7a65a5",
    "chr16": "#e4af69", "chr17": "#3b1b53", "chr18": "#cddeb7",
    "chr19": "#612a79", "chr20": "#ae1f63", "chr21": "#e7c76f",
    "chr22": "#5a655e", "chrX": "#cc9900", "chrY": "#99cc00",
    "chrUn": "#404040",
    "chr23": "#cc9900", "chr24": "#99cc00", "chr25": "#33cc00",
    "chr26": "#00cc33", "chr27": "#00cc99", "chr28": "#0099cc",
    "chr29": "#0a47ff", "chr30": "#4775ff", "chr31": "#ffc20a",
    "chr32": "#ffd147", "chr33": "#990033", "chr34": "#991a00",
    "chr35": "#996600", "chr36": "#809900", "chr37": "#339900",
    "chr38": "#00991a", "chr39": "#009966", "chr40": "#008099",
    "chr41": "#003399", "chr42": "#1a0099", "chr43": "#660099",
    "chr44": "#990080", "chr45": "#d60047", "chr46": "#ff1463",
    "chr47": "#00d68f", "chr48": "#14ffb1",
    "chrI": "#8b9bbb", "chrII": "#ce3d32", "chrIII": "#749b58",
    "chrIV": "#f0e685", "chr2a": "#ff5747", "chr2b": "#ff7c65",
}
DEFAULT_STYLES = {
    "row_height_in": 0.06,
    "squish_row_height_in": 0.015,
    "row_margin": 0.08,
    "squish_row_margin": 0.02,
    "annotation_row_height_in": 0.30,
    "cnv_track_height_in": 1.15,
    "baf_track_height_in": 1.05,
    "peak_track_height_in": 1.05,
    "coverage_track_height_in": 1.40,
    "ideogram_height_in": 0.34,
    "panel_header_height_in": 0.30,
    "reference_height_in": 0.38,
    "alignment_alpha": 0.90,
    "secondary_alignment_alpha": 0.50,
    "mapq_alpha_floor": 0.15,
    "coverage_alpha": 0.85,
    "coverage_bins_per_pixel": 1.00,
    "reference_base_alpha": 0.20,
    "cnv_fill_alpha": 0.28,
    "baf_alpha": 0.90,
    "peak_fill_alpha": 0.55,
    "signal_fill_alpha": 0.82,
    "density_fill_alpha": 0.45,
    "haplotype_lane_alpha": 0.055,
    "alignment_edge_width": 0.00,
    "squish_alignment_edge_width": 0.00,
    "annotation_edge_width": 0.30,
    "pair_link_width": 0.85,
    "squish_pair_link_width": 0.40,
    "deletion_line_width": 0.65,
    "squish_deletion_line_width": 0.30,
    "split_read_line_width": 0.55,
    "squish_split_read_line_width": 0.30,
    "grid_line_width": 0.60,
    "gene_line_width": 0.55,
    "primary_gene_line_width": 0.85,
    "gene_arrow_size": 2.70,
    "gene_arrow_spacing_px": 24.0,
    "peak_summit_width": 0.85,
    "signal_line_width": 0.75,
    "signal_y_max": 0.0,
    "center_guide_alpha": 0.65,
    "center_guide_width": 0.80,
    "center_guide_line_style": "--",
    "legend_font_size": 6.80,
    "legend_title_size": 7.30,
    "legend_compartment_gap": 0.010,
    "sashimi_track_height_in": 1.25,
    "sashimi_arc_alpha": 0.85,
    "sashimi_min_line_width": 0.80,
    "sashimi_max_line_width": 4.00,
    "sashimi_arc_height": 0.78,
    "sashimi_label_size": 6.00,
}

DEFAULT_CONFIG: Dict[str, Any] = {
    "alignment_colors": DEFAULT_ALIGNMENT_COLORS,
    "base_colors": DEFAULT_BASE_COLORS,
    "track_colors": DEFAULT_TRACK_COLORS,
    "visual_colors": DEFAULT_VISUAL_COLORS,
    "haplotype_colors": DEFAULT_HAPLOTYPE_COLORS,
    "cytoband_colors": DEFAULT_CYTOBAND_COLORS,
    "chromosome_colors": DEFAULT_CHROMOSOME_COLORS,
    "chromosome_palette": DEFAULT_CHROMOSOME_PALETTE,
    "styles": DEFAULT_STYLES,
    "preferences": {},
}

CONFIG_SECTIONS = set(DEFAULT_CONFIG)
COLOR_SECTIONS = (
    "alignment_colors", "base_colors", "track_colors", "visual_colors",
    "haplotype_colors", "cytoband_colors", "chromosome_colors",
)
ALPHA_STYLE_KEYS = {
    "alignment_alpha", "secondary_alignment_alpha", "mapq_alpha_floor",
    "coverage_alpha", "reference_base_alpha", "cnv_fill_alpha", "baf_alpha",
    "peak_fill_alpha", "signal_fill_alpha", "density_fill_alpha",
    "haplotype_lane_alpha",
    "center_guide_alpha",
    "sashimi_arc_alpha",
}
STRING_STYLE_CHOICES = {
    "center_guide_line_style": ("-", "--", ":", "-."),
}
NONNEGATIVE_STYLE_KEYS = {
    "alignment_edge_width", "squish_alignment_edge_width",
    "annotation_edge_width", "signal_y_max",
}


def read_yaml(path: str) -> Dict[str, Any]:
    config_path = Path(path)
    if not config_path.is_file():
        raise ValueError(f"Config file not found: {path}")
    try:
        with config_path.open("r", encoding="utf-8") as handle:
            data = yaml.safe_load(handle)
    except yaml.YAMLError as exc:
        raise ValueError(f"Invalid YAML in config file '{path}': {exc}") from exc
    if data is None:
        return {}
    if not isinstance(data, dict):
        raise ValueError("Config root must be a YAML mapping.")
    return data


def validate_color(section: str, key: str, value: Any) -> None:
    if section == "alignment_colors" and key == "interchrom" and value is None:
        return
    if not isinstance(value, str) or not is_color_like(value):
        raise ValueError(f"Invalid color for {section}.{key}: {value!r}")


def load_config(path: Optional[str] = None) -> Dict[str, Any]:
    """Load a partial configuration and merge it over every built-in default."""
    config = deepcopy(DEFAULT_CONFIG)
    if not path:
        return config
    data = read_yaml(path)
    unknown_sections = sorted(set(data) - CONFIG_SECTIONS)
    if unknown_sections:
        raise ValueError(f"Unknown config section(s): {', '.join(unknown_sections)}.")

    for section in COLOR_SECTIONS:
        configured = data.get(section, {})
        if not isinstance(configured, dict):
            raise ValueError(f"'{section}' must be a YAML mapping.")
        unknown = (
            [] if section in ("haplotype_colors", "chromosome_colors")
            else sorted(set(configured) - set(config[section]))
        )
        if unknown:
            noun = "alignment color" if section == "alignment_colors" else f"{section} key"
            raise ValueError(f"Unknown {noun}(s): {', '.join(unknown)}.")
        for key, value in configured.items():
            validate_color(section, key, value)
            config[section][key] = value

    palette = data.get("chromosome_palette", config["chromosome_palette"])
    if not isinstance(palette, list) or not palette:
        raise ValueError("'chromosome_palette' must be a non-empty YAML list.")
    for index, value in enumerate(palette):
        validate_color("chromosome_palette", str(index), value)
    config["chromosome_palette"] = list(palette)

    styles = data.get("styles", {})
    if not isinstance(styles, dict):
        raise ValueError("'styles' must be a YAML mapping.")
    unknown_styles = sorted(set(styles) - set(DEFAULT_STYLES))
    if unknown_styles:
        raise ValueError(f"Unknown style key(s): {', '.join(unknown_styles)}.")
    for key, value in styles.items():
        if key in STRING_STYLE_CHOICES:
            if value not in STRING_STYLE_CHOICES[key]:
                choices = ", ".join(STRING_STYLE_CHOICES[key])
                raise ValueError(f"Style {key} must be one of: {choices}.")
            config["styles"][key] = value
            continue
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise ValueError(f"Style {key} must be numeric, got {value!r}.")
        if key in ALPHA_STYLE_KEYS and not 0 <= value <= 1:
            raise ValueError(f"Style {key} must be between 0 and 1.")
        if key in NONNEGATIVE_STYLE_KEYS and value < 0:
            raise ValueError(f"Style {key} must be zero or greater.")
        if key not in ALPHA_STYLE_KEYS and key not in NONNEGATIVE_STYLE_KEYS and value <= 0:
            raise ValueError(f"Style {key} must be greater than zero.")
        config["styles"][key] = float(value)

    preferences = data.get("preferences", {})
    if not isinstance(preferences, dict):
        raise ValueError("'preferences' must be a YAML mapping.")
    config["preferences"] = dict(preferences)
    return config


def load_alignment_colors(path: str) -> Dict[str, Optional[str]]:
    """Backward-compatible loader for callers that only need alignment colours."""
    return load_config(path)["alignment_colors"]
