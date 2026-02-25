from dataclasses import dataclass
from typing import Tuple

@dataclass
class PlotConfig:
    """Configuration for smudgeplot visualization."""
    palette: str = 'viridis'
    invert_colours: bool = False
    show_histograms: bool = True
    figsize: Tuple[int, int] = (20, 20)
    dpi: int = 100
    fontsize: int = 32 # Font size for axis labels.
    label_fontsize: int = 28
    legend_fontsize: int = 20
    hist_bins_minor: int = 60
    adjust: bool = True # Adjust text alignment for structures at 0.5 minor coverage (default: False).
    xmax: int = 0.49 # Maximum x-coordinate for labels (default: 0.49).

@dataclass
class AnalysisConfig:
    """Configuration for analysis parameters."""
    max_ploidy: int = 16
    min_ploidy: int = 1
    max_structure_range: int = 17
    min_smudge_size_default: float = 0.02
    min_report_size_default: float = 0.03
    rel_coverage_tolerance: float = 0.5
    local_min_multiplier: float = 1.25
    min_coverage_floor: int = 10
    quantile_998: float = 0.998
    quantile_95: int = 95
    max_manhattan_distance: float = 0.5
    task_all_noise_filter: int = 1000
    error_limit: float = 0.7
    smudge_size_cutoff: float = 0 # 0.01  # this is % of all k-mer pairs smudge needs to have to be considered a valid smudge