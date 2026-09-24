"""TLR reporter assay data loading and analysis."""

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt


def load_tlr_data(
    data_dir: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, float]]:
    """Load TLR2 (Pam3) and TLR4 (LPS) data from supplementary tables."""
    tlr4_raw = pd.read_csv(data_dir / "Supplementary_Table_5.csv")
    tlr4_lps = tlr4_raw[tlr4_raw["OD630nm_LPS_Replicate1"].notna()]
    tlr4_df = pd.DataFrame(
        {
            "Concentration_EU_mL": tlr4_lps["Concentration_(EU_mL)"],
            "Average": tlr4_lps[
                ["OD630nm_LPS_Replicate1", "OD630nm_LPS_Replicate2"]
            ].mean(axis=1),
        }
    )
    tlr4_fla = tlr4_raw[tlr4_raw["OD630nm_Fla-PA_Replicate1"].notna()].iloc[0]

    tlr2_raw = pd.read_csv(data_dir / "Supplementary_Table_6.csv")
    tlr2_pam = tlr2_raw[tlr2_raw["OD630nm_Pam3_Replicate1"].notna()]
    tlr2_df = pd.DataFrame(
        {
            "Concentration_ng_mL": tlr2_pam["Concentration_(ng_mL)"],
            "Average": tlr2_pam[
                ["OD630nm_Pam3_Replicate1", "OD630nm_Pam3_Replicate2"]
            ].mean(axis=1),
        }
    )
    tlr2_fla = tlr2_raw[tlr2_raw["OD630nm_Fla-PA_Replicate1"].notna()].iloc[0]

    replicates = ["OD630nm_Fla-PA_Replicate1", "OD630nm_Fla-PA_Replicate2"]
    flapa_data = {
        "tlr4": tlr4_fla[replicates].mean(),
        "tlr2": tlr2_fla[replicates].mean(),
    }

    return tlr2_df, tlr4_df, flapa_data


def plot_tlr_panel(
    ax_main,
    ax_bar,
    df: pd.DataFrame,
    conc_col: str,
    fla_pa_val: float | None = None,
    xlabel: str = "Concentration",
    title: str = "TLR",
    label: str = "Ligand",
    color: str = "#1f77b4",
    xlim: tuple[float, float] = (0.01, 100),
):
    """Plot a single TLR dose-response panel with optional Fla-PA bar."""
    df_plot = df[df[conc_col] > 0].sort_values(conc_col)

    ax_main.plot(
        df_plot[conc_col], df_plot["Average"], "o-", linewidth=2, markersize=8,
        color=color, label=label,
    )

    ax_main.set_xscale("log")
    ax_main.set_xlabel(xlabel, fontsize=12)
    ax_main.set_ylabel("OD (630 nm)", fontsize=12)
    ax_main.set_title(title, fontsize=14)
    ax_main.grid(True, alpha=0.3, linestyle="-", linewidth=0.5)
    ax_main.legend(fontsize=10)
    ax_main.set_xlim(xlim)

    if fla_pa_val is not None:
        ax_bar.bar(
            ["Fla-PA"],
            [fla_pa_val],
            color="lightblue",
            alpha=0.4,
            width=0.5,
            edgecolor="black",
            linewidth=1.2,
        )
        ax_bar.set_ylim(ax_main.get_ylim())
        ax_bar.set_ylabel("")
        ax_bar.tick_params(left=False, labelleft=False)
        ax_bar.grid(True, alpha=0.3, linestyle="-", linewidth=0.5, axis="y")
        ax_main.axhline(y=fla_pa_val, color="black", linestyle=":", linewidth=1.5)
    else:
        ax_bar.axis("off")


def plot_tlr_hek_blue(
    tlr2_df: pd.DataFrame,
    tlr4_df: pd.DataFrame,
    fla_pa_data: dict[str, float],
    output_path: Path,
    output_filename: str,
) -> Path:
    """Create TLR2/TLR4 dose-response plots with Fla-PA bar."""
    fig, axes = plt.subplots(
        2, 2, figsize=(10, 10), gridspec_kw={"width_ratios": [4, 1]}
    )
    ax1, ax1_bar = axes[0]
    ax2, ax2_bar = axes[1]

    plot_tlr_panel(
        ax_main=ax1,
        ax_bar=ax1_bar,
        df=tlr4_df,
        conc_col="Concentration_EU_mL",
        fla_pa_val=fla_pa_data["tlr4"],
        xlabel="Concentration (EU/ml)",
        title="HEK-Blue™ Reporter Line TLR4 LPS",
        label="LPS",
    )

    plot_tlr_panel(
        ax_main=ax2,
        ax_bar=ax2_bar,
        df=tlr2_df,
        conc_col="Concentration_ng_mL",
        fla_pa_val=fla_pa_data["tlr2"],
        xlabel="Concentration (ng/mL)",
        title="HEK-Blue™ Reporter Line TLR2 Pam3",
        label="Pam3",
    )

    plt.tight_layout()

    output_path.mkdir(parents=True, exist_ok=True)

    save_path = output_path / output_filename
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    print(f"Plot saved as '{save_path}'")
    plt.close()

    return save_path
