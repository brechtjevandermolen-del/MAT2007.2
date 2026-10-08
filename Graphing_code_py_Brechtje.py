#Cytokine graphs using Python.

#All statistical analysis is performed in R.
#Python only reads the exported R results and makes the graphs.

#Required R output files:
#    significant_cytokine_summary.csv
#    fold_change_results.csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch


#MAKE sure you have first run the R script and checked that the significant_cytokine_summary and fold_change_results have been generated!
#1Update the filepath to match the location of significant_cytokine_summary and fold_change_results. 
# ----------------------------------------------------------------------
BASE_DIR = Path(r"C:\Users\brech\OneDrive - Maastricht University\MSP\PRA\PRA2026")
SUMMARY_PATH = BASE_DIR / "significant_cytokine_summary.csv"
FOLD_CHANGE_PATH = BASE_DIR / "fold_change_results.csv"

# the graphs (PNG) are saved in the same folder
OUT_DIR = BASE_DIR

GROUPS = ["No Dex", "Dex"]
GROUP_LABELS = ["Untreated", "Treated"]

# Colours for each significant cytokine: untreated, treated, untreated error, treated error
BAR_STYLE = {
    "IL-6":      ("#f7b257", "#e36f60", "#000000", "#000000"),
    "IL-10":     ("#d085df", "#cb4f73", "#000000", "#000000"),
    "IFN-gamma": ("#89eed0", "#a0a4d1", "#000000", "#000000"),
}

#establish the colours for the different cytokine function types
GROUP_COLORS = {
    "Vascular Injury markers": "#ffdae0",
    "Pro-inflammatory cyto-/chemokines": "#EE799F",
    "Anti-inflammatory cytokines": "#FFDEAD",
    "Tissue damage markers": "#DDA0DD",
}



# Load results already calculated in R
# ----------------------------------------------------------------------
summary = pd.read_csv(SUMMARY_PATH)
fold_change = pd.read_csv(FOLD_CHANGE_PATH)


# Graph 1: mean +/- SEM for significant cytokines
# ----------------------------------------------------------------------
def graph1(summary):
    cytokines = ["IL-6", "IL-10", "IFN-gamma"]
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.5))

    for ax, cyt in zip(axes, cytokines):
        bar_u, bar_t, err_u, err_t = BAR_STYLE[cyt]

        rows = (
            summary.loc[summary["Cytokine"] == cyt]
            .set_index("Treatment")
            .reindex(GROUPS)
        )

        means = rows["Mean"].to_numpy()
        sems = rows["SEM"].to_numpy()

        x = np.array([0, 1])

        ax.bar(x, means, width=0.5, color=[bar_u, bar_t])

        for xi, mean, sem, error_color in zip(
            x, means, sems, [err_u, err_t]
        ):
            ax.errorbar(
                xi,
                mean,
                yerr=sem,
                fmt="none",
                ecolor=error_color,
                capsize=4,
                lw=1.5,
            )

        ax.set_xticks(x)
        ax.set_xticklabels(GROUP_LABELS)
        ax.set_xlim(-0.6, 1.6)
        ax.set_ylabel("Concentration (pg/mL)")
        ax.set_title(f"{cyt} Concentration\nby Dexamethasone treatment")
        ax.spines[["top", "right"]].set_visible(False)

        ymax = np.max(means + sems)
        ax.set_ylim(0, ymax * 1.25)

        # p-value and FDR were calculated in R and exported with the summary
        p_value = rows["p.value"].iloc[0]
        fdr = rows["FDR"].iloc[0]

        ax.text(
            0.7,
            ymax * 1.12,
            f"p = {p_value:.2e}\nFDR = {fdr:.2e}",
            ha="center",
            va="bottom",
            fontsize=8,
        )

        ax.legend(
            handles=[
                Patch(color=bar_u, label="w/o Dexamethasone"),
                Patch(color=bar_t, label="Dexamethasone"),
            ],
            fontsize=7,
            loc="upper left",
            frameon=False,
        )

    fig.suptitle("Significant Cytokines following Dexamethasone treatment")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "significant_cytokines_py.png", dpi=150)
    plt.show()


# Graph 2: log2 fold change with SEM already calculated in R
# ----------------------------------------------------------------------
def graph2(fold_change):
    cytokines = fold_change["Cytokine"].tolist()
    log2_fc = fold_change["log2_FC"].to_numpy()
    sem_log2 = fold_change["SEM_log2_FC"].to_numpy()
    upper = fold_change["Upper"].to_numpy()
    lower = fold_change["Lower"].to_numpy()

    colors = [
        GROUP_COLORS[group]
        for group in fold_change["Group"]
    ]

    x = np.arange(len(cytokines))

    fig, ax = plt.subplots(figsize=(8, 6))

    ax.bar(
        x,
        log2_fc,
        width=0.6,
        color=colors,
        zorder=2,
    )

    ax.errorbar(
        x,
        log2_fc,
        yerr=sem_log2,
        fmt="none",
        ecolor="#8B475D",
        capsize=3,
        lw=1.5,
        zorder=3,
    )

    ax.axhline(0, color="black", lw=2, zorder=4)

    # Dotted grey guide lines
    for xi in x:
        ax.axvline(
            xi,
            color="lightgrey",
            ls=":",
            lw=2,
            zorder=1,
        )

    # Significance was determined in R
    pad = 0.015 * (upper.max() - lower.min())

    for xi, row in fold_change.iterrows():
        if row["Significant"]:
            if row["log2_FC"] >= 0:
                ax.text(
                    xi,
                    row["Upper"] + pad,
                    "*",
                    ha="center",
                    va="center",
                    fontsize=13,
                )
            else:
                ax.text(
                    xi,
                    row["Lower"] - pad,
                    "*",
                    ha="center",
                    va="top",
                    fontsize=13,
                )

    ax.set_xticks(x)
    ax.set_xticklabels(cytokines, rotation=90)
    ax.set_xlabel("Immune-system cytokines", labelpad=10)
    ax.set_ylabel("Log2 Fold Change (Dex / No Dex)")
    ax.set_title(
        "Log2 fold change in Cytokine concentration\n"
        "after Dexamethasone treatment"
    )

    ax.spines[["top", "right"]].set_visible(False)
    ax.margins(x=0.03)

    ax.legend(
        handles=[
            Patch(color=color, label=label)
            for label, color in GROUP_COLORS.items()
        ],
        fontsize=8,
        loc="lower right",
        frameon=True,
        facecolor="white",
    )

    fig.tight_layout()
    fig.savefig(OUT_DIR / "fold_change_py.png", dpi=150)
    plt.show()



# Make graphs
# ----------------------------------------------------------------------
graph1(summary)
graph2(fold_change)
