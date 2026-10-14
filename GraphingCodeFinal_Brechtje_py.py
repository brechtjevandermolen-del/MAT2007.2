#Cytokine graphs using Python.

#All statistical analysis is performed in R.
#Python only reads the exported R results and makes the graphs.

#Required R output files:
#    significant_cytokine_summary.csv
#    fold_change_results.csv
from pathlib import Path #handles file locations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch


#MAKE sure you have first run the R script and checked that the significant_cytokine_summary and fold_change_results have been generated!
#1Update the filepath to match the location of significant_cytokine_summary and fold_change_results. 
# ----------------------------------------------------------------------
FOLDER = Path(r"C:\Users\brech\OneDrive - Maastricht University\MSP\PRA\PRA2026")
SUMMARY_PATH = FOLDER / "significant_cytokine_summary.csv"
FOLD_CHANGE_PATH = FOLDER / "fold_change_results.csv"

# the graphs (PNG) are saved in the same folder
OUTPUT = FOLDER

GROUPS = ["No Dex", "Dex"]
GROUP_LABELS = ["Untreated", "Treated"]

# Colors for each significant cytokine: untreated, treated, untreated error, treated error
BAR_STYLE = {
    "IL-6":      ("#f7b257", "#e36f60", "#000000", "#000000"),
    "IL-10":     ("#d085df", "#cb4f73", "#000000", "#000000"),
    "IFN-gamma": ("#89eed0", "#a0a4d1", "#000000", "#000000"),
}

#establish the colors for the different cytokine function types
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
    #specify the significant cytokines
    cytokines = ["IL-6", "IL-10", "IFN-gamma"]
    #we want them to be plotted next to one another
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.5))

    #zip pairs each panel with a cytokine
    for ax, cyt in zip(axes, cytokines):    
        #assign the colours of a cytokine to each variable (bar_u = untreated bar color, bar_t = treated color, and then error is both black.)
        bar_u, bar_t, err_u, err_t = BAR_STYLE[cyt]

        #Keep only the rows of each cytokine (summary.loc), rownames are the Treatments, and sort (reindex) by NO dex, then dex
        rows = (
            summary.loc[summary["Cytokine"] == cyt]
            .set_index("Treatment")
            .reindex(GROUPS)
        )
        #pull the means and standard deviation from the R files
        means = rows["Mean"].to_numpy()
        sems = rows["SEM"].to_numpy()

        #bars are placed at 0 and 1 with its own color.
        x = np.array([0, 1])
        ax.bar(x, means, width=0.5, color=[bar_u, bar_t])

        for xi, mean, sem, error_color in zip(
            x, means, sems, [err_u, err_t]
        ):
            ax.errorbar(
                xi,
                mean,
                yerr=sem,
                fmt="none",     #No marker, only error bar
                ecolor=error_color,
                capsize=4,      #width error bar feet
                lw=1.5,
            )

        #have untreated and treated on the x axis
        ax.set_xticks(x)
        ax.set_xticklabels(GROUP_LABELS)
        ax.set_xlim(-0.6, 1.6)

        #set y axis label
        ax.set_ylabel("Concentration (pg/mL)")

        #title per significant cytokine
        ax.set_title(f"{cyt} Concentration\nby Dexamethasone treatment")
        ax.spines[["top", "right"]].set_visible(False)

        #establish y max (max mean + max sem * 1.25 so theres a little room)
        ymax = np.max(means + sems)
        ax.set_ylim(0, ymax * 1.25)

        # p-value and FDR were calculated in R and exported with the summary
        p_value = rows["p.value"].iloc[0]
        fdr = rows["FDR"].iloc[0]

        #add the p and fdr values to the graphs
        ax.text(
            0.7,
            ymax * 1.12,
            f"p = {p_value:.2e}\nFDR = {fdr:.2e}",
            ha="center",
            va="bottom",
            fontsize=8,
        )

        #add legend
        ax.legend(
            handles=[
                Patch(color=bar_u, label="w/o Dexamethasone"),
                Patch(color=bar_t, label="Dexamethasone"),
            ],
            fontsize=7,
            loc="upper left",
            frameon=False,
        )

    #title for entire graph
    fig.suptitle("Significant Cytokines following Dexamethasone treatment")
    fig.tight_layout()
    fig.savefig(OUTPUT / "significant_cytokines_py.png", dpi=150)
    plt.show()


# Graph 2: log2 fold change with SEM already calculated in R
# ----------------------------------------------------------------------
def graph2(fold_change):
    #pull columns out of R
    cytokines = fold_change["Cytokine"].tolist()
    log2_fc = fold_change["log2_FC"].to_numpy()
    sem_log2 = fold_change["SEM_log2_FC"].to_numpy()
    upper = fold_change["Upper"].to_numpy()
    lower = fold_change["Lower"].to_numpy()

    #establish the colors for each cytokine group
    colors = [
        GROUP_COLORS[group]
        for group in fold_change["Group"]
    ]

    #establish X positions
    x = np.arange(len(cytokines))

    fig, ax = plt.subplots(figsize=(8, 6))

    #plot colored bars
    ax.bar(
        x,
        log2_fc,
        width=0.6,
        color=colors,
        zorder=2,
    )

    #Plot error bars
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

    #Plot the black line at 0
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

    #iterrows() goes through the table row by row. xi is the row number, which equals the x position of that bar
    #grab the significant cytokines and add asterisks either below or above bar depending on fold change. 
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

    #set x ticks and their labels
    ax.set_xticks(x)
    ax.set_xticklabels(cytokines, rotation=90)
    #set x axis label
    ax.set_xlabel("Immune-system cytokines", labelpad=10)

    #set y label
    ax.set_ylabel("Log2 Fold Change (Dex / No Dex)")
    #title
    ax.set_title(
        "Log2 fold change in Cytokine concentration\n"
        "after Dexamethasone treatment"
    )

    ax.spines[["top", "right"]].set_visible(False)
    ax.margins(x=0.03)

    #add a legend with the functions
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

    #save the figure in the output
    fig.tight_layout()
    fig.savefig(OUTPUT / "fold_change_py.png", dpi=150)
    plt.show()



# Make graphs
# ----------------------------------------------------------------------
graph1(summary)
graph2(fold_change)
