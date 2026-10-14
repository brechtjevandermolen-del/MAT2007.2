# MAT2007.2
**Dexamethasone and Cytokine Concentrations**  
  Analysis of how dexamethasone treatment affects the concentration of 18 immune-system cytokines and markers. All statistics are done in R and the final graphs are made in Python, which only reads the results exported by R.

**Download** dex_cytokine.tsv via the following link: https://datadryad.org/downloads/file_stream/3137209 
  
**Save the downloaded file in the same folder as the scripts**

**Files**:
The files uploaded in this Github include _MAT2007_Rcode_Brechtje_ (Statistical analysis, R bar graphs, and export of results for Python), _GraphingCodeFinal_Brechtje_py_ (Python graphs made from the exported R results), _MAT2007_Report_Brechtje_ (Report on the context of the project, methodology and results), and lastly _MAT2007_Overleaf_Files_ (Containing the overleaf final report as well as the generated figures and bibliogrpahy). 

**Created by the R script (saved in the working directory):**
- **significant_cytokine_summary.csv**: Mean, SEM, p-value and FDR per group for IL-6, IL-10 and IFN-gamma (read by Python)
- **fold_change_results.csv**:	Log2 fold change, SEM, error bar limits, cytokine group and significance for all 18 cytokines (read by Python)

**Created by the Python script:**
- **significant_cytokines_py.png**:	Bar graphs (mean ± SEM) for the three significant cytokines, with p-value and FDR
- **fold_change_py.png**:	Log2 fold change of all 18 cytokines, coloured by group, with asterisks for significant ones

**Setup**:
For this code to work you are required to have installed R (any recent version; only base R is used, no extra packages), as well as Python 3 with pandas, numpy and matplotlib (3.4 or newer).
You can install the python add-ons via "pip install pandas numpy matplotlib"

**Folder**:
Put dex_cytokine.tsv, MAT2007_Rcode_Brechtje and GraphingCodeFinal_Brechtje_py in one folder. All outputs are written to that same folder.

**Update the paths**
**The scripts contain the author's folder path. Change it to your own folder in both places:**
- In MAT2007_Rcode_Brechtje: setwd("C:/Users/brech/OneDrive - Maastricht University/MSP/PRA/PRA2026")
- In GraphingCodeFinal_Brechtje_py: FOLDER = Path(r"C:\Users\brech\OneDrive - Maastricht University\MSP\PRA\PRA2026")

**How to run**:
Run Project_2_1.R first (the whole script). This performs the analysis and writes the two CSV files that Python needs.
Check that significant_cytokine_summary.csv and fold_change_results.csv now exist in the folder.
Run Graphing_code_py_Brechtje.py. It reads the two CSV files and saves the two PNG figures.
