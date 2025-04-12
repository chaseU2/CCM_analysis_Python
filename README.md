# CCM Analysis Package
[![PyPI version](https://img.shields.io/pypi/v/ccm-analysis.svg)](https://pypi.org/project/ccm-analysis/)

Package implementing Convergent Cross Mapping for causality inference in dynamical systems as defined by ["Sugihara et al (2012)"](https://www.science.org/doi/10.1126/science.1227079).
> **Note:** This package is optimized for Python scripts (`.py`)scripts and standalone Python usage. If you need to perform CCM analysis in a Python notebook (`.ipynb`), please use [CCM_analysis_Python_Notebook](https://github.com/chaseU2/CCM_analysis_Python_Notebook).


## Features
- Accepts **tab-separated files (`.tsv`)** or **pandas DataFrames** as input.
- Computes **CCM scores** for all variable pairs.
- Performs **convergence analysis** to validate causal relationships.
- Saves **results, plots, and logs** in an output directory.
- Generates **heatmaps** of significant causal relationships.

## Installation
To install the package, run:

```bash
pip install ccm_analysis
```

## Example Usage

To use the `run_ccm_analysis` function in your project, follow this simple example:

### Example 1: CCM Analysis with a DataFrame

```python
# Import the ccm_analysis package
from ccm_analysis import run_ccm_analysis
import pandas as pd

# Load example data (replace this with your own file path or DataFrame)
# Assume you have a file "example_data.tsv" with time series data in tab-separated format
data = pd.read_csv("example_data.tsv", delimiter="\t")

# Run the CCM analysis
results = run_ccm_analysis(
    data_input=data,        # Pass the loaded DataFrame here
    L=110,                  # Time series length
    E=2,                    # Embedding dimension
    tau=1,                  # Time delay
    THRESHOLD=0.8,          # Significance threshold
    save_output=True,       # Save results
    output_dir="output"     # Output directory for results
)

# The results will be stored in the `results` variable.
# You can further use or save the results.

```

## Test Dataset

You can test the `ccm_analysis` package using the sample dataset provided in this repository. You can download the test dataset from the following link:

- [Test Dataset (ccm_test_data.txt)](https://github.com/chaseU2/ccm-analysis-tool/blob/master/ccm_test_data.txt)

### Data Format

The dataset is structured as follows:

- **Columns** represent the variables you want to analyze.
- **Rows** represent the time points for each variable.


---

## Interpreting Convergence Plots

The package generates diagnostic convergence plots that reveal causal relationships between variables. Below are the three characteristic patterns to analyze:


When analyzing convergence plots, you'll be prompted:

```python
Analyzing: A ↔ B
Final ρ: A→B: 0.82, B→A: 0.41

Convergence in? (1=both, 2=A→B only, 3=B→A only, 0=none): 
```
p is in this case the calculated crossmap score 


Here you have to take a look at the shown plot and decide in which directions you can see a convergence of the crossmap skill with increasing library size.
You have to enter a number according to this table:

### Options Table
| Key | Action                  |
|-----|-------------------------|
| 1   | Keep both directions    |
| 2   | Keep only A→B           |
| 3   | Keep only B→A           |
| 0   | Discard both            |

---

### 1. No Significant Causality

![No Causal Relationship](https://raw.githubusercontent.com/chaseU2/CCM_analysis_Python/main/ccm_analysis/Screenshot%204.png)

- Neither directional curve (X→Y in blue, Y→X in red) show a clear convergence to a final crossmap score
- Example use case: Independent systems

---

### 2. Unidirectional Causality
![Unidirectional Causality](https://raw.githubusercontent.com/chaseU2/CCM_analysis_Python/main/ccm_analysis/Screenshot%203.png)

- One direction (X→Y) converges to a final crossmap score
- Reverse direction (Y→X) does ot show a clear convegence
- Interpretation: X drives Y but not vice versa
- Pay attention to whether convergence is present in the X→Y or Y→X plot, and enter 2 or 3 accordingly

---

### 3. Bidirectional Causality
![Bidirectional Causality](https://raw.githubusercontent.com/chaseU2/CCM_analysis_Python/main/ccm_analysis/Screenshot%202.png)

- Both directions show positive convergence
- Typical of feedback systems
- Convergence rates may differ (e.g., X→Y stronger than Y→X)


---

## Dependencies and Acknowledgements

Parts of this project are based on the Convergent Cross Mapping (CCM) implementation from the following repository from Prince Javier :

- [Convergent Cross Mapping GitHub Repository Prince Javier ](https://github.com/PrinceJavier/causal_ccm.git)

I have utilized parts of the CCM algorithm from this repository to help analyze causality in time series data in my own project.
