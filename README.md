# CCM Analysis Package
[![PyPI version](https://img.shields.io/pypi/v/ccm-analysis.svg)](https://pypi.org/project/ccm-analysis/)

This package performs **Convergent Cross Mapping (CCM) analysis** to infer causal relationships in time series data.  
It is optimized for `.py` scripts and standalone Python usage.

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

