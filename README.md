# Just Two CT Scans: Detecting Nodule Progression

**_`NPDS4Clib`_** is an R package designed to calculate NPDS statistics and assess. This repository contains the code for [**A Novel Statistic Guided by Clinical Experience for Assessing the Progression of Lung Nodules**](url). 
the progression of lung nodules using baseline and follow-up CT images. It offers 
a simple and efficient approach for lung nodule analysis, inspired by clinical practices. 
This package is an Rcpp-based version of the Python version **_`NPDS4Clib`_**, originally developed by Hang Yu.

Check:
- [@Hang Yu](https://github.com/hangyustat)'s [NPDS calculation in Python](https://github.com/hangyustat/NPDS)
For python version.

The **_NPDS_** statistic, introduced in the paper **_“A Novel Statistic Guided 
by Clinical Experience for Assessing the Progression of Lung Nodules”_**, provides 
a simple way to determine if a lung nodule has progressed using hypothesis testing.

In this package, you can easily use your baseline and follow-up CT images to calculate the **_NPDS_** statistic, and assess whether the lung nodule has shown any progression.

## Before starting to calculate your own **_NPDS_**, **make sure**:

1. Your CT images are in **NIfTI (.nii)** or compressed **NIfTI (.nii.gz)** format.
If they are in DICOM format, please convert them to NIfTI first.

2. The two CT scans should have similar slice thickness (usually listed in scan 
settings or metadata), ideally differing by no more than 1–2 mm and the lung images 
are sharp and free of distortions, graininess, or visual artifacts that could affect analysis.

3. The nodule’s center coordinates **(X, Y)**, **Z-axis range**, and **maximum diameter** (mm), 
all obtained from the follow-up CT scan., have already been measured and recorded.

If all the data are ready, replace the paths and variables in **Step 2** below with 
your own, then rerun the code and check your own **_NPDS_**!

## Step 1 : Import the necessary package for detecting Nodule Progression.

```r
# If devtools is not installed, install it first
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# Load devtools package
library(devtools)

# Install the NPDS4Clib package from GitHub
devtools::install_github("JingyuGUi/NPDS4Clib")

# Load the NPDS4Clib package
library(NPDS4Clib)
```

## Step 2 : Replace the nii CT data paths, nodule coordinates, and nodule diameter with your own.

```r
# Initialize the nodule progress detector with specified parameters
nodule_progress_detector <- initialization(
  X = 209, 
  Y = 356, 
  range_Z = '325-347', 
  diameter = 12, 
  baseline_CT_nii_path = "Data/0002358111-20180516.nii.gz", 
  followup_CT_nii_path = "Data/0002358111-20220707.nii.gz"  
)

# Perform image registration using elastix and update the nodule progress detector
nodule_progress_detector <- registration_by_elastix(nodule_progress_detector)

# Segment the lungs from the CT scan and update the nodule progress detector
nodule_progress_detector <- get_segmented_lungs(nodule_progress_detector)
```

Below is a detailed explanation of the parameters. You can also find this information 
by typing _"?initialization"_ in the R console.

- **X**:
The X-coordinate of the center of the nodule in the follow-up CT scan. This coordinate 
should correspond to the slice where the nodule’s maximum diameter is largest.

- **Y**:
The Y-coordinate of the center of the nodule in the follow-up CT scan, similarly 
based on the slice with the largest maximum diameter.

- **range_Z**:
A string specifying the range of Z slices that include the nodule in the follow-up 
CT scan. For example, '325-347' indicates that the nodule spans slices 325 to 347.

- **diameter**:
The largest value of the nodule’s maximum diameter across slices in the follow-up 
CT scan, measured in millimeters.

- **baseline_CT_nii_path**:
The file path to the baseline CT scan. This file should be in .nii format.

- **followup_CT_nii_path**:
The file path to the follow-up CT scan, also in .nii format.

## Step 3 : Run the program to perform calculations and hypothesis testing.

```r
# NPDS_calculate_based_on_C++
nodule_progress_detectorC <- NPDS_calculateC(nodule_progress_detector)
nodule_progress_detectorC$NPDS

# Hypothesis test
test <- hypothesis_test_by_ClinvNod_sample(nodule_progress_detectorC)
```

After you run this function, it provides the following key results to help you 
understand and evaluate nodule progression:

1. **Nodule Progression Detection Score** (**_NPDS_**):
This numeric score quantifies the degree of detected progression in the nodule. 
Higher scores suggest more pronounced progression. It serves as a critical indicator 
for tracking changes over time.

2. **Progression Prediction** (**_Progression_**):
A boolean value (TRUE or FALSE) that directly answers whether the nodule is likely progressing:
  -  _TRUE_ : The analysis suggests significant nodule progression.
  - _FALSE_ : No meaningful progression is detected based on the data.

3. **Statistical Significance** (**_p_value_**):
The p-value evaluates the reliability of the NPDS compared to a standard clinical 
reference distribution.
  - _A low p-value_ (e.g., < 0.05) indicates that the progression is statistically 
  significant and unlikely to be due to chance.
  - _A higher p-value_ suggests weaker evidence for progression, implying that the 
  detected NPDS might not reflect clinically meaningful change.

These outputs are designed to work together: the **_NPDS_** gives a measurable 
progression score, the **_progression_** provides a straightforward interpretation, 
and the **_p-value_** ensures the statistical rigor of the result. 

By combining these elements, the function helps you make data-driven decisions with 
both clarity and confidence.

## Acknowledgments
Portions of this package include code adapted from the `EBImage` package 
(https://github.com/aoles/EBImage), which is licensed under LGPL.

## License
This package is licensed under LGPL-3.0. For more details, please refer to the LICENSE file.


## Getting Help
To learn more about specific functions, use the following command in R:

```R
?function_name
```

## Contributors
We'd like to acknowledge the following contributors to this project:

- [@Hang Yu](https://github.com/hangyustat): Developed the original Python version of NPDS calculation.
- [@Jingyu Gui](https://github.com/JingyuGUi): Developed the R version of NPDS calculation.
- [@Yan Lan](https://github.com/lyannnisme): Documentation and testing.



For example, to learn more about the `registration_by_elastix` function, use:

```R
?registration_by_elastix
```
