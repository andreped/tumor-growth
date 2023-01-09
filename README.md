# tumor-growth
Statistical analysis to assess meningioma tumour growth.
Analysis was performed using Python 3.7.9 on macOS. See requirements.txt for info regarding dependencies.

## Implementation setup
The initial statistical analysis was performed in Python 3.7.9 mainly using pandas, numpy, and rpy2.
Two-point regression analysis was conducted in R-4.2.2 through the rpy2 tool in Python.
The time-series regression curve analysis was performed in Stata/MP 17.
The analysis was conducted on a MacBook Pro 2016, with Intel Core i7 CPU, 16 GB RAM, and macOS Monterey (v12.6) operating system.

## Info regarding variables:
* OP_ID: Patients with IDs starting with "O" means that they got surgery (one tumor cluster will have disappeared between timestamps), while those starting with "P" have solely been followed at the outpatient clinic over time.
* Timestamp: Sequential timestamp integer.
* Date: Calendar date (YYYY-MM-DD) for the specified timestamp.
* Clusters total: Total amount of tumor elements identified for the patient, using standard connected components analysis.
* Cluster index: For each line, specifying which of the potential multiple clusters is considered (matching the annotation label value), to track volume evolution over time properly.
* Volume: Volume of the specified tumor cluster, reported in ml, and obtained in patient space.
* Birth_Year: Year of birth (YYYY).
* Gender: Man/woman.
* Dim1/Dim2/Dim3: Dimensions in pixels along each MRI volume axis, order does not necessarily mean X/Y/Z because of potential shenanigans due to volume orientation.
* Spacing1/Spacing2/Spacing3: Same as above but for spacing values.

## Info regarding the inter-variability study:
A total of 30 volumes have been (pseudo-randomly) picked at different timestamps, and Per has been asked to segment them all manually again. In some cases, Per was unable to find the tumor again, hence a 0 ml volume.
* Raw volume: Ground truth volume in ml, obtained after automatic segmentation and manual refinement (from Johanna/Kathrine).
* Per volume: Same as above but from Per as annotator.
* Dice: Dice score between the two annotated volumes.

## Perform experiments
1. Setup Python virtual environment and activate it:
```
virtualenv -ppython3 venv --clear
source venv/bin/activate
```

2. Install Python dependencies:
```
pip install -r requirements.txt
```

3. Install R-4.2.2 and from there install the R dependencies by:
```
install.packages(c("nlme", "car", "outliers", "minpack.lm"))
```

4. Given that the data lies in the `data/` directory, generate summary statistics by:
```
python src/python/main.py
```

5. Finally, perform growth curve modelling in Stata using the DO-file that lies [here](src/stata/curve_fitting.do).
