# tumor-growth
Analysis was performed using Python 3.7.9 on macOS. See requirements.txt for info regarding dependencies.

## Info regarding variables:
* OP_ID: Patients with IDs starting with "O" means that they got surgery (one tumor cluster will have disappeared between timestamps), while those starting with "P" have solely been followed at the outpatient clinic over time.
* Timestamp: Sequential timestamp integer.
* Date: calendar date (YYYY-MM-DD) for the specified timestamp.
* Clusters total: total amount of tumor elements identified for the patient, using standard connected components analysis.
* Cluster index: for each line, specifying which of the potential multiple clusters is considered (matching the annotation label value), to track volume evolution over time properly.
* Volume: volume of the specified tumor cluster, reported in ml, and obtained in patient space.
* Birth_Year: year of birth (YYYY).
* Gender: man/woman.
* Dim1/Dim2/Dim3: dimensions in pixels along each MRI volume axis, order does not necessarily mean X/Y/Z because of potential shenanigans due to volume orientation.
* Spacing1/Spacing2/Spacing3: same as above but for spacing values.

## Info regarding the inter-variability study:
A total of 30 volumes have been (pseudo-randomly) picked at different timestamps, and Per has been asked to segment them all manually again. In some cases, Per was unable to find the tumor again, hence a 0 ml volume.
* Raw volume: ground truth volume in ml, obtained after automatic segmentation and manual refinement (from Johanna/Kathrine).
* Per volume: same as above but from Per as annotator.
* Dice: Dice score between the two annotated volumes.

## Implementation setup
To perform the statistical analysis data was preprocessed in Python 3.7.9 using pandas and numpy.
Regression analysis were conducted in R-4.2.1 through the rpy2 tool in Python.
