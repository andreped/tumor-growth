# tumor-growth
Statistical analysis to assess meningioma tumour growth.
Analysis was performed using Python 3.7.9 and Stata/MP 17 on macOS.

## Implementation setup
The initial statistical analysis was performed in Python 3.7.9 mainly using SciPy, Pandas, and NumPy.

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
