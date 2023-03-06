# tumor-growth
This project contains the source code relevant for the study titled _"Growth patterns of meningiomas followed by repeated MRI scans"_.

## Implementation setup

The initial statistical analysis was performed in Python 3.7.9 on macOS (12.6 Monterey) using the following libraries:
* pandas==1.3.5
* scipy==1.7.3

The growth analysis was performed using Stata/MP 17 using the [menl](https://www.stata.com/manuals/memenl.pdf) library.

## Project structure

The source code in this project expects some structure on the data, and was tailored for this application and not meant to generalize to new datasets and applications.

    └── tumor-growth/
        ├── src/
        │   ├── stata/
        |   |   └── curve_fitting.do
        │   └──  python/
        |       ├── main.py
        |       └── utils.py
        └── data/
            ├── cohort_personal_info.csv
            ├── cohort_volumes_quality-filtered.csv
            ├── T2_and_peritumorial_oedema.csv
            ├── scanners_info.csv
            └── volumes.csv

Note that the CSV files under `data/` are not provided as this dataset is not made public.

## Perform experiments

1. Setup Python virtual environment and activate it:
```
virtualenv -ppython3 venv --clear
source venv/bin/activate
```

2. Install Python dependencies:
```
pip install -r assets/requirements.txt
```

3. Given that the data lies in the `data/` directory, generate summary statistics by:
```
python src/python/main.py
```

4. Finally, perform growth curve modelling in Stata using the DO-file that lies [here](src/stata/curve_fitting.do).
