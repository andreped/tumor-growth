# [tumor-growth](https://github.com/andreped/tumor-growth#tumor-growth)

[![license](https://img.shields.io/github/license/DAVFoundation/captain-n3m0.svg?style=flat-square)](https://github.com/DAVFoundation/captain-n3m0/blob/master/LICENSE)

This project contains the source code relevant for the study titled _"Growth dynamics of untreated meningiomas"_ accepted for publication in [Neuro-Oncology Advances](https://academic.oup.com/noa).

## [Setup](https://github.com/andreped/tumor-growth#setup)
The initial statistical analysis was performed in Python 3.7.9 on macOS (12.6 Monterey) using the following libraries:
* [pandas==1.3.5](https://pypi.org/project/pandas/1.3.5/)
* [scipy==1.7.3](https://pypi.org/project/scipy/1.7.3/)

The growth analysis was performed using [Stata/MP 17](https://www.stata.com/statamp/) using the [menl](https://www.stata.com/manuals/memenl.pdf) library.

## [Project structure](https://github.com/andreped/tumor-growth#project-structure)
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

## [Analysis](https://github.com/andreped/tumor-growth#analysis)

1. Setup Python virtual environment and activate it:
```
virtualenv -ppython3 venv --clear
source venv/bin/activate
```

2. Install Python dependencies:
```
pip install -r requirements.txt
```

3. Given that the data lies in the `data/` directory, generate summary statistics by:
```
python src/python/main.py --remove-missing --export-csv
```

4. Finally, perform growth curve modelling in Stata using the DO-file that lies [here](src/stata/curve_fitting.do).

Note that the `main.py` script support various arguments. Run `python src/python/main.py --help` to which arguments are available.

To activate the virtual environment on Windows, instead of `source venv/bin/activate` run `./venv/Scripts/activate`.

## [License](https://github.com/andreped/tumor-growth#license)

The code in this repository is released under [MIT license](https://github.com/andreped/tumor-growth/blob/main/LICENSE).

## [Citation](https://github.com/andreped/tumor-growth#citation)

A bibtex will be added here when the paper is published.
