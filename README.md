# [tumor-growth](https://github.com/andreped/tumor-growth#tumor-growth)

[![license](https://img.shields.io/github/license/DAVFoundation/captain-n3m0.svg?style=flat-square)](https://github.com/DAVFoundation/captain-n3m0/blob/master/LICENSE)

This project contains the source code relevant for the study titled [_"Growth dynamics of untreated meningiomas"_](https://academic.oup.com/noa/advance-article/doi/10.1093/noajnl/vdad157/7484549) published in [Neuro-Oncology Advances](https://academic.oup.com/noa).

<details open>
<summary>

## [Abstract](https://github.com/andreped/tumor-growth#abstract)</summary>

    Background: Knowledge about meningioma growth characteristics is needed for
    developing biologically rational follow-up routines. In this study of
    untreated meningiomas followed with repeated MRIs, we studied growth
    dynamics and explored potential factors associated with tumor growth.
    
    Methods: In a single-center cohort study, we included 235 adult patients
    with a radiologically suspected intracranial meningioma and at least three
    MRI scans during follow-up. Tumors were segmented using an automatic
    algorithm from contrast enhanced T1-series, and if needed manually
    corrected. Potential meningioma growth curves were statistically compared;
    linear, exponential, linear radial, or Gompertzian. Factors associated with
    growth were explored.

    Results: In 235 patients, 1394 MRI scans were carried out in the median
    five-year observational period. Of the models tested, a Gompertzian growth
    curve best described growth dynamics of meningiomas on group level. 59 % of
    the tumors grew, 27 % remained stable, and 14 % shrunk. Only 13 patients(5%)
    underwent surgery during the observational period and were excluded after
    surgery. Tumor size at time of diagnosis, multifocality, and length of
    follow-up were associated with tumor growth, whereas age, sex, presence of
    peritumoral edema or hyperintense T2-signal were not significant factors.

    Conclusion: Untreated meningiomas follow a Gompertzian growth curve,
    indicating that increasing and potentially doubling of subsequent follow-up
    intervals between MRIs seems biologically reasonable, instead of fixed time
    intervals. Tumor size at diagnosis is the strongest predictor of future
    growth, indicating a potential for longer follow up intervals for smaller
    tumors. Although most untreated meningiomas grow, few require surgery.

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

If you found our research article or this repository relevant in your research, consider citing our paper:

```
@article{10.1093/noajnl/vdad157,
    title = {{Growth dynamics of untreated meningiomas}},
    author = {Strand, Per Sveino and Wågø, Kathrine Jørgensen and Pedersen, André and Reinertsen, Ingerid and Nälsund, Olivia and Jakola, Asgeir Store and Bouget, David and Hosainey, Sayied Abdol Mohieb and Sagberg, Lisa Millgård and Vanel, Johanna and Solheim, Ole},
    journal = {Neuro-Oncology Advances},
    pages = {vdad157},
    year = {2023},
    month = {12},
    issn = {2632-2498},
    doi = {10.1093/noajnl/vdad157},
    url = {https://doi.org/10.1093/noajnl/vdad157},
}
```
