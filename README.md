## Overview

This page provides the codes for generating random variables, estimation, and data fitting used in 
Vila, Quintino and Bourguignon (2026). Modeling double bounded data based on correlated gamma random variables. Preprint.

## Author Information

- **Roberto Vila**: [Scopus](https://www.scopus.com/authid/detail.uri?authorId=56924856000&origin=AuthorProfile), [Orcid](https://orcid.org/0000-0003-1073-0114),
  [Google Scholar](https://scholar.google.com/citations?hl=en&user=-u38Si8AAAAJ&view_op=list_works&sortby=pubdate)
- **Felipe Quintino** [Scopus](https://www.scopus.com/authid/detail.uri?authorId=57604163400), [Orcid](https://orcid.org/0000-0003-0286-0541), [Google Scholar](https://scholar.google.com/citations?user=BY-KytEAAAAJ&hl=en&oi=ao)
- **Marcelo Bourguignon**: [Scopus](https://www.scopus.com/authid/detail.uri?authorId=55372187300&origin=AuthorProfile), [Orcid](https://orcid.org/0000-0002-1182-5193), [Google Scholar](https://scholar.google.com/citations?user=px50wRYAAAAJ&hl=en&oi=ao)


## Repository Structure

📄 [01_plots](./01-PDF_plots) # Basic codes for PDF plots

📄 [02_estimation](./02-Estimation) # Parameter estimation methods for EB model

📄 [03_modeling](./03-Modelling) # Data modeling scripts (requires 01 and 02)

📄 [04_data](./04-Dataset) # rfam08 and risfam08 datasets

## Requirements

- R (version ≥ 4.5 recommended)
- Required R packages: gsl, MASS, latex2exp, AdequacyModel, goftest

## Citation

If you use these codes in your work, please cite the original paper: 
Vila, R., Quintino, F. and Bourguignon, M., 2026. Modeling double bounded data based on correlated gamma random variables. arXiv preprint arXiv:2603.02566.
