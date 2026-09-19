# Two-Part-SEM-Template-R

A reusable, generalized R template for two-part (hurdle) semi-continuous
moderated-mediation structural equation models on zero-inflated
continuous/count outcomes (e.g., NSSI frequency, problem-behavior
severity, screen-time overuse, dissociation-episode counts).

## Contents
- `template_twopart_sem.R` — the script. Edit only the CONFIGURATION
  block at the top (your X, mediator(s), zero-inflated Y, and
  moderator(s)); the rest generalizes automatically.
- `toy_data.csv` — 10 rows of illustrative data matching the default
  configuration, so you can test the script end-to-end before using
  your own data.
- `Method_Guide.pdf` — a 2-page explanation of the modeling logic and
  how to adapt the template.

## Quick start
```r
install.packages(c("haven", "psych", "lavaan", "ggplot2", "pscl"))
source("template_twopart_sem.R")   # runs against toy_data.csv by default
```
Then edit the CONFIGURATION block and `data_path` to point at your own
study before re-running.

## Citation
If you use this template, please cite it via its Zenodo DOI (add DOI
here after archiving) and/or acknowledge it in your Methods section.
See `Method_Guide.pdf` for a note on attribution and authorship.

## License
Choose a license before publishing (e.g., MIT for the code). State it
here and add a LICENSE file.
