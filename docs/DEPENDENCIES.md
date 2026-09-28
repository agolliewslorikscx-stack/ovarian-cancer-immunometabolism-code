# Software dependencies

| Workflow | Packages |
| --- | --- |
| MR | `ieugwasr`, `TwoSampleMR`; optional `httr` for timeout configuration |
| Module association | `ggplot2`, `svglite`, `ragg` |
| Core genes | `readr`, `dplyr`, `tidyr`, `stringr`, `ggplot2`, `gridExtra`, `svglite`, `ragg`, `scales`; base `grid` |
| Enrichment and composition | `ggplot2`, `svglite`, `ragg`; base `grid` |

Obtain packages using their official distribution instructions. This release does not install packages automatically or change the user's R library. `check_dependencies.R` reports installed versions and missing packages. Plot rendering requires usable fonts and Cairo support. The existing scripts request Arial; font substitution can alter layout.

R 4.6.1 parsed all included R files and ran the offline smoke checks during packaging. No full research-data execution or live OpenGWAS query was performed. No claim is made that this R version was used for every original analysis. Package versions are not pinned; versions and API defaults can affect exact reproduction.
