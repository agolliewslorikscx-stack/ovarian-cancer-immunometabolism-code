# Ovarian cancer immunometabolism analysis code

Selected existing R analysis scripts for the study of diet-related metabolic contexts and myeloid immune remodeling in ovarian cancer.

## Release scope

This is a **partial code release**, not an end-to-end reproduction of all manuscript results. It contains five analysis scripts from the research project, with project-root paths made portable. Statistical routines, thresholds and historical output names are retained. `SOURCE_MANIFEST.json` records their source filenames, checksums and packaging changes.

| Command | Included analysis | Required inputs |
| --- | --- | --- |
| `mr` | One exposure/outcome pair: instrument retrieval, harmonization and two-sample MR | Public OpenGWAS IDs and the user's own OpenGWAS authentication |
| `module-association` | Gene-wise z-score programme means, correlations with existing module eigengenes, BH correction and overlap-excluded sensitivity | TCGA expression, precomputed WGCNA assignments/eigengenes and fixed programme definitions |
| `core-genes` | MEblue/MEbrown core ranking and bootstrap rank stability | The same inputs plus programme scores from `module-association` |
| `enrichment` | Module over-representation analysis and pathway redundancy reduction | WGCNA assignments and separately obtained GMT gene sets |
| `composition` | Donor-level broad cell-compartment proportions and group comparisons | Annotated public single-cell metadata |

The repository contains **no original experimental data**, participant-level datasets, expression matrices, microscopy images, immunoblots, qPCR/ELISA measurements, simulated manuscript results, manuscript documents, or credentials. It also does not include the complete upstream WGCNA construction, cell annotation, trajectory, survival or NHANES workflows. Input data and frozen intermediate objects must be obtained separately. Historical figure numbers in filenames are provenance labels and do not certify correspondence to final manuscript panels.

## Requirements

R 4.1 or later is required by the native pipe syntax. Syntax and isolated base-R function checks were performed with R 4.6.1. These checks are not a full rerun against research data. Required packages are listed by workflow in [docs/DEPENDENCIES.md](docs/DEPENDENCIES.md); an exact historical package lockfile is not available.

Start from the repository root. Use `Rscript --vanilla` so that project-local startup files are not loaded. Set `OPENGWAS_JWT` in your own environment for MR; never commit its value.

```sh
# List supported commands and inspect requirements
Rscript --vanilla run_analysis.R --help
Rscript --vanilla check_dependencies.R

# Check an input tree without running analysis or downloading data
Rscript --vanilla run_analysis.R module-association /path/to/input-project --check-only

# Run modules in order, after supplying inputs
Rscript --vanilla run_analysis.R module-association /path/to/input-project
Rscript --vanilla run_analysis.R core-genes /path/to/input-project
Rscript --vanilla run_analysis.R enrichment /path/to/input-project
Rscript --vanilla run_analysis.R composition /path/to/input-project

# One public-data MR pair; ./work is a local output directory
Rscript --vanilla run_analysis.R mr ./work met-d-DHA ieu-a-1120

# Offline code checks; no research data or network access
Rscript --vanilla tests/smoke.R
python tools/audit_release.py
```

See [docs/INPUTS.md](docs/INPUTS.md) for the input layout, fields and execution order, and [docs/METHOD_NOTES.md](docs/METHOD_NOTES.md) for the retained implementation details and limitations. Analyses write generated outputs into the supplied project directory; use a separate working copy of inputs to avoid replacing earlier outputs. MR requires network access and is subject to the database's access conditions and quota.

## Data and code availability wording

This repository supports a statement that **selected analysis scripts** are available. It does not support a statement that all analysis code, original experimental data, or all figure source data have been deposited here.

## 中文说明

本仓库整理并公开五部分现有分析代码：单暴露–结局MR、WGCNA模块与髓系程序关联、核心基因排序、模块通路富集、供者水平细胞组成比较。仅对项目路径做可移植化处理，并增加运行入口、依赖说明与代码检查。

这是部分代码公开版本，不是整篇论文的完整复现包。仓库不附原始实验数据、患者或样本数据、表达矩阵、实验图片、模拟结果和访问凭据。读者需另行取得相应公开数据及预处理结果；本次整理未重新运行研究分析，也未认证最终稿全部图表与这些历史脚本逐项一致。
