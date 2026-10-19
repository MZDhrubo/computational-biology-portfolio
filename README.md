# Computational Biology Portfolio

**Md. Morshed Zaman | Bioinformatics, Structural Biology and Omics Workflows**

This repository complements my protein-phylogeny work by documenting other computational areas connected to my research and training: transcriptomics quality control, protein-structure file inspection, and transparent review of docking and molecular-descriptor tables.

I want a reader to be able to open the code, understand what it does, run it on a small example, and see its limits. The included datasets are synthetic and are used only to demonstrate software behavior. They are not results from my research projects.

## What is included

| Area | What this repository demonstrates | Important boundary |
|---|---|---|
| RNA-seq / transcriptomics | Count-matrix validation, sample-level QC, log-CPM transformation, and exploratory PCA | This is exploratory QC, not differential-expression analysis or a completed production RNA-seq study. |
| Structural bioinformatics | A lightweight inventory of protein PDB files: chains, residues, C-alpha coverage, occupancy, and B-factor summaries | This is file inspection, not a full structure-quality validation or molecular-dynamics result. |
| Docking and compound triage | Transparent ranking and descriptor-rule flags for a docking results table | The example does not predict ADMET, prove binding, or establish a drug candidate. |
| Immunoinformatics | A documented workflow checklist for vaccine-candidate research | The original project sequences, construct, simulation outputs, and unpublished findings are not included. |
| Molecular dynamics | A reproducibility checklist reflecting GROMACS training context | No claim is made that the included code performs or validates a full MD simulation. |

## Repository layout

```text
morshed-zaman-computational-biology-portfolio/
├── portfolio_tools/
│   ├── rnaseq_qc.py
│   ├── docking_descriptor_triage.py
│   └── pdb_inventory.py
├── examples/
│   ├── synthetic_gene_counts.csv
│   ├── synthetic_docking_descriptors.csv
│   └── toy_protein.pdb
├── tests/
├── docs/
│   ├── CV_PROJECT_MAP.md
│   ├── IMMUNOINFORMATICS_WORKFLOW.md
│   ├── GROMACS_REPRODUCIBILITY.md
│   └── RESEARCH_DATA_POLICY.md
├── CONTRIBUTIONS.md
├── requirements.txt
└── requirements-dev.txt
```

## Quick start

Python 3.10 or newer is recommended.

```bash
python -m venv .venv
```

Activate the environment on Windows PowerShell:

```powershell
.venv\Scripts\Activate.ps1
```

On macOS or Linux:

```bash
source .venv/bin/activate
```

Install dependencies:

```bash
python -m pip install -r requirements.txt
```

## 1. Run the RNA-seq quality-control demonstration

```bash
python -m portfolio_tools.rnaseq_qc --counts examples/synthetic_gene_counts.csv --outdir results/rnaseq_demo
```

The input is a gene-by-sample count matrix: the first column contains gene identifiers and each remaining column is a sample. The tool writes sample QC metrics, exploratory PCA coordinates, and a run summary. It checks for missing values, duplicate identifiers, negative counts, and non-integer counts before analysis.

**Interpretation note:** log-CPM and PCA are exploratory summaries. They do not replace a design-aware differential-expression workflow such as DESeq2 or edgeR, and the synthetic example should not be used for biological inference.

## 2. Run the docking/descriptors table demonstration

```bash
python -m portfolio_tools.docking_descriptor_triage --input examples/synthetic_docking_descriptors.csv --outdir results/docking_demo
```

The tool sorts records by docking score and marks compounds against a small set of explicitly stated, illustrative descriptor thresholds. It does not calculate docking scores, predict ADMET properties, or determine whether a compound binds a target. The output is a transparent table-processing example only.

## 3. Inspect a protein PDB file

```bash
python -m portfolio_tools.pdb_inventory --input examples/toy_protein.pdb --out results/pdb_inventory.json
```

The report summarizes chains, residue counts, C-alpha atom coverage, occupancy, and B-factors. These checks help describe a structure file; they are not sufficient to establish that a model is biologically correct or structurally validated.

## Tests

Install the development dependency and run the tests:

```bash
python -m pip install -r requirements-dev.txt
python -m pytest -q
```

The tests cover input validation and key output behavior for the demonstration utilities. Passing these tests does not validate external docking, structure-prediction, immune-simulation, or molecular-dynamics tools.

## How this relates to my CV

My CV describes experience and training in functional genomics, vaccine immunoinformatics, protein structure modelling, docking and virtual screening, ADMET analysis, and bulk/single-cell RNA-seq. `docs/CV_PROJECT_MAP.md` explains how this public code portfolio relates to those areas without releasing unpublished project data or claiming that every CV activity is implemented here as software.



## A note to reviewers

I built this portfolio so reviewers can inspect the actual code, reproduce the examples, and see where I distinguish computational output from biological interpretation. The examples are intentionally modest: reliable input checks, traceable processing, and honest limitations matter as much as producing a result table.
