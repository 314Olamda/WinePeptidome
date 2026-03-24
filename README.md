# 🍷 Wine Peptidome — UniProt Small Proteins & Large Peptides Pipeline

> *Automated retrieval and cross-validation of the yeast and grapevine peptidome from dual UniProt APIs — built for wine lees research, SO₂-alternative discovery, and antifungal peptidomics.*

[![Python](https://img.shields.io/badge/Python-3.8%2B-blue?logo=python&logoColor=white)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![UniProt REST](https://img.shields.io/badge/UniProt-REST%20API-8B0000)](https://rest.uniprot.org)
[![Proteins API](https://img.shields.io/badge/EBI-Proteins%20API-003087)](https://www.ebi.ac.uk/proteins/api)
[![ORCID](https://img.shields.io/badge/ORCID-0000--0002--7720--3733-a6ce39?logo=orcid)](https://orcid.org/0000-0002-7720-3733)

---

## 🔬 What is this?

Wine lees — the yeast-rich by-product of fermentation — are a largely untapped reservoir of bioactive peptides. Released through yeast autolysis during *sur lies* aging, these small proteins and large peptides (500 Da → 100 kDa) exhibit antioxidant, antifungal, and antimicrobial properties that make them compelling natural alternatives to SO₂ in winemaking.

This pipeline performs **systematic, reproducible retrieval** of the full peptidome landscape for the two organisms at the heart of wine biochemistry:

| Organism | TaxID | Role |
|---|---|---|
| *Saccharomyces cerevisiae* | 559292 | Primary winemaking yeast; source of lees peptides via autolysis |
| *Vitis vinifera* | 29760 | Grapevine; source of grape-derived macromolecules and defence peptides |

It combines two complementary UniProt data sources to give you both **canonical annotation** and **experimental mass-spectrometry evidence** in a single run.

---

## 🧬 Pipeline Architecture

```mermaid
graph TD
    A[🚀 Start pipeline] --> B[UniProt REST API<br/>rest.uniprot.org]
    A --> C[Proteins API<br/>ebi.ac.uk/proteins/api]

    B --> D[KB Entry Query<br/>mass:501–100000 + taxon]
    D --> E[Link-header pagination<br/>no memory overflow]
    E --> F[UniParc fallback<br/>for deleted accessions]
    E --> G[FASTA stream download]

    C --> H[Coverage metadata check<br/>before bulk download]
    H --> I[/proteomics/non-ptm<br/>Experimental peptide evidence]
    H --> J[/proteomics/ptm<br/>PTM-modified peptides<br/>entry_pos = pep_start + pos − 1]

    F --> K[(output/)]
    G --> K
    I --> K
    J --> K

    K --> L[uniprot_kb_entries.tsv]
    K --> M[sequences.fasta]
    K --> N[proteomics_nonptm_peptides.tsv]
    K --> O[proteomics_ptm_peptides.tsv]
    K --> P[pipeline_summary.txt]
```

---

## ⚡ Quick start

```bash
# 1. Clone
git clone https://github.com/314Olamda/WinePeptidome.git
cd WinePeptidome

# 2. Install dependency (single package)
pip install requests

# 3. Run
python uniprot_wine_peptidome.py
```

Results land in `./output/` automatically.

---

## 📦 Output files

| File | Content | Key columns |
|---|---|---|
| `uniprot_kb_entries.tsv` | UniProt KB canonical entries | accession, gene, mass, annotation_score, GO terms, subcellular location |
| `sequences.fasta` | FASTA sequences for all retrieved entries | standard FASTA format |
| `proteomics_nonptm_peptides.tsv` | Experimental peptide evidence (non-PTM) | peptide, pep_start, pep_end, psm_count, dataset_id |
| `proteomics_ptm_peptides.tsv` | PTM-modified peptides with canonical site mapping | modification, ptm_pos_in_entry, confidence, final_site_prob |
| `pipeline_summary.txt` | Run statistics and coverage report | totals per organism and globally |

---

## ⚙️ Configuration

All adjustable parameters live in one block at the top of the script:

```python
ORGANISMS = {
    "Saccharomyces_cerevisiae": 559292,
    "Vitis_vinifera":           29760,
}

MASS_MIN =     501      # Da — strictly > 500 Da
MASS_MAX = 100_000      # Da — 100 kDa ceiling
PAGE_SIZE =    500      # entries per paginated request
```

To add a third organism (e.g. *Hanseniaspora uvarum*, TaxID 5483) for non-*Saccharomyces* lees research, simply extend the `ORGANISMS` dict.

---

## 🧪 Scientific context

This pipeline was developed within the **reLees** project (H.F.R.I. / NextGenerationEU) and ongoing work at **ISVV – Université de Bordeaux** on antioxidant and antifungal peptides from white wine lees.

The mass window (501 Da → 100 kDa) was chosen to:
- Exclude free amino acids and small dipeptides (already well characterised)
- Capture the biologically active mid-range peptides released during yeast autolysis
- Include the HMW protein fraction with known antifungal properties against *Phaeoacremonium minimum* and *Phaeomoniella chlamydospora*

The **dual-API strategy** cross-validates canonical annotation (UniProt KB) with experimental mass-spectrometry evidence (Proteins API), producing a confirmed peptidome — only accessions present in both layers have genuine biological support.

---

## 🔗 Related resources

- [UniProt REST API documentation](https://rest.uniprot.org/docs)
- [EBI Proteins API documentation](https://www.ebi.ac.uk/proteins/api/doc)
- [GNPS molecular networking](https://gnps.ucsd.edu/) — for downstream peptidomics annotation
- [reLees project](https://relees.uniwa.gr) — wine lees circular economy research

---

## 📄 Citation

If you use this pipeline in your research, please cite:

```bibtex
@software{gimenez_gil_wine_peptidome_2025,
  author  = {Giménez-Gil, Pol},
  title   = {Wine Peptidome: UniProt Small Proteins and Large Peptides Pipeline},
  year    = {2025},
  url     = {https://github.com/314Olamda/WinePeptidome},
  orcid   = {0000-0002-7720-3733}
}
```

---

## 👤 Author

**Pol Giménez-Gil**, PhD  
Postdoctoral Researcher — ISVV, Université de Bordeaux  
Scopus ID: 57219336109 · ORCID: [0000-0002-7720-3733](https://orcid.org/0000-0002-7720-3733)  
ResearchGate: [Pol_Gimenez2](https://www.researchgate.net/profile/Pol_Gimenez2)

---

## 📜 License

MIT — see [LICENSE](LICENSE)
