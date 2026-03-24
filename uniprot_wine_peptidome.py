#!/usr/bin/env python3
"""
╔══════════════════════════════════════════════════════════════════════════════╗
║          UniProt Wine Peptidome Pipeline  –  Small Proteins & Large Peptides                        ║
║          Organisms : Saccharomyces cerevisiae (TaxID 559292)                                        ║
║                      Vitis vinifera          (TaxID 29760)                                          ║
║          Mass range: 501 Da  →  100 kDa                                                             ║
╚══════════════════════════════════════════════════════════════════════════════╝

Two UniProt APIs combined into a single reproducible pipeline
─────────────────────────────────────────────────────────────
1. UniProt REST API  (rest.uniprot.org)
   • Paginated KB entry retrieval with mass:[501 TO 100000] filter
   • Link-header pagination – no memory overflow for large result sets
   • UniParc fallback for inactive / deleted accessions
   • FASTA stream download

2. Proteins API  (www.ebi.ac.uk/proteins/api)
   • Coverage metadata check before downloading (avoids empty calls)
   • /proteomics/non-ptm  – experimental peptide evidence (23 species)
   • /proteomics/ptm      – PTM-modified peptides (yeast phosphorylation)
   • PTM site mapping: entry_pos = pep_start + pos_in_pep − 1

Output files (written to ./output/)
────────────────────────────────────
  uniprot_kb_entries.tsv           UniProt KB canonical entries + annotation
  sequences.fasta                  FASTA for all retrieved KB accessions
  proteomics_nonptm_peptides.tsv   Experimental peptide evidence (non-PTM)
  proteomics_ptm_peptides.tsv      PTM-modified peptides with site positions
  pipeline_summary.txt             Run statistics and coverage report

Requirements
────────────
  pip install requests
  Python >= 3.8

Author  : Pol Giménez-Gil
ORCID   : 0000-0002-7720-3733
Affil.  : ISVV – Université de Bordeaux / AUTH / UNIWA
GitHub  : github.com/314Olamda
"""

import csv
import json
import sys
import time
from pathlib import Path

try:
    import requests
except ImportError:
    sys.exit("ERROR: 'requests' is not installed.  Run: pip install requests")


# ═══════════════════════════════════════════════════════════════════════════════
#  CONFIGURATION  –  edit only this block if needed
# ═══════════════════════════════════════════════════════════════════════════════

UNIPROT_API  = "https://rest.uniprot.org"
PROTEINS_API = "https://www.ebi.ac.uk/proteins/api"

ORGANISMS = {
    "Saccharomyces_cerevisiae": 559292,   # S288C reference strain
    "Vitis_vinifera":           29760,
}

MASS_MIN   =     501      # Da  (strictly > 500 Da)
MASS_MAX   = 100_000      # Da  (100 kDa ceiling)

# Approximate AA length window (avg 111 Da/residue + 18 Da water)
LEN_MIN    = max(1, int((MASS_MIN - 18) / 111))
LEN_MAX    = int((MASS_MAX - 18) / 111)

PAGE_SIZE  = 500          # KB entries per page
SLEEP_S    = 0.35         # polite delay between requests (s)
OUTPUT_DIR = Path(__file__).parent / "output"

# UniProt KB fields to retrieve
KB_FIELDS = ",".join([
    "accession", "id", "gene_names", "protein_name",
    "organism_name", "organism_id",
    "length", "mass",
    "sequence_version", "reviewed",
    "annotation_score",
    "go_p", "go_f", "go_c",
    "cc_subcellular_location",
    "ft_signal", "ft_propep",
    "protein_existence",
    "xref_pdb",
])

# ═══════════════════════════════════════════════════════════════════════════════
#  UTILITY FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════════════

def _get(url: str, params: dict = None, headers: dict = None,
         retries: int = 4) -> requests.Response:
    """GET with exponential back-off and informative error messages."""
    h = {"Accept": "application/json"}
    if headers:
        h.update(headers)
    for attempt in range(retries):
        try:
            r = requests.get(url, params=params, headers=h, timeout=60)
            if r.status_code == 200:
                return r
            if r.status_code == 429:          # rate-limited
                wait = 2 ** attempt
                print(f"  [429 rate-limit] waiting {wait}s …")
                time.sleep(wait)
                continue
            if r.status_code == 400:
                print(f"  [400 bad request] {url}")
                return r
            print(f"  [HTTP {r.status_code}] attempt {attempt+1}/{retries}")
            time.sleep(2 ** attempt)
        except requests.RequestException as exc:
            print(f"  [connection error] {exc}")
            time.sleep(2 ** attempt)
    raise RuntimeError(f"Failed after {retries} attempts: {url}")


def save_tsv(rows: list, fields: list, path: Path):
    """Write list-of-dicts to a TSV file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields,
                                delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    print(f"  ✓ Saved {len(rows):,} rows → {path.name}")


def write_summary(lines: list, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")
    print(f"  ✓ Summary → {path.name}")


# ═══════════════════════════════════════════════════════════════════════════════
#  1. UniProt KB  –  REST API
# ═══════════════════════════════════════════════════════════════════════════════

def fetch_uniprot_kb(org_name: str, tax_id: int) -> tuple[list, list]:
    """
    Retrieve all UniProt KB entries for one organism within the mass window.
    Uses Link-header pagination (recommended for large result sets).

    Returns (list_of_row_dicts, list_of_field_names).
    """
    query = (
        f"(taxonomy_id:{tax_id}) "
        f"AND (mass:[{MASS_MIN} TO {MASS_MAX}])"
    )
    url    = f"{UNIPROT_API}/uniprotkb/search"
    params = {
        "query":    query,
        "fields":   KB_FIELDS,
        "format":   "tsv",
        "size":     PAGE_SIZE,
    }

    rows, header, page_num = [], [], 1
    print(f"\n  Querying UniProt KB: {org_name} …")
    print(f"  Query: {query}")

    while url:
        r = _get(url, params=params if page_num == 1 else None)
        time.sleep(SLEEP_S)

        lines = r.text.strip().splitlines()
        if not lines:
            break

        if page_num == 1 and lines[0].startswith("Entry"):
            header = lines[0].split("\t")
            data_lines = lines[1:]
        else:
            data_lines = lines[1:] if lines[0].startswith("Entry") else lines

        for line in data_lines:
            if line.strip():
                values = line.split("\t")
                row = dict(zip(header, values))
                rows.append(row)

        total_results = int(r.headers.get("X-Total-Results", 0))
        link_header   = r.headers.get("Link", "")
        url           = None
        if 'rel="next"' in link_header:
            import re
            m = re.search(r'<([^>]+)>;\s*rel="next"', link_header)
            if m:
                url = m.group(1)
        params = None   # params embedded in next URL

        print(f"    Page {page_num:3d} | "
              f"{len(rows):,}/{total_results:,} entries retrieved")
        page_num += 1

    print(f"  → {len(rows):,} KB entries for {org_name}")
    return rows, header


def check_uniparc_fallback(accession: str) -> str | None:
    """
    Check if an accession is active in UniProt KB.
    If not (deleted / merged), fall back to UniParc and return the UPI.
    Mirrors the pattern demonstrated in the UniProt programmatic access webinar.
    """
    r = _get(f"{UNIPROT_API}/uniprotkb/{accession}", params={"format": "json"})
    if r.status_code == 200:
        return accession   # still active

    print(f"    [UniParc fallback] {accession} not in KB – querying UniParc …")
    r2 = _get(f"{UNIPROT_API}/uniparc/accession/{accession}",
              params={"format": "json"})
    if r2.status_code == 200:
        data = r2.json()
        upi  = data.get("uniParcId", "NOT_FOUND")
        print(f"    → UniParc UPI: {upi}")
        return upi
    print(f"    → {accession} not found in UniParc either")
    return None


def fetch_fasta(org_name: str, tax_id: int) -> str:
    """Download FASTA for the entire mass-filtered set via the stream endpoint."""
    query = (
        f"(taxonomy_id:{tax_id}) "
        f"AND (mass:[{MASS_MIN} TO {MASS_MAX}])"
    )
    url = f"{UNIPROT_API}/uniprotkb/stream"
    params = {"query": query, "format": "fasta"}
    print(f"  Downloading FASTA for {org_name} …")
    r = _get(url, params=params, headers={"Accept": "text/plain"})
    time.sleep(SLEEP_S)
    seq_count = r.text.count(">")
    print(f"  → {seq_count:,} FASTA sequences for {org_name}")
    return r.text


# ═══════════════════════════════════════════════════════════════════════════════
#  2. Proteins API  –  Experimental Proteomics Evidence
# ═══════════════════════════════════════════════════════════════════════════════

def check_proteomics_coverage(tax_id: int) -> dict:
    """
    Query the Proteins API metadata endpoint to discover which pipelines
    (non-ptm, ptm, hpp) have data for this taxon BEFORE downloading.
    Returns dict of {pipeline: bool}.
    """
    pipelines = {"nonptm": False, "ptm": False}
    for pipe in pipelines:
        endpoint = "non-ptm" if pipe == "nonptm" else "ptm"
        r = _get(
            f"{PROTEINS_API}/proteomics/{endpoint}",
            params={"taxid": tax_id, "size": 1},
            headers={"Accept": "application/json"},
        )
        time.sleep(SLEEP_S)
        if r.status_code == 200:
            data = r.json()
            if isinstance(data, list) and len(data) > 0:
                pipelines[pipe] = True
            elif isinstance(data, dict) and data.get("results"):
                pipelines[pipe] = True
    return pipelines


def fetch_proteomics_nonptm(org_name: str, tax_id: int) -> list:
    """
    Retrieve experimental non-PTM peptide evidence from Proteins API.
    Applies a peptide length filter as a proxy for the mass window
    (the Proteins API does not expose mass directly).
    """
    print(f"  Fetching non-PTM proteomics evidence for {org_name} …")
    url    = f"{PROTEINS_API}/proteomics/non-ptm"
    params = {"taxid": tax_id, "size": 100, "offset": 0}
    records, page_n = [], 0

    while True:
        r = _get(url, params=params, headers={"Accept": "application/json"})
        time.sleep(SLEEP_S)
        if r.status_code != 200:
            print(f"  [non-PTM] HTTP {r.status_code} — stopping")
            break
        batch = r.json()
        if not batch:
            break

        for entry in batch:
            acc     = entry.get("accession", "")
            org     = entry.get("organism", {}).get("names", [{}])[0].get("value", "")
            taxid_e = entry.get("organism", {}).get("taxonomy", tax_id)
            for feat in entry.get("features", []):
                pep   = feat.get("peptide", "")
                start = feat.get("begin", 0)
                end   = feat.get("end", 0)
                ln    = len(pep)
                if not (LEN_MIN <= ln <= LEN_MAX):
                    continue   # mass-window proxy filter
                for src in feat.get("sources", [{}]):
                    records.append({
                        "organism":      org,
                        "tax_id":        taxid_e,
                        "accession":     acc,
                        "peptide":       pep,
                        "length_aa":     ln,
                        "approx_mass_Da": round(ln * 111 + 18),
                        "pep_start":     start,
                        "pep_end":       end,
                        "unique_to_isoform": feat.get("unique", False),
                        "dataset_id":    src.get("id", ""),
                        "psm_count":     src.get("count", ""),
                        "pipeline":      "non-ptm",
                    })

        page_n += 1
        print(f"    Page {page_n:3d} | {len(records):,} peptide records")
        params["offset"] += 100
        if len(batch) < 100:
            break

    print(f"  → {len(records):,} non-PTM records for {org_name}")
    return records


def fetch_proteomics_ptm(org_name: str, tax_id: int) -> list:
    """
    Retrieve PTM-modified peptide evidence from Proteins API.
    Maps each modification to the canonical entry position using:
        entry_pos = pep_start + pos_in_pep − 1
    """
    print(f"  Fetching PTM proteomics evidence for {org_name} …")
    url    = f"{PROTEINS_API}/proteomics/ptm"
    params = {"taxid": tax_id, "size": 100, "offset": 0}
    records, page_n = [], 0

    while True:
        r = _get(url, params=params, headers={"Accept": "application/json"})
        time.sleep(SLEEP_S)
        if r.status_code != 200:
            print(f"  [PTM] HTTP {r.status_code} — stopping")
            break
        batch = r.json()
        if not batch:
            break

        for entry in batch:
            acc     = entry.get("accession", "")
            org     = entry.get("organism", {}).get("names", [{}])[0].get("value", "")
            taxid_e = entry.get("organism", {}).get("taxonomy", tax_id)
            for feat in entry.get("features", []):
                pep       = feat.get("peptide", "")
                pep_start = feat.get("begin", 0)
                pep_end   = feat.get("end", 0)
                ln        = len(pep)
                if not (LEN_MIN <= ln <= LEN_MAX):
                    continue
                for mod in feat.get("ptms", []):
                    pos_in_pep  = mod.get("position", 0)
                    entry_pos   = pep_start + pos_in_pep - 1   # canonical mapping
                    records.append({
                        "organism":       org,
                        "tax_id":         taxid_e,
                        "accession":      acc,
                        "peptide":        pep,
                        "length_aa":      ln,
                        "approx_mass_Da": round(ln * 111 + 18),
                        "pep_start":      pep_start,
                        "pep_end":        pep_end,
                        "ptm_pos_in_entry": entry_pos,
                        "modification":   mod.get("name", ""),
                        "confidence":     mod.get("confidence", ""),
                        "final_site_prob":mod.get("probability", ""),
                        "psm_count":      mod.get("count", ""),
                        "dataset_id":     mod.get("sources", [{}])[0].get("id", "")
                                          if mod.get("sources") else "",
                        "pipeline":       "ptm",
                    })

        page_n += 1
        print(f"    Page {page_n:3d} | {len(records):,} PTM records")
        params["offset"] += 100
        if len(batch) < 100:
            break

    print(f"  → {len(records):,} PTM records for {org_name}")
    return records


# ═══════════════════════════════════════════════════════════════════════════════
#  MAIN
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    print("\n" + "═" * 65)
    print("  UniProt Wine Peptidome Pipeline")
    print(f"  Mass window  : {MASS_MIN:,} – {MASS_MAX:,} Da  (~{LEN_MIN}–{LEN_MAX} aa)")
    print(f"  Organisms    : {', '.join(ORGANISMS)}")
    print("═" * 65)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    all_kb_rows   = []
    all_kb_header = []
    all_fasta     = ""
    all_nonptm    = []
    all_ptm       = []

    summary_lines = [
        "UniProt Wine Peptidome Pipeline – Run Summary",
        f"Mass window  : {MASS_MIN} – {MASS_MAX} Da",
        f"AA length    : ~{LEN_MIN} – {LEN_MAX} residues",
        "",
    ]

    for org_name, tax_id in ORGANISMS.items():
        print(f"\n{'═'*65}")
        print(f"  ORGANISM: {org_name}  (TaxID {tax_id})")
        print(f"{'═'*65}")

        # ── A. UniProt KB entries ──────────────────────────────────────────
        kb_rows, kb_header = fetch_uniprot_kb(org_name, tax_id)
        all_kb_rows.extend(kb_rows)
        if kb_header:
            all_kb_header = kb_header

        # Spot-check for inactive accessions (UniParc fallback demo)
        sample_accs = [r.get("Entry", "") for r in kb_rows if r.get("Entry")][:5]
        print(f"\n  UniParc fallback check on first {len(sample_accs)} accessions …")
        for acc in sample_accs:
            check_uniparc_fallback(acc)

        # ── B. FASTA ───────────────────────────────────────────────────────
        fasta_block = fetch_fasta(org_name, tax_id)
        all_fasta  += fasta_block

        # ── C. Proteomics metadata ─────────────────────────────────────────
        print(f"\n  Checking Proteins API coverage for TaxID {tax_id} …")
        pipelines = check_proteomics_coverage(tax_id)
        print(f"  Coverage: {pipelines}")

        org_nonptm, org_ptm = [], []
        if pipelines.get("nonptm"):
            org_nonptm = fetch_proteomics_nonptm(org_name, tax_id)
        else:
            print(f"  No non-PTM data available for {org_name}")

        if pipelines.get("ptm"):
            org_ptm = fetch_proteomics_ptm(org_name, tax_id)
        else:
            print(f"  No PTM data available for {org_name}")

        all_nonptm.extend(org_nonptm)
        all_ptm.extend(org_ptm)

        summary_lines += [
            f"── {org_name} (TaxID {tax_id})",
            f"   KB entries       : {len(kb_rows):,}",
            f"   FASTA sequences  : {fasta_block.count('>'):,}",
            f"   Non-PTM records  : {len(org_nonptm):,}",
            f"   PTM records      : {len(org_ptm):,}",
            "",
        ]

    # ── Save all outputs ───────────────────────────────────────────────────
    print(f"\n{'═'*65}")
    print("  Saving outputs …")

    if all_kb_rows and all_kb_header:
        save_tsv(all_kb_rows, all_kb_header,
                 OUTPUT_DIR / "uniprot_kb_entries.tsv")

    fasta_path = OUTPUT_DIR / "sequences.fasta"
    fasta_path.write_text(all_fasta, encoding="utf-8")
    print(f"  ✓ FASTA → sequences.fasta ({all_fasta.count('>'):,} sequences)")

    nonptm_fields = [
        "organism", "tax_id", "accession", "peptide",
        "length_aa", "approx_mass_Da",
        "pep_start", "pep_end", "unique_to_isoform",
        "dataset_id", "psm_count", "pipeline",
    ]
    save_tsv(all_nonptm, nonptm_fields,
             OUTPUT_DIR / "proteomics_nonptm_peptides.tsv")

    ptm_fields = [
        "organism", "tax_id", "accession", "peptide",
        "length_aa", "approx_mass_Da",
        "pep_start", "pep_end", "ptm_pos_in_entry",
        "modification", "confidence", "final_site_prob",
        "psm_count", "dataset_id", "pipeline",
    ]
    save_tsv(all_ptm, ptm_fields,
             OUTPUT_DIR / "proteomics_ptm_peptides.tsv")

    # ── Global summary ─────────────────────────────────────────────────────
    all_kb_acc  = {r.get("Entry", "") for r in all_kb_rows}
    all_exp_acc = {r["accession"] for r in all_nonptm} | \
                  {r["accession"] for r in all_ptm}
    confirmed   = all_kb_acc & all_exp_acc

    summary_lines += [
        "══ GLOBAL TOTALS ══",
        f"  KB entries total             : {len(all_kb_rows):,}",
        f"  FASTA sequences              : {all_fasta.count('>'):,}",
        f"  Non-PTM peptide records      : {len(all_nonptm):,}",
        f"  PTM records                  : {len(all_ptm):,}",
        f"  KB accessions with any",
        f"    proteomics evidence        : {len(confirmed):,}",
        "",
        "Output files",
        "  output/uniprot_kb_entries.tsv",
        "  output/sequences.fasta",
        "  output/proteomics_nonptm_peptides.tsv",
        "  output/proteomics_ptm_peptides.tsv",
        "  output/pipeline_summary.txt",
    ]
    write_summary(summary_lines, OUTPUT_DIR / "pipeline_summary.txt")

    print("\n" + "═" * 65)
    print("  DONE  –  Global totals")
    print(f"  KB entries total             : {len(all_kb_rows):,}")
    print(f"  FASTA sequences              : {all_fasta.count('>'):,}")
    print(f"  Non-PTM peptide records      : {len(all_nonptm):,}")
    print(f"  PTM records                  : {len(all_ptm):,}")
    print(f"  Accessions w/ prot. evidence : {len(confirmed):,}")
    print("═" * 65)
    print(f"\n  All files written to: {OUTPUT_DIR.resolve()}\n")


if __name__ == "__main__":
    main()
