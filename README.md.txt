# 🧬 CNV HGVS Annotator with Cytoband

A **Streamlit web app** for generating HGVS-style annotations for **Copy Number Variants (CNVs)** — such as deletions and duplications — enriched with **cytoband information**.

---

## 🔍 Overview

This app allows users to input genomic coordinates and CNV event types (deletion or duplication), and it automatically generates **HGVS-like notations**.  
It also integrates **cytoband data** (from a user-uploaded `cytoBand.txt` file) to add cytogenetic-level annotation.

---

## ⚙️ Features

- 🧩 Parse genomic coordinates (e.g., `chr16:15489724-16367962`)
- 🧬 Identify cytoband regions using uploaded `cytoBand.txt`
- 🧾 Generate HGVS-like CNV notation (duplication/deletion)
- 💻 Streamlit-based interactive web interface

---

## 🧠 Example Input

| Input Type | Example |
|-------------|----------|
| **Coordinate** | `chr16:15489724-16367962` |
| **Event Type** | `duplication` or `deletion` |
| **Cytoband File** | `cytoBand.txt` (tab-separated: chrom, start, end, cytoband, stain) |

---

## 🧪 Example Output

```text
chr16:(?_15489724)_(16367962_?) [3]
(chr16p13.11-p12.3 partial duplication)
