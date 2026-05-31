# Biosight-genomic-pipeline
This application is  aimed at enabling medical laboratory technicians and clinical researchers to input genomic FASTA / FNA data, extract the information in an efficient manner, perform nucleotide/k-mer analysis, execute an optimal substring search/pattern match algorithm, and create publication quality visualizations.


# 🧬 BioSight: Clinical Variant & Sequence Motif Analytics Engine

BioSight is a high-performance, production-grade clinical genomics analytics dashboard built natively in Python. It bridges the gap between raw sequencing outputs and clinical interpretation, providing medical laboratory technicians, researchers, and pathologists with a zero-configuration web interface to analyze genetic assets, track nucleotide distributions, and discover structural sequence motifs.


---

## 🛠️ Clinical Use Cases & Solved Roadblocks

Unlike general-purpose string parsers, BioSight is engineered around the strict file processing demands of clinical bioinformatics:

1.  **Bypassing Web Engine Stream Restrictions:** Built-in integration with `io.StringIO` decodes inbound binary data streams on-the-fly, completely neutralizing the traditional `StreamModeError` common in native web-based sequence parsing.
2.  **Multi-Contig/Chromosome Diagnostics:** Includes a dynamic record filter that handles complex, multi-record genomic files seamlessly, allowing users to isolate specific chromosomes from a single dashboard.
3. **Molecular Diagnostics (PCR Primer Design):** Deploys an accurate global $GC\%$ mathematical formula to calculate sequence stability, directly assisting in determining optimum melting temperatures for clinical assay design.
4.  **Overlapping Motif Identification:** Implements a high-performance regex engine applying lookahead assertions `(?=(...))` to capture overlapping sequence patterns (such as transcription initiation points or restriction enzyme cut sites) that traditional iteration loops miss.

## 🧬 Data Flow Pipeline Architecture

The application processes data across five isolated analytical modules:

[ User Browser ] ──► Uploads RAW Genomic Asset (.fna / .fasta)
│
▼
[ Ingestion & Memory Layer ] ──► Decodes stream & handles multi-contig selection
│
▼
[ Base Profiling Engine ] ──► Extracts global counts & runs strict GC% formula
│
▼
[ Sliding-Window K-mer Visualizer ] ──► Extracts top 10 dominant structural motifs
│
▼
[ Lookahead Motif Target Engine ] ──► Vectorized re.finditer flags absolute variants

## 💻 Tech Stack & Dependencies

This asset is written completely within the Python data science and structural biology ecosystem:
* **Core GUI UI Framework:** Streamlit
* **Bioinformatics Parser:** Biopython (SeqIO)
* **Data Structures & Arrays:** NumPy & Pandas
* **Visualization Engines:** Matplotlib

---

## 📈 Sample Dashboard Outputs

### 1. Primary Metadata & Global Composition View
When a file is loaded, the engine isolates contig headers, formats sequence lengths with proper localized typography (e.g., `4,639,675 bp`), calculates approximate molecular weight, and displays a normalized bar chart of sequence abundances.

### 2. Flanking Sequence Context Window
When a sequence motif is isolated, the engine doesn't just print indices; it isolates and extracts a $\pm 10\text{ bp}$ upstream and downstream flanking sequence boundary around the target locus to provide researchers with vital biological context.

---

## 📄 License

Distributed under the MIT License. See `LICENSE` for more information.

---
**Developed by Sreehari Nandanan** *Clinical Laboratory Technology & Computational Automation Specialist* [Connect with me on LinkedIn](YOUR_LINKEDIN_PROFILE_URL)

