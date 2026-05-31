import streamlit as st
from Bio import SeqIO
from collections import Counter
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import io
import re

# -----------------------------------------------------------------------------
# 1. PAGE CONFIGURATION & THEME INITIALIZATION
# -----------------------------------------------------------------------------
st.set_page_config(
    page_title="BioSight Analytics Engine",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS for polished typography and metric cards
st.markdown("""
    <style>
    .main-header { font-size: 2.2rem; font-weight: 800; color: #1E3A8A; margin-bottom: 0.5rem; }
    .sub-header { font-size: 1.2rem; color: #4B5563; margin-bottom: 2rem; }
    .metric-card { 
        background-color: #F3F4F6; 
        padding: 1rem; 
        border-radius: 0.5rem; 
        border-left: 5px solid #2563EB;
        margin-bottom: 1rem;
    }
    .metric-title { font-size: 0.85rem; font-weight: 600; color: #4B5563; text-transform: uppercase; }
    .metric-value { font-size: 1.6rem; font-weight: 700; color: #1F2937; }
    </style>
""", unsafe_allow_html=True)

st.markdown('<div class="main-header">🧬 BioSight: Clinical Variant & Sequence Motif Analytics Engine</div>', unsafe_allow_html=True)
st.markdown('<div class="sub-header">Diagnostic sequence processing and visualization pipeline</div>', unsafe_allow_html=True)

# -----------------------------------------------------------------------------
# 2. SIDEBAR NAVIGATION & CONTROLS
# -----------------------------------------------------------------------------
st.sidebar.header("🛠️ Pipeline Configurations")

# K-mer slider configuration
k_length = st.sidebar.slider(
    label="K-mer Window Length",
    min_value=2,
    max_value=6,
    value=3,
    step=1,
    help="Select the size of overlapping nucleotide chunks to analyze."
)

# Complement toggle configuration
enable_antisense = st.sidebar.checkbox(
    label="Analyze Antisense Strand",
    value=False,
    help="Generate complementary sequence statistics alongside sense data."
)

# -----------------------------------------------------------------------------
# 3. FILE INGESTION & ROBUST DATA DECODING (MODULE A)
# -----------------------------------------------------------------------------
uploaded_file = st.file_uploader(
    label="Upload Clinical Genomic Asset (.fasta, .fna)", 
    type=["fasta", "fna"],
    help="Upload a standardized genomic sequence file to initiate the diagnostic pipeline."
)

if uploaded_file is not None:
    try:
        # Step A1: Bypass BioPython StreamModeError via memory decoding translation
        file_bytes = uploaded_file.read().decode("utf-8")
        text_stream = io.StringIO(file_bytes)
        
        # Step A2: Parse all available records within the sequence file
        all_records = list(SeqIO.parse(text_stream, "fasta"))
        
        if not all_records:
            st.error("Data Parsing Failure: The uploaded asset contains no verifiable FASTA records.")
            st.stop()
            
        # Step A3: Dynamic Record Selection handling multiple contigs/chromosomes
        record_dict = {f"{rec.id} | {rec.description[:40]}...": rec for rec in all_records}
        
        if len(all_records) > 1:
            st.sidebar.info(f"Multi-record asset detected. Total contigs: {len(all_records)}")
            selected_key = st.sidebar.selectbox("Active Sequence Record Target", list(record_dict.keys()))
            active_record = record_dict[selected_key]
        else:
            active_record = all_records[0]
            
        # Step A4: Strict validation, case stabilization, and molecular assessment
        raw_sequence = str(active_record.seq).upper()
        sequence_length = len(raw_sequence)
        
        if sequence_length == 0:
            st.error("Sequence Boundary Error: The selected genomic record contains zero bases.")
            st.stop()
            
        # Calculate approximate molecular mass (Average single-stranded DNA nucleotide ~330.79 g/mol)
        approx_mw_g_mol = sequence_length * 330.79
        
        # Display Core High-Level Metadata Layout Dashboard
        st.markdown("### 📊 Primary Metadata Overview")
        col1, col2, col3, col4 = st.columns(4)
        
        with col1:
            st.markdown(f'<div class="metric-card"><div class="metric-title">Target Sequence ID</div><div class="metric-value">{active_record.id}</div></div>', unsafe_allow_html=True)
        with col2:
            st.markdown(f'<div class="metric-card"><div class="metric-title">Total Sequence Length</div><div class="metric-value">{sequence_length:,} bp</div></div>', unsafe_allow_html=True)
        with col3:
            # Placeholder calculated globally down in Module B
            gc_metric_container = st.empty()
        with col4:
            st.markdown(f'<div class="metric-card"><div class="metric-title">Est. Molecular Weight</div><div class="metric-value">{approx_mw_g_mol/1e6:.2f} M g/mol</div></div>', unsafe_allow_html=True)
            
        st.text_area("Record Header Information", value=active_record.description, disabled=True)
        
        # -----------------------------------------------------------------------------
        # 4. NUCLEOTIDE DATA SCIENCE LAYERS & FORMULAS (MODULE B)
        # -----------------------------------------------------------------------------
        st.markdown("---")
        st.markdown("### 🔢 Nucleotide Profiling & Global Composition Statistics")
        
        # Comprehensive counting including ambiguous degenerate biological base codes
        base_counts = Counter(raw_sequence)
        standard_bases = ['A', 'T', 'C', 'G']
        total_standard = sum(base_counts[base] for base in standard_bases)
        total_ambiguous = sum(base_counts[base] for base in base_counts if base not in standard_bases)
        
        # Accurate Global GC Calculation Formula deployment
        g_count = base_counts.get('G', 0)
        c_count = base_counts.get('C', 0)
        total_bases = sum(base_counts.values())
        
        global_gc_percentage = ((g_count + c_count) / total_bases) * 100 if total_bases > 0 else 0.0
        
        # Deliver global metric up into card container
        gc_metric_container.markdown(f'<div class="metric-card"><div class="metric-title">Global GC Content</div><div class="metric-value">{global_gc_percentage:.2f} %</div></div>', unsafe_allow_html=True)
        
        # Build pandas interactive matrix
        base_data_matrix = []
        for base, count in base_counts.items():
            frequency = (count / total_bases) * 100
            base_data_matrix.append({"Nucleotide Base": base, "Absolute Count": count, "Frequency (%)": round(frequency, 4)})
            
        df_bases = pd.DataFrame(base_data_matrix).sort_values(by="Absolute Count", ascending=False).reset_index(drop=True)
        
        col_b1, col_b2 = st.columns([1, 1])
        with col_b1:
            st.markdown("**Composition Matrix Table**")
            st.dataframe(df_bases, use_container_width=True)
            st.metric("Total Ambiguous/Degenerate Bases Found (e.g., N)", f"{total_ambiguous:,} bp")
            
        with col_b2:
            st.markdown("**Relative Base Abundance Vector View**")
            fig_base, ax_base = plt.subplots(figsize=(6, 4))
            ax_base.bar(
                df_bases["Nucleotide Base"], 
                df_bases["Absolute Count"], 
                color=['#2563EB', '#DC2626', '#16A34A', '#D97706'] + ['#9CA3AF']*len(df_bases),
                edgecolor='black',
                linewidth=0.7
            )
            ax_base.set_ylabel("Absolute Counts", fontsize=10)
            ax_base.set_xlabel("Identified Bases", fontsize=10)
            ax_base.grid(axis='y', linestyle='--', alpha=0.5)
            plt.tight_layout()
            st.pyplot(fig_base)
            plt.close(fig_base)

        # -----------------------------------------------------------------------------
        # 5. OPTIMIZED SLIDING-WINDOW K-MER ANALYSIS ENGINE (MODULE C)
        # -----------------------------------------------------------------------------
        st.markdown("---")
        st.markdown(f"### 📈 Top 10 Dominant {k_length}-mer Local Distributions")
        
        # Optimization: Restrict profiling array memory load on extreme sequence vectors
        if sequence_length > 5000000:
            st.warning("Computational Warning: Sequence exceeds 5Mbp boundary. Sampling first 5Mbp for sliding-window k-mer calculation efficiency.")
            kmer_target_seq = raw_sequence[:5000000]
        else:
            kmer_target_seq = raw_sequence
            
        # Fast iterative sliding window extraction allocation
        kmer_dict = {}
        for i in range(len(kmer_target_seq) - k_length + 1):
            kmer = kmer_target_seq[i:i + k_length]
            kmer_dict[kmer] = kmer_dict.get(kmer, 0) + 1
            
        # Extract and sort the top 10 dominant structural motifs
        df_kmers = pd.DataFrame(list(kmer_dict.items()), columns=["K-mer Substring", "Observation Count"])
        df_top_kmers = df_kmers.sort_values(by="Observation Count", ascending=False).head(10).reset_index(drop=True)
        
        col_c1, col_c2 = st.columns([1, 2])
        with col_c1:
            st.write(f"Top 10 Detected Frequency Rankings ({k_length}-length chunks)")
            st.dataframe(df_top_kmers, use_container_width=True)
            
        with col_c2:
            fig_kmer, ax_kmer = plt.subplots(figsize=(8, 3.8))
            ax_kmer.bar(
                df_top_kmers["K-mer Substring"], 
                df_top_kmers["Observation Count"], 
                color='#0D9488', 
                edgecolor='black',
                linewidth=0.8
            )
            ax_kmer.set_xlabel(f"Extracted {k_length}-mer Structure", fontsize=10)
            ax_kmer.set_ylabel("Frequency Counts", fontsize=10)
            ax_kmer.grid(axis='y', linestyle=':', alpha=0.6)
            
            # Clean tick mark label adjustments prevent text collision strings
            if k_length >= 4:
                plt.xticks(rotation=45, ha='right')
                
            plt.tight_layout()
            st.pyplot(fig_kmer)
            plt.close(fig_kmer)

        # -----------------------------------------------------------------------------
        # 6. VECTORIZED ADVANCED PATTERN & MOTIF DISCOVERY (MODULE D)
        # -----------------------------------------------------------------------------
        st.markdown("---")
        st.markdown("### 🔍 High-Performance Target Variant & Motif Search Engine")
        
        search_pattern = st.text_input(
            label="Input Target Structural Motif Sequence Pattern (IUPAC Strict Match)", 
            value="ATG",
            help="Enter a sequence motif pattern to look up. Example: Enter ATG to isolate transcription init sites."
        ).upper().strip()
        
        if search_pattern:
            # Strict regex structural input validation verification check
            if not re.match(r"^[ATGCNRYKMSWBDHV]+$", search_pattern):
                st.error("IUPAC Validation Refusal: Search input array query contains illegal character flags. Restrict syntax strings explicitly to standard or degenerate nucleotide labels.")
            else:
                st.info(f"Initiating parsing scanner string array matching for target motif: '{search_pattern}'")
                
                # Regex execution engine applying lookahead to extract overlapping sub-string coordinates safely
                # Expression handles overlapping patterns securely via lookahead grouping syntax assertions
                regex_compiled_query = f"(?=({search_pattern}))"
                
                matches = [m.start() for m in re.finditer(regex_compiled_query, raw_sequence)]
                total_matches = len(matches)
                
                st.metric("Total Motifs Discovered Across Active Strand", f"{total_matches:,} matches")
                
                if total_matches > 0:
                    results_list = []
                    
                    # Optimization: Cap reporting views array allocation limits inside display tables
                    max_render_limit = min(total_matches, 1000)
                    
                    for idx in matches[:max_render_limit]:
                        # Dynamically slice contextual local flanking coordinates
                        start_upstream = max(0, idx - 10)
                        end_downstream = min(sequence_length, idx + len(search_pattern) + 10)
                        
                        upstream_segment = raw_sequence[start_upstream:idx]
                        match_segment = raw_sequence[idx:idx + len(search_pattern)]
                        downstream_segment = raw_sequence[idx + len(search_pattern):end_downstream]
                        
                        flanking_visualization = f"... {upstream_segment} -> [{match_segment}] <- {downstream_segment} ..."
                        
                        results_list.append({
                            "Absolute Position Index (0-based)": idx,
                            "Chromosomal Position Display (1-based)": idx + 1,
                            "Local Flanking Sequence Architecture Window (±10 bp)": flanking_visualization
                        })
                        
                    df_search_results = pd.DataFrame(results_list)
                    st.dataframe(df_search_results, use_container_width=True)
                    
                    if total_matches > 1000:
                        st.caption("⚠️ Layout Optimization: Table rendering output limits truncated at 1,000 matches to protect memory pipelines.")
                else:
                    st.warning(f"Diagnostic Scan Negative: No structural instances of motif '{search_pattern}' discovered inside the sequence index.")

        # -----------------------------------------------------------------------------
        # 7. ANTISENSE STRAND STRUCTURAL TRANSFORMATION LAYERS (MODULE E)
        # -----------------------------------------------------------------------------
        if enable_antisense:
            st.markdown("---")
            st.markdown("### 🔄 Antisense Coordinate Mapping & Complementary Architecture Matrix")
            
            # Unshakeable Watson-Crick structural conversion map
            complement_map = str.maketrans("ATCGN", "TAGCN")
            
            st.markdown("**Antisense Real-time Molecular Processing Sequence Information**")
            
            # Step E1: Generate the reverse complement strand safely string flipped
            with st.spinner("Reversing and calculating global matrix transformations..."):
                reverse_complement_string = raw_sequence.translate(complement_map)[::-1]
                rev_counts = Counter(reverse_complement_string)
                
                col_e1, col_e2 = st.columns([1, 2])
                with col_e1:
                    st.markdown("**Complementary Stand Properties**")
                    df_rev = pd.DataFrame(list(rev_counts.items()), columns=["Antisense Base", "Count"]).sort_values(by="Count", ascending=False)
                    st.dataframe(df_rev, use_container_width=True)
                
                with col_e2:
                    st.markdown("**Antisense Sequence Vector (Terminal Slice View)**")
                    st.text_area(
                        label="First 500bp Antisense Sequence Stream (5' ──► 3')", 
                        value=reverse_complement_string[:500] + ("..." if len(reverse_complement_string) > 500 else ""),
                        height=150,
                        disabled=True
                    )
                    
    except Exception as pipeline_exception:
        st.error(f"Critical System Execution Exception Trapped: An unhandled data exception halted processing operations.")
        st.error(f"System Message Context: {str(pipeline_exception)}")
        st.caption("Ensure your input asset coordinates contain legal formatting and standard uncorrupted text buffers.")

else:
    # Diagnostic standby screen execution layer shown when file missing context
    st.info("System Standby: Awaiting inbound file pipeline connection streams. Upload a validated genomic file to continue.")
