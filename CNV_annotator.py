import pandas as pd
import re
import streamlit as st


class HGVSCNVAnnotator:
    def __init__(self, cytoband_file=None, genelist_file=None):
        self.coordinate_pattern = r'^(chr)?(\d+|X|Y|M|MT):(\d+)-(\d+)$'
        self.cytoband_df = None
        self.genelist_df = None

        if cytoband_file:
            self.load_cytoband(cytoband_file)
        if genelist_file:
            self.load_genelist(genelist_file)

    def load_cytoband(self, cytoband_file):
        try:
            self.cytoband_df = pd.read_csv(cytoband_file, header=None)
            self.cytoband_df.columns = ['chrom', 'start', 'end', 'cytoband', 'stain']
            print(f"Loaded {len(self.cytoband_df)} cytoband entries")
        except Exception as e:
            print(f"Error loading cytoband file: {e}")

    def load_genelist(self, genelist_file):
        try:
            self.genelist_df = pd.read_csv(genelist_file, header=None)
            self.genelist_df.columns = ['chrom', 'start', 'end', 'gene_name']
            self.genelist_df['chrom'] = self.genelist_df['chrom'].astype(str)
            print(f"Loaded {len(self.genelist_df)} gene entries")
        except Exception as e:
            print(f"Error loading genelist file: {e}")

    def parse_coordinate(self, coordinate):
        match = re.match(self.coordinate_pattern, coordinate.strip())
        if not match:
            return None
        chrom = match.group(2)
        start = int(match.group(3))
        end = int(match.group(4))
        if start >= end:
            print(f"Error: Start position ({start}) must be less than end position ({end})")
            return None
        return chrom, start, end

    def get_overlapping_genes(self, chrom, start, end):
        if self.genelist_df is None:
            return []
        chrom_normalize = str(chrom).replace('chr', '')
        chr_genes = self.genelist_df[
            self.genelist_df['chrom'].str.replace('chr', '') == chrom_normalize
        ].copy()
        if chr_genes.empty:
            return []
        overlapping = chr_genes[
            (chr_genes['start'] <= end) & (chr_genes['end'] >= start)
        ].copy()
        if overlapping.empty:
            return []
        genes = sorted(overlapping['gene_name'].unique().tolist())
        return genes

    
    def get_cytoband(self, chrom, start, end):
        if self.cytoband_df is None:
            return None
        chrom_str = f"chr{chrom}"
        chr_data = self.cytoband_df[self.cytoband_df['chrom'] == chrom_str].copy()
        if chr_data.empty:
            return None
        overlapping = chr_data[
            (chr_data['start'] <= end) & (chr_data['end'] >= start)
        ].copy()
        if overlapping.empty:
            return None
        cytobands = overlapping['cytoband'].unique()
        if len(cytobands) == 1:
            return cytobands[0]
        first_band = cytobands[0]
        last_band = cytobands[-1]
        arm_prefix = first_band[0] if first_band[0] in ['p', 'q'] else ""
        if arm_prefix and last_band.startswith(arm_prefix):
            last_band_clean = last_band[1:]
            return f"{first_band}-{last_band_clean}"
        else:
            return f"{first_band}-{last_band}"
            
        cytoband = self.get_cytoband(chrom, start, end)

    def generate_hgvs(self, coordinate, event_type, zygosity=None):
        parsed = self.parse_coordinate(coordinate)
        if not parsed:
            return None, None, []
        chrom, start, end = parsed
        event_type = event_type.lower().strip()
        base_notation = f"chr{chrom}:(?_{start})_({end}_?)"
        if event_type in ['duplication', 'dup']:
            if zygosity:
                zygosity=zygosity.lower().strip()
                if zygosity in ['homozygous', 'hom']:
                    copy_number="[4]"
                elif zygosity in ['heterozygous', 'het']:
                    copy_number="[3]"
                else:
                    copy_number = "[3]"
            else:
                copy_number = "[3]"            
            hgvs_notation = f"{base_notation} {copy_number}"
            event_name = "duplication"
        elif event_type in ['deletion', 'del']:
            hgvs_notation = f"{base_notation}del"
            event_name = "deletion"
        else:
            print(f"Error: Invalid event type '{event_type}'.")
            return None, None, []
        cytoband = self.get_cytoband(chrom, start, end)
        genes = self.get_overlapping_genes(chrom, start, end)
        if cytoband:
            full_annotation = f"{hgvs_notation}\n(chr{chrom}{cytoband} partial {event_name})"
        else:
            full_annotation = hgvs_notation
        return hgvs_notation, full_annotation, genes


# ---------------------------
# ✅ Streamlit UI starts here
# ---------------------------
st.markdown(
    """
    <style>
    .stApp {
        background: radial-gradient(ellipse at 80% 1200%, #1de9f6 0%, #0d47a1 50%, #020924 100%);
        background-blend-mode: lighten;
    }
    /* Annotate button */
    div.stButton > button {
        background-color: #28a745 !important;
        color: white !important;
        border: none !important;
        padding: 0.55rem 1rem;
        font-weight: 700;
    }
    /* Custom labels for inputs */
    .custom-label {
        font-size: 18px;
        font-weight: 600;
        color: #FFFFFF;

        margin-bottom: 6px;
        display: block;
    }
    /* Placeholder style inside the text input */
    div[data-testid="stTextInput"] input::placeholder {
        font-size: 15px;
        color: #cfeff6;

        opacity: 0.95;
    }
    /* Warning card */
    .warning-card {
        background: linear-gradient(90deg,#fffbe6,#fff3cd);

        color: #665100;

        padding: 10px 12px;
        border-radius: 8px;
        font-weight: 600;
        margin-bottom: 10px;
    }
    /* Result badge */
    .result-badge {
        padding:10px 14px;
        border-radius:10px;
        color:white;
        background: linear-gradient(90deg, #8900fa 0%, #188afc 100%);

        font-weight:700;
        display:inline-block;
        margin-bottom:8px;
    }
    /* Result area card */
    .result-card {
        background: rgba(255,255,255,0.03);
        padding:12px;
        border-radius:8px;
        color:#e6f7ff;
    }
    </style>
    """,
    unsafe_allow_html=True
)

st.title("GRCh37 CNV Annotation Tool")
st.markdown(
     """
    **⚠️ Important:**  
    This tool is designed exclusively for **multigene CNV annotation**.  
    It should **not** be used for **exonic events** or **whole-chromosome duplications/deletions**.
    """,
    unsafe_allow_html=True
)
annotator = HGVSCNVAnnotator("cytoBand.csv", "Genelist.csv")
st.markdown("<label class='custom-label'>Enter coordinate</label>", unsafe_allow_html=True)
coordinate = st.text_input("Enter coordinate (e.g., chr16:15489724-16367962)", label_visibility="collapsed")
st.markdown("<label class='custom-label'>Select zygosity</label>", unsafe_allow_html=True)
zygosity=st.selectbox("Select zygosity", ["Homozygous", "Heterozygous"], label_visibility="collapsed")
st.markdown("<label class='custom-label'>Select event type</label>", unsafe_allow_html=True)
event_type = st.selectbox("Select event type", ["duplication", "deletion"], label_visibility="collapsed")

if st.button("Annotate"):
    hgvs, full, genes = annotator.generate_hgvs(coordinate, event_type)
    if full:
        lines = full.split("\n")
        st.text(lines[0])
        if len(lines) > 1:
            st.text(lines[1])
        st.subheader("Overlapping Genes:")
        if genes and len(genes) > 0:
            st.markdown(
            f"<div style='padding:12px; border-radius:10px; color:white; background: linear-gradient(90deg, #8900fa 0%, #188afc 100%); font-weight:bold; font-size:1.2em; width:fit-content;'>"
            f"Found {len(genes)} gene(s)"
            "</div>",
            unsafe_allow_html=True,
        )
        st.markdown("&nbsp;")    
        st.write(", ".join(genes))
        if len(genes) > 1:
            gene_df = pd.DataFrame({'Gene Name': genes})
            st.dataframe(gene_df, use_container_width=True)
        else:
            st.info("No overlapping genes found in this region.")
    else:
         st.error("Invalid input or no cytoband found.")

st.markdown(
    "<div style='color:#dff6ff; margin-top:14px; font-size:13px;'>Tip: Use the format <code>chrN:start-end</code>. Example: <code>chr16:15489724-16367962</code></div>",
    unsafe_allow_html=True,
)

