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
            self.cytoband_df = pd.read_table(cytoband_file, header=None)
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

    def generate_hgvs(self, coordinate, event_type):
        parsed = self.parse_coordinate(coordinate)
        if not parsed:
            return None, None, []
        chrom, start, end = parsed
        event_type = event_type.lower().strip()
        base_notation = f"chr{chrom}:(?_{start})_({end}_?)"
        if event_type in ['duplication', 'dup']:
            hgvs_notation = f"{base_notation} [3]"
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
st.title("HGVS CNV Annotator")
st.markdown(
     """
    **⚠️ Important:**  
    This tool is designed exclusively for **multigene CNV analysis**.  
    It should **not** be used for **single-exon events** or **whole-chromosome duplications/deletions**.
    """,
    unsafe_allow_html=True
)
annotator = HGVSCNVAnnotator("cytoBand.txt", "Genelist.csv")

coordinate = st.text_input("Enter coordinate (e.g., chr16:15489724-16367962)")
event_type = st.selectbox("Select event type", ["duplication", "deletion"])

if st.button("Annotate"):
    hgvs, full, genes = annotator.generate_hgvs(coordinate, event_type)

    if full:
        lines = full.split("\n")
        st.text(lines[0])
        if len(lines) > 1:
            st.text(lines[1])

        st.subheader("Overlapping Genes:")
        if genes and len(genes) > 0:
            st.success(f"Found {len(genes)} gene(s)")
            st.write(", ".join(genes))
            if len(genes) > 5:
                gene_df = pd.DataFrame({'Gene Name': genes})
                st.dataframe(gene_df, use_container_width=True)
        else:
            st.info("No overlapping genes found in this region.")
    else:
        st.error("Invalid input or no cytoband found.")
