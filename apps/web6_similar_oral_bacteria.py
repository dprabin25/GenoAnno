"""
6_similar_oral_bacteria.py

Streamlit tab: find oral bacteria with similar gene composition by combining
FOUR processed summary tables derived from three raw annotation files:

  Raw file                          -> Processed table(s)
  --------------------------------------------------------------------------
  KBASE annotations (.tsv)          -> Table 2: KEGG pathway gene counts
                                    -> Table 3: Protein function family gene counts
  Proksee/Bakta annotations (.tsv)  -> Table 1: Functional category gene counts
  Products.tsv (KEGG-Decoder style) -> Table 4: Metabolic pathway completeness

"""

import io
import re
import pandas as pd
import streamlit as st

# ---------------------------------------------------------------------------
# Combined prompt template
# ---------------------------------------------------------------------------

DEFAULT_INPUT_DESCRIPTION = """I annotated genes from a bacterial genome (an oral bacterium) using multiple annotation pipelines and summarized the results into four processed tables describing its functional gene composition.

Table 1 - Functional category gene counts (Bakta):
I selected the Product column from the Bakta annotation and counted how many genes map to each repeated product/function term.
Column 1: Functional category
Column 2: Gene count

Table 2 - KEGG pathway gene counts (KBASE):
I assigned each gene to its corresponding KEGG pathway and counted the number of genes mapped to each pathway.
Column 1: KEGG pathway
Column 2: Gene count

Table 3 - Protein function family gene counts (KBASE):
I assigned each gene to its corresponding protein function family and counted the number of genes mapped to each family.
Column 1: Protein function family
Column 2: Gene count

Table 4 - Metabolic pathway completeness (KEGG-Decoder style):
For each metabolic pathway, I calculated its completeness and, where applicable, determined whether the pathway is present (TRUE/FALSE).
Column 1: Metabolic pathway name
Column 2: Pathway status (completeness score or TRUE/FALSE)
Filtering rule already applied: pathways with 0, FALSE, empty, missing, NA, or NaN values were removed before analysis.

Only these four processed tables (not the raw annotation files) should be used for phenotype grouping and comparison.

Question: Are there other oral bacteria that have similar gene composition to my bacteria? List name of well characterized bacteria with the most consistent shared functional categories. Describe the phenotype of these bacteria. Provide only a table with columns containing a list of closely related bacteria with their shared functions and phenotypes and a summary. 
Rank the species and provide the reasoning. Sort it based on the confidence, with higher confidence listed on the top. 

"""

# ---------------------------------------------------------------------------
# Parsers: raw file -> processed table
# ---------------------------------------------------------------------------

def parse_bakta_functional_categories(bakta_file) -> pd.DataFrame:
    """Proksee/Bakta annotation TSV -> Table 1 (Functional category, Gene count).

    Bakta TSVs have a handful of '#'-prefixed metadata lines before the header
    row (which itself starts with '#Sequence Id'). We skip the metadata lines
    but keep the header.
    """
    raw = bakta_file.read()
    if isinstance(raw, bytes):
        raw = raw.decode("utf-8", errors="replace")
    lines = raw.splitlines()
    header_idx = next(i for i, l in enumerate(lines) if l.startswith("#Sequence Id"))
    data_str = "\n".join(lines[header_idx:])
    df = pd.read_csv(io.StringIO(data_str), sep="\t")
    df.columns = [c.lstrip("#").strip() for c in df.columns]

    products = df["Product"].dropna()
    products = products[products.str.strip() != ""]
    counts = (
        products.value_counts()
        .rename_axis("Functional category")
        .reset_index(name="Gene count")
    )
    return counts


def parse_kbase_kegg_pathways(kbase_file) -> pd.DataFrame:
    """KBASE annotation TSV -> Table 2 (KEGG pathway, Gene count).

    Uses the 'kegg_hit' description tied to each gene's ko_id. Swap this for
    a ko_id -> KEGG pathway lookup if you have one; as shipped it counts genes
    per unique KEGG hit description, which is the closest pathway-level label
    available directly in this file.
    """
    df = pd.read_csv(kbase_file, sep="\t")
    hits = df["kegg_hit"].dropna()
    hits = hits[hits.str.strip() != ""]
    counts = (
        hits.value_counts()
        .rename_axis("KEGG pathway")
        .reset_index(name="Gene count")
    )
    return counts


def parse_kbase_protein_families(kbase_file) -> pd.DataFrame:
    """KBASE annotation TSV -> Table 3 (Protein function family, Gene count)."""
    df = pd.read_csv(kbase_file, sep="\t")
    fams = df["peptidase_family"].dropna()
    fams = fams[fams.str.strip() != ""]
    counts = (
        fams.value_counts()
        .rename_axis("Protein function family")
        .reset_index(name="Gene count")
    )
    return counts


def parse_products_completeness(products_file) -> pd.DataFrame:
    """Products.tsv (one row per genome, one column per pathway) -> Table 4
    (Metabolic pathway name, Pathway status), with the FALSE/0/empty/NA rows
    already filtered out.
    """
    df = pd.read_csv(products_file, sep="\t")
    row = df.iloc[0]
    row = row.drop(labels=["genome"], errors="ignore")
    long = row.rename_axis("Metabolic pathway name").reset_index(name="Pathway status")

    def is_present(v):
        if pd.isna(v):
            return False
        if isinstance(v, str):
            v_clean = v.strip().lower()
            if v_clean in ("", "na", "nan", "false", "0"):
                return False
            return True
        if isinstance(v, bool):
            return v
        try:
            return float(v) != 0
        except (TypeError, ValueError):
            return bool(v)

    long = long[long["Pathway status"].apply(is_present)].reset_index(drop=True)
    return long


# ---------------------------------------------------------------------------
# Prompt assembly
# ---------------------------------------------------------------------------

def _to_markdown_table(df: pd.DataFrame) -> str:
    try:
        return df.to_markdown(index=False)
    except ImportError:
        # tabulate not installed - fall back to CSV-ish text
        return df.to_csv(index=False, sep="|")


def build_combined_prompt(
    functional_category_df: pd.DataFrame,
    kegg_pathway_df: pd.DataFrame,
    protein_family_df: pd.DataFrame,
    pathway_completeness_df: pd.DataFrame,
    extra_instructions: str = "",
) -> str:
    prompt = DEFAULT_INPUT_DESCRIPTION.strip() + "\n\n"
    prompt += "Table 1 - Functional category gene counts:\n"
    prompt += _to_markdown_table(functional_category_df) + "\n\n"
    prompt += "Table 2 - KEGG pathway gene counts:\n"
    prompt += _to_markdown_table(kegg_pathway_df) + "\n\n"
    prompt += "Table 3 - Protein function family gene counts:\n"
    prompt += _to_markdown_table(protein_family_df) + "\n\n"
    prompt += "Table 4 - Metabolic pathway completeness:\n"
    prompt += _to_markdown_table(pathway_completeness_df) + "\n\n"
    if extra_instructions.strip():
        prompt += "Additional instructions:\n" + extra_instructions.strip() + "\n"
    return prompt


# ---------------------------------------------------------------------------
# Streamlit page
# ---------------------------------------------------------------------------

def render():
    st.header("Find Similar Oral Bacteria")
    st.caption(
        "Upload your three raw annotation files. This tab derives four "
        "summary tables from them (functional categories, KEGG pathways, "
        "protein families, pathway completeness) and combines all four into "
        "one prompt to search for bacteria with a similar gene composition."
    )

    col1, col2, col3 = st.columns(3)
    with col1:
        bakta_file = st.file_uploader(
            "Bakta/Proksee annotation (.tsv)", type=["tsv"], key="bakta_upload"
        )
    with col2:
        kbase_file = st.file_uploader(
            "KBASE annotation (.tsv)", type=["tsv"], key="kbase_upload"
        )
    with col3:
        products_file = st.file_uploader(
            "Pathway completeness table, e.g. Products.tsv (.tsv)",
            type=["tsv"],
            key="products_upload",
        )

    input_description = st.text_area(
        "Prompt sent to the model",
        value=DEFAULT_INPUT_DESCRIPTION,
        height=380,
    )

    extra_instructions = st.text_area(
        "Any additional instructions (optional)", value="", height=80
    )

    if st.button("Find similar bacteria", type="primary"):
        if not (bakta_file and kbase_file and products_file):
            st.warning("Please upload all three files first.")
            return

        with st.spinner("Processing annotation files..."):
            functional_category_df = parse_bakta_functional_categories(bakta_file)
            kbase_file.seek(0)
            kegg_pathway_df = parse_kbase_kegg_pathways(kbase_file)
            kbase_file.seek(0)
            protein_family_df = parse_kbase_protein_families(kbase_file)
            pathway_completeness_df = parse_products_completeness(products_file)

        with st.expander("Processed tables (review before sending to the model)"):
            st.subheader("Table 1 - Functional category gene counts")
            st.dataframe(functional_category_df)
            st.subheader("Table 2 - KEGG pathway gene counts")
            st.dataframe(kegg_pathway_df)
            st.subheader("Table 3 - Protein function family gene counts")
            st.dataframe(protein_family_df)
            st.subheader("Table 4 - Metabolic pathway completeness")
            st.dataframe(pathway_completeness_df)

        # Use the (possibly edited) text from the textbox as the description,
        # but still append the freshly parsed tables.
        base_description = input_description if input_description.strip() else DEFAULT_INPUT_DESCRIPTION
        prompt = base_description.strip() + "\n\n"
        prompt += "Table 1 - Functional category gene counts:\n" + _to_markdown_table(functional_category_df) + "\n\n"
        prompt += "Table 2 - KEGG pathway gene counts:\n" + _to_markdown_table(kegg_pathway_df) + "\n\n"
        prompt += "Table 3 - Protein function family gene counts:\n" + _to_markdown_table(protein_family_df) + "\n\n"
        prompt += "Table 4 - Metabolic pathway completeness:\n" + _to_markdown_table(pathway_completeness_df) + "\n\n"
        if extra_instructions.strip():
            prompt += "Additional instructions:\n" + extra_instructions.strip() + "\n"

        with st.spinner("Querying the model..."):
            # TODO: replace this with the same LLM-calling function your other
            # tabs (1-5) already use, e.g.:
            #   response = call_llm(prompt)
            # Left as a placeholder so this file runs standalone.
            response = "TODO: wire up your existing LLM call here.\n\nPrompt built:\n\n" + prompt

        st.subheader("Result")
        st.markdown(response)


if __name__ == "__main__":
    render()
