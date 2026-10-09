# -*- coding: utf-8 -*-
"""
Created on Wed May 13 15:38:41 2026

@author: newfaculty
"""

import os
from pathlib import Path

import pandas as pd
import streamlit as st
from openai import OpenAI


# ============================================================
# Read config.txt from the same folder as this app
# ============================================================
HERE = Path(__file__).resolve().parent
CONFIG_TXT = HERE / "config.txt"


def read_config(path):
    config = {}

    # config.txt is optional in the website version.
    # OpenAI settings are normally provided from the main website form
    # and stored only in st.session_state for the current browser session.
    if not path.exists():
        return config

    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()

        if not line or line.startswith("#"):
            continue

        if "=" not in line:
            continue

        key, value = line.split("=", 1)
        config[key.strip().upper()] = value.strip()

    return config


def clean_filename(name):
    name = name.strip()

    if not name:
        return "protein_family_count_output"

    invalid_chars = '<>:"/\\|?*'

    for char in invalid_chars:
        name = name.replace(char, "_")

    for ext in [".txt", ".tsv", ".csv"]:
        if name.lower().endswith(ext):
            name = name[: -len(ext)]

    return name


def clean_model_output(text):
    if not text:
        return ""

    replacements = {
        "<br>": "; ",
        "<br/>": "; ",
        "<br />": "; ",
        "<BR>": "; ",
        "<BR/>": "; ",
        "<BR />": "; ",
    }

    cleaned = text

    for old, new in replacements.items():
        cleaned = cleaned.replace(old, new)

    return cleaned


config = read_config(CONFIG_TXT)


# ============================================================
# OpenAI API key setup
# Priority:
# 1. User-entered API key from the website session
# 2. Streamlit Secrets
# 3. Streamlit Secrets, optional fallback
# 4. OPENAI_API_KEY environment variable, optional fallback
# ============================================================
def get_secret_or_config(key, default=""):
    # Website-provided settings are preferred and are stored only in the current Streamlit session.
    try:
        if key == "KEY":
            value = st.session_state.get("user_openai_api_key", "")
            if value:
                return str(value).strip()

        if key == "DEFAULT_MODEL":
            value = st.session_state.get("user_selected_model", "")
            if value:
                return str(value).strip()

        if key == "TEMPERATURE":
            value = st.session_state.get("user_temperature", "")
            if value != "":
                return str(value).strip()

        if key == "MAX_TOKENS":
            value = st.session_state.get("user_max_tokens", "")
            if value != "":
                return str(value).strip()
    except Exception:
        pass

    try:
        value = st.secrets.get(key, default)
    except Exception:
        value = default

    if value in [None, ""]:
        value = config.get(key, default)

    return str(value).strip()


API_KEY = get_secret_or_config("KEY", "")
DEFAULT_MODEL = get_secret_or_config("DEFAULT_MODEL", "gpt-4o-mini")

try:
    TEMPERATURE = float(get_secret_or_config("TEMPERATURE", "0.5"))
except ValueError:
    TEMPERATURE = 0.5

try:
    MAX_TOKENS = int(get_secret_or_config("MAX_TOKENS", "2000"))
except ValueError:
    MAX_TOKENS = 2000

if not API_KEY:
    API_KEY = os.getenv("OPENAI_API_KEY", "").strip()

if not API_KEY:
    st.error("OpenAI API key not found. Enter your API key on the main GenoAnno page before opening a tool.")
    st.stop()

client = OpenAI(api_key=API_KEY)


# ============================================================
# Exact protein-family interpretation prompt
# ============================================================
EXACT_PROMPT_TEMPLATE = """Prompt for protein-family interpretation

1. AI role

You are an expert in oral microbiology, bacterial proteins, and microbial physiology. Your task is to interpret protein family information from a bacterial genome and generate evidence-based candidate phenotypes.

2. Input description

The bacterial genome (an oral bacterium) has been annotated, and genes have been assigned to protein function families. The input table contains two columns:

Column 1: Protein function family — the name of the protein function family.

Column 2: Gene count — the number of genes within the protein function family.

The table below has already been processed, and generic or uninformative protein families were removed.

{processed table text}

3. Analysis instructions

Using only the protein function families provided in the input table:

1. Group related protein function families into biologically meaningful candidate phenotype categories.
2. Explain the biological relationship between the associated protein function families and the candidate phenotype.
3. Give more weight to protein function families with higher gene counts, but do not overclaim phenotype prediction from count alone.
4. Do not use protein function families that are absent from the input table as evidence.
5. Do not introduce phenotypes that cannot be reasonably supported by the provided protein function family information.
6. Do not force every protein function family into a phenotype category. If a protein function family does not provide sufficient information to support a meaningful phenotype, it may be left ungrouped.
7. Prefer biologically specific phenotype categories over overly broad categories such as “metabolism” or “energy production.”

4. Reporting instructions

Generate a table with the following three columns:

Column 1: Candidate phenotype — Provide a concise description of a potential bacterial characteristic.

Column 2: Associated protein families — List the protein families supporting the phenotype. If multiple families belong to the same phenotype, separate them using semicolons.

Column 3: Evidence-based explanation — Explain how the associated families collectively support the candidate phenotype. Do not provide evidence beyond the information available from the input protein families.
"""


# ============================================================
# Helper functions
# ============================================================
def guess_protein_family_column(df):
    columns = list(df.columns)

    priority_terms = [
        "protein function family",
        "protein_family",
        "protein family",
        "pfam",
        "family",
        "function",
        "product",
        "description",
        "protein",
    ]

    for term in priority_terms:
        for col in columns:
            if term.lower() in col.lower():
                return col

    return columns[0]


def split_protein_family_terms(value):
    if pd.isna(value):
        return []

    text = str(value).strip()

    if not text:
        return []

    separators = [";", "|"]

    terms = [text]

    for sep in separators:
        new_terms = []
        for item in terms:
            new_terms.extend(item.split(sep))
        terms = new_terms

    cleaned_terms = []

    for term in terms:
        clean_term = term.strip()

        if not clean_term:
            continue

        # The prompt tells the model that "generic or uninformative protein
        # families were removed", so the uninformative placeholder annotations
        # have to be dropped here for that statement to hold. Otherwise
        # placeholder rows compete for the top-N slots with real families and
        # the model reports them as a phenotype category.
        if clean_term.lower() in [
            "nan",
            "na",
            "n/a",
            "none",
            "null",
            "-",
            "false",
            "0",
            "uncharacterized protein",
            "hypothetical protein",
            "unknown",
            "unknown function",
            "putative protein",
            "predicted protein",
        ]:
            continue

        cleaned_terms.append(clean_term)

    return cleaned_terms


def build_protein_family_count_table(df, protein_family_column):
    all_terms = []

    for value in df[protein_family_column]:
        terms = split_protein_family_terms(value)

        # Count each family at most once per gene. The source column lists
        # every domain hit on a protein, so a protein carrying the same
        # family twice (repeat domains) would otherwise contribute 2 to a
        # column the prompt labels "Gene count". De-duplicating per row keeps
        # the count a gene count rather than a domain-occurrence count.
        # dict.fromkeys preserves first-seen order, so output stays stable.
        all_terms.extend(dict.fromkeys(terms))

    if not all_terms:
        return pd.DataFrame(columns=["Protein function family", "Gene count"])

    count_df = pd.Series(all_terms).value_counts().reset_index()
    count_df.columns = ["Protein function family", "Gene count"]

    count_df = count_df.sort_values(
        by=["Gene count", "Protein function family"],
        ascending=[False, True],
    ).reset_index(drop=True)

    count_df.index = count_df.index + 1

    return count_df


def get_top_n_with_ties(count_df, n=10):
    """
    Returns the top n protein families, including all protein families tied with the nth row.

    If the 10th protein family has Gene count = 7, every protein family with Gene count >= 7
    will be included. Therefore, the final table can contain 10, 11, 12, etc. rows.
    """
    if count_df.empty:
        return pd.DataFrame(columns=["Protein function family", "Gene count"])

    sorted_df = count_df.sort_values(
        by=["Gene count", "Protein function family"],
        ascending=[False, True],
    ).reset_index(drop=True)

    if len(sorted_df) <= n:
        top_df = sorted_df.copy()
    else:
        cutoff_value = sorted_df.iloc[n - 1]["Gene count"]
        top_df = sorted_df[sorted_df["Gene count"] >= cutoff_value].copy()

    top_df = top_df.reset_index(drop=True)
    top_df.index = top_df.index + 1

    return top_df


def dataframe_to_tsv_text(df):
    if df.empty:
        return "No protein function families remained after processing."

    return df.to_csv(sep="\t", index=False)


def build_prompt(processed_table_text):
    return EXACT_PROMPT_TEMPLATE.replace(
        "{processed table text}",
        processed_table_text.strip(),
    )


def generate_output(final_prompt):
    response = client.chat.completions.create(
        model=DEFAULT_MODEL,
        messages=[
            {
                "role": "user",
                "content": final_prompt,
            },
        ],
        temperature=TEMPERATURE,
        max_tokens=MAX_TOKENS,
    )

    return response.choices[0].message.content


# ============================================================
# Streamlit website
# ============================================================
st.set_page_config(
    page_title="GenoAnno Protein Family Count App",
    layout="wide",
    initial_sidebar_state="collapsed",
)

st.markdown(
    """
    <style>
    .stApp {
        background: #f3f6fb;
        color: #111827;
    }

    header[data-testid="stHeader"] {
        background: rgba(243, 246, 251, 0.95);
    }

    #MainMenu {
        visibility: hidden;
    }

    footer {
        visibility: hidden;
    }

    .block-container {
        max-width: 96%;
        padding-top: 1.5rem;
        padding-left: 2.5rem;
        padding-right: 2.5rem;
        padding-bottom: 4rem;
    }

    .main-header {
        width: 100%;
        background: #111827;
        color: white;
        padding: 2.2rem 2.6rem;
        border-radius: 18px;
        margin-bottom: 2rem;
        box-shadow: 0 18px 40px rgba(17, 24, 39, 0.18);
    }

    .main-header h1 {
        font-size: 2.6rem;
        font-weight: 800;
        margin: 0;
        letter-spacing: -0.03em;
    }

    .main-header p {
        font-size: 1.15rem;
        margin-top: 0.8rem;
        margin-bottom: 0;
        color: #d1d5db;
        line-height: 1.65;
        max-width: 1250px;
    }

    .section-title {
        font-size: 1.35rem;
        font-weight: 800;
        color: #111827;
        margin-top: 1.5rem;
        margin-bottom: 0.3rem;
    }

    .section-note {
        font-size: 1rem;
        color: #4b5563;
        margin-bottom: 0.75rem;
        line-height: 1.5;
    }

    div[data-testid="stTextArea"] textarea {
        background-color: #ffffff;
        border: 1px solid #cfd8e3;
        border-radius: 12px;
        color: #111827;
        font-size: 1.05rem;
        line-height: 1.65;
        padding: 1rem;
        box-shadow: inset 0 1px 2px rgba(15, 23, 42, 0.04);
    }

    div[data-testid="stTextInput"] input {
        background-color: #ffffff;
        border: 1px solid #cfd8e3;
        border-radius: 12px;
        color: #111827;
        font-size: 1.05rem;
        padding: 0.85rem 1rem;
        box-shadow: inset 0 1px 2px rgba(15, 23, 42, 0.04);
    }

    div[data-testid="stFileUploader"] {
        background: #ffffff;
        border: 1px solid #d9e1ec;
        border-radius: 16px;
        padding: 1rem;
        box-shadow: 0 8px 20px rgba(15, 23, 42, 0.04);
    }

    label {
        color: #1f2937 !important;
        font-size: 1rem !important;
        font-weight: 700 !important;
    }

    .stButton > button {
        background: #111827;
        color: #ffffff;
        border: 1px solid #111827;
        border-radius: 12px;
        padding: 0.85rem 1.8rem;
        font-size: 1.05rem;
        font-weight: 800;
        box-shadow: 0 12px 26px rgba(17, 24, 39, 0.18);
        transition: all 0.15s ease-in-out;
    }

    .stButton > button:hover {
        background: #374151;
        border-color: #374151;
        color: #ffffff;
        transform: translateY(-1px);
    }

    .stExpander {
        background: #ffffff;
        border-radius: 14px;
        border: 1px solid #e5e7eb;
        box-shadow: 0 6px 18px rgba(15, 23, 42, 0.04);
        margin-bottom: 1.3rem;
    }

    .output-panel {
        background: #ffffff;
        border: 1px solid #e1e7ef;
        border-radius: 18px;
        padding: 1.7rem 1.9rem;
        margin-top: 1.6rem;
        box-shadow: 0 12px 32px rgba(15, 23, 42, 0.08);
    }

    .output-panel h2 {
        font-size: 2rem;
        font-weight: 800;
        color: #111827;
        margin-top: 0;
        margin-bottom: 1rem;
    }

    .save-panel {
        background: #ffffff;
        border: 1px solid #e1e7ef;
        border-radius: 18px;
        padding: 1.5rem 1.7rem;
        margin-top: 1.4rem;
        box-shadow: 0 12px 32px rgba(15, 23, 42, 0.06);
    }

    .save-panel h3 {
        font-size: 1.55rem;
        font-weight: 800;
        color: #111827;
        margin-top: 0;
        margin-bottom: 0.4rem;
    }

    .small-muted {
        color: #4b5563;
        font-size: 1rem;
        margin-bottom: 1rem;
        line-height: 1.5;
    }

    div[data-testid="stMarkdownContainer"] {
        font-size: 1.05rem;
        line-height: 1.65;
    }

    div[data-testid="stAlert"] {
        border-radius: 14px;
        font-size: 1rem;
    }

    .divider {
        height: 1px;
        background: #dce3ed;
        margin: 1.8rem 0;
    }

    .metric-card {
        background: #ffffff;
        border: 1px solid #e1e7ef;
        border-radius: 16px;
        padding: 1.1rem 1.3rem;
        box-shadow: 0 8px 22px rgba(15, 23, 42, 0.05);
    }

    .metric-card .label {
        color: #6b7280;
        font-size: 0.95rem;
        font-weight: 700;
        margin-bottom: 0.3rem;
    }

    .metric-card .value {
        color: #111827;
        font-size: 1.7rem;
        font-weight: 850;
    }

    .stDownloadButton > button {
        background: #111827;
        color: #ffffff;
        border: 1px solid #111827;
        border-radius: 12px;
        padding: 0.75rem 1.35rem;
        font-size: 0.98rem;
        font-weight: 800;
        box-shadow: 0 10px 22px rgba(17, 24, 39, 0.14);
        transition: all 0.15s ease-in-out;
        margin-top: 0.55rem;
        margin-bottom: 0.8rem;
    }

    .stDownloadButton > button:hover {
        background: #374151;
        border-color: #374151;
        color: #ffffff;
        transform: translateY(-1px);
    }

    </style>
    """,
    unsafe_allow_html=True,
)

st.markdown(
    """
    <div class="main-header">
        <h1>Protein Family Interpreter</h1>
    </div>
    """,
    unsafe_allow_html=True,
)


# ============================================================
# Fixed prompt
# ============================================================
st.markdown('<div class="section-title">Prompt</div>', unsafe_allow_html=True)
st.markdown(
    '<div class="section-note">The protein-family interpretation prompt is fixed exactly as specified. The processed protein-family count table is inserted automatically.</div>',
    unsafe_allow_html=True,
)


# ============================================================
# 2. Input
# ============================================================
st.markdown('<div class="section-title">2. Input</div>', unsafe_allow_html=True)
st.markdown(
    '<div class="section-note">Upload a TSV annotation file and select the column containing protein function family terms.</div>',
    unsafe_allow_html=True,
)

uploaded_file = st.file_uploader(
    "Upload protein family annotation TSV file",
    type=["tsv", "txt", "csv"],
)

count_df = None
top10_df = None
raw_df = None
processed_table_text = None
top10_table_text = None
include_top10 = True

if uploaded_file is not None:
    try:
        file_name = uploaded_file.name.lower()

        if file_name.endswith(".csv"):
            raw_df = pd.read_csv(uploaded_file)
        else:
            raw_df = pd.read_csv(uploaded_file, sep="\t")

        st.markdown('<div class="section-title">2a. Original uploaded input</div>', unsafe_allow_html=True)
        st.markdown(
            '<div class="section-note">This is the raw uploaded annotation table before protein family counting.</div>',
            unsafe_allow_html=True,
        )

        st.dataframe(
            raw_df,
            use_container_width=True,
            height=300,
        )

        guessed_column = guess_protein_family_column(raw_df)
        column_options = list(raw_df.columns)
        default_index = column_options.index(guessed_column)

        protein_family_column = st.selectbox(
            "Select column containing protein function family terms",
            options=column_options,
            index=default_index,
        )

        count_df = build_protein_family_count_table(
            df=raw_df,
            protein_family_column=protein_family_column,
        )

        total_rows = len(raw_df)
        unique_families = len(count_df)
        total_mapped_terms = int(count_df["Gene count"].sum()) if not count_df.empty else 0

        metric_col1, metric_col2, metric_col3 = st.columns(3)

        with metric_col1:
            st.markdown(
                f"""
                <div class="metric-card">
                    <div class="label">Uploaded rows</div>
                    <div class="value">{total_rows}</div>
                </div>
                """,
                unsafe_allow_html=True,
            )

        with metric_col2:
            st.markdown(
                f"""
                <div class="metric-card">
                    <div class="label">Unique protein families</div>
                    <div class="value">{unique_families}</div>
                </div>
                """,
                unsafe_allow_html=True,
            )

        with metric_col3:
            st.markdown(
                f"""
                <div class="metric-card">
                    <div class="label">Total mapped terms</div>
                    <div class="value">{total_mapped_terms}</div>
                </div>
                """,
                unsafe_allow_html=True,
            )

        st.markdown('<div class="section-title">2b. Protein family count table</div>', unsafe_allow_html=True)
        st.markdown(
            '<div class="section-note">The table below groups repeated protein function family terms and counts how many times each term appears.</div>',
            unsafe_allow_html=True,
        )

        st.dataframe(
            count_df,
            use_container_width=True,
            height=420,
        )

        include_top10 = st.checkbox(
            "Generate top 10 pivot/count table",
            value=True,
        )

        if include_top10:
            top10_df = get_top_n_with_ties(count_df, n=10)

            st.markdown(
                '<div class="section-title">2c. Top 10 protein family pivot/count table</div>',
                unsafe_allow_html=True,
            )

            st.markdown(
                '<div class="section-note">This table includes the top 10 gene-count level and all protein families tied at the 10th position. Therefore, it may contain more than 10 protein families when counts are tied.</div>',
                unsafe_allow_html=True,
            )

            st.dataframe(
                top10_df,
                use_container_width=True,
                height=380,
            )

            top10_download_text = top10_df.to_csv(
                sep="\t",
                index=False,
            )

            st.download_button(
                label="Download Top 10 protein family table",
                data=top10_download_text,
                file_name="protein_family_top10_with_ties_table.tsv",
                mime="text/tab-separated-values",
                use_container_width=False,
            )
        else:
            top10_df = pd.DataFrame(columns=["Protein function family", "Gene count"])

        processed_table_text = dataframe_to_tsv_text(count_df)
        top10_table_text = dataframe_to_tsv_text(top10_df)

    except Exception as upload_error:
        st.error(f"Could not process uploaded file:\n\n{upload_error}")

else:
    st.info("Upload a protein family annotation TSV file to begin analysis.")


# ============================================================
# Final prompt preview and generation
# ============================================================
if processed_table_text is not None:
    final_prompt = build_prompt(
        processed_table_text=processed_table_text,
    )

    with st.expander("Preview final prompt sent to model", expanded=False):
        st.text_area(
            "Final prompt",
            value=final_prompt,
            height=420,
            label_visibility="collapsed",
        )

    st.markdown('<div class="divider"></div>', unsafe_allow_html=True)

    generate_button = st.button("Generate phenotype grouping")

    if generate_button:
        if count_df is None or count_df.empty:
            st.warning("No protein function family terms were found. Please check the selected column.")
        else:
            with st.spinner("Generating phenotype grouping..."):
                try:
                    output = generate_output(
                        final_prompt=final_prompt,
                    )

                    output = clean_model_output(output)

                    st.session_state["output"] = output

                except Exception as error:
                    st.error(f"Generation failed:\n\n{error}")


# ============================================================
# Output
# ============================================================
if "output" in st.session_state:
    output = st.session_state["output"]

    st.markdown('<div class="output-panel">', unsafe_allow_html=True)
    st.markdown("<h2>Output</h2>", unsafe_allow_html=True)
    st.markdown(output)
    st.markdown("</div>", unsafe_allow_html=True)

    st.download_button(
        label="Download AI interpretation",
        data=output,
        file_name="protein_family_interpretation.txt",
        mime="text/plain",
        use_container_width=False,
    )
