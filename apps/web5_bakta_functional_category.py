import os
from io import StringIO
from pathlib import Path

import pandas as pd
import streamlit as st
from openai import OpenAI


# ============================================================
# Read config.txt
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
        return "bakta_functional_category_output"

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
# Exact gene-function interpretation prompt
# ============================================================
EXACT_PROMPT_TEMPLATE = """Prompt for gene function interpretation

1. AI role

You are an expert in oral microbiology, bacterial genes, and microbial physiology. Your task is to interpret gene function information from a bacterial genome and generate evidence-based candidate phenotypes.

2. Input description

The bacterial genome (an oral bacterium) has been annotated, and molecular function has been described for each gene. The input table contains two columns:

Column 1: Functional category — the name of the gene function.

Column 2: Gene count — the number of genes within the functional category.

The table below has already been processed, and generic or uninformative categories were removed. The table contains top 40 most abundant categories

{processed table text}

3. Analysis instructions

Using only the gene functions provided in the input table:

1. Group related functional categories into biologically meaningful candidate phenotype categories.
2. Explain the biological relationship between the associated functional categories and the candidate phenotype.
3. Give more weight to functional categories with higher gene counts, but do not overclaim phenotype prediction from count alone.
4. Do not use functional categories that are absent from the input table as evidence.
5. Do not introduce phenotypes that cannot be reasonably supported by the provided gene function information.
6. Do not force every gene function into a phenotype category. If a gene function does not provide sufficient information to support a meaningful phenotype, it may be left ungrouped.
7. Prefer biologically specific phenotype categories over overly broad categories such as “metabolism” or “energy production.”

4. Reporting instructions

Generate a table with the following three columns:

Column 1: Candidate phenotype — Provide a concise description of a potential bacterial characteristic.

Column 2: Associated functional categories — List the functional categories supporting the phenotype. If multiple categories belong to the same phenotype, separate them using semicolons.

Column 3: Evidence-based explanation — Explain how the associated categories collectively support the candidate phenotype. Do not provide evidence beyond the information available from the input functional categories.
"""


# ============================================================
# Bakta / table parsing
# ============================================================
def read_uploaded_annotation_file(uploaded_file):
    file_name = uploaded_file.name.lower()
    raw_text = uploaded_file.getvalue().decode("utf-8", errors="replace")
    lines = raw_text.splitlines()

    # Bakta TSV usually has metadata lines beginning with "#"
    # and the real header starts with "#Sequence Id".
    header_index = None

    for i, line in enumerate(lines):
        if line.startswith("#Sequence Id"):
            header_index = i
            break

    if header_index is not None:
        header_line = lines[header_index].lstrip("#")
        data_lines = lines[header_index + 1 :]
        table_text = "\n".join([header_line] + data_lines)
        df = pd.read_csv(StringIO(table_text), sep="\t")
        return df, "bakta"

    # Fallback for normal CSV/TSV.
    if file_name.endswith(".csv"):
        df = pd.read_csv(StringIO(raw_text))
        return df, "csv"

    df = pd.read_csv(StringIO(raw_text), sep="\t")
    return df, "tsv"


def guess_functional_column(df):
    columns = list(df.columns)

    priority_terms = [
        "Product",
        "product",
        "functional category",
        "function",
        "description",
        "gene",
    ]

    for term in priority_terms:
        for col in columns:
            if term.lower() == col.lower():
                return col

    for term in priority_terms:
        for col in columns:
            if term.lower() in col.lower():
                return col

    return columns[0]


def is_removed_function(value, remove_hypothetical=True):
    if pd.isna(value):
        return True

    text = str(value).strip()

    if not text:
        return True

    bad_values = [
        "nan",
        "na",
        "n/a",
        "none",
        "null",
        "-",
        "false",
        "0",
        "unknown",
        "uncharacterized protein",
    ]

    if text.lower() in bad_values:
        return True

    if remove_hypothetical and text.lower() == "hypothetical protein":
        return True

    return False


def split_function_terms(value):
    if pd.isna(value):
        return []

    text = str(value).strip()

    if not text:
        return []

    # For Bakta Product column, usually each row has one product.
    # Semicolon splitting is kept optional for compound annotations.
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

        if clean_term:
            cleaned_terms.append(clean_term)

    return cleaned_terms


def split_cds_and_non_cds(df):
    """Separate protein-coding rows from everything else in a Bakta table.

    A Bakta TSV lists every annotated feature, not just proteins: tRNA, rRNA,
    ncRNA, ncRNA-region, ncRNA riboswitches, tmRNA and CRISPR arrays all carry
    a Product string and are counted as "genes" if they are not excluded.
    Those products are not protein functional categories, and because some
    ncRNA families are multi-copy they can outrank real protein products and
    reach the top-N table handed to the model.

    Returns (cds_df, non_cds_df). If the table has no "Type" column -- e.g. a
    pre-made two-column count table, or a non-Bakta source -- every row is
    treated as coding so behaviour is unchanged.
    """
    if "Type" not in df.columns:
        return df, df.iloc[0:0]

    is_cds = df["Type"].astype(str).str.strip().str.lower() == "cds"

    return df[is_cds], df[~is_cds]


def build_function_count_table(df, selected_column, remove_hypothetical=True):
    retained_terms = []
    removed_records = []

    # Drop non-coding features before counting, and log them in the removal
    # table so the row arithmetic shown in the UI still reconciles.
    df, non_cds_df = split_cds_and_non_cds(df)

    for _, non_cds_row in non_cds_df.iterrows():
        removed_records.append(
            {
                "Removed functional category": non_cds_row[selected_column],
                "Removal reason": "non-coding feature ({})".format(
                    str(non_cds_row["Type"]).strip()
                ),
            }
        )

    for value in df[selected_column]:
        if is_removed_function(value, remove_hypothetical=remove_hypothetical):
            removed_records.append(
                {
                    "Removed functional category": value,
                    "Removal reason": "empty / unknown / hypothetical / invalid",
                }
            )
            continue

        terms = split_function_terms(value)

        for term in terms:
            if not is_removed_function(term, remove_hypothetical=remove_hypothetical):
                retained_terms.append(term)

    if retained_terms:
        count_df = pd.Series(retained_terms).value_counts().reset_index()
        count_df.columns = ["Functional category", "Gene count"]

        count_df = count_df.sort_values(
            by=["Gene count", "Functional category"],
            ascending=[False, True],
        ).reset_index(drop=True)

        count_df.index = count_df.index + 1
    else:
        count_df = pd.DataFrame(columns=["Functional category", "Gene count"])

    removed_df = pd.DataFrame(removed_records)

    return count_df, removed_df


def get_top_n_with_ties(count_df, n=10):
    if count_df.empty:
        return pd.DataFrame(columns=["Functional category", "Gene count"])

    sorted_df = count_df.sort_values(
        by=["Gene count", "Functional category"],
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



def get_top_40(count_df):
    if count_df.empty:
        return pd.DataFrame(columns=["Functional category", "Gene count"])

    top40_df = count_df.sort_values(
        by=["Gene count", "Functional category"],
        ascending=[False, True],
    ).head(40).reset_index(drop=True)

    top40_df.index = top40_df.index + 1

    return top40_df


def dataframe_to_tsv_text(df):
    if df.empty:
        return "No functional categories remained after processing."

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
    page_title="Functional Interpreter",
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
    }

    div[data-testid="stTextInput"] input {
        background-color: #ffffff;
        border: 1px solid #cfd8e3;
        border-radius: 12px;
        color: #111827;
        font-size: 1.05rem;
        padding: 0.85rem 1rem;
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
    }

    .stButton > button:hover {
        background: #374151;
        border-color: #374151;
        color: #ffffff;
    }

    .output-panel {
        background: #ffffff;
        border: 1px solid #e1e7ef;
        border-radius: 18px;
        padding: 1.7rem 1.9rem;
        margin-top: 1.6rem;
        box-shadow: 0 12px 32px rgba(15, 23, 42, 0.08);
    }

    .save-panel {
        background: #ffffff;
        border: 1px solid #e1e7ef;
        border-radius: 18px;
        padding: 1.5rem 1.7rem;
        margin-top: 1.4rem;
    }

    .small-muted {
        color: #4b5563;
        font-size: 1rem;
        margin-bottom: 1rem;
        line-height: 1.5;
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
        <h1>Functional Interpreter</h1>
    </div>
    """,
    unsafe_allow_html=True,
)


# ============================================================
# Fixed prompt
# ============================================================
st.markdown('<div class="section-title">Prompt</div>', unsafe_allow_html=True)
st.markdown(
    '<div class="section-note">The gene-function interpretation prompt is fixed exactly as specified. The top 40 most abundant processed functional categories are inserted automatically.</div>',
    unsafe_allow_html=True,
)


# ============================================================
# 2. Input
# ============================================================
st.markdown('<div class="section-title">2. Input</div>', unsafe_allow_html=True)
st.markdown(
    '<div class="section-note"> </div>',
    unsafe_allow_html=True,
)

uploaded_file = st.file_uploader(
    "Upload Bakta TSV file",
    type=["tsv", "txt", "csv"],
)

count_df = None
top40_df = None
raw_df = None
removed_df = None
processed_table_text = None

if uploaded_file is not None:
    try:
        raw_df, file_format = read_uploaded_annotation_file(uploaded_file)

        st.markdown('<div class="section-title">2a. Original uploaded input</div>', unsafe_allow_html=True)
        st.markdown(
            '<div class="section-note">This is the parsed Bakta table before functional category counting.</div>',
            unsafe_allow_html=True,
        )

        st.dataframe(
            raw_df,
            use_container_width=True,
            height=300,
        )

        guessed_column = guess_functional_column(raw_df)
        column_options = list(raw_df.columns)
        default_index = column_options.index(guessed_column)

        selected_column = st.selectbox(
            "Select column containing functional category terms",
            options=column_options,
            index=default_index,
        )

        remove_hypothetical = st.checkbox(
            "Remove hypothetical protein / unknown annotations",
            value=True,
        )

        count_df, removed_df = build_function_count_table(
            df=raw_df,
            selected_column=selected_column,
            remove_hypothetical=remove_hypothetical,
        )

        total_rows = len(raw_df)
        unique_categories = len(count_df)
        total_mapped_terms = int(count_df["Gene count"].sum()) if not count_df.empty else 0
        removed_terms = len(removed_df)

        metric_col1, metric_col2, metric_col3, metric_col4 = st.columns(4)

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
                    <div class="label">Unique retained categories</div>
                    <div class="value">{unique_categories}</div>
                </div>
                """,
                unsafe_allow_html=True,
            )

        with metric_col3:
            st.markdown(
                f"""
                <div class="metric-card">
                    <div class="label">Mapped retained terms</div>
                    <div class="value">{total_mapped_terms}</div>
                </div>
                """,
                unsafe_allow_html=True,
            )

        with metric_col4:
            st.markdown(
                f"""
                <div class="metric-card">
                    <div class="label">Removed rows</div>
                    <div class="value">{removed_terms}</div>
                </div>
                """,
                unsafe_allow_html=True,
            )

        st.markdown('<div class="section-title">2b. Removed annotations</div>', unsafe_allow_html=True)
        st.markdown(
            '<div class="section-note">Rows removed because they were empty, unknown, hypothetical, or invalid.</div>',
            unsafe_allow_html=True,
        )

        if removed_df.empty:
            st.success("No annotations were removed.")
        else:
            st.dataframe(
                removed_df,
                use_container_width=True,
                height=250,
            )

        st.markdown('<div class="section-title">2c. Functional category count table</div>', unsafe_allow_html=True)
        st.markdown(
            '<div class="section-note">Repeated functional categories are grouped and counted here.</div>',
            unsafe_allow_html=True,
        )

        st.dataframe(
            count_df,
            use_container_width=True,
            height=420,
        )

        top40_df = get_top_40(count_df)

        st.markdown(
            '<div class="section-title">2d. Top 40 functional categories used for AI analysis</div>',
            unsafe_allow_html=True,
        )

        st.markdown(
            '<div class="section-note">These are the 40 most abundant retained functional categories and are the only categories inserted into the interpretation prompt.</div>',
            unsafe_allow_html=True,
        )

        st.dataframe(
            top40_df,
            use_container_width=True,
            height=500,
        )

        top40_download_text = top40_df.to_csv(
            sep="\t",
            index=False,
        )

        st.download_button(
            label="Download Top 40 functional categories",
            data=top40_download_text,
            file_name="bakta_functional_category_top40_table.tsv",
            mime="text/tab-separated-values",
            use_container_width=False,
        )

        processed_table_text = dataframe_to_tsv_text(top40_df)

    except Exception as upload_error:
        st.error(f"Could not process uploaded file:\n\n{upload_error}")

else:
    st.info("Upload a Bakta TSV file to begin analysis.")


# ============================================================
# Generate
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
            st.warning("No functional category terms were found. Please check the selected column.")
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
        file_name="functional_category_interpretation.txt",
        mime="text/plain",
        use_container_width=False,
    )
