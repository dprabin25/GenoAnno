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

    if not path.exists():
        st.error(f"config.txt was not found here:\n\n{path}")
        st.stop()

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
        return "similar_oral_bacteria_output"

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
# OpenAI / Streamlit Cloud secret setup
# Priority: Streamlit Secrets > config.txt > environment variable
# ============================================================
def get_secret_or_config(key, default=""):
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
    st.error("OpenAI API key not found. Add KEY in Streamlit Secrets or set OPENAI_API_KEY.")
    st.stop()

client = OpenAI(api_key=API_KEY)


# ============================================================
# Default editable chunks
# ============================================================
DEFAULT_AI_ROLE = """A Veteran Professor in Biological Science"""

DEFAULT_INPUT_DESCRIPTION = """I annotated genes from a bacterial genome using Bakta. This is an oral bacterium. The uploaded Bakta TSV contains gene-level annotations. I selected the Product column and counted how many genes map to each repeated product/function term.

The processed input table contains two columns:
Column 1: Functional category
Column 2: Gene count

Only the processed functional category count table should be used for oral-bacteria functional similarity comparison."""

DEFAULT_ANALYSIS = """Are there other oral bacteria that have similar gene composition to my bacteria?

Using the processed functional category count table from this oral bacterium, identify well-characterized oral bacteria that may have similar gene composition or similar functional-category patterns.

Focus on shared functional categories such as:
- TonB-dependent receptors
- SusC/SusD/RagA/RagB nutrient uptake systems
- outer membrane proteins
- lipoproteins
- ABC transporters
- MFS transporters
- glycoside hydrolases
- glycosyltransferases
- peptidases/proteases
- hemolysin or virulence-related proteins
- secretion-related proteins
- stress response and redox proteins
- DNA mobility/transposase systems
- biofilm-associated functions

Important interpretation rules:
- Do not claim exact genomic similarity unless the input strongly supports it.
- Describe similarity as functional similarity or phenotype-level similarity when appropriate.
- Include only well-characterized oral bacteria.
- Be conservative and scientifically accurate.
- Do not invent shared functions that are not represented in the processed input table."""

DEFAULT_REPORTING = """Provide only the following output:

1. A markdown table with these columns:
   - Closely related or functionally similar oral bacterium
   - Shared functional categories
   - Shared phenotype
   - Reason for similarity
   - Confidence

2. A short summary paragraph explaining which bacterial phenotypes are most consistent with the uploaded functional category profile.

Formatting rules:
- Use a clean markdown table.
- Do not use HTML tags.
- Do not use <br>, <br/>, or <br />.
- Separate multiple shared functions using semicolons.
- Confidence must be High, Moderate, or Low.
- If the similarity is functional rather than taxonomic, say so clearly.
- Do not provide extra sections outside the table and summary."""


# ============================================================
# Bakta / table parsing
# ============================================================
def read_uploaded_annotation_file(uploaded_file):
    file_name = uploaded_file.name.lower()
    raw_text = uploaded_file.getvalue().decode("utf-8", errors="replace")
    lines = raw_text.splitlines()

    header_index = None

    for i, line in enumerate(lines):
        if line.startswith("#Sequence Id"):
            header_index = i
            break

    if header_index is not None:
        header_line = lines[header_index].lstrip("#")
        data_lines = lines[header_index + 1:]
        table_text = "\n".join([header_line] + data_lines)
        df = pd.read_csv(StringIO(table_text), sep="\t")
        return df, "bakta"

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


def build_function_count_table(df, selected_column, remove_hypothetical=True):
    retained_terms = []
    removed_records = []

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


def dataframe_to_tsv_text(df):
    if df.empty:
        return "No functional categories remained after processing."

    return df.to_csv(sep="\t", index=False)


def build_prompt(
    ai_role,
    input_description,
    processed_table_text,
    top10_table_text,
    include_top10,
    analysis,
    reporting,
):
    if include_top10:
        top10_section = f"""
Top 10 functional category count table, including all functional categories tied at the 10th position:
{top10_table_text}
""".strip()
    else:
        top10_section = "Top 10 functional category count table was not requested."

    final_prompt = f"""
1. AI role:
{ai_role}

2. Input:
{input_description}

Processed functional category count table:
{processed_table_text}

{top10_section}

3. Analysis:
{analysis}

4. Reporting instructions:
{reporting}
""".strip()

    return final_prompt


def generate_output(final_prompt, ai_role):
    response = client.chat.completions.create(
        model=DEFAULT_MODEL,
        messages=[
            {
                "role": "system",
                "content": ai_role,
            },
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
    page_title="GenoAnno Similar Oral Bacteria App",
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
    </style>
    """,
    unsafe_allow_html=True,
)

st.markdown(
    """
    <div class="main-header">
        <h1>GenoAnno Similar Oral Bacteria App</h1>
        <p>
            Upload a Bakta TSV file, count repeated Product annotations, then identify well-characterized oral bacteria
            with similar functional-category composition and phenotype-level similarity.
        </p>
    </div>
    """,
    unsafe_allow_html=True,
)


# ============================================================
# 1. AI role
# ============================================================
st.markdown('<div class="section-title">1. AI role</div>', unsafe_allow_html=True)
st.markdown(
    '<div class="section-note">Define the biological expertise and reasoning style.</div>',
    unsafe_allow_html=True,
)

ai_role = st.text_area(
    "AI role",
    value=DEFAULT_AI_ROLE,
    height=110,
    label_visibility="collapsed",
)


# ============================================================
# 2. Input
# ============================================================
st.markdown('<div class="section-title">2. Input</div>', unsafe_allow_html=True)
st.markdown(
    '<div class="section-note">Upload Bakta TSV and select the annotation column to count.</div>',
    unsafe_allow_html=True,
)

input_description = st.text_area(
    "Input description",
    value=DEFAULT_INPUT_DESCRIPTION,
    height=190,
    label_visibility="collapsed",
)

uploaded_file = st.file_uploader(
    "Upload Bakta TSV file",
    type=["tsv", "txt", "csv"],
)

count_df = None
top10_df = None
raw_df = None
removed_df = None
processed_table_text = None
top10_table_text = None
include_top10 = True

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

        include_top10 = st.checkbox(
            "Generate top 10 count table",
            value=True,
        )

        if include_top10:
            top10_df = get_top_n_with_ties(count_df, n=10)

            st.markdown(
                '<div class="section-title">2d. Top 10 functional category count table</div>',
                unsafe_allow_html=True,
            )

            st.markdown(
                '<div class="section-note">Includes the top 10 count level and all categories tied at the 10th position.</div>',
                unsafe_allow_html=True,
            )

            st.dataframe(
                top10_df,
                use_container_width=True,
                height=380,
            )
        else:
            top10_df = pd.DataFrame(columns=["Functional category", "Gene count"])

        processed_table_text = dataframe_to_tsv_text(count_df)
        top10_table_text = dataframe_to_tsv_text(top10_df)

    except Exception as upload_error:
        st.error(f"Could not process uploaded file:\n\n{upload_error}")

else:
    st.info("Upload a Bakta TSV file to begin analysis.")


# ============================================================
# 3. Analysis
# ============================================================
st.markdown('<div class="section-title">3. Analysis</div>', unsafe_allow_html=True)
st.markdown(
    '<div class="section-note">Define the oral-bacteria comparison question.</div>',
    unsafe_allow_html=True,
)

analysis = st.text_area(
    "Analysis",
    value=DEFAULT_ANALYSIS,
    height=260,
    label_visibility="collapsed",
)


# ============================================================
# 4. Reporting instructions
# ============================================================
st.markdown('<div class="section-title">4. Reporting instructions</div>', unsafe_allow_html=True)
st.markdown(
    '<div class="section-note">Define the output table and summary format.</div>',
    unsafe_allow_html=True,
)

reporting = st.text_area(
    "Reporting instructions",
    value=DEFAULT_REPORTING,
    height=300,
    label_visibility="collapsed",
)


# ============================================================
# Generate
# ============================================================
if processed_table_text is not None:
    final_prompt = build_prompt(
        ai_role=ai_role.strip(),
        input_description=input_description.strip(),
        processed_table_text=processed_table_text,
        top10_table_text=top10_table_text,
        include_top10=include_top10,
        analysis=analysis.strip(),
        reporting=reporting.strip(),
    )

    with st.expander("Preview final prompt sent to model", expanded=False):
        st.text_area(
            "Final prompt",
            value=final_prompt,
            height=420,
            label_visibility="collapsed",
        )

    st.markdown('<div class="divider"></div>', unsafe_allow_html=True)

    generate_button = st.button(
        "Generate similar oral bacteria comparison",
        key="generate_similar_oral_bacteria",
    )

    if generate_button:
        if not ai_role.strip():
            st.warning("Please enter the AI role.")
        elif not input_description.strip():
            st.warning("Please enter the input description.")
        elif count_df is None or count_df.empty:
            st.warning("No functional category terms were found. Please check the selected column.")
        elif not analysis.strip():
            st.warning("Please enter the analysis instructions.")
        elif not reporting.strip():
            st.warning("Please enter the reporting instructions.")
        else:
            with st.spinner("Generating similar oral bacteria comparison..."):
                try:
                    output = generate_output(
                        final_prompt=final_prompt,
                        ai_role=ai_role.strip(),
                    )

                    output = clean_model_output(output)
                    st.session_state["output"] = output

                except Exception as error:
                    st.error(f"Generation failed:\n\n{error}")


# ============================================================
# Save GPT output only
# ============================================================
if "output" in st.session_state:
    output = st.session_state["output"]

    st.markdown('<div class="output-panel">', unsafe_allow_html=True)
    st.markdown("<h2>Output</h2>", unsafe_allow_html=True)
    st.markdown(output)
    st.markdown("</div>", unsafe_allow_html=True)

    st.markdown(
        """
        <div class="save-panel">
            <h3>Save output as TXT</h3>
            <p class="small-muted">
                Choose the file name and local directory where the text output should be saved.
            </p>
        </div>
        """,
        unsafe_allow_html=True,
    )

    col1, col2 = st.columns(2, gap="large")

    with col1:
        output_file_name = st.text_input(
            "Output file name",
            value="similar_oral_bacteria_output",
        )

    with col2:
        output_directory = st.text_input(
            "Output directory",
            value=str(HERE),
        )

    safe_file_name = clean_filename(output_file_name)
    output_path = Path(output_directory) / f"{safe_file_name}.txt"

    if st.button("Save TXT to selected directory"):
        try:
            Path(output_directory).mkdir(parents=True, exist_ok=True)
            output_path.write_text(output, encoding="utf-8")

            st.success(f"Saved successfully:\n\n{output_path}")

        except Exception as save_error:
            st.error(f"Could not save file:\n\n{save_error}")