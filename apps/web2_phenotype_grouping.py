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
        return "phenotype_grouping_output"

    invalid_chars = '<>:"/\\|?*'

    for char in invalid_chars:
        name = name.replace(char, "_")

    if name.lower().endswith(".txt"):
        name = name[:-4]

    return name


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
# Default editable chunks
# ============================================================
DEFAULT_AI_ROLE = """A Veteran Professor in Biological Science"""

DEFAULT_ANALYSIS = """Using the processed input table, group the retained metabolic pathways into biologically meaningful phenotypes. Interpret only pathways that remain after filtering."""

DEFAULT_REPORTING = """Generate an output table with the following three columns:

1. Phenotype
2. Associated pathways
3. Explanation describing how the listed pathways support or define the phenotype

Finally, provide a brief overall summary describing the predicted phenotype(s) of the bacterium based on these pathway groupings."""

DEFAULT_INPUT_DESCRIPTION = """I annotated genes from a bacterial genome and assigned each gene to its corresponding metabolic pathway. For each pathway, I calculated its completeness and, where applicable, determined whether the pathway is present (TRUE/FALSE).

The input table contains two columns after processing:
Column 1: Metabolic pathway name
Column 2: Pathway status, either completeness score or status indicator

Filtering rule:
Pathways with 0, FALSE, empty, missing, NA, or NaN values are removed before analysis.

Only the processed input table should be used for phenotype grouping."""


# ============================================================
# TSV processing
# ============================================================
def is_removed_value(value):
    if pd.isna(value):
        return True

    value_str = str(value).strip()

    if value_str == "":
        return True

    if value_str.lower() in [
        "false",
        "f",
        "no",
        "none",
        "nan",
        "na",
        "n/a",
        "missing",
    ]:
        return True

    try:
        numeric_value = float(value_str)
        if numeric_value == 0:
            return True
    except ValueError:
        pass

    return False


def process_products_tsv(uploaded_file, selected_index=0):
    df = pd.read_csv(uploaded_file, sep="\t")

    if df.empty:
        empty_processed = pd.DataFrame(
            columns=["Metabolic pathway name", "Pathway status"]
        )
        empty_removed = pd.DataFrame(
            columns=["Metabolic pathway name", "Removed value", "Removal reason"]
        )
        return empty_processed, empty_removed, df, "empty"

    # Wide format:
    # genome | pathway1 | pathway2 | pathway3 ...
    if "genome" in df.columns and df.shape[1] > 2:
        selected_index = min(selected_index, len(df) - 1)
        row = df.iloc[selected_index]

        retained_records = []
        removed_records = []

        for col in df.columns:
            if col == "genome":
                continue

            value = row[col]

            if is_removed_value(value):
                removed_records.append(
                    {
                        "Metabolic pathway name": col,
                        "Removed value": value,
                        "Removal reason": "0, FALSE, empty, missing, NA, or NaN",
                    }
                )
            else:
                retained_records.append(
                    {
                        "Metabolic pathway name": col,
                        "Pathway status": value,
                    }
                )

        processed_df = pd.DataFrame(retained_records)
        removed_df = pd.DataFrame(removed_records)

        return processed_df, removed_df, df, "wide"

    # Long format:
    # pathway | value
    if df.shape[1] >= 2:
        first_col = df.columns[0]
        second_col = df.columns[1]

        temp_df = df[[first_col, second_col]].copy()
        temp_df.columns = ["Metabolic pathway name", "Pathway status"]

        removed_df = temp_df[
            temp_df["Pathway status"].apply(is_removed_value)
        ].copy()

        removed_df = removed_df.rename(
            columns={"Pathway status": "Removed value"}
        )

        if not removed_df.empty:
            removed_df["Removal reason"] = "0, FALSE, empty, missing, NA, or NaN"

        processed_df = temp_df[
            ~temp_df["Pathway status"].apply(is_removed_value)
        ].copy()

        processed_df = processed_df.reset_index(drop=True)
        removed_df = removed_df.reset_index(drop=True)

        return processed_df, removed_df, df, "long"

    empty_processed = pd.DataFrame(
        columns=["Metabolic pathway name", "Pathway status"]
    )
    empty_removed = pd.DataFrame(
        columns=["Metabolic pathway name", "Removed value", "Removal reason"]
    )

    return empty_processed, empty_removed, df, "unsupported"


def dataframe_to_tsv_text(df):
    if df.empty:
        return "No pathways remained after filtering."

    return df.to_csv(sep="\t", index=False)


# ============================================================
# Build final prompt
# ============================================================
def build_prompt(ai_role, input_description, processed_table_text, analysis, reporting):
    final_prompt = f"""
1. AI role:
{ai_role}

2. Input:
{input_description}

Processed input table after removing pathways with 0, FALSE, empty, missing, NA, or NaN values:
{processed_table_text}

3. Analysis:
{analysis}

4. Reporting instructions:
{reporting}
""".strip()

    return final_prompt


# ============================================================
# OpenAI call
# ============================================================
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
    page_title="GenoAnno Phenotype Grouping App",
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
    </style>
    """,
    unsafe_allow_html=True,
)

st.markdown(
    """
    <div class="main-header">
        <h1>GenoAnno Phenotype Grouping App</h1>
        <p>
            Upload a products.tsv file, inspect the original input, remove inactive pathways,
            review removed and retained pathways, then generate phenotype groupings using the processed input.
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
    '<div class="section-note">Define the biological expertise and reasoning style for phenotype grouping.</div>',
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
    '<div class="section-note">Upload products.tsv and define how the processed input should be interpreted.</div>',
    unsafe_allow_html=True,
)

input_description = st.text_area(
    "Input description",
    value=DEFAULT_INPUT_DESCRIPTION,
    height=230,
    label_visibility="collapsed",
)

uploaded_tsv = st.file_uploader(
    "Upload products.tsv",
    type=["tsv", "txt"],
)

processed_df = None
removed_df = None
raw_df = None
file_format = None
processed_table_text = None

if uploaded_tsv is not None:
    try:
        raw_df_preview = pd.read_csv(uploaded_tsv, sep="\t")
        uploaded_tsv.seek(0)

        selected_index = 0

        if "genome" in raw_df_preview.columns and raw_df_preview.shape[1] > 2:
            genome_options = raw_df_preview["genome"].astype(str).tolist()

            selected_index = st.selectbox(
                "Select genome row for analysis",
                options=list(range(len(genome_options))),
                format_func=lambda i: genome_options[i],
            )

        uploaded_tsv.seek(0)

        processed_df, removed_df, raw_df, file_format = process_products_tsv(
            uploaded_file=uploaded_tsv,
            selected_index=selected_index,
        )

        if file_format == "wide":
            total_raw_pathways = max(raw_df.shape[1] - 1, 0)
        elif file_format == "long":
            total_raw_pathways = len(raw_df)
        else:
            total_raw_pathways = 0

        retained_pathways = len(processed_df)
        removed_pathways = len(removed_df)

        metric_col1, metric_col2, metric_col3 = st.columns(3)

        with metric_col1:
            st.markdown(
                f"""
                <div class="metric-card">
                    <div class="label">Raw pathways</div>
                    <div class="value">{total_raw_pathways}</div>
                </div>
                """,
                unsafe_allow_html=True,
            )

        with metric_col2:
            st.markdown(
                f"""
                <div class="metric-card">
                    <div class="label">Retained after filtering</div>
                    <div class="value">{retained_pathways}</div>
                </div>
                """,
                unsafe_allow_html=True,
            )

        with metric_col3:
            st.markdown(
                f"""
                <div class="metric-card">
                    <div class="label">Removed FALSE / 0 / empty</div>
                    <div class="value">{removed_pathways}</div>
                </div>
                """,
                unsafe_allow_html=True,
            )

        st.markdown('<div class="section-title">2a. Original uploaded input</div>', unsafe_allow_html=True)
        st.markdown(
            '<div class="section-note">This is the raw uploaded TSV before filtering.</div>',
            unsafe_allow_html=True,
        )

        st.dataframe(
            raw_df,
            use_container_width=True,
            height=300,
        )

        st.markdown('<div class="section-title">2b. Removed pathways</div>', unsafe_allow_html=True)
        st.markdown(
            '<div class="section-note">These pathways were removed because their values were 0, FALSE, empty, missing, NA, or NaN.</div>',
            unsafe_allow_html=True,
        )

        if removed_df.empty:
            st.success("No pathways were removed.")
        else:
            st.dataframe(
                removed_df,
                use_container_width=True,
                height=300,
            )

        st.markdown('<div class="section-title">2c. Processed input used for AI analysis</div>', unsafe_allow_html=True)
        st.markdown(
            '<div class="section-note">Only this processed table will be sent to the model.</div>',
            unsafe_allow_html=True,
        )

        st.dataframe(
            processed_df,
            use_container_width=True,
            height=420,
        )

        processed_table_text = dataframe_to_tsv_text(processed_df)

    except Exception as upload_error:
        st.error(f"Could not process uploaded TSV file:\n\n{upload_error}")

else:
    st.info("Upload a products.tsv file to begin analysis.")


# ============================================================
# 3. Analysis
# ============================================================
st.markdown('<div class="section-title">3. Analysis</div>', unsafe_allow_html=True)
st.markdown(
    '<div class="section-note">Define the biological question for the processed input.</div>',
    unsafe_allow_html=True,
)

analysis = st.text_area(
    "Analysis",
    value=DEFAULT_ANALYSIS,
    height=140,
    label_visibility="collapsed",
)


# ============================================================
# 4. Reporting instructions
# ============================================================
st.markdown('<div class="section-title">4. Reporting instructions</div>', unsafe_allow_html=True)
st.markdown(
    '<div class="section-note">Define the expected output structure.</div>',
    unsafe_allow_html=True,
)

reporting = st.text_area(
    "Reporting instructions",
    value=DEFAULT_REPORTING,
    height=210,
    label_visibility="collapsed",
)


# ============================================================
# Final prompt preview and generation
# ============================================================
if processed_table_text is not None:
    final_prompt = build_prompt(
        ai_role=ai_role.strip(),
        input_description=input_description.strip(),
        processed_table_text=processed_table_text,
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

    generate_button = st.button("Generate phenotype grouping")

    if generate_button:
        if not ai_role.strip():
            st.warning("Please enter the AI role.")
        elif not input_description.strip():
            st.warning("Please enter the input description.")
        elif processed_df is None or processed_df.empty:
            st.warning("No pathways remain after filtering. Please check your uploaded TSV file.")
        elif not analysis.strip():
            st.warning("Please enter the analysis instructions.")
        elif not reporting.strip():
            st.warning("Please enter the reporting instructions.")
        else:
            with st.spinner("Generating phenotype grouping..."):
                try:
                    output = generate_output(
                        final_prompt=final_prompt,
                        ai_role=ai_role.strip(),
                    )

                    st.session_state["output"] = output

                except Exception as error:
                    st.error(f"Generation failed:\n\n{error}")


# ============================================================
# Output and save
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
            value="phenotype_grouping_output",
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