import os
from pathlib import Path

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
        return "pathway_expectation_output"

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
# Editable default chunks
# ============================================================
DEFAULT_AI_ROLE = """A Veteran Professor in Biological Science"""

DEFAULT_PATHWAYS = """3-Hydroxypropionate bi-cycle	Acetyl-CoA pathway, CO2 => acetyl-CoA	Citrate cycle (TCA cycle, Krebs cycle)	Dicarboxylate-hydroxybutyrate cycle	Entner-Doudoroff pathway, glucose-6P => glyceraldehyde-3P + pyruvate	Glycolysis (Embden-Meyerhof pathway), glucose => pyruvate	Glyoxylate cycle	Hydroxypropionate-hydroxybutylate cycle	Methanogenesis, CO2 => methane	Pentose phosphate pathway (Pentose phosphate cycle)	Reductive acetyl-CoA pathway (Wood-Ljungdahl pathway)	Reductive citrate cycle (Arnon-Buchanan cycle)	Reductive pentose phosphate cycle (Calvin cycle)	Complex I: NAD(P)H:quinone oxidoreductase, chloroplasts and cyanobacteria	Complex I: NADH dehydrogenase (ubiquinone) 1 alpha subcomplex	Complex I: NADH:quinone oxidoreductase, prokaryotes	Complex II: Fumarate reductase, prokaryotes	Complex II: Succinate dehydrogenase (ubiquinone)	Complex II: Succinate dehydrogenase, prokaryotes	Complex III: Cytochrome bc1 complex	Complex III: Cytochrome bc1 complex respiratory unit	Complex III: Cytochrome bd ubiquinol oxidase	Complex IV High affinity: Cytochrome bd ubiquinol oxidase	Complex IV High affinity: Cytochrome c oxidase, cbb3-type	Complex IV Low affinity: Cytochrome aa3-600 menaquinol oxidase	Complex IV Low affinity: Cytochrome c oxidase	Complex IV Low affinity: Cytochrome c oxidase, prokaryotes	Complex IV Low affinity: Cytochrome o ubiquinol oxidase	Complex V: F-type ATPase, eukaryotes	Complex V: F-type ATPase, prokaryotes and chloroplasts	Complex V: V-type ATPase, eukaryotes	Complex V: V/A-type ATPase, prokaryotes	CAZy: Alpha-galactans	CAZy: Alpha-mannan	CAZy: Amorphous Cellulose	CAZy: Arabinan	CAZy: Arabinose cleavage	CAZy: Beta-galactan (pectic galactan)	CAZy: Beta-mannan	CAZy: Chitin	CAZy: Crystalline Cellulose	CAZy: Fucose Cleavage	CAZy: Mixed-Linkage glucans	CAZy: Mucin	CAZy: Pectin	CAZy: Polyphenolics	CAZy: Rhamnose cleavage	CAZy: Starch	CAZy: Sulf-Polysachharides	CAZy: Xylans	CAZy: Xyloglucan	Methanogenesis and methanotrophy: Key functional gene	Methanogenesis and methanotrophy: acetate => methane, pt 1	Methanogenesis and methanotrophy: acetate => methane, pt 2	Methanogenesis and methanotrophy: acetate => methane, pt 3	Methanogenesis and methanotrophy: dimethylamine => monomethylamine	Methanogenesis and methanotrophy: methane => methanol, with oxygen (mmo)	Methanogenesis and methanotrophy: methane => methanol, with oxygen (pmo)	Methanogenesis and methanotrophy: methanol => methane	Methanogenesis and methanotrophy: monomethylamine => ammonia	Methanogenesis and methanotrophy: putative but not defining CO2 => methane	Methanogenesis and methanotrophy: trimethylamine => dimethylamine	Nitrogen metabolism: Bacterial (aerobic-specific) ammonia oxidation	Nitrogen metabolism: Bacterial (anaerobic-specific) ammonia oxidation	Nitrogen metabolism: Bacterial/Archaeal ammonia oxidation	Nitrogen metabolism: Dissimilatory nitrite reduction to ammonia (DNRA)	Nitrogen metabolism: Nitrogen fixation altennative	Nitrogen metabolism: ammonia => nitrite	Nitrogen metabolism: nitrate => nitrite	Nitrogen metabolism: nitric oxide => nitrous oxide	Nitrogen metabolism: nitrite => nitrate	Nitrogen metabolism: nitrite => nitric oxide	Nitrogen metabolism: nitrogen => ammonia	Nitrogen metabolism: nitrous oxide => nitrogen	Other Reductases: TMAO reductase	Other Reductases: arsenate reduction, pt 1	Other Reductases: arsenate reduction, pt 2	Other Reductases: mercury reduction	Other Reductases: selenate/Chlorate reduction	Photosynthesis: Photosystem I	Photosynthesis: Photosystem II	SCFA and alcohol conversions: Alcohol production	SCFA and alcohol conversions: Butyrate, pt 1	SCFA and alcohol conversions: Butyrate, pt 2	SCFA and alcohol conversions: Propionate, pt 1	SCFA and alcohol conversions: Propionate, pt 2	SCFA and alcohol conversions: acetate, pt 1	SCFA and alcohol conversions: acetate, pt 2	SCFA and alcohol conversions: acetate, pt 3	SCFA and alcohol conversions: lactate D	SCFA and alcohol conversions: lactate L	SCFA and alcohol conversions: pyruvate => acetyl CoA v1	SCFA and alcohol conversions: pyruvate => acetyl CoA v2	SCFA and alcohol conversions: pyruvate => acetylCoA f+ formate v3	Sulfur metabolism: Thiosulfate oxidation by SOX complex, thiosulfate => sulfate	Sulfur metabolism: dissimilatory sulfate reduction (and oxidation) sulfate => sulfide	Sulfur metabolism: tetrathionate => thiosulfate	Sulfur metabolism: thiosulfate => sulfite."""

DEFAULT_BACTERIA_PHENOTYPE = """The name of my bacteria is “Tannerella forsythia” and it utilizes potent proteolysis as a primary virulence mechanism to drive tissue damage and periodontal disease. The bacteria secretes KLIKK proteases—including karilysin and mirolysin—which break down host connective tissue (collagen, fibrinogen) and inactivate immune defenses like the complement system."""

DEFAULT_ANALYSIS = """Which pathways are active for 2(b). Also please provide reasons."""

DEFAULT_REPORTING = """Please provide a table with where first column is metabolic pathway term and second column is link to the phenotype. Also provide a summary paragraph."""


# ============================================================
# Build final prompt
# ============================================================
def build_prompt(ai_role, pathways, bacteria_phenotype, analysis, reporting):
    final_prompt = f"""
1. AI role: {ai_role}

2. Input type:
2a List of pathways {pathways}

2b {bacteria_phenotype}

3. Analysis: {analysis}

4. Reporting instructions: {reporting}
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
    page_title="GenoAnno Pathway Expectation App",
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
    </style>
    """,
    unsafe_allow_html=True,
)

st.markdown(
    """
    <div class="main-header">
        <h1>Making Expectation</h1>
    </div>
    """,
    unsafe_allow_html=True,
)

st.markdown('<div class="section-title">1. AI role</div>', unsafe_allow_html=True)
st.markdown(
    '<div class="section-note">Define the scientific expertise and perspective the model should use.</div>',
    unsafe_allow_html=True,
)

ai_role = st.text_area(
    "AI role",
    value=DEFAULT_AI_ROLE,
    height=110,
    label_visibility="collapsed",
)

st.markdown('<div class="section-title">2a. List of pathways</div>', unsafe_allow_html=True)
st.markdown(
    '<div class="section-note">Paste or edit the pathway terms that should be evaluated.</div>',
    unsafe_allow_html=True,
)

pathways = st.text_area(
    "Pathway list",
    value=DEFAULT_PATHWAYS,
    height=360,
    label_visibility="collapsed",
)

st.markdown('<div class="section-title">2b. Bacteria and phenotype</div>', unsafe_allow_html=True)
st.markdown(
    '<div class="section-note">Describe the organism and phenotype or biological property being analyzed.</div>',
    unsafe_allow_html=True,
)

bacteria_phenotype = st.text_area(
    "Bacteria and phenotype",
    value=DEFAULT_BACTERIA_PHENOTYPE,
    height=190,
    label_visibility="collapsed",
)

left_col, right_col = st.columns(2, gap="large")

with left_col:
    st.markdown('<div class="section-title">3. Analysis</div>', unsafe_allow_html=True)
    st.markdown(
        '<div class="section-note">State the specific biological question.</div>',
        unsafe_allow_html=True,
    )

    analysis = st.text_area(
        "Analysis question",
        value=DEFAULT_ANALYSIS,
        height=160,
        label_visibility="collapsed",
    )

with right_col:
    st.markdown('<div class="section-title">4. Reporting instructions</div>', unsafe_allow_html=True)
    st.markdown(
        '<div class="section-note">Define the expected output structure.</div>',
        unsafe_allow_html=True,
    )

    reporting = st.text_area(
        "Reporting instructions",
        value=DEFAULT_REPORTING,
        height=160,
        label_visibility="collapsed",
    )

final_prompt = build_prompt(
    ai_role=ai_role.strip(),
    pathways=pathways.strip(),
    bacteria_phenotype=bacteria_phenotype.strip(),
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

generate_button = st.button("Generate analysis")

if generate_button:
    if not ai_role.strip():
        st.warning("Please enter the AI role.")
    elif not pathways.strip():
        st.warning("Please enter the pathway list.")
    elif not bacteria_phenotype.strip():
        st.warning("Please enter the bacteria and phenotype.")
    elif not analysis.strip():
        st.warning("Please enter the analysis question.")
    elif not reporting.strip():
        st.warning("Please enter the reporting instructions.")
    else:
        with st.spinner("Generating expectation table..."):
            try:
                output = generate_output(
                    final_prompt=final_prompt,
                    ai_role=ai_role.strip(),
                )

                st.session_state["output"] = output

            except Exception as error:
                st.error(f"Generation failed:\n\n{error}")


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
            value="pathway_expectation_output",
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
