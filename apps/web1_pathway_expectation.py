# ============================================================
# GenoAnno Master Streamlit App - Clickable Tool Dashboard
# Runs six separate Streamlit apps from one professional home screen.
# ============================================================

from __future__ import annotations

import hashlib
import inspect
import runpy
from contextlib import contextmanager
from pathlib import Path
from typing import Callable, Dict, List

import streamlit as st


# ============================================================
# Page configuration: Streamlit allows this only once.
# Child app calls to st.set_page_config(...) are ignored safely.
# ============================================================
st.set_page_config(
    page_title="GenoAnno",
    layout="wide",
    initial_sidebar_state="collapsed",
)


BASE_DIR = Path(__file__).resolve().parent
APPS_DIR = BASE_DIR / "apps"


APP_FILES: Dict[str, str] = {
    "Making Expectation": "web1_pathway_expectation.py",
    "Metabolic Pathway Interpreter": "web2_phenotype_grouping.py",
    "KEGG Pathway Interpreter": "web3_kegg_pathway_count.py",
    "Protein Family Interpreter": "web4_protein_family_count.py",
    "Functional Interpreter": "web5_bakta_functional_category.py",
    "Similar Oral Bacteria Interpreter": "web6_similar_oral_bacteria.py",
}


APP_DESCRIPTIONS: Dict[str, str] = {
    "Making Expectation": "Predict phenotype-relevant active pathways from known bacterial behavior and pathway lists.",
    "Metabolic Pathway Interpreter": "Filter annotated metabolic pathway tables and organize retained pathways into biological phenotype groups. Use products.tsv for your bacteria from KBASE-Dram",
    "KEGG Pathway Interpreter": "Count repeated KEGG pathway terms and interpret dominant KEGG-derived functional pathway signals. Use annotation.tsv for your bacteria from KBASE-DRAM",
    "Protein Family Interpreter": "Summarize repeated protein-family terms into phenotype-level functional interpretations. Use annotation.tsv for your bacteria from KBASE-DRAM",
    "Functional Interpreter": "Parse product annotations, count functional categories, and generate phenotype summaries. Use annotation.tsv from Proksee/Bacta",
    "Similar Oral Bacteria Interpreter": "Compare functional-category patterns against well-characterized oral bacteria. Use annotation.tsv from Proksee/Bacta"
}


APP_NUMBERS: Dict[str, str] = {
    "Making Expectation": "01",
    "Metabolic Pathway Interpreter": "02",
    "KEGG Pathway Interpreter": "03",
    "Protein Family Interpreter": "04",
    "Functional Interpreter": "05",
    "Similar Oral Bacteria Interpreter": "06",
}


APP_GROUPS = [
    ("Prediction", [
        "Making Expectation",
    ]),
    ("Pathway interpretation", [
        "Metabolic Pathway Interpreter",
        "KEGG Pathway Interpreter",
    ]),
    ("Gene / protein interpretation", [
        "Protein Family Interpreter",
        "Functional Interpreter",
    ]),
    ("Comparative analysis", [
        "Similar Oral Bacteria Interpreter",
    ]),
]

APP_SHORT_NAMES: Dict[str, str] = {
    "Making Expectation": "Making Expectation",
    "Metabolic Pathway Interpreter": "Metabolic Pathway",
    "KEGG Pathway Interpreter": "KEGG Pathway",
    "Protein Family Interpreter": "Protein Family",
    "Functional Interpreter": "Gene Function",
    "Similar Oral Bacteria Interpreter": "Similar Oral Bacteria",
}


MODEL_OPTIONS: List[str] = [
    "gpt-4o-mini",
    "gpt-4.1-mini",
    "gpt-4.1",
    "gpt-4o",
    "gpt-5-mini",
    "gpt-5",
    "Custom model name",
]


WIDGET_FUNCTIONS_TO_NAMESPACE: List[str] = [
    "button",
    "checkbox",
    "color_picker",
    "date_input",
    "download_button",
    "file_uploader",
    "multiselect",
    "number_input",
    "radio",
    "select_slider",
    "selectbox",
    "slider",
    "text_area",
    "text_input",
    "time_input",
    "toggle",
]


CONTAINER_FUNCTIONS_TO_NAMESPACE: List[str] = [
    "expander",
    "form",
]


# ============================================================
# Visual system
# ============================================================
def inject_css() -> None:
    st.markdown(
        """
<style>
:root {
    --ga-bg: #f6f8fb;
    --ga-surface: #ffffff;
    --ga-surface-2: #f9fbfd;
    --ga-text: #102033;
    --ga-muted: #617083;
    --ga-muted-2: #94a3b8;
    --ga-border: #dfe6ee;
    --ga-border-strong: #c9d3df;

    /* Core scientific palette */
    --ga-navy: #15345b;
    --ga-navy-deep: #0f2744;
    --ga-teal: #147d78;
    --ga-teal-deep: #0f625e;

    /* Premium accent: used sparingly */
    --ga-gold: #b58a2a;
    --ga-gold-deep: #8f6b1f;
    --ga-gold-soft: #fbf7ea;

    --ga-primary: var(--ga-navy);
    --ga-primary-700: var(--ga-navy);
    --ga-primary-800: var(--ga-navy-deep);
    --ga-primary-soft: #eef4fa;
    --ga-accent-soft: #eef9f7;

    --ga-success: #16735d;
    --ga-success-soft: #effaf6;
    --ga-shadow-xs: 0 1px 2px rgba(15, 23, 42, 0.025);
    --ga-shadow-sm: 0 7px 22px rgba(15, 23, 42, 0.05);
    --ga-shadow-md: 0 18px 48px rgba(15, 23, 42, 0.07);
}

html {
    scroll-behavior: smooth;
}

.stApp {
    background: var(--ga-bg) !important;
    color: var(--ga-text) !important;
}

header[data-testid="stHeader"] {
    background: rgba(246, 248, 251, 0.97) !important;
    backdrop-filter: blur(14px) !important;
    border-bottom: 1px solid rgba(226, 232, 240, 0.95) !important;
    z-index: 999 !important;
}

#MainMenu, footer {
    visibility: hidden !important;
}

.block-container {
    max-width: 1540px !important;
    padding-top: 4.55rem !important;
    padding-left: 1.15rem !important;
    padding-right: 1.15rem !important;
    padding-bottom: 3rem !important;
}

/* ============================================================
   TOP APP BAR
   ============================================================ */
.ga-topbar {
    position: relative;
    overflow: hidden;
    display: flex;
    align-items: center;
    justify-content: space-between;
    gap: 1rem;
    padding: 0.95rem 1.15rem;
    margin-bottom: 0.85rem;
    border: 1px solid var(--ga-border);
    border-radius: 14px;
    background: var(--ga-surface);
    box-shadow: var(--ga-shadow-xs);
}

.ga-topbar::before {
    content: "";
    position: absolute;
    left: 0;
    top: 0;
    bottom: 0;
    width: 4px;
    background: linear-gradient(
        180deg,
        var(--ga-navy) 0%,
        var(--ga-teal) 58%,
        var(--ga-gold) 100%
    );
}

.ga-brand-kicker {
    color: var(--ga-gold-deep);
    font-size: 0.66rem;
    font-weight: 900;
    letter-spacing: 0.12em;
    text-transform: uppercase;
    margin-bottom: 0.14rem;
}

.ga-brand-title {
    margin: 0;
    color: var(--ga-text);
    font-size: 1.78rem;
    line-height: 1.02;
    font-weight: 950;
    letter-spacing: -0.055em;
}

.ga-brand-subtitle {
    margin-top: 0.28rem;
    color: var(--ga-muted);
    font-size: 0.86rem;
    line-height: 1.42;
    max-width: 980px;
}

.ga-status-pill {
    flex: 0 0 auto;
    display: inline-flex;
    align-items: center;
    gap: 0.38rem;
    padding: 0.4rem 0.66rem;
    border-radius: 999px;
    border: 1px solid #cbe8df;
    background: var(--ga-success-soft);
    color: var(--ga-success);
    font-size: 0.72rem;
    font-weight: 850;
}

.ga-status-dot {
    width: 0.43rem;
    height: 0.43rem;
    border-radius: 50%;
    background: var(--ga-teal);
}

/* ============================================================
   LEFT RAIL
   ============================================================ */
.ga-nav-shell {
    position: sticky;
    top: 5.25rem;
    padding: 0.78rem;
    border-radius: 14px;
    background: var(--ga-surface);
    border: 1px solid var(--ga-border);
    border-top: 3px solid var(--ga-gold);
    box-shadow: var(--ga-shadow-sm);
}

.ga-nav-head {
    padding: 0.1rem 0.12rem 0.55rem 0.12rem;
}

.ga-nav-title {
    color: var(--ga-text);
    font-size: 0.88rem;
    font-weight: 900;
    margin-bottom: 0.1rem;
}

.ga-nav-subtitle {
    color: var(--ga-muted);
    font-size: 0.73rem;
    line-height: 1.4;
}

.ga-nav-group {
    margin-top: 0.65rem;
    margin-bottom: 0.22rem;
    color: var(--ga-gold-deep);
    font-size: 0.63rem;
    font-weight: 900;
    letter-spacing: 0.095em;
    text-transform: uppercase;
}

.ga-nav-divider {
    height: 1px;
    background: var(--ga-border);
    margin: 0.68rem 0;
}

.ga-connection {
    padding: 0.7rem 0.72rem;
    border-radius: 10px;
    background: linear-gradient(135deg, #f7fbfa 0%, #fbfaf5 100%);
    border: 1px solid var(--ga-border);
    margin-top: 0.65rem;
}

.ga-connection-row {
    display: flex;
    align-items: center;
    gap: 0.4rem;
    margin-bottom: 0.22rem;
}

.ga-connection-dot {
    width: 0.42rem;
    height: 0.42rem;
    border-radius: 50%;
    background: var(--ga-teal);
}

.ga-connection-title {
    color: var(--ga-text);
    font-size: 0.7rem;
    font-weight: 900;
}

.ga-connection-meta {
    color: var(--ga-muted);
    font-size: 0.66rem;
    line-height: 1.42;
}

/* ============================================================
   WORKSPACE HEADER + STATUS
   ============================================================ */
.ga-workspace-head {
    position: relative;
    overflow: hidden;
    padding: 0.9rem 1rem 0.9rem 1.08rem;
    border-radius: 13px;
    background:
        linear-gradient(135deg, rgba(21, 52, 91, 0.025), rgba(20, 125, 120, 0.035)),
        var(--ga-surface);
    border: 1px solid var(--ga-border);
    box-shadow: var(--ga-shadow-xs);
    margin-bottom: 0.6rem;
}

.ga-workspace-head::before {
    content: "";
    position: absolute;
    top: 0;
    left: 0;
    bottom: 0;
    width: 3px;
    background: linear-gradient(
        180deg,
        var(--ga-navy) 0%,
        var(--ga-teal) 65%,
        var(--ga-gold) 100%
    );
}

.ga-workspace-statusline {
    display: flex;
    align-items: center;
    gap: 0.45rem;
    color: var(--ga-muted);
    font-size: 0.69rem;
    font-weight: 750;
    margin-bottom: 0.22rem;
}

.ga-workspace-statusline .ga-dot {
    width: 0.39rem;
    height: 0.39rem;
    border-radius: 50%;
    background: var(--ga-teal);
}

.ga-workspace-title {
    margin-top: 0.04rem;
    color: var(--ga-text);
    font-size: 1.3rem;
    font-weight: 950;
    letter-spacing: -0.042em;
}

.ga-workspace-desc {
    margin-top: 0.18rem;
    color: var(--ga-muted);
    font-size: 0.84rem;
    line-height: 1.48;
}

/* ============================================================
   HORIZONTAL WORKFLOW
   ============================================================ */
.ga-flow {
    display: grid;
    grid-template-columns: repeat(4, minmax(0, 1fr));
    gap: 0.5rem;
    margin: 0.55rem 0 0.7rem 0;
}

.ga-flow-step {
    position: relative;
    padding: 0.68rem 0.72rem;
    min-height: 3.25rem;
    border-radius: 10px;
    background: var(--ga-surface-2);
    border: 1px solid var(--ga-border);
}

.ga-flow-num {
    color: var(--ga-gold-deep);
    font-size: 0.62rem;
    font-weight: 900;
    letter-spacing: 0.08em;
    text-transform: uppercase;
}

.ga-flow-label {
    margin-top: 0.13rem;
    color: var(--ga-text);
    font-size: 0.75rem;
    font-weight: 850;
    line-height: 1.28;
}

/* ============================================================
   EMPTY STATE
   ============================================================ */
.ga-empty-state {
    padding: 3.2rem 1.5rem;
    border: 1px dashed var(--ga-border-strong);
    border-radius: 13px;
    background: rgba(255, 255, 255, 0.74);
    text-align: center;
}

.ga-empty-title {
    color: var(--ga-text);
    font-weight: 900;
    font-size: 1.02rem;
}

.ga-empty-desc {
    color: var(--ga-muted);
    margin-top: 0.25rem;
    font-size: 0.84rem;
    line-height: 1.5;
}

/* ============================================================
   CHILD APP WRAPPER
   ============================================================ */
.ga-app-frame {
    padding: 0.72rem;
    border-radius: 14px;
    border: 1px solid var(--ga-border);
    background: var(--ga-surface);
    box-shadow: var(--ga-shadow-sm);
}

/* Suppress duplicate child hero/header blocks */
.main-header {
    display: none !important;
}

/* Normalize child section hierarchy */
.ga-app-frame .section-title {
    color: var(--ga-text) !important;
    font-size: 1.02rem !important;
    font-weight: 900 !important;
    letter-spacing: -0.025em !important;
    margin-top: 1rem !important;
    margin-bottom: 0.14rem !important;
}

.ga-app-frame .section-note {
    color: var(--ga-muted) !important;
    font-size: 0.78rem !important;
    line-height: 1.42 !important;
    margin-bottom: 0.55rem !important;
}

.ga-app-frame .metric-card {
    border-radius: 10px !important;
    border: 1px solid var(--ga-border) !important;
    background: var(--ga-surface-2) !important;
    box-shadow: none !important;
    padding: 0.78rem 0.82rem !important;
}

.ga-app-frame .metric-card .label {
    color: var(--ga-muted) !important;
    font-size: 0.68rem !important;
    font-weight: 800 !important;
}

.ga-app-frame .metric-card .value {
    color: var(--ga-text) !important;
    font-size: 1.25rem !important;
    font-weight: 950 !important;
}

.ga-app-frame .output-panel {
    border-radius: 12px !important;
    border: 1px solid #d4dfeb !important;
    border-top: 3px solid var(--ga-gold) !important;
    background: linear-gradient(180deg, #ffffff 0%, #fbfcfe 100%) !important;
    box-shadow: none !important;
    padding: 1rem 1.05rem !important;
}

.ga-app-frame .output-panel h2 {
    font-size: 1.14rem !important;
    color: var(--ga-text) !important;
    letter-spacing: -0.025em !important;
}

.ga-app-frame .save-panel {
    border-radius: 12px !important;
    border: 1px solid var(--ga-border) !important;
    background: var(--ga-surface-2) !important;
    box-shadow: none !important;
    padding: 0.95rem 1rem !important;
}

.ga-app-frame .save-panel h3 {
    font-size: 0.96rem !important;
    color: var(--ga-text) !important;
    letter-spacing: -0.02em !important;
}

.ga-app-frame .small-muted {
    color: var(--ga-muted) !important;
    font-size: 0.76rem !important;
}

.ga-app-frame .divider {
    background: var(--ga-border) !important;
    margin: 1rem 0 !important;
}

/* ============================================================
   STREAMLIT COMPONENTS
   ============================================================ */
.stButton > button,
.stDownloadButton > button {
    border-radius: 9px !important;
    border: 1px solid transparent !important;
    background: linear-gradient(135deg, var(--ga-navy), var(--ga-teal-deep)) !important;
    color: white !important;
    font-weight: 800 !important;
    box-shadow: none !important;
    min-height: 2.35rem;
    transition: background 0.12s ease, border-color 0.12s ease !important;
}

.stButton > button:hover,
.stDownloadButton > button:hover {
    background: linear-gradient(135deg, var(--ga-navy-deep), var(--ga-teal-deep)) !important;
    border-color: var(--ga-navy-deep) !important;
}

button[kind="primary"],
button[data-testid="baseButton-primary"] {
    background: linear-gradient(135deg, var(--ga-navy), var(--ga-teal-deep)) !important;
    color: #ffffff !important;
    border: 1px solid transparent !important;
    box-shadow: 0 5px 14px rgba(21, 52, 91, 0.13) !important;
}

button[kind="primary"]:hover,
button[data-testid="baseButton-primary"]:hover {
    background: linear-gradient(135deg, var(--ga-navy-deep), var(--ga-teal-deep)) !important;
}

button[kind="secondary"],
button[data-testid="baseButton-secondary"] {
    background: #ffffff !important;
    color: var(--ga-text) !important;
    border: 1px solid var(--ga-border-strong) !important;
}

button[kind="secondary"]:hover,
button[data-testid="baseButton-secondary"]:hover {
    background: var(--ga-surface-2) !important;
    color: var(--ga-navy-deep) !important;
    border-color: #94a3b8 !important;
}

div[data-testid="stTextArea"] textarea,
div[data-testid="stTextInput"] input,
div[data-testid="stNumberInput"] input,
div[data-testid="stSelectbox"] > div > div,
div[data-testid="stFileUploader"] {
    border-radius: 9px !important;
    border-color: var(--ga-border) !important;
}

div[data-testid="stFileUploader"] {
    background: linear-gradient(180deg, #fbfcfe 0%, #f7fafc 100%) !important;
    border: 1px dashed var(--ga-border-strong) !important;
}

div[data-testid="stExpander"] {
    border: 1px solid var(--ga-border) !important;
    border-radius: 10px !important;
    background: var(--ga-surface) !important;
    box-shadow: none !important;
}

div[data-testid="stDataFrame"] {
    border: 1px solid var(--ga-border) !important;
    border-top: 2px solid rgba(21, 52, 91, 0.22) !important;
    border-radius: 10px !important;
    overflow: hidden !important;
}

.ga-app-frame div[data-testid="stMarkdownContainer"] p {
    line-height: 1.52;
}

/* Make labels more compact */
label[data-testid="stWidgetLabel"] p {
    font-size: 0.8rem !important;
    font-weight: 800 !important;
}

/* ============================================================
   RESPONSIVE
   ============================================================ */
@media (max-width: 1050px) {
    .ga-flow {
        grid-template-columns: repeat(2, minmax(0, 1fr));
    }
}

@media (max-width: 980px) {
    .ga-topbar {
        display: block;
    }

    .ga-status-pill {
        display: inline-flex;
        margin-top: 0.7rem;
    }

    .block-container {
        padding-top: 4.35rem !important;
        padding-left: 0.75rem !important;
        padding-right: 0.75rem !important;
    }

    .ga-nav-shell {
        position: relative;
        top: auto;
    }
}

@media (max-width: 680px) {
    .ga-flow {
        grid-template-columns: 1fr;
    }
}

/* Premium focus states */
div[data-testid="stTextArea"] textarea:focus,
div[data-testid="stTextInput"] input:focus,
div[data-testid="stNumberInput"] input:focus {
    border-color: var(--ga-teal) !important;
    box-shadow: 0 0 0 2px rgba(20, 125, 120, 0.10) !important;
}

/* High-value headings use navy, never pure black */
.ga-brand-title,
.ga-workspace-title,
.ga-nav-title,
.ga-empty-title,
.ga-app-frame .section-title {
    color: var(--ga-navy-deep) !important;
}

/* Gold is reserved for orientation and premium cues */
.ga-brand-kicker,
.ga-nav-group,
.ga-flow-num {
    color: var(--ga-gold-deep) !important;
}

</style>
        """.strip(),
        unsafe_allow_html=True,
    )


# ============================================================
# Child-app isolation
# ============================================================
def _stable_auto_key(prefix: str, function_name: str, args: tuple, kwargs: dict) -> str:
    """Create a stable widget key using app prefix, function, label, caller file, and caller line."""
    label = ""
    if args:
        label = str(args[0])
    elif "label" in kwargs:
        label = str(kwargs["label"])

    caller_bits = []
    for frame in inspect.stack()[2:8]:
        filename = Path(frame.filename).name
        if filename.startswith("web") or filename == "master_app.py":
            caller_bits.append(f"{filename}:{frame.lineno}")
            break

    raw = f"{prefix}|{function_name}|{label}|{'|'.join(caller_bits)}"
    digest = hashlib.md5(raw.encode("utf-8")).hexdigest()[:12]
    return f"{prefix}_{function_name}_{digest}"


@contextmanager
def streamlit_app_isolation(prefix: str):
    """
    Temporarily patches Streamlit while one child app is executed.

    Solves:
    1. Duplicate st.set_page_config calls from independent apps.
    2. Duplicate widget keys/labels across independent app files.
    """
    originals: Dict[str, Callable] = {}

    originals["set_page_config"] = st.set_page_config
    st.set_page_config = lambda *args, **kwargs: None

    def patch_widget_function(function_name: str) -> None:
        if not hasattr(st, function_name):
            return

        original_function = getattr(st, function_name)
        originals[function_name] = original_function

        def wrapped(*args, **kwargs):
            if "key" not in kwargs or kwargs["key"] is None:
                kwargs["key"] = _stable_auto_key(prefix, function_name, args, kwargs)
            else:
                kwargs["key"] = f"{prefix}_{kwargs['key']}"
            return original_function(*args, **kwargs)

        setattr(st, function_name, wrapped)

    def patch_container_function(function_name: str) -> None:
        if not hasattr(st, function_name):
            return

        original_function = getattr(st, function_name)
        originals[function_name] = original_function

        def wrapped(*args, **kwargs):
            if function_name == "form":
                if "key" not in kwargs or kwargs["key"] is None:
                    kwargs["key"] = _stable_auto_key(prefix, function_name, args, kwargs)
                else:
                    kwargs["key"] = f"{prefix}_{kwargs['key']}"
            return original_function(*args, **kwargs)

        setattr(st, function_name, wrapped)

    for widget_name in WIDGET_FUNCTIONS_TO_NAMESPACE:
        patch_widget_function(widget_name)

    for container_name in CONTAINER_FUNCTIONS_TO_NAMESPACE:
        patch_container_function(container_name)

    try:
        yield
    finally:
        for name, original in originals.items():
            setattr(st, name, original)


# ============================================================
# App utilities
# ============================================================
def ensure_config_available() -> None:
    """
    Website-only configuration version.

    No config.txt file is required. OpenAI settings are collected from the
    main website form and stored in st.session_state for the current browser
    session. This function is kept as a harmless no-op for compatibility.
    """
    return None


def config_has_key() -> bool:
    if st.session_state.get("user_openai_api_key", ""):
        return True

    try:
        if st.secrets.get("KEY", ""):
            return True
    except Exception:
        pass

    import os
    return bool(os.getenv("OPENAI_API_KEY", "").strip())


def count_available_child_apps() -> int:
    return sum(1 for filename in APP_FILES.values() if (APPS_DIR / filename).exists())


def app_prefix(app_name: str) -> str:
    names = list(APP_FILES.keys())
    return f"app{names.index(app_name) + 1}"


def clear_child_transient_state(previous_app: str | None = None) -> None:
    """
    Clear analysis-only state when moving between tools while preserving the
    user's OpenAI connection settings.
    """
    transient_exact_keys = {
        "output",
        "analysis_output",
        "generated_output",
        "model_output",
        "result",
        "results",
        "processed_table_text",
        "top10_table_text",
        "top40_table_text",
    }

    for key in list(st.session_state.keys()):
        if key in transient_exact_keys:
            st.session_state.pop(key, None)

    # Remove namespaced child-widget values belonging to the previous tool.
    # This prevents stale file selections / controls from visually carrying
    # into a fresh visit to that tool after switching away.
    if previous_app in APP_FILES:
        prefix = app_prefix(previous_app) + "_"
        protected = {
            "user_openai_api_key",
            "user_selected_model",
            "user_temperature",
            "user_max_tokens",
            "active_app",
        }

        for key in list(st.session_state.keys()):
            if key in protected:
                continue
            if str(key).startswith(prefix):
                st.session_state.pop(key, None)


def switch_active_app(app_name: str) -> None:
    previous_app = st.session_state.get("active_app")

    if previous_app != app_name:
        clear_child_transient_state(previous_app=previous_app)
        st.session_state["active_app"] = app_name


def run_child_app(app_name: str) -> None:
    filename = APP_FILES[app_name]
    app_path = APPS_DIR / filename

    if not app_path.exists():
        st.error(f"Missing app file: {app_path}")
        return

    prefix = app_prefix(app_name)

    with streamlit_app_isolation(prefix):
        runpy.run_path(str(app_path), run_name=f"__genoanno_{prefix}__")



def api_key_setup() -> bool:
    """
    Collect OpenAI settings for the current browser session.

    When already connected, this function stays visually silent. Connection
    status and settings are shown compactly in the left navigation rail.
    """

    if st.session_state.get("user_openai_api_key", ""):
        return True

    st.markdown(
        """
<div class="ga-topbar">
    <div>
        <div class="ga-brand-kicker">GenoAnno setup</div>
        <h1 class="ga-brand-title">Connect OpenAI</h1>
        <div class="ga-brand-subtitle">
            Enter your API key and model settings once. They are stored only in
            the current Streamlit browser session.
        </div>
    </div>
</div>
        """.strip(),
        unsafe_allow_html=True,
    )

    with st.container(border=True):
        user_key = st.text_input(
            "OpenAI API key",
            type="password",
            placeholder="sk-...",
            key="user_openai_api_key_input",
            help="Stored only in the current Streamlit session.",
        )

        col_model, col_temp, col_tokens = st.columns([2, 1, 1])

        with col_model:
            selected_model_option = st.selectbox(
                "Model",
                MODEL_OPTIONS,
                index=0,
                key="user_model_select",
            )

            custom_model = ""
            if selected_model_option == "Custom model name":
                custom_model = st.text_input(
                    "Custom model name",
                    placeholder="Example: gpt-4o-mini",
                    key="user_custom_model_name",
                )

        with col_temp:
            selected_temperature = st.number_input(
                "Temperature",
                min_value=0.0,
                max_value=2.0,
                value=0.5,
                step=0.1,
                key="user_temperature_input",
            )

        with col_tokens:
            selected_max_tokens = st.number_input(
                "Max tokens",
                min_value=256,
                max_value=16000,
                value=2000,
                step=256,
                key="user_max_tokens_input",
            )

        save_clicked = st.button(
            "Connect and open GenoAnno",
            key="save_user_api_settings",
            use_container_width=True,
        )

        if save_clicked:
            cleaned_key = user_key.strip()

            if not cleaned_key:
                st.error("Please enter an OpenAI API key.")
                return False

            final_model = (
                custom_model.strip()
                if selected_model_option == "Custom model name"
                else selected_model_option
            )

            if not final_model:
                st.error("Please choose a model or enter a custom model name.")
                return False

            st.session_state["user_openai_api_key"] = cleaned_key
            st.session_state["user_selected_model"] = final_model
            st.session_state["user_temperature"] = float(selected_temperature)
            st.session_state["user_max_tokens"] = int(selected_max_tokens)
            st.rerun()

    return False

# ============================================================
# Page sections
# ============================================================
def workflow_for_app(app_name: str) -> list[str]:
    if app_name == "Making Expectation":
        return [
            "Select tool",
            "Provide bacterial behavior and pathway information",
            "Review provided inputs",
            "Generate pathway expectations",
        ]

    return [
        "Select tool",
        "Upload annotation",
        "Review processed data",
        "Generate interpretation",
    ]


def render_workflow(app_name: str) -> None:
    steps = workflow_for_app(app_name)

    cards = []
    for idx, label in enumerate(steps, start=1):
        cards.append(
            f"""
<div class="ga-flow-step">
    <div class="ga-flow-num">Step {idx}</div>
    <div class="ga-flow-label">{label}</div>
</div>
            """.strip()
        )

    st.markdown(
        '<div class="ga-flow">' + "".join(cards) + "</div>",
        unsafe_allow_html=True,
    )


def render_hero(config_ready: bool, child_apps_ready: int) -> None:
    st.markdown(
        f"""
<div class="ga-topbar">
    <div>
        <div class="ga-brand-kicker">Genome annotation interpretation suite</div>
        <h1 class="ga-brand-title">GenoAnno</h1>
        <div class="ga-brand-subtitle">
            AI-assisted interpretation of bacterial pathways, protein families,
            gene functions, phenotypes, and oral-bacteria similarity.
        </div>
    </div>
    <div class="ga-status-pill">
        <span class="ga-status-dot"></span>
        {child_apps_ready}/{len(APP_FILES)} tools ready
    </div>
</div>
        """.strip(),
        unsafe_allow_html=True,
    )


def render_left_navigation(active_app: str | None) -> None:
    st.markdown('<div class="ga-nav-shell">', unsafe_allow_html=True)

    st.markdown(
        """
<div class="ga-nav-head">
    <div class="ga-nav-title">Analysis modules</div>
    <div class="ga-nav-subtitle">Choose a module to open it in the workspace.</div>
</div>
        """.strip(),
        unsafe_allow_html=True,
    )

    for group_name, app_names in APP_GROUPS:
        st.markdown(
            f'<div class="ga-nav-group">{group_name}</div>',
            unsafe_allow_html=True,
        )

        for app_name in app_names:
            is_active = active_app == app_name
            label = f"{APP_NUMBERS[app_name]}  {APP_SHORT_NAMES[app_name]}"

            if st.button(
                label,
                key=f"nav_{app_prefix(app_name)}",
                use_container_width=True,
                type="primary" if is_active else "secondary",
            ):
                switch_active_app(app_name)
                st.rerun()

    st.markdown('<div class="ga-nav-divider"></div>', unsafe_allow_html=True)

    active_model = st.session_state.get("user_selected_model", "gpt-4o-mini")
    active_temperature = st.session_state.get("user_temperature", 0.5)
    active_max_tokens = st.session_state.get("user_max_tokens", 2000)

    st.markdown(
        f"""
<div class="ga-connection">
    <div class="ga-connection-row">
        <span class="ga-connection-dot"></span>
        <span class="ga-connection-title">OpenAI connected</span>
    </div>
    <div class="ga-connection-meta">
        {active_model}<br>
        T={active_temperature} · {active_max_tokens} tokens
    </div>
</div>
        """.strip(),
        unsafe_allow_html=True,
    )

    with st.expander("Settings", expanded=False):
        st.caption("Stored only for this browser session.")

        if st.button(
            "Clear API settings",
            key="clear_user_api_settings",
            use_container_width=True,
            type="secondary",
        ):
            for session_key in [
                "user_openai_api_key",
                "user_selected_model",
                "user_temperature",
                "user_max_tokens",
            ]:
                st.session_state.pop(session_key, None)

            clear_child_transient_state(
                previous_app=st.session_state.get("active_app")
            )
            st.session_state["active_app"] = None
            st.rerun()

    if active_app is not None:
        if st.button(
            "Close current tool",
            key="nav_close_current",
            use_container_width=True,
            type="secondary",
        ):
            clear_child_transient_state(previous_app=active_app)
            st.session_state["active_app"] = None
            st.rerun()

    st.markdown('</div>', unsafe_allow_html=True)


def render_selected_app_header(app_name: str) -> None:
    st.markdown(
        f"""
<div class="ga-workspace-head">
    <div class="ga-workspace-statusline">
        <span class="ga-dot"></span>
        Tool {APP_NUMBERS[app_name]} · Ready
    </div>
    <div class="ga-workspace-title">{app_name}</div>
    <div class="ga-workspace-desc">{APP_DESCRIPTIONS[app_name]}</div>
</div>
        """.strip(),
        unsafe_allow_html=True,
    )

    render_workflow(app_name)

    col_clear, col_space = st.columns([1, 5])

    with col_clear:
        if st.button(
            "Reset",
            key="reset_active_tool",
            use_container_width=True,
            type="secondary",
            help=(
                "Reset the current tool to its original state. "
                "Uploaded files, generated outputs, and tool-specific temporary "
                "inputs are cleared. OpenAI settings are preserved."
            ),
        ):
            clear_child_transient_state(previous_app=app_name)
            st.rerun()


def main() -> None:
    inject_css()

    # Website-based API-key setup. No password is required.
    if not api_key_setup():
        st.stop()

    config_ready = config_has_key()
    child_apps_ready = count_available_child_apps()

    if "active_app" not in st.session_state:
        st.session_state["active_app"] = None

    active_app = st.session_state.get("active_app")

    if active_app is not None and active_app not in APP_FILES:
        clear_child_transient_state(previous_app=active_app)
        st.session_state["active_app"] = None
        active_app = None

    render_hero(
        config_ready=config_ready,
        child_apps_ready=child_apps_ready,
    )

    nav_col, workspace_col = st.columns([1.2, 4.8], gap="large")

    with nav_col:
        render_left_navigation(active_app=active_app)

    with workspace_col:
        if active_app is None:
            st.session_state.pop("output", None)

            st.markdown(
                """
<div class="ga-empty-state">
    <div class="ga-empty-title">Select a module from the left navigation</div>
    <div class="ga-empty-desc">
        GenoAnno will open the selected workflow here. Each module keeps its
        analysis state isolated so results do not carry into another tool.
    </div>
</div>
                """.strip(),
                unsafe_allow_html=True,
            )
            return

        render_selected_app_header(active_app)

        st.markdown('<div class="ga-app-frame">', unsafe_allow_html=True)
        run_child_app(active_app)
        st.markdown('</div>', unsafe_allow_html=True)


if __name__ == "__main__":
    main()
