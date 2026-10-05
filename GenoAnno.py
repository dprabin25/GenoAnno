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
    --ga-bg: #f4f6f8;
    --ga-surface: #ffffff;
    --ga-surface-soft: #f8fafc;
    --ga-text: #0f172a;
    --ga-muted: #64748b;
    --ga-muted-2: #94a3b8;
    --ga-border: #e2e8f0;
    --ga-border-strong: #cbd5e1;
    --ga-primary: #1d4ed8;
    --ga-primary-dark: #1e3a8a;
    --ga-primary-soft: #eff6ff;
    --ga-accent: #0f766e;
    --ga-success: #166534;
    --ga-shadow-xs: 0 1px 2px rgba(15, 23, 42, 0.03);
    --ga-shadow-sm: 0 6px 20px rgba(15, 23, 42, 0.055);
    --ga-shadow-md: 0 18px 48px rgba(15, 23, 42, 0.075);
}

.stApp {
    background: var(--ga-bg) !important;
    color: var(--ga-text) !important;
}

header[data-testid="stHeader"] {
    background: rgba(244, 246, 248, 0.97) !important;
    backdrop-filter: blur(14px) !important;
    border-bottom: 1px solid rgba(226, 232, 240, 0.95) !important;
    z-index: 999 !important;
}

#MainMenu, footer {
    visibility: hidden !important;
}

.block-container {
    max-width: 1500px !important;
    padding-top: 4.6rem !important;
    padding-left: 1.3rem !important;
    padding-right: 1.3rem !important;
    padding-bottom: 3rem !important;
}

/* ============================================================
   MASTER HEADER
   ============================================================ */
.ga-topbar {
    display: flex;
    align-items: center;
    justify-content: space-between;
    gap: 1rem;
    padding: 1rem 1.25rem;
    margin-bottom: 0.9rem;
    border: 1px solid var(--ga-border);
    border-radius: 16px;
    background: var(--ga-surface);
    box-shadow: var(--ga-shadow-xs);
}

.ga-brand-kicker {
    color: var(--ga-primary);
    font-size: 0.7rem;
    font-weight: 900;
    letter-spacing: 0.12em;
    text-transform: uppercase;
    margin-bottom: 0.18rem;
}

.ga-brand-title {
    margin: 0;
    color: var(--ga-text);
    font-size: 1.85rem;
    line-height: 1.02;
    font-weight: 950;
    letter-spacing: -0.055em;
}

.ga-brand-subtitle {
    margin-top: 0.34rem;
    color: var(--ga-muted);
    font-size: 0.9rem;
    line-height: 1.45;
    max-width: 900px;
}

.ga-status-pill {
    flex: 0 0 auto;
    display: inline-flex;
    align-items: center;
    gap: 0.38rem;
    padding: 0.42rem 0.7rem;
    border-radius: 999px;
    border: 1px solid #bbf7d0;
    background: #f0fdf4;
    color: #166534;
    font-size: 0.76rem;
    font-weight: 850;
}

.ga-status-dot {
    width: 0.48rem;
    height: 0.48rem;
    border-radius: 999px;
    background: #22c55e;
}

/* ============================================================
   LEFT NAVIGATION RAIL
   ============================================================ */
.ga-nav-shell {
    position: sticky;
    top: 5.35rem;
    padding: 0.9rem;
    border-radius: 16px;
    background: var(--ga-surface);
    border: 1px solid var(--ga-border);
    box-shadow: var(--ga-shadow-sm);
}

.ga-nav-brand {
    padding: 0.15rem 0.15rem 0.75rem 0.15rem;
}

.ga-nav-title {
    color: var(--ga-text);
    font-size: 0.92rem;
    font-weight: 900;
    margin-bottom: 0.12rem;
}

.ga-nav-subtitle {
    color: var(--ga-muted);
    font-size: 0.76rem;
    line-height: 1.4;
}

.ga-nav-divider {
    height: 1px;
    background: var(--ga-border);
    margin: 0.7rem 0;
}

.ga-workflow {
    padding: 0.72rem 0.75rem;
    border-radius: 11px;
    background: var(--ga-surface-soft);
    border: 1px solid var(--ga-border);
    margin-bottom: 0.75rem;
}

.ga-workflow-title {
    color: var(--ga-text);
    font-size: 0.72rem;
    font-weight: 900;
    text-transform: uppercase;
    letter-spacing: 0.075em;
    margin-bottom: 0.45rem;
}

.ga-workflow-step {
    color: var(--ga-muted);
    font-size: 0.72rem;
    line-height: 1.55;
}

.ga-workflow-step b {
    color: var(--ga-text);
}

.ga-nav-active {
    margin-top: 0.72rem;
    padding: 0.72rem 0.8rem;
    border-radius: 11px;
    background: var(--ga-primary-soft);
    border-left: 3px solid var(--ga-primary);
    border-top: 1px solid #dbeafe;
    border-right: 1px solid #dbeafe;
    border-bottom: 1px solid #dbeafe;
}

.ga-nav-active-label {
    color: var(--ga-primary);
    font-size: 0.63rem;
    font-weight: 900;
    text-transform: uppercase;
    letter-spacing: 0.09em;
}

.ga-nav-active-title {
    color: var(--ga-text);
    font-size: 0.82rem;
    font-weight: 900;
    margin-top: 0.16rem;
    line-height: 1.3;
}

.ga-connection {
    padding: 0.72rem 0.75rem;
    border-radius: 11px;
    background: #f8fafc;
    border: 1px solid var(--ga-border);
    margin-top: 0.72rem;
}

.ga-connection-row {
    display: flex;
    align-items: center;
    gap: 0.42rem;
    margin-bottom: 0.3rem;
}

.ga-connection-dot {
    width: 0.46rem;
    height: 0.46rem;
    border-radius: 50%;
    background: #22c55e;
}

.ga-connection-title {
    color: var(--ga-text);
    font-size: 0.73rem;
    font-weight: 900;
}

.ga-connection-meta {
    color: var(--ga-muted);
    font-size: 0.69rem;
    line-height: 1.45;
}

/* ============================================================
   WORKSPACE
   ============================================================ */
.ga-workspace-head {
    padding: 0.95rem 1.1rem;
    border-radius: 14px;
    background: var(--ga-surface);
    border: 1px solid var(--ga-border);
    box-shadow: var(--ga-shadow-xs);
    margin-bottom: 0.7rem;
}

.ga-workspace-kicker {
    color: var(--ga-primary);
    font-size: 0.66rem;
    font-weight: 900;
    letter-spacing: 0.1em;
    text-transform: uppercase;
}

.ga-workspace-title {
    margin-top: 0.16rem;
    color: var(--ga-text);
    font-size: 1.32rem;
    font-weight: 950;
    letter-spacing: -0.04em;
}

.ga-workspace-desc {
    margin-top: 0.22rem;
    color: var(--ga-muted);
    font-size: 0.86rem;
    line-height: 1.48;
}

.ga-empty-state {
    padding: 3rem 1.5rem;
    border: 1px dashed var(--ga-border-strong);
    border-radius: 14px;
    background: rgba(255, 255, 255, 0.72);
    text-align: center;
}

.ga-empty-title {
    color: var(--ga-text);
    font-weight: 900;
    font-size: 1.03rem;
}

.ga-empty-desc {
    color: var(--ga-muted);
    margin-top: 0.28rem;
    font-size: 0.86rem;
    line-height: 1.5;
}

.ga-app-frame {
    padding: 0.8rem;
    border-radius: 16px;
    border: 1px solid var(--ga-border);
    background: var(--ga-surface);
    box-shadow: var(--ga-shadow-sm);
}

/* ============================================================
   NORMALIZE CHILD APPS
   ============================================================ */

/* Hide duplicate child hero/header blocks. */
.main-header {
    display: none !important;
}

/* Normalize child app page spacing/background. */
.ga-app-frame .block-container {
    padding: 0 !important;
    max-width: 100% !important;
}

.ga-app-frame .section-title {
    color: var(--ga-text) !important;
    font-size: 1.08rem !important;
    font-weight: 900 !important;
    letter-spacing: -0.025em !important;
    margin-top: 1.15rem !important;
    margin-bottom: 0.18rem !important;
}

.ga-app-frame .section-note {
    color: var(--ga-muted) !important;
    font-size: 0.82rem !important;
    line-height: 1.45 !important;
    margin-bottom: 0.65rem !important;
}

.ga-app-frame .metric-card {
    border-radius: 12px !important;
    border: 1px solid var(--ga-border) !important;
    background: var(--ga-surface-soft) !important;
    box-shadow: none !important;
    padding: 0.85rem 0.9rem !important;
}

.ga-app-frame .metric-card .label {
    color: var(--ga-muted) !important;
    font-size: 0.72rem !important;
    font-weight: 800 !important;
}

.ga-app-frame .metric-card .value {
    color: var(--ga-text) !important;
    font-size: 1.35rem !important;
    font-weight: 950 !important;
}

.ga-app-frame .output-panel,
.ga-app-frame .save-panel {
    border-radius: 13px !important;
    border: 1px solid var(--ga-border) !important;
    background: var(--ga-surface) !important;
    box-shadow: none !important;
    padding: 1rem 1.1rem !important;
}

.ga-app-frame .output-panel h2,
.ga-app-frame .save-panel h3 {
    color: var(--ga-text) !important;
    letter-spacing: -0.03em !important;
}

.ga-app-frame .output-panel h2 {
    font-size: 1.25rem !important;
}

.ga-app-frame .save-panel h3 {
    font-size: 1rem !important;
}

.ga-app-frame .divider {
    background: var(--ga-border) !important;
    margin: 1.1rem 0 !important;
}

/* ============================================================
   STREAMLIT CONTROLS
   ============================================================ */
.stButton > button,
.stDownloadButton > button {
    border-radius: 10px !important;
    border: 1px solid transparent !important;
    background: var(--ga-primary) !important;
    color: white !important;
    font-weight: 800 !important;
    box-shadow: none !important;
    min-height: 2.45rem;
    transition: background 0.12s ease, border-color 0.12s ease !important;
}

.stButton > button:hover,
.stDownloadButton > button:hover {
    background: var(--ga-primary-dark) !important;
    border-color: var(--ga-primary-dark) !important;
}

/* Make secondary Streamlit buttons actually look secondary when type="secondary". */
button[kind="secondary"],
button[data-testid="baseButton-secondary"] {
    background: #ffffff !important;
    color: var(--ga-text) !important;
    border: 1px solid var(--ga-border-strong) !important;
}

button[kind="secondary"]:hover,
button[data-testid="baseButton-secondary"]:hover {
    background: var(--ga-surface-soft) !important;
    color: var(--ga-primary-dark) !important;
    border-color: #94a3b8 !important;
}

div[data-testid="stTextArea"] textarea,
div[data-testid="stTextInput"] input,
div[data-testid="stNumberInput"] input,
div[data-testid="stSelectbox"] > div > div,
div[data-testid="stFileUploader"] {
    border-radius: 10px !important;
    border-color: var(--ga-border) !important;
}

div[data-testid="stFileUploader"] {
    background: var(--ga-surface-soft) !important;
}

div[data-testid="stExpander"] {
    border: 1px solid var(--ga-border) !important;
    border-radius: 11px !important;
    background: var(--ga-surface) !important;
    box-shadow: none !important;
}

/* Cleaner dataframes */
div[data-testid="stDataFrame"] {
    border: 1px solid var(--ga-border) !important;
    border-radius: 11px !important;
    overflow: hidden !important;
}

/* Reduce excessive markdown spacing inside child apps */
.ga-app-frame div[data-testid="stMarkdownContainer"] p {
    line-height: 1.55;
}

/* ============================================================
   RESPONSIVE
   ============================================================ */
@media (max-width: 980px) {
    .ga-topbar {
        display: block;
    }

    .ga-status-pill {
        display: inline-flex;
        margin-top: 0.75rem;
    }

    .block-container {
        padding-top: 4.4rem !important;
        padding-left: 0.8rem !important;
        padding-right: 0.8rem !important;
    }

    .ga-nav-shell {
        position: relative;
        top: auto;
    }
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
<div class="ga-nav-brand">
    <div class="ga-nav-title">Analysis tools</div>
    <div class="ga-nav-subtitle">Select a module to open it in the workspace.</div>
</div>

<div class="ga-workflow">
    <div class="ga-workflow-title">Workflow</div>
    <div class="ga-workflow-step"><b>1.</b> Select tool</div>
    <div class="ga-workflow-step"><b>2.</b> Upload annotation</div>
    <div class="ga-workflow-step"><b>3.</b> Review processed data</div>
    <div class="ga-workflow-step"><b>4.</b> Generate interpretation</div>
</div>
        """.strip(),
        unsafe_allow_html=True,
    )

    for app_name in APP_FILES:
        is_active = active_app == app_name
        label = f"{APP_NUMBERS[app_name]}  {app_name}"

        if st.button(
            label,
            key=f"nav_{app_prefix(app_name)}",
            use_container_width=True,
            type="primary" if is_active else "secondary",
        ):
            switch_active_app(app_name)
            st.rerun()

    if active_app is not None:
        st.markdown(
            f"""
<div class="ga-nav-active">
    <div class="ga-nav-active-label">Current analysis</div>
    <div class="ga-nav-active-title">{active_app}</div>
</div>
            """.strip(),
            unsafe_allow_html=True,
        )

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
        Temperature {active_temperature} · {active_max_tokens} tokens
    </div>
</div>
        """.strip(),
        unsafe_allow_html=True,
    )

    with st.expander("Settings", expanded=False):
        st.caption("Connection settings are stored only for this browser session.")

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
    <div class="ga-workspace-kicker">Analysis workspace</div>
    <div class="ga-workspace-title">{app_name}</div>
    <div class="ga-workspace-desc">{APP_DESCRIPTIONS[app_name]}</div>
</div>
        """.strip(),
        unsafe_allow_html=True,
    )

    col_clear, col_space = st.columns([1, 5])

    with col_clear:
        if st.button(
            "Clear result",
            key="clear_active_output",
            use_container_width=True,
            type="secondary",
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
    <div class="ga-empty-title">Choose an analysis module from the left</div>
    <div class="ga-empty-desc">
        The selected tool will open in this workspace. Results and temporary
        analysis state are isolated between tools to keep each workflow clean.
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
