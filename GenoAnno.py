# ============================================================
# GenoAnno Master Streamlit App - Clickable Tool Dashboard
# Runs six separate Streamlit apps from one professional home screen.
#
# How to run:
#   streamlit run master_app.py
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
    --ga-bg-1: #f7fbff;
    --ga-bg-2: #edf4ff;
    --ga-text: #0f172a;
    --ga-muted: #5f6b7a;
    --ga-border: rgba(148, 163, 184, 0.32);
    --ga-primary: #1e40af;
    --ga-primary-2: #0f766e;
    --ga-card: rgba(255, 255, 255, 0.86);
    --ga-shadow: 0 24px 70px rgba(15, 23, 42, 0.12);
    --ga-soft-shadow: 0 14px 34px rgba(15, 23, 42, 0.08);
}

.stApp {
    background:
        radial-gradient(circle at 10% 5%, rgba(30, 64, 175, 0.13), transparent 30%),
        radial-gradient(circle at 90% 4%, rgba(15, 118, 110, 0.14), transparent 32%),
        linear-gradient(135deg, var(--ga-bg-1) 0%, var(--ga-bg-2) 52%, #f9fbff 100%) !important;
    color: var(--ga-text) !important;
}

header[data-testid="stHeader"] {
    background: rgba(248, 251, 255, 0.76) !important;
    backdrop-filter: blur(18px) !important;
    border-bottom: 1px solid rgba(148, 163, 184, 0.18) !important;
}

#MainMenu, footer { visibility: hidden !important; }

.block-container {
    max-width: 94% !important;
    padding-top: 1.25rem !important;
    padding-left: 2.2rem !important;
    padding-right: 2.2rem !important;
    padding-bottom: 4rem !important;
}

.ga-hero {
    position: relative;
    overflow: hidden;
    border-radius: 30px;
    padding: 2.8rem 3rem;
    margin-bottom: 1.25rem;
    background:
        linear-gradient(135deg, rgba(15, 23, 42, 0.98), rgba(30, 58, 138, 0.95) 52%, rgba(15, 118, 110, 0.92));
    color: white;
    box-shadow: var(--ga-shadow);
    border: 1px solid rgba(255, 255, 255, 0.14);
}

.ga-hero::after {
    content: "";
    position: absolute;
    right: -80px;
    top: -100px;
    width: 420px;
    height: 420px;
    border-radius: 999px;
    background: rgba(255, 255, 255, 0.10);
}

.ga-hero-content {
    position: relative;
    z-index: 2;
    max-width: 1120px;
}

.ga-badge-row {
    display: flex;
    gap: 0.65rem;
    flex-wrap: wrap;
    margin-bottom: 1.15rem;
}

.ga-badge {
    display: inline-flex;
    align-items: center;
    padding: 0.42rem 0.82rem;
    border-radius: 999px;
    background: rgba(255, 255, 255, 0.12);
    border: 1px solid rgba(255, 255, 255, 0.22);
    color: #e0f2fe;
    font-size: 0.86rem;
    font-weight: 800;
    letter-spacing: 0.02em;
}

.ga-hero h1 {
    margin: 0;
    font-size: clamp(2.5rem, 5vw, 4.8rem);
    line-height: 0.96;
    font-weight: 950;
    letter-spacing: -0.075em;
}

.ga-hero h1 span {
    background: linear-gradient(90deg, #ffffff 0%, #c7f9ef 48%, #dbeafe 100%);
    -webkit-background-clip: text;
    color: transparent;
}

.ga-hero p {
    margin-top: 1.2rem;
    margin-bottom: 0;
    max-width: 980px;
    color: #e5efff;
    font-size: 1.14rem;
    line-height: 1.72;
}

.ga-stats {
    display: grid;
    grid-template-columns: repeat(4, minmax(0, 1fr));
    gap: 1rem;
    margin: 1.15rem 0 1.45rem 0;
}

.ga-stat-card {
    border-radius: 22px;
    padding: 1.18rem 1.22rem;
    background: var(--ga-card);
    backdrop-filter: blur(18px);
    border: 1px solid var(--ga-border);
    box-shadow: var(--ga-soft-shadow);
}

.ga-stat-label {
    color: var(--ga-muted);
    font-size: 0.82rem;
    font-weight: 850;
    text-transform: uppercase;
    letter-spacing: 0.075em;
}

.ga-stat-value {
    color: var(--ga-text);
    font-size: 1.55rem;
    font-weight: 950;
    margin-top: 0.34rem;
    letter-spacing: -0.04em;
}

.ga-section-title {
    margin-top: 0.25rem;
    margin-bottom: 0.2rem;
    color: var(--ga-text);
    font-size: 1.55rem;
    font-weight: 950;
    letter-spacing: -0.045em;
}

.ga-section-subtitle {
    color: var(--ga-muted);
    margin-bottom: 1rem;
    font-size: 1rem;
}

.ga-tool-number {
    display: inline-block;
    min-width: 2.4rem;
    text-align: center;
    padding: 0.2rem 0.55rem;
    border-radius: 999px;
    background: rgba(30, 64, 175, 0.10);
    color: #1e3a8a;
    font-weight: 950;
    font-size: 0.78rem;
    letter-spacing: 0.04em;
    margin-bottom: 0.65rem;
}

.ga-tool-title {
    color: var(--ga-text);
    font-size: 1.16rem;
    font-weight: 950;
    letter-spacing: -0.035em;
    margin-bottom: 0.4rem;
}

.ga-tool-desc {
    color: var(--ga-muted);
    font-size: 0.95rem;
    line-height: 1.55;
    min-height: 4.4rem;
}

.ga-selected-bar {
    margin-top: 1.2rem;
    margin-bottom: 0.8rem;
    border-radius: 24px;
    padding: 1.1rem 1.25rem;
    background: linear-gradient(135deg, rgba(30, 64, 175, 0.10), rgba(15, 118, 110, 0.10));
    border: 1px solid rgba(30, 64, 175, 0.16);
}

.ga-selected-label {
    color: var(--ga-muted);
    font-size: 0.78rem;
    text-transform: uppercase;
    font-weight: 900;
    letter-spacing: 0.08em;
}

.ga-selected-title {
    color: var(--ga-text);
    font-size: 1.45rem;
    font-weight: 950;
    letter-spacing: -0.045em;
    margin-top: 0.15rem;
}

.ga-selected-desc {
    color: var(--ga-muted);
    margin-top: 0.25rem;
    font-size: 0.98rem;
}



.ga-api-box {
    margin-top: 1.1rem;
    margin-bottom: 1.2rem;
    border-radius: 26px;
    padding: 1.4rem 1.5rem;
    background: rgba(255, 255, 255, 0.86);
    border: 1px solid rgba(148, 163, 184, 0.32);
    box-shadow: var(--ga-soft-shadow);
}

.ga-api-title {
    color: var(--ga-text);
    font-size: 1.45rem;
    font-weight: 950;
    letter-spacing: -0.045em;
    margin-bottom: 0.35rem;
}

.ga-api-desc {
    color: var(--ga-muted);
    font-size: 0.98rem;
    line-height: 1.6;
    margin-bottom: 0.7rem;
}

.ga-api-good {
    color: #166534;
    font-weight: 900;
}

.ga-api-warning {
    color: #92400e;
    font-weight: 900;
}

.ga-app-frame {
    background: rgba(255, 255, 255, 0.72);
    border: 1px solid rgba(148, 163, 184, 0.24);
    border-radius: 28px;
    padding: 1rem;
    box-shadow: var(--ga-soft-shadow);
    margin-top: 0.7rem;
}

.stButton > button,
.stDownloadButton > button {
    border-radius: 14px !important;
    border: 0 !important;
    background: linear-gradient(135deg, #1e40af, #0f766e) !important;
    color: white !important;
    font-weight: 900 !important;
    box-shadow: 0 13px 28px rgba(30, 64, 175, 0.20) !important;
    transition: transform 0.15s ease, box-shadow 0.15s ease !important;
}

.stButton > button:hover,
.stDownloadButton > button:hover {
    transform: translateY(-1px) !important;
    box-shadow: 0 18px 36px rgba(30, 64, 175, 0.28) !important;
}

div[data-testid="stVerticalBlockBorderWrapper"] {
    border-radius: 24px !important;
    background: rgba(255, 255, 255, 0.72) !important;
    box-shadow: 0 12px 30px rgba(15, 23, 42, 0.06) !important;
}

div[data-testid="stTextArea"] textarea,
div[data-testid="stTextInput"] input,
div[data-testid="stNumberInput"] input,
div[data-testid="stSelectbox"] div,
div[data-testid="stFileUploader"] {
    border-radius: 16px !important;
    border-color: rgba(148, 163, 184, 0.38) !important;
}

div[data-testid="stFileUploader"] {
    background: rgba(255, 255, 255, 0.78) !important;
    box-shadow: 0 10px 26px rgba(15, 23, 42, 0.05) !important;
}

@media (max-width: 1100px) {
    .ga-stats { grid-template-columns: repeat(2, minmax(0, 1fr)); }
    .ga-hero { padding: 2rem; }
}

@media (max-width: 740px) {
    .ga-stats { grid-template-columns: 1fr; }
    .block-container { padding-left: 1rem !important; padding-right: 1rem !important; }
    .ga-hero { border-radius: 22px; padding: 1.6rem; }
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
    Collect OpenAI settings from the website user.

    Stored only in st.session_state for the current browser session:
    - user_openai_api_key
    - user_selected_model
    - user_temperature
    - user_max_tokens
    Nothing is written to GitHub or any local app file.
    """

    if st.session_state.get("user_openai_api_key", ""):
        active_model = st.session_state.get("user_selected_model", "gpt-4o-mini")
        active_temperature = st.session_state.get("user_temperature", 0.5)
        active_max_tokens = st.session_state.get("user_max_tokens", 2000)

        st.markdown(
            f"""
<div class="ga-api-box">
    <div class="ga-api-title">OpenAI settings active for this session</div>
    <div class="ga-api-desc">
        Your API key is stored only in this browser session. Current model: <b>{active_model}</b>;
        temperature: <b>{active_temperature}</b>; max tokens: <b>{active_max_tokens}</b>.
        You can clear these settings below before sharing the screen or closing work.
    </div>
</div>
            """.strip(),
            unsafe_allow_html=True,
        )

        if st.button("Clear API settings from this session", key="clear_user_api_settings"):
            for session_key in [
                "user_openai_api_key",
                "user_selected_model",
                "user_temperature",
                "user_max_tokens",
            ]:
                st.session_state.pop(session_key, None)
            st.success("API settings cleared from this session.")
            st.rerun()

        return True

    st.markdown(
        """
<div class="ga-api-box">
    <div class="ga-api-title">Enter your OpenAI API settings</div>
    <div class="ga-api-desc">
        This app uses OpenAI only when you run an analysis. Enter your own API key and choose the model.
        The settings are stored only for this browser session and are not saved to GitHub or any app file.
    </div>
</div>
        """.strip(),
        unsafe_allow_html=True,
    )

    user_key = st.text_input(
        "OpenAI API key",
        type="password",
        placeholder="sk-...",
        key="user_openai_api_key_input",
        help="Paste your own OpenAI API key. It is kept only in the current Streamlit session.",
    )

    col_model, col_temp, col_tokens = st.columns([2, 1, 1])

    with col_model:
        selected_model_option = st.selectbox(
            "Model",
            MODEL_OPTIONS,
            index=0,
            key="user_model_select",
            help="Choose the OpenAI model used by all six GenoAnno tools.",
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
            help="Lower values are more deterministic; higher values are more creative.",
        )

    with col_tokens:
        selected_max_tokens = st.number_input(
            "Max tokens",
            min_value=256,
            max_value=16000,
            value=2000,
            step=256,
            key="user_max_tokens_input",
            help="Maximum output length requested from the model.",
        )

    col1, col2 = st.columns([1, 3])

    with col1:
        save_clicked = st.button("Use these settings", key="save_user_api_settings", use_container_width=True)

    with col2:
        st.caption("No password is required. Users provide their own API key and model choice before opening the analysis tools.")

    if save_clicked:
        cleaned_key = user_key.strip()

        if not cleaned_key:
            st.error("Please enter an OpenAI API key.")
            return False

        if not cleaned_key.startswith("sk-"):
            st.warning("This does not look like a standard OpenAI API key. If it is correct, the app will still try to use it.")

        final_model = custom_model.strip() if selected_model_option == "Custom model name" else selected_model_option

        if not final_model:
            st.error("Please choose a model or enter a custom model name.")
            return False

        st.session_state["user_openai_api_key"] = cleaned_key
        st.session_state["user_selected_model"] = final_model
        st.session_state["user_temperature"] = float(selected_temperature)
        st.session_state["user_max_tokens"] = int(selected_max_tokens)

        st.success("API settings saved for this session.")
        st.rerun()

    return False

# ============================================================
# Page sections
# ============================================================
def render_hero(config_ready: bool, child_apps_ready: int) -> None:
    status_text = "Ready" if config_ready else "Needs key"

    st.markdown(
        f"""
<div class="ga-hero">
    <div class="ga-hero-content">
        <div class="ga-badge-row">
            <div class="ga-badge">Genome annotation suite</div>
            <div class="ga-badge">Oral microbiology focused</div>
            <div class="ga-badge">AI-assisted interpretation</div>
        </div>
        <h1>GenoAnno <span>Master Suite</span></h1>
        <p>
            A professional workspace for translating bacterial genome annotations into pathway signals,
            phenotype groups, functional-category summaries, and oral-bacteria similarity insights.
            Select a tool below to open it directly on this page.
        </p>
    </div>
</div>

<div class="ga-stats">
    <div class="ga-stat-card">
        <div class="ga-stat-label">Available tools</div>
        <div class="ga-stat-value">{len(APP_FILES)}</div>
    </div>
    <div class="ga-stat-card">
        <div class="ga-stat-label">Detected app files</div>
        <div class="ga-stat-value">{child_apps_ready}/{len(APP_FILES)}</div>
    </div>
    <div class="ga-stat-card">
        <div class="ga-stat-label">OpenAI setup</div>
        <div class="ga-stat-value">{status_text}</div>
    </div>
    <div class="ga-stat-card">
        <div class="ga-stat-label">Interface</div>
        <div class="ga-stat-value">Tool launcher</div>
    </div>
</div>
        """.strip(),
        unsafe_allow_html=True,
    )


def render_tool_launcher() -> None:
    st.markdown('<div class="ga-section-title">Analysis tools</div>', unsafe_allow_html=True)
    st.markdown(
        '<div class="ga-section-subtitle">Choose a tool, then click Open this tool. No analysis app loads until you click Open.</div>',
        unsafe_allow_html=True,
    )

    app_items = list(APP_FILES.keys())
    for row_start in range(0, len(app_items), 3):
        cols = st.columns(3)
        for col, app_name in zip(cols, app_items[row_start: row_start + 3]):
            with col:
                with st.container(border=True):
                    st.markdown(
                        f"""
<div class="ga-tool-number">{APP_NUMBERS[app_name]}</div>
<div class="ga-tool-title">{app_name}</div>
<div class="ga-tool-desc">{APP_DESCRIPTIONS[app_name]}</div>
                        """.strip(),
                        unsafe_allow_html=True,
                    )
                    if st.button("Open this tool", key=f"open_{app_prefix(app_name)}", use_container_width=True):
                        st.session_state["active_app"] = app_name
                        st.rerun()


def render_selected_app_header(app_name: str) -> None:
    st.markdown(
        f"""
<div class="ga-selected-bar">
    <div class="ga-selected-label">Selected analysis tool</div>
    <div class="ga-selected-title">{app_name}</div>
    <div class="ga-selected-desc">{APP_DESCRIPTIONS[app_name]}</div>
</div>
        """.strip(),
        unsafe_allow_html=True,
    )


def main() -> None:
    inject_css()

    # Website-based API-key setup. No password is required.
    if not api_key_setup():
        st.stop()

    config_ready = config_has_key()
    child_apps_ready = count_available_child_apps()

    if "active_app" not in st.session_state:
        st.session_state["active_app"] = None

    render_hero(config_ready=config_ready, child_apps_ready=child_apps_ready)
    render_tool_launcher()

    active_app = st.session_state.get("active_app")

    if active_app is None:
        st.markdown(
            """
<div class="ga-selected-bar">
    <div class="ga-selected-label">No tool opened</div>
    <div class="ga-selected-title">Select a tool and click Open this tool</div>
    <div class="ga-selected-desc">The six analysis apps stay unloaded until you intentionally open one.</div>
</div>
            """.strip(),
            unsafe_allow_html=True,
        )
        return

    render_selected_app_header(active_app)

    with st.container():
        st.markdown('<div class="ga-app-frame">', unsafe_allow_html=True)
        run_child_app(active_app)
        st.markdown('</div>', unsafe_allow_html=True)


if __name__ == "__main__":
    main()
