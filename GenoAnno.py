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
    --ga-bg: #f5f7fb;
    --ga-surface: #ffffff;
    --ga-surface-soft: #f8fafc;
    --ga-text: #111827;
    --ga-muted: #64748b;
    --ga-border: #e2e8f0;
    --ga-primary: #1d4ed8;
    --ga-primary-dark: #1e3a8a;
    --ga-accent: #0f766e;
    --ga-shadow-sm: 0 6px 18px rgba(15, 23, 42, 0.06);
    --ga-shadow-md: 0 18px 45px rgba(15, 23, 42, 0.09);
}

.stApp {
    background: var(--ga-bg) !important;
    color: var(--ga-text) !important;
}

header[data-testid="stHeader"] {
    background: rgba(245, 247, 251, 0.94) !important;
    backdrop-filter: blur(14px) !important;
    border-bottom: 1px solid rgba(226, 232, 240, 0.95) !important;
}

#MainMenu, footer {
    visibility: hidden !important;
}

.block-container {
    max-width: 1320px !important;
    padding-top: 1.25rem !important;
    padding-left: 2rem !important;
    padding-right: 2rem !important;
    padding-bottom: 4rem !important;
}

/* ---------- Header ---------- */
.ga-topbar {
    display: flex;
    align-items: center;
    justify-content: space-between;
    gap: 1.25rem;
    padding: 1.35rem 1.55rem;
    margin-bottom: 1rem;
    border: 1px solid var(--ga-border);
    border-radius: 22px;
    background: var(--ga-surface);
    box-shadow: var(--ga-shadow-sm);
}

.ga-brand-kicker {
    color: var(--ga-primary);
    font-size: 0.76rem;
    font-weight: 900;
    letter-spacing: 0.12em;
    text-transform: uppercase;
    margin-bottom: 0.2rem;
}

.ga-brand-title {
    margin: 0;
    color: var(--ga-text);
    font-size: 2rem;
    line-height: 1.05;
    font-weight: 950;
    letter-spacing: -0.055em;
}

.ga-brand-subtitle {
    margin-top: 0.45rem;
    color: var(--ga-muted);
    font-size: 0.98rem;
    line-height: 1.55;
    max-width: 820px;
}

.ga-status-pill {
    flex: 0 0 auto;
    display: inline-flex;
    align-items: center;
    padding: 0.5rem 0.8rem;
    border-radius: 999px;
    border: 1px solid #bfdbfe;
    background: #eff6ff;
    color: #1e40af;
    font-size: 0.82rem;
    font-weight: 850;
}

/* ---------- Section headers ---------- */
.ga-section {
    margin-top: 1.2rem;
}

.ga-section-title {
    color: var(--ga-text);
    font-size: 1.35rem;
    font-weight: 950;
    letter-spacing: -0.04em;
    margin-bottom: 0.22rem;
}

.ga-section-subtitle {
    color: var(--ga-muted);
    font-size: 0.96rem;
    line-height: 1.55;
    margin-bottom: 0.9rem;
}

/* ---------- Tool cards ---------- */
div[data-testid="stVerticalBlockBorderWrapper"] {
    border-radius: 18px !important;
    border: 1px solid var(--ga-border) !important;
    background: var(--ga-surface) !important;
    box-shadow: var(--ga-shadow-sm) !important;
}

.ga-tool-number {
    display: inline-flex;
    align-items: center;
    justify-content: center;
    min-width: 2rem;
    height: 2rem;
    padding: 0 0.55rem;
    border-radius: 10px;
    background: #eff6ff;
    color: var(--ga-primary-dark);
    font-size: 0.76rem;
    font-weight: 950;
    letter-spacing: 0.05em;
    margin-bottom: 0.7rem;
}

.ga-tool-title {
    color: var(--ga-text);
    font-size: 1.06rem;
    font-weight: 900;
    line-height: 1.3;
    letter-spacing: -0.025em;
    margin-bottom: 0.4rem;
}

.ga-tool-desc {
    color: var(--ga-muted);
    font-size: 0.91rem;
    line-height: 1.55;
    min-height: 5.7rem;
}

/* ---------- Active tool workspace ---------- */
.ga-workspace-head {
    margin-top: 1.3rem;
    padding: 1.2rem 1.35rem;
    border-radius: 18px;
    background: linear-gradient(135deg, #eff6ff 0%, #f0fdfa 100%);
    border: 1px solid #dbeafe;
}

.ga-workspace-kicker {
    color: var(--ga-primary);
    font-size: 0.75rem;
    font-weight: 900;
    letter-spacing: 0.1em;
    text-transform: uppercase;
}

.ga-workspace-title {
    margin-top: 0.22rem;
    color: var(--ga-text);
    font-size: 1.45rem;
    font-weight: 950;
    letter-spacing: -0.04em;
}

.ga-workspace-desc {
    margin-top: 0.3rem;
    color: var(--ga-muted);
    font-size: 0.95rem;
    line-height: 1.55;
}

.ga-empty-state {
    margin-top: 1.2rem;
    padding: 1.5rem;
    border: 1px dashed #cbd5e1;
    border-radius: 18px;
    background: rgba(255, 255, 255, 0.72);
    text-align: center;
}

.ga-empty-title {
    color: var(--ga-text);
    font-weight: 900;
    font-size: 1.05rem;
}

.ga-empty-desc {
    color: var(--ga-muted);
    margin-top: 0.3rem;
    font-size: 0.92rem;
}

/* ---------- API settings ---------- */
.ga-api-summary {
    padding: 1rem 1.1rem;
    border-radius: 16px;
    border: 1px solid var(--ga-border);
    background: var(--ga-surface);
    box-shadow: var(--ga-shadow-sm);
    margin-bottom: 0.75rem;
}

.ga-api-title {
    color: var(--ga-text);
    font-size: 1rem;
    font-weight: 900;
}

.ga-api-desc {
    color: var(--ga-muted);
    font-size: 0.9rem;
    line-height: 1.5;
    margin-top: 0.25rem;
}

/* ---------- Child app frame ---------- */
.ga-app-frame {
    margin-top: 0.8rem;
    padding: 0.8rem;
    border-radius: 20px;
    border: 1px solid var(--ga-border);
    background: var(--ga-surface);
    box-shadow: var(--ga-shadow-md);
}

/* ---------- Streamlit controls ---------- */
.stButton > button,
.stDownloadButton > button {
    border-radius: 12px !important;
    border: 1px solid transparent !important;
    background: var(--ga-primary) !important;
    color: white !important;
    font-weight: 850 !important;
    box-shadow: none !important;
    min-height: 2.65rem;
}

.stButton > button:hover,
.stDownloadButton > button:hover {
    background: var(--ga-primary-dark) !important;
    border-color: var(--ga-primary-dark) !important;
}

div[data-testid="stTextArea"] textarea,
div[data-testid="stTextInput"] input,
div[data-testid="stNumberInput"] input,
div[data-testid="stSelectbox"] > div > div,
div[data-testid="stFileUploader"] {
    border-radius: 12px !important;
    border-color: var(--ga-border) !important;
}

div[data-testid="stFileUploader"] {
    background: var(--ga-surface-soft) !important;
}

div[data-testid="stExpander"] {
    border: 1px solid var(--ga-border) !important;
    border-radius: 14px !important;
    background: var(--ga-surface) !important;
}

hr {
    border-color: var(--ga-border) !important;
}

@media (max-width: 900px) {
    .ga-topbar {
        display: block;
    }

    .ga-status-pill {
        margin-top: 0.9rem;
    }

    .block-container {
        padding-left: 1rem !important;
        padding-right: 1rem !important;
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
    Clear transient analysis results when moving between tools.

    Child apps currently store generated model text in the shared session key
    "output". Without clearing it, the previous tool's result can appear in the
    next tool before a new analysis is run.
    """
    st.session_state.pop("output", None)

    # Also clear common transient result-like keys if child apps add them later.
    transient_exact_keys = {
        "analysis_output",
        "generated_output",
        "model_output",
        "result",
        "results",
    }
    for key in list(st.session_state.keys()):
        if key in transient_exact_keys:
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
    Settings remain in st.session_state and are never written to app files.
    """

    if st.session_state.get("user_openai_api_key", ""):
        active_model = st.session_state.get("user_selected_model", "gpt-4o-mini")
        active_temperature = st.session_state.get("user_temperature", 0.5)
        active_max_tokens = st.session_state.get("user_max_tokens", 2000)

        st.markdown(
            f"""
<div class="ga-api-summary">
    <div class="ga-api-title">OpenAI connection ready</div>
    <div class="ga-api-desc">
        Model: <b>{active_model}</b> &nbsp;·&nbsp;
        Temperature: <b>{active_temperature}</b> &nbsp;·&nbsp;
        Max tokens: <b>{active_max_tokens}</b>
    </div>
</div>
            """.strip(),
            unsafe_allow_html=True,
        )

        with st.expander("OpenAI settings", expanded=False):
            st.caption(
                "Your API key is stored only in the current Streamlit browser session."
            )

            col_a, col_b = st.columns([1, 3])
            with col_a:
                if st.button(
                    "Clear session settings",
                    key="clear_user_api_settings",
                    use_container_width=True,
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

            with col_b:
                st.caption(
                    "Clear the session before handing the browser to another person."
                )

        return True

    st.markdown(
        """
<div class="ga-api-summary">
    <div class="ga-api-title">Connect OpenAI to begin</div>
    <div class="ga-api-desc">
        Enter your API key and choose the model once. The settings are kept only
        for this browser session.
    </div>
</div>
        """.strip(),
        unsafe_allow_html=True,
    )

    with st.expander("OpenAI settings", expanded=True):
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
            "Connect and continue",
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
    status_text = "OpenAI ready" if config_ready else "OpenAI setup required"

    st.markdown(
        f"""
<div class="ga-topbar">
    <div>
        <div class="ga-brand-kicker">Genome annotation interpretation suite</div>
        <h1 class="ga-brand-title">GenoAnno</h1>
        <div class="ga-brand-subtitle">
            A focused workspace for pathway, protein-family, gene-function,
            phenotype, and oral-bacteria interpretation.
            {child_apps_ready}/{len(APP_FILES)} analysis tools detected.
        </div>
    </div>
    <div class="ga-status-pill">{status_text}</div>
</div>
        """.strip(),
        unsafe_allow_html=True,
    )


def render_tool_launcher() -> None:
    st.markdown('<div class="ga-section">', unsafe_allow_html=True)
    st.markdown(
        '<div class="ga-section-title">Choose an analysis tool</div>',
        unsafe_allow_html=True,
    )
    st.markdown(
        '<div class="ga-section-subtitle">Each tool opens in the workspace below. Switching tools clears the previous generated output so results never carry over between analyses.</div>',
        unsafe_allow_html=True,
    )

    app_items = list(APP_FILES.keys())

    for row_start in range(0, len(app_items), 3):
        cols = st.columns(3, gap="medium")

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

                    if st.button(
                        "Open tool",
                        key=f"open_{app_prefix(app_name)}",
                        use_container_width=True,
                    ):
                        switch_active_app(app_name)
                        st.rerun()

    st.markdown('</div>', unsafe_allow_html=True)


def render_selected_app_header(app_name: str) -> None:
    st.markdown(
        f"""
<div class="ga-workspace-head">
    <div class="ga-workspace-kicker">Active workspace</div>
    <div class="ga-workspace-title">{app_name}</div>
    <div class="ga-workspace-desc">{APP_DESCRIPTIONS[app_name]}</div>
</div>
        """.strip(),
        unsafe_allow_html=True,
    )

    col_back, col_clear, col_space = st.columns([1, 1, 5])

    with col_back:
        if st.button(
            "Close tool",
            key="close_active_tool",
            use_container_width=True,
        ):
            clear_child_transient_state(previous_app=app_name)
            st.session_state["active_app"] = None
            st.rerun()

    with col_clear:
        if st.button(
            "Clear output",
            key="clear_active_output",
            use_container_width=True,
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

    render_hero(config_ready=config_ready, child_apps_ready=child_apps_ready)
    render_tool_launcher()

    active_app = st.session_state.get("active_app")

    if active_app is not None and active_app not in APP_FILES:
        clear_child_transient_state(previous_app=active_app)
        st.session_state["active_app"] = None
        active_app = None

    if active_app is None:
        st.session_state.pop("output", None)

        st.markdown(
            """
<div class="ga-empty-state">
    <div class="ga-empty-title">No analysis tool is open</div>
    <div class="ga-empty-desc">
        Choose one of the six tools above. Only the selected tool is loaded,
        keeping the workspace clean and preventing outputs from other tools
        from appearing here.
    </div>
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
