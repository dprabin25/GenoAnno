# ============================================================
# GenoAnno Master Streamlit App - Password Protected Tool Dashboard
# Runs six separate Streamlit apps from one professional home screen.
#
# How to run:
#   streamlit run master_app.py
# ============================================================

from __future__ import annotations

import hashlib
import hmac
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
    page_title="GenoAnno Analyzer | Master Suite",
    layout="wide",
    initial_sidebar_state="collapsed",
)


BASE_DIR = Path(__file__).resolve().parent
APPS_DIR = BASE_DIR / "apps"
CONFIG_FILE = BASE_DIR / "config.txt"


APP_FILES: Dict[str, str] = {
    "Pathway Expectation": "web1_pathway_expectation.py",
    "Phenotype Grouping": "web2_phenotype_grouping.py",
    "KEGG Pathway Count": "web3_kegg_pathway_count.py",
    "Protein Family Count": "web4_protein_family_count.py",
    "Bakta Functional Category": "web5_bakta_functional_category.py",
    "Similar Oral Bacteria": "web6_similar_oral_bacteria.py",
}


APP_DESCRIPTIONS: Dict[str, str] = {
    "Pathway Expectation": "Predict phenotype-relevant active pathways from known bacterial behavior and pathway lists.",
    "Phenotype Grouping": "Filter annotated pathway tables and organize retained pathways into biological phenotype groups.",
    "KEGG Pathway Count": "Count repeated KEGG pathway terms and interpret dominant functional pathway signals.",
    "Protein Family Count": "Summarize repeated protein-family terms into phenotype-level functional interpretations.",
    "Bakta Functional Category": "Parse Bakta product annotations, count functional categories, and generate phenotype summaries.",
    "Similar Oral Bacteria": "Compare functional-category patterns against well-characterized oral bacteria.",
}


APP_NUMBERS: Dict[str, str] = {
    "Pathway Expectation": "01",
    "Phenotype Grouping": "02",
    "KEGG Pathway Count": "03",
    "Protein Family Count": "04",
    "Bakta Functional Category": "05",
    "Similar Oral Bacteria": "06",
}


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

.ga-app-frame {
    background: rgba(255, 255, 255, 0.72);
    border: 1px solid rgba(148, 163, 184, 0.24);
    border-radius: 28px;
    padding: 1rem;
    box-shadow: var(--ga-soft-shadow);
    margin-top: 0.7rem;
}

.ga-login-box {
    max-width: 560px;
    margin: 4rem auto 0 auto;
    padding: 2rem;
    border-radius: 28px;
    background: rgba(255, 255, 255, 0.86);
    border: 1px solid rgba(148, 163, 184, 0.32);
    box-shadow: var(--ga-shadow);
}

.ga-login-title {
    font-size: 2rem;
    font-weight: 950;
    letter-spacing: -0.05em;
    color: var(--ga-text);
    margin-bottom: 0.4rem;
}

.ga-login-desc {
    color: var(--ga-muted);
    font-size: 1rem;
    line-height: 1.6;
    margin-bottom: 1rem;
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
# Password protection
# ============================================================
def check_password() -> bool:
    """
    Password gate for the public Streamlit app.

    Required Streamlit Secret:
        APP_PASSWORD = "your_private_password_here"
    """

    app_password = st.secrets.get("APP_PASSWORD", "")

    if not app_password:
        st.markdown(
            """
<div class="ga-login-box">
    <div class="ga-login-title">GenoAnno is protected</div>
    <div class="ga-login-desc">
        APP_PASSWORD is not configured yet. Add APP_PASSWORD in Streamlit Secrets,
        then reboot the app.
    </div>
</div>
            """.strip(),
            unsafe_allow_html=True,
        )
        st.error("Missing APP_PASSWORD in Streamlit Secrets.")
        return False

    if st.session_state.get("password_correct", False):
        return True

    st.markdown(
        """
<div class="ga-login-box">
    <div class="ga-login-title">Protected GenoAnno App</div>
    <div class="ga-login-desc">
        Enter the private app password to access the analysis tools.
    </div>
</div>
        """.strip(),
        unsafe_allow_html=True,
    )

    password = st.text_input(
        "App password",
        type="password",
        key="app_password_input",
    )

    if st.button("Unlock app", key="unlock_app_button"):
        if hmac.compare_digest(password, app_password):
            st.session_state["password_correct"] = True
            st.rerun()
        else:
            st.error("Incorrect password.")

    return False


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
    The uploaded child apps read config.txt from their own folder.
    This copies the master config.txt into apps/config.txt when available.
    """
    child_config = APPS_DIR / "config.txt"

    if CONFIG_FILE.exists():
        child_config.write_text(CONFIG_FILE.read_text(encoding="utf-8"), encoding="utf-8")
    elif not child_config.exists():
        st.warning(
            "config.txt was not found beside master_app.py. Create one with KEY=your_openai_api_key before running OpenAI-based analyses."
        )


def config_has_key() -> bool:
    """
    Checks whether OpenAI credentials are available.

    Priority:
    1. Streamlit Secrets KEY
    2. Local config.txt KEY
    """

    try:
        secret_key = st.secrets.get("KEY", "")
        if secret_key:
            return True
    except Exception:
        pass

    if not CONFIG_FILE.exists():
        return False

    text = CONFIG_FILE.read_text(encoding="utf-8", errors="ignore")
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if line.startswith("KEY="):
            value = line.split("=", 1)[1].strip()
            return bool(value and "your_openai_api_key" not in value.lower())

    return False


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

                    if st.button(
                        "Open this tool",
                        key=f"open_{app_prefix(app_name)}",
                        use_container_width=True,
                    ):
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


# ============================================================
# Main app
# ============================================================
def main() -> None:
    inject_css()

    # Password gate: nothing below this line is shown unless password is correct.
    if not check_password():
        st.stop()

    ensure_config_available()

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
