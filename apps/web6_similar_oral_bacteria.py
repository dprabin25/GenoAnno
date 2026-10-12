# ============================================================
# web6_similar_oral_bacteria.py
#
# GenoAnno tab: find oral bacteria with similar gene composition by
# combining FOUR processed summary tables derived from three raw
# annotation files:
#
#   Raw file                           -> Processed table(s)
#   -------------------------------------------------------------
#   Proksee/Bakta annotation (.tsv)    -> Table 1: Functional category gene counts
#   KBASE annotation (.tsv)            -> Table 2: KEGG pathway gene counts
#                                      -> Table 3: Protein function family gene counts
#   Products.tsv (KEGG-Decoder style)  -> Table 4: Metabolic pathway completeness
#
# This file is run by master_app.py via runpy.run_path(..., run_name=
# f"__genoanno_{prefix}__"), NOT as __main__. That means any code guarded
# by `if __name__ == "__main__":` never executes when opened from the
# dashboard. All UI code below therefore runs at module level.
# ============================================================

import functools
import io
import re

import pandas as pd
import streamlit as st

# ---------------------------------------------------------------------------
# Combined prompt template
# ---------------------------------------------------------------------------

DEFAULT_INPUT_DESCRIPTION = """I annotated genes from a bacterial genome (an oral bacterium) using multiple annotation pipelines and summarized the results into four processed tables describing its functional gene composition.

Table 1 - Functional category gene counts (Bakta):
I selected the Product column from the Bakta annotation and counted how many genes map to each repeated product/function term. Ultra-generic, non-informative labels ("hypothetical protein", "uncharacterized protein") have been removed. Every remaining category is listed, ordered by gene count, including categories represented by only one gene - most of this table is single-copy categories, and a low count is not a reason to disregard a row.
Column 1: Functional category
Column 2: Gene count

Table 2 - KEGG-annotated gene function descriptions (KBASE):
Each row is the free-text function description ("kegg_hit") tied to a gene's best KO match, counted by how many genes share that exact description. These are gene-level descriptions, NOT a curated pathway hierarchy - there is no row literally named "glycolysis" or "TCA cycle" unless it appears verbatim below. Do not invent or assume pathway-level labels that are not literal rows in this table. Generic/uninformative rows have been removed; every remaining entry is listed, ordered by gene count, including single-gene entries.
Column 1: KEGG-hit description
Column 2: Gene count

Table 3 - Protein function family gene counts (KBASE):
I assigned each gene to its corresponding protein function family (peptidase family) and counted the number of genes mapped to each family.
Column 1: Protein function family
Column 2: Gene count

Table 4 - Metabolic pathway completeness (KEGG-Decoder style):
For each metabolic pathway, I calculated its completeness and, where applicable, determined whether the pathway is present (TRUE/FALSE).
Column 1: Metabolic pathway name
Column 2: Pathway status (completeness score or TRUE/FALSE)
Filtering rule already applied: pathways with 0, FALSE, empty, missing, NA, or NaN values were removed before analysis - so Table 4 below ONLY lists pathways that are actually present/complete. A separate list of pathway GROUPS that are entirely absent (zero present rows in that group) is also provided below Table 4 - treat total absence of a group as meaningful negative evidence about this organism's capabilities, not as missing data. For example, if every "CAZy:" row is absent, this genome shows no detectable capacity to degrade the corresponding dietary polysaccharides, which is a strong asaccharolytic signal - do not describe such an organism as "broadly saccharolytic."

Only these four processed tables (not the raw annotation files) should be used for phenotype grouping and comparison.

GROUNDING RULE (important): Every category, pathway, or family you cite as "shared" or as evidence must be a value that appears verbatim as a row in Tables 1-4 below. Do not cite a category unless you can point to its literal row. Do not describe capabilities (e.g., "broadly saccharolytic," "diverse CAZy families") unless multiple distinct, clearly-named carbohydrate-degradation rows are actually present in Table 4 - a handful of generic "glycosyltransferase" or "glycosyl transferase" entries in Tables 1/3 usually reflect structural cell-envelope biosynthesis (peptidoglycan, lipopolysaccharide, lipid A), not dietary carbohydrate degradation, unless the label explicitly names an extracellular substrate (e.g., cellulase, amylase, xylanase, chitinase, pectinase).

DETERMINISM RULE (important): This is a grounded lookup/classification task, not creative writing. Given the same four tables, your answer should be reproducible - do not vary your conclusion between runs on identical data, and do not introduce candidate bacteria that aren't anchored to STEP 0's checklist below. If two candidates are close, say so explicitly rather than picking one arbitrarily.

STEP 0 - Before any narrative reasoning, extract a fixed evidence checklist as a short bullet list, pulling only discriminating (non-universal) rows verbatim from Tables 1-4: the specific CAZy/carbohydrate categories present (if any) and the fact of total CAZy absence (if that's the case), the specific SCFA/fermentation end-product pathways present, the specific respiratory complexes present vs. entirely absent, and any unusual protein families (Table 3) or KEGG-hit descriptions (Table 2) that are not generic housekeeping. This checklist is the sole basis for every later step - do not add or reason from anything outside it.

STEP 1 - Using only the STEP 0 checklist, derive an internal phenotype profile of this genome. Explicitly reason about:
- Respiratory/fermentative capacity: which respiratory chain complexes (Table 1/Table 4) are present or absent, and whether the profile looks aerobic, microaerophilic, anaerobic, or strictly fermentative.
- Carbon source usage: which CAZy/carbohydrate-active categories (Table 1) are present vs. largely absent - i.e., whether this organism looks broadly saccharolytic (degrades many sugars) or asaccharolytic/limited (relies on other substrates such as amino acids, lactate, or organic acids).
- Fermentation end products: which SCFA/alcohol conversion pathways (Table 4) are complete (e.g., lactate, propionate, acetate, butyrate). End-product profile is one of the strongest, most specific phenotype signals for oral bacteria and often more diagnostic than raw gene overlap.
- Distinctive protein families (Table 3) and KEGG pathways (Table 2) beyond universal housekeeping genes.

STEP 2 - When weighing evidence across all four tables, explicitly downweight categories that are near-universal across nearly all bacteria and carry little taxonomic signal (e.g., "hypothetical protein", generic ABC transporters, ribosomal proteins, DNA replication/repair machinery, chaperones, generic transport systems). These will dominate the raw gene counts but should NOT drive the match. Instead prioritize categories that vary meaningfully between oral genera: CAZy/carbohydrate degradation breadth, SCFA/fermentation end products, respiratory chain complex composition, nitrogen/sulfur metabolism, and any pathway present in some oral taxa but absent in most others.

STEP 3 - Only after deriving this phenotype profile, identify well-characterized oral bacteria whose published phenotype (fermentation strategy, carbon source range, oxygen tolerance, cell/colony morphology) is genuinely consistent with the derived profile - not simply bacteria that happen to share the largest raw count of generic annotated genes. If a well-known genus overlaps only on universal housekeeping categories, exclude it even if its raw overlap count is high.

Question: Identify other oral bacteria that have similar gene composition to mine. Use all four tables together - a strong match should show consistent overlap across functional categories (Table 1), KEGG pathways (Table 2), protein function families (Table 3), and completed metabolic pathways (Table 4), weighted as described above, not just raw agreement in one table.

Output, in this order:
1. The STEP 0 evidence checklist (bullet list, verbatim rows only).
2. A short (2-4 sentence) summary of the phenotype profile you derived in Step 1, so the reasoning can be sanity-checked.
3. A table with columns: Bacterium name | Phylum | Shared functional categories, pathways, and families (cite which of Tables 1-4 support the match, and note whether each is a discriminating or universal category) | Phenotype description | Confidence.

For the Phylum column, give the candidate's phylum, using the current name with the older name in brackets where they differ - for example "Bacteroidota (Bacteroidetes)" or "Bacillota (Firmicutes)".
4. A short summary paragraph.

List well-characterized bacteria with the most consistent shared, discriminating (not universal) functional categories/pathways/families across the four tables. Describe the phenotype of each of these bacteria.

Rank the species and provide your reasoning for each rank. Sort by confidence, with the highest-confidence match listed first.
"""

# ---------------------------------------------------------------------------
# Parsers: raw file -> processed table
# ---------------------------------------------------------------------------

# Labels that appear in nearly every bacterial genome and carry ~no
# taxonomic signal. These get dropped from Tables 1 and 2 so the model's
# attention isn't spent on rows that can't discriminate between species.
GENERIC_LABELS = {
    "hypothetical protein",
    "uncharacterized protein",
}

# Row cap for Tables 1 and 2. None = send every category.
#
# Previously 40. Ranking by gene count selects for gene-family expansion,
# which is the least species-specific part of a genome, while the categories
# that distinguish related species are overwhelmingly single-copy: 86% of
# distinct CDS products in the Tannerella forsythia input occur exactly once,
# and because ties break alphabetically those singletons are scattered from
# rank ~220 to 1526. Any intermediate cap is therefore an arbitrary slice
# through alphabetically-ordered ties - at 500 rows the S-layer subunits are
# included but the secreted protease is not; at 1000 the protease is included
# but the sialidase is not. Sending every row costs ~0.3 cents more per run at
# gpt-4o-mini list prices (3,268 -> 26,858 input tokens, well inside a 128k
# context) and needs no cutoff to be justified in a methods section.
TOP_N_ROWS = None


def _clean_and_cap(counts: pd.DataFrame, label_col: str, top_n=TOP_N_ROWS):
    """Drop generic/uninformative labels, sort by count, keep the top N.

    top_n=None keeps every row. Returns (trimmed_df, note) where note states
    what the table does and does not contain, so the prompt stays accurate
    whichever setting is in force.
    """
    total_rows = len(counts)
    filtered = counts[~counts[label_col].str.strip().str.lower().isin(GENERIC_LABELS)]
    dropped_generic = total_rows - len(filtered)
    filtered = filtered.sort_values("Gene count", ascending=False).reset_index(drop=True)
    trimmed = filtered if top_n is None else filtered.head(top_n)
    omitted = len(filtered) - len(trimmed)

    if omitted > 0:
        note = (
            f"(showing top {len(trimmed)} of {total_rows} total categories by gene count; "
            f"{dropped_generic} generic/non-discriminating rows removed; "
            f"{omitted} additional lower-count categories omitted for brevity)"
        )
    else:
        note = (
            f"(all {len(trimmed)} categories listed, ordered by gene count, "
            f"including single-copy categories; "
            f"{dropped_generic} generic/non-discriminating rows removed)"
        )

    return trimmed, note


def parse_bakta_functional_categories(bakta_file):
    """Proksee/Bakta annotation TSV -> Table 1 (Functional category, Gene count)."""
    raw = bakta_file.read()
    if isinstance(raw, bytes):
        raw = raw.decode("utf-8", errors="replace")
    lines = raw.splitlines()
    header_idx = next(i for i, l in enumerate(lines) if l.startswith("#Sequence Id"))
    data_str = "\n".join(lines[header_idx:])
    df = pd.read_csv(io.StringIO(data_str), sep="\t")
    df.columns = [c.lstrip("#").strip() for c in df.columns]

    # Keep only protein-coding features. A Bakta TSV annotates every feature
    # type -- tRNA, rRNA, ncRNA, ncRNA-region (riboswitches), tmRNA, CRISPR --
    # and each carries a Product string, so they are otherwise counted as
    # "genes" in a table the prompt presents as protein functional categories.
    # Multi-copy RNA families can then outrank real protein products and claim
    # top-N slots: in the Tannerella forsythia input, "Acido-Lenti-1 RNA" (12)
    # and "Cobalamin riboswitch" (6) both reach Table 1 untreated. Guarded so
    # a table lacking a Type column behaves exactly as before.
    if "Type" in df.columns:
        df = df[df["Type"].astype(str).str.strip().str.lower() == "cds"]

    products = df["Product"].dropna()
    products = products[products.str.strip() != ""]
    counts = (
        products.value_counts()
        .rename_axis("Functional category")
        .reset_index(name="Gene count")
    )
    return _clean_and_cap(counts, "Functional category")


def parse_kbase_kegg_pathways(kbase_file):
    """KBASE annotation TSV -> Table 2 (KEGG-hit description, Gene count)."""
    df = pd.read_csv(kbase_file, sep="\t")
    hits = df["kegg_hit"].dropna()
    hits = hits[hits.str.strip() != ""]
    counts = (
        hits.value_counts()
        .rename_axis("KEGG-hit description")
        .reset_index(name="Gene count")
    )
    return _clean_and_cap(counts, "KEGG-hit description")


def parse_kbase_protein_families(kbase_file) -> pd.DataFrame:
    """KBASE annotation TSV -> Table 3 (Protein function family, Gene count)."""
    df = pd.read_csv(kbase_file, sep="\t")
    fams = df["peptidase_family"].dropna()
    fams = fams[fams.str.strip() != ""]
    counts = (
        fams.value_counts()
        .rename_axis("Protein function family")
        .reset_index(name="Gene count")
    )
    return counts


def parse_products_completeness(products_file):
    """Products.tsv (one row per genome, one column per pathway) -> Table 4
    (Metabolic pathway name, Pathway status), FALSE/0/empty/NA rows removed.

    Also returns a list of pathway GROUPS (text before the first ":") that
    have ZERO surviving rows - i.e. every member of that group was filtered
    out as absent. This is the key fix for the "broadly saccharolytic"
    misread: when every "CAZy: ..." row is absent, that's strong negative
    evidence, not silence, and the model needs it spelled out explicitly
    rather than inferring it from a row that simply isn't there.
    """
    df = pd.read_csv(products_file, sep="\t")
    row = df.iloc[0]
    row = row.drop(labels=["genome"], errors="ignore")
    long = row.rename_axis("Metabolic pathway name").reset_index(name="Pathway status")

    def is_present(v):
        if pd.isna(v):
            return False
        if isinstance(v, str):
            v_clean = v.strip().lower()
            if v_clean in ("", "na", "nan", "false", "0"):
                return False
            return True
        if isinstance(v, bool):
            return v
        try:
            return float(v) != 0
        except (TypeError, ValueError):
            return bool(v)

    def group_of(name: str) -> str:
        return name.split(":")[0].strip() if ":" in name else name

    long["Group"] = long["Metabolic pathway name"].apply(group_of)
    present_mask = long["Pathway status"].apply(is_present)
    kept = long[present_mask].reset_index(drop=True)

    all_groups = long["Group"].unique().tolist()
    present_groups = set(kept["Group"].unique().tolist())
    absent_groups = [g for g in all_groups if g not in present_groups]

    kept = kept.drop(columns=["Group"])
    return kept, absent_groups


# ---------------------------------------------------------------------------
# Prompt assembly
# ---------------------------------------------------------------------------

# The model writes multi-item evidence cells in its comparison table as
# "- item A <br> - item B", because an HTML <br> is the usual way to force a
# line break inside a markdown table cell. st.markdown() does not render raw
# HTML unless unsafe_allow_html=True, so those tags show up as literal "<br>"
# text in the Result panel. Rather than enabling raw HTML -- which would let
# arbitrary model output inject markup into the page -- normalise the tags
# into plain "; " separators, matching how the other GenoAnno tabs clean their
# output. Handles the bullet prefixes the model pairs them with so cells read
# as "item A; item B" and not "- item A; - item B".

_BR_PATTERN = re.compile(r"\s*<\s*br\s*/?\s*>\s*", re.IGNORECASE)


def clean_model_output(text):
    if not text:
        return ""

    cleaned = _BR_PATTERN.sub("; ", str(text))

    # Drop the bullet markers that followed each <br>, now stranded after a
    # separator ("; - item" -> "; item").
    cleaned = re.sub(r";\s*[-*•]\s+", "; ", cleaned)

    # A cell that began with a bullet keeps it only if it is a real list item
    # at line start; inside a table cell it is noise ("| - item" -> "| item").
    cleaned = re.sub(r"\|\s*[-*•]\s+", "| ", cleaned)

    # Tidy artefacts from the substitutions.
    cleaned = re.sub(r"(;\s*){2,}", "; ", cleaned)
    cleaned = re.sub(r";\s*\|", " |", cleaned)
    cleaned = re.sub(r"[ \t]{2,}", " ", cleaned)
    cleaned = re.sub(r";\s*$", "", cleaned, flags=re.MULTILINE)

    return cleaned.strip()


# The prompt asks the model to write out a STEP 0 evidence checklist and calls
# it "the sole basis for every later step". That checklist is the grounding
# mechanism, not decoration: a non-reasoning model has no hidden scratchpad, so
# the written rows are what hold the later steps to verbatim table values.
# It is therefore still requested and still generated -- it is only hidden from
# the Result panel, which readers want to start at the phenotype summary. The
# full response, checklist included, is what the download button writes, so the
# grounding stays auditable.

_STEP0_PATTERN = re.compile(
    r"""(?:^|\n)            # start of a line
        [^\S\n]*            # leading spaces
        (?:[#>*_\-\d.)\s]*) # markdown heading / bold / list decoration
        STEP\s*0\b          # the heading itself
        .*?                 # the checklist body
        (?=\n[^\S\n]*(?:[#>*_\-\d.)\s]*)STEP\s*1\b)  # stop before STEP 1
    """,
    re.IGNORECASE | re.DOTALL | re.VERBOSE,
)


def strip_step0_section(text):
    """Remove the STEP 0 checklist from text shown in the Result panel.

    Only strips when a STEP 1 heading follows, so a response that is shaped
    differently than expected is passed through untouched rather than being
    truncated.
    """
    if not text:
        return ""

    return _STEP0_PATTERN.sub("\n", str(text)).strip()


# "STEP n" is prompt scaffolding -- it tells the model what order to reason in
# and is meaningless to someone reading the finished result. The model is also
# inconsistent about it, labelling some sections ("STEP 2 Bacterial Comparison
# Table") and not others ("Comparison Table"), which makes the output look
# half-numbered. Strip the prefix and keep the descriptive title, so headings
# read uniformly however the model chose to label them. Markdown heading and
# bold markers are preserved; a step heading with no title of its own is given
# a sensible default rather than being left blank.

_STEP_LABEL_PATTERN = re.compile(
    r"""^(?P<lead>[^\S\n]*(?:\#{1,6}\s*)?(?:\*\*|__)?)   # heading / bold open
        STEP\s*\d+\b                                      # the step label
        [^\S\n]*(?:[:.–—-]+[^\S\n]*)?           # separator, if any
        (?P<title>.*?)                                    # descriptive title
        (?P<trail>(?:\*\*|__)?[^\S\n]*)$                  # bold close
    """,
    re.IGNORECASE | re.MULTILINE | re.VERBOSE,
)

_STEP_FALLBACK_TITLES = {
    "1": "Phenotype profile summary",
    "2": "Comparison table",
    "3": "Summary",
}


def strip_step_labels(text):
    """Drop "STEP n" prefixes from headings, keeping each section's title."""
    if not text:
        return ""

    def replace(match):
        title = match.group("title").strip()

        if not title:
            step_number = re.search(r"STEP\s*(\d+)", match.group(0), re.IGNORECASE)
            key = step_number.group(1) if step_number else ""
            title = _STEP_FALLBACK_TITLES.get(key, "Result")

        return match.group("lead") + title + match.group("trail")

    return _STEP_LABEL_PATTERN.sub(replace, str(text))


# ---------------------------------------------------------------------------
# Candidate table parsing
# ---------------------------------------------------------------------------

# A binomial: capitalised genus, lower-case species epithet. Accepts an
# abbreviated genus ("P. intermedia") so those rows are not silently lost.
_BINOMIAL = re.compile(r"^(?:[A-Z][a-z]+|[A-Z]\.)\s+[a-z][a-z-]{2,}")

_NOT_A_SPECIES = {
    "bacterium name", "bacterium", "name", "species", "candidate",
    "summary", "summary paragraph", "phylum", "confidence",
}


# ---------------------------------------------------------------------------
# Venn diagram of the evidence each candidate shares with the query genome
# ---------------------------------------------------------------------------
#
# Each row of the candidate table cites the categories, pathways and families
# it shares with the query genome. Those cells are free text written by the
# model, so they are not compared as raw strings: every cited item is matched
# back to its literal row in Tables 1-4 and that row label becomes the set
# element. Two candidates citing "Lactate dehydrogenase (Table 4,
# discriminating)" and "lactate dehydrogenase" therefore count as one shared
# element, and the source table of each element is known.
#
# Items that match no row are kept as written: every category the candidate
# table cites is shown.
#
# Two views of the same items, side by side:
#   heatmap - WHICH categories each candidate shares, every one named, any
#             number of candidates.
#   Venn    - HOW MANY categories each combination of candidates shares or
#             holds alone, as counts, for up to five candidates.

# Brown, blue and orange/gold are reserved in this project's figures for the
# query bacterium, immune cells and proteins; greens and red carry audit
# meaning. The set colours avoid all of them.
#
# VENN_COLORS carry identity (fills, filled markers), in this fixed order:
# deep purple, teal, light rose, olive, lavender. Validated as a categorical
# palette over all pairs: lightness band, chroma floor (each >= 0.1 OKLCH, so
# none reads as grey), CVD separation (worst dE 8.1, protan, rose vs teal) and
# normal-vision separation (worst dE 16.3, lavender vs rose) all pass. Teal
# and rose stay apart under protanopia because they differ in LIGHTNESS; a
# darker rose collapses onto the teal (dE 2.1).
#
# VENN_STROKES draw edges only (outlines, marker rims). They equal the fills
# except the rose, whose fill is below 3:1 against white (2.34); its stroke is
# a darker step of the same hue (3.45). Edges never carry identity -- the
# adjacent species titles do.
#
# Species names are set in ink, never in the set colour.
VENN_COLORS = ("#5d4a9e", "#048a9c", "#e88ab4", "#8a8f2e", "#9e6fd6")
VENN_STROKES = ("#5d4a9e", "#048a9c", "#d0679b", "#8a8f2e", "#9e6fd6")

# Secondary encoding. Under protanopia the teal and rose fills stay apart on
# large areas, but small markers converge. Each candidate therefore also has
# a shape (circle, square, triangle, diamond, inverted triangle), drawn beside
# its species title, so identity never rests on colour alone.
VENN_MARKERS = ("o", "s", "^", "D", "v")

# Neutral inks shared by every evidence figure.
INK = "#222222"
INK_MUTED = "#555555"
MARK_DARK = "#3d3d3d"
MARK_EMPTY = "#ececec"          # absent heatmap cell (present vs absent 9.2:1)
RULE = "#b8b8b8"

FIGURE_STYLE = {
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8,
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
}


TABLE_LABELS = {
    1: "Table 1 - Functional categories",
    2: "Table 2 - KEGG-hit descriptions",
    3: "Table 3 - Protein function families",
    4: "Table 4 - Metabolic pathways",
}

# Shortest row label allowed to match inside a longer cited item. Shorter
# labels ("Sec", "ATPase") would match inside unrelated text.
_MIN_CONTAINED_LABEL = 8

_TABLE_CITATION = re.compile(
    r"\(?\s*\bTables?\s*\d(?:\s*(?:[-–,/&+]|and)\s*\d)*\s*\)?\s*[:\-–—]?\s*",
    re.IGNORECASE,
)
_QUALIFIER_PAREN = re.compile(
    r"\(\s*(?:discriminating|universal|near[- ]universal|housekeeping|generic)\b[^)]*\)",
    re.IGNORECASE,
)
_QUALIFIER_TAIL = re.compile(
    r"[\s,;:\-–—]*\b(?:discriminating|universal|near[- ]universal)\b\s*$",
    re.IGNORECASE,
)


def _norm_label(text):
    text = re.sub(r"[*_`\"“”'‘’]", "", str(text))
    text = " ".join(text.split()).lower()
    return text.strip(" .,;:")


def collect_table_rows(functional_category_df, kegg_pathway_df,
                       protein_family_df, pathway_completeness_df):
    """{normalised label: (verbatim label, table number)} for Tables 1-4.

    Table 4 is indexed first, then 3, 2, 1, so a label that occurs in more
    than one table is attributed to the pathway-level table.
    """
    rows = {}
    sources = (
        (4, pathway_completeness_df, "Metabolic pathway name"),
        (3, protein_family_df, "Protein function family"),
        (2, kegg_pathway_df, "KEGG-hit description"),
        (1, functional_category_df, "Functional category"),
    )

    for table_number, frame, column in sources:
        if frame is None or column not in frame:
            continue

        for label in frame[column].dropna().astype(str):
            key = _norm_label(label)

            if key and key not in rows:
                rows[key] = (label.strip(), table_number)

    return rows


def _ground_item(item, table_rows):
    """The table row an item refers to, as (label, table), or None."""
    key = _norm_label(item)

    if not key:
        return None
    if key in table_rows:
        return table_rows[key]

    best = None

    for row_key, value in table_rows.items():
        if len(row_key) < _MIN_CONTAINED_LABEL or row_key not in key:
            continue
        if best is not None and len(row_key) <= len(best[0]):
            continue
        if re.search(r"(?<!\w)" + re.escape(row_key) + r"(?!\w)", key):
            best = (row_key, value)

    return best[1] if best else None


def _drop_unbalanced_parens(text):
    while text.endswith(")") and text.count(")") > text.count("("):
        text = text[:-1].rstrip()
    while text.startswith("(") and text.count("(") > text.count(")"):
        text = text[1:].lstrip()
    return text


# Words the model attaches to an item to qualify it, and placeholders. They
# are never items themselves, however a cell is formatted.
_NOT_AN_ITEM = {
    "none", "n/a", "na", "-", "discriminating", "universal",
    "near-universal", "near universal", "housekeeping", "generic",
    "non-discriminating", "discriminating category", "universal category",
}


def split_evidence_cell(cell):
    """Cited items from one evidence cell, with table citations and
    discriminating/universal qualifiers removed."""
    text = re.sub(r"[*_`]", "", str(cell))

    # "A (Table 1), B (Table 4)": a citation followed by a comma ends an
    # item, so treat it as a separator. Commas inside a name are untouched.
    text = re.sub(r"(\(\s*Tables?\s*\d[^)]*\))\s*,", r"\1;", text, flags=re.IGNORECASE)

    # Remove every parenthesised citation whole, qualifiers included:
    # "(Table 1, discriminating)" must not leave ", discriminating" behind
    # to be split off as an item of its own.
    text = re.sub(r"\(\s*Tables?\b[^)]*\)", " ", text, flags=re.IGNORECASE)
    items = []

    for piece in re.split(r"\s*(?:;|\n|•)\s*", text):
        piece = _TABLE_CITATION.sub(" ", piece)
        piece = _QUALIFIER_PAREN.sub(" ", piece)
        piece = _QUALIFIER_TAIL.sub("", piece)
        piece = " ".join(piece.split()).strip(" .,:-–—\"'")
        piece = _drop_unbalanced_parens(piece).strip(" .,:-–—\"'")

        if piece and _norm_label(piece) not in _NOT_AN_ITEM:
            items.append(piece)

    return items


def ground_evidence(items, table_rows):
    """{normalised label: (label, table or None)} for a candidate's items.

    Every cited item is kept: matched items carry their verbatim row label and
    table, unmatched ones are kept as written with table None, so nothing the
    model cited disappears silently.

    The model lists several items behind one citation, separated by commas
    ("Glycolysis, Pentose phosphate pathway (Table 4)"), so an item that is
    not a row as a whole is split on commas. The parts are rebuilt greedily
    from the left, joining adjacent parts whenever the join is itself a row,
    so a table row whose name contains a comma ("Mixed acid: Ethanol, Acetate
    to Acetylaldehyde") is never broken up.
    """
    grounded = {}

    def add(hit):
        grounded.setdefault(_norm_label(hit[0]), hit)

    for item in items:
        exact = table_rows.get(_norm_label(item))

        if exact:
            add(exact)
            continue

        parts = [p.strip() for p in item.split(",") if p.strip()]
        i = 0

        while i < len(parts):
            for j in range(len(parts), i, -1):
                joined = table_rows.get(_norm_label(", ".join(parts[i:j])))
                if joined:
                    add(joined)
                    i = j
                    break
            else:
                if _norm_label(parts[i]) not in _NOT_AN_ITEM:
                    add(_ground_item(parts[i], table_rows) or (parts[i], None))
                i += 1

    return grounded


def collect_candidate_evidence(text, table_rows, limit=12):
    """{species: grounded evidence} from the model's candidate table.

    The evidence column is found by its header ("Shared ..."); without a
    recognisable header it is taken to be the third column, as the prompt
    specifies.
    """
    evidence = {}
    evidence_index = 2

    for line in str(text or "").splitlines():
        stripped = line.strip()

        if stripped.count("|") >= 3:
            cells = [c.strip() for c in stripped.strip("|").split("|")]
        elif "\t" in stripped:
            cells = [c.strip() for c in stripped.split("\t")]
        else:
            continue

        if len(cells) < 3 or set("".join(cells)) <= set("-: "):
            continue

        lowered = [_norm_label(c) for c in cells]
        header_hits = [i for i, c in enumerate(lowered) if c.startswith("shared")]

        if header_hits and not _BINOMIAL.match(re.sub(r"[*_`]", "", cells[0]).strip()):
            evidence_index = header_hits[0]
            continue

        name = re.sub(r"[*_`]", "", cells[0]).strip().strip(".")

        if not name or name.lower() in _NOT_A_SPECIES or not _BINOMIAL.match(name):
            continue
        if name in evidence or evidence_index >= len(cells):
            continue

        evidence[name] = ground_evidence(split_evidence_cell(cells[evidence_index]), table_rows)

        if len(evidence) >= limit:
            break

    return evidence


def filter_evidence(evidence, tables):
    """Restrict each candidate's categories to the chosen source tables.

    Categories that are not a verbatim table row have no source table and
    are always kept, so every category the candidate table cites is shown.
    """
    return {
        name: {key: value for key, value in items.items()
               if value[1] is None or value[1] in tables}
        for name, items in evidence.items()
    }


def abbreviate_species(name):
    parts = str(name).split()

    if len(parts) >= 2 and len(parts[0]) > 2 and not parts[0].endswith("."):
        return f"{parts[0][0]}. " + " ".join(parts[1:])

    return str(name)


def collect_candidate_phyla(text):
    """{species: current phylum name} from the candidate table's Phylum
    column, located by its header (default: second column, as the prompt
    specifies). The prompt asks for the current name with the older one in
    brackets, but the model does not always keep that order ("Firmicutes
    (Bacillota)"), so the name ending in -ota -- the suffix the 2021
    nomenclature rule gives every current phylum name -- is preferred, and
    the first name is used only when neither ends in -ota."""
    phyla = {}
    phylum_index = 1

    for line in str(text or "").splitlines():
        stripped = line.strip()

        if stripped.count("|") >= 3:
            cells = [c.strip() for c in stripped.strip("|").split("|")]
        elif "\t" in stripped:
            cells = [c.strip() for c in stripped.split("\t")]
        else:
            continue

        lowered = [_norm_label(c) for c in cells]
        if "phylum" in lowered:
            phylum_index = lowered.index("phylum")
            continue

        name = re.sub(r"[*_`]", "", cells[0]).strip().strip(".")
        if not _BINOMIAL.match(name) or phylum_index >= len(cells):
            continue

        names_given = [n.strip() for n in re.split(r"[()\[\]/,;]", re.sub(
            r"[*_`]", "", cells[phylum_index])) if n.strip()]
        current = [n for n in names_given if n.lower().endswith("ota")]
        phylum = (current or names_given or [""])[0]
        if phylum and name not in phyla:
            phyla[name] = phylum

    return phyla


def draw_organism_label(axes, x, y, name, phyla, fontsize, ha="left",
                        va="baseline", rotation=0, color=INK, zorder=4):
    """Draw "P. intermedia (Bacteroidota)": species italic, phylum upright.

    The two pieces are separate text objects -- so each is one clean,
    editable text element in the SVG -- with the phylum placed exactly after
    the measured width of the species name, along the label's rotation.
    The axes limits must already be set: positions are mapped through the
    axes transform.
    """
    import math
    from matplotlib.font_manager import FontProperties

    renderer = axes.figure.canvas.get_renderer()
    pixels_per_point = axes.figure.dpi / 72
    space = 0.3 * fontsize * pixels_per_point

    phylum = (phyla or {}).get(name)
    pieces = [(abbreviate_species(name), "italic")]
    if phylum:
        pieces.append((f"({phylum})", "normal"))

    widths = [
        renderer.get_text_width_height_descent(
            text, FontProperties(style=style, size=fontsize), ismath=False)[0]
        for text, style in pieces
    ]

    # Positions are worked out in display pixels, where "along the label"
    # is simply the rotation angle whatever the axes' direction or scale
    # (the heatmap's y-axis runs downward), then mapped back to data units.
    total = sum(widths) + space * (len(pieces) - 1)
    offset = -{"left": 0.0, "center": total / 2, "right": total}[ha]
    dx, dy = math.cos(math.radians(rotation)), math.sin(math.radians(rotation))
    anchor = axes.transData.transform((x, y))
    to_data = axes.transData.inverted()

    for (text, style), width in zip(pieces, widths):
        px, py = to_data.transform((anchor[0] + dx * offset, anchor[1] + dy * offset))
        axes.text(px, py, text, ha="left", va=va, rotation=rotation,
                  rotation_mode="anchor", fontstyle=style, fontsize=fontsize,
                  color=color, zorder=zorder)
        offset += width + space


def organism_label_chars(name, phyla=None):
    """Visible length of an organism label, for layout estimates."""
    phylum = (phyla or {}).get(name)
    return len(abbreviate_species(name)) + (len(phylum) + 3 if phylum else 0)


def row_id(index):
    """0 -> A, 25 -> Z, 26 -> AA, 27 -> AB, ..."""
    letters = ""
    index += 1
    while index:
        index, remainder = divmod(index - 1, 26)
        letters = chr(65 + remainder) + letters
    return letters


def venn_regions(named_sets):
    """{bitmask: sorted element keys}; bit i set = element is in set i."""
    names = list(named_sets)
    regions = {}

    for element in set().union(*named_sets.values()):
        mask = sum(1 << i for i, n in enumerate(names) if element in named_sets[n])
        regions.setdefault(mask, []).append(element)

    return {mask: sorted(members) for mask, members in regions.items()}


def region_name(mask, names):
    members = [abbreviate_species(n) for i, n in enumerate(names) if mask >> i & 1]

    if len(members) == len(names):
        return "All " + str(len(names)) if len(names) > 2 else "Both"

    return " ∩ ".join(members) + " only"


# Venn geometry for 2-5 sets, in data units (circle radius = 1).
#
# Each set is an ellipse (cx, cy, width, height, angle in degrees). 2 and 3
# sets are circles. 4 and 5 sets use the standard ellipse arrangements (as in
# pyvenn), rescaled to these units. They are true Venn diagrams: every one of
# the 2^n - 1 intersections exists as a region, which _venn_geometry checks
# on every call rather than trusting the constants. Six or more sets have no
# legible Venn layout; beyond five, the heatmap and UpSet views carry the
# comparison.
VENN_MAX_SETS = 5

_VENN_ELLIPSES = {
    2: [(-0.55, 0.0, 2.0, 2.0, 0.0), (0.55, 0.0, 2.0, 2.0, 0.0)],
    3: [(-0.55, 0.32, 2.0, 2.0, 0.0), (0.55, 0.32, 2.0, 2.0, 0.0),
        (0.0, -0.63, 2.0, 2.0, 0.0)],
    4: [(-0.600, -0.400, 2.88, 1.80, 140.0), (-0.200, 0.000, 2.88, 1.80, 140.0),
        (0.176, 0.000, 2.88, 1.80, 40.0), (0.576, -0.400, 2.88, 1.80, 40.0)],
    5: [(-0.288, -0.204, 3.48, 2.00, 155.0), (-0.124, 0.172, 3.48, 2.00, 82.0),
        (0.232, 0.092, 3.48, 2.00, 10.0), (0.312, -0.272, 3.48, 2.00, 118.0),
        (-0.044, -0.468, 3.48, 2.00, 46.0)],
}

# For 2-3 circles the titles sit above (and, for the third circle, below)
# the diagram. The two top titles are anchored outward from the centre line
# -- the left one ends there, the right one starts there -- so labels that
# carry a phylum never run into each other. Count anchors are computed for
# every n.
_VENN_HAND_PLACED = {
    2: {"titles": [(-0.12, 1.10, "right", "bottom"),
                   (0.12, 1.10, "left", "bottom")]},
    3: {"titles": [(-0.12, 1.42, "right", "bottom"),
                   (0.12, 1.42, "left", "bottom"),
                   (0.0, -1.72, "center", "top")]},
}


def _ellipse_half_extent(width, height, angle):
    import math
    a, b, t = width / 2, height / 2, math.radians(angle)
    return (math.hypot(a * math.cos(t), b * math.sin(t)),
            math.hypot(a * math.sin(t), b * math.cos(t)))


def _inside_ellipse(ellipse, x, y):
    import numpy as np
    cx, cy, width, height, angle = ellipse
    t = np.radians(angle)
    dx, dy = x - cx, y - cy
    u = dx * np.cos(t) + dy * np.sin(t)
    v = -dx * np.sin(t) + dy * np.cos(t)
    return (u / (width / 2)) ** 2 + (v / (height / 2)) ** 2 <= 1.0


@functools.lru_cache(maxsize=None)
def _venn_geometry(n_sets):
    """Region anchors and title anchors for an n-set Venn.

    Rasterises the ellipses, labels every pixel with the bitmask of the sets
    containing it, and anchors each region's count at its deepest point (the
    pixel farthest from the region's edge, found by eroding all regions at
    once). The depth is returned as the region's inscribed radius, which
    sizes the count so it fits small regions. Raises if any of the 2^n - 1
    regions is missing, so a wrong layout is never drawn.
    """
    import numpy as np

    ellipses = _VENN_ELLIPSES[n_sets]
    extents = [_ellipse_half_extent(w, h, a) for _, _, w, h, a in ellipses]
    x_low = min(e[0] - ext[0] for e, ext in zip(ellipses, extents)) - 0.05
    x_high = max(e[0] + ext[0] for e, ext in zip(ellipses, extents)) + 0.05
    y_low = min(e[1] - ext[1] for e, ext in zip(ellipses, extents)) - 0.05
    y_high = max(e[1] + ext[1] for e, ext in zip(ellipses, extents)) + 0.05

    pixel = max(x_high - x_low, y_high - y_low) / 480
    xs = np.arange(x_low, x_high, pixel)
    ys = np.arange(y_low, y_high, pixel)
    grid_x, grid_y = np.meshgrid(xs, ys)

    label = np.zeros(grid_x.shape, dtype=np.int64)
    for index, ellipse in enumerate(ellipses):
        label |= _inside_ellipse(ellipse, grid_x, grid_y).astype(np.int64) << index

    # Erode every region together: a pixel survives a round while all four
    # neighbours are still alive and in the same region.
    depth = np.zeros(label.shape, dtype=np.int32)
    alive = label > 0
    while alive.any():
        depth[alive] += 1
        same = np.ones(label.shape, dtype=bool)
        for shift, axis in ((1, 0), (-1, 0), (1, 1), (-1, 1)):
            same &= np.roll(alive, shift, axis) & (np.roll(label, shift, axis) == label)
        alive &= same

    regions = {}
    for mask in range(1, 2 ** n_sets):
        selected = label == mask
        if not selected.any():
            raise ValueError(f"{n_sets}-set layout has no region for {mask:0{n_sets}b}")
        deepest = depth[selected].max()
        rows, cols = np.nonzero(selected & (depth == deepest))
        pick = np.argmin((rows - rows.mean()) ** 2 + (cols - cols.mean()) ** 2)
        row, col = int(rows[pick]), int(cols[pick])

        regions[mask] = {
            "x": float(xs[col]), "y": float(ys[row]),
            "radius": float(deepest * pixel),
        }

    if n_sets in _VENN_HAND_PLACED:
        titles = list(_VENN_HAND_PLACED[n_sets]["titles"])
    else:
        # Title just outside the diagram, on the ray from the diagram's centre
        # through the set's exclusive region.
        centre_x = sum(e[0] for e in ellipses) / n_sets
        centre_y = sum(e[1] for e in ellipses) / n_sets
        titles = []

        for index in range(n_sets):
            x, y = regions[1 << index]["x"], regions[1 << index]["y"]
            dx, dy = x - centre_x, y - centre_y
            length = (dx ** 2 + dy ** 2) ** 0.5 or 1.0
            dx, dy = dx / length, dy / length

            while any(_inside_ellipse(e, x, y) for e in ellipses):
                x, y = x + dx * pixel, y + dy * pixel

            x, y = x + dx * 0.08, y + dy * 0.08
            ha = "left" if dx > 0.35 else "right" if dx < -0.35 else "center"
            va = "bottom" if dy > 0.35 else "top" if dy < -0.35 else "center"
            titles.append((x, y, ha, va))

    return {
        "ellipses": ellipses,
        "regions": regions,
        "titles": titles,
        "bounds": ((x_low + 0.05, x_high - 0.05), (y_low + 0.05, y_high - 0.05)),
    }


_VENN_LETTER_SIZE = 10.5     # pt, letters in a region with room for them
_VENN_MIN_LETTER_SIZE = 5.5  # pt, smallest letters in a cramped region


def _region_letter_lines(ids, radius_pt):
    """Wrap a region's heatmap letters ("A, B, C") to fit inside it.

    The region's inscribed circle bounds the text block. Tries font sizes
    from the full size down to the minimum and keeps the largest that fits;
    returns (lines, font size)."""
    tokens = [f"{i}," for i in ids[:-1]] + [ids[-1]]
    width_pt, height_pt = 1.7 * radius_pt, 1.5 * radius_pt
    size = _VENN_LETTER_SIZE

    while True:
        chars = max(1, int(width_pt / (0.62 * size)))
        lines, current = [], ""
        for token in tokens:
            candidate = f"{current} {token}".strip()
            if current and len(candidate) > chars:
                lines.append(current)
                current = token
            else:
                current = candidate
        lines.append(current)

        # Same number of lines, letters spread evenly across them, so a list
        # never ends on a lone letter ("A, B, C, D," / "E").
        per_line = -(-len(tokens) // len(lines))
        lines = [" ".join(tokens[i:i + per_line])
                 for i in range(0, len(tokens), per_line)]

        fits = len(lines) * size * 1.2 <= height_pt and \
            max(len(l) for l in lines) <= chars
        if fits or size <= _VENN_MIN_LETTER_SIZE:
            return lines, size
        size = max(_VENN_MIN_LETTER_SIZE, size - 0.5)


def _letter_sort_key(letter):
    return (len(letter), letter)


def draw_evidence_venn(named_sets, phyla=None, ids=None):
    """Draw an n-set Venn (2-5) labelled with the heatmap's category letters.

    Each region lists the letters of the categories held by exactly that
    combination of candidates (ids maps category key -> letter), so the
    diagram reads straight across to the heatmap rows: overlaps are what
    candidates share, single-candidate regions what each holds alone. Empty
    regions are left blank. Titles give each candidate's marker, italic
    species name and phylum. Returns (figure, {"svg", "pdf", "png": bytes}).
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Ellipse

    names = list(named_sets)
    n_sets = len(names)
    geometry = _venn_geometry(n_sets)
    regions = venn_regions(named_sets)
    ids = ids or {}

    (bx_low, bx_high), (by_low, by_high) = geometry["bounds"]
    target_width_in = {2: 4.4, 3: 4.4, 4: 5.6, 5: 6.6}[n_sets]
    points_per_unit = 72 * target_width_in / (bx_high - bx_low + 1.6)
    unit = 1.0 / points_per_unit

    name_size = 10
    marker_width, marker_gap = 7 * unit, 4 * unit
    title_height = name_size * 1.35 * unit

    title_blocks = []
    for name, (x, y, ha, va) in zip(names, geometry["titles"]):
        width = marker_width + marker_gap \
            + organism_label_chars(name, phyla) * 0.6 * name_size * unit
        left = x if ha == "left" else x - width if ha == "right" else x - width / 2
        bottom = y if va == "bottom" else y - title_height if va == "top" \
            else y - title_height / 2
        title_blocks.append((left, bottom, width))

    x_low = min([bx_low] + [b[0] for b in title_blocks]) - 0.08
    x_high = max([bx_high] + [b[0] + b[2] for b in title_blocks]) + 0.08
    y_low = min([by_low] + [b[1] for b in title_blocks]) - 0.08
    y_high = max([by_high] + [b[1] + title_height for b in title_blocks]) + 0.08

    fill_alpha = 0.22 if n_sets <= 3 else 0.16

    with plt.rc_context(FIGURE_STYLE):
        figure = plt.figure(figsize=((x_high - x_low) * points_per_unit / 72,
                                     (y_high - y_low) * points_per_unit / 72))
        axes = figure.add_axes([0, 0, 1, 1])
        axes.set_xlim(x_low, x_high)
        axes.set_ylim(y_low, y_high)
        axes.set_aspect("equal")
        axes.axis("off")

        for index, (cx, cy, w, h, angle) in enumerate(geometry["ellipses"]):
            axes.add_patch(Ellipse((cx, cy), w, h, angle=angle,
                                   facecolor=VENN_COLORS[index], alpha=fill_alpha,
                                   edgecolor="none", zorder=1))
            axes.add_patch(Ellipse((cx, cy), w, h, angle=angle, facecolor="none",
                                   edgecolor=VENN_STROKES[index],
                                   linewidth=1.4 if n_sets <= 3 else 1.1,
                                   zorder=2))

        for mask, room in geometry["regions"].items():
            members = regions.get(mask, [])
            if not members:
                continue

            letters = sorted((ids.get(m, "?") for m in members), key=_letter_sort_key)
            lines, size = _region_letter_lines(letters, room["radius"] * points_per_unit)
            axes.text(room["x"], room["y"], "\n".join(lines), ha="center",
                      va="center", fontsize=size, fontweight="bold",
                      linespacing=1.2, color=INK, zorder=4)

        # Title: marker (shape + colour), then italic species and phylum.
        for index, (name, (left, bottom, width)) in enumerate(
                zip(names, title_blocks)):
            middle = bottom + title_height / 2
            axes.plot(left + marker_width / 2, middle,
                      marker=VENN_MARKERS[index], markersize=6,
                      markerfacecolor=VENN_COLORS[index],
                      markeredgecolor=VENN_STROKES[index],
                      markeredgewidth=0.9, linestyle="none", zorder=4)
            draw_organism_label(axes, left + marker_width + marker_gap, middle,
                                name, phyla, name_size, ha="left",
                                va="center_baseline")

        exports = export_figure(figure)

    return figure, exports


def export_figure(figure):
    """SVG (editable text), PDF (TrueType) and 600-dpi PNG of a figure.

    Every format has an opaque white background: the inks are chosen for
    white, and a transparent export placed on a dark slide would lose them.

    Call inside the FIGURE_STYLE context: the font-embedding settings are
    read at save time.
    """
    exports = {}

    for fmt in ("svg", "pdf", "png"):
        buffer = io.BytesIO()
        figure.savefig(buffer, format=fmt, dpi=600 if fmt == "png" else None,
                       bbox_inches="tight", pad_inches=0.08,
                       facecolor="white")
        exports[fmt] = buffer.getvalue()

    return exports


def evidence_matrix(filtered, names):
    """Presence/absence of every cited category across candidates.

    One row per category, ordered for display: by how many candidates cite
    it, then by presence pattern in candidate order -- so categories shared
    by the same candidates sit together -- then alphabetically. Rows are
    lettered A, B, C, ... in that order.
    """
    labels = {}
    for items in filtered.values():
        for key, value in items.items():
            labels.setdefault(key, value)

    rows = []

    for key, (label, table_number) in labels.items():
        present = [key in filtered[name] for name in names]
        rows.append({
            "key": key,
            "Category": label,
            "Table": table_number,
            "Source": _table_tag(table_number),
            **dict(zip(names, present)),
            "Cited by": sum(present),
            "_pattern": tuple(not p for p in present),
        })

    if not rows:
        return pd.DataFrame()

    frame = pd.DataFrame(rows)
    frame["_label"] = frame["Category"].str.lower()
    frame = frame.sort_values(["Cited by", "_pattern", "_label"],
                              ascending=[False, True, True])
    frame = frame.drop(columns=["_pattern", "_label"]).reset_index(drop=True)
    frame.insert(0, "ID", [row_id(i) for i in range(len(frame))])
    return frame


def _fit_text(text, max_chars):
    text = " ".join(str(text).split())
    return text if len(text) <= max_chars else text[: max_chars - 1].rstrip() + "…"


def draw_evidence_heatmap(matrix, names):
    """Binary presence/absence heatmap: categories down, candidates across.

    Each row reads "A  Category" then one cell per candidate and the
    "Shared by" fraction. Column heads are the italic species names; the
    phylum is given only in the Venn titles beside it. Drawn in cell units (1 cell = 1 data unit) on an
    axes that fills the figure, with the figure sized from the row/column
    counts, so cells stay square and text sizes stay fixed however many rows
    there are.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    cell_in = 0.19
    ppu = cell_in * 72                          # points per cell
    label_chars = 72
    row_font, head_font = 7, 7.5
    gap = 2 / ppu                               # 2 pt surface gap between cells
    heading_rows = 1.4

    n_rows, n_cols = len(matrix), len(names)
    labels = [_fit_text(c, label_chars) for c in matrix["Category"]]
    number_units = max(len(i) for i in matrix["ID"]) * 0.68 * row_font / ppu
    label_units = max(len(l) for l in labels) * 0.52 * row_font / ppu
    left_units = number_units + 0.45 + label_units + 1.0
    column_head_units = max(organism_label_chars(n) for n in names) \
        * 0.52 * head_font / ppu * 0.87 + 1.0
    count_units = 2.2

    x_low, x_high = -left_units, n_cols + count_units
    y_top, y_bottom = -column_head_units, heading_rows + n_rows + 0.4
    number_right = x_low + 0.15 + number_units

    with plt.rc_context(FIGURE_STYLE):
        figure = plt.figure(figsize=((x_high - x_low) * cell_in,
                                     (y_bottom - y_top) * cell_in))
        axes = figure.add_axes([0, 0, 1, 1])
        axes.set_xlim(x_low, x_high)
        axes.set_ylim(y_bottom, y_top)
        axes.axis("off")

        for col, name in enumerate(names):
            draw_organism_label(axes, col + 0.5, -0.35, name, None,
                                head_font, va="bottom", rotation=60)

        axes.text(n_cols + 0.6, -0.35, "Shared by", rotation=60,
                  rotation_mode="anchor", ha="left", va="bottom",
                  fontsize=row_font, color=INK_MUTED)

        axes.text(x_low + 0.15, 0.62, "Categories", ha="left", va="center",
                  fontsize=head_font, fontweight="bold", color=INK)
        axes.plot([x_low + 0.15, n_cols + count_units - 0.2], [1.0, 1.0],
                  color=RULE, linewidth=0.5)

        for i, (_, row) in enumerate(matrix.iterrows()):
            y = heading_rows + i

            axes.text(number_right, y + 0.5, row["ID"], ha="right",
                      va="center", fontsize=row_font, fontweight="bold",
                      color=INK_MUTED)
            axes.text(number_right + 0.45, y + 0.5, labels[i], ha="left",
                      va="center", fontsize=row_font, color=INK)

            for col, name in enumerate(names):
                axes.add_patch(Rectangle(
                    (col + gap / 2, y + gap / 2), 1 - gap, 1 - gap,
                    facecolor=MARK_DARK if row[name] else MARK_EMPTY,
                    edgecolor="none",
                ))

            axes.text(n_cols + 0.45, y + 0.5, f"{row['Cited by']}/{n_cols}",
                      ha="left", va="center", fontsize=row_font,
                      color=INK_MUTED)

        exports = export_figure(figure)

    return figure, exports


def _figure_downloads(exports, stem, csv_label, csv_frame, csv_name):
    """Row of SVG / PDF / PNG / CSV download buttons."""
    columns = st.columns(4)
    items = (
        ("SVG", exports["svg"], f"{stem}.svg", "image/svg+xml"),
        ("PDF", exports["pdf"], f"{stem}.pdf", "application/pdf"),
        ("PNG, 600 dpi", exports["png"], f"{stem}.png", "image/png"),
        (csv_label, csv_frame.to_csv(index=False), csv_name, "text/csv"),
    )

    for column, (label, data, file_name, mime) in zip(columns, items):
        with column:
            st.download_button(label, data=data, file_name=file_name, mime=mime,
                               use_container_width=True,
                               key=f"web6_dl_{file_name}")


def _table_tag(table_number):
    return f"Table {table_number}" if table_number else "not in Tables 1-4"


def render_shared_evidence(response_text, table_rows):
    """Heatmap and Venn, side by side, of the categories candidates share."""
    import matplotlib.pyplot as plt

    evidence = collect_candidate_evidence(response_text, table_rows)
    phyla = collect_candidate_phyla(response_text)

    if len(evidence) < 2:
        st.info(
            "Fewer than two candidates with a readable evidence column were "
            "found in the candidate table, so there is nothing to compare."
        )
        return

    if not table_rows:
        st.warning(
            "The processed tables from this run are not available, so cited "
            "categories cannot be matched to their rows. Run the analysis again."
        )
        return

    names = list(evidence)

    tables = st.multiselect(
        "Evidence from",
        options=[1, 2, 3, 4],
        default=[1, 2, 3, 4],
        format_func=lambda t: TABLE_LABELS[t],
        key="web6_venn_tables",
    )

    filtered = filter_evidence(evidence, set(tables))
    matrix = evidence_matrix(filtered, names)

    if matrix.empty:
        st.info("No cited categories remain under the current filters.")
        return

    # Venn candidates: all of them when they fit; otherwise the top five by
    # the model's confidence order, with a choice. The key carries this
    # result's candidate list, so a new analysis starts from its own
    # candidates instead of a selection remembered from a previous run.
    venn_names = names[:VENN_MAX_SETS]
    if len(names) > VENN_MAX_SETS:
        venn_names = st.multiselect(
            f"Candidates in the Venn diagram (2-{VENN_MAX_SETS}; the heatmap "
            f"shows all {len(names)})",
            options=names,
            default=names[:VENN_MAX_SETS],
            max_selections=VENN_MAX_SETS,
            key="web6_venn_candidates_" + str(abs(hash(tuple(names)))),
        )

    heatmap_column, venn_column = st.columns([1, 1])

    with heatmap_column:
        figure, heatmap_exports = draw_evidence_heatmap(matrix, names)
        st.pyplot(figure, use_container_width=True)
        plt.close(figure)

    with venn_column:
        named_sets = {name: set(filtered[name]) for name in venn_names}

        if len(venn_names) < 2:
            st.info("Select at least two candidates for the Venn diagram.")
            venn_exports = None
        else:
            ids = dict(zip(matrix["key"], matrix["ID"]))
            figure, venn_exports = draw_evidence_venn(named_sets, phyla, ids)
            st.pyplot(figure, use_container_width=True)
            plt.close(figure)

    st.caption(
        "Left: each filled cell is a category the candidate's evidence cell "
        "cites as shared with the query genome; rows are ordered by how many "
        "candidates share them, and candidates by the model's confidence. "
        "Right: the same categories by their heatmap letters - overlaps are "
        "categories the candidates have in common, single-candidate regions "
        "are categories only that candidate cites; empty regions are blank."
    )

    csv = matrix[["ID", "Category", "Source"] + names + ["Cited by"]].copy()
    csv[names] = csv[names].astype(int)

    st.markdown("**Heatmap**")
    _figure_downloads(heatmap_exports, "candidate_categories_heatmap",
                      "Matrix (CSV)", csv, "candidate_categories_matrix.csv")

    if venn_exports:
        region_rows = []
        for mask, members in sorted(
            venn_regions(named_sets).items(),
            key=lambda kv: (-bin(kv[0]).count("1"), kv[0]),
        ):
            title = region_name(mask, venn_names)
            region_rows.extend(
                {"Region": title, "ID": ids[m],
                 "Category": next(filtered[n][m][0] for n in venn_names
                                  if m in filtered[n]),
                 "Source": _table_tag(next(filtered[n][m][1] for n in venn_names
                                           if m in filtered[n]))}
                for m in sorted(members, key=lambda k: _letter_sort_key(ids[k]))
            )

        st.markdown("**Venn diagram**")
        _figure_downloads(venn_exports, "candidate_categories_venn",
                          "Regions (CSV)", pd.DataFrame(region_rows),
                          "candidate_categories_venn_regions.csv")


def _to_markdown_table(df: pd.DataFrame) -> str:
    try:
        return df.to_markdown(index=False)
    except ImportError:
        return df.to_csv(index=False, sep="|")


def build_combined_prompt(
    functional_category_df: pd.DataFrame,
    functional_category_note: str,
    kegg_pathway_df: pd.DataFrame,
    kegg_pathway_note: str,
    protein_family_df: pd.DataFrame,
    pathway_completeness_df: pd.DataFrame,
    absent_pathway_groups: list,
    base_description: str,
    extra_instructions: str = "",
) -> str:
    prompt = base_description.strip() + "\n\n"
    prompt += "Table 1 - Functional category gene counts " + functional_category_note + ":\n"
    prompt += _to_markdown_table(functional_category_df) + "\n\n"
    prompt += "Table 2 - KEGG-hit description gene counts " + kegg_pathway_note + ":\n"
    prompt += _to_markdown_table(kegg_pathway_df) + "\n\n"
    prompt += "Table 3 - Protein function family gene counts:\n" + _to_markdown_table(protein_family_df) + "\n\n"
    prompt += "Table 4 - Metabolic pathway completeness (present/complete pathways only):\n"
    prompt += _to_markdown_table(pathway_completeness_df) + "\n\n"
    if absent_pathway_groups:
        prompt += (
            "Pathway GROUPS with ZERO present/complete rows in Table 4 (entirely absent - "
            "treat as meaningful negative evidence, not missing data):\n"
            + ", ".join(sorted(absent_pathway_groups)) + "\n\n"
        )
    if extra_instructions.strip():
        prompt += "Additional instructions:\n" + extra_instructions.strip() + "\n"
    return prompt


# ---------------------------------------------------------------------------
# OpenAI call - uses the same session-state keys master_app.py's
# api_key_setup() already populates (user_openai_api_key, user_selected_model,
# user_temperature, user_max_tokens). No new config screen needed.
# ---------------------------------------------------------------------------

DEFAULT_TEMPERATURE = 0.5  # this tab uses the dashboard temperature setting
# unchanged, like every other GenoAnno tab, so all modules are queried under
# identical sampling conditions and cross-module differences reflect the
# annotation source rather than a per-tab parameter. Run-to-run stability of
# this tab is assessed separately rather than enforced here; the prompt's own
# DETERMINISM RULE still asks the model for a reproducible answer.


def call_openai(prompt: str) -> str:
    from openai import OpenAI

    api_key = st.session_state.get("user_openai_api_key", "")
    model = st.session_state.get("user_selected_model", "gpt-4o-mini")

    temperature = st.session_state.get("user_temperature", DEFAULT_TEMPERATURE)
    try:
        temperature = float(temperature)
    except (TypeError, ValueError):
        temperature = DEFAULT_TEMPERATURE

    max_tokens = st.session_state.get("user_max_tokens", 2000)

    if not api_key:
        return "No OpenAI API key found in this session. Please enter one on the dashboard screen."

    client = OpenAI(api_key=api_key)
    messages = [{"role": "user", "content": prompt}]

    try:
        completion = client.chat.completions.create(
            model=model,
            messages=messages,
            temperature=temperature,
            max_tokens=max_tokens,
        )
    except Exception:
        # Newer reasoning-style models (e.g. gpt-5 family) reject `max_tokens`
        # and fixed `temperature` on the chat.completions endpoint - retry
        # with the parameters they do accept.
        try:
            completion = client.chat.completions.create(
                model=model,
                messages=messages,
                max_completion_tokens=max_tokens,
            )
        except Exception as exc:  # surface the real error in the UI
            return f"OpenAI request failed: {exc}"

    return completion.choices[0].message.content


# ---------------------------------------------------------------------------
# UI - runs at module (import) level, since master_app.py executes this file
# with runpy.run_path rather than as __main__.
# ---------------------------------------------------------------------------

# Every session-state key that belongs to one analysis run. Cleared together
# when the inputs disappear or the user clears the result.
_RUN_STATE_KEYS = (
    "web6_response",
    "web6_table_rows",
)

st.header("Find Similar Oral Bacteria")
st.caption(
    "Upload your three raw annotation files. This tab derives four summary "
    "tables from them (functional categories, KEGG pathways, protein "
    "families, pathway completeness) and combines all four into one prompt "
    "to search for bacteria with a similar gene composition."
)

col1, col2, col3 = st.columns(3)
with col1:
    bakta_file = st.file_uploader(
        "Bakta/Proksee annotation (.tsv)", type=["tsv"], key="bakta_upload"
    )
with col2:
    kbase_file = st.file_uploader(
        "KBASE annotation (.tsv)", type=["tsv"], key="kbase_upload"
    )
with col3:
    products_file = st.file_uploader(
        "Pathway completeness table, e.g. Products.tsv (.tsv)",
        type=["tsv"],
        key="products_upload",
    )

input_description = st.text_area(
    "Prompt sent to the model",
    value=DEFAULT_INPUT_DESCRIPTION,
    height=380,
    key="input_description",
)

extra_instructions = st.text_area(
    "Any additional instructions (optional)", value="", height=80, key="extra_instructions"
)

if st.button("Find similar bacteria", type="primary", key="run_similar_bacteria"):
    if not (bakta_file and kbase_file and products_file):
        st.warning("Please upload all three files first.")
    else:
        with st.spinner("Processing annotation files..."):
            functional_category_df, functional_category_note = parse_bakta_functional_categories(bakta_file)
            kbase_file.seek(0)
            kegg_pathway_df, kegg_pathway_note = parse_kbase_kegg_pathways(kbase_file)
            kbase_file.seek(0)
            protein_family_df = parse_kbase_protein_families(kbase_file)
            pathway_completeness_df, absent_pathway_groups = parse_products_completeness(products_file)

        with st.expander("Processed tables (review before sending to the model)"):
            st.subheader("Table 1 - Functional category gene counts")
            st.caption(functional_category_note)
            st.dataframe(functional_category_df)
            st.subheader("Table 2 - KEGG-hit description gene counts")
            st.caption(kegg_pathway_note)
            st.dataframe(kegg_pathway_df)
            st.subheader("Table 3 - Protein function family gene counts")
            st.dataframe(protein_family_df)
            st.subheader("Table 4 - Metabolic pathway completeness")
            st.dataframe(pathway_completeness_df)
            if absent_pathway_groups:
                st.caption("Entirely absent pathway groups: " + ", ".join(sorted(absent_pathway_groups)))

        base_description = input_description if input_description.strip() else DEFAULT_INPUT_DESCRIPTION
        prompt = build_combined_prompt(
            functional_category_df,
            functional_category_note,
            kegg_pathway_df,
            kegg_pathway_note,
            protein_family_df,
            pathway_completeness_df,
            absent_pathway_groups,
            base_description,
            extra_instructions,
        )

        with st.spinner("Querying the model..."):
            full_response = clean_model_output(call_openai(prompt))

        # Held in session state so the result and its figures survive the
        # reruns Streamlit performs when a widget below is used. Without this,
        # changing a figure control would re-trigger the model call.
        st.session_state["web6_response"] = full_response
        # Row labels of Tables 1-4, so the Venn section can match the model's
        # cited items back to verbatim rows without re-parsing the uploads.
        st.session_state["web6_table_rows"] = collect_table_rows(
            functional_category_df,
            kegg_pathway_df,
            protein_family_df,
            pathway_completeness_df,
        )


# ---------------------------------------------------------------------------
# Result and shared-evidence figures
#
# Rendered outside the button block and driven from session state, so the
# figure controls can rerun on their own without calling the model again.
# ---------------------------------------------------------------------------

# A stored result belongs to the files that produced it. Closing the tool and
# reopening it, or clearing the uploaders, leaves no files but -- because the
# result lives in session state so the figure controls can rerun without
# re-querying the model -- would otherwise leave the previous run's output on
# screen as if it described the new (empty) input. Drop it instead.
if not (bakta_file and kbase_file and products_file):
    for stale_key in _RUN_STATE_KEYS:
        st.session_state.pop(stale_key, None)


if st.session_state.get("web6_response"):
    stored_response = st.session_state["web6_response"]

    result_header, clear_column = st.columns([4, 1])

    with result_header:
        st.subheader("Result")

    with clear_column:
        if st.button("Clear result", key="web6_clear_result",
                     use_container_width=True):
            for stale_key in _RUN_STATE_KEYS:
                st.session_state.pop(stale_key, None)

            st.rerun()

    st.markdown(strip_step_labels(strip_step0_section(stored_response)))

    # The STEP 0 checklist is still requested and generated (it grounds the
    # later steps) but not shown; it remains in the downloaded interpretation.

    st.download_button(
        label="Download AI interpretation",
        data=stored_response,
        file_name="similar_oral_bacteria_interpretation.txt",
        mime="text/plain",
    )

    # Heatmap and Venn of the categories each candidate shares with the query
    # genome, read off the candidate table directly.
    st.divider()
    st.subheader("Shared evidence among candidate bacteria")
    render_shared_evidence(stored_response, st.session_state.get("web6_table_rows", {}))
