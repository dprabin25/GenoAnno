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

import io
import re

import pandas as pd
import streamlit as st

# ---------------------------------------------------------------------------
# Combined prompt template
# ---------------------------------------------------------------------------

DEFAULT_INPUT_DESCRIPTION = """I annotated genes from a bacterial genome (an oral bacterium) using multiple annotation pipelines and summarized the results into four processed tables describing its functional gene composition.

Table 1 - Functional category gene counts (Bakta):
I selected the Product column from the Bakta annotation and counted how many genes map to each repeated product/function term. Ultra-generic, non-informative labels ("hypothetical protein", "uncharacterized protein") have been removed, and the table has been capped to the top most frequent categories - it is a partial, curated list, not the full gene set.
Column 1: Functional category
Column 2: Gene count

Table 2 - KEGG-annotated gene function descriptions (KBASE):
Each row is the free-text function description ("kegg_hit") tied to a gene's best KO match, counted by how many genes share that exact description. These are gene-level descriptions, NOT a curated pathway hierarchy - there is no row literally named "glycolysis" or "TCA cycle" unless it appears verbatim below. Do not invent or assume pathway-level labels that are not literal rows in this table. Generic/uninformative rows have been removed and the table capped to the top most frequent entries.
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

STEP 0 - Before any narrative reasoning, extract a fixed evidence checklist as a short bullet list, pulling only discriminating (non-universal) rows verbatim from Tables 1-4: the specific CAZy/carbohydrate categories present (if any) and the fact of total CAZy absence (if that's the case), the specific SCFA/fermentation end-product pathways present, the specific respiratory complexes present vs. entirely absent, any unusual protein families (Table 3) or KEGG-hit descriptions (Table 2) that are not generic housekeeping, and every cell-surface, secretion-system and virulence-associated category in Table 1 (surface/S-layer proteins, secretion systems, motility apparatus, adhesins, capsule, toxins, host-substrate-degrading enzymes) - include these even when their gene count is 1, because such genes are typically single-copy yet are among the most taxonomically specific evidence available. This checklist is the sole basis for every later step - do not add or reason from anything outside it.

STEP 1 - Using only the STEP 0 checklist, derive an internal phenotype profile of this genome. Explicitly reason about:
- Respiratory/fermentative capacity: which respiratory chain complexes (Table 1/Table 4) are present or absent, and whether the profile looks aerobic, microaerophilic, anaerobic, or strictly fermentative.
- Carbon source usage: which CAZy/carbohydrate-active categories (Table 1) are present vs. largely absent - i.e., whether this organism looks broadly saccharolytic (degrades many sugars) or asaccharolytic/limited (relies on other substrates such as amino acids, lactate, or organic acids).
- Fermentation end products: which SCFA/alcohol conversion pathways (Table 4) are complete (e.g., lactate, propionate, acetate, butyrate). End-product profile is one of the strongest, most specific phenotype signals for oral bacteria and often more diagnostic than raw gene overlap.
- Distinctive protein families (Table 3) and KEGG pathways (Table 2) beyond universal housekeeping genes.
- Cell surface architecture and host interaction: which surface-layer or surface-glycoprotein components, secretion systems, motility/gliding apparatus, adhesins, capsule machinery, toxins and host-substrate-degrading enzymes (Table 1) are present. Treat these as high-weight evidence despite low gene counts - two organisms with near-identical fermentation strategies can still be different genera, and when that is the case the surface and secretion machinery is usually what separates them. State explicitly which of these the genome does and does not encode.

STEP 2 - When weighing evidence across all four tables, explicitly downweight categories that are near-universal across nearly all bacteria and carry little taxonomic signal (e.g., "hypothetical protein", generic ABC transporters, ribosomal proteins, DNA replication/repair machinery, chaperones, generic transport systems). These will dominate the raw gene counts but should NOT drive the match. Instead prioritize categories that vary meaningfully between oral genera: CAZy/carbohydrate degradation breadth, SCFA/fermentation end products, respiratory chain complex composition, nitrogen/sulfur metabolism, and any pathway present in some oral taxa but absent in most others.

STEP 3 - Only after deriving this phenotype profile, identify well-characterized oral bacteria whose published phenotype (fermentation strategy, carbon source range, oxygen tolerance, cell/colony morphology, cell surface architecture, secretion systems and characterised virulence mechanisms) is genuinely consistent with the derived profile - not simply bacteria that happen to share the largest raw count of generic annotated genes. If a well-known genus overlaps only on universal housekeeping categories, exclude it even if its raw overlap count is high. Where several candidates are metabolically indistinguishable, let the surface, secretion and virulence evidence decide the ranking, and name the specific feature that separates your top candidate from the runner-up. A candidate that is metabolically plausible but cannot account for the genome's surface and secretion machinery is a weaker match than one that explains both, and should be ranked accordingly or excluded.

Question: Identify other oral bacteria that have similar gene composition to mine. Use all four tables together - a strong match should show consistent overlap across functional categories (Table 1), KEGG pathways (Table 2), protein function families (Table 3), and completed metabolic pathways (Table 4), weighted as described above, not just raw agreement in one table.

Output, in this order:
1. The STEP 0 evidence checklist (bullet list, verbatim rows only).
2. A short (2-4 sentence) summary of the phenotype profile you derived in Step 1, so the reasoning can be sanity-checked.
3. A table with columns: Bacterium name | Phylum | Shared functional categories, pathways, and families (cite which of Tables 1-4 support the match, and note whether each is a discriminating or universal category) | Phenotype description | Confidence.

For the Phylum column, give the candidate's phylum and, in brackets, whether it is Gram-positive or Gram-negative - for example "Bacillota [Gram-positive]" or "Bacteroidota [Gram-negative]". Before citing any row as shared evidence, check that the candidate's cell envelope can actually possess it: outer-membrane systems (TonB-dependent receptors and their plug domains, SusC/SusD/RagA/RagB nutrient uptake, LPS export, outer-membrane beta-barrel proteins, TolC-family channels, type IX secretion) exist only in Gram-negative organisms and cannot be shared with a Gram-positive candidate. If a row is incompatible with a candidate's envelope, do not list it as shared evidence for that candidate, and say in the Shared column which evidence had to be excluded on envelope grounds.
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

TOP_N_ROWS = 40  # cap on how many rows go into Tables 1 and 2


# ---------------------------------------------------------------------------
# Discriminating-row rescue for Table 1
# ---------------------------------------------------------------------------
#
# Ranking Table 1 by gene count selects for gene-family EXPANSION, which is the
# least species-specific part of a genome: in Bacteroidota the top rows are
# always SusC/SusD/sigma-factor PUL machinery, shared across the whole phylum.
# The genes that distinguish one species from its relatives -- secreted
# proteases, surface layers, secretion systems, adhesins -- are almost all
# single-copy. In the Tannerella forsythia input 85.6% of distinct CDS products
# occur exactly once and the rank-40 cutoff sits at a count of 4, so every such
# gene is excluded before the model sees the table.
#
# These two selectors rescue those rows. Both are properties of the ANNOTATION,
# not of any particular organism, so the same code applied to any Bakta table
# surfaces whatever that genome's characterised, virulence-associated genes are
# -- gingipains and fimbrillins in Porphyromonas, FadA in Fusobacterium,
# leukotoxin in Aggregatibacter. No species, gene or locus name appears here,
# so the model is not handed a lookup table for an expected answer.
#
#   Selector A -- a BlastRules cross-reference. NCBI BlastRules are curated,
#   high-specificity protein-naming rules that fire only for well-studied
#   proteins, so the xref marks a gene the reference databases treat as
#   diagnostically named rather than generically matched. Typically a handful
#   of genes per genome.
#
#   Selector C -- a generic functional-class keyword. Terms name the CLASS of
#   virulence-associated function (secretion system, surface layer, adhesin,
#   toxin, motility, capsule, host-substrate-degrading enzyme), as any
#   a-priori virulence panel would. Matching is on the product description, so
#   organism-specific protein names are matched incidentally by their class
#   words and never enumerated.
#
# Selector B from the design discussion -- "has an assigned gene symbol" -- is
# deliberately NOT used: it admits ~10% of all CDS, mostly housekeeping genes
# (uvrA, greA, pfkB), adding bulk without discriminating power.

DISCRIMINATING_XREF_RULES = ("BlastRules",)

DISCRIMINATING_KEYWORDS = (
    "secretion system",
    "s-layer",
    "surface layer",
    "adhesin",
    "adhesion",
    "virulence",
    "toxin",
    "motility",
    "capsul",
    "fimbri",
    "pilus",
    "pili",
    "hemolysin",
    "haemolysin",
    "invasin",
    "autotransporter",
    "sialid",
    "neuraminidase",
    "mucinase",
    "collagenase",
    "hyaluronidase",
    "gingipain",
    "leukotoxin",
)


def select_discriminating_products(df: pd.DataFrame) -> pd.DataFrame:
    """Rows worth keeping in Table 1 regardless of how rare they are.

    Returns a (Functional category, Gene count) frame, counted over the same
    CDS rows Table 1 is built from. Returns an empty frame when the table
    carries neither of the columns the selectors read, so a non-Bakta source
    simply contributes nothing.
    """
    if "Product" not in df.columns:
        return pd.DataFrame(columns=["Functional category", "Gene count"])

    product = df["Product"].fillna("").astype(str)
    keep = pd.Series(False, index=df.index)

    # Selector A: curated protein-naming rule fired for this gene.
    if "DbXrefs" in df.columns:
        xrefs = df["DbXrefs"].fillna("").astype(str)
        for rule in DISCRIMINATING_XREF_RULES:
            keep = keep | xrefs.str.contains(rule, case=False, regex=False)

    # Selector C: product description names a virulence-associated class.
    for keyword in DISCRIMINATING_KEYWORDS:
        keep = keep | product.str.contains(keyword, case=False, regex=False)

    selected = product[keep & (product.str.strip() != "")]

    if selected.empty:
        return pd.DataFrame(columns=["Functional category", "Gene count"])

    return (
        selected.value_counts()
        .rename_axis("Functional category")
        .reset_index(name="Gene count")
    )


def _clean_and_cap(counts: pd.DataFrame, label_col: str, top_n: int = TOP_N_ROWS):
    """Drop generic/uninformative labels, sort by count, keep the top N.

    Returns (trimmed_df, note) where note documents what was omitted so the
    prompt can be transparent about the fact that this is a partial table.
    """
    total_rows = len(counts)
    filtered = counts[~counts[label_col].str.strip().str.lower().isin(GENERIC_LABELS)]
    dropped_generic = total_rows - len(filtered)
    filtered = filtered.sort_values("Gene count", ascending=False).reset_index(drop=True)
    trimmed = filtered.head(top_n)
    omitted = len(filtered) - len(trimmed)
    note = (
        f"(showing top {len(trimmed)} of {total_rows} total categories by gene count; "
        f"{dropped_generic} generic/non-discriminating rows removed"
        + (f"; {omitted} additional lower-count categories omitted for brevity)" if omitted > 0 else ")")
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

    trimmed, note = _clean_and_cap(counts, "Functional category")

    # Add back the discriminating rows the count cutoff would have discarded.
    rescued = select_discriminating_products(df)

    if not rescued.empty:
        already_present = set(trimmed["Functional category"])
        rescued = rescued[
            ~rescued["Functional category"].isin(already_present)
            & ~rescued["Functional category"].str.strip().str.lower().isin(GENERIC_LABELS)
        ]

    if not rescued.empty:
        trimmed = (
            pd.concat([trimmed, rescued], ignore_index=True)
            .sort_values(
                by=["Gene count", "Functional category"],
                ascending=[False, True],
            )
            .reset_index(drop=True)
        )

        # _clean_and_cap's note ends in ")"; extend it rather than replace it,
        # so the table stays self-describing in the prompt.
        note = note.rstrip()
        note = (
            note[:-1] if note.endswith(")") else note
        ) + (
            f"; plus {len(rescued)} virulence-associated categories retained "
            "regardless of gene count, selected by curated-annotation flag or "
            "functional-class keyword)"
        )

    return trimmed, note


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
# Model output cleanup
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
# Prompt assembly
# ---------------------------------------------------------------------------

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

        # Shown: phenotype summary onward. Downloaded: the whole response,
        # including the STEP 0 checklist that grounds it.
        st.subheader("Result")
        st.markdown(strip_step_labels(strip_step0_section(full_response)))

        with st.expander("Show STEP 0 evidence checklist", expanded=False):
            st.markdown(full_response)

        st.download_button(
            label="Download AI interpretation",
            data=full_response,
            file_name="similar_oral_bacteria_interpretation.txt",
            mime="text/plain",
        )
