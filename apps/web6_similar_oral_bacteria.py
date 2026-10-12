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
import json
import re
import time
import urllib.parse
import urllib.request

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
# 16S phylogeny of the candidate species
# ---------------------------------------------------------------------------
#
# The candidate table names species; this section fetches one 16S rRNA
# sequence per species from NCBI Nucleotide and draws a neighbour-joining tree
# so the candidates can be seen in relation to each other rather than only
# read as a list. A species with no usable 16S record is reported as missing
# rather than silently dropped.
#
# SCOPE: distances here are 1 - pairwise identity from a Biopython global
# alignment, not a multiple-sequence alignment, and no model selection or
# bootstrapping is performed. That is adequate for a quick look at whether the
# candidates cluster sensibly, and NOT adequate for publication. The FASTA
# download exists for that: realign it (MAFFT/MUSCLE) and infer the tree with
# IQ-TREE or RAxML with support values if the figure is going into a paper.

NCBI_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# NCBI asks that programmatic users identify themselves and stay under 3
# requests/second without an API key (10/second with one).
NCBI_TOOL_NAME = "GenoAnno"
NCBI_REQUEST_INTERVAL = 0.4

# 16S in bacteria is ~1540 bp. The length bound keeps the hit list to
# full-length gene records and rejects short partial-sequence submissions and
# whole genomes that merely mention 16S.
SIXTEEN_S_MIN_LEN = 1200
SIXTEEN_S_MAX_LEN = 1700


def _eutils_get(endpoint, params, timeout=30):
    url = f"{NCBI_EUTILS}/{endpoint}?" + urllib.parse.urlencode(params)
    request = urllib.request.Request(url, headers={"User-Agent": NCBI_TOOL_NAME})

    with urllib.request.urlopen(request, timeout=timeout) as response:
        return response.read().decode("utf-8", errors="replace")


@st.cache_data(show_spinner=False, ttl=60 * 60 * 24)
def fetch_16s(organism, email="", api_key=""):
    """One 16S rRNA sequence for `organism` from NCBI Nucleotide.

    Returns (accession, sequence) or (None, reason). Cached for a day so
    repeated runs on the same candidate list do not re-query NCBI.
    """
    organism = str(organism).strip()

    if not organism:
        return None, "no name given"

    shared = {"tool": NCBI_TOOL_NAME}

    if email:
        shared["email"] = email
    if api_key:
        shared["api_key"] = api_key

    term = (
        f'"{organism}"[Organism] AND '
        f'("16S ribosomal RNA"[Title] OR "16S rRNA"[Title]) AND '
        f"{SIXTEEN_S_MIN_LEN}:{SIXTEEN_S_MAX_LEN}[SLEN]"
    )

    try:
        payload = json.loads(
            _eutils_get(
                "esearch.fcgi",
                {
                    **shared,
                    "db": "nucleotide",
                    "term": term,
                    "retmode": "json",
                    "retmax": "1",
                    "sort": "relevance",
                },
            )
        )
    except Exception as error:
        return None, f"NCBI search failed ({type(error).__name__})"

    ids = payload.get("esearchresult", {}).get("idlist", [])

    if not ids:
        return None, "no 16S record found"

    time.sleep(NCBI_REQUEST_INTERVAL)

    try:
        fasta = _eutils_get(
            "efetch.fcgi",
            {**shared, "db": "nucleotide", "id": ids[0], "rettype": "fasta", "retmode": "text"},
        )
    except Exception as error:
        return None, f"NCBI fetch failed ({type(error).__name__})"

    lines = [line.strip() for line in fasta.splitlines() if line.strip()]

    if not lines or not lines[0].startswith(">"):
        return None, "unreadable FASTA from NCBI"

    accession = lines[0][1:].split()[0]
    sequence = "".join(lines[1:]).upper()

    if len(sequence) < SIXTEEN_S_MIN_LEN:
        return None, "record shorter than expected for 16S"

    return accession, sequence


# A binomial: capitalised genus, lower-case species epithet. Accepts an
# abbreviated genus ("P. intermedia") so those rows are not silently lost.
_BINOMIAL = re.compile(r"^(?:[A-Z][a-z]+|[A-Z]\.)\s+[a-z][a-z-]{2,}")

# Stricter form for the plain-text fallback: the whole field must be a
# binomial, optionally with one subspecies/strain qualifier.
_BINOMIAL_EXACT = re.compile(
    r"(?:[A-Z][a-z]+|[A-Z]\.)\s+[a-z][a-z-]{2,}"
    r"(?:\s+(?:subsp\.|var\.|str\.|sp\.)?\s*[A-Za-z0-9-]+)?"
)

_NOT_A_SPECIES = {
    "bacterium name", "bacterium", "name", "species", "candidate",
    "summary", "summary paragraph", "phylum", "confidence",
}


def collect_candidate_names(text, limit=12):
    """Species names from the first column of the model's candidate table.

    The model is inconsistent about table format -- sometimes pipe-delimited
    markdown, sometimes tab-separated -- so both are accepted. If neither
    parses, falls back to scanning for binomials at line starts, which still
    recovers the candidates from a plain-text table.
    """
    text = str(text or "")
    names = []

    def consider(cell):
        name = re.sub(r"[*_`]", "", str(cell)).strip().strip(".")

        if not name or name.lower() in _NOT_A_SPECIES:
            return
        if not _BINOMIAL.match(name):
            return
        if name not in names:
            names.append(name)

    for line in text.splitlines():
        stripped = line.strip()

        if not stripped:
            continue

        # Markdown pipe table.
        if stripped.count("|") >= 3:
            cells = [c.strip() for c in stripped.strip("|").split("|")]

            if len(cells) >= 3 and not set("".join(cells)) <= set("-: "):
                consider(cells[0])

            continue

        # Tab-separated table (what Streamlit's markdown renderer emits, and
        # what the model sometimes produces directly).
        if "\t" in stripped:
            cells = [c.strip() for c in stripped.split("\t")]

            if len(cells) >= 3:
                consider(cells[0])

            continue

    # Last resort: a table flattened to text with runs of spaces instead of
    # tabs. Requires both a column separator on the line and a first field
    # that is nothing but a binomial, so ordinary prose ("The analysis shows
    # that...") cannot be mistaken for a species.
    if not names:
        for line in text.splitlines():
            stripped = re.sub(r"[*_`]", "", line).strip()

            if not re.search(r"\s{2,}|\t|\|", stripped):
                continue

            first_field = re.split(r"\s{2,}|\t|\|", stripped)[0].strip()

            if _BINOMIAL_EXACT.fullmatch(first_field):
                consider(first_field)

    return names[:limit]


def build_16s_tree(sequences):
    """Neighbour-joining tree from pairwise identity. Returns (tree, newick)."""
    try:
        from Bio import Phylo
        from Bio.Align import PairwiseAligner
        from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
    except ImportError as error:
        raise RuntimeError(
            "Biopython is not installed in this environment. Add "
            "'biopython' to requirements.txt and reboot the app "
            f"({error})."
        ) from error

    labels = list(sequences)

    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 1
    aligner.mismatch_score = 0
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -0.5

    # Lower triangle, as DistanceMatrix expects.
    matrix = []

    for i, a in enumerate(labels):
        row = []

        for j, b in enumerate(labels[: i + 1]):
            if i == j:
                row.append(0.0)
                continue

            score = aligner.score(sequences[a], sequences[b])
            identity = score / max(len(sequences[a]), len(sequences[b]))
            row.append(round(max(0.0, 1.0 - identity), 6))

        matrix.append(row)

    tree = DistanceTreeConstructor().nj(DistanceMatrix(labels, matrix))
    tree.root_at_midpoint()

    # NJ can return small negative branch lengths; they are meaningless on a
    # drawn tree and make it unreadable, so clamp them.
    for clade in tree.find_clades():
        if clade.branch_length and clade.branch_length < 0:
            clade.branch_length = 0.0

    handle = io.StringIO()
    Phylo.write(tree, handle, "newick")

    return tree, handle.getvalue()


def draw_16s_tree(tree, n_taxa):
    from Bio import Phylo
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, axes = plt.subplots(figsize=(8, max(2.5, 0.55 * n_taxa)))
    Phylo.draw(
        tree,
        axes=axes,
        do_show=False,
        show_confidence=False,
        label_func=lambda c: c.name if c.is_terminal() else "",
    )
    axes.set_xlabel("substitutions per site (pairwise identity)")
    axes.spines[["top", "right", "left"]].set_visible(False)
    axes.set_ylabel("")
    axes.set_yticks([])
    figure.tight_layout()

    return figure


def render_16s_phylogeny(response_text, extra_organisms=(), email="", api_key=""):
    """Fetch 16S for each candidate, draw the tree, report what was missing.

    extra_organisms are names the user added after the analysis ran. They are
    treated exactly like candidates for fetching and tree-building, and are
    marked in the status table so it stays clear which came from the model.
    """
    names = collect_candidate_names(response_text)
    added = [str(o).strip() for o in extra_organisms if str(o).strip()]

    # Added organisms lead, so they are easy to find in the status table.
    names = added + [n for n in names if n.lower() not in (a.lower() for a in added)]

    if not names:
        st.warning(
            "No species names could be read from the candidate table, so no "
            "tree was built. The parser accepts pipe-delimited or "
            "tab-separated tables and lines beginning with a binomial; the "
            "response matched none of these."
        )

        with st.expander("Show the response the parser received", expanded=False):
            st.code(str(response_text)[:4000] or "(empty)")

        return


    sequences = {}
    status = []

    progress = st.progress(0.0, text="Retrieving 16S sequences from NCBI...")

    for index, name in enumerate(names, start=1):
        accession, result = fetch_16s(name, email=email, api_key=api_key)

        source = "added by you" if name in added else "model candidate"

        if accession:
            sequences[name] = result
            status.append({"Species": name, "Source": source,
                           "16S accession": accession,
                           "Length (bp)": len(result), "Status": "retrieved"})
        else:
            status.append({"Species": name, "Source": source,
                           "16S accession": "-",
                           "Length (bp)": "-", "Status": f"missing - {result}"})

        progress.progress(index / len(names), text=f"Retrieving 16S sequences... ({index}/{len(names)})")
        time.sleep(NCBI_REQUEST_INTERVAL)

    progress.empty()

    status_df = pd.DataFrame(status)
    st.dataframe(status_df, use_container_width=True)

    missing = status_df[status_df["Status"].str.startswith("missing")]

    if not missing.empty:
        st.warning(
            "16S missing for: " + ", ".join(missing["Species"]) +
            ". These are left out of the tree."
        )

    if len(sequences) < 3:
        st.info(
            f"{len(sequences)} sequence(s) retrieved - a tree needs at least 3 taxa, "
            "so none is drawn."
        )
        return

    with st.spinner(f"Aligning {len(sequences)} sequences and inferring tree..."):
        try:
            tree, newick = build_16s_tree(sequences)
        except Exception as error:
            st.error(f"Tree construction failed: {type(error).__name__}: {error}")
            return

        figure = draw_16s_tree(tree, len(sequences))

    st.pyplot(figure, use_container_width=True)

    st.caption(
        "Neighbour-joining tree from pairwise 16S identity, midpoint-rooted. "
        "Quick-look topology only - no multiple-sequence alignment, no model "
        "selection, no bootstrap support. For a publication figure, download "
        "the FASTA below and infer the tree with MAFFT plus IQ-TREE or RAxML."
    )

    fasta = "".join(f">{name.replace(' ', '_')}\n{seq}\n" for name, seq in sequences.items())

    column_left, column_right = st.columns(2)

    with column_left:
        st.download_button(
            "Download 16S FASTA",
            data=fasta,
            file_name="candidate_16S.fasta",
            mime="text/plain",
        )

    with column_right:
        st.download_button(
            "Download tree (Newick)",
            data=newick,
            file_name="candidate_16S_nj.nwk",
            mime="text/plain",
        )


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


def split_evidence_cell(cell):
    """Cited items from one evidence cell, with table citations and
    discriminating/universal qualifiers removed."""
    text = re.sub(r"[*_`]", "", str(cell))

    # "A (Table 1), B (Table 4)": a citation followed by a comma ends an
    # item, so treat it as a separator. Commas inside a name are untouched.
    text = re.sub(r"(\(\s*Tables?\s*\d[^)]*\))\s*,", r"\1;", text, flags=re.IGNORECASE)
    items = []

    for piece in re.split(r"\s*(?:;|\n|•)\s*", text):
        piece = _TABLE_CITATION.sub(" ", piece)
        piece = _QUALIFIER_PAREN.sub(" ", piece)
        piece = _QUALIFIER_TAIL.sub("", piece)
        piece = " ".join(piece.split()).strip(" .,:-–—\"'")
        piece = _drop_unbalanced_parens(piece).strip(" .,:-–—\"'")

        if piece and _norm_label(piece) not in ("none", "n/a", "na", "-"):
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
# the diagram, the familiar layout. Count anchors are computed for every n.
_VENN_HAND_PLACED = {
    2: {"titles": [(-0.85, 1.10, "center", "bottom"),
                   (0.85, 1.10, "center", "bottom")]},
    3: {"titles": [(-0.85, 1.42, "center", "bottom"),
                   (0.85, 1.42, "center", "bottom"),
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


_VENN_COUNT_SIZE = 10.5      # pt, count in a region with room for it
_VENN_MIN_COUNT_SIZE = 5.5   # pt, smallest count in a cramped region
INK_ZERO = "#a3a3a3"         # empty regions recede; populated ones carry the eye


def draw_evidence_venn(named_sets):
    """Draw an n-set Venn (2-5): the count of categories held by exactly each
    combination of candidates.

    Overlaps show what candidates share, the single-candidate regions what
    each holds alone. Empty regions are drawn in a light grey so the
    populated ones carry the eye. Returns (figure, {"svg", "pdf", "png"}).
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Ellipse

    names = list(named_sets)
    n_sets = len(names)
    geometry = _venn_geometry(n_sets)
    regions = venn_regions(named_sets)

    (bx_low, bx_high), (by_low, by_high) = geometry["bounds"]
    target_width_in = {2: 4.4, 3: 4.4, 4: 5.6, 5: 6.6}[n_sets]
    points_per_unit = 72 * target_width_in / (bx_high - bx_low + 1.6)
    unit = 1.0 / points_per_unit

    name_size, count_line_size = 10, 9
    title_height = (name_size * 1.25 + count_line_size * 1.3) * unit

    title_blocks = []
    for name, (x, y, ha, va) in zip(names, geometry["titles"]):
        count_text = f"n = {len(named_sets[name])}"
        width = max(len(abbreviate_species(name)) * 0.55 * name_size,
                    10 + len(count_text) * 0.55 * count_line_size) * unit
        left = x if ha == "left" else x - width if ha == "right" else x - width / 2
        bottom = y if va == "bottom" else y - title_height if va == "top" \
            else y - title_height / 2
        title_blocks.append((left, bottom, width, ha))

    x_low = min([bx_low] + [b[0] for b in title_blocks]) - 0.08
    x_high = max([bx_high] + [b[0] + b[2] for b in title_blocks]) + 0.08
    y_low = min([by_low] + [b[1] for b in title_blocks]) - 0.08
    y_high = max([by_high] + [b[1] + title_height for b in title_blocks]) + 0.08

    fill_alpha = 0.22 if n_sets <= 3 else 0.16

    with plt.rc_context(FIGURE_STYLE):
        figure = plt.figure(figsize=((x_high - x_low) * points_per_unit / 72,
                                     (y_high - y_low) * points_per_unit / 72))
        axes = figure.add_axes([0, 0, 1, 1])

        for index, (cx, cy, w, h, angle) in enumerate(geometry["ellipses"]):
            axes.add_patch(Ellipse((cx, cy), w, h, angle=angle,
                                   facecolor=VENN_COLORS[index], alpha=fill_alpha,
                                   edgecolor="none", zorder=1))
            axes.add_patch(Ellipse((cx, cy), w, h, angle=angle, facecolor="none",
                                   edgecolor=VENN_STROKES[index],
                                   linewidth=1.4 if n_sets <= 3 else 1.1,
                                   zorder=2))

        for mask, room in geometry["regions"].items():
            count = len(regions.get(mask, []))
            size = max(_VENN_MIN_COUNT_SIZE,
                       min(_VENN_COUNT_SIZE, room["radius"] * points_per_unit * 0.8))
            axes.text(room["x"], room["y"], str(count), ha="center",
                      va="center", fontsize=size,
                      fontweight="bold" if count else "normal",
                      color=INK if count else INK_ZERO, zorder=4)

        # Title block: italic name, then marker (shape + colour) and n.
        for index, (name, (left, bottom, width, ha)) in enumerate(
                zip(names, title_blocks)):
            top = bottom + title_height
            anchor_x = left if ha == "left" else left + width if ha == "right" \
                else left + width / 2
            axes.text(anchor_x, top, abbreviate_species(name), ha=ha, va="top",
                      fontstyle="italic", fontsize=name_size, color=INK, zorder=4)

            count_text = f"n = {len(named_sets[name])}"
            pair_width = (10 + len(count_text) * 0.55 * count_line_size) * unit
            pair_left = left if ha == "left" else left + width - pair_width \
                if ha == "right" else anchor_x - pair_width / 2
            line_top = top - name_size * 1.25 * unit

            axes.plot(pair_left + 3.5 * unit, line_top - 4.6 * unit,
                      marker=VENN_MARKERS[index], markersize=6,
                      markerfacecolor=VENN_COLORS[index],
                      markeredgecolor=VENN_STROKES[index],
                      markeredgewidth=0.9, linestyle="none", zorder=4)
            axes.text(pair_left + 10 * unit, line_top, count_text, ha="left",
                      va="top", fontsize=count_line_size, color=INK_MUTED,
                      zorder=4)

        axes.set_xlim(x_low, x_high)
        axes.set_ylim(y_low, y_high)
        axes.set_aspect("equal")
        axes.axis("off")

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
    numbered in that order.
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
    frame.insert(0, "No.", range(1, len(frame) + 1))
    return frame


def _fit_text(text, max_chars):
    text = " ".join(str(text).split())
    return text if len(text) <= max_chars else text[: max_chars - 1].rstrip() + "…"


def draw_evidence_heatmap(matrix, names):
    """Binary presence/absence heatmap: categories down, candidates across.

    Each row reads "No.  Category" then one cell per candidate and the
    "Shared by" fraction. Drawn in cell units (1 cell = 1 data unit) on an
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
    number_units = len(str(n_rows)) * 0.6 * row_font / ppu
    label_units = max(len(l) for l in labels) * 0.52 * row_font / ppu
    left_units = number_units + 0.45 + label_units + 1.0
    column_head_units = max(len(abbreviate_species(n)) for n in names) \
        * 0.52 * head_font / ppu * 0.87 + 1.0
    count_units = 2.2

    x_low, x_high = -left_units, n_cols + count_units
    y_top, y_bottom = -column_head_units, heading_rows + n_rows + 0.4
    number_right = x_low + 0.15 + number_units

    with plt.rc_context(FIGURE_STYLE):
        figure = plt.figure(figsize=((x_high - x_low) * cell_in,
                                     (y_bottom - y_top) * cell_in))
        axes = figure.add_axes([0, 0, 1, 1])

        for col, name in enumerate(names):
            axes.text(col + 0.5, -0.35, abbreviate_species(name), rotation=60,
                      rotation_mode="anchor", ha="left", va="bottom",
                      fontstyle="italic", fontsize=head_font, color=INK)

        axes.text(n_cols + 0.6, -0.35, "Shared by", rotation=60,
                  rotation_mode="anchor", ha="left", va="bottom",
                  fontsize=row_font, color=INK_MUTED)

        axes.text(x_low + 0.15, 0.62, "Categories", ha="left", va="center",
                  fontsize=head_font, fontweight="bold", color=INK)
        axes.plot([x_low + 0.15, n_cols + count_units - 0.2], [1.0, 1.0],
                  color=RULE, linewidth=0.5)

        for i, (_, row) in enumerate(matrix.iterrows()):
            y = heading_rows + i

            axes.text(number_right, y + 0.5, str(row["No."]), ha="right",
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

        axes.set_xlim(x_low, x_high)
        axes.set_ylim(y_bottom, y_top)
        axes.axis("off")

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
            figure, venn_exports = draw_evidence_venn(named_sets)
            st.pyplot(figure, use_container_width=True)
            plt.close(figure)

    st.caption(
        "Left: each filled cell is a category the candidate's evidence cell "
        "cites as shared with the query genome; rows are ordered by how many "
        "candidates share them, and candidates by the model's confidence. "
        "Right: the same categories as counts - overlaps are categories the "
        "candidates have in common, single-candidate regions are categories "
        "only that candidate cites. Shape size is fixed and not proportional "
        "to n."
    )

    csv = matrix[["No.", "Category", "Source"] + names + ["Cited by"]].copy()
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
                {"Region": title, "Category": filtered_label,
                 "Source": _table_tag(table_number)}
                for filtered_label, table_number in sorted(
                    (next(filtered[n][m] for n in venn_names if m in filtered[n])
                     for m in members),
                    key=lambda v: v[0].lower(),
                )
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
    "web6_added_organisms",
    "web6_tree_built",
    "web6_add_counter",
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

        # Held in session state so the result and its tree survive the reruns
        # Streamlit performs when a widget below is used. Without this, adding
        # an organism to the tree would re-trigger the model call.
        st.session_state["web6_response"] = full_response
        # Row labels of Tables 1-4, so the Venn section can match the model's
        # cited items back to verbatim rows without re-parsing the uploads.
        st.session_state["web6_table_rows"] = collect_table_rows(
            functional_category_df,
            kegg_pathway_df,
            protein_family_df,
            pathway_completeness_df,
        )
        st.session_state.setdefault("web6_added_organisms", [])
        # A new analysis means new candidates, so any tree already on screen
        # belongs to the previous run; require Generate again.
        st.session_state["web6_tree_built"] = False


# ---------------------------------------------------------------------------
# Result, Venn, tree, and post-hoc organism addition
#
# Rendered outside the button block and driven from session state, so the
# Venn filters and the "add an organism" control below can rerun on their
# own. Neither calls the model again; the tree step only queries NCBI and the
# already-fetched 16S sequences come from cache.
# ---------------------------------------------------------------------------

# A stored result belongs to the files that produced it. Closing the tool and
# reopening it, or clearing the uploaders, leaves no files but -- because the
# result lives in session state so the tree controls can rerun without
# re-querying the model -- would otherwise leave the previous run's output on
# screen as if it described the new (empty) input. Drop it instead.
if not (bakta_file and kbase_file and products_file):
    for stale_key in _RUN_STATE_KEYS:
        st.session_state.pop(stale_key, None)


if st.session_state.get("web6_response"):
    stored_response = st.session_state["web6_response"]
    added_organisms = st.session_state.get("web6_added_organisms", [])

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

    # Heatmap / UpSet / Venn of the evidence each candidate shares with the
    # query genome. Sits between the result and the tree: it is read off the
    # candidate table directly, needs no network call, and is shown as soon as
    # there is a result.
    st.divider()
    st.subheader("Shared evidence among candidate bacteria")
    render_shared_evidence(stored_response, st.session_state.get("web6_table_rows", {}))

    # The tree and everything that configures it live here, after the result,
    # because none of it is meaningful until there are candidate species to
    # tree. Nothing about the phylogeny appears in the pre-run form.
    st.divider()
    st.subheader("16S phylogeny of the candidate species")

    # The tree is built on request, not automatically. Species come from the
    # candidate table; the user can add their own organisms to that list; only
    # then does pressing Generate query NCBI. Changing the list clears the
    # tree, so what is on screen always matches the list above it and no
    # network call happens without an explicit click.
    candidate_names = collect_candidate_names(stored_response)
    tree_taxa = added_organisms + [
        n for n in candidate_names
        if n.lower() not in (a.lower() for a in added_organisms)
    ]

    if candidate_names:
        st.markdown("**Species read from the candidate table**")
        st.write(", ".join(candidate_names))
    else:
        st.warning(
            "No species names could be read from the candidate table. The "
            "parser accepts pipe-delimited or tab-separated tables and lines "
            "beginning with a binomial."
        )

    st.markdown("**Add your own organisms (optional)**")
    st.caption(
        "Added organisms are treated like candidates: one 16S sequence each "
        "from NCBI. The model is not queried again and the candidate table "
        "above does not change."
    )

    # The text box key carries a counter so that incrementing it after an add
    # gives a brand-new, empty widget. Setting a widget's value directly after
    # it has been created is not allowed, and st.form is not usable here
    # because the host app wraps Streamlit's callables.
    add_counter = st.session_state.setdefault("web6_add_counter", 0)
    add_column, add_button_column = st.columns([4, 1])

    with add_column:
        # Use the widget's RETURN VALUE, never a session_state lookup by key:
        # the host app wraps Streamlit's callables and may rewrite widget
        # keys, in which case st.session_state[<our key>] is empty even
        # though the box visibly contains text. The return value is correct
        # either way, and on the rerun the Add click triggers it still holds
        # what the user typed.
        typed_organism = st.text_input(
            "Organism name",
            placeholder="e.g. Tannerella forsythia",
            key=f"web6_organism_to_add_{add_counter}",
            label_visibility="collapsed",
        )

    with add_button_column:
        add_clicked = st.button(
            "Add", key="web6_add_organism", use_container_width=True
        )

    if add_clicked:
        name = " ".join(str(typed_organism or "").split())

        if not name:
            st.warning("Enter an organism name first.")
        elif name.lower() in (o.lower() for o in tree_taxa):
            st.info(f"{name} is already in the list.")
        else:
            st.session_state["web6_added_organisms"] = added_organisms + [name]
            st.session_state["web6_add_counter"] = add_counter + 1
            # The list changed, so any tree on screen is now out of date.
            st.session_state["web6_tree_built"] = False
            st.rerun()

    if added_organisms:
        st.success("Added: " + ", ".join(added_organisms))

        if st.button("Clear added organisms", key="web6_clear_added"):
            st.session_state["web6_added_organisms"] = []
            st.session_state["web6_tree_built"] = False
            st.rerun()

    with st.expander("NCBI options", expanded=False):
        ncbi_email = st.text_input(
            "Contact email (recommended)",
            value="",
            key="ncbi_email",
            help="NCBI asks programmatic users to identify themselves. "
                 "Sent only to NCBI, as the E-utilities 'email' parameter.",
        )
        ncbi_api_key = st.text_input(
            "API key (optional, raises the rate limit)",
            value="",
            type="password",
            key="ncbi_api_key",
        )

    tree_ready = st.session_state.get("web6_tree_built", False)

    st.caption(
        f"{len(tree_taxa)} species will be placed in the tree"
        + (" - at least 3 are needed." if len(tree_taxa) < 3 else ".")
    )

    if st.button(
        "Regenerate tree" if tree_ready else "Generate tree",
        type="primary",
        key="web6_generate_tree",
        disabled=len(tree_taxa) < 3,
    ):
        st.session_state["web6_tree_built"] = True
        st.rerun()

    if tree_ready:
        render_16s_phylogeny(
            stored_response,
            extra_organisms=added_organisms,
            email=ncbi_email,
            api_key=ncbi_api_key,
        )
    else:
        st.info(
            "Tree not generated yet. Add any organisms you want included, "
            "then press Generate tree."
        )
