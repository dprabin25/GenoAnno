# GenoAnno

[![Live app](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://genoanno.streamlit.app/)

GenoAnno is a Streamlit application that turns bacterial genome annotation output into phenotype-level biological interpretation. It targets researchers who already have annotation tables from pipelines such as **KBase-DRAM**, **Proksee**, or **Bakta** and want a fast, no-code read on what those genes imply about a bacterium's phenotype — pathway completeness, functional-category composition, and similarity to well-characterized oral bacteria — without writing analysis scripts by hand.

The app itself does no biological inference on its own: it parses and summarizes the uploaded annotation table, then sends that structured summary to an OpenAI model (chosen and paid for by the user) to generate the biological interpretation.

**Live app:** [genoanno.streamlit.app](https://genoanno.streamlit.app/)

## How it works

`GenoAnno.py` is a single entry point — a dashboard that lists all six tools. Selecting a tool and clicking **Open this tool** loads that tool's module from `apps/` in an isolated namespace, so switching tools mid-session never leaks widget state or cached results between them. Each tool follows the same underlying pattern:

1. Take in the annotation data (typed/pasted for tool 01, uploaded as a file for tools 02–06 — see the input method below).
2. Count and group the relevant terms (pathways, KEGG IDs, protein families, or functional categories).
3. Build a structured prompt from those counts and send it to the configured OpenAI model.
4. Render the model's interpretation alongside the underlying counts table.

## Tools

| # | Tool | Purpose | Input method | Expected input |
|---|------|---------|---------------|-----------------|
| 01 | Making Expectation | Predicts phenotype-relevant active pathways from known bacterial behavior and pathway lists. | Type or paste — form fields come pre-filled with an editable example; replace the sample text with your own data. | Pathway list, bacteria/phenotype description |
| 02 | Metabolic Pathway Interpreter | Filters annotated metabolic pathway tables and groups retained pathways into biological phenotypes. | Upload only — no paste option. | `products.tsv` (KBase-DRAM) |
| 03 | KEGG Pathway Interpreter | Counts repeated KEGG pathway terms and interprets dominant pathway-level functional signals. | Upload only — no paste option. | `annotation.tsv` (KBase-DRAM) |
| 04 | Protein Family Interpreter | Summarizes repeated protein-family terms into phenotype-level functional interpretations. | Upload only — no paste option. | `annotation.tsv` (KBase-DRAM) |
| 05 | Functional Interpreter | Parses product annotations, counts functional categories, and generates phenotype summaries. | Upload only — no paste option. | `annotation.tsv` (Proksee/Bakta) |
| 06 | Similar Oral Bacteria Interpreter | Compares functional-category patterns against well-characterized oral bacteria, ranked by confidence. | Upload only — no paste option. | `annotation.tsv` (Proksee/Bakta) |

## Using the web app

1. Open the [GenoAnno app](https://genoanno.streamlit.app/).
2. Enter your OpenAI API key in the API settings section.

   *(screenshot)*

3. Select an OpenAI model, temperature, and maximum token limit, then click **Use these settings**.

   *(screenshot)*

4. Select a tool and click **Open this tool**.

   *(screenshot)*

5. Provide your data:
   - **Making Expectation (01):** the form loads with an editable example already filled in — type over it or paste your own pathway list and phenotype description.
   - **All other tools (02–06):** upload your annotation file. There is no paste option for these — a file is required.
6. Run the analysis, then review, copy, or download the result.

   *(screenshot)*

## Model settings

The dashboard offers `gpt-4o-mini`, `gpt-4.1-mini`, `gpt-4.1`, `gpt-4o`, `gpt-5-mini`, `gpt-5`, or a custom model name. Reasoning-style models (the `gpt-5` family) don't accept a fixed `temperature` or `max_tokens` on the chat-completions endpoint — GenoAnno detects this automatically and retries the request with the parameters those models do accept, so no manual adjustment is needed when switching models.

You are billed directly by OpenAI for any requests GenoAnno makes on your behalf; GenoAnno itself does not add a markup or fee.

## OpenAI API key handling

Each user supplies their own OpenAI API key before running an analysis.

- The key lives only in the current Streamlit session state — it is never written to disk, saved permanently, or shared across sessions or users.
- It is cleared automatically when the browser session ends, or manually via **Clear API settings from this session**.

## Outputs

Each tool renders the model's biological interpretation alongside the underlying counts table, with options to copy the text or download the result.

*(screenshot)*

## Troubleshooting

- **"No OpenAI API key found in this session"** — enter a key and click **Use these settings** before opening a tool; the key doesn't carry over between tools until it's saved at the dashboard level.
- **Analysis fails immediately after running** — the OpenAI request likely errored (invalid/expired key, no access to the selected model, or a rate limit); the error message returned by OpenAI is shown in place of the interpretation.
- **Empty or all-zero counts** — check that the uploaded file matches the expected input for that tool (see the table above) and uses the expected column layout and delimiter (tab-separated `.tsv` for most tools).
- **"Upload" doesn't accept pasted text** — that's expected for tools 02–06; only Making Expectation (01) accepts typed/pasted input.

## Project structure

```
GenoAnno/
├── GenoAnno.py          # Master app — dashboard and tool launcher
└── apps/
    ├── web1_pathway_expectation.py
    ├── web2_phenotype_grouping.py
    ├── web3_kegg_pathway_count.py
    ├── web4_protein_family_count.py
    ├── web5_bakta_functional_category.py
    └── web6_similar_oral_bacteria.py
```
