# GenoAnno

[![Live app](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://genoanno.streamlit.app/)

GenoAnno is a Streamlit application that turns bacterial genome annotation output into phenotype-level biological interpretation. It targets researchers who already have annotation tables from pipelines such as **KBase-DRAM**, **Proksee**, or **Bakta** and want a fast, no-code read on what those genes imply about a bacterium's phenotype — pathway completeness, functional-category composition, and similarity to well-characterized oral bacteria — without writing analysis scripts by hand.

The app itself does no biological inference on its own: it parses and summarizes the uploaded annotation table, then sends that structured summary to an OpenAI model (chosen and paid for by the user) to generate the biological interpretation.

**Live app:** [genoanno.streamlit.app](https://genoanno.streamlit.app/)

## How it works

`GenoAnno.py` is a single entry point — a dashboard that lists all six tools. Selecting a tool and clicking **Open this tool** loads that tool's module from `apps/` in an isolated namespace, so switching tools mid-session never leaks widget state or cached results between them. Each tool follows the same underlying pattern:

1. Parse the uploaded annotation table (tab- or comma-separated).
2. Count and group the relevant terms (pathways, KEGG IDs, protein families, or functional categories).
3. Build a structured prompt from those counts and send it to the configured OpenAI model.
4. Render the model's interpretation alongside the underlying counts table.

## Tools

| # | Tool | Purpose | Expected input |
|---|------|---------|-----------------|
| 01 | Making Expectation | Predicts phenotype-relevant active pathways from known bacterial behavior and pathway lists. | Pathway/behavior list |
| 02 | Metabolic Pathway Interpreter | Filters annotated metabolic pathway tables and groups retained pathways into biological phenotypes. | `products.tsv` (KBase-DRAM) |
| 03 | KEGG Pathway Interpreter | Counts repeated KEGG pathway terms and interprets dominant pathway-level functional signals. | `annotation.tsv` (KBase-DRAM) |
| 04 | Protein Family Interpreter | Summarizes repeated protein-family terms into phenotype-level functional interpretations. | `annotation.tsv` (KBase-DRAM) |
| 05 | Functional Interpreter | Parses product annotations, counts functional categories, and generates phenotype summaries. | `annotation.tsv` (Proksee/Bakta) |
| 06 | Similar Oral Bacteria Interpreter | Compares functional-category patterns against well-characterized oral bacteria, ranked by confidence. | `annotation.tsv` (Proksee/Bakta) |

## Using the web app

1. Open the [GenoAnno app](https://genoanno.streamlit.app/).
2. Enter your OpenAI API key in the API settings section.

<img width="348" height="81" alt="image" src="https://github.com/user-attachments/assets/931501ef-48a2-447d-ad68-3ab60663a3bf" />



3. Select an OpenAI model, temperature, and maximum token limit.
4. Click **Use these settings**.
   
<img width="589" height="58" alt="image" src="https://github.com/user-attachments/assets/ca3e2791-d045-4d0e-96ca-effd0ab8e82d" />

5. Select a tool and click **Open this tool**.

E.g. 

<img width="749" height="197" alt="image" src="https://github.com/user-attachments/assets/eb14cf83-7fee-4515-b8fd-11a852436922" />

7. Upload or paste the required annotation data.
8. Run the analysis, then review, copy, or download the result.

E.g. 
<img width="383" height="97" alt="image" src="https://github.com/user-attachments/assets/32e2ea5d-e538-4c1c-b859-1c7ca87590af" />


### Model settings

The dashboard offers `gpt-4o-mini`, `gpt-4.1-mini`, `gpt-4.1`, `gpt-4o`, `gpt-5-mini`, `gpt-5`, or a custom model name. Reasoning-style models (the `gpt-5` family) don't accept a fixed `temperature` or `max_tokens` on the chat-completions endpoint — GenoAnno detects this automatically and retries the request with the parameters those models do accept, so no manual adjustment is needed when switching models.

You are billed directly by OpenAI for any requests GenoAnno makes on your behalf; GenoAnno itself does not add a markup or fee.

## OpenAI API key handling

Each user supplies their own OpenAI API key before running an analysis.

- The key lives only in the current Streamlit session state — it is never written to disk, saved permanently, or shared across sessions or users.
- It is cleared automatically when the browser session ends, or manually via **Clear API settings from this session**.



## Outputs

E.g.

<img width="2408" height="1066" alt="image" src="https://github.com/user-attachments/assets/1008d3c1-0917-4550-b8a0-0143edebd143" />




