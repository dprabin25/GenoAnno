# GenoAnno

GenoAnno is a Streamlit application for interpreting bacterial genome annotation output. It turns annotation tables from tools such as KBase-DRAM, Proksee, and Bakta into phenotype-level biological summaries, using an OpenAI model to do the interpretation.

**Live app:** [genoanno.streamlit.app](https://genoanno.streamlit.app/)

## Tools

| # | Tool | Purpose | Expected input |
|---|------|---------|-----------------|
| 01 | Making Expectation | Predicts phenotype-relevant active pathways from known bacterial behavior and pathway lists. | Pathway/behavior list |
| 02 | Metabolic Pathway Interpreter | Filters annotated metabolic pathway tables and groups retained pathways into biological phenotypes. | `products.tsv` (KBase-DRAM) |
| 03 | KEGG Pathway Interpreter | Counts repeated KEGG pathway terms and interprets dominant pathway-level functional signals. | `annotation.tsv` (KBase-DRAM) |
| 04 | Protein Family Interpreter | Summarizes repeated protein-family terms into phenotype-level functional interpretations. | `annotation.tsv` (KBase-DRAM) |
| 05 | Functional Interpreter | Parses product annotations, counts functional categories, and generates phenotype summaries. | `annotation.tsv` (Proksee/Bakta) |
| 06 | Similar Oral Bacteria Interpreter | Compares functional-category patterns against well-characterized oral bacteria. | `annotation.tsv` (Proksee/Bakta) |

## Using the web app

1. Open the [GenoAnno app](https://genoanno.streamlit.app/).
2. Enter your OpenAI API key in the API settings section.
3. Select the OpenAI model, temperature, and maximum token limit.
4. Click **Use these settings**.
5. Select a tool and click **Open this tool**.
6. Upload or paste the required annotation data.
7. Run the analysis, then review, copy, or download the result.

## OpenAI API key handling

Each user supplies their own OpenAI API key before running an analysis.

- The key lives only in the current Streamlit session state — it is never written to disk, saved permanently, or shared across sessions.
- It is cleared automatically when the browser session ends, or manually via **Clear API settings from this session**.

## Running locally

```bash
git clone https://github.com/dprabin25/GenoAnno.git
cd GenoAnno
pip install -r requirements.txt
streamlit run GenoAnno.py
```

Open the local Streamlit URL printed in the terminal and follow the same steps as the web version.

To avoid re-entering an API key every run, set defaults in `config.txt`:

```
KEY=
DEFAULT_MODEL=gpt-4o-mini
TEMPERATURE=0.5
MAX_TOKENS=2000
```

Leave `KEY` blank for local use and enter it in the app instead, or paste it here for convenience. On Streamlit Community Cloud, keep `KEY` blank and set it under **App settings → Secrets**.

## Project structure

```
GenoAnno/
├── GenoAnno.py          # Master app — dashboard and tool launcher
├── config.txt           # Default model/temperature/token settings
├── requirements.txt
└── apps/
    ├── web1_pathway_expectation.py
    ├── web2_phenotype_grouping.py
    ├── web3_kegg_pathway_count.py
    ├── web4_protein_family_count.py
    ├── web5_bakta_functional_category.py
    └── web6_similar_oral_bacteria.py
```

## Requirements

- Python 3.9+
- [streamlit](https://streamlit.io/)
- [openai](https://pypi.org/project/openai/)
- [pandas](https://pandas.pydata.org/)

See `requirements.txt` for exact packages.
