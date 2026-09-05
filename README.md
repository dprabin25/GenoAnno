# GenoAnno



**GenoAnno **is a Streamlit application that turns bacterial genome annotation output into phenotype-level biological interpretation. It targets researchers who already have annotation tables from pipelines such as **KBase-DRAM**, **Proksee**, or **Bakta** and want a fast, no-code read on what those genes imply about a bacterium's phenotype: pathway completeness, functional-category composition, and similarity to well-characterized oral bacteria without writing analysis scripts.

The app itself does no biological inference on its own: it parses and summarizes the uploaded annotation table, then sends that structured summary to an OpenAI model (chosen and paid for by the user) to generate the biological interpretation.

**Live app:** [![Live app](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://genoanno.streamlit.app/)

## How it works

`GenoAnno.py` is a single entry point — a dashboard that lists all six tools. Selecting a tool and clicking **Open this tool** loads that tool's module from [`apps/`](https://github.com/dprabin25/GenoAnno/tree/main/apps) in an isolated namespace, so switching tools mid-session never leaks widget state or cached results between them. Each tool follows the same underlying pattern:

1. Take in the annotation data (typed/pasted for tool 01, uploaded as a file for tools 02–06 — see the input method below).
2. Count and group the relevant terms (pathways, KEGG IDs, protein families, or functional categories).
3. Build a structured prompt from those counts and send it to the configured OpenAI model.
4. Render the model's interpretation alongside the underlying counts table.

See the six tool implementations directly in [`apps/`](https://github.com/dprabin25/GenoAnno/tree/main/apps).



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

<img width="777" height="99" alt="image" src="https://github.com/user-attachments/assets/1bacf740-ad4d-49b4-b0f4-a7e5536ac4a4" />


3. Select an OpenAI model, temperature, and maximum token limit, then click **Use these settings**.

 <img width="614" height="93" alt="image" src="https://github.com/user-attachments/assets/b56812a0-7720-41fa-b69b-2f5fe04b3687" />


4. Select a tool and click **Open this tool**.


E.g. 
   <img width="750" height="199" alt="image" src="https://github.com/user-attachments/assets/5cd1887b-8937-4a61-89fc-d7c5205352d7" />


5. Scroll down and provide your data:
   - **Making Expectation (01):** the form loads with an editable example already filled in — type over it or paste your own pathway list and phenotype description.
   - **All other tools (02–06):** upload your annotation file. There is no paste option for these — a file is required.

E.g. 

<img width="379" height="117" alt="image" src="https://github.com/user-attachments/assets/9dec0f22-7b11-46e8-9cb2-dfa17d8a48c2" />


6. Run the analysis, then review, copy, or download the result.

E.g.

<img width="379" height="120" alt="image" src="https://github.com/user-attachments/assets/7075dd43-344c-4a06-b3ae-860bfc4d6046" />

   
## Notes

### Model settings

The dashboard offers `gpt-4o-mini`, `gpt-4.1-mini`, `gpt-4.1`, `gpt-4o`, `gpt-5-mini`, `gpt-5`, or a custom model name. Reasoning-style models (the `gpt-5` family) don't accept a fixed `temperature` or `max_tokens` on the chat-completions endpoint — GenoAnno detects this automatically and retries the request with the parameters those models do accept, so no manual adjustment is needed when switching models.

You are billed directly by OpenAI for any requests GenoAnno makes on your behalf; GenoAnno itself does not add a markup or fee.

### OpenAI API key handling

Each user supplies their own OpenAI API key before running an analysis.

- The key lives only in the current Streamlit session state — it is never written to disk, saved permanently, or shared across sessions or users.
- It is cleared automatically when the browser session ends, or manually via **Clear API settings from this session**.


## Outputs

Each tool renders the model's biological interpretation alongside the underlying counts table, with options to copy the text or download the result.

E.g. 

<img width="2351" height="689" alt="image" src="https://github.com/user-attachments/assets/5751fc71-7eec-4080-b81c-e42b8fd0022d" />

## Example Inputs for 02-06 apps

- Provided in this Github repo for _Veillonella parvula_

## Reference

[1] Connor Pasvantis, Preston Wgunn,.... Prabin Dawadi, Sayaka Miura. Large language model assistance for interpreting oral microbiota phenotype (2026). Under Review.

## License

Copyright 2026, Authors and University of Mississippi

BSD 3-Clause "New" or "Revised" License, which is a permissive license similar to the BSD 2-Clause License except that it prohibits others from using the name of the project or its contributors to promote derived products without written consent.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
