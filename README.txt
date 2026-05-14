# GenoAnno

GenoAnno is a Streamlit web application for interpreting bacterial genome annotation outputs. It provides six tools for pathway expectation, metabolic pathway interpretation, KEGG pathway interpretation, protein-family interpretation, Bakta functional annotation interpretation, and oral-bacteria similarity analysis.

## Tools Included

1. Making Expectation  
   Predicts phenotype-relevant pathways from known bacterial behavior and pathway information.

2. Metabolic Pathway Interpreter  
   Interprets metabolic pathway tables and groups retained pathways into biologically meaningful phenotypes.

3. KEGG Pathway Interpreter  
   Counts KEGG pathway terms and summarizes dominant pathway-level functional signals.

4. Protein Family Interpreter  
   Counts protein-family terms and interprets their phenotype-level relevance.

5. Functional Interpreter  
   Processes Bakta product annotations, counts repeated functional categories, and generates biological summaries.

6. Similar Oral Bacteria Interpreter  
   Compares functional-category patterns with well-characterized oral bacteria.

## Required Project Files

The project should contain:

apps/
GenoAnno.py
README.txt
requirements.txt

## Streamlit Cloud Deployment

When deploying on Streamlit Community Cloud, set the main file path as:

GenoAnno.py

## How to Use the Website

1. Open the GenoAnno website.
2. Enter your OpenAI API key in the API settings section.
3. Select the OpenAI model.
4. Set the temperature and maximum token limit.
5. Click "Use these settings".
6. Select one analysis tool.
7. Click "Open this tool".
8. Upload or paste the required biological data.
9. Run the analysis.
10. Review, copy, or download the result.

## OpenAI API Key Handling

Each user enters their own OpenAI API key before running an analysis.

The key is stored only in the current Streamlit browser session using st.session_state. It is not written to project files, saved permanently, or shared across sessions.

The key is removed when the browser session ends or when the user clicks "Clear API settings from this session".

## Running Locally

Open a terminal in the project folder and run:

pip install -r requirements.txt
streamlit run GenoAnno.py

Then open the local Streamlit link in your browser and follow the same steps described above.
