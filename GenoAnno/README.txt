GenoAnno Master App - Ready for Streamlit Cloud

How to run locally
------------------
1. Open a terminal in this folder.
2. Install dependencies:

   pip install -r requirements.txt

3. For local testing, either add your OpenAI key to config.txt:

   KEY=your_openai_api_key_here

   OR set an environment variable named OPENAI_API_KEY.

4. Run:

   streamlit run master_app.py


How to publish on Streamlit Community Cloud
-------------------------------------------
1. Upload this full folder to GitHub.
2. Keep config.txt with KEY blank. Do not commit your real OpenAI key.
3. In Streamlit Cloud, deploy with:

   Main file path: master_app.py

4. In Streamlit Cloud > App settings > Secrets, add:

   KEY = "your_openai_api_key_here"
   DEFAULT_MODEL = "gpt-4o-mini"
   TEMPERATURE = "0.5"
   MAX_TOKENS = "2000"

The child apps now read keys in this order:
1. Streamlit Secrets
2. config.txt
3. OPENAI_API_KEY environment variable

Website behavior
----------------
- No app opens automatically.
- Click Open this tool to launch one selected analysis app.
- The six tools are kept as independent modules inside the apps folder.
