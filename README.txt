GenoAnno Master App - Website API Settings Version

This version does NOT require config.txt.
Users enter their OpenAI API key and model settings directly on the website.

Folder contents to upload to GitHub:

apps/
master_app.py
README.txt
requirements.txt

Do not upload a config.txt file. It is not needed.

How the website works:

1. Open the Streamlit website.
2. Enter your own OpenAI API key in the API settings section.
3. Choose the model.
4. Choose temperature and max tokens.
5. Click "Use these settings".
6. Select a GenoAnno tool.
7. Click "Open this tool".
8. Upload or paste your biological data inside the selected tool and run the analysis.

Security note:

The OpenAI API key entered on the website is stored only in the current Streamlit browser session using st.session_state. It is not written to GitHub, not written to config.txt, and not saved permanently by this app. When the browser session ends or the user clicks "Clear API settings from this session", the key is removed from that session.

Optional fallback:

You may still use Streamlit Secrets or an OPENAI_API_KEY environment variable as a fallback, but this is not recommended for a public app because other users could consume your credits. The safest public setup is to require each user to enter their own API key on the website.

Run locally:

pip install -r requirements.txt
streamlit run master_app.py

Deploy on Streamlit Community Cloud:

1. Upload the unzipped app contents to GitHub.
2. Go to Streamlit Community Cloud.
3. Create a new app from your GitHub repository.
4. Set the main file path to:

master_app.py

5. Deploy.

No Streamlit Secrets are required for the normal website API-key workflow.
