import streamlit as st

st.set_page_config(
    page_title="Drug Discovery Web",
    layout="wide"
)
 
home_page = st.Page("pages/home.py", title = "Home")

pg = st.navigation(
    [
        home_page,
    ]
)
pg.run()