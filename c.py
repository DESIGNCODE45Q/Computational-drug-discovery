import streamlit as st
from streamlit_option_menu import option_menu


# Import necessary modules
import subprocess
import os
import base64
import pickle
import pandas as pd
from PIL import Image
import requests
from stmol import showmol
import py3Dmol
import biotite.structure.io as bsio
from langchain.chains.summarize import load_summarize_chain
from langchain.document_loaders import PyPDFLoader
from langchain import OpenAI
import tempfile

# Set up OpenAI API
os.environ["OPENAI_API_KEY"] = "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"
llm = OpenAI(temperature=0)

# Define the function for rendering protein structure
def render_mol(pdb):
    pdbview = py3Dmol.view()
    pdbview.addModel(pdb, 'pdb')
    pdbview.setStyle({'cartoon': {'color': 'spectrum'}})
    pdbview.setBackgroundColor('white')
    pdbview.zoomTo()
    pdbview.zoom(2, 800)
    pdbview.spin(True)
    showmol(pdbview, height=500, width=800)
background_image_path = "C:\Projeects\Drug discovery\Platform-designed-to-accelerate-drug-discovery-development.jpg"
st.markdown(
    f"""
    <style>
    .reportview-container {{
        background: url('{background_image_path}');
        background-size: cover;
        background-repeat: no-repeat;
        opacity: 0.85;
    }}
    </style>
    """,
    unsafe_allow_html=True
)

# Define page selection
options = ["Protein Structure Prediction", "Chemical Bioactivity Prediction", "Document Summarizer", "Chatbot","About"]
selected_page = option_menu(None, options, default_index=0, orientation="horizontal")

st.markdown("""
    # Computational Machine Learning for predicting a Compound Bioactivity :A Drug Discovery Approach
    """)

# Increase the size of the options using CSS
st.markdown("""
    <style>
    .option-menu-container button {
        width: 150%;
    }
    </style>
""", unsafe_allow_html=True)


# Page 1: Protein Structure Prediction
if selected_page == "Protein Structure Prediction":
    # Protein sequence input
    DEFAULT_SEQ = "MGSSHHHHHHSSGLVPRGSHMRGPNPTAASLEASAGPFTVRSFTVSRPSGYGAGTVYYPTNAGGTVGAIAIVPGYTARQSSIKWWGPRLASHGFVVITIDTNSTLDQPSSRSSQQMAALRQVASLNGTSSSPIYGKVDTARMGVMGWSMGGGGSLISAANNPSLKAAAPQAPWDSSTNFSSVTVPTLIFACENDSIAPVNSSALPIYDSMSRNAKQFLEINGGSHSCANSGNSNQALIGKKGVAWMKRFMDNDTRYSTFACENPNSTRVSDFRTANCSLEDPAANKARKEAELAAATAEQ"
    txt = st.sidebar.text_area('Input sequence', DEFAULT_SEQ, height=275)
    
    # ESMfold
    def update(sequence=txt):
        headers = {
            'Content-Type': 'application/x-www-form-urlencoded',
        }
        response = requests.post('https://api.esmatlas.com/foldSequence/v1/pdb/', headers=headers, data=sequence, verify=False)

        name = sequence[:3] + sequence[-3:]
        pdb_string = response.content.decode('utf-8')

        with open('predicted.pdb', 'w') as f:
            f.write(pdb_string)

        struct = bsio.load_structure('predicted.pdb', extra_fields=["b_factor"])
        b_value = round(struct.b_factor.mean(), 4)

        # Display protein structure
        st.subheader('Visualization of predicted protein structure')
        render_mol(pdb_string)

        # plDDT value is stored in the B-factor field
        st.subheader('plDDT')
        st.write('plDDT is a per-residue estimate of the confidence in prediction on a scale from 0-100.')
        st.info(f'plDDT: {b_value}')

        st.download_button(
            label="Download PDB",
            data=pdb_string,
            file_name='predicted.pdb',
            mime='text/plain',
        )
    image = Image.open('3.png')

    st.image(image, use_column_width=True)
    predict = st.sidebar.button('Predict', on_click=update)

    if not predict:
        st.warning('👈 Enter protein sequence data')




        

# Page 2: Chemical Bioactivity Prediction
elif selected_page == "Chemical Bioactivity Prediction":
    # Function for descriptor calculation
    def desc_calc():
        # Perform the descriptor calculation
        bashCommand = "java -Xms2G -Xmx2G -Djava.awt.headless=true -jar ./PaDEL-Descriptor/PaDEL-Descriptor.jar -removesalt -standardizenitro -fingerprints -descriptortypes ./PaDEL-Descriptor/PubchemFingerprinter.xml -dir ./ -file descriptors_output.csv"
        process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
        output, error = process.communicate()
        os.remove('molecule.smi')

    # Function to download a CSV file
    def filedownload(df):
        csv = df.to_csv(index=False)
        b64 = base64.b64encode(csv.encode()).decode()
        href = f'<a href="data:file/csv;base64,{b64}" download="prediction.csv">Download Predictions</a>'
        return href

    # Function to build the model
    def build_model(input_data):
        # Read the saved regression model
        load_model = pickle.load(open('acetylcholinesterase_model.pkl', 'rb'))
        # Apply the model to make predictions
        prediction = load_model.predict(input_data)
        st.header('**Prediction output**')
        prediction_output = pd.Series(prediction, name='pIC50')
        molecule_name = pd.Series(load_data[1], name='molecule_name')
        df = pd.concat([molecule_name, prediction_output], axis=1)
        st.write(df)
        st.markdown(filedownload(df), unsafe_allow_html=True)

    # Logo image
    image = Image.open('1.png')

    st.image(image, use_column_width=True)

    # Page title
  
    # Sidebar
    with st.sidebar.header('1. Upload your CSV data'):
        uploaded_file = st.sidebar.file_uploader("Upload your input file", type=['txt'])
        st.sidebar.markdown("""
    [Example input file](https://raw.githubusercontent.com/dataprofessor/bioactivity-prediction-app/main/example_acetylcholinesterase.txt)
    """)

    if st.sidebar.button('Predict'):
        load_data = pd.read_table(uploaded_file, sep=' ', header=None)
        load_data.to_csv('molecule.smi', sep='\t', header=False, index=False)

        st.header('**Original input data**')
        st.write(load_data)

        with st.spinner("Calculating descriptors..."):
            desc_calc()

        # Read in calculated descriptors and display the dataframe
        st.header('**Calculated molecular descriptors**')
        desc = pd.read_csv('descriptors_output.csv')
        st.write(desc)
        st.write(desc.shape)

        # Read descriptor list used in previously built model
        st.header('**Subset of descriptors from previously built models**')
        Xlist = list(pd.read_csv('descriptor_list.csv').columns)
        desc_subset = desc[Xlist]
        st.write(desc_subset)
        st.write(desc_subset.shape)

        # Apply trained model to make predictions on query compounds
        build_model(desc_subset)
    else:
        st.info('Upload input data in the sidebar to start!')

# Page 3: PDF Summarizer
elif selected_page == "Document Summarizer":
    def summarize_pdfs_from_folder(pdfs_folder):
        summaries = []
        for pdf_file in pdfs_folder:
            with tempfile.NamedTemporaryFile(delete=False) as temp_file:
                temp_path = temp_file.name
                temp_file.write(pdf_file.read())

            loader = PyPDFLoader(temp_path)
            docs = loader.load_and_split()
            chain = load_summarize_chain(llm, chain_type="map_reduce")
            summary = chain.run(docs)
            summaries.append(summary)

            # Delete the temporary file
            os.remove(temp_path)

        return summaries

    # Streamlit App
    st.title("Document Summarizer")

    # Allow user to upload PDF files
     # Logo image
    image = Image.open('2.png')

    st.image(image, use_column_width=True)

    pdf_files = st.file_uploader("Upload PDF files", type="pdf", accept_multiple_files=True)

    if pdf_files:
        # Generate summaries when the "Generate Summary" button is clicked
        if st.button("Generate Summary"):
            st.write("Summaries:")
            summaries = summarize_pdfs_from_folder(pdf_files)
            for i, summary in enumerate(summaries):
                st.write(f"Summary for PDF {i+1}:")
                st.write(summary)

# You can add more pages or features as needed
elif selected_page == "Chatbot":
    # Page 4: Chatbot
    # Chatbot code
    from openai import OpenAI
    import streamlit as st

    # Provided OpenAI API key
    openai_api_key = "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"

    st.title("💬 Chatbot")
    st.caption("🚀 A Streamlit chatbot powered by OpenAI LLM")

    # Initialize chat history if not present
    if "messages" not in st.session_state:
        st.session_state["messages"] = [{"role": "assistant", "content": "How can I help you?"}]

    # Display chat history
    for msg in st.session_state.messages:
        st.chat_message(msg["role"]).write(msg["content"])

    # Accept user input
    if prompt := st.chat_input():
        # Initialize OpenAI client
        client = OpenAI(api_key=openai_api_key)

        # Add user input to chat history
        st.session_state.messages.append({"role": "user", "content": prompt})
        st.chat_message("user").write(prompt)

        # Generate response from OpenAI
        response = client.chat.completions.create(model="gpt-3.5-turbo", messages=st.session_state.messages)
        msg = response.choices[0].message.content

        # Add response to chat history
        st.session_state.messages.append({"role": "assistant", "content": msg})
        st.chat_message("assistant").write(msg)

elif selected_page == "About":
    st.markdown("""
        <h1 style='font-size: 36px;'>About</h1>
        <p style='font-size: 28px;'>
            "Computational Machine Learning for Predicting Compound Bioactivity: 
            A Drug Discovery Approach" is an integrated bioinformatics project focused on data collection, 
            preprocessing, exploratory analysis, and model building using bioactivity data from the ChEMBL database. It involves building regression models with random forest algorithms, comparing regressors for performance, and developing a Streamlit web application for bioactivity prediction, protein structure prediction, document summarization, and chatbot integration. The project emphasizes practical application and user interaction, providing a comprehensive platform for researchers in drug discovery.
        </p>
    """, unsafe_allow_html=True)