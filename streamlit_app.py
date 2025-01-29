import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw
from stmol import showmol
import py3Dmol
from pathlib import Path
import pandas as pd
import os
import numpy as np



# Load data
@st.cache_data()
def download_data():
    "Download the ChEMBL database"
    current_file = Path(os.path.abspath(''))
    #csv_file = current_file.parent / "CF3_seed_stats.csv"
    csv_file = "CF3_seed_stats.csv"
    df = pd.read_csv(csv_file).sort_values("pooled_CF3_De", ascending = True).reset_index(drop = True)
    # Sort by CF3_De
    return df


# Visualize molecule 
def visualize_3D(xyzstring, w, h, ts=False, no2 = False):
    "Visualize the molecule in 3D using stmol"
    xyzview = py3Dmol.view(width=w,height=w)
    xyzview.addModel(xyzstring,'xyz')
    xyzview.setStyle(
        {'sphere': {'colorscheme': 'cyanCarbon', 'scale': 0.25}, 'stick': {'colorscheme': 'cyanCarbon'}}
    )

    # Set style for atoms 1-35 (gray carbons with alpha of 0.5)
    # Default style for all atoms (cyan)
    xyzview.setStyle(
        {},  # Empty selector applies to all atoms
        {'sphere': {'colorscheme': 'cyanCarbon', 'scale': 0.25}, 'stick': {'colorscheme': 'cyanCarbon'}}
    )
    if no2:
        end = 35
    else:
        end = 36
    if ts:    # Custom style for atoms 1-35 (gray with alpha of 0.5)
        for atom_idx in range(1, end):  # Atoms 1-35 (py3Dmol uses 1-based indexing)
            xyzview.setStyle(
                {'serial': atom_idx},  # Target specific atom by serial number
                {
                    'sphere': {'colorscheme': 'grayCarbon', 'scale': 0.25, 'alpha': 0.5},
                    'stick': {'colorscheme': 'grayCarbon', 'alpha': 0.5},
                }
            )
    xyzview.zoomTo()
    xyzview.setBackgroundColor('white')
    showmol(xyzview, height = w,width=h)


# APP
# Title
st.set_page_config(layout="wide")
st.title('Saturn Results')
st.caption("Visualize top generated molecules")



# Download data, calculate descriptors and fingerprints
df = download_data()
# Sidebar

molecule_index = st.sidebar.selectbox("Molecule Index", range(len(df)))
smiles = df.iloc[molecule_index, :]["smiles"]



cola, colb = st.columns([0.5,1])

with cola:
    st.header(f"Molecule {molecule_index}")
    de = df.loc[molecule_index, "pooled_CF3_De"]
    NO2_de = df.loc[molecule_index, "pooled_NO2_De"]
    df['CF3_De_rank'] = df['pooled_CF3_De'].rank(ascending=True)
    df['NO2_De_rank'] = df['pooled_NO2_De'].rank(ascending=True)
    natoms = Chem.MolFromSmiles(smiles).GetNumHeavyAtoms()
    sa_score = df.loc[molecule_index, "mbh_sa_raw_values"]

    def rank_color(rank, total_count):
        if rank <= 10:
            return 'green'
        elif rank <= 50:
            return 'orange'
        elif rank > total_count - 50:
            return 'red'
        else:
            return 'black'


    # Display metrics in rows
    metric_row1 = st.columns(2)
    metric_row2 = st.columns(2)

    # Row 1: CF3 ΔE and NO2 ΔE
    with metric_row1[0]:
        rank_color_de = rank_color(int(df.loc[molecule_index, 'CF3_De_rank']), len(df))
        st.metric("CF3 ΔE", round(de, 2))
        # Custom rank color using markdown
        st.markdown(f"<span style='color:{rank_color_de}; font-weight:bold;'>Rank: {int(df.loc[molecule_index, 'CF3_De_rank'])}</span>", unsafe_allow_html=True)

    with metric_row1[1]:
        rank_color_NO2_de = rank_color(int(df.loc[molecule_index, 'NO2_De_rank']), len(df))
        st.metric("NO2 ΔE", round(NO2_de, 2))
        # Custom rank color using markdown
        st.markdown(f"<span style='color:{rank_color_NO2_de}; font-weight:bold;'>Rank: {int(df.loc[molecule_index, 'NO2_De_rank'])}</span>", unsafe_allow_html=True)

    # Row 2
    with metric_row2[0]:
        st.metric("SA Score", round(sa_score, 2))
    with metric_row2[1]:
        st.metric("Number of heavy atoms", natoms)



with colb:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        st.error('INVALID MOLECULE', icon="🚨")
    st.image(Draw.MolToImage(mol, size=(500,500)), width = 400)

with st.expander("CF3 Geometries", expanded = True):
    # 3D coordinates
    col1, col2 = st.columns([1,1])

    with col1:
        st.header("Best TS conformer")
        #st.caption("Best TS conformer")
        tsxyz = df.loc[molecule_index, "CF3_BestTSXYZ"]
        if tsxyz is not None:
            visualize_3D(tsxyz, 500, 500, True)
        else:
            st.error('INVALID MOLECULE', icon="🚨")


    with col2:
        st.header("Best catalyst conformer")
        #st.caption("Best catalyst conformer")
        catxyz = df.loc[molecule_index, "CF3_BestCATXYZ"]
        if catxyz is not None:
            visualize_3D(catxyz, 500, 500)
        else:
            st.error('INVALID MOLECULE', icon="🚨")

# Add a tab below for NO2
with st.expander("NO2 Geometries", expanded = False):
    tab1, tab2 = st.columns([1, 1])

    with tab1:
        st.header("Best TS conformer (NO2)")
        tsxyz_no2 = df.loc[molecule_index, "NO2_BestTSXYZ"]
        if tsxyz_no2 is not None:
            visualize_3D(tsxyz_no2, 500, 500, True, True)
        else:
            st.error('INVALID MOLECULE', icon="🚨")

    # Best Catalyst Conformer for NO2
    with tab2:
        st.header("Best catalyst conformer (NO2)")
        catxyz_no2 = df.loc[molecule_index, "NO2_BestCATXYZ"]
        if catxyz_no2 is not None:
            visualize_3D(catxyz_no2, 500, 500)
        else:
            st.error('INVALID MOLECULE', icon="🚨")
