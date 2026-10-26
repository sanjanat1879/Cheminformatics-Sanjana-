"""
Chemical Compound Dashboard
---------------------------
Compounds: Benzene, Naphthalene, Salicylic acid, Ethanol, Pyrimidine, Aspirin

Install:   pip install streamlit rdkit pillow
Run:       streamlit run compound_dashboard.py

(The 3D viewer loads 3Dmol.js from a CDN, so an internet connection is needed for it.)
"""

import json

import streamlit as st
import streamlit.components.v1 as components
from rdkit import Chem
from rdkit.Chem import AllChem, Draw

# ----------------------------------------------------------------------------
# DATABASE
# ----------------------------------------------------------------------------
# Each question: (question text, [options], index of correct option)
COMPOUNDS = {
    "Benzene": {
        "smiles": "c1ccccc1",
        "iupac": "Benzene",
        "physical": {
            "Molecular formula": "C6H6",
            "Molar mass": "78.11 g/mol",
            "Appearance": "Colorless, highly flammable liquid",
            "Odor": "Sweet, aromatic (gasoline-like)",
            "Melting point": "5.5 °C",
            "Boiling point": "80.1 °C",
            "Density": "0.8765 g/cm³ (20 °C)",
            "Solubility in water": "1.79 g/L (25 °C); miscible with most organic solvents",
            "Vapor pressure": "12.7 kPa (25 °C)",
            "Refractive index (nD)": "1.5011",
            "Flash point": "−11 °C",
            "Dipole moment": "0 D (non-polar, symmetrical)",
            "log P (octanol/water)": "2.13",
            "Structure": "Planar hexagon, all C–C bonds 139 pm",
        },
        "chemical": [
            "Aromatic: planar, cyclic, fully conjugated with 6 π electrons (Hückel 4n+2, n = 1), giving exceptional stability (~150 kJ/mol resonance energy).",
            "Characteristic reaction is electrophilic aromatic substitution (EAS) rather than addition.",
            "Nitration (conc. HNO3/H2SO4) gives nitrobenzene; sulfonation (fuming H2SO4) gives benzenesulfonic acid.",
            "Halogenation (Cl2 or Br2 with FeX3/AlX3) gives halobenzenes.",
            "Friedel–Crafts alkylation and acylation (RCl or RCOCl with AlCl3) give alkyl- and acyl-benzenes.",
            "Catalytic hydrogenation (H2, Ni/Pt, high T and P) gives cyclohexane.",
            "Burns with a sooty, luminous flame: 2 C6H6 + 15 O2 → 12 CO2 + 6 H2O.",
            "Does not decolorize bromine water or alkaline KMnO4 (unlike alkenes).",
            "Neutral; very weakly acidic C–H (pKa ≈ 43).",
            "Toxic and carcinogenic (long-term exposure is linked to leukemia).",
        ],
        "questions": [
            ("Which SMILES string represents benzene?",
             ["c1ccccc1", "C1CCCCC1", "CCO", "c1ccc2ccccc2c1"], 0),
            ("How many delocalised π electrons does benzene have?",
             ["4", "6", "8", "10"], 1),
            ("What is the approximate boiling point of benzene?",
             ["35 °C", "80 °C", "118 °C", "218 °C"], 1),
            ("Which reaction type is characteristic of benzene?",
             ["Electrophilic substitution", "Nucleophilic addition",
              "Free-radical polymerisation", "Elimination"], 0),
            ("What is the dipole moment of benzene?",
             ["0 D", "1.5 D", "2.3 D", "3.0 D"], 0),
        ],
    },
    "Naphthalene": {
        "smiles": "c1ccc2ccccc2c1",
        "iupac": "Naphthalene",
        "physical": {
            "Molecular formula": "C10H8",
            "Molar mass": "128.17 g/mol",
            "Appearance": "White crystalline solid (flakes or balls)",
            "Odor": "Strong, characteristic mothball odor",
            "Melting point": "80.3 °C",
            "Boiling point": "218 °C",
            "Density": "1.145 g/cm³ (solid, 20 °C)",
            "Solubility in water": "≈ 31 mg/L (25 °C); soluble in ether, benzene, chloroform",
            "Vapor pressure": "≈ 8.5 Pa (25 °C)",
            "Flash point": "80 °C",
            "log P (octanol/water)": "3.30",
            "Special property": "Sublimes readily at room temperature",
            "Structure": "Planar, two fused benzene rings (bicyclic, 10 π electrons)",
        },
        "chemical": [
            "Aromatic (10 π electrons) but less aromatic per ring than benzene, so it is more reactive than benzene.",
            "Undergoes electrophilic substitution preferentially at the α-position (C1), the kinetically favoured site.",
            "Nitration gives 1-nitronaphthalene; halogenation gives 1-halonaphthalenes.",
            "Sulfonation is temperature dependent: ~80 °C gives the 1-sulfonic acid (kinetic), ~160 °C gives the 2-sulfonic acid (thermodynamic).",
            "Oxidation with O2/V2O5 at high temperature gives phthalic anhydride (industrially important).",
            "Reduction: Na/ethanol gives 1,4-dihydronaphthalene; H2/Ni gives tetralin and then decalin.",
            "Forms a red-orange picrate (m.p. 149 °C) with picric acid, used for identification.",
            "Combustion gives CO2 and H2O with a sooty flame.",
            "Used as a moth repellent, and as a precursor of phthalic anhydride, dyes and naphthols.",
        ],
        "questions": [
            ("Which SMILES string represents naphthalene?",
             ["c1ccccc1", "c1ccc2ccccc2c1", "c1cncnc1", "CCO"], 1),
            ("What is the molecular formula of naphthalene?",
             ["C6H6", "C8H8", "C10H8", "C10H12"], 2),
            ("Which physical behaviour is characteristic of naphthalene?",
             ["It sublimes at room temperature", "It is a gas at room temperature",
              "It is miscible with water", "It is a liquid at room temperature"], 0),
            ("At which position does electrophilic substitution mainly occur (kinetic control)?",
             ["α (C1)", "β (C2)", "C9", "C10"], 0),
            ("Oxidation of naphthalene with O2/V2O5 gives:",
             ["Benzoic acid", "Phthalic anhydride", "Phenol", "Acetic acid"], 1),
        ],
    },
    "Salicylic acid": {
        "smiles": "OC(=O)c1ccccc1O",
        "iupac": "2-Hydroxybenzoic acid",
        "physical": {
            "Molecular formula": "C7H6O3",
            "Molar mass": "138.12 g/mol",
            "Appearance": "White crystalline powder / colorless needles",
            "Odor": "Odorless",
            "Melting point": "158.6 °C",
            "Boiling point": "≈ 211 °C at 20 mmHg (sublimes; decomposes at atmospheric pressure)",
            "Density": "1.443 g/cm³",
            "Solubility in water": "2.48 g/L (25 °C); soluble in ethanol, ether, acetone",
            "pKa": "2.97 (COOH); ≈ 13.6 (phenolic OH)",
            "log P (octanol/water)": "2.26",
            "Flash point": "157 °C",
            "Structure": "Benzene ring with –COOH and –OH ortho to each other (intramolecular H-bond)",
        },
        "chemical": [
            "Bifunctional: contains a carboxylic acid and a phenolic hydroxyl group.",
            "Stronger acid than benzoic acid (pKa 4.2) because the ortho –OH stabilises the carboxylate through intramolecular hydrogen bonding.",
            "Acetylation with acetic anhydride (H+ catalyst) gives acetylsalicylic acid (aspirin).",
            "Esterification with methanol/H2SO4 gives methyl salicylate (oil of wintergreen).",
            "Gives a violet/purple colour with neutral FeCl3 (test for phenolic –OH).",
            "Reacts with NaHCO3 to liberate CO2 (carboxylic acid test).",
            "Heating (or with soda lime) causes decarboxylation to phenol.",
            "Undergoes EAS (e.g. bromination, nitration) at positions ortho/para to the –OH group.",
            "Prepared industrially by the Kolbe–Schmitt reaction (sodium phenoxide + CO2, then H+).",
            "Keratolytic action: used in acne and wart treatments; precursor of aspirin.",
        ],
        "questions": [
            ("Which SMILES string represents salicylic acid?",
             ["CC(=O)Oc1ccccc1C(=O)O", "OC(=O)c1ccccc1O", "CCO", "c1ccccc1"], 1),
            ("Which functional groups does salicylic acid contain?",
             ["Carboxylic acid and phenol", "Ester and ether", "Aldehyde and amine",
              "Ketone and alcohol"], 0),
            ("Which colour does salicylic acid give with neutral FeCl3?",
             ["Green", "Red-brown", "Violet", "No colour"], 2),
            ("Acetylation of salicylic acid gives:",
             ["Methyl salicylate", "Aspirin", "Phenol", "Benzoic acid"], 1),
            ("Which is the approximate pKa of the carboxylic acid group?",
             ["2.97", "4.76", "9.9", "15.9"], 0),
        ],
    },
    "Ethanol": {
        "smiles": "CCO",
        "iupac": "Ethanol",
        "physical": {
            "Molecular formula": "C2H6O",
            "Molar mass": "46.07 g/mol",
            "Appearance": "Colorless, clear, volatile liquid",
            "Odor": "Characteristic, pleasant alcoholic odor",
            "Melting point": "−114.1 °C",
            "Boiling point": "78.37 °C",
            "Density": "0.789 g/cm³ (20 °C)",
            "Solubility in water": "Miscible in all proportions",
            "Vapor pressure": "5.95 kPa (20 °C)",
            "Refractive index (nD)": "1.3611",
            "Viscosity": "1.074 mPa·s (25 °C)",
            "Flash point": "13 °C",
            "Dipole moment": "1.69 D",
            "log P (octanol/water)": "−0.31",
            "Structure": "CH3–CH2–OH; hydrogen-bonding polar molecule",
        },
        "chemical": [
            "Very weakly acidic (pKa ≈ 15.9), about as acidic as water.",
            "Reacts with active metals (Na, K): 2 C2H5OH + 2 Na → 2 C2H5ONa + H2↑.",
            "Oxidation: with PCC gives acetaldehyde; with K2Cr2O7/H+ or KMnO4 gives acetic acid.",
            "Dehydration: conc. H2SO4 at ~170 °C gives ethene; at ~140 °C gives diethyl ether.",
            "Esterification (Fischer): with a carboxylic acid and H+ gives an ester, e.g. ethyl acetate.",
            "Reaction with HX, PCl5, PCl3 or SOCl2 gives ethyl halides.",
            "Iodoform test: I2/NaOH gives yellow CHI3 (positive for CH3CH(OH)– groups).",
            "Combustion: C2H5OH + 3 O2 → 2 CO2 + 3 H2O (clean blue flame, used as fuel).",
            "Fermentation of sugars by yeast gives ethanol (C6H12O6 → 2 C2H5OH + 2 CO2).",
            "Acts as a solvent and a mild nucleophile; denatured forms are used industrially.",
        ],
        "questions": [
            ("Which SMILES string represents ethanol?",
             ["CCO", "CC(=O)O", "C1CCCCC1", "COC"], 0),
            ("What is the approximate boiling point of ethanol?",
             ["34.6 °C", "78.4 °C", "100 °C", "118 °C"], 1),
            ("Dehydration of ethanol with conc. H2SO4 at ~170 °C gives:",
             ["Ethene", "Diethyl ether", "Acetaldehyde", "Ethane"], 0),
            ("Ethanol reacts with sodium metal to give which gas?",
             ["O2", "CO2", "H2", "Cl2"], 2),
            ("How does ethanol behave in water?",
             ["Immiscible", "Slightly soluble", "Miscible in all proportions",
              "Reacts violently"], 2),
        ],
    },
    "Pyrimidine": {
        "smiles": "c1cncnc1",
        "iupac": "Pyrimidine (1,3-diazine)",
        "physical": {
            "Molecular formula": "C4H4N2",
            "Molar mass": "80.09 g/mol",
            "Appearance": "Colorless liquid or low-melting crystalline solid",
            "Odor": "Pyridine-like",
            "Melting point": "20–22 °C",
            "Boiling point": "123–124 °C",
            "Density": "1.016 g/cm³",
            "Solubility in water": "Miscible (very soluble); soluble in most organic solvents",
            "pKa (conjugate acid)": "1.23",
            "Dipole moment": "≈ 2.3 D",
            "log P (octanol/water)": "−0.40",
            "Structure": "Planar six-membered aromatic ring with N at positions 1 and 3",
        },
        "chemical": [
            "Aromatic heterocycle with 6 π electrons; each nitrogen lone pair lies in the ring plane (not part of the π system).",
            "Much weaker base than pyridine (pKa 5.2) because the second nitrogen withdraws electron density; protonates on one N only.",
            "π-Deficient ring: resistant to electrophilic substitution, which occurs (if at all) at C5, the least electron-poor carbon.",
            "Nucleophilic substitution is easy at C2, C4 and C6 (e.g. halopyrimidines react with amines, alkoxides).",
            "N-alkylation (with alkyl halides) gives pyrimidinium salts; peracids give N-oxides.",
            "Can be reduced (e.g. with NaBH4 or catalytic hydrogenation) to di- and tetrahydropyrimidines.",
            "Ring is stable to oxidation, but alkyl side chains can be oxidised.",
            "Core of the nucleobases cytosine, thymine and uracil, and of vitamin B1 (thiamine), barbiturates and many drugs.",
            "Classic synthesis: Pinner synthesis (1,3-dicarbonyl + amidine) and Biginelli reaction.",
        ],
        "questions": [
            ("Which SMILES string represents pyrimidine?",
             ["c1cncnc1", "c1ccncc1", "c1ccccc1", "c1ccc2ccccc2c1"], 0),
            ("How many nitrogen atoms are in the pyrimidine ring?",
             ["1", "2", "3", "4"], 1),
            ("Nucleophilic substitution on pyrimidine occurs mainly at:",
             ["C2, C4, C6", "C5 only", "N1 only", "No position"], 0),
            ("Which of these is a pyrimidine-based nucleobase?",
             ["Adenine", "Guanine", "Cytosine", "Caffeine"], 2),
            ("What is the molecular formula of pyrimidine?",
             ["C5H5N", "C4H4N2", "C4H6N2", "C6H6"], 1),
        ],
    },
    "Aspirin": {
        "smiles": "CC(=O)Oc1ccccc1C(=O)O",
        "iupac": "2-(Acetyloxy)benzoic acid (acetylsalicylic acid)",
        "physical": {
            "Molecular formula": "C9H8O4",
            "Molar mass": "180.16 g/mol",
            "Appearance": "White crystalline powder or colorless crystals",
            "Odor": "Odorless (vinegar smell when hydrolysed)",
            "Melting point": "135 °C (decomposes)",
            "Boiling point": "≈ 140 °C (decomposes)",
            "Density": "1.40 g/cm³",
            "Solubility in water": "≈ 3 g/L (20 °C); soluble in ethanol, ether, chloroform",
            "pKa": "3.5 (25 °C)",
            "log P (octanol/water)": "1.19",
            "Flash point": "250 °C (approx.)",
            "Structure": "Benzene ring with ortho –COOH and –OCOCH3 (ester) groups",
        },
        "chemical": [
            "Contains a carboxylic acid and a phenolic ester (acetate) group.",
            "Hydrolysis (moisture, acid or base) gives salicylic acid + acetic acid, hence the vinegar smell of old tablets.",
            "Acidic (pKa ≈ 3.5); reacts with NaHCO3 to give CO2 and sodium acetylsalicylate (soluble salt).",
            "Gives no colour with FeCl3 when pure (the phenol is masked); a violet colour appears after hydrolysis.",
            "Synthesised by acetylation of salicylic acid with acetic anhydride (H2SO4 or H3PO4 catalyst).",
            "Acetylates a serine residue of cyclooxygenase (COX-1/COX-2), irreversibly inhibiting prostaglandin synthesis.",
            "Pharmacology: analgesic, antipyretic, anti-inflammatory and antiplatelet (low-dose cardioprotection).",
            "Decomposes on heating above its melting point; incompatible with strong bases and moisture.",
        ],
        "questions": [
            ("Which SMILES string represents aspirin?",
             ["OC(=O)c1ccccc1O", "CC(=O)Oc1ccccc1C(=O)O", "CCO", "c1cncnc1"], 1),
            ("Hydrolysis of aspirin gives:",
             ["Salicylic acid and acetic acid", "Phenol and CO2",
              "Benzoic acid and methanol", "Ethanol and benzene"], 0),
            ("Aspirin irreversibly inhibits which enzyme?",
             ["Cyclooxygenase (COX)", "Catalase", "Amylase", "DNA polymerase"], 0),
            ("What is the molar mass of aspirin?",
             ["138.12 g/mol", "180.16 g/mol", "46.07 g/mol", "128.17 g/mol"], 1),
            ("Which functional groups are present in aspirin?",
             ["Ester and carboxylic acid", "Aldehyde and ketone", "Amine and amide",
              "Ether and alcohol"], 0),
        ],
    },
}


# ----------------------------------------------------------------------------
# HELPERS
# ----------------------------------------------------------------------------
def draw_2d(smiles: str, size=(450, 350)):
    mol = Chem.MolFromSmiles(smiles)
    return Draw.MolToImage(mol, size=size)


def viewer_3d_html(smiles: str, spin: bool, width=600, height=420) -> str:
    """Build a 3D structure (RDKit) and render it with 3Dmol.js."""
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(mol, randomSeed=42)
    try:
        AllChem.MMFFOptimizeMolecule(mol)
    except Exception:
        pass
    molblock = Chem.MolToMolBlock(mol)
    return f"""
    <div id="viewer" style="width:{width}px; max-width:100%; height:{height}px;
         position:relative; border:1px solid #ccc; border-radius:8px;"></div>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.0.4/3Dmol-min.js"></script>
    <script>
      const viewer = $3Dmol.createViewer("viewer", {{backgroundColor: "white"}});
      viewer.addModel({json.dumps(molblock)}, "mol");
      viewer.setStyle({{}}, {{stick: {{radius: 0.15}}, sphere: {{scale: 0.28}}}});
      viewer.zoomTo();
      viewer.render();
      {"viewer.spin(true);" if spin else ""}
    </script>
    """


# ----------------------------------------------------------------------------
# APP
# ----------------------------------------------------------------------------
st.set_page_config(page_title="Compound Dashboard", page_icon="🧪", layout="wide")
st.title("🧪 Chemical Compound Dashboard")

tabs = st.tabs([
    "1. Compound",
    "2. 3D Structure",
    "3. Physical Properties",
    "4. Chemical Properties",
    "5. SMILES",
    "6. Assessment",
])

# ---- 1. Compound -----------------------------------------------------------
with tabs[0]:
    st.header("1. Compound")
    st.selectbox("Select a compound", list(COMPOUNDS.keys()), key="compound")

name = st.session_state.get("compound", list(COMPOUNDS.keys())[0])
data = COMPOUNDS[name]

with tabs[0]:
    st.subheader(name)
    st.write(f"**IUPAC name:** {data['iupac']}")
    st.write(f"**Molecular formula:** {data['physical']['Molecular formula']}")
    st.write(f"**Molar mass:** {data['physical']['Molar mass']}")
    st.image(draw_2d(data["smiles"]), caption=f"2D structure of {name}")

# ---- 2. 3D structure -------------------------------------------------------
with tabs[1]:
    st.header(f"2. 3D Structure of {name}")
    spin = st.checkbox("Auto-rotate", value=True, key="spin")
    components.html(viewer_3d_html(data["smiles"], spin), height=450)
    st.caption("Drag to rotate, scroll to zoom. Ball-and-stick model (hydrogens included).")

# ---- 3. Physical properties ------------------------------------------------
with tabs[2]:
    st.header(f"3. Physical Properties of {name}")
    col1, col2 = st.columns([1, 1.4])
    with col1:
        st.subheader("2D structure")
        st.image(draw_2d(data["smiles"]), caption=name)
    with col2:
        st.subheader("Properties")
        st.table({"Property": list(data["physical"].keys()),
                  "Value": list(data["physical"].values())})

# ---- 4. Chemical properties ------------------------------------------------
with tabs[3]:
    st.header(f"4. Chemical Properties of {name}")
    for item in data["chemical"]:
        st.markdown(f"- {item}")

# ---- 5. SMILES -------------------------------------------------------------
with tabs[4]:
    st.header(f"5. SMILES Representation of {name}")
    st.code(data["smiles"], language="text")
    mol = Chem.MolFromSmiles(data["smiles"])
    st.write(f"**Canonical SMILES (RDKit):** `{Chem.MolToSmiles(mol)}`")
    st.write(f"**Atoms (heavy):** {mol.GetNumAtoms()}   |   **Bonds:** {mol.GetNumBonds()}")
    st.markdown(
        "**How to read it:** lowercase letters (`c`, `n`) denote aromatic atoms, "
        "digits mark ring-closure points, `=` is a double bond, parentheses show branches, "
        "and `O` / `C` are oxygen / aliphatic carbon."
    )
    st.image(draw_2d(data["smiles"]), caption=name)

# ---- 6. Assessment ---------------------------------------------------------
with tabs[5]:
    st.header(f"6. Assessment: {name}")
    st.write("Answer all 5 questions, then press **Submit**.")
    with st.form(key=f"quiz_{name}"):
        answers = []
        for i, (q, options, _) in enumerate(data["questions"]):
            answers.append(
                st.radio(f"Q{i + 1}. {q}", options, index=None, key=f"{name}_q{i}")
            )
        submitted = st.form_submit_button("Submit")

    if submitted:
        score = 0
        for i, (q, options, correct) in enumerate(data["questions"]):
            if answers[i] == options[correct]:
                score += 1
                st.success(f"Q{i + 1}: Correct ✅ ({options[correct]})")
            else:
                st.error(f"Q{i + 1}: Incorrect ❌. Correct answer: {options[correct]}")
        st.subheader(f"Score: {score} / {len(data['questions'])}")
        if score == len(data["questions"]):
            st.balloons()
