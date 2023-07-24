import streamlit as st
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.Seq import Seq
from streamlit_tags import st_tags
import pandas as pd
import json


# Declare chemical functions on top
def calc_mol_weight_prot_param(seq):
    pp = ProteinAnalysis(str(seq))
    return pp.molecular_weight()


# Logic for the prediciton of cleavage sites
def predict_cleavage_sites(seq_list, amino_sequence, default_mass, margin):
    # Difference is the default mass
    if margin == "" or margin is None:
        error_rate = 5
    else:
        error_rate = margin

    margin_of_error = error_rate / 100 * default_mass
    saved_dict = {}
    prob_dict = {}
    for ik in range(0, len(seq_list)):
        if ik == 0:
            prob_dict[amino_sequence[0]] = calc_mol_weight_prot_param(amino_sequence[0])
        else:
            prob_dict[amino_sequence[0:ik]] = calc_mol_weight_prot_param(amino_sequence[0:ik])

    for keys, value in prob_dict.items():
        diff = abs(default_mass - value)
        if diff <= margin_of_error:
            saved_dict[keys] = [value, diff]

    saved_dict = dict(sorted(saved_dict.items(), key=lambda item: item[1][1]))
    return saved_dict


# Mathematical formula for prediction cleavage sites based on molecular weight difference
def calc_mass_diff(diff, disul_bonds, nq, ck, water_value, n_sugar=None, o_sugar=None, mod=None):
    mass_diff = diff + 2 * int(disul_bonds) - nq + ck - water_value - float(n_sugar if n_sugar else 0) - \
                float(o_sugar if o_sugar else 0) - float(mod if mod else 0)
    return mass_diff


# Get values of PTMs stored in a dictionary
def get_values_from_dict(d, val):
    if val is None:
        return None
    else:
        try:
            position = list(d.values()).index(float(val))
            return list(d.keys())[position]
        except ValueError:
            return None


# Creating rows on Streamlit GUI
def create_row(inputs, diff, mass_diff, pot_cleav_sites, keys):
    new_row = {"Inputs": inputs, "Initial Mass Difference": abs(diff), "Modified Mass Difference": abs(mass_diff),
               "Cleavage Sites": keys}
    return new_row


# Declare GUI Functions
def highlight_sequence(seq, target_aas):
    target_aas_count_dict = {}
    colors = ["yellow", "cyan", "navajowhite", "lime", "orange", "pink"]
    result_str = ""
    with open('mod.json') as s:
        high_dict = json.load(s)['Highlights']

    for aa in seq:
        if aa in target_aas:
            color_index = target_aas.index(aa)
            result_str += f'<span style="background-color: {colors[color_index]};">{aa}</span>'
            if aa not in target_aas_count_dict:
                target_aas_count_dict[aa] = 1
            else:
                target_aas_count_dict[aa] += 1
        else:
            result_str += aa

    # Display the count of the number of amino acids that the user wants to find
    st.markdown(result_str, unsafe_allow_html=True)
    for target_aa, value in sorted(target_aas_count_dict.items()):
        new_string = f"Count of amino acid {target_aa} {high_dict[target_aa]}: " + str(value)
        st.markdown(new_string, unsafe_allow_html=True)

    # Display the count of the entire length of user's input amino acid sequence
    total_length = len(seq)
    st.markdown(f'Total length of amino acid sequence: {total_length}', unsafe_allow_html=True)


# Constant
water = 18
st.session_state['error_check'] = 0
st.session_state['error_message'] = ''


# Error Checking

# Error Checking for Highlighting Sequence
def check_input_highlight(a_seq, a_st):
    st.session_state['error_message'] = ''
    # Added a new condition such that if the string is empty, it'll also trigger this
    if a_seq.isdigit() or a_seq == "":
        st.session_state['error_check'] = 1
        st.session_state['error_message'] += 'Full amino acid sequence should only be alphabets and cannot be empty.'
    if all(item.isdigit() for item in a_st):
        st.session_state['error_check'] = 1
        # Added a check because due to separators of \n problems
        if st.session_state["error_message"]:
            st.session_state['error_message'] += '\nAmino acid sequence search filter should only be alphabets.'
        else:
            st.session_state['error_message'] += 'Amino Acid Sequence search filter should only be alphabets.'


# Define GUI Script
st.set_page_config(layout="wide")
with open("styles/styles.css") as source_des:
    st.markdown(f"<style>{source_des.read()} </style>", unsafe_allow_html=True)

st.title("Amino Acid Sequence Analyzer for QC Department")
st.markdown('<hr class=hr-1></hr>', unsafe_allow_html=True)

tab1, tab2 = st.tabs(["Amino Acid Finder & Highlighter", "Cleavage Site Predictor"])

# Tab 1 - amino acid highlighting function
with tab1:
    amino_seq = st.text_input(label="Enter target amino acid sequence")
    amino_short = st_tags(label="Enter the amino acid(s) to search (e.g, 'C', 'G', 'N')", maxtags=6)
    analyse = st.button(label="Analyze Sequence")

    st.markdown('<hr class=hr-1></hr>', unsafe_allow_html=True)

    if analyse:
        check_input_highlight(amino_seq, amino_short)
        col1, _, col2 = st.columns([1, 0.5, 1])
        if st.session_state['error_check'] == 1:
            with col1:
                st.subheader("Output Results")
                error_list = st.session_state['error_message'].split("\n")
                for i in error_list:
                    st.error(i)
            with col2:
                st.subheader("Molecular Weight")
                st.markdown("NIL", unsafe_allow_html=True)
        else:
            analyzed_seq = ProteinAnalysis(amino_seq)
            molecular_weight_mw = analyzed_seq.molecular_weight()

            with col1:
                st.subheader("Output Results")
                highlight_sequence(amino_seq, amino_short)

            with col2:
                st.subheader("Molecular Weight")
                molecular_weight = "{:.2f}".format(molecular_weight_mw) + " Dalton (Da)"
                st.markdown(molecular_weight, unsafe_allow_html=True)

# Tab 2 - cleavage site predictor (work in progress)
with tab2:
    with st.expander("Minimise"):
        amino_seq_input = st.text_input(label="Amino Acid Sequence")
        col3, col4, col5 = st.columns([1, 1, 1])
        col6, col7, col8 = st.columns([1, 1, 1])
        col9, col10, _ = st.columns([1, 1, 1])

        with col3:
            theorectical_mass = st.text_input(label="Theoretical Mass")
        with col4:
            observed_mass = st.text_input(label="Observed Mass")
        with col5:
            num_disul_bonds = st.text_input(label="Number of Disulfide Bonds")

        with col6:
            n_sugar_mass = st.text_input(label="N Sugar Mass (If present, seperate it by comma) e.g 1000, 2000")
        with col7:
            o_sugar_mass = st.text_input(label="O Sugar Mass (If present, seperate it by comma) e.g 1000, 2000")
        with col8:
            mod_mass = st.text_input(label="Modified Mass (If present, seperate it by comma) e.g 1000, 2000")

        with col9:
            error_margin = st.text_input(
                label="Error of Margin For Cleavage Site Prediction e.g (5 for 5% out of 100%)")
        with col10:
            terminus_selection = st.radio(label="Terminus Selection", options=("N-Terminus", "C-Terminus"), index=0)

        predict = st.button(label="Predict Cleavage Sites")

    st.markdown('<hr class=hr-1></hr>', unsafe_allow_html=True)

    if predict:
        sequence = Seq(amino_seq_input)

        # Opens the mod.json files
        with open('mod.json') as f:
            all_dict = json.load(f)

        n_sugar_dict = all_dict["N_Sugars"]
        o_sugar_dict = all_dict["O_Sugars"]
        mod_dict = all_dict["Other_Mods"]

        # This gets the input from the user, and obtains all inputs in a list based on a filter.
        n_sugar_list = list(filter(None, n_sugar_mass.split(", ")))
        o_sugar_list = list(filter(None, o_sugar_mass.split(", ")))
        mod_list = list(filter(None, mod_mass.split(", ")))

        # Switched to this location

        # Normal
        if terminus_selection == "N-Terminus":
            amino_seq_list = list(amino_seq_input)

        # Reverse
        elif terminus_selection == "C-Terminus":
            amino_seq_list = list(amino_seq_input)
            amino_seq_list.reverse()
            amino_seq_input = amino_seq_input[::-1]

        n_terminus_q = 17 if amino_seq_input[0] == "Q" else 0
        c_terminus_k = 128 if amino_seq_input[len(amino_seq_input) - 1] == "K" else 0

        difference = float(observed_mass) - float(theorectical_mass)

        df = pd.DataFrame(columns=["Inputs", "Initial Mass Difference", "Modified Mass Difference", "Cleavage Sites"])

        if difference < 0:
            for i in (n_sugar_list if len(n_sugar_list) > 0 else [None]):
                for j in (o_sugar_list if len(o_sugar_list) > 0 else [None]):
                    for k in (mod_list if len(mod_list) > 0 else [None]):
                        mass_diff_with_mod = calc_mass_diff(difference, num_disul_bonds, n_terminus_q, c_terminus_k,
                                                            water, i, j, k)

                        n_sugar_type = get_values_from_dict(n_sugar_dict, i)
                        o_sugar_type = get_values_from_dict(o_sugar_dict, j)
                        mod_type = get_values_from_dict(mod_dict, k)

                        n_sugar_input = f"N Sugar: {n_sugar_type} || " if n_sugar_type else ""
                        o_sugar_input = f"O Sugar: {o_sugar_type} || " if o_sugar_type else ""
                        mod_input = f"Modifications: {mod_type} || " if mod_type else ""
                        di_input = f"Disulfide Bonds: {str(num_disul_bonds)}"

                        potential_cleavage_sites = predict_cleavage_sites(amino_seq_list, amino_seq_input,
                                                                          abs(mass_diff_with_mod), float(error_margin))
                        all_input = n_sugar_input + o_sugar_input + mod_input + di_input

                        if not potential_cleavage_sites:
                            key = "No Potential Cleavages Found"
                        else:
                            key = list(potential_cleavage_sites.keys())[0]

                        row = create_row(all_input, difference, mass_diff_with_mod, potential_cleavage_sites, key)
                        df = pd.concat([df, pd.DataFrame([row])], ignore_index=True)

            st.dataframe(df, use_container_width=True)
