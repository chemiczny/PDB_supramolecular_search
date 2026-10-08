"""
Created on Fri Aug  3 10:37:23 2018

@author: michal
"""

import pandas as pd

CATIONIC_AA = ["LYS", "ARG"]
ACIDIC_AA = ["ASP", "GLU"]
REST_AA = [
    "ALA",
    "CYS",
    "GLY",
    "ILE",
    "LEU",
    "MET",
    "ASN",
    "PRO",
    "GLN",
    "SER",
    "THR",
    "VAL",
]
ALL_AA = REST_AA + ACIDIC_AA
RING_AA = ["PHE", "HIS", "TRP", "TYR"]
NU = ["A", "G", "T", "C", "U", "I", "DA", "DC", "DG", "DT", "DI"]


def no_aa_in_pi_acids(actual_data):
    return actual_data[~actual_data["Pi acid Code"].isin(RING_AA)]


def no_aa_in_pi_res(actual_data):
    return actual_data[~actual_data["Pi res code"].isin(RING_AA)]


def no_aa_in_anions(actual_data):
    return actual_data[~actual_data["Anion code"].isin(ALL_AA)]


def no_nu_in_anions(actual_data):
    return actual_data[~actual_data["Anion code"].isin(NU)]


def no_nu_in_pi_acids(actual_data):
    return actual_data[~actual_data["Pi acid Code"].isin(NU)]


def no_nu_in_pi_res(actual_data):
    return actual_data[~actual_data["Pi res code"].isin(NU)]


def no_aa_in_h_acceptors(actual_data):
    return actual_data[~actual_data["Acceptor code"].isin(ACIDIC_AA)]


def no_nu_in_h_acceptors(actual_data):
    return actual_data[~actual_data["Acceptor code"].isin(NU)]


def no_aa_in_h_donors(actual_data):
    return actual_data[
        ~actual_data["Donor code"].isin(["TYR", "PHE", "HIS", "TRP", "LYS", "ARG"])
    ]


def no_nu_in_h_donors(actual_data):
    return actual_data[~actual_data["Donor code"].isin(NU)]


def no_aa_in_cations(actual_data):
    return actual_data[~actual_data["Cation code"].isin(CATIONIC_AA)]


def only_anions(actual_data):
    return actual_data[actual_data["isAnion"].eq(True)]


def only_complexes(actual_data):
    return actual_data[actual_data["Complex"].eq(True)]


def simple_merge(
    data_frames_to_merge,
    data_frame_merge_headers,
    data_frames_to_exclude,
    data_frame_exclude_headers,
):
    if len(data_frames_to_merge) + len(data_frames_to_exclude) < 2:
        return

    unique_data = []
    data_excluded = []

    actual_keys = []
    excluded_keys = []

    for df, headers in zip(data_frames_to_merge, data_frame_merge_headers):
        if len(unique_data) == 0:
            unique_data = df[headers].drop_duplicates()
        elif len(df) > 0:
            unique_data = pd.merge(
                unique_data, df[headers], on=list(set(actual_keys) & set(headers))
            )
            unique_data = unique_data.drop_duplicates()
        actual_keys = list(set(actual_keys + headers))

    for df, headers in zip(data_frames_to_exclude, data_frame_exclude_headers):
        if len(data_excluded) == 0:
            data_excluded = df[headers].drop_duplicates()
        elif len(df) > 0:
            data_excluded = pd.merge(
                data_excluded, df[headers], on=list(set(actual_keys) & set(headers))
            )
            data_excluded = data_excluded.drop_duplicates()
        excluded_keys = list(set(excluded_keys + headers))

    if len(data_excluded) > 0:
        merging_keys = list(set(actual_keys) & set(excluded_keys))
        sub_merged = pd.merge(
            unique_data, data_excluded, on=merging_keys, how="left", indicator=True
        )
        unique_data = sub_merged[sub_merged["_merge"] == "left_only"]

    all_df = data_frames_to_merge + data_frames_to_exclude
    all_headers = data_frame_merge_headers + data_frame_exclude_headers

    new_df = []
    for df, headers in zip(all_df, all_headers):
        if len(unique_data) == 0:
            break

        merging_keys = list(set(actual_keys) & set(headers))
        temp_data_frame = unique_data[merging_keys].drop_duplicates()
        new_df.append(pd.merge(df, temp_data_frame, on=merging_keys))

    return new_df
