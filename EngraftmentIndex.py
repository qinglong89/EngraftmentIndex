#!/usr/bin/env python3

# Written by Qinglong Wu, PhD (Qinglong.Wu@bcm.edu or qinglong@connect.hku.hk)

import os
import pandas as pd
import math
from math import atan2, cos, radians, degrees, hypot

#specify the working directory
os.chdir('C:/Users/qinglong/PycharmProjects/EngraftmentIndex')

print(os.getcwd())

#import data
pcoa = pd.read_csv("pcoa_abund-jaccard_4EngraftmentIndexCalculation.csv")

#extract patient IDs and its corresponding donor info.
patient_donor = pcoa[["PatientID", "Donor"]][pcoa["PatientID"] != "Donor"].drop_duplicates(subset=['PatientID'])

print(patient_donor)

output = []

for row in patient_donor.index:
    #print(patient_donor["PatientID"][row], patient_donor["Donor"][row])
    pcoa_subset_patient = pcoa[pcoa["PatientID"] == patient_donor["PatientID"][row]]
    pcoa_subset_donor = pcoa[(pcoa["PatientID"] == "Donor") & (pcoa["Donor"] == patient_donor["Donor"][row])]

    pcoa_work = pd.concat([pcoa_subset_donor, pcoa_subset_patient], ignore_index = True, sort = False)

    #append a column with degree based on coordinate (x, y); use lambda for a Series
    pcoa_work["angle"] = pcoa_work[["PCo1", "PCo2"]].apply(lambda x: (degrees(atan2(x.PCo1, x.PCo2))+360)%360, axis=1)

    #append a column with relative degree to Donor
    pcoa_work["angle_relative_to_Donor"] = pcoa_work[["angle"]].apply(lambda x: abs(x.angle - pcoa_work[pcoa_work["Group"] == "Donor"].angle), axis=1)

    #append a column with cosine(rel_degree); use radians() to convert degrees into radians
    pcoa_work["cosine_angle_relative_to_Donor"] = pcoa_work[["angle_relative_to_Donor"]].apply(lambda x: cos(radians(x.angle_relative_to_Donor)), axis=1)

    #append a column with distance relative to PreFMT
    pcoa_work["distance_relative_to_PreFMT"] = pcoa_work[["PCo1", "PCo2"]].apply(lambda x: hypot(x.PCo1 - pcoa_work[pcoa_work["Group"] == "PreFMT"].PCo1, x.PCo2 - pcoa_work[pcoa_work["Group"] == "PreFMT"].PCo2), axis=1)

    #append a column with the ration of relative distance between PostFMTs and Donor
    pcoa_work["distance_ratio_PostFMT_vs_Donor"] = pcoa_work[["distance_relative_to_PreFMT"]].apply(lambda x: x.distance_relative_to_PreFMT / pcoa_work[pcoa_work["Group"] == "Donor"].distance_relative_to_PreFMT, axis=1)

    #append a column with engraftment index
    pcoa_work["relative_engraftment_index"] = pcoa_work[["distance_ratio_PostFMT_vs_Donor", "cosine_angle_relative_to_Donor"]].apply(lambda x: x.distance_ratio_PostFMT_vs_Donor * x.cosine_angle_relative_to_Donor, axis=1)

    pcoa_work["absolute_engraftment_index"] = pcoa_work[["distance_ratio_PostFMT_vs_Donor", "cosine_angle_relative_to_Donor", "distance_relative_to_PreFMT"]].apply(lambda x: x.distance_ratio_PostFMT_vs_Donor * x.cosine_angle_relative_to_Donor * x.distance_relative_to_PreFMT, axis=1)

    pcoa_PostFMTs = pcoa_work[pcoa_work["Group"] == "PostFMT"]

    output.append(pcoa_PostFMTs)

output_df = pd.concat(output, ignore_index = True, sort = False)

output_df.to_csv("EngraftmentIndex_calculated.csv")

