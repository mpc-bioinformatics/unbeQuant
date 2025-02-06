#!/bin/env python


import argparse
import numpy as np
import pyopenms
import os 
import pandas as pd
import json
import tqdm
from ast import literal_eval


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--feature_tsv", help="The Feature-TSV file from biosaur2")
    parser.add_argument("--hills_tsv", help="The Hills-TSV file from biosaur2 (retriavable via 'write_hills')")
    parser.add_argument("--mzml", help="The original mzML file to create a connection of MS2-events to the features")
    parser.add_argument("--mz_tolerance", help="The PPM-Error Tolerance, as has been set in biosaur2", type=int)
    parser.add_argument("--rt_enlarge_factor", help="The Retention Time enlargement factor. The higher, the further away MS2 spectra (in RT dimension) can match to a (or more) features", type=float, default=0.0)
    parser.add_argument("--output_featurexml", help="The output featureXML file, inwhich it stores the information. WARNING: This file (if it exists) will be overwritten!")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    # Read all features, but only necessary columns via pandas
    features_df = pd.read_csv(args.feature_tsv, sep="\t", usecols=["charge", "nIsotopes", "intensitySum", "rtApex", "nScans", "mz", "rtStart", "rtEnd", "isoerror", "isoerror2", "isoerror", "isotopes", "monoisotope hill idx"])

    # Read all hills, but only neccessary columns via pandas
    hills_df = pd.read_csv(args.hills_tsv, sep="\t", index_col="hill_idx", usecols=["rtStart", "rtEnd", "mz", "hill_idx", "hills_scan_lists", "hills_mz_array"])

    # Load each feature and isotope (hill) into custom dict:
    print("Loading all features from biosaur2")
    all_features = []
    for feature_idx, feature_row in tqdm.tqdm(features_df.iterrows(), total=len(features_df)):
        feat_d = dict()
        feat_d["charge"] = feature_row["charge"]
        feat_d["intensitySum"] = feature_row["intensitySum"]
        feat_d["rtApex"] = feature_row["rtApex"]
        feat_d["nIsotopes"] = feature_row["nIsotopes"]
        feat_d["isoerror"] = feature_row["isoerror"]
        feat_d["isoerror2"] = feature_row["isoerror2"]

        isos = []
        for iso in [dict(isotope_hill_idx=feature_row["monoisotope hill idx"])] + json.loads(feature_row["isotopes"].replace("'", '"')):
            hill_d = dict()
            hill = hills_df.loc[iso["isotope_hill_idx"]]
            hill_d["rtStart"] = hill["rtStart"]
            hill_d["rtEnd"] = hill["rtEnd"]
            hill_d["mz"] = hill["mz"]
            hill_d["hills_mz_array"] = hill["hills_mz_array"]
            isos.append(hill_d)

        feat_d["isotopes"] = isos
        all_features.append(feat_d)

    ### Read mzML file  (but only we only retrieve the MS2 scans, specifically the m/z, charge, RT and scan number)
    ### and save it in the dataframe
    print("Loading MS2 scans")
    exp = pyopenms.MSExperiment()
    pyopenms.MzMLFile().load(args.mzml, exp)
    ms2_rts = []
    ms2_mzs = []
    ms2_cs = []
    ms2_scans = []
    for spec in tqdm.tqdm(exp.getSpectra()):

        if spec.getMSLevel() == 2:
            p = spec.getPrecursors()[0]
            ms2_cs.append(p.getCharge())
            ms2_mzs.append(p.getMZ())
            ms2_scans.append(int(spec.getNativeID()[spec.getNativeID().index("scan")+5:].split(" ", 1)[0]))  # Get Scan "scan=XXXX" from the native ID information
            ms2_rts.append(spec.getRT())  # In Seconds!!!

    ms2_df = pd.DataFrame.from_dict(dict(rt=ms2_rts, mz=ms2_mzs, charge=ms2_cs, scan=ms2_scans))
    ms2_df["mz_tol"] = (ms2_df["mz"]/1000000)*args.mz_tolerance



    # Write features in featureXML
    print("Generating PyOpenMS Features")
    features = pyopenms.FeatureMap()
    for f_d in tqdm.tqdm(all_features):
        f = pyopenms.Feature()
        f.setCharge(f_d["charge"])
        f.setIntensity(f_d["intensitySum"])
        f.setOverallQuality(0.9)

        # Save infos from the very first isotope into the feature
        iso0 = f_d["isotopes"][0]
        f.setMZ(np.median(literal_eval(iso0["hills_mz_array"])))
        f.setRT(f_d["rtApex"]*60)

        # Create convex hull for each isotope
        chulls = []
        for iso in f_d["isotopes"]:
            hull = pyopenms.ConvexHull2D()

            iso_mz = np.median(literal_eval(iso["hills_mz_array"]))
            iso_mz_tol = (iso_mz/1000000)*5
            upper_mz = iso_mz + iso_mz_tol
            lower_mz = iso_mz - iso_mz_tol

            hull.addPointXY(iso["rtStart"]*60, lower_mz)
            hull.addPointXY(iso["rtEnd"]*60, lower_mz)
            hull.addPointXY(iso["rtEnd"]*60, upper_mz)
            hull.addPointXY(iso["rtStart"]*60, upper_mz)

            chulls.append(hull)

        # Set isotopes and id, add it to the list of features
        f.setConvexHulls(chulls)
        f.ensureUniqueId()
        features.push_back(f)


    ### Post processing, afterwards we match alls the ms2 spectra with the features and add it as a CV param
    print("Matching MS2 to Features")
    ms2_matching_to_feature = []
    for f in tqdm.tqdm(features, total=features.size()):
        # Match each feature if it corresponds to a ms2
        poss_ms2 = ms2_df[ms2_df["charge"] == f.getCharge()]
        min_box = [x.getBoundingBox().minPosition() for x in f.getConvexHulls()]
        max_box = [x.getBoundingBox().maxPosition() for x in f.getConvexHulls()]
        min_rt = min([x[0] for x in min_box])
        max_rt = max([x[0] for x in max_box])
        rt_tolerance = (max_rt-min_rt)*args.rt_enlarge_factor

        ms2_matching = []
        if len(poss_ms2) > 0:
            overlapping_rt = np.logical_and(poss_ms2["rt"] >= (min_rt - rt_tolerance), poss_ms2["rt"] <= (max_rt + rt_tolerance))
            remianing_ms = poss_ms2[overlapping_rt]

            if len(remianing_ms) > 0:
                for (_, min_mz), (_, max_mz) in zip(min_box, max_box):
                    # Check mz for overlapping intervals
                    # poss_ms2 = a
                    # min_mz, max_mz = b
                    overlapping_mz = np.logical_and(   # b2>=a1   and  a2>=b1
                        max_mz >= (remianing_ms["mz"] - remianing_ms["mz_tol"]),
                        (remianing_ms["mz"] + remianing_ms["mz_tol"]) >= min_mz
                    )

                    remaining_after_rt_and_mz_filter = remianing_ms[overlapping_mz]

                    if len(remaining_after_rt_and_mz_filter) > 0:
                        ms2_matching += [x for x in remaining_after_rt_and_mz_filter["scan"]]

        ms2_matching_to_feature.append(ms2_matching)


    ### Match Features with MS2s and write the final results featureXML
    final_features = pyopenms.FeatureMap()
    final_features.ensureUniqueId()  # Needed for the ConsensusMap generation via OpenMs
    final_features.setMetaValue("spectra_data", [os.path.basename(args.mzml)])  # Needed for the ConsensusMap generation
    for ms2s, f  in zip(ms2_matching_to_feature, features):
        f.setMetaValue("unbeQuant_MS2_Scan_Map", ms2s)
        final_features.push_back(f)

    fh = pyopenms.FeatureXMLFile()
    fh.store(args.output_featurexml, final_features)
