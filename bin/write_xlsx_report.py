#!/bin/env python

import xlsxwriter
import argparse
import csv
import tqdm

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", "-o", help="Output xlsx file")
    parser.add_argument("--input", "-i", help="Input tsv file")
    return parser.parse_args()


def merge_cells(y1, y2, header, workbook, worksheet):
    merge_format = workbook.add_format({"bold": True, "align": "center", "right": 2})
    if y1 == y2:
        worksheet.write(0, y1, header, merge_format)
    else:
        worksheet.merge_range(0, y1, 0, y2, header, merge_format)



if __name__ == "__main__":
    args = parse_args()

    # Get number of lines
    with open(args.input, "r") as in_file:
        no_lines = sum(1 for _ in in_file)

    with open(args.input, "r") as in_file:
        csv_reader = csv.reader(in_file, delimiter="\t")

        RAW_header = next(csv_reader)

        # Parse Logic, get the indices of the columns we are interested in. Yield None or empty list if the column is not found
        ceid_idx = RAW_header.index("openms_ceid")
        intensity_idcs = [i for i, h in enumerate(RAW_header) if h.endswith("_____sum_intensity")]
        pep_ident_idcs = [i for i, h in enumerate(RAW_header) if h.endswith("_____l_pep_ident")]
        prot_ident_idcs = [i for i, h in enumerate(RAW_header) if h.endswith("_____l_prot_ident")]
        charge_idcs = [i for i, h in enumerate(RAW_header) if h.endswith("_____charge")]
        ms2_scans_idcs = [i for i, h in enumerate(RAW_header) if h.endswith("_____l_ms2_scans")]
        global_min_mz_idx = RAW_header.index("feature_global_min_mz") if "feature_global_min_mz" in RAW_header else None
        global_max_mz_idx = RAW_header.index("feature_global_max_mz") if "feature_global_max_mz" in RAW_header else None
        global_min_rt_idx = RAW_header.index("feature_global_min_rt") if "feature_global_min_rt" in RAW_header else None
        global_max_rt_idx = RAW_header.index("feature_global_max_rt") if "feature_global_max_rt" in RAW_header else None
        first_iso_min_mz_idx = RAW_header.index("first_iso_global_min_mz") if "first_iso_global_min_mz" in RAW_header else None
        first_iso_max_mz_idx = RAW_header.index("first_iso_global_max_mz") if "first_iso_global_max_mz" in RAW_header else None
        first_iso_min_rt_idx = RAW_header.index("first_iso_global_min_rt") if "first_iso_global_min_rt" in RAW_header else None
        first_iso_max_rt_idx = RAW_header.index("first_iso_global_max_rt") if "first_iso_global_max_rt" in RAW_header else None
        fid_idcs = [i for i, h in enumerate(RAW_header) if h.endswith("_____openms_fid")]
        single_iso_bound_mz_start_idcs = [i for i, h in enumerate(RAW_header) if h.endswith("_____l_mz_start")]
        single_iso_bound_mz_end_idcs = [i for i, h in enumerate(RAW_header) if h.endswith("_____l_mz_end")]
        single_iso_bound_rt_start_idcs = [i for i, h in enumerate(RAW_header) if h.endswith("_____l_rt_start")]
        single_iso_bound_rt_end_idcs = [i for i, h in enumerate(RAW_header) if h.endswith("_____l_rt_end")]
        xic_intenst_idcs = [i for i, h in enumerate(RAW_header) if h.endswith("_____l_intensities")]
        xic_rt_idcs = [i for i, h in enumerate(RAW_header) if h.endswith("_____l_retention_times")]


        final_output_list = [
            (ceid_idx, "Identifier", 128, False),
            (intensity_idcs, "Intensities", 96, True),
            (pep_ident_idcs, "Peptide identification", 96, True),
            (prot_ident_idcs, "Protein identification", 96, True),
            (charge_idcs, "Charge", 96, True),
            (ms2_scans_idcs, "MS2_Scans", 96, True),
            (global_min_mz_idx, "Feature_Global_mz_start", None, False),
            (global_max_mz_idx, "Feature_Global_mz_end", None, False),
            (global_min_rt_idx, "Feature_Global_rt_start", None, False),
            (global_max_rt_idx, "Feature_Global_rt_end", None, False),
            (first_iso_min_mz_idx, "Feature_First_Isotope_mz_start", None, False),
            (first_iso_max_mz_idx, "Feature_First_Isotope_mz_end", None, False),
            (first_iso_min_rt_idx, "Feature_First_Isotope_rt_start", None, False),
            (first_iso_max_rt_idx, "Feature_First_Isotope_rt_end", None, False),
            (fid_idcs, "Feauture_IDs", None, True),
            (single_iso_bound_mz_start_idcs, "Single_Isotope_mz_start", None, True),
            (single_iso_bound_mz_end_idcs, "Single_Isotope_mz_end", None, True),
            (single_iso_bound_rt_start_idcs, "Single_Isotope_rt_start", None, True),
            (single_iso_bound_rt_end_idcs, "Single_Isotope_rt_end", None, True),
            (xic_intenst_idcs, "XICs_intensities", None, True),
            (xic_rt_idcs, "XICs_retention_times", None, True),
        ]

        # Create workbook and set the toplevel header with merged cells
        workbook = xlsxwriter.Workbook(args.output)
        workbook.use_zip64()
        worksheet = workbook.add_worksheet()

        # Write concatenated Headers (merged Header)
        offset = 0
        for idcs, header, _, _ in final_output_list:
            if idcs is None:
                continue
            if isinstance(idcs, int):        
                merge_cells(offset, offset, header, workbook, worksheet)
                offset += 1
                continue
            else: # its a list
                if len(idcs) == 0:
                    continue
                merge_cells(offset, offset+len(idcs)-1, header, workbook, worksheet)
                offset += len(idcs)
                

        # Write column headers:
        offset = 0
        header_format = workbook.add_format({"bold": True, "align": "center", "bottom": 2})
        header_last_format = workbook.add_format({"bold": True, "align": "center", "bottom": 2, "right": 2})

        for idcs, header, width, remove_suffix in final_output_list:
            if idcs is None:
                continue
            if isinstance(idcs, int):
                worksheet.write(1, offset, RAW_header[idcs] if not remove_suffix else RAW_header[idcs].split("_____", 1)[0], header_last_format)
                if width is not None:
                    worksheet.set_column_pixels(offset, offset, width)
                offset += 1
                continue
            else:
                if len(idcs) == 0:  # Redundant, but to keep the code structure the same
                    continue
                for offoff, i in enumerate(idcs):
                    if len(idcs) -1 == offoff:
                        worksheet.write(1, offset + offoff, RAW_header[i] if not remove_suffix else RAW_header[i].split("_____")[0], header_last_format)
                    else:
                        worksheet.write(1, offset + offoff, RAW_header[i] if not remove_suffix else RAW_header[i].split("_____")[0], header_format)
                    if width is not None:
                        worksheet.set_column_pixels(offset + offoff, offset + offoff, width)
                offset += len(idcs)


        data_last_format = workbook.add_format({"right": 2})
        row_offset = 2
        for row in tqdm.tqdm(csv_reader, unit="lines", total=no_lines-1):
            offset = 0
            for idcs, _, _, _ in final_output_list:
                if idcs is None:
                    continue
                if isinstance(idcs, int):
                    worksheet.write(row_offset, offset, row[idcs], data_last_format)
                    offset += 1
                    continue
                else:
                    if len(idcs) == 0:  # Redundant, but to keep the code structure the same
                        continue
                    for offoff, i in enumerate(idcs):
                        if len(idcs) -1 == offoff:
                            worksheet.write(row_offset, offset + offoff, row[i], data_last_format)
                        else:
                            worksheet.write(row_offset, offset + offoff, row[i])

                    offset += len(idcs)
            row_offset += 1

    workbook.close()
