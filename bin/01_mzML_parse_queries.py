#!/bin/env python

from xml.etree import ElementTree
import sys

input_mzML = sys.argv[1]
output_csv = sys.argv[2]
limit_da =  int(sys.argv[3])
ppm = int(sys.argv[4])

HYDROGEN_MONO_MASS = 1.007825035

offset = -18.0105647 # subtract water, since proteingraphs do not encode this information!


if __name__ == "__main__":
    with open(output_csv, "w") as out_file:
        entries = []
        tree = ElementTree.parse(input_mzML)
        root = tree.getroot()

        mzml_node = root.find(".//{http://psi.hupo.org/ms/mzml}mzML")
        mzml_run = mzml_node.find(".//{http://psi.hupo.org/ms/mzml}run")
        mzml_spectrumlist = mzml_run.find(".//{http://psi.hupo.org/ms/mzml}spectrumList")

        for spec in mzml_spectrumlist:


            prec_list = spec.find(".//{http://psi.hupo.org/ms/mzml}precursorList")
            if prec_list is None:
                continue

            for prec in prec_list.findall(".//{http://psi.hupo.org/ms/mzml}precursor"):
                ions = prec.find(".//{http://psi.hupo.org/ms/mzml}selectedIonList")

                mz_got, charge_got = False, False
                for cvParam in ions[0].findall(".//{http://psi.hupo.org/ms/mzml}cvParam"):
                    if cvParam.get("name") == "selected ion m/z":
                        mz = cvParam.get("value")
                        mz_got = True

                    if cvParam.get("name") == "charge state":
                        charge = cvParam.get("value")
                        charge_got = True

                if mz_got and charge_got:
                    # Do here the calculation of lower and upper limit
                    da = (float(mz) * float(charge)) - (HYDROGEN_MONO_MASS * float(charge))
                    da = da + offset

                    lower = da - (da / 1000000) * ppm
                    upper = da + (da / 1000000) * ppm

                    entries.append((lower, upper))


        for l, u in sorted(entries):
            if l < limit_da:
                out_file.write(str(l) + "," + str(u) + "\n")
