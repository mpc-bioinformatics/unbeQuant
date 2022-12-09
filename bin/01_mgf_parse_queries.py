#!/bin/env python

import sys

input_mgf = sys.argv[1]
output_csv = sys.argv[2]
limit_da =  int(sys.argv[3])
ppm = int(sys.argv[4])

HYDROGEN_MONO_MASS = 1.007825035
offset = -18.0105647 # subtract water, since proteingraphs do not encode this information!


if __name__ == "__main__":
    with open(output_csv, "w") as out_file, open(input_mgf, "r") as in_file:

        entries = []
        in_entry = False
        retrieved_pep_mass = False
        retrieved_charge = False
        for line in in_file:

            if line.startswith("BEGIN IONS"):
                in_entry = True
                continue

            if in_entry:
                if line.startswith("PEPMASS="):
                    pepmass = float(line[len("PEPMASS="):-1].split(" ")[0])
                    retrieved_pep_mass = True
                    continue
                if line.startswith("CHARGE="):
                    charge = int(line[len("CHARGE="):-1].replace("+", ""))
                    retrieved_charge = True
                    continue

            if line.startswith("END IONS"):
                # Special check if pepmass was retrieved
                if retrieved_charge and retrieved_charge:
                    # Do here the calculation of lower and upper limit
                    da = (float(pepmass) * float(charge)) - (HYDROGEN_MONO_MASS * float(charge))
                    da = da + offset

                    lower = da - (da / 1000000) * ppm
                    upper = da + (da / 1000000) * ppm

                    entries.append((lower, upper))

                in_entry = False
                retrieved_pep_mass = False
                retrieved_charge = False


        # Write only the sorted and limited queries
        for l, u in sorted(entries):
            if l < limit_da:
                out_file.write(str(l) + "," + str(u) + "\n")
