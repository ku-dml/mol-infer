# -*- coding: utf-8 -*-

import sys
from pathlib import Path

if len(sys.argv) < 3:
    print("Usage:")
    print("python First_n_SDF_selection.py input_file.sdf number_of_molecules")
    print("Example:")
    print("python First_n_SDF_selection.py qm9_eli.sdf 500")
    sys.exit(1)

input_sdf = Path(sys.argv[1])

try:
    number_of_molecules = int(sys.argv[2])
except ValueError:
    print("Error: Number of molecules must be an integer.")
    sys.exit(1)

if number_of_molecules <= 0:
    print("Error: Number of molecules must be greater than zero.")
    sys.exit(1)

if not input_sdf.exists():
    print(f"Error: File not found: {input_sdf}")
    sys.exit(1)

output_sdf = input_sdf.parent / (
    f"{input_sdf.stem}_first_{number_of_molecules}.sdf"
)

count = 0
molecule_record = []

with open(
    input_sdf,
    "r",
    encoding="utf-8",
    errors="replace"
) as infile, open(
    output_sdf,
    "w",
    encoding="utf-8"
) as outfile:

    for line in infile:
        molecule_record.append(line)

        # Each molecule ends with $$$$
        if line.strip() == "$$$$":
            outfile.writelines(molecule_record)
            molecule_record = []
            count += 1

            if count >= number_of_molecules:
                break

print(f"Completed: {count} molecules saved.")
print(f"Output file: {output_sdf}")