# RNAnneal_LiLab_2026-04-24

This package contains RNAnneal representative structures for customer delivery.

## Contents

- `sequence/input.fa`
  Input RNA sequences.

- `sequence/reps/`
  Representative PDB models using the FASTA/workdir IDs.

- `manifest.csv`
  Per-sequence packaging summary including sequence ID, sequence length,
  number of copied representatives, packaging status, source directory, and
  source zip path.

- `chimerax/color_ie_bins.cxc`
  ChimeraX coloring script copied from the example package.

- `chimerax/open_and_color_all.cxc`
  Convenience script to open all representative PDBs and then load the
  coloring script.

## Notes

Only representative PDB files from each `results.zip` archive were included
when results archives were present.
