# RNAnneal_LiLab_2026-03-27

This package contains RNAnneal representative structures for customer delivery.

## Contents

- `sequence/input.fa`
  Input RNA sequences.

- `sequence/reps/`
  Representative PDB models using the original FASTA/workdir IDs, for example:
  - `id_4_row5_rank01.pdb`
  - `id_4_row5_rank02.pdb`
  - `id_10_row11_rank01.pdb`

- `manifest.csv`
  Per-sequence packaging summary including:
  - sequence ID
  - sequence length
  - number of copied representatives
  - packaging status
  - original source directory
  - source zip path

- `chimerax/color_ie_bins.cxc`
  ChimeraX coloring script copied from the RNAnneal repo.

- `chimerax/open_and_color_all.cxc`
  Convenience script to open all representative PDBs and then load the coloring script.

## Notes

Only representative PDB files from each `results.zip` archive were included.
