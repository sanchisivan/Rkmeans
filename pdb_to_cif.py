from Bio.PDB import PDBParser, MMCIFIO
import sys
import os

def convert_pdb_to_cif(pdb_path, cif_path):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('X', pdb_path)
    io = MMCIFIO()
    io.set_structure(structure)
    io.save(cif_path)

if __name__ == "__main__":
    pdb_file = sys.argv[1]
    cif_file = sys.argv[2]
    convert_pdb_to_cif(pdb_file, cif_file)
