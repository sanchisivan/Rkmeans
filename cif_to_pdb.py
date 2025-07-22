from Bio.PDB import MMCIFParser, PDBIO
import sys
import os

def convert_cif_to_pdb(cif_path, pdb_path):
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure('X', cif_path)
    io = PDBIO()
    io.set_structure(structure)
    io.save(pdb_path)

if __name__ == "__main__":
    cif_file = sys.argv[1]
    pdb_file = sys.argv[2]
    convert_cif_to_pdb(cif_file, pdb_file)
