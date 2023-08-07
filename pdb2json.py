%%time
from Bio.PDB import PDBParser
import json
import yaml

def pdb_to_json(pdb_file, json_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)
    
    data = {
        "atoms": [],
        "residues": [],
        "chains": []
    }
    
    for model in structure:
        for chain in model:
            data["chains"].append({"id": chain.id})
            for residue in chain:
                data["residues"].append({
                    "id": residue.id[1],
                    "name": residue.resname,
                    "chain": chain.id
                })
                for atom in residue:
                    data["atoms"].append({
                        "name": atom.name,
                        "element": atom.element,
                        "residue_id": residue.id[1],
                        "chain": chain.id
                    })
    
    with open(json_file, "w") as outfile:
        json.dump(data, outfile, indent=4)

def pdb_to_yaml(pdb_file, yaml_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)
    
    data = {
        "atoms": [],
        "residues": [],
        "chains": []
    }
    
    for model in structure:
        for chain in model:
            data["chains"].append({"id": chain.id})
            for residue in chain:
                data["residues"].append({
                    "id": residue.id[1],
                    "name": residue.resname,
                    "chain": chain.id
                })
                for atom in residue:
                    data["atoms"].append({
                        "name": atom.name,
                        "element": atom.element,
                        "residue_id": residue.id[1],
                        "chain": chain.id
                    })
    
    with open(yaml_file, "w") as outfile:
        yaml.dump(data, outfile, default_flow_style=False)

if __name__ == "__main__":
    pdb_file = "2R_80.pdb"
    json_file = "2R_80ns.json"
    yaml_file = "2R_80n.yaml"
    
    pdb_to_json(pdb_file, json_file)
    pdb_to_yaml(pdb_file, yaml_file)
