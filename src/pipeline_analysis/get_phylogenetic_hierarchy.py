import os
from dataclasses import dataclass
from typing import List
from enum import IntEnum

class Level(IntEnum):
    U = 0
    R = 1
    D = 2
    K = 3
    P = 4
    C = 5
    O = 6
    F = 7
    G = 8
    S = 9

@dataclass
class TreeNode:
    name: str
    taxid: str
    level: str
    level_enum : 'Level'
    abundance: int = 0
    acumulated_abundance: int = 0
    children: List['TreeNode'] = None
    parent: 'TreeNode' = None

    def __init__(self, name: str, taxid: str, level: str):
        self.name = name
        self.taxid = taxid
        self.level = level

    def __str__(self):
        parent =  self.parent.name if self.parent else 'None'
        return f"{self.name}, {self.taxid}, {self.level}, {parent}"

    def add_child(self, child_node):
        if self.children is None:
            self.children = []
        self.children.append(child_node)
        child_node.parent = self

    def get_level_enum(self):
        return Level[self.level[0]] if self.level[0] in Level.__members__ else None
    
    def get_sublevel(self):
        return  int(self.level[1:]) if len(self.level) > 1 else 0
    
    def clear_abundance(self):
        self.abundance = 0
        self.acumulated_abundance = 0

    def set_abundance(self, abundance):
        self.abundance = abundance
        self.acumulated_abundance += abundance
        parent = self.parent
        while parent is not None:
            parent.acumulated_abundance += abundance
            parent = parent.parent



def get_files_in_folder(input_path, input_extension):
    print("Start process")
    files_fullpath = []
    for root, dirs, files in os.walk(input_path):
        for file_name in files:
            if file_name.endswith(input_extension):
                file_path = os.path.join(root, file_name)
                files_fullpath.append(file_path)
    return files_fullpath


def read_report_tree(input_file):
    tree_by_name = {}
    tree_by_taxid = {}
    last_parent_node = TreeNode("unclassified", "0", "U")
    last_node = last_parent_node

    with open(input_file, 'r') as file:
        for line in file:
            line = line.strip()
            line_splits = line.split(maxsplit=5)
            if len(line_splits) >= 6:
                level = line_splits[3] 
                taxid = line_splits[4]
                name = line_splits[5]
                new_node = TreeNode(name, taxid, level)

                if taxid == '0':
                    print(f"Invalid taxid: {line}")
                    continue
                if not new_node.get_level_enum():                    
                    print(f"Invalid level: {line}")
                    continue
                if name in tree_by_name:                    
                    print(f"Duplicate name: new='{line}' existing='{tree_by_name[name]}'")
                    continue
                if taxid in tree_by_taxid:                    
                    print(f"Duplicate taxid: new='{line}' existing='{tree_by_taxid[taxid]}")
                    continue

                while (new_node.get_level_enum() < last_node.get_level_enum() or
                       (new_node.get_level_enum() == last_node.get_level_enum() 
                        and new_node.get_sublevel() < last_node.get_sublevel())):
                    last_node = last_node.parent

                if new_node.get_level_enum() > last_node.get_level_enum():
                    parent_node = last_node
                elif new_node.get_level_enum() == last_node.get_level_enum():
                    if new_node.get_sublevel() == last_node.get_sublevel():
                        parent_node = last_node.parent
                    else:
                        parent_node = last_node

                parent_node.add_child(new_node)
                tree_by_taxid[taxid] = new_node
                tree_by_name[name] = new_node
                last_node = new_node
                #print(f"New node: {new_node}")
            else:
                print(f"Invalid line: {line}")
    print(f"Length report taxid tree: {len(tree_by_name)}")
    return tree_by_taxid

def clear_abundance_from_tree(dict_tree):
    # Loop throught all tree nodes and clear ir
    for value in dict_tree.values():
        value.clear_abundance()



def count_abundance_by_species(accession_taxid_counts, accession_abundance, tree_by_taxid, output_file):
    with open(output_file, 'w') as file:
        file.write("parent_taxid,level,taxid,name,kraken_classifed_reads,NT_rPM\n")

    clear_abundance_from_tree(tree_by_taxid)
    for taxid_counts in accession_taxid_counts.values():
        for taxid, counts in taxid_counts.items():
            if taxid not in tree_by_taxid:
                print(f"Invalid taxid: {taxid}")
                continue
            tree_by_taxid[taxid].set_abundance(counts)

    total_reads = float(sum(accession_abundance.values()))
    for taxid, node in tree_by_taxid.items():
        if node.level == 'S' or node.level == 'G':
            level = node.get_level_enum() - Level.G + 1
            abundance = node.acumulated_abundance
            nt_rpm = int((abundance/total_reads)*1000000)
            with open(output_file, 'a') as file:
                file.write(f"{node.parent.taxid},{level},{taxid},{node.name},{abundance},{nt_rpm}\n")



def count_abundance_by_level(accession_taxid_counts, accession_abundance, accession_metadata, tree_by_taxid, output_file):
    with open(output_file, 'w') as file:
        file.write("read_accession_id,read_accession_total_count,root,kraken_classified_root")
        file.write(",superkingdom,kraken_classified_superkingdom,phylum,kraken_classified_phylum")
        file.write(",class,kraken_classified_class,order,kraken_classified_order,family,kraken_classified_family")
        file.write(",genus,kraken_classified_genus,species,kraken_classified_species\n")

    for accession_id, taxid_counts in accession_taxid_counts.items():
        clear_abundance_from_tree(tree_by_taxid)
        for taxid, counts in taxid_counts.items():
            if taxid not in tree_by_taxid:
                print(f"Invalid taxid: {taxid}")
                continue
            tree_by_taxid[taxid].set_abundance(counts)

        accession_lineage = accession_metadata[accession_id]
        lineage_output = f"{accession_id},{accession_abundance[accession_id]}"
        for taxid, name in accession_lineage:
            if taxid not in tree_by_taxid:
                lineage_output += f",{name},0"
            else:
                lineage_output += f",{name},{tree_by_taxid[taxid].acumulated_abundance}"                
        with open(output_file, 'a') as file:
            file.write(f"{lineage_output}\n")



def get_accession_abundance(filename):
    accession_abundance = {}

    with open(filename, 'r') as file:
        header = True
        for line in file:
            if header:
                header = False
                continue

            line_splits = line.strip().split(',')
            accession_id = line_splits[0]
            abundance = int(line_splits[1])
            if accession_id in accession_abundance:
                print(f"Duplicated accession id {accession_id}")
                continue

            accession_abundance[accession_id] = abundance

    return accession_abundance


def get_accession_taxid_abundance(filename):
    accession_taxid_abundance = {}

    with open(filename, 'r') as file:
        header = True
        for line in file:
            if header:
                header = False
                continue

            line_splits = line.strip().split(',')
            accession_id = line_splits[0]
            taxid = line_splits[1]
            abundance = int(line_splits[2])

            if accession_id not in accession_taxid_abundance:
               accession_taxid_abundance[accession_id] = {}
            
            accession_taxid_abundance[accession_id][taxid] = abundance

    return accession_taxid_abundance


def get_accession_taxid_by_level(filename):
    accession_taxid_by_level = {}

    with open(filename, 'r') as file:
        header = True
        for line in file:
            if header:
                header = False
                continue

            line_splits = line.strip().split(',')
            accession_id = line_splits[0]
            lineage = line_splits[2:9]
            taxids = line_splits[9:]

            if accession_id in accession_taxid_by_level:
                print(f"Duplicated accession id {accession_id}")
                continue

            lineage_list = []
            lineage_list.append(('1', 'root'))
            lineage_list.extend( [(taxids[i], lineage[i]) for i in range(0, 7)])
            accession_taxid_by_level[accession_id] = lineage_list

    return accession_taxid_by_level


def main():
    #folder_path = sys.argv[1] if len(sys.argv) > 1 else ''
    input_extension = 'mock01.report'
    input_path = r"C:\Users\pablo\Documents\github\taxonomy_analysis\data"
    output_path = r"C:\Users\pablo\Documents\github\taxonomy_analysis\results"

    if not os.path.isdir(input_path):
        print("Invalid input path")
        return

    all_files = get_files_in_folder(input_path, input_extension)
    print(all_files)

    for file in all_files:
        print(f"Analyzing file: {file}")
        tree_by_taxid = read_report_tree(file)
        filename = os.path.basename(file).split(".")[0]
        
        metadata_file = os.path.join(output_path, filename + "_metadata.csv")
        accession_metadata = get_accession_taxid_by_level(metadata_file)

        accession_taxid_abundance_file = os.path.join(output_path, filename + "_out.csv")
        accession_taxid_counts = get_accession_taxid_abundance(accession_taxid_abundance_file)

        accession_abundance_file = os.path.join(output_path, filename + "_accession_abundance.csv")
        accession_abundance = get_accession_abundance(accession_abundance_file)

        abundance_by_level_file = os.path.join(output_path, filename + "_level_abundance.csv")
        count_abundance_by_level(accession_taxid_counts, accession_abundance, accession_metadata, tree_by_taxid, abundance_by_level_file)

        abundance_by_species_file = os.path.join(output_path, filename + "_species_abundance.csv")
        count_abundance_by_species(accession_taxid_counts, accession_abundance, tree_by_taxid, abundance_by_species_file)

        #output_file = os.path.join(output_path, filename + "_tree_out.csv")
        #write_output_file(read_counts, output_file)


if __name__ == '__main__':
    main()