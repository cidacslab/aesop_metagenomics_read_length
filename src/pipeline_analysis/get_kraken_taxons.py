import os, sys, csv
import kraken_report_parser as KrakenParser


def main():
    report_file = "../functional_redundancy/data/kraken_inspect_ids.txt"
    # report_file = "results/pluspfp_20230605_inspect.txt"
    taxid_file = "results/new_compositions_taxids.txt"
    composition_file = "data/AESOP_AMVB_MOCKS_with_accession_final.csv"
    # output_file = "results/complete_kraken_summary_pfpdb.csv"
    output_file = "results/include_in_krakendb.txt"

    root_node, tree_by_taxid = KrakenParser.load_kraken_report_tree(report_file)

    # all_nodes = root_node.get_all_nodes()

    # output_content = "taxid,level,level_enum,level_count,name,parent_taxid,parent_level,parent_name\n"
    # output_content += "\n".join([str(node) for node in all_nodes])

    # with open(output_file, 'w') as file:
    #     file.write(output_content)

    count = 0
    taxid_not_found = set()
    with open(taxid_file, "r") as file:
        for line in file:
            taxid = line.strip()
            #print(line)
            if taxid not in tree_by_taxid:
                taxid_not_found.add(taxid)
                #print(f"Node not found for {taxid}")
                count += 1
    print(f"\nNot found {count} taxids!")


    accession_by_taxid = {}
    with open(composition_file, "r") as file:
        csv_reader = csv.reader(file, delimiter=",")
        row = next(csv_reader)
        for row in csv_reader:
            accession = row[0].strip()
            taxid = row[2].strip()
            accession_by_taxid[taxid] = accession

    output_content = ""
    for id in taxid_not_found:
        output_content += f"{accession_by_taxid[id]}\n"
        print(f"Accession for taxid {id}: {accession_by_taxid[id]}")

    with open(output_file, 'w') as file:
        file.write(output_content)



if __name__ == '__main__':
    main()