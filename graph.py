from pyorient.ogm import Graph, Config
from schema import *
from variant import *

def add_clinvar_data(graph, Involves, num_variants=500):
    batch = graph.batch()
    statement_count = 0
    variant_count = 0
    snv_count = 0
    copy_gain_count = 0
    infile = open("/Users/ninaxiong/projects/orientdb/clinvar/variant_summary.txt")
    for i in range(num_variants):
        line = infile.readline()
        line = line.strip().split("\t")
        # extract data
        type, name, geneid, gene = line[1], line[2], line[3], line[4]
        chr, start, stop = line[18], line[19], line[20]
        clinical_significance = line[6]
        phenotypes = line[13].split(", ")
        build = line[16]
        chr, start, stop = line[18], line[19], line[20]
        # add variants
        if build != "GRCh38":
            continue
        if type == "single nucleotide variant":
            if add_snv(batch, name, gene, build, start, snv_count):
                statement_count += add_statement(batch, name, clinical_significance, phenotypes, Involves)
                variant_count += 1
                snv_count += 1
        if type == "copy number gain":
            if geneid == "-1":
                gene = ""
            if add_copy_gain(batch, name, gene, build, copy_gain_count):
                statement_count += add_statement(batch, name, clinical_significance, phenotypes, Involves)
                variant_count += 1
                copy_gain_count += 1
    infile.close()
    batch.commit()
    print("\nFrom ClinVar:")
    print("Added " + str(variant_count) + " variants")
    print("   " + str(snv_count) + " single nucleotide variants")
    print("   " + str(copy_gain_count) + " copy gains")
    print("Added " + str(statement_count) + " statements")

def get_parent_name(records, parent_code):
    parent = [r for r in records if r["code"] == parent_code]
    if parent:
        return parent[0]["name"]

def add_disease_tree(graph, SubclassOf):
    batch = graph.batch()
    infile = open("/Users/ninaxiong/Downloads/Thesaurus.txt")
    records = []
    for line in infile:
        if "lymphoma" not in line or "Lymphoma" not in line:
            continue
        line = line.strip().split("\t")
        if len(line) < 8:
            continue
        code = line[0]
        parents = line[2].split("|")
        synonyms, name, type = line[3], line[5], line[7]
        if type != "Neoplastic Process":
            continue
        if not name:
            name = synonyms.split("|")[0]
        records.append({"code": code, "parents": parents, "name": name})
    infile.close()
    for r in records:
        child_name = r["name"]
        try:
            batch[:clean_string(child_name)]
        except:
            batch[clean_string(child_name)] = batch.diseases.create(
                source = "ncit",
                code = r["code"],
                name = child_name
            )
        for parent_code in r["parents"]:
            parent_name = get_parent_name(records, parent_code)
            if parent_name:
                try:
                    batch[:clean_string(parent_name)]
                except:
                    batch[clean_string(parent_name)] = batch.diseases.create(
                        source = "ncit",
                        code = parent_code,
                        name = parent_name
                    )
                batch[:] = batch[:clean_string(child_name)](SubclassOf) > batch[:clean_string(parent_name)]
    batch.commit()
    print("\nFrom NCIt:")
    print("Added {} diseases".format(len(records)))

def main():
    # when beginning, set initial drop to false
    graph = Graph(Config("localhost", 2424, "root", "root", "test", initial_drop=True))
    print("Created graph!")
    Variant, SNV, CopyGain, Statement, Disease, Drug, Involves, SubclassOf = initialize_graph_schema(graph)
    add_clinvar_data(graph, Involves)
    # add_disease_tree(graph, SubclassOf)

if __name__ == "__main__":
    main()