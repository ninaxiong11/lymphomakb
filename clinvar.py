from variant import *

def add_clinvar_data(graph, Involves, num_variants=500):
    batch = graph.batch()
    disease_dict = {}
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
        phenotypes = line[13].strip().split("|")
        build = line[16]
        chr, start, stop = line[18], line[19], line[20]
        # add variants
        if build != "GRCh38":
            continue
        if type == "single nucleotide variant":
            if add_snv(batch, name, gene, build, chr, start):
                count, disease_dict = add_statement(batch, name, clinical_significance, phenotypes, Involves, disease_dict)
                statement_count += count
                variant_count += 1
                snv_count += 1
        if type == "copy number gain":
            if geneid == "-1":
                gene = ""
            if add_copy_gain(batch, name, gene, build, chr):
                count, disease_dict = add_statement(batch, name, clinical_significance, phenotypes, Involves, disease_dict)
                statement_count += count
                variant_count += 1
                copy_gain_count += 1
    infile.close()
    batch.commit()
    print("\nFrom ClinVar:")
    print("Added " + str(variant_count) + " variants")
    print("   " + str(snv_count) + " single nucleotide variants")
    print("   " + str(copy_gain_count) + " copy gains")
    print("Added " + str(statement_count) + " statements")