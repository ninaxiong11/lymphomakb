from variant import *

def add_clinvar_data(graph, Involves, filepath, num_variants=500):
    batch = graph.batch()
    disease_dict = {}
    statement_count = 0
    variant_count = {"single nucleotide variants": 0, 
                    "copy number gains": 0, 
                    "copy number losses": 0,
                    "indels": 0,
                    "deletions": 0,
                    "insertions": 0
                    }
    infile = open(filepath) # /Users/ninaxiong/projects/orientdb/clinvar/variant_summary.txt
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
                variant_count["single nucleotide variants"] += 1
                statement_count += count
        if type == "copy number gain":
            if geneid == "-1":
                gene = "-"
            if add_copy_gain(batch, name, gene, build, chr):
                count, disease_dict = add_statement(batch, name, clinical_significance, phenotypes, Involves, disease_dict)
                variant_count["copy number gains"] += 1
                statement_count += count
        if type == "copy number loss":
            if geneid == "-1":
                gene = "-"
            if add_copy_loss(batch, name, gene, build, chr):
                count, disease_dict = add_statement(batch, name, clinical_significance, phenotypes, Involves, disease_dict)
                variant_count["copy number losses"] += 1
                statement_count += count
        if type == "Indel":
            if add_indel(batch, name, gene, build, chr, start, stop):
                count, disease_dict = add_statement(batch, name, clinical_significance, phenotypes, Involves, disease_dict)
                variant_count["indels"] += 1
                statement_count += count
        if type == "Deletion":
            if add_deletion(batch, name, gene, build, chr, start, stop):
                count, disease_dict = add_statement(batch, name, clinical_significance, phenotypes, Involves, disease_dict)
                variant_count["deletions"] += 1
                statement_count += count
        if type == "Insertion":
            if add_insertion(batch, name, gene, build, chr, start, stop):
                count, disease_dict = add_statement(batch, name, clinical_significance, phenotypes, Involves, disease_dict)
                variant_count["insertions"] += 1
                statement_count += count
    infile.close()
    batch.commit()
    print("\nFrom ClinVar:")
    print("Added " + str(sum(variant_count.values())) + " variants")
    for k,v in variant_count.items():
        print("   " + str(v) + " " + k)
    print("Added " + str(statement_count) + " statements")