from pyorient.ogm import Graph, Config
from pyorient.ogm.declarative import declarative_node, declarative_relationship
from pyorient.ogm.property import String, Integer
from types import *
import re

graph = Graph(Config("localhost", 2424, "root", "root", "test", initial_drop=True))

Node = declarative_node()
Relationship = declarative_relationship()

class Variant(Node):
    element_type = "variant"
    element_plural = "variants"
    name = String()
    source = String()

class SNV(Variant):
    element_type = "snv"
    element_plural = "snvs"
    build = String()
    gene = String()
    transcript = String()
    genomic_pos = Integer()
    coding_pos = Integer()
    ref_nt = String()
    alt_nt = String()
    protein_pos = Integer()
    ref_aa = String()
    alt_aa = String()

class CopyGain(Variant):
    element_type = "copygain"
    element_plural = "copygains"
    build = String()
    gene = String()
    copies = Integer()
    genomic_loc = String()
    cytogenetic_loc = String()

class Statement(Node):
    element_type = "statement"
    element_plural = "statements"
    source = String()
    significance = String()

class Disease(Node):
    element_type = "disease"
    element_plural = "diseases"
    name = String(nullable=False)

class Drug(Node):
    element_type = "drug"
    element_plural = "drugs"
    name = String(nullable=False, unique=True)

class Involves(Relationship):
    label = "involves"

# class Property(Node):
#     pass

# class HasProperty(Relationship):
#    pass

# class SubclassOf(Relationship):
#     pass

# class Infers(Relationship):
#     pass

graph.create_all(Node.registry)
graph.create_all(Relationship.registry)

batch = graph.batch()

count = 0
snv_count = 0
copy_gain_count = 0

def get_mut_info(change, type):
    if type == "nt":
        match = re.match("(.+)\.(\d+)(.+)>(.+)", change)
        if not match:
            return None
        match = match.groups()
        pos = match[1]
        ref = match[2]
        alt = match[3]
    elif type == "protein":
        match = re.match("(.+)\.(\D+)(\d+)(\D+)", change)
        if not match:
            return None
        match = match.groups()
        pos = match[2]
        ref = match[1]
        alt = match[3]  
    return pos, ref, alt

def add_snv(name, gene, build, genomic_pos):
    global count, snv_count
    match = re.match("(.+)\((.+)\)\:(.+).\((.+)\)", name)
    if not match:
        return
    match = match.groups()
    transcript = match[0]
    gene = match[1]
    if match[2].startswith("c"):
        coding_change = match[2]
        if not get_mut_info(coding_change, "nt"):
            return
        else:
            coding_pos, ref_nt, alt_nt = get_mut_info(coding_change, "nt")
    else:
        return
    protein_change = match[3]
    if not get_mut_info(protein_change, "protein"):
        return
    else:
        protein_pos, ref_aa, alt_aa = get_mut_info(protein_change, "protein")
    variant_name = "snv" + str(snv_count)
    batch[variant_name] = batch.snvs.create(
        name = gene + coding_change,
        source = "clinvar",
        build = build,
        gene = gene,
        transcript = transcript,
        genomic_pos = int(genomic_pos),
        coding_pos = int(coding_pos),
        ref_nt = ref_nt,
        alt_nt = alt_nt,
        protein_pos = int(protein_pos),
        ref_aa = ref_aa,
        alt_aa = alt_aa,
    )
    count += 1
    snv_count += 1

def add_copy_gain(name, gene, build):
    global count, copy_gain_count
    match = re.match("(.+)\s(.+)\((.+)\)x(\d+)", name)
    if not match:
        return
    match = match.groups()
    cytogenetic_loc = match[1]
    genomic_loc = match[2]
    copies = match[3]
    variant_name = "copygain" + str(copy_gain_count)
    batch[variant_name] = batch.copygains.create(
        source = "clinvar",
        build = build,
        gene = gene,
        copies = copies,
        genomic_loc = genomic_loc,
        cytogenetic_loc = cytogenetic_loc
    )
    count += 1
    copy_gain_count += 1

def add_statement(clinsig, phenotypes, variant_type, variant_count):
    variant_name = variant_type + str(variant_count - 1)
    statement_count = 0
    disease_count = 0
    for pheno in phenotypes:
        statement_name = variant_name + "statement" + str(statement_count)
        batch[statement_name] = batch.statements.create(
            source = "clinvar",
            significance = clinsig
        )
        disease_name = variant_name + "disease" + str(disease_count)
        batch[disease_name] = batch.diseases.create(
            name = pheno
        )
        batch[:] = batch[:statement_name](Involves) > batch[:disease_name]
        batch[:] = batch[:statement_name](Involves) > batch[:variant_name]
        statement_count += 1
        disease_count += 1

infile = open("../clinvar/variant_summary.txt")
for i in range(500): #3049146
    line = infile.readline()
    line = line.strip().split("\t")
    #
    type, name, geneid, gene = line[1], line[2], line[3], line[4]
    chr, start, stop = line[18], line[19], line[20]
    clinical_significance = line[6]
    phenotypes = line[13].split(", ")
    build = line[16]
    chr, start, stop = line[18], line[19], line[20]
    #
    if build != "GRCh38":
        continue
    if type == "single nucleotide variant":
        add_snv(name, gene, build, start)
        add_statement(clinical_significance, phenotypes, "snv", snv_count)
    # if type == "copy number gain":
    #     if geneid == "-1":
    #         gene = ""
    #     add_copy_gain(name, gene, build)
    #     add_statement(clinical_significance, phenotypes, "copygain", copy_gain_count)
infile.close()

batch.commit()

print("Finished!")
print("Added " + str(count) + " variants:")
print("   " + str(snv_count) + " single nucleotide variants")
print("   " + str(copy_gain_count) + " copy gains")




def main():
    pass

if __name__ == "__main__":
    main()

exit()

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
    synonyms = line[3]
    definition = line[4]
    name = line[5]
    type = line[7]
    if type != "Neoplastic Process":
        continue
    if not name:
        name = synonyms.split("|")[0]
    records.append({"code": code, "parents": parents, "name": name})
infile.close()

def get_parent_name(parent_code):
    parent = [r for r in records if r["code"] == parent_code]
    if parent:
        return parent[0]["name"]

for r in records:
    child_name = r["name"]
    result = client.command("SELECT FROM Disease WHERE name='{}'".format(child_name)) # result.oRecordData
    if not result:
        client.command("CREATE VERTEX Disease SET name='{}', code='{}'".format(child_name, r["code"]))
    for parent_code in r["parents"]:
        parent_name = get_parent_name(parent_code)
        if parent_name:
            result = client.command("SELECT FROM Disease WHERE name='{}'".format(parent_name))
            if not result:
                client.command("CREATE VERTEX Disease SET name='{}', code='{}'".format(parent_name, parent_code))
            try:
                client.command("CREATE EDGE subclassOf FROM (SELECT FROM Disease WHERE name='{}') TO (SELECT FROM Disease WHERE name='{}')".format(child_name, parent_name))
            except:
                continue
            print(child_name + " is a subclass of " + parent_name)
    
print("Finished!")