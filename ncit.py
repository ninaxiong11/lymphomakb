from variant import *

def get_synonyms(records, code):
    record = [r for r in records if r["code"] == code]
    if record:
        return " | ".join(record[0]["synonyms"]) # lower()

def get_parent_name(records, parent_code):
    parent = [r for r in records if r["code"] == parent_code]
    if parent:
        return parent[0]["name"]

def add_disease_tree(graph, SubclassOf):
    batch = graph.batch()
    infile = open("/Users/ninaxiong/Downloads/Thesaurus.txt") # change filepath
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
        synonyms = synonyms.split("|")
        if type != "Neoplastic Process":
            continue
        if not name:
            name = synonyms[0]
        records.append({"code": code, "parents": parents, "name": name, "synonyms": synonyms})
    infile.close()
    for r in records:
        child_name = r["name"]
        try:
            batch[:clean_string(child_name)]
        except:
            batch[clean_string(child_name)] = batch.diseases.create(
                source = "ncit",
                code = r["code"],
                name = child_name,
                synonyms = get_synonyms(records, r["code"])
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
                        name = parent_name,
                        synonyms = get_synonyms(records, parent_code)
                    )
                batch[:] = batch[:clean_string(child_name)](SubclassOf) > batch[:clean_string(parent_name)]
    batch.commit()
    print("\nFrom NCIt:")
    print("Added {} diseases".format(len(records)))