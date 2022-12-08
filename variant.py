import re

def clean_string(string):
    return re.sub("\W|_", "", string)

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

def add_snv(batch, name, gene, build, chr, genomic_pos):
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
    batch[clean_string(name)] = batch.snvs.create(
        source = "clinvar",
        name = gene + coding_change,
        build = build,
        gene = gene,
        transcript = transcript,
        chromosome = int(chr),
        genomic_pos = int(genomic_pos),
        coding_pos = int(coding_pos),
        ref_nt = ref_nt,
        alt_nt = alt_nt,
        protein_pos = int(protein_pos),
        ref_aa = ref_aa,
        alt_aa = alt_aa,
    )
    return True

def add_copy_gain(batch, name, gene, build, chr):
    match = re.match("(.+)\s(.+)\((.+)\)x(\d+)", name)
    if not match:
        return
    match = match.groups()
    cytogenetic_loc = match[1]
    genomic_loc = match[2]
    copies = match[3]
    batch[clean_string(name)] = batch.copygains.create(
        source = "clinvar",
        build = build,
        gene = gene,
        copies = copies,
        chromosome = int(chr),
        genomic_loc = genomic_loc,
        cytogenetic_loc = cytogenetic_loc
    )
    return True

def add_statement(batch, variant_name, clinsig, phenotypes, statement_type, disease_dict):
    statement_count = 0
    for pheno in phenotypes:
        if pheno not in disease_dict.keys():
            disease_dict[pheno] = "disease" + str(len(disease_dict))
        statement_name = clean_string(variant_name) + clean_string(pheno)
        disease_name = disease_dict[pheno] # clean_string(pheno)
        batch[statement_name] = batch.statements.create(
            source = "clinvar",
            significance = clinsig
        )
        try:
            batch[:disease_name]
        except:
            batch[disease_name] = batch.diseases.create(
                source = "clinvar",
                name = pheno
            )
        batch[:] = batch[:statement_name](statement_type) > batch[:disease_name]
        batch[:] = batch[:statement_name](statement_type) > batch[:clean_string(variant_name)]
        statement_count += 1
    return statement_count, disease_dict