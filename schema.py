from pyorient.ogm.declarative import declarative_node, declarative_relationship
from pyorient.ogm.property import String, Integer

def initialize_graph_schema(graph):
    
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
        chromosome = String()
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
        chromosome = String()
        genomic_start = Integer()
        genomic_end = Integer()
        cytogenetic_loc = String()

    class CopyLoss(Variant):
        element_type = "copyloss"
        element_plural = "copylosses"
        build = String()
        gene = String()
        chromosome = String()
        genomic_start = Integer()
        genomic_end = Integer()
        cytogenetic_loc = String()

    class Indel(Variant):
        element_type = "indel"
        element_plural = "indels"
        build = String()
        gene = String()
        transcript = String()
        chromosome = String()
        genomic_start = Integer()
        genomic_end = Integer()
        coding_start = String()
        coding_end = String()
        ins_nt = String()
    
    class Deletion(Variant):
        element_type = "deletion"
        element_plural = "deletions"
        build = String()
        gene = String()
        transcript = String()
        chromosome = String()
        genomic_start = Integer()
        genomic_end = Integer()
        coding_start = String()
        coding_end = String()

    class Insertion(Variant):
        element_type = "insertion"
        element_plural = "insertions"
        build = String()
        gene = String()
        transcript = String()
        chromosome = String()
        genomic_start = Integer()
        genomic_end = Integer()
        coding_start = String()
        coding_end = String()
        ins_nt = String()
        
    class Statement(Node):
        element_type = "statement"
        element_plural = "statements"
        source = String()
        significance = String()

    class Disease(Node):
        element_type = "disease"
        element_plural = "diseases"
        name = String(nullable=False)
        source = String()
        synonyms = String()

    class Drug(Node):
        element_type = "drug"
        element_plural = "drugs"
        name = String(nullable=False, unique=True)

    class Involves(Relationship):
        label = "involves"

    class SubclassOf(Relationship):
        label = "subclassOf"

    # class Property(Node):
    #     pass

    # class HasProperty(Relationship):
    #    pass

    # class Infers(Relationship):
    #     pass

    graph.create_all(Node.registry)
    graph.create_all(Relationship.registry)

    return Involves, SubclassOf