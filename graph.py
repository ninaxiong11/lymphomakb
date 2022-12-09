from pyorient.ogm import Graph, Config
from schema import *
from variant import *
from clinvar import add_clinvar_data
from ncit import add_disease_tree

def main():
    # when beginning, set initial drop to false
    graph = Graph(Config("localhost", 2424, "root", "root", "test", initial_drop=True))
    print("Created graph!")
    Variant, SNV, CopyGain, CopyLoss, Statement, Disease, Drug, Involves, SubclassOf = initialize_graph_schema(graph)
    # add_disease_tree(graph, SubclassOf)
    add_clinvar_data(graph, Involves, 500)

if __name__ == "__main__":
    main()