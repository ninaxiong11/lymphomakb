from pyorient.ogm import Graph, Config
from schema import *
from variant import *
from clinvar import add_clinvar_data
from ncit import add_disease_tree
import argparse

parser = argparse.ArgumentParser(description="Setup graph KB")
parser.add_argument("-n", '--num_variants', type=int, default=500)
args = parser.parse_args()

def main():
    # when beginning, set initial drop to false
    graph = Graph(Config("localhost", 2424, "root", "root", "test", initial_drop=True))
    print("Created graph!")
    Involves, SubclassOf = initialize_graph_schema(graph)
    # add_disease_tree(graph, SubclassOf)
    add_clinvar_data(graph, Involves, args.num_variants)

if __name__ == "__main__":
    main()