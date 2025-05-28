import argparse

from analysis import cds, clustering, variants

''''
Entry point for the package 
Handles all commands from command line
'''

def main(): 

    
    # Initialize parser object
    parser = argparse.ArgumentParser(prog="apggr", description="A package with a set of modules for the analysis of population genetics of genomic regions targeting bacterial species")

    subparsers = parser.add_subparsers(title="commands", dest="command")
    
    # Coding DNA sequence command and arguments
    cds_command = subparsers.add_parser("cds", help="Extract coding DNA sequence from a genomic region in annotation files")
    cds_command.add_argument("--annotation", required=True, type=str, help="Annotation file path")
    cds_command.add_argument("--region", required=True, type=str, help="Genomic region")
    cds_command.add_argument("--output", required=True, type=str, help="Output directory path")
    
    # Cluster command and arguments
    cluster_command = subparsers.add_parser("cluster", help="Cluster samples in a population using single nucleotide polymorphism (SNP) data")
    cluster_command.add_argument("--params", required=True, type=str, help="DBSCAN model parameters in format 'epsilon:minimum_samples:distance_metric'")
    cluster_command.add_argument("--snpdir", required=True, type=str, help="Directory path for variants in plain text format subset to genomic region")
    cluster_command.add_argument("--outdir", required=True, type=str, help="Output directory path")
        
    # Variant command and arguments
    variant_command = subparsers.add_parser("variants", help="Analysis of variant effect predictions")
    variant_command.add_argument("--effectsdir", required=True, type=str, help="Effects directory")
    variant_command.add_argument("--annotation", required=True, type=str, help="Genome nnotation file")
    variant_command.add_argument("--region", required=True, type=str, help="Genomic region")
    variant_command.add_argument("--outdir", required=True, type=str, help="Output directory path")
                
    args = parser.parse_args()

    # Dispatch to the appropriate command
    if args.command == "cds":
        cds.run(args)
        
    elif args.command == "cluster":
        clustering.run(args)
        
    elif args.command == "variants":
        variants.run(args)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()