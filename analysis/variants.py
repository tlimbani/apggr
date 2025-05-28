
"""
Gets all the predicted variant effects for a region across all the sequnces
Summary tables for various pivots
These pivots can later be visualized with the visualize module
"""

from collections import Counter
from pandas import DataFrame
from fnmatch import fnmatch
from os import listdir
from os.path import join, splitext

from utils import save_to_csv
from analysis.cds import CDS

class Variants:
    def __init__(self, effects_dir, annotation, region, out_dir):
        self.effectsdir = effects_dir
        self.annotation = annotation
        self.region = region
        self.outdir = out_dir
        
        
    def merge_variant_effects(self):
        """
        Frequency of SNPs consequence overall
        """
        variant_effects_list = []

        # For each variant text file, add to the variants list (CHROM, POS, REF, ALT)
        for effects_file in listdir(self.effectsdir):
            
            # Get text files for varients
            if fnmatch(effects_file, "*.txt") and not fnmatch(effects_file, "*warnings.txt"):
                
                # Full path
                annotation_path = join(self.effectsdir, effects_file)
                
                # Get file name only without the extension
                variant_name = splitext(effects_file)[0]
                
                # Open file and read content
                with open(annotation_path, 'r') as file:
                    rows = file.readlines()[1:]
                    
                    
                    for row in rows:
                        if not row.startswith("#"):
                            effects_list_row = []
                            
                            fields = row.strip().split('\t')
                            
                            # SNP position
                            snp_location=fields[1].split(':')[1]
                            
                            # conseuqnce
                            consequence=fields[6]
                            
                            # Sequence name
                            sequence_name=variant_name
                                                       
                            effects_list_row.append(sequence_name)
                            effects_list_row.append(snp_location)
                            effects_list_row.append(consequence)
                    
                            variant_effects_list.append(effects_list_row)
        return variant_effects_list
                            
    
    def get_effects_overall(self, merged_effects):
        """
        Across all samples in the population

        Parameters
        ----------
        merged_effects : TYPE
            DESCRIPTION.

        Returns
        -------
        effects_frequencies : TYPE
            DESCRIPTION.

        """
        consequence_list = [effects_[2] for effects_ in merged_effects]
        effects_frequencies = Counter(consequence_list)
        
        return effects_frequencies
                              
    
    def get_effects_by_gene(self, merged_effects, cds):
        """
        Frequency of SNPs consequence by coding DNA sequence
        Tag each effect for the CDS it belongs to
        """
              
        effects_merge_gene_product = []
        
        
        for effect_row in merged_effects:
        
            snp_position_lookup = int(effect_row[1])
                       
            matcher = cds[(cds['cds_start'] <= snp_position_lookup) & (cds['cds_end'] >= snp_position_lookup)]
            
            if not matcher.empty:
                product_description = matcher['product_description'].iloc[0]
            
            else: 
                product_description = ""
            
            
            effect_row.append(product_description)
            
            effects_merge_gene_product.append(effect_row)
                
        # Get all unique possible effects
        unique_effects = {consequences[2] for consequences in effects_merge_gene_product}
                
        gene_products = cds['product_description'].replace('', None).dropna().tolist()
        
        effects_summary = {}
        
        # For each consequence count how many in each gene product
        for consequence in unique_effects:
            
            consequence_dict = {key: 0 for key in gene_products}
            
            # Retrieves all similar consequences but different gene products
            filtered_effects = [row for row in effects_merge_gene_product if row[2] == consequence]
            
            gene_product_column = 3
            
            filtered_effects_no_blanks = [row for row in filtered_effects if row[gene_product_column] not in ("", None)]
                                        
            # Get all gene products in a list
            gene_product_values = [row[3] for row in filtered_effects_no_blanks if row[3] not in [None, ''] and row[3] in gene_products]
            
            # Get the frequencies for each gene product 
            frequencies = Counter(gene_product_values)
            
            for key in consequence_dict:
                if key in frequencies:
                    consequence_dict[key] = frequencies[key]
             
            effects_summary[consequence] = consequence_dict
         
        summary_effects_df = DataFrame.from_dict(effects_summary)
                
        return summary_effects_df            

def run(args):
    # Create variants object
    variants_obj = Variants(args.effectsdir, args.annotation, args.region, args.outdir)
    
    # Create CDS object
    cds_obj = CDS(args.annotation, args.region)
    
    # Data frame for all CDS features
    cds_df = cds_obj.get_cds_features()
    
    # Merge all predicted variant effects
    effects_merged_list = variants_obj.merge_variant_effects()
    
    # Summary of overall predicted effects
    overall_effects = variants_obj.get_effects_overall(effects_merged_list)
    
    # Predicted effects by gene Data Frame
    effects_by_gene = variants_obj.get_effects_by_gene(effects_merged_list, cds_df)
    
    save_to_csv(effects_by_gene, "effects_summary_by_gene_", variants_obj.outdir)