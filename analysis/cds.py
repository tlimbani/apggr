from Bio import SeqIO
from pandas import DataFrame
from utils import detect_filetype, save_to_csv, strip_non_alphanumeric, valid_region

class CDS:
    """
    Uses The DDBJ/ENA/GenBank Feature Table definition version 11.3 October 2024
    """
   
    def __init__(self, annotation, region):
        """
        Parameters
        ----------
        annotation : str
            Genome annotation file path (genbank, embl)
        region : str
            Genomic region to extract CDS features
        """
        self.annotation = annotation
        self.region = region
        
    def get_cds_features(self):
        """
        Get all CDS features from genomic region
        """
        
        # get annotation file type
        annotation_file_type = detect_filetype(self.annotation)
            
        # Get annotation
        record = SeqIO.read(self.annotation, annotation_file_type)
                  
        # Get start and end for genomic region
        start, end = map(int, self.region.split(":"))
       
        # Subset annotation to region
        sub_record = record[start:end]
                
        # CDS features list
        cds_collection = []
        
        for feature in sub_record.features:
            
            if feature.type == 'CDS':
                
                # UniProt ID
                protein_id = strip_non_alphanumeric(str(feature.qualifiers.get("protein_id")))
                
                # CDS description / gene
                gene_product = strip_non_alphanumeric(str(feature.qualifiers.get("product")))
                                    
                # CDS feature record list
                feature_record = []
                                    
                # Add UniProt ID
                feature_record.append(protein_id)
                
                # Add CDS description / gene
                feature_record.append(gene_product)
                
                # Length of CDS
                feature_record.append(len(feature.location))
                
                # CDS location start, adding start offset
                feature_record.append(int(feature.location.start + start))
                
                # CDS location end, adding start offset
                feature_record.append(int(feature.location.end + start))
              
                # Strand where CDS is located
                feature_record.append(feature.location.strand)
                
                # Add CDS feature record to CDS collection list
                cds_collection.append(feature_record)
        
        # Coding DNA sequence (CDS) DataFrame
        cds_collection_df = DataFrame(cds_collection, columns = ['protein_id', 'product_description', 'cds_length', 'cds_start', 'cds_end', 'strand'])
                
        return cds_collection_df

def run(args): 
    """
    Get CDS from passed arguments
    Parameters
    ----------
    args : Namespace
        Passed arguments for CDS command
    """
    
    # New CDS object
    cds_obj = CDS(args.annotation, args.region)
    
    file_type = detect_filetype(args.annotation)
    
    validated_region = valid_region(args.region)
    
    # Extract features if passed annotation file is genbank or embl format and the region is valid and save to CSV file
    if file_type != "unknown" and validated_region:
    
        # Get all CDS features into a DataFrame
        cds_collection = cds_obj.get_cds_features()
        
        file_name = "CDS_"
        
        # Save CDS features to CSV file
        save_to_csv(cds_collection, file_name, args.output)
    else: 
        raise Exception("error: wrong arguments passed for CDS analysis")