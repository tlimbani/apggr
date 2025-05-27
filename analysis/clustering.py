from fnmatch import fnmatch
from os import listdir
from os.path import join, splitext
from pandas import DataFrame
from sklearn.cluster import DBSCAN
from sklearn.decomposition import PCA
from utils import save_to_csv
from visualize import visualize
import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics import silhouette_score

class Clustering:
    """
    Creates a Cluster object that can be used to run a model for DBSCAN 
    Adds metrics to assess model performance
    """
    def __init__(self, model_params, variant_directory, output_directory):
        self.model_params = model_params
        self.variant_directory = variant_directory
        self.output_directory = output_directory
        self.performance_metrics = {}
   
    # List of lists of variant rows
    def merge_variants(self):
        variants_list = []
        
        # For each variant text file, add to the variants list (CHROM, POS, REF, ALT)
        for variant_file in listdir(self.variant_directory):
            
            # Get text files for varients
            if fnmatch(variant_file, "*.txt") and not fnmatch(variant_file, "*txt_warnings.txt"):
          
                
                # Full path
                variant_path = join(self.variant_directory, variant_file)
                
                # Get file name only without the extension
                variant_name = splitext(variant_file)[0]
                
                # Open file and read content
                with open(variant_path, 'r') as file:
                    
                    rows = file.readlines()[1:]
                  
                    # For each variant row create a list split by tab and add to the main list
                    for row in rows:
                        
                        
                        if not row.startswith("#"):
                            variant_list_row = []
                            
                            fields = row.strip().split('\t')
                            
                            # SNP position
                            snp_location=fields[1].split(':')[1]
                            
                            # conseuqnce
                            consequence=fields[6]
                            
                            # Sequence name
                            sequence_name=variant_name
                          
                            variant_list_row.append(sequence_name)
                            variant_list_row.append(snp_location)
                            variant_list_row.append(consequence)
                            
                            fields = row.strip().split('\t')
                            
                            # Add sequence ID 
                            fields.append(variant_name)
                           
                            # Add complete row for varient to master list
                            variants_list.append(variant_list_row)
                    
        return variants_list
    
    def get_snp_list(self, variants_list):
        
       
        variants_df = DataFrame(variants_list, columns = ['sequence_id', 'POS', 'consequence'])
       
        # Frequency of SNP positions
        snp_position_frequency = variants_df['POS'].value_counts()
      
        snp_list = list(snp_position_frequency.index)
    
        # Return SNP list
        return snp_list
    
    
    def get_variants_in_directory(self, variants_directory):
        # Sequence names list
        variants_list_tmp = []
        
        # Loop through all sequences in study data directory
        for sequence in listdir(variants_directory):
            
            # Get text fasta for varients
            if fnmatch(sequence, "*.txt") and not fnmatch(sequence, "*txt_warnings.txt"):
                
                # Get file name only without the extension
                sequence_name = splitext(sequence)[0]
                
                # Append sequence name to sequences list
                variants_list_tmp.append(sequence_name)
                
        # Return list of all sequences in 
        return variants_list_tmp
    

    def get_matrix(self, variants_list, snp_list):
       
        variants_df = DataFrame(variants_list, columns = ['sequence_id', 'POS', 'consequence'])
        
        snp_positions_sequencies_matrix_list = []

        study_sequence_list = self.get_variants_in_directory(self.variant_directory)
      
        # For each sequence get SNPs and lookup positions, populating the matrix
        for sequence_id in study_sequence_list:
            
            # Sequence row in matrix list
            sequence_row = []
            
            # Get a subset of variants data frame for positions in sequence
            subset_variants_df = variants_df[variants_df['sequence_id'] == sequence_id]
            
            # Look if SNP present at position in sequence
            for position in snp_list:
                
                # Match found in sequence for position
                if position in subset_variants_df['POS'].values:
                    sequence_row.append(1)
                
                # Match not found for position in sequence
                else:
                    sequence_row.append(0)
            
            # Add row for this sequence to matrix list
            snp_positions_sequencies_matrix_list.append(sequence_row)
                
        # Create a matrix data frame for SNP positions and sequences
        snp_positions_sequencies_matrix = DataFrame(snp_positions_sequencies_matrix_list, columns = snp_list, index = study_sequence_list)
        
        save_to_csv(snp_positions_sequencies_matrix, "-------csv----.", self.output_directory)
        
        return snp_positions_sequencies_matrix
    
    def run_DBSCAN_model(self, matrix_df):
      
        # Indexes for saving
        indexes = matrix_df.index
        
        matrix_df.reset_index(drop=True, inplace=True)
        
        matrix_np_array = matrix_df.values
              
        pca = PCA(n_components=2)
             
        X_pca = pca.fit_transform(matrix_np_array)
                
        epsilon, min_samples_, distance = map(lambda func, val: func(val), [float, int, str], self.model_params.split(":"))
        
        # DBSCAN object
        db = DBSCAN(eps=epsilon, min_samples=min_samples_, metric=distance)
       
        # Fit model and predict cluster labels
        labels = db.fit_predict(X_pca)
        
        # Combine PCA components and cluster labels into a DataFrame for visualization in R
        cluster_data = DataFrame(X_pca, columns=['PC1', 'PC2'])
        cluster_data['Cluster'] = labels
        cluster_data['row_index'] = indexes
        cluster_data.to_csv("/path/to/csv/file", index=False)
                
        unique_labels = set(labels)
        
        # Number of clusters
        n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
        n_noise_ = list(labels).count(-1)
        silhouette_coefficient = silhouette_score(X_pca, labels)
        
       
        # Performance metrics
        self.performance_metrics['number_of_clusters'] = n_clusters_
        self.performance_metrics['number_of_noise'] = n_noise_
        self.performance_metrics['silhoutte_score'] = silhouette_coefficient
        
        plt.figure(figsize=(8, 6))
        
        for label in unique_labels:
            cluster_points = X_pca[labels == label]
            cluster_size = len(cluster_points)  # Number of points in the cluster
            if label == -1:
                # Plot noise points in black
                plt.scatter(cluster_points[:, 0], cluster_points[:, 1], color='black', label='Outlier')
            else:
                # Plot clusters in different colors
                plt.scatter(cluster_points[:, 0], cluster_points[:, 1], label=f'Cluster {label}')
                
                # Add the number of points in each cluster as text
                # The position of the label is the mean of the cluster's points
                cluster_center = np.mean(cluster_points, axis=0)
                plt.text(cluster_center[0], cluster_center[1], f'{cluster_size}', color='red', fontsize=15, ha='center', va='center')
        
        plt.title(r'S. Typhimurium O-Antigen SNP clusters of sequences')
        plt.legend()
              
        plt.xlabel('Principal Component 1')
        plt.ylabel('Principal Component 2')
      
        plt.savefig('/path/to/plot/image', dpi=300)
        
        plt.close()
            
    def get_indexes_with_cluster_assignment(self):
        pass
    
    def get_model_performance_metrics(self):
        return self.performance_metrics
        
        
    
    def visualize_cluster(DBSCAN_model_data):
        """
        Scatter plot
        """
      
        visualize.visualize_cluster(DBSCAN_model_data)
 
def run(args):
    # Create cluster object
    clustering_obj = Clustering(args.params, args.snpdir, args.outdir)
    
    # Get list of variants
    variants_ = clustering_obj.merge_variants()
    
    # List of SNP positions present in one percent of the population
    one_percent_snp = clustering_obj.get_snp_list(variants_)
        
    df_d = DataFrame(one_percent_snp)
    
    file_name = "ome_pecent_SNPs"
    
    save_to_csv(df_d, file_name, clustering_obj.output_directory)
    
    # Create matrix from list of variants and SNP positions present in one percent of population
    matrix = clustering_obj.get_matrix(variants_, one_percent_snp)
    
    DBSCAN_model_data = clustering_obj.run_DBSCAN_model(matrix)