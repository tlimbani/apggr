import datetime
from pandas import DataFrame
import re

def detect_filetype(annotation_file_path):
    """
    Gets file and detects whether it is in embl or genbank format

    Parameters
    ----------
    file : File
        Gets an annotation file and detects if embl or genbank
    """
    with open(annotation_file_path, 'r') as file:
        lines = [file.readline().strip() for _ in range(5)]

    # GenBank annotation
    if lines[0].startswith("LOCUS"):
        return "genbank"
    
    # EMBL annotation
    elif lines[0].startswith("ID"):
        return "embl"
    
    return "unknown"
    

def get_datenow(date_format):
    """
    Get data now in format "ddmmYYYY"
    """
    # Get today's date
    date_today = datetime.datetime.now()

    # Format today's date 
    fmt_date_today = date_today.strftime(date_format)
 
    return fmt_date_today 

def get_directory(path_str):
    pass

def save_to_csv(data_frame: DataFrame, file_name, file_path):
    """
    Save data frame as a CSV file to specified file path
    
    Parameters
    ----------
    data_frame : Pandas.DataFrame
        A Pandas data frame
    file_path : String
        File path to save CSV
    """
    file = file_path + "/" + file_name + get_datenow("%d%m%Y") + "_.csv"
    
    data_frame.to_csv(file, index=True)

def strip_non_alphanumeric(string):
    """
    Strip all non-alphanumeric characters from beginning and end of string
    
    Parameters
    ----------
    string : String
        A regular string
    """
    return re.sub('^[^A-Za-z0-9]+|[^A-Za-z0-9]+$', '', string)

def valid_region(genomic_region):
    """
    Validates genomic region to be in the correct format
    From:To
    
    Parameters
    ----------
    region : String
        Genomic region
    """
    pattern = r"^[0-9]+:[0-9]+$"
    
    return bool(re.match(pattern, genomic_region))