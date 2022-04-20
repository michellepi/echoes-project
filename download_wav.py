# Imports
from urllib import request
from urllib.error import URLError, HTTPError
import os
import logging
import pandas as pd
import numpy as np
from tqdm import tqdm
logging.basicConfig(filename='download_wav.log', filemode='w', format='%(levelname)s:%(message)s', level=logging.DEBUG)

# Reading the required excel files, change to appropriate file path
# META - This is the file containing the metadata of the recordings
# URLS - This is the file of URLs with access tokens
META_DF = pd.read_excel(r'../ref/echoes-meta.xlsx')
URLS_DF = pd.read_excel(r'../ref/dl_urls.xls')

# Getting the URLS without any duplicates
URLS = URLS_DF.values
URLS = np.unique(URLS)
# Getting the filenames
patient_files = META_DF['FILENAME'].values

# directory - The name of the folder you want to download the files to
# No need to create a new folder, it will create one for you if it doesn't exist yet
# It will go in the same directory this script is stored in
root = '../echoes-patient-data' 


# Here we are going through each url 
# (Unroll the loop if number of files getting too large and takes too long)
i = 0
pbar = tqdm(URLS)
for url in pbar:
    try:
        # Sending a GET request to the url
        response = request.urlopen(request.Request(url))
        # Getting the file name from the response 
        name = response.info().get_filename()
        name_only = name[:-4]
        # Updating progress bar text
        pbar.set_description("Processing %s" % name)
        
        # Check if the URL file is the file we want
        if name_only in patient_files:
            file_meta = META_DF.loc[META_DF['FILENAME'] == name_only]
            patient_num = file_meta['PATIENT'].values[0]
            position = file_meta['POS_CODE'].values[0]
            directory = root + "/Patient " + str(patient_num) + "/" + str(position)
            file_path =  directory + "/" + name

            if not os.path.exists(file_path):
                if not os.path.exists(directory):
                    os.makedirs(directory)
                try:
                    with open(file_path, "wb") as file:
                        file.write(response.read())
                        logging.info("Downloaded " + name)
                        i += 1
                except OSError as e:
                    logging.error("Failed to download" + name)
                    logging.error(f"{type(e)}: {e}")
            else:
                logging.warning(name + " already exists in directiry, skipped download")
    except HTTPError as e:
        logging.error("The server couldn't fulfill the request.")
        logging.error("Error code: " + e.code)
    except URLError as e:
        logging.error("We failed to reach a server.")
        logging.error("Reason: " + e.reason)

print("Successfuly downloaded", i, "out of", len(patient_files), " required files,", "check download_wav.log for more info")
