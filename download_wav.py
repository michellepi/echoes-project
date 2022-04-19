# Imports
from urllib import request
from urllib.error import URLError, HTTPError
import os
import logging
import pandas as pd
import numpy as np
from tqdm import tqdm
logging.basicConfig(filename='download_wav.log', filemode='w', format='%(levelname)s:%(message)s', level=logging.DEBUG)

# URLS - This is the list of URLs with access tokens
# Use the appropriate path to the excel file containing the links
URLS = pd.read_excel(r'dl_urls.xls')
URLS = URLS.values
URLS = np.unique(URLS)
# directory - The name of the folder you want to download the files to
# No need to create a new folder, it will create one for you if it doesn't exist yet
# It will go in the same directory this script is stored in
directory = 'echoes-patient-data' 
if not os.path.exists(directory):
    os.makedirs(directory)

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
        # Updating progress bar text
        pbar.set_description("Processing %s" % name)

        full_path = directory + "/" + name
        # Download the file if it is not a duplicate
        if not os.path.exists(full_path):
            try:
                with open(full_path, "wb") as file:
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

print("Successfuly downloaded", i, "out of", len(URLS), "files,", "check download_wav.log for more info")
