#!/usr/bin/env python
# -*- coding: utf-8 -*-


import pandas as pd


def get_swiss_index(download=False):
    swiss = pd.read_csv("databases/INDEX", comment = "#", sep = "\t")
    unique_prot_list = swiss["UniProtKB_ac"].unique()
    if download:
        with open("swiss.csv", "w") as f:
            for i in range(0, len(unique_prot_list)):
                get_unique = swiss[swiss["UniProtKB_ac"] == unique_prot_list[i]]
                try:
                    max_qmean = max(get_unique["qmeandisco_global"])
                    download_link = get_unique[get_unique["qmeandisco_global"]== max_qmean]["url"].tolist()[0]
                    print(download_link, file=f)
                except:
                    download_link = get_unique["url"].tolist()[0]
                    print(download_link,file=f)
        f.close()
    return unique_prot_list




def request(file_database):
    import requests
    import os
    prot_name = get_swiss_index()
    os.mkdir(file_database)
    database = pd.read_csv(str(file_database+".csv"), header = None)
    for i in range(0, len(database)):
        url = database[0][i]
        response = requests.get(url)
        open(str(file_database+"/%s.pdb")%prot_name[i], "wb").write(response.content)





#get_swiss_index()
request("swiss")






