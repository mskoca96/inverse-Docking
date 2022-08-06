#!/usr/bin/env python
# -*- coding: utf-8 -*-
import requests
import json
import pandas as pd
import urllib.parse
import urllib.request
from collections import defaultdict
import numpy as np

def prep_json_uniprot(uniprot_id):
    """ This function for the preparation of json format from uniprot list """
    uniprot_json = '{"query": {"type": "group","logical_operator": "and","nodes": [{"type": "terminal","service": "text","parameters": {"operator": "exact_match","value": "uniprot_id","attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession"}},{"type": "terminal","service": "text","parameters": {"operator": "exact_match","value": "UniProt","attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_name"}}]},"request_options": {"results_verbosity": "compact"},"return_type": "polymer_entity"}'.replace("uniprot_id",uniprot_id)
    return uniprot_json

def get_rscb(json):
    """ For the get PDBID from RSCB """
    url = "https://search.rcsb.org/rcsbsearch/v2/query?json="
    response = requests.post(url,data=json)
    #print(response.text)
    return response.text

def get_rscb_tabular(query):
    """For the get content of pdbid"""
    url = "https://data.rcsb.org/graphql?query="
    response = requests.post(url,json={'query':query})
    print(response.text)
    return response.text

def prep_organism(organism_name, structure_number):
    """ Preparation of json spesific organism name """
    organism = '{"query":{"type":"group","nodes":[{"type":"terminal","service":"full_text","parameters":{"value":"protein"}},{"type":"group","nodes":[{"type":"terminal","service":"text","parameters":{"attribute":"rcsb_entity_source_organism.ncbi_scientific_name","value":"Homo sapiens","operator":"exact_match"}},{"type":"terminal","service":"text","parameters":{"attribute":"exptl.method","value":"X-RAY DIFFRACTION","operator":"exact_match"}},{"type":"terminal","service":"text","parameters":{"attribute":"entity_poly.rcsb_entity_polymer_type","value":"Protein","operator":"exact_match"}},{"type":"group","nodes":[{"type":"terminal","service":"text","parameters":{"attribute":"rcsb_entry_info.resolution_combined","value":0.5,"operator":"less"}},{"type":"terminal","service":"text","parameters":{"attribute":"rcsb_entry_info.resolution_combined","value":{"from":0.5,"to":1,"include_lower":true,"include_upper":false},"operator":"range"}},{"type":"terminal","service":"text","parameters":{"attribute":"rcsb_entry_info.resolution_combined","value":{"from":1,"to":1.5,"include_lower":true,"include_upper":false},"operator":"range"}},{"type":"terminal","service":"text","parameters":{"attribute":"rcsb_entry_info.resolution_combined","value":{"from":1.5,"to":2,"include_lower":true,"include_upper":false},"operator":"range"}},{"type":"terminal","service":"text","parameters":{"attribute":"rcsb_entry_info.resolution_combined","value":{"from":2,"to":2.5,"include_lower":true,"include_upper":false},"operator":"range"}},{"type":"terminal","service":"text","parameters":{"attribute":"rcsb_entry_info.resolution_combined","value":{"from":2.5,"to":3,"include_lower":true,"include_upper":false},"operator":"range"}},{"type":"terminal","service":"text","parameters":{"attribute":"rcsb_entry_info.resolution_combined","value":{"from":3,"to":3.5,"include_lower":true,"include_upper":false},"operator":"range"}}],"logical_operator":"or","label":"rcsb_entry_info.resolution_combined"}],"logical_operator":"and"}],"logical_operator":"and"},"return_type":"entry","request_options":{"paginate":{"start":0,"rows":number_of_structure},"scoring_strategy":"combined","sort":[{"sort_by":"score","direction":"desc"}]}}'.replace("Homo Sapiens",organism_name).replace("number_of_structure",str(structure_number))
    return organism


def prep_tabular():
    """Preparation of json for entries knowledge"""
    df = pd.read_csv("metadata.csv")
    pdb_list = df["identifier"].tolist()
    #tabular = '{ entry(entry_id:"4HHB") { exptl { method } } }'
    sub_pdb_list = np.array_split(pdb_list,100)
    df_res = pd.DataFrame()
    for i in sub_pdb_list:
        df_list = []
        pdb_list = i.tolist()
        pdb_list = json.dumps(pdb_list)
        tabular = '{entries(entry_ids: pdb_list){rcsb_id exptl_crystal_grow { pH } rcsb_binding_affinity { comp_id value type unit } rcsb_entry_container_identifiers { entry_id } polymer_entities { rcsb_entity_source_organism { rcsb_gene_name { value } } } nonpolymer_entities { nonpolymer_comp { chem_comp { formula_weight id name } } } } }'.replace("pdb_list",str(pdb_list))
        res = get_rscb_tabular(tabular)
        res = prep_metadata(res,"rscb_tabular","tabular")
        df_list.append(res)
        df_list.append(df_res)
        print(df_list)
        df_res = pd.concat(df_list)
        print(df_res)
    return df_res


def prep_metadata(api_result,output_name,typ = "rscb"):
    """ This function prepare metadata """
    from pandas.io.json import json_normalize
    import numpy as np
    df = pd.DataFrame()
    print(api_result)
    alist = json.loads(api_result)
    if typ == "rscb":
        for i in range(0,len(alist["result_set"])):
            df = df.append(alist["result_set"][i], ignore_index = True)
        df.to_csv(output_name,index=False)
    else :
        df = json_normalize(alist)
        df = pd.concat((json_normalize(d) for d in df["data.entries"][0]), axis=0)
        col_names = ['exptl_crystal_grow', 'rcsb_binding_affinity','polymer_entities', 'nonpolymer_entities']
        col_names = ['rcsb_id', 'exptl_crystal_grow', 'rcsb_binding_affinity', 'rcsb_entry_container_identifiers', 'polymer_entities', 'nonpolymer_entities']
        df_last = pd.DataFrame()
        for i in range(0,len(df)):
            try :
                rcsb_id = json_normalize(alist["data"]["entries"][i])["rcsb_id"]
            except:
                rcsb_id = pd.DataFrame(columns=["rscb_id"])
            try:    
                exptl_crystal_grow = json_normalize(alist["data"]["entries"][i]["exptl_crystal_grow"][0])
            except:
                exptl_crystal_grow = pd.DataFrame(columns=["pH"])
            try:
                rcsb_binding_affinity = json_normalize(alist["data"]["entries"][i]["rcsb_binding_affinity"])
            except:
                rcsb_binding_affinity = pd.DataFrame(columns=["comp_id", "value", "type", "unit"])
            try:
                rcsb_entry_container_identifiers = json_normalize(alist["data"]["entries"][i]["rcsb_entry_container_identifiers"])
            except:
                rcsb_entry_container_identifiers = pd.DataFrame(columns=["entry_id"])
            try:
                polymer_entities = json_normalize(alist["data"]["entries"][i]["polymer_entities"][0]["rcsb_entity_source_organism"][0]["rcsb_gene_name"])
            except:
                polymer_entities = pd.DataFrame(columns=["value"])
            try:
                nonpolymer_entities = json_normalize(alist["data"]["entries"][i]["nonpolymer_entities"][0]["nonpolymer_comp"])
            except:
                nonpolymer_entities = pd.DataFrame(columns=["chem_comp.formula_weight","chem_comp.id","chem_comp.name"])
            print(nonpolymer_entities) 
            df2 = [rcsb_id, exptl_crystal_grow["pH"], rcsb_binding_affinity["comp_id"],rcsb_binding_affinity["value"],rcsb_binding_affinity["type"],rcsb_binding_affinity["unit"] , rcsb_entry_container_identifiers["entry_id"], polymer_entities["value"], nonpolymer_entities["chem_comp.formula_weight"],nonpolymer_entities["chem_comp.id"],nonpolymer_entities["chem_comp.name"]]
            result_1 = pd.concat(df2, join='outer', axis=1).reset_index(drop=True)
            df_last = df_last.append(result_1,ignore_index=True)
        df = df_last.reset_index(drop=True)
    return df 


def get_uniprotid_from_pdbid():
    import Bio.SwissProt as sp
    df = pd.read_csv("metadata.csv")
    pdb_list = df["identifier"].tolist()
    pdb_str = "ids="+ ",".join(pdb_list)
     



def main():
    #a = prep_organism("Homo Sapiens",60000)
    #b = get_rscb(a)
    #c = prep_tabular()
    #prep_metadata(b,"metadata.csv")
    #c.to_csv("rscb_tabular.csv")
    get_uniprotid_from_pdbid()
main()

