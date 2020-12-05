# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 13:05:34 2019

Data warehouse control class

Organizer o-- DWHControl o-- StoredRef <|-- Msigdb

@author: tadahaya
"""
import pandas as pd
import numpy as np
import os
import csv

class Msigdb():
    """
    Reference gene set data from MSigDB by Broad Institute
    https://www.gsea-msigdb.org/gsea/index.jsp
    entrez type data like "c2.all.v7.1.entrez.gmt" should be stored in ~\geneutil\dwh\reference\msigdb
    
    """
    def __init__(self):
        self.ref = dict()
        self.__base = os.path.dirname(__file__) # ~\geneutil\dwh\reference
        self.__library = []
        self.__state = {"database":"msigdb","name":""}
        self.__fileinfo = {"h":"hallmark gene sets",
                           "c1":"positional gene sets",
                           "c2":"curated gene sets",
                           "c3":"regulatory target gene sets",
                           "c4":"computational gene sets",
                           "c5":"GO gene sets",
                           "c6":"oncogenic signatures gene sets",
                           "c7":"immunologic signatures gene sets"}
        available = os.listdir(self.__base + "\\msigdb")
        available_key = [v.split(".")[0] for v in available]
        self.__available = dict(zip(available_key,available))
        print("Reference database: MsigDB")
        print("--- all libraries currently available ---")
        for k in available_key:
            print(k,"({})".format(self.__fileinfo[k]))
        print("-----------------------------------------")

    def get_(self):
        return self.ref

    def load_(self,library="c2"):
        """
        load gmt files
        both datasets and whole genes in gmt files are loaded
        
        Parameters
        ----------
        library: str
            indicate the name of library of interest
            refer to "all libraries currently available" shown in initialization
                
        """
        if library not in set(self.__available.keys()):
            raise KeyError("!! Wrong key for indicating MsigDB data !!")
        self.__library = library
        print("library='{}'".format(self.__library))
        url = self.__base + "\\msigdb\\{}".format(self.__available[self.__library])
        temp = self.__load(url)
        self.ref = dict(zip(temp[0],temp[1]))
        self.__state["name"] = self.__library

    def prep_(self):
        raise NotImplementedError

    def __load(self,url):
        """ convert gmt file to feature set """
        with open(url) as f:
            reader = csv.reader(f,delimiter="\t")
            members = []
            terms = []
            ap = members.append
            ap2 = terms.append
            for l in reader:
                ap2(l[0])
                del l[0]
                del l[0]
                ap(set(map(lambda x: int(x),l)))            
        return (terms,members)

    def get_state(self):
        return self.__state