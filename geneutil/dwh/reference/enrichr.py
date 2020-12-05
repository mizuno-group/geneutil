# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 13:05:34 2019

Data warehouse control class

Organizer o-- DWHControl o-- StoredRef <|-- Enrichr

@author: tadahaya
"""
import pandas as pd
import numpy as np
import os
import csv

class Enrichr():
    """
    Reference gene set data from Enrichr by Ma'ayan Laboratory
    https://amp.pharm.mssm.edu/Enrichr/
    
    """
    def __init__(self):
        self.ref = dict()
        self.__base = os.path.dirname(__file__) # ~\geneutil\dwh\reference
        self.__base_en = self.__base + "\\enrichr"
        self.__library = ""
        self.__species = ""
        self.__state = {"database":"enrichr","name":""}
        h_ready = self.__get_filenames(self.__base_en + "\\human")
        m_ready = self.__get_filenames(self.__base_en + "\\mouse")
        self.__available = {"human":h_ready,"mouse":m_ready}
        print("Reference database: Enrichr")
        print("--- all libraries currently available ---")
        print("- human")
        for l in h_ready:
            print("    ",l)
        print("- mouse")
        for l in m_ready:
            print("    ",l)
        print("-----------------------------------------")


    def get_(self):
        return self.ref


    def get_species(self):
        return self.__species


    def load_(self,library="GO_Biological_Process_2018"):
        """
        load a data set from Enrichr library
        
        Parameters
        ----------
        library: str
            indicate the name of library of interest
                    
        """
        if len(library)==0:
            raise ValueError("!! Give library:str !!")
        print("library='{}'".format(library))
        self.__library = library
        if library in self.__available["human"]:
            self.prep_(library=library,species="human")
        elif library in self.__available["mouse"]:
            self.prep_(library=library,species="mouse")
        else:
            raise ValueError("!! No indicated file. check enrichr directory !!")


    def prep_(self,library="",species=""):
        """
        prepare a reference data set by converting txt
        
        Parameters
        ----------
        library: str
            indicate the name of library of interest

        dic: SynoDict
            SynoDict object for conversion        
        
        """
        if len(species) > 0:
            self.__species = species
        else:
            raise ValueError("!! No species: indicate species !!")
        if len(library) > 0:
            self.__library = library
        url = self.__base + "\\enrichr\\{0}\\{1}.txt".format(self.__species,self.__library)
        with open(url,encoding="utf_8") as f:
            reader = csv.reader(f,delimiter='\t')
            data = [row for row in reader]
        terms = []
        ap = terms.append
        members = []
        ap2 = members.append
        for v in data:
            ap(v[0])
            del v[:2]
            temp = set(v) - {""}
            temp = {x.lower() for x in temp}
            ap2(temp)
        self.ref = dict(zip(terms,members))
        self.__state["name"] = self.__library
        # 200930 Human_Gene_Atlas data is wrong and should be treated by hard coding
        if library=="Human_Gene_Atlas":
            for k,v in self.ref.items():
                temp = {w.split(",")[0] for w in v}
                self.ref[k] = temp
                

    def __get_filenames(self,url):
        """ get filenames """
        n_all = os.listdir(url)
        return [v.replace(".txt","") for v in n_all if ".txt" in v]


    def get_state(self):
        return self.__state