# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 13:05:34 2019

Data warehouse control class

Organizer o-- DWHControl

@author: tadahaya
"""
import pandas as pd
import numpy as np

from .dictionary import *
from .reference import *

class DWHControl():
    def __init__(self,dic:str="",ref:str=""):
        self.__state = {"dict":None,"ref":None}
        if dic=="biomart":
            self.dic = biomart.Biomart()
        else:
            self.dic = None
        if ref=="enrichr":
            self.ref = enrichr.Enrichr()
        elif ref=="msigdb":
            self.ref = msigdb.Msigdb()
        else:
            self.ref = None
        if self.dic is not None:
            self.__state["dict"] = self.dic.get_state()
        if self.ref is not None:
            self.__state["ref"] = self.ref.get_state()

    # change dict
    def to_biomart(self):
        """ switch to BioMart """
        self.dic = biomart.Biomart()
        self.__state["dict"] = "biomart"

    # change ref
    def to_enrichr(self):
        """ switch to Enrichr """
        self.ref = enrichr.Enrichr()
        self.__state["ref"] = "enrichr"

    def to_msigdb(self):
        """ switch to MsigDB """
        self.ref = msigdb.Msigdb()
        self.__state["ref"] = "msigdb"

    # realization
    def load_dict(self,**kwargs): # delegation
        """ load a gene dictionary from data warehouse """
        self.dic.load_(**kwargs)
        self.__state["dict"] = self.dic.get_state()
        
    def get_dict(self): # delegation
        return self.dic.get_()

    def prep_dict(self,**kwargs): # delegation
        """ prepare a gene dictionary from data warehouse """
        self.dic.prep_(**kwargs)

    def load_ref(self,**kwargs): # delegation
        """ load a gene reference from data warehouse """
        try:
            self.ref.load_(**kwargs)
        except ValueError:
            raise ValueError("!! Sth wrong in reference data loading !!")
        self.__state["ref"] = self.ref.get_state()            

    def get_ref(self): # delegation
        return self.ref.get_()

    def prep_ref(self,**kwargs): # delegation
        """ prepare a gene reference from data warehouse """
        self.ref.prep_(**kwargs)

    def get_state(self):
        """ get states of objects """
        return self.__state
