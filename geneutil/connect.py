# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 13:05:34 2019

organizer of geneutil package

@author: tadahaya
"""
import pandas as pd
import numpy as np
from itertools import chain

from .enan.connect import Connect
from .dwh.dwh_control import DWHControl
from .identifier.identifier import Identifier

class Analysis():
    def __init__(self,species:str="human",dic:str="biomart"):
        self.__identification = False
        self.__analyzer = Connect()
        self.__identifier = Identifier()
        self.whole = set()
        self.whole_key = set()
        self.dic = None
        self.ref = dict()
        self.obj = set()
        self.res = pd.DataFrame()
        if dic not in ["biomart"]:
            raise KeyError("!! Wrong key for dic: Choose biomart !!")
        self.__dwh = DWHControl(dic,ref=None)
        try:
            self.set_species(species=species)
        except KeyError("!! Wrong key for species: Choose human, mouse, or rat !!"):
            pass


    ### main calculation ###
    def set_species(self,species:str="human"):
        """
        specify species to be analyzed
        based on HGNC(human), MGI(mouse), or RGD(rat) data
        
        """
        self._load_dic(species=species)
        temp = dict(zip(self.dic.keys,self.dic.values))
        self.whole_key = set(self.dic.keys)
        self.define_whole(set(temp.values()),conversion=False)
        self.set_identifier(self.dic)
        self.__identifier.extend_encoder(dic=temp,name="whole")

    def prep_ref(self,data,fold:float=3.0,
                   nmin:int=None,nmax:int=None,**kwargs):
        """
        convert dataframe to the outlier set for reference data

        Parameters
        ----------
        data: dataframe
            feature x sample dataframe
        
        fold: float
            indicates fold change determining the outliers
        
        nmin,nmax: int
            indicate the minimum/maximum number of each set 
        
        """
        temp = data.copy()
        idx = list(temp.index) # care for upper case
        try:
            idx = [v.lower() for v in idx]
            temp.index = idx
        except AttributeError:
            pass
        return self.__analyzer.vector2set(data=temp,fold=fold,nmin=nmin,nmax=nmax,**kwargs)

    def load_ref(self,ref:dict,nmin=10): # cf. MsigDB is given as ID
        """
        load a reference data
        
        Parameters
        ----------
        data: dict of up-/down-tags
            keys: tag name
            values: tuple of up-/down-gene set
            {tag_name:(up-tag set,down-tag set)}        
        
        nmin: int
            indicates the number of features necessary for each set

        """
        temp0 = ref.copy()
        temp = [v[0]|v[1] for v in temp0.values()]
        whole_ref = set(chain.from_iterable(temp))
        subtract = whole_ref - self.__identifier.decode_set(self.whole,key="whole")
        if len(subtract) > 0:
            self.__identifier.register_set(subtract,name="temp")
            self.__identifier.merge_decoder(keys=["temp","whole"],new_name="whole")
        self.ref = self.__identifier.encode_dic_2set(temp0)
        self._fit(self.ref,nmin=nmin)

    def define_whole(self,whole:set,nmin:int=10,conversion:bool=True):
        """
        define whole features of enrichment analysis

        Parameters
        ----------
        whole: set
            indicates whole features

        nmin: int
            indicates the minimum number of each set

        conversion: boolean
            whether conversion to ID is done or not

        """
        if conversion:
            whole = self.__identifier.encode_set(whole)
        self.set_whole(whole)
        self.__analyzer.set_whole(whole)
        if len(self.ref)!=0:
            self.__analyzer.fit(self.ref,keep_whole=True,nmin=nmin) # adjust ref to whole

    def calc(self,data):
        """
        conduct connectivity analysis
        
        Parameters
        ----------
        data: dataframe
            feature x sample dataframe
        
        Returns res
        -------
        res: df
            gene set enrichment score
        
        """
        idx = list(data.index)
        idx = {x.lower() for x in idx}
        data2 = data.copy()
        data2.index = idx
        idx2 = idx & self.whole_key
        subtract = idx - idx2
        if len(subtract) > 10:
            print(list(subtract)[:10],"etc. were not found ({} genes were analyzed)".format(len(idx2)))
            data2 = data2.loc[idx2,:]
        elif len(subtract) > 0:
            print(subtract," were not found ({} genes were analyzed)".format(len(idx2)))
            data2 = data2.loc[idx2,:]
        self.obj = data2
        temp = self.__identifier.encode_df(data2)
        self.res = self.__analyzer.calc(data=temp)
        return self.res


    ### setter & getter ###
    def set_dic(self,dic):
        """ set SynoDict """
        self.dic = dic

    def set_whole(self,whole:set):
        """ set whole features """
        self.whole = whole

    def set_ref(self,ref:dict):
        """ set a reference for enrichment analysis """
        self.ref = ref

    def get_dic(self):
        return self.dic

    def get_whole(self):
        return self.whole

    def get_ref(self):
        return self.ref

    def get_obj(self):
        return self.obj


    ### stored data handling ###
    def _load_dic(self,species:str="human"):
        """ load a stored dictionary """
        self.__dwh.load_dict(species=species)
        self.dic = self.__dwh.get_dict()


    ### analysis setting ###
    def _fit(self,data:dict,nmin=10):
        """
        set a reference data instance to analyzer
        
        Parameters
        ----------
        data: dict
            a dictionary of sets like {"XXXX":{"aa","bb"},"YYYY":{"cc","dd","ee"},...}

        nmin: int
            indicates the minimum number of each set

        """
        self.__analyzer.fit(data,keep_whole=True,nmin=nmin) # adjust to whole by default


    ### visualization ###
    def set_res(self,res):
        """ load result data for visualization """
        self.__analyzer.set_res(res)

    def plot(self,sample_name:str=None,highlight:list=[],ylabel:str="connectivity score",**kwargs):
        """
        visualize a result of connectivity score

        Parameters
        ----------
        sample_name: str
            indicate the sample name to be visualized

        highlight: list
            indicate the plots to be highlightened

        fileout: str
            indicate the path for the output image

        dpi: int
            indicate dpi of the output image
            
        ylabel: str
            indicate the name of y axis

        title: str
            indicate the title of the plot

        color: str
            indicate the color of the bars

        fontsize: float
            indicate the fontsize in the plot

        size: float
            indicate the size of the plot

        figsize: tuple
            indicate the size of the plot

        """
        self.__analyzer.plot(highlight=highlight,sample_name=sample_name,ylabel=ylabel,**kwargs)


    ### other ###
    def check_ref(self,keyword:str):
        """ check contents of reference data """
        if len(self.ref)==0:
            raise ValueError("!! load_ref() before this process !!")
        try:
            temp = self.ref[keyword]
            temp0 = self.__identifier.decode_set(temp[0],key="whole")
            temp1 = self.__identifier.decode_set(temp[1],key="whole")
            print("{0}: {1}".format(keyword,temp0|temp1))
            return temp
        except KeyError:
            print("!! Wrong keyword !!")
            hit = {v for v in self.ref.keys() if keyword in v}
            print("perhaps: {}".format(hit))
            return tuple()

    def set_identifier(self,dic):
        """ set SynoDict for identification """
        self.__identifier.set_dic(dic)