# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 13:05:34 2019

organizer of geneutil package

@author: tadahaya
"""
import pandas as pd
import numpy as np
from itertools import chain

from .enan.fet import FET
from .dwh.dwh_control import DWHControl
from .identifier.identifier import Identifier

class Analysis():
    def __init__(self,species:str="human",dic:str="biomart",ref:str="enrichr"):
        self.__identification = False
        self.__analyzer = FET()
        self.__identifier = Identifier()
        self.whole = set()
        self.whole_key = set()
        self.dic = None
        self.ref = dict()
        self.obj = set()
        self.res = pd.DataFrame()
        if dic not in ["biomart"]:
            raise KeyError("!! Wrong key for dic: Choose biomart !!")
        if ref not in ["enrichr","msigdb"]:
            raise KeyError("!! Wrong key for ref: Choose enrichr or msigdb !!")
        self.__dwh = DWHControl(dic,ref)
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

    def load_ref(self,library:str="",nmin=10):
        """
        load a stored reference
        
        Parameters
        ----------
        library: str
            indicate a library name such as GO_Biological_Process_2018

        nmin: int
            indicate the minimum number of the members of a gene group
        
        """
        self._load_ref(library=library)
        temp = self.ref.copy()
        state = self.__dwh.get_state()
        if state["ref"]["database"]=="enrichr": # cf. MsigDB is given as ID
            for k,v in temp.items():
                temp[k] = {w.lower() for w in v}
            whole_ref = set(chain.from_iterable(temp.values()))
            subtract = whole_ref - self.__identifier.decode_set(self.whole,key="whole")
            if len(subtract) > 0:
                self.__identifier.register_set(subtract,name="temp")
                self.__identifier.merge_decoder(keys=["temp","whole"],new_name="whole")
            self.ref = self.__identifier.encode_dic_set(temp)
        self._fit(self.ref,nmin=nmin)
        print(self.__dwh.get_state())

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

    def calc(self,data:set,**kwargs):
        """
        conduct enrichment analysis

        Parameters
        ----------            
        correction: str
            indicate method for correcting multiple tests
            depend on "statsmodels.stats.multitest.multipletests"
            
        focus: int
            export results by XX th lowest p value
        
        mode: int
            indicate the type of test
            "two-sided", "less", "greater"
            in many case of ontology analysis, "greater" is appropriate
    
        """
        data = {x.lower() for x in data}
        self.obj = data & self.whole_key
        subtract = data - self.obj
        if len(subtract) > 10:
            print(list(subtract)[:10],"etc. were not found ({} genes were analyzed)".format(len(self.obj)))
        elif len(subtract) > 0:
            print(subtract," were not found ({} genes were analyzed)".format(len(self.obj)))
        data2 = self.__identifier.encode_set(self.obj)
        temp = self.__analyzer.calc(data2,**kwargs)
        temp["overlap"] = temp["overlap"].map(lambda x: self.__identifier.decode_set(x,key="whole"))
        self.res = temp
        return self.res

    def calc_binary(self,data:set,alpha:float=0.01,correction:str="holm",**kwargs):
        """
        return whether the subjected data shows significant enrichment or not

        Parameters
        ----------            
        alpha: float
            indicate the threshold for judging the significance of the test

        correction: str
            indicate method for correcting multiple tests
            depend on "statsmodels.stats.multitest.multipletests"
            
        focus: int
            export results by XX th lowest p value
        
        mode: str
            indicate the type of test
            "two-sided", "less", "greater"
            in many case of ontology analysis, "greater" is appropriate
    
        """
        data = {x.lower() for x in data}
        self.obj = data & self.whole_key
        subtract = data - self.obj
        data2 = self.__identifier.encode_set(self.obj)
        self.res = self.__analyzer.calc(data2,correction=correction,**kwargs)
        if self.res.shape[0]==0:
            return False
        else:
            return np.min(self.res["adjusted p value"]) < alpha


    ### setter & getter ###
    def set_dic(self,dic):
        """ set SynoDict """
        self.dic = dic

    def set_whole(self,whole:set):
        """ set whole features """
        self.whole = whole

    def set_ref(self,ref:dict,conversion:bool=True,nmin:int=10):
        """
        set a reference for enrichment analysis
        
        Parameters
        ----------
        ref: dict
            a dict for reference {'key':{gene set},...}

        conversion: bool
            whether conversion into ID is necessary or not

        nmin: int
            indicate the minimum number of the members of a gene group

        """
        self.ref = ref
        temp = self.ref.copy()
        if conversion:
            for k,v in temp.items():
                temp[k] = {w.lower() for w in v}
            whole_ref = set(chain.from_iterable(temp.values()))
            subtract = whole_ref - self.__identifier.decode_set(self.whole,key="whole")
            if len(subtract) > 0:
                self.__identifier.register_set(subtract,name="temp")
                self.__identifier.merge_decoder(keys=["temp","whole"],new_name="whole")
            self.ref = self.__identifier.encode_dic_set(temp)
        self._fit(self.ref,nmin=nmin)

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

    def _load_ref(self,library:str=""):
        """
        load a stored reference
        ref converted to lower case
        
        """
        self.__dwh.load_ref(library=library)
        self.ref = self.__dwh.get_ref()


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


    def plot(self,display_num:int=5,**kwargs): # realization
        """
        visualize a result of enrichment analysis

        Parameters
        ----------
        display_num: int
            determine how many groups are visualized

        fileout: str
            indicate the path for the output image

        dpi: int
            indicate dpi of the output image
            
        thresh: float
            indicate the threshold value of adjusted p value to be colored

        xlabel,ylabel: str
            indicate the name of x and y axes

        title: str
            indicate the title of the plot

        color: str
            indicate the color of the bars

        alpha: float
            indicate transparency of the bars: (0,1)

        height: float
            indicate the height of the bars

        fontsize: float
            indicate the fontsize in the plot

        textsize: float
            indicate the fontsize of the texts in the bars

        figsize: tuple
            indicate the size of the plot

        """
        self.__analyzer.plot(display_num=display_num,**kwargs)


    ### other ###
    def check_ref(self,keyword:str):
        """ check contents of reference data """
        if len(self.ref)==0:
            raise ValueError("!! load_ref() before this process !!")
        try:
            temp = self.ref[keyword]    
            temp = self.__identifier.decode_set(temp,key="whole")
            print("{0}: {1}".format(keyword,temp))
            return temp
        except KeyError:
            print("!! Wrong keyword !!")
            hit = {v for v in self.ref.keys() if keyword in v}
            print("perhaps: {}".format(hit))
            return set()

    def set_identifier(self,dic):
        """ set SynoDict for identification """
        self.__identifier.set_dic(dic)