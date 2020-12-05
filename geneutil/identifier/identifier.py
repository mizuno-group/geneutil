# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 13:05:34 2019

Converter

@author: tadahaya
"""
import pandas as pd
import numpy as np
from  itertools import chain

from .._utils.converter import SynoDict,Integrator

class Identifier():
    def __init__(self):
        self.integrator = Integrator()
        
    def set_dic(self,dic):
        """
        set SynoDict for name identification

        """
        self.integrator.load_ref(dic)

    def generate_dic(self,keys=[],values=[],synonyms=[]):
        """ generate SynoDict for identification """
        return self.integrator.make_ref(keys,values,synonyms)

    def get_keys(self):
        """ return keys of decoders """
        return self.integrator.keys

    def extend_encoder(self,dic:dict,name="EXTENDED"):
        """
        extend encoder by adding dic

        Parameters
        ----------
        name: str
            indicates the key for the data
        
        """
        self.integrator.extend_encoder(dic=dic,name=name)

    def merge_decoder(self,keys:list=[],new_name:str="MERGED"):
        """
        merge decoders with the indicated keys

        Parameters
        ----------
        keys: list
            indicates the keys to be concatenated

        new_name: str
            indicates the key for the concatenated data
        
        """
        self.integrator.merge_decoder(keys=keys,new_name=new_name)


    def register_df(self,data,name:str=""):
        """
        registration of dataframe for integration
        index is employed for integration

        Parameters
        ----------
        name: str
            indicates the key for calling the decoder of the input
        
        """
        lst = list(data.index)
        self.integrator.register(keys=lst,name=name)

    def encode_df(self,data):
        """ convert indices of df from keys into values """
        idx = list(data.index)
        new_idx = self.integrator.enc_list(idx)
        temp = data.copy()
        temp.index = new_idx
        return temp

    def decode_df(self,data,key):
        """ convert indices of df from values into keys """
        idx = list(data.index)
        new_idx = self.integrator.dec_list(idx,key)
        temp = data.copy()
        temp.index = new_idx
        return temp


    def register_set(self,data,name:str=""):
        """
        registration of set for integration

        Parameters
        ----------
        name: str
            indicates the key for calling the decoder of the input
        
        """
        lst = list(data)
        self.integrator.register(keys=lst,name=name)

    def encode_set(self,data):
        """ convert set from keys into values """
        return self.integrator.enc_set(data)

    def decode_set(self,data,key):
        """ convert set from values into keys """
        return self.integrator.dec_set(data,key)


    def register_list(self,data,name:str=""):
        """
        registration of list for integration

        Parameters
        ----------
        name: str
            indicates the key for calling the decoder of the input
        
        """
        lst = list(data)
        self.integrator.register(keys=lst,name=name)

    def encode_list(self,data):
        """ convert list from keys into values """
        return self.integrator.enc_list(data)

    def decode_list(self,data,key):
        """ convert list from values into keys """
        return self.integrator.dec_list(data,key)


    def register_dic_set(self,data,name:str=""):
        """
        registration of dictionary of sets for integration
        {"xxxx":{"aa","bb"},"yyyy":{"cc","dd"},....}

        Parameters
        ----------
        name: str
            indicates the key for calling the decoder of the input
        
        """
        values = set(chain.from_iterable(data.values()))
        lst = list(values)
        self.integrator.register(keys=lst,name=name)

    def encode_dic_set(self,data):
        """
        convert dictionary of sets from keys into values
        {"xxxx":{"aa","bb"},"yyyy":{"cc","dd"},....}

        """
        temp = data.copy()
        for k,v in temp.items():
            temp[k] = self.integrator.enc_set(v)
        return temp

    def decode_dic_set(self,data,key):
        """
        convert dictionary of sets from values into keys
        {"xxxx":{"aa","bb"},"yyyy":{"cc","dd"},....}

        """
        temp = data.copy()
        for k,v in temp.items():
            temp[k] = self.integrator.dec_set(v,key)
        return temp


    def register_dic_2set(self,data,name:str=""):
        """
        registration of dictionary of tuples of 2 sets for integration
        {"xxxx":{{"aa","bb"},{"cc","dd"}},"yyyy":{{"ee","ff"},{"gg","hh"}},....}

        Parameters
        ----------
        name: str
            indicates the key for calling the decoder of the input
        
        """
        temp = [v[0]|v[1] for v in data.values()]
        values = set(chain.from_iterable(temp))
        lst = list(values)
        self.integrator.register(keys=lst,name=name)

    def encode_dic_2set(self,data):
        """
        convert dictionary of tuples of 2 sets from keys into values
        {"xxxx":{{"aa","bb"},{"cc","dd"}},"yyyy":{{"ee","ff"},{"gg","hh"}},....}

        """
        temp = data.copy()
        for k,v in temp.items():
            val1 = self.integrator.enc_set(v[0])
            val2 = self.integrator.enc_set(v[1])
            temp[k] = (val1,val2)
        return temp

    def decode_dic_2set(self,data,key):
        """
        convert dictionary of tuples of 2 sets from values into keys
        {"xxxx":{{"aa","bb"},{"cc","dd"}},"yyyy":{{"ee","ff"},{"gg","hh"}},....}

        """
        temp = data.copy()
        for k,v in temp.items():
            val1 = self.integrator.dec_set(v[0],key)
            val2 = self.integrator.dec_set(v[1],key)
            temp[k] = (val1,val2)
        return temp

