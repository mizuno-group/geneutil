# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 13:05:34 2019

PreProcessor class

Analyzer o-- PreProcess

@author: tadahaya
"""
import pandas as pd
import numpy as np

from .._utils import statistics as st
from .._utils.converter import SynoDict,FixedDict


class PreProcessor():
    def __init__(self):
        self.__dic = DictControl()
        self.__mod = DataModifier()

    ### data modification
    def vector2set(self,data,**kwargs):
        """
        convert dataframe into dict of tags

        Parameters
        ----------
        mtx: dataframe
            feature x sample matrix

        fold: float
            determine threshold of outliers

        two_sided: boolean
            determine whether up/down is kept

        method: str
            "std" or "iqr"

        nmin: int
            indicate the minimum number for set

        nmax: int
            indicate the maximum number for set

        """
        return self.__mod.vec2set(data,**kwargs)


    ### dictionary preparation
    def generate_dict(self,names,ids,synonyms,processing=False):
        """
        generate dict considering synonyms
        !! takes a long time !!

        Parameters
        ----------
        names: list
            representative names for decoder
            
        synonyms: list
            a list of synonym sets [{},{},...]
            !! nan should be replaced with "" before construction !!
            
        ids: list
            list of ID

        processing: boolean
            whether converter generation from 0 or lists

        """
        self.__dic.generate_dict(names,ids,synonyms,processing=processing)
        return self.__dic.get_dict()

    def call_DictControl(self):
        """ call DictControl object for detailed handling """
        return self.__dic


class DataModifier():
    def __init__(self):
        pass

    def __vec2tpl(self,mtx,fold=2.0,method="iqr",nmin=None,nmax=None):
        """ convert dataframe into tuple of tags """
        sample_name = list(mtx.columns)
        n_sample = len(sample_name)
        n_feature = len(mtx.index)
        if nmin is None:
            nmin = 15
        if nmax is None:
            nmax = int(0.01*n_feature)
        if method=="std":
            upper,lower = st.outlier_std(mtx=mtx,fold=fold,axis=0)
        else:
            upper,lower = st.outlier_iqr(mtx=mtx,fold=fold,axis=0)
        res = []
        ap = res.append
        for i in range(n_sample):
            temp = mtx.iloc[:,i].sort_values(ascending=False)
            up_val = upper[i]
            low_val = lower[i]
            temp_l = list(temp.index)
            upper_tag = set(temp[temp > up_val].index)
            lower_tag = set(temp[temp < low_val].index)
            n_up = len(upper_tag)
            n_low = len(lower_tag)
            if n_up > nmax:
                upper_tag = set(temp_l[:nmax])
            elif n_up < nmin:
                upper_tag = set(temp_l[:nmin])
            if n_low > nmax:
                lower_tag = set(temp_l[-nmax:])
            elif n_low < nmin:
                lower_tag = set(temp_l[-nmin:])
            ap((upper_tag,lower_tag))
        return res


    def vec2set(self,mtx,fold=2.0,two_sided=True,method="iqr",nmin=None,nmax=None):
        """
        convert dataframe into dict of tags

        Parameters
        ----------
        mtx: dataframe
            feature x sample matrix

        fold: float
            determine threshold of outliers

        two_sided: boolean
            determine whether up/down is kept

        method: str
            "std" or "iqr"

        """
        temp = self.__vec2tpl(mtx,fold,method,nmin,nmax)
        if two_sided:
            dic = dict(zip(list(mtx.columns),temp))
        else:
            temp = [v[0] | v[1] for v in temp]
            dic = dict(zip(list(mtx.columns),temp))
        return dic


class DictControl():
    def __init__(self,names=[],ids=[],synonyms=[]):
        self.__dic = SynoDict(names,ids,synonyms,processing=False)

    def get_dict(self):
        return self.__dic

    def to_pickle(self,url):
        """ to save Synonym Dict """
        self.__dic.to_pickle(url)

    def read_pickle(self,url):
        """ to load Synonym Dict """
        self.__dic.read_pickle(url)

    def generate_dict(self,names,ids,synonyms,processing=False):
        """
        generate dict considering synonyms
        !! takes a long time !!

        Parameters
        ----------
        names: list
            representative names for decoder
            
        synonyms: list
            a list of synonym sets [{},{},...]
            !! nan should be replaced with "" before construction !!
            
        ids: list
            list of ID

        processing: boolean
            whether converter generation from 0 or lists

        """
        self.__dic = SynoDict(names,ids,synonyms,processing=processing)