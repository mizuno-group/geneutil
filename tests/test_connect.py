# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 13:05:34 2019

@author: tadahaya
"""
import unittest
import pandas as pd
import os
import sys
import math

from enan.connect import Connect

BASEPATH = os.path.dirname(os.path.abspath(__file__))

class SampleTest(unittest.TestCase):
    CLS_VAL = 'none'

    # called when test class initialization
    @classmethod
    def setUpClass(cls):
        if sys.flags.debug:
            print('> setUpClass method is called.')
        cls.CLS_VAL = '> setUpClass : initialized!'
        if sys.flags.debug:
            print(cls.CLS_VAL)

    # called when test class end
    @classmethod
    def tearDownClass(cls):
        if sys.flags.debug:
            print('> tearDownClass method is called.')
        cls.CLS_VAL = '> tearDownClass : released!'
        if sys.flags.debug:
            print(cls.CLS_VAL)

    # called when a test method runs
    def setUp(self):
        if sys.flags.debug:
            print(os.linesep + '> setUp method is called.')
        self.smpl = Connect()

    # called when a test method ends
    def tearDown(self):
        if sys.flags.debug:
            print(os.linesep + '> tearDown method is called.')

    def _df_checker(self,df):
        if type(df)!=pd.core.frame.DataFrame:
            return False
        elif df.shape[0]==0:
            return False
        else:
            head = df.head(1)
            judge = math.isnan(head.iat[0,0])
            return not judge

    def _sr_checker(self,sr):
        if type(sr)!=pd.core.series.Series:
            return False
        if sr.shape[0]==0:
            return False
        else:
            head = sr.head(1)
            judge = math.isnan(head.iat[0])
            return not judge

    def test_calc(self):
        ### preparation
        df = pd.read_csv(BASEPATH.replace("tests","enan\\_test\\") + "grm.csv",index_col=0)
        df_ref = pd.read_csv(BASEPATH.replace("tests","enan\\_test\\") + "large_grm.csv",index_col=0)
        self.smpl.gene(ref="enrichr",species="human")
        prep = self.smpl.call_PreProcessor()
        ref = prep.vector2set(df_ref,fold=4,two_sided=True)
        self.smpl.set_ref(ref)
        self.smpl.set_obj(df)

        ### test
        self.assertTrue(self._df_checker(self.smpl.calc()))