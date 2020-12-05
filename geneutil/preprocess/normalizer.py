# -*- coding: utf-8 -*-
"""
Created on Mon Dec 10 10:39:56 2018

normalization code for data matrix

ver 200319

@author: tadahaya
"""
import pandas as pd
import numpy as np
from scipy.stats import rankdata
import time


def z_array(x):
    """
    to calculate z scores
    
    Parameters
    ----------
    x: a numpy array
        a numpy array to be analyzed
        
    """
    myu = np.mean(x,axis=0)
    sigma = np.std(x,axis=0,ddof=1)
    return (x - myu)/sigma    


def z_pop(x,axis=0,drop=True):
    """
    to calculate z scores from dataframe
    the scores employ population control
    
    Parameters
    ----------
    x: a dataframe
        a dataframe to be analyzed
    
    axis: 0 or 1
        whether z scores are calculate in column or row

    drop: boolean
        whether drop inf and nan
        
    """
    if axis==0:
        myu = np.mean(x.values,axis=0)
        sigma = np.std(x.values,axis=0,ddof=1)
    else:
        myu = np.c_[np.mean(x.values,axis=1)]
        sigma = np.c_[np.std(x.values,axis=1,ddof=1)]
    df = pd.DataFrame((x.values - myu)/sigma)    
    df.index = x.index
    df.columns = x.columns
    if drop:
        df = df.replace(np.inf,np.nan)
        df = df.replace(-np.inf,np.nan)
        df = df.dropna()
    return df


def z(x,control="",drop=True):
    """
    to calculate z scores based on control data from dataframe
    
    Parameters
    ----------
    x: a dataframe
        a dataframe to be analyzed
    
    control: string, default ""
        indicates the control column name

    drop: boolean
        whether drop inf and nan
        
    """
    if len(control) > 0:
        print("control column name: {0}".format(control))
        con = x.loc[:,x.columns.str.contains(control)]        
        n = len(con.columns)
        print("control column No.: {0}".format(n))
        if n < 3:
            print("<< CAUTION >> control columns are too few: population control was employed")
            return z_pop(x,axis=1,drop=drop)
        else:
            myu = np.c_[np.mean(con.values,axis=1)]
            sigma = np.c_[np.std(con.values,axis=1,ddof=1)]
            x = x.loc[:,~x.columns.str.contains(control)]
            df = pd.DataFrame((x.values - myu)/sigma)    
            df.index = x.index
            df.columns = x.columns
            if drop:
                df = df.replace(np.inf,np.nan)
                df = df.replace(-np.inf,np.nan)
                df = df.dropna()
            return df
    else:
        print("<< CAUTION >> no control columns: population control was employed")
        return z_pop(x,axis=1,drop=drop)


def madz_array(x):
    """
    to calculate MAD Z
    
    Parameters
    ----------
    x: a numpy array
        a numpy array to be analyzed
        
    """
    med = np.median(x,axis=0)
    mad = np.median(np.abs(x - med),axis=0)
    return (x - med)/(1.4826*mad)


def madz_pop(x,axis=0,drop=True):
    """
    to calculate MAD Z from dataframe
    the scores employ population control
    
    Parameters
    ----------
    x: a dataframe
        a dataframe to be analyzed
    
    axis: 0 or 1
        whether MAD Z are calculate in column or row

    drop: boolean
        whether drop inf and nan
        
    """
    if axis==0:
        med = np.median(x.values,axis=0)
        mad = np.median(np.abs(x.values - med),axis=0)
    else:
        med = np.c_[np.median(x.values,axis=1)]
        mad = np.c_[np.median(np.abs(x.values - med),axis=1)]
    df = pd.DataFrame((x.values - med)/(1.4826*mad))
    df.index = x.index
    df.columns = x.columns
    if drop:
        df = df.replace(np.inf,np.nan)
        df = df.replace(-np.inf,np.nan)
        df = df.dropna()
    return df


def madz(x,control="",drop=True):
    """
    to calculate MAD Z based on control data from dataframe
    
    Parameters
    ----------
    x: a dataframe
        a dataframe to be analyzed
    
    control: string, default ""
        indicates the control column name

    drop: boolean
        whether drop inf and nan
        
    """
    if len(control) > 0:
        print("control column name: {0}".format(control))
        con = x.loc[:,x.columns.str.contains(control)]        
        n = len(con.columns)
        print("control column No.: {0}".format(n))
        if n < 3:
            print("<< CAUTION >> control columns are too few: population control was employed")
            return madz_pop(x,axis=1,drop=drop)
        else:
            med = np.c_[np.median(con.values,axis=1)]
            mad = np.c_[np.median(np.abs(con.values - med),axis=1)]
            x = x.loc[:,~x.columns.str.contains(control)]
            df = pd.DataFrame((x.values - med)/(1.4826*mad))
            df.index = x.index
            df.columns = x.columns
            if drop:
                df = df.replace(np.inf,np.nan)
                df = df.replace(-np.inf,np.nan)
                df = df.dropna()
            return df
    else:
        print("<< CAUTION >> no control columns: population control was employed")
        return madz_pop(x,axis=1,drop=drop)


def robz_array(x):
    """
    to calculate robust z scores
    
    Parameters
    ----------
    x: a numpy array
        a numpy array to be analyzed
        
    """
    med = np.median(x,axis=0)
    q1,q3 = np.percentile(x,[25,75],axis=0)
    niqr = (q3-q1)*0.7413
    return (x - med)/niqr


def robz_pop(x,axis=0,drop=True):
    """
    to calculate robust z scores from dataframe
    the scores employ population control
    
    Parameters
    ----------
    x: a dataframe
        a dataframe to be analyzed
    
    axis: 0 or 1
        whether robust z scores are calculate in rows or columns

    drop: boolean
        whether drop inf and nan
        
    """
    if axis==0:
        med = np.median(x.values,axis=0)
        q1,q3 = np.percentile(x.values,[25,75],axis=0)
    else:
        med = np.c_[np.median(x.values,axis=1)]
        q1 = np.c_[np.percentile(x.values,25,axis=1)]
        q3 = np.c_[np.percentile(x.values,75,axis=1)]
    niqr = (q3-q1)*0.7413
    df = pd.DataFrame((x.values - med)/niqr)
    df.index = x.index
    df.columns = x.columns
    if drop:
        df = df.replace(np.inf,np.nan)
        df = df.replace(-np.inf,np.nan)
        df = df.dropna()
    return df


def robz(x,control="",drop=True):
    """
    to calculate robust z score based on control data from dataframe
    
    Parameters
    ----------
    x: a dataframe
        a dataframe to be analyzed
    
    control: string, default ""
        indicates the control column name

    drop: boolean
        whether drop inf and nan
        
    """
    if len(control) > 0:
        print("control column name: {0}".format(control))
        con = x.loc[:,x.columns.str.contains(control)]        
        n = len(con.columns)
        print("control column No.: {0}".format(n))
        if n < 3:
            print("<< CAUTION >> control columns are too few: population control was employed")
            return robz_pop(x,axis=1)
        else:
            con = x.loc[:,x.columns.str.contains(control)]  
            med = np.c_[np.median(con.values,axis=1)]
            q1 = np.c_[np.percentile(con.values,25,axis=1)]
            q3 = np.c_[np.percentile(con.values,75,axis=1)]
            niqr = (q3-q1)*0.7413
            x = x.loc[:,~x.columns.str.contains(control)]
            df = pd.DataFrame((x.values - med)/niqr)
            df.index = x.index
            df.columns = x.columns
            if drop:
                df = df.replace(np.inf,np.nan)
                df = df.replace(-np.inf,np.nan)
                df = df.dropna()
            return df
    else:
        print("<< CAUTION >> no control columns: population control was employed")
        return robz_pop(x,axis=1,drop=drop)


def quantile(df,method="median"):
    """
    quantile normalization of dataframe (variable x sample)
    
    Parameters
    ----------
    df: dataframe
        a dataframe subjected to QN
    
    method: str, default "median"
        determine median or mean values are employed as the template    

    """
    print("quantile normalization (QN)")
    df_c = df.copy() # deep copy
    lst_index = list(df_c.index)
    lst_col = list(df_c.columns)
    n_ind = len(lst_index)
    n_col = len(lst_col)

    ### prepare mean/median distribution
    x_sorted = np.sort(df_c.values,axis=0)[::-1]
    if method=="median":
        temp = np.median(x_sorted,axis=1)
    else:
        temp = np.mean(x_sorted,axis=1)
    temp_sorted = np.sort(temp)[::-1]

    ### prepare reference rank list
    x_rank_T = np.array([rankdata(v,method="ordinal") for v in df_c.T.values])

    ### conversion
    rank = sorted([v + 1 for v in range(n_ind)],reverse=True)
    converter = dict(list(zip(rank,temp_sorted)))
    converted = []
    converted_ap = converted.append
    print("remaining conversion count")    
    for i in range(n_col):
        transient = [converter[v] for v in list(x_rank_T[i])]
        converted_ap(transient)
        rem = n_col - i
        print("\r"+str(rem),end="")
    np_data = np.matrix(converted).T
    df2 = pd.DataFrame(np_data)
    df2.index = lst_index
    df2.columns = lst_col
    print("")
    print("use {}".format(method))
    return df2


def ts_norm(df,axis=0,ts=False):
    """
    normalization with total strength of each sample (columns)
    
    Parameters
    ----------
    df: dataframe
        a dataframe subjected to ts normalization
    
    axis: int, default 0
        determine direction of normalization, row or column
        0: normalization in column vector
        1: normalization in row vector
    
    ts: boolean, default "False"
        whether total strength array is exported or not    

    """
    if axis==0:
        norms = np.linalg.norm(df,axis=0)
        df2 = df/norms
    else:
        df = df.T
        norms = np.linalg.norm(df,axis=0)
        df2 = df/norms
        df2 = df2.T
    if ts:
        return df2,norms
    else:
        return df2


def consensus_sig(data,sep="_",position=0):
    """
    to generate consensus signature
    by linear combination with weightning Spearman correlation
    
    Parameters
    ----------
    data: a dataframe
        a dataframe to be analyzed
    
    sep: str, default "_"
        separator for sample name
        
    position: int, default 0
        indicates position of sample name such as drug    
    
    """
    print("generate consensus signature (time-consuming)")
    start = time.time()
    col = list(data.columns)
    ind = list(data.index)
    samples = list(set([v.split(sep)[position] for v in col]))
    samples.sort()
    rank = data.rank()
    df2 = pd.DataFrame()
    for i,v in enumerate(samples):
        lst_temp = [y.split(sep)[position]==v for y in col]
        temp = data.loc[:,lst_temp]
        col_temp = list(temp.columns)
        temp_rank = rank.loc[:,lst_temp]
        if len(col_temp)<=1:
            if i==0:
                df2 = data.loc[:,lst_temp]
            else:
                df2 = pd.concat([df2,data.loc[:,lst_temp]],axis=1,join="inner")            
        else:
            corr = np.corrcoef(temp_rank.T)
            corr_sum = np.sum(corr,axis=1) - 1
            corr = corr/np.c_[corr_sum]
            x = temp.values.T
            lst = []
            ap = lst.append
            for j in range(len(col_temp)):
                temp2 = x*(np.c_[corr[j]])
                temp2 = np.delete(temp2,j,axis=0)
                ap(np.sum(temp2,axis=0))
            x = pd.DataFrame(np.array(lst).T,index=ind,columns=col_temp)
            if i==0:
                df2 = x
            else:
                df2 = pd.concat([df2,x],axis=1,join="inner")
    end = time.time()
    h,mod = divmod(end - start,3600)
    m,s = divmod(mod,60)
    print("{0} hr {1} min {2} sec".format(int(h),int(m),round(s,4)))
    return df2


def undo_autoconv(data):
    """ correct wrong indices by excel """
    DEC = {}
    MARCH = {}
    OCT = {}
    SEPT = {}
    for i in range(20):
        DEC["{}-Dec".format(i)] = "DEC{}".format(i)
        DEC["{}-DEC".format(i)] = "DEC{}".format(i)
        MARCH["{}-Mar".format(i)] = "MARCH{}".format(i)
        MARCH["{}-MAR".format(i)] = "MARCH{}".format(i)
        OCT["{}-Oct".format(i)] = "OCT{}".format(i)
        OCT["{}-OCT".format(i)] = "OCT{}".format(i)
        SEPT["{}-Sep".format(i)] = "SEPT{}".format(i)    
        SEPT["{}-SEP".format(i)] = "SEPT{}".format(i)    
    converter = {}
    converter.update(DEC)
    converter.update(MARCH)
    converter.update(OCT)
    converter.update(SEPT)
    if type(data)==list:
        return [converter[v] if v in converter.keys() else v for v in lst]
    else:
        ind = list(data.index)
        new = [converter[v] if v in converter.keys() else v for v in ind]
        data.index = new
        data = data.loc[~data.index.str.contains("---"),:]
        return data