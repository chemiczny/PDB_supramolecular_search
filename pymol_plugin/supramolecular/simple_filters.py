#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  3 10:37:23 2018

@author: michal
"""
import pandas as pd

cationicAA = ["LYS","ARG"]
acidicAA = ["ASP","GLU"]
restAA = ["ALA", "CYS","GLY","ILE","LEU","MET","ASN","PRO","GLN","SER","THR","VAL"]
allAA = restAA + acidicAA
ringAA = ["PHE", "HIS", "TRP", "TYR"]
NU = ["A","G","T","C","U","I","DA", "DC", "DG", "DT", "DI" ]

def noAAinPiAcids(actualData):
    return actualData[ ~actualData["Pi acid Code"].isin(ringAA ) ]

def noAAinPiRes(actualData):
    return actualData[ ~actualData["Pi res code"].isin(ringAA) ]

def noAAinAnions(actualData):
    return actualData[~actualData["Anion code"].isin( allAA ) ]
   
def noNUinAnions(actualData):
    return actualData[ ~actualData["Anion code"].isin( NU) ]
   
def noNUinPiAcids(actualData):
    return actualData[~actualData["Pi acid Code"].isin(NU )  ]

def noNUinPiRes(actualData):
    return actualData[~actualData["Pi res code"].isin( NU ) ]

def noAAinHAcceptors(actualData):
    return actualData[~actualData["Acceptor code"].isin( acidicAA) ]

def noNUinHAcceptors(actualData):
    return actualData[~actualData["Acceptor code"].isin( NU ) ]

def noAAinHDonors(actualData):
    return actualData[~actualData["Donor code"].isin([ "TYR" , "PHE" , "HIS" ,"TRP", "LYS","ARG" ]) ]

def noNUinHDonors(actualData):
    return actualData[~actualData["Donor code"].isin( NU ) ]

def noAAinCations(actualData):
    return actualData[~actualData["Cation code"].isin(cationicAA)]

def simpleMerge( dataFramesToMerge , dataFrameMergeHeaders, dataFramesToExclude, dataFrameExcludeHeaders ):
    if len(dataFramesToMerge) + len(dataFramesToExclude) < 2:
        return
        
    uniqueData = []
    dataExcluded =[]
    
    actualKeys = []
    excludedKeys = []
    
    for df, headers in zip(dataFramesToMerge, dataFrameMergeHeaders):               
        if len(uniqueData) == 0:
            uniqueData = df[ headers ].drop_duplicates()
        elif len(df) > 0:
            uniqueData = pd.merge( uniqueData,  df[ headers ], on = list( set(actualKeys) & set(headers) ))
            uniqueData = uniqueData.drop_duplicates()
        actualKeys = list(set( actualKeys + headers ))
            
            
    for df, headers in zip(dataFramesToExclude, dataFrameExcludeHeaders):   
            if len(dataExcluded) == 0:
                dataExcluded = df[ headers ].drop_duplicates()
            elif len(df) > 0:
                dataExcluded = pd.merge( dataExcluded,  df[ headers ], on = list( set(actualKeys) & set(headers) ))
                dataExcluded = dataExcluded.drop_duplicates()
            excludedKeys = list(set( excludedKeys + headers ))
            
    if len(dataExcluded) > 0:
        mergingKeys = list(set(actualKeys) & set(excludedKeys) )
        subMerged = pd.merge( uniqueData,  dataExcluded , on = mergingKeys, how='left', indicator=True )
        uniqueData = subMerged[ subMerged['_merge'] == 'left_only' ]
        
    allDf = dataFramesToMerge + dataFramesToExclude
    allHeaders = dataFrameMergeHeaders + dataFrameExcludeHeaders
    
    newDf = []
    for df, headers in zip(allDf, allHeaders):        
        if len(uniqueData) == 0 :
            break
        
        mergingKeys = list(set(actualKeys) & set(headers) )
        tempDataFrame = uniqueData[ mergingKeys   ].drop_duplicates()
        newDf.append(pd.merge( df, tempDataFrame, on = mergingKeys ))
        
    return newDf
        
        
        
        