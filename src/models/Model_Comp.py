# -*- coding: utf-8 -*-
###################################################################
#
# Author: Matthew W. Boyce, PhD, boyce.matthew@epa.gov
#
# Version: 1.0 3-26-20202
#
# Description:  This script includes a number of functions that clean
#               data generated from metabolite prediction softwares:
#               OCED Toolbox, Nexus Meteor, Oasis TIMES, and Biotransformer.
#               Each function returns a data frame consisting of the DTXSID 
#               of the parent molecule and the metabolite's InChI key.
#               Of note, the data require SMILES to be included in each data set
#               and for cleaning OCED Toolbox data, a dictionary object 
#               with the parent molecules' QSAR ready InChI keys as keys 
#               and the DTXSID as the values.
#
# Notes: This script uses pandas, numpy, and rdkit and their dependecies
#
# Potential issues: None known
#
###################################################################

"""
Created on Tue Mar 17 10:30:53 2020

@author: MBOYCE
"""
import os as os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import PandasTools
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def TIMES_cleanup (file, Model_Module):
    """Cleans data genrated by Oasis TIMES and returns a dataframe witht he DTXSID of the parent compound and InChI key of each metabolite"""
    """The Model_Module argument should be a string to designate the model used for metabolism (e.g., TIMES_RatLiver S9, TIMES_RatInVivo"""
    df = pd.read_csv(file, delimiter = ",")    
    df.rename(columns={'Chem. Name':'DTXSID'},inplace = True) 
    df = df[['DTXSID', 'Smiles']]
    df = df[:-2]                                                                      
    df[Model_Module] = 1                                                            #Adds Times_(Model_Module) to designate model generating metabolite
    df['DTXSID'].replace({' ': np.NaN}, inplace = True)                            #Cleans DTXSID to list NaN in empty rows
    df['Metabolite_INCHIKEY'] = np.NaN                                             #Initialized column for metabolite InChI key
    metabList = df.Smiles[df['DTXSID'].isnull()]                                    #Establishes boolean list to desgiante indecies with metabolite smiles
    df['Metabolite_INCHIKEY'] = SMILES_to_InchiKey(metabList)                       #Converts metabolie SMILES to InChI keys
    df['DTXSID'] = df['DTXSID'].fillna(method = 'ffill')                            #Fills empty spaces with DTXSID, empty spaces are filled with the preceeding DTXSID
    df = df[df['Metabolite_INCHIKEY'].notnull()]                                    #Removes any all parent entries, whcih are represented as nulls in the metabolite INCHIKEY list
    df = df.drop_duplicates()  
    df[['Formula','[M+H]']] = SMILES_to_MW(df.Smiles)
    df['Clean_SMILES'] = clean_SMILES(df['Smiles'])
    return df[['DTXSID','Metabolite_INCHIKEY','Clean_SMILES','Formula','[M+H]', Model_Module]];

def Meteor_cleanup (file):
    """Cleans and returns  a dataframe for the imported OCED Toolbox metabolite data."""
    df = []
    df = pd.read_csv(file,                                                          #Reads Meteor data after it is merged into a single file
                     header = 0, 
                     usecols = ['SMILES','Name','Query Name', 'Parent'])
    df['Metabolite_INCHIKEY'] = np.NaN
    df['Meteor'] = 1                                                                #Added column to designate model that generated metabolite
    metabList = df.Parent.notnull()                                                 #Establishes boolean list of metabolites SMILES strings
    df['Metabolite_INCHIKEY'] = SMILES_to_InchiKey(df.SMILES[metabList])            #Converts metabolite SMILES to InChI keys
    df['Query Name'] = df['Query Name'].str.replace(' \([^()]*\)',"")               #Removes ' (Query)' from each DTXSID in column
    df = df.rename(columns={'Query Name':'DTXSID'})
    df = df[df.Metabolite_INCHIKEY.notnull()]                                       #Removes parent SMILES for the dataframe
    df = df.drop_duplicates()
    df[['Formula','[M+H]']] = SMILES_to_MW(df.SMILES)
    df['Clean_SMILES'] = clean_SMILES(df['SMILES'])
    return df[['DTXSID','Metabolite_INCHIKEY','Clean_SMILES','Formula','[M+H]','Meteor']];

def ToolBox_cleanup(file, DTXSIDdict, coding = 'UTF-8', delimiter = ','):
    """Cleans and returns  a dataframe for the exported OCED Toolbox metabolite data.""" 
    """"Input file requires that SMILES by exported as part of the .csv file. The DTXSDIdict argument should be a dictionary with the QSAR Ready InChI keys as the key and the DSTXID as teh value."""
    """If issues occur reading the file, try coding = 'UTF-16' and delimiter = '\t' """
    df = pd.read_csv(file, sep = delimiter,                                                      #Reads toolbox data as a tab-delimited filewith UTF-16 encoding
                     encoding = coding,                                                   #Using UTF-16 encoding due to errors in most recent file saves,
                     header = 0, usecols =                                                  #but have been able to use UTF-8 prior
                     ['SMILES','Metabolite'])
    df = df[:-2]                                                                            #Removes empty bottom row
    df = df[df.Metabolite.notnull()]                                                        #Establishes boolean list indicating indecies with metabolite
    df['Metabolite_INCHIKEY'] = SMILES_to_InchiKey(df['Metabolite'])                        #Converts metabolite SMILES to InChI keys
    df['Parent_INCHIKEY'] = SMILES_to_InchiKey(df['SMILES'],stereoisomer=False)             #Converts parent SMILES to QSAR Ready InChI keys (removes stereoisomer features during conversion)
    df['DTXSID'] = [DTXSIDdict.get(e) for e in df['Parent_INCHIKEY']]                       #Uses dictionary of parent molecules to extract 
    df['ToolBox'] = 1                                                                       #Generate column indicating the model source of the metabolite
    df = df.drop_duplicates()    
    df[['Formula','[M+H]']] = SMILES_to_MW(df.Metabolite)
    df['Clean_SMILES'] = clean_SMILES(df['Metabolite'])                                                   
    return df[['DTXSID','Metabolite_INCHIKEY','Clean_SMILES','Formula','[M+H]','ToolBox']];

def BioTransformer_cleanup(file, DTXSIDdict):
    """Cleans and returns a dataframe for results of BioTransformer data."""
    """Input dictionary should have the normal InChI keys (i.e., not QSAR ready) as the key and DTXSID as the value"""
    df = pd.read_csv(file, header = 0, usecols = ['InChIKey','Precursor InChIKey', 'Molecular formula','Major Isotope Mass','SMILES', 'InChI'])         #Reads metabolite InChI key and DTXSID in the file
    df = df.rename(columns = {'InChIKey':'Metabolite_INCHIKEY','Precursor InChIKey':'Parent_INCHIKEY','Molecular formula':'Formula','Major Isotope Mass':'[M+H]'})    #Renames columns
    df['DTXSID'] = [DTXSIDdict.get(e) for e in df['Parent_INCHIKEY']]
    df['DTXSID'] = df['DTXSID'].fillna(method = 'ffill')                       #Reads dictionary with Parent InChI keys are the key, and DTXSIDs as the value
    df['BioTransformer'] = 1                                                                #Generate column indicating the model source of the metabolite
    df['[M+H]'] = df['[M+H]'].apply(lambda x: x + Descriptors.ExactMolWt(Chem.MolFromSmiles('[H+]')))
    df = df.drop_duplicates()
    df['Clean_SMILES'] = clean_SMILES(df['InChI'], source = 'InChI')
    return df[['DTXSID','Metabolite_INCHIKEY','Clean_SMILES','Formula','[M+H]','BioTransformer']];

def SMILES_to_InchiKey (smile_List, stereoisomer = True):
    """Uses RDKit to convert a lsit of SMILES to a list of InChI keys"""
    """If stereochemistry is not wanted (e.g., to generate QSAR Ready InChI keys, the stereoisomer argument should be set to false"""
    molList = []                                                                            #initializes a series of lists for the mols, smiles, and inchi keys
    clean_SMILES = []                                                                 
    InchiList = []
    if stereoisomer == False:                                                               #If QSAR Ready InChI Keys are needed, a SMILES without specified stereochemistry
        molList = smile_List.apply(lambda x: Chem.MolFromSmiles(x))                         #is generated from the initial set of mols. New mols are generated from the
        clean_SMILES = molList.apply(lambda x: Chem.MolToSmiles(x, isomericSmiles=False))   #stereochemistry-free SMILES, and QSAR ready InChI keys are genreated from
        molList = clean_SMILES.apply(lambda x: Chem.MolFromSmiles(x))                       #the cleaned mols
        InchiList = molList.apply(lambda x: Chem.MolToInchi(x))
        InchiList = InchiList.apply(lambda x : Chem.InchiToInchiKey(x))
    else:
        molList = smile_List.apply(lambda x: Chem.MolFromSmiles(x))                         #Uses RDKit to convert SMILES to mols,
        InchiList = molList.apply(lambda x: Chem.MolToInchi(x))
        InchiList = InchiList.apply(lambda x : Chem.InchiToInchiKey(x))                     #then mols to InChI keys
    return InchiList;

def clean_SMILES (series_List, source = 'SMILES'):
    cleanSMILES = []
    mol_list = []
    if source == 'SMILES':
        mol_list = series_List.apply(lambda x: Chem.MolFromSmiles(x))
    elif source == 'InChI':
        mol_list = series_List.apply(lambda x: Chem.MolFromInchi(x))
    cleanSMILES = mol_list.apply(lambda x: Chem.MolToSmiles(x, isomericSmiles=False))
    return cleanSMILES;

def SMILES_to_MW (smile_List):
    molList = []
    MH_list = []
    molForm_list = []
    mass_H = Descriptors.ExactMolWt(Chem.MolFromSmiles('[H+]'))
    molList = smile_List.apply(lambda x: Chem.MolFromSmiles(x))
    MH_list = molList.apply(lambda x: Descriptors.ExactMolWt(x) + mass_H)
    molForm_list = molList.apply(lambda x: rdMolDescriptors.CalcMolFormula(x))
    data = pd.DataFrame({'Formula':molForm_list, '[M+H]':MH_list})
    return data;

#Aggregate all data single dataframe by iteratively combining dataframes
def aggregate_DFs(DF_list, arg_on = ['DTXSID','Metabolite_INCHIKEY','Formula','[M+H]'], arg_how = 'outer'):
    numDF = len(DF_list)
    if numDF < 2:
        return
    agg_DF = DF_list[0]
    comp_DF = DF_list[1:]
    for DF in range(len(comp_DF)):
        agg_DF = pd.merge(agg_DF, comp_DF[DF], on = arg_on, how = arg_how)
    agg_DF = agg_DF.replace(np.NaN, 0)
    all_smiles = agg_DF.loc[:,('Clean_SMILES','Clean_SMILES_x','Clean_SMILES_y')]
    all_smiles.replace(0,np.NaN, inplace = True)
    mode_smiles = all_smiles.mode(axis = 1, dropna = True)
    agg_DF['SMILES'] = mode_smiles[0]
    agg_DF.drop(columns = ['Clean_SMILES','Clean_SMILES_x','Clean_SMILES_y'], inplace = True)
    return agg_DF
        

#####Example of run sequenc to process the a set of data using original CopTox file and model outputs
###Set directory of files exported from each of the models
#os.chdir('C:\\Users\\MBOYCE\\Documents\\ExpoCast_39CompData\\All Data\\')
    
#####Import starting DSSToxID data and intializes dictionary for importing ToolBox data
#DSSToxList = pd.read_csv("CompToxList.csv", header = 0)
#DSSToxList = DSSToxList.rename(columns={'INCHIKEY':'Parent_INCHIKEY'})
#DSSToxList['QSAR_READY_INCHIKEY'] = SMILES_to_InchiKey(DSSToxList['QSAR_READY_SMILES'],stereoisomer = False)
#Norm_DTXSID_dict = dict(zip(DSSToxList['Parent_INCHIKEY'],DSSToxList['DTXSID']))
#QSAR_DTXSID_dict = dict(zip(DSSToxList['QSAR_READY_INCHIKEY'],DSSToxList['DTXSID']))

#####Import and clean up each DF
# toolBoxDF = ToolBox_cleanup('ToolBox_Report.csv', QSAR_DTXSID_dict)
# meteorDF = Meteor_cleanup('Meteor_Report.csv')
# bioTransformerDF = BioTransformer_cleanup('BioTransformer_Report.csv', Norm_DTXSID_dict)
# times_inVivoDF = TIMES_cleanup('TIMES_invivo.txt', 'TIMES_InVivo')
# times_inVitroDF = TIMES_cleanup('TIMES_invitro.txt', 'TIMES_InVitro')

######Combine DFs into a list, and aggregate DFs into single file
# dfList = [toolBoxDF, meteorDF, bioTransformerDF,times_inVivoDF, times_inVitroDF]
# agg_Data = aggregate_DFs(dfList)












