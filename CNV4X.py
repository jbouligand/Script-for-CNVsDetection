#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = 'Christophe Habib'
__copyright__ = 'Copyright (C) 2017 GMPH'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'christophe.habib@gmail.com'
__status__ = 'prod'

"""
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
"""

import sys,os,argparse
import scipy as sp
import numpy as np



def cov2dico(filin):
    #Format tabule : ChrX Start Stop Gene 0 strand Pos Cov
    fichier = open(filin)
    lignes  = fichier.readlines()
    fichier.close()
    dico = {}
    for ligne in lignes:
        col = ligne.split("\t")
        exon = ";".join(col[:4])
        if dico.has_key(exon):
            dico[exon].append(int(col[7]))
        else:
            dico[exon] = [int(col[7])]

    cles = dico.keys()
    cles.sort()
    for cle in cles:
        dico[cle]=np.asarray(dico[cle])
        #print "Moyenne : ",np.mean(dico[cle])
        #print "Ecart-type : ",np.std(dico[cle])

    return dico,cles

def patientMoyen(dico):
    patients = dico.keys()
    regions = dico[patients[0]].keys()
    #print len(regions)
    PatientMoyen={}
    PatientSTD={}
    for region in regions:
            chrom = region.split(";")[0].upper()
            if chrom.count("Y")==0 :
                PatientMoyen[region]=np.asarray([])
                for patient in patients :
                    PatientMoyen[region]=np.concatenate((PatientMoyen[region],dico[patient][region]))
                PatientSTD[region] = np.std(PatientMoyen[region])
                PatientMoyen[region]=np.mean(PatientMoyen[region])

    return PatientMoyen,PatientSTD

def getCov(path):
    allFiles=os.listdir(path)
    if path[-1]!="/":
        path+="/"
    dico={}
    boul = False
    for filin in allFiles:
        dico[filin.split("_")[0]],cles = cov2dico(path+filin)
    return dico

def covAjuste(dicoPat,DPmeanRef,DPmean):
    cov_ajuste={}
    patients = dicoPat.keys()
    regions = dicoPat[patients[0]].keys()
    for patient in patients :
            cov_ajuste[patient]={}
            for region in regions:
                chrom = region.split(";")[0].upper()
                if chrom.count("Y")==0  :
                    cov = np.mean(dicoPat[patient][region])
                    cov_ajuste[patient][region]=(DPmeanRef*cov)/DPmean[patient]
    return cov_ajuste

def calcDPmeanRef(temoin):
    DPmeanRef=[]
    for reg in temoin:
        chrom = reg.split(";")[0].upper()
        if chrom.count("Y")==0  :
            DPmeanRef.append(temoin[reg])
    DPmeanRef=np.mean(np.asarray(DPmeanRef))

    return DPmeanRef

def calcDPmean(dicoPat):
    DPmean={}
    patients = dicoPat.keys()
    regions = dicoPat[patients[0]].keys()
    for patient in patients :
        DPmean[patient]=np.asarray([])
        for reg in regions:
            chrom = reg.split(";")[0].upper()
            if chrom.count("Y")==0  :
                DPmean[patient]=np.concatenate((DPmean[patient],dicoPat[patient][reg]))
        DPmean[patient]=np.mean(DPmean[patient])
    return DPmean

def calcRapport(cov_ajuste,temoin):
    rapport={}
    patients = cov_ajuste.keys()
    regions = cov_ajuste[patients[0]].keys()
    for patient in patients:
        rapport[patient]={}
        for reg in regions:
            chrom = reg.split(";")[0].upper()
            if chrom.count("Y")==0   :
                rapport[patient][reg]=cov_ajuste[patient][reg]/temoin[reg]
    return rapport


if __name__ == '__main__':
    parser=argparse.ArgumentParser(usage="prog -p patCov -t temoCov")
    parser.add_argument("-p","--path",dest="path",help="Fichier de couverture des patients",default=None)
    parser.add_argument("-t","--temoins",dest="temoins",help="Repertoire des temoins")
    args=parser.parse_args()

    if args.path is None :
        parser.error("Path to cov files (-p) is missing.")
    else:
        path = args.path

    if args.temoins is None :
        parser.error("Path to temoins files (-p) is missing.")
    else:
        temFile = args.temoins

    #couvertures pour chaque exon chez TOUS les patients

    temoinDico = getCov(temFile)
    temoin,temoinstd = patientMoyen(temoinDico)

    dicoPat = getCov(path)
    dpmean,dpmeanstd = patientMoyen(dicoPat)

    DPmean=calcDPmean(dicoPat)
    DPmeanRef=calcDPmeanRef(temoin)
    cov_ajuste=covAjuste(dicoPat,DPmeanRef,DPmean)
    rapport=calcRapport(cov_ajuste,temoin)

    patients = dicoPat.keys()
    regions = dicoPat[patients[0]].keys()
    regions.sort()
    patients.sort()

    os.system("mkdir -p results")

    for patient in patients:
        fichier=open("results/patient"+patient+"_CNV.txt","w")
        fichier.write("chrom\tstart\tend\tGene\tMeanCov\tStdCov\tCV\tRapport\tAdjCov\n")
        for region in regions:
            chrom = region.split(";")[0].upper()
            if chrom.count("Y")==0  :
                mean=np.mean(dicoPat[patient][region])
                std =np.std(dicoPat[patient][region])
                CV = std/mean*100.0
                results = region.split(";")+[mean,std,CV,rapport[patient][region],cov_ajuste[patient][region]]
                results = [str(i) for i in results]
                ligne = "\t".join(results)+"\n"
                ligne = ligne.replace(".",",")
                fichier.write(ligne)
        fichier.close()
