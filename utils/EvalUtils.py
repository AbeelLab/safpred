#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import numpy as np
from copy import deepcopy

def CalculateFmax(realpredictions, ytrue, thresholds, testproteins, goclasses, \
                  evalgoclass, removeroot=False, evalmode='partial'):
    n = len(testproteins)
    predictions = deepcopy(realpredictions)
    # Remove the root term for evaluations
    obsoleteterms = {'GO:0052312', 'GO:1902586', 'GO:2000775'}
    if removeroot:
        rootgo = ['GO:0008150','GO:0005575','GO:0003674']
    else:
        rootgo = []
    rootgo = obsoleteterms.union(rootgo)  
    pr = []
    rc = []
    ff = []
    covlist = []

    fmax = -np.inf
    numprot = 0

    for t in thresholds:
        precision = 0.0
        recall = 0.0

        # Total number of test proteins for which any prediction was made and 
        # a prediction with probability >= t
        numpredany = 0
        numpredthreshold = 0

        for testprot in testproteins:
            predflag = 0 # variable to check if there were any predictions for this query
            
            gopred = deepcopy(predictions[testprot]) # predictions
            gotrue = set(deepcopy(ytrue[testprot])) # real values
            gotrueclass = set(filter(lambda x: (goclasses[x] == evalgoclass) and \
                                     (x not in rootgo), gotrue))
            if len(gotrueclass) == 0:
                continue
            intersectionsize = 0
            numprot += 1

            # Keep track of the number of terms that belong to the desired class
            numclassterms = 0
            for goterm,goprob in gopred.items():
                if (goclasses[goterm] == evalgoclass) and (goterm not in rootgo):
                    numclassterms = numclassterms + 1

                    # Check threshold condition
                    if goprob >= t:
                        predflag = 1
                        if goterm in gotrueclass:
                            intersectionsize = intersectionsize + 1

            # If t = 0, override precision -- all terms in the ontology are predicted with probability 0
            if t == 0:
                precision = precision + (float(len(gotrueclass)) / len(goclasses))
            elif (numclassterms > 0):
                precision = precision + (float(intersectionsize) / numclassterms)
                
            if len(gotrueclass) > 0:
                recall = recall + (float(intersectionsize) / len(gotrueclass))
            
            numpredthreshold = numpredthreshold + predflag
            if numclassterms > 0:
                numpredany += 1

        if numpredthreshold > 0:
            precision = float(precision) / numpredthreshold
        else:
            precision = 0.0
        if evalmode == 'full':
            recall = float(recall) / numprot
        elif evalmode == 'partial':
            # Saw this option in the CAFA assessment tool
            recall = float(recall) / numpredany

        # If t is 0, override recall as 1
        if t == 0:
            recall = 1.0    
        
        pr.append(precision)
        rc.append(recall)
        covlist.append(float(numpredany) / float(numprot))

        if precision + recall > 0:
            ffcurr = (2 * precision * recall) / (precision + recall)
            ff.append((2 * precision * recall) / (precision + recall))
            fmax = max(fmax, ffcurr)
        else:
            ffcurr = 0.0
            ff.append(ffcurr)
            fmax = max(fmax, ffcurr)
    fmaxidx = np.argmax(ff)
    thmax = thresholds[fmaxidx]
    coverage = covlist[fmaxidx]
    
    return fmax, pr, rc, coverage, thmax


def CalculateSmin(realpredictions, ytrue, thresholds, testproteins, goclasses, \
                  goics, evalgoclass, removeroot=False):
    import math
    n = len(testproteins)
    predictions = deepcopy(realpredictions)
    # Remove the root term for evaluations
    obsoleteterms = {'GO:0052312', 'GO:1902586', 'GO:2000775'}
    if removeroot:
        rootgo = ['GO:0008150','GO:0005575','GO:0003674']
    else:
        rootgo = []
    rootgo = obsoleteterms.union(rootgo)  
    numprot = 0 # number of proteins with a GO term in this ontology    
    rulist = []
    milist = []
    smin = np.inf

    for t in thresholds:
        ru = 0.0
        mi = 0.0

        for testprot in testproteins:
            
            gotrue = set(deepcopy(ytrue[testprot])) # real values
            gotrueclass = set(filter(lambda x: (goclasses[x] == evalgoclass) and \
                                     (x not in rootgo), gotrue))
            if len(gotrueclass) == 0:
                continue
            numprot = numprot + 1
            # Fetch predictions for this test protein
            gopred = deepcopy(predictions[testprot]) # predictions
            gopredthreshold = set()
            numclassterms = 0
            for goterm,goprob in gopred.items():
                if (goclasses[goterm] == evalgoclass) and (goterm not in rootgo):
                    numclassterms = numclassterms + 1
                    # Check threshold condition
                    if goprob >= t:
                        gopredthreshold.add(goterm)
            for trueterm in gotrueclass:
                if trueterm not in gopredthreshold:
                    ru = ru + goics[trueterm]
            for predterm in gopredthreshold:
                if predterm not in gotrueclass:
                    mi = mi + goics[predterm]
        ru = float(ru) / numprot
        mi = float(ru) / numprot
        
        if (ru > 0) or (mi > 0):
            scur = math.sqrt((ru ** 2) + (mi ** 2))
            smin = min(smin, scur)

        rulist.append(ru)
        milist.append(mi)

    return smin, rulist, milist

def CalculateROCAUC(pr, rc):
    prarr = np.array(pr)
    rcarr = np.array(rc)

    sortedidx = np.argsort(rcarr)
    rcarrsorted = rcarr[sortedidx]
    prarrsorted = prarr[sortedidx]
    aupr = np.trapz(prarrsorted, rcarrsorted)
    return aupr
