#!/usr/bin/env python3

import matplotlib
matplotlib.use('Agg')
import pandas as pd
import argparse
import matplotlib.pyplot as plt
import sys

TAB=chr(9)





def parseArguments():
    parser=argparse.ArgumentParser()
    parser.add_argument("--phe",  dest="phe", action="store", nargs=3, default=[])
    parser.add_argument("--tree", dest="tree", action="store", default="mt.csv")
    parser.add_argument("--cut-off", dest="cut_off", action="store", type=float,\
                        default=0.04)
    parser.add_argument("--labels",  dest="labels", action="store", \
                     default="",type=str,\
                     help="relabelling rules",\
                     metavar="str")    
    parser.add_argument("--normalise",  dest="normalise", action="store_true", default=False)
    parser.add_argument("--terminal-nodes",  dest="terminals", action="store", \
                     default="",type=str,\
                     help="Terminal nodes (comma-separated)",\
                     metavar="str")    
    parser.add_argument("haplos", action="store")
    parser.add_argument("out",   action="store")
    args = parser.parse_args()
    return args




args = parseArguments()

terminals = args.terminals.split(",")
spec_parent = {}
parent      = {}

def getSpecParent(fname):
    f = open(fname)
    for line in f:
        line=line.strip().split()
        spec_parent[line[0]]=line[1]
    f.close()


def getParent(node):
    if node in spec_parent:
        return spec_parent[node]
    if node in parent:
        return parent[node]
    if node in ["L","Root"]: return "Root"
    the_parent = node[:-1]
    if the_parent[-1]=="'" or len(the_parent)==0:
        sys.exit("Invalid parent for <%s>"%node)
    parent[node]=the_parent
    return the_parent




def  getClassification(want):
    def hg_sub(x):
        gname = hg.iloc[x]['Haplogroup']
        if "'" in gname:
            return gname
        else:
            return gname[:9]
    cats = hg[(hg['Quality']>0.5) &\
              (hg['SampleID'].isin(want))].groupby(hg_sub).count()['Range']
    return cats


def draw(fname,score):
    labels = sorted(score.keys())
    for label in labels:
        if label not in cdict:
            next_c = len(cdict) % len(colour_hex)
            cdict[label]=colour_hex[next_c]

    count  = [score[lab]  for lab in labels]
    fig1, ax1 = plt.subplots()
    these_colours = [cdict[l] for l in labels]
    ax1.pie(count, labels=labels, colors=these_colours,autopct='%1.1f%%', startangle=90)
    ax1.axis('equal')  
    plt.savefig("%s.pdf"%fname)

def completeTree(node,children):
    if node=="Root":return
    the_parent=getParent(node)
    if the_parent in children:
        if node not in children[the_parent]:
            children[the_parent].append(node)
        return
    children[the_parent]=[node]
    completeTree(the_parent,children)
    return


def makeTree(hgs):
    children={}
    for hg in hgs:
        p = getParent(hg)
        if p in children:
            children[p].append(hg)
        else:
            children[p]=[hg]
        if hg not in children: children[hg]=[]
    curr_children = list(children.keys())
    for c in curr_children:
        completeTree(c,children)
    return children

    
def scoreNodes(cats,children):
    labels=cats.index
    counts = cats.values
    cut_off = args.cut_off * sum(counts)
    score = {}
    #T2e1a1b T2e1a1b1
    def bubble(node, below_terminal):
        if node in terminals: below_terminal=True
        if node not in score: score[node]=0
        for c in children[node]:
            bubble(c,below_terminal)
            if below_terminal or score[c]<cut_off:
                score[node]=score[node]+score[c]
                score[c]=0
    for label,count in zip(labels,counts):
        score[label]=int(count)
    bubble("Root",False) # get score 
    root=score["Root"]
    score = { k:v for k,v in score.items() if v>=cut_off }
    if root>cut_off/4:
        score["L"]=root
    return score


colour_hex = ['#1f77b4', '#ff1f8e', '#2ca02c', '#d62728', '#9467bd', '#8c564b',\
              '#e377c2', '#7f7f7f', '#bcbd22', '#17becf', '#4c24b4',\
              '#ef8f8e','#bcfde2','#1070cf','#20af3c', '#a4478d']
cdict = {}

rdict={}
if args.labels:
    rules = args.labels.split(",")
    for r in rules:
        old,new=r.split("=")
        rdict[old]=new

getSpecParent(args.tree)
labels=["ALL"]
hg = pd.read_csv(args.haplos,delim_whitespace=True)
if args.phe:
  [phefname,column,cond]=args.phe
  phe=pd.read_csv(phefname,usecols=["FID",column],delim_whitespace=True,dtype=str)
  conditions = cond.split(",")
  interesting = [phe[phe[column].isin(conditions)]['FID']]
  for c in conditions:
     labels.append(rdict.get(c,c))
     interesting.append(phe[phe[column].isin([c])]['FID'])
else:
  interesting=[hg['SampleID']]

scores = {}
for i,want in enumerate(interesting):
    cats = getClassification(want)
    children = makeTree(cats.index)
    score =scoreNodes(cats,children)
    scores[labels[i]]=score
    draw("%s-%s"%(args.out,rdict.get(labels[i],labels[i])),score)

overall = pd.DataFrame.from_dict(scores,dtype=int).fillna(0).astype(int)
srt_col = sorted(overall.columns[1:])
if args.normalise:
    for col in overall.columns:
        overall[col]=overall[col]*100/sum(overall[col])
overall = overall.rename(rdict,axis=1).round(1).sort_index().reindex(['ALL']+srt_col,axis=1)
overall.to_csv(args.out+".csv",sep="\t")



