#!/usr/bin/env python314

import re
import sys
import sqlite3 as _sqlite3
from os.path import splitext,join
from ete3 import Tree, faces, AttrFace, TreeStyle, NodeStyle
import argparse
import math
from collections import defaultdict
import pandas as pd
from pandas_plink import read_plink
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


TAB=chr(9)
count = "___count___"
group_all = "___all___"

#12345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
def parseArguments():
    if len(sys.argv)<=1:
        sys.argv="mafplot.py $input $output".split()
    parser=argparse.ArgumentParser()
    parser.add_argument("--format",  dest="format", action="store", default="plink")
    parser.add_argument("--depth",  dest="depth", action="store", default=10000, type=int)
    parser.add_argument("--scan",  dest="scan", action="store_true", default=False)
    parser.add_argument("--table",  dest="table", action="store_true", default=False)
    parser.add_argument("--pdf",  dest="pdf", action="store_true", default=False,\
                                    help="produce PDF of the tree")    
    parser.add_argument("--result",  dest="result", action="store_true", default=False,\
                                    help="show text file of result")    
    parser.add_argument("--best",  dest="best", action="store_true", default=False,\
                                    help="use internal node instead of leaves if better")
    parser.add_argument("--group",  dest="group", action="store", default="",type=str,\
                        help="file with group info, name of id and group column",\
                        nargs=3, metavar="FILE COL")    
    parser.add_argument("--mingroup",  dest="mingroup", action="store", \
                                    default=0,type=int,\
                                    help="min size group should be",\
                                    metavar="INT")    
    parser.add_argument("--labels",  dest="labels", action="store", \
                                    default="",type=str,\
                                    help="relabelling rules",\
                                    metavar="str")    
    parser.add_argument("--points",  dest="points", action="store", default=15,type=int,\
                                    help="font size of labels")    
    parser.add_argument("--terminal-nodes",  dest="terminals", action="store", \
                        default="",type=str,\
                        help="comma separated list of terminal nodes")    
    parser.add_argument("--pie",  dest="pie", action="store", default=0, type=float,\
                                    help="Produce pie chart",metavar='min-wedge')    
    parser.add_argument("--normalise", dest="norm", action="store_true", default=False,\
                                    help="normalise summary files")    
    parser.add_argument("--prune",  dest="prune", action="store_true", default=False,\
                                    help="prune tree for PDF")    
    parser.add_argument("--show-probes",  dest="show_probes",\
                    action="store_true", default=False,\
                    help="show relevant SNPs in output file")    
    parser.add_argument("--show-bad-probes",  dest="show_bad_probes",\
                    action="store_true", default=False,\
                    help="show  SNPs inconsistent with call in output file")    
    parser.add_argument("--ignore", dest="ignore",action="store_true",default=False)
    parser.add_argument("--overall", dest="overall",action="store_true",default=False,help="Show overall tree")
    parser.add_argument("--cut", dest="cut",action="store",default=0,type=float)
    parser.add_argument("--num", dest="num",action="store",default=0,type=int,help="number of samples")
    parser.add_argument("--out", dest="out",action="store",default="output.txt")    
    parser.add_argument("--out-dir", dest="out_dir",action="store",default=".",help="output directory")    
    parser.add_argument("--score-node-only", dest="score_node_only",action="store_true",default=False)    
    parser.add_argument('tree', type=str, metavar='tree')
    parser.add_argument('mute', type=str, metavar='mute', help="list of mutations")
    parser.add_argument('sample', type=str, metavar='sample', help="sample",\
                                    nargs = '+')
    args = parser.parse_args()
    return args


def my_layout(node):
  faces.add_face_to_node(AttrFace("name"), node, column=0, position="branch-right")

args  = parseArguments()
ts = TreeStyle()
ts.show_leaf_name = False
ts.show_branch_support = False
ts.show_branch_length = False
ts.layout_fn = my_layout
ts.show_scale = False

def readTree(treef):
    marker = {}
    parent = {}
    treef.readline()
    n=0
    for line in treef:
        data    = line.strip().split()
        haplo   = data[0]
        alt     = data[1]
        the_par = data[2]
        snps    = data[3:]
        parent[haplo]=the_par
        for snp in snps:
            marker[snp]=haplo
        n=n+1
    return marker, parent


def recCvt(node, children):
    if node not in children:
        return node
    curr = []
    for c in children[node]:
        curr.append(recCvt(c,children))
    if len(curr)==1: curr.append("_")
    this_level = ",".join(curr)
    return "(%s)%s"%(this_level,node)

def cvtNewick(parent):
    children = {}
    root = []
    for k in parent.keys():
        siblings = children.get(parent[k],False)
        if siblings:
            siblings.append(k)
        else:
            children[parent[k]]=[k]
    newick = recCvt("Root",children)
    return children, newick+";"
    
def getMarkers(treef):
    treef.readline()
    marker_pos  = {}
    mute_alleles = {}
    for line in treef:
        data = line.strip().split("\t")
        name = data[0]
        pos  = data[2]
        mutation = data[3]
        ignore = data[5]
        if args.ignore and  ignore=="yes":
            continue
        if pos in marker_pos:
            marker_pos[pos].append(name)
        else:
            marker_pos[pos]=[name]
        mute_alleles[name]=mutation.split("->")
    return marker_pos, mute_alleles
        

def drawTree(outname,tree):
    tree.render(join(args.out_dir,outname+".pdf"), w=1383, units="mm", tree_style=ts)


def getColNums():
    if args.format == "amy":
        pos = 1;
        anc = 2;
        alt  = 3
    elif args.format == "bim":
        pos = 3
        anc = 5
        alt = 4
    elif args.format == "plink":
        pos=0
        anc=1
        alt=2
    else:
        sys.exit("Format <%s> unknown"%args.format)
    return (pos,anc,alt)


def getInd(p,bed,bim):
  gt = bed[:,p].compute()
  def cvt(s):
        name = "x" if len(bim.snp[s])==0 else bim.snp[s]
        if gt[s]==0:
            return (bim.pos[s],bim.a1[s],bim.a0[s],name)
        elif gt[s]==1:
            return (bim.pos[s],".",".",name)
        elif gt[s]==2:
            return (bim.pos[s],bim.a1[s],".",name)
        else:
            return (bim.pos[s],".",".",name)            
  return  list(map(cvt,range(len(bim))))



    
def classify(samplef):
    probes  = {}
    not_m   = {}
    pos_col, anc_col, call_col = getColNums()
    for line in samplef:
        if type(line)==str:
            data=line.strip().split()
        else:
            data=line
        pos = str(data[pos_col])
        anc = data[anc_col]
        alt  = data[call_col]
        if anc == "." or pos not in marker_pos: continue
        if alt == ".":
            call = anc
        else:
            call = alt;
        for the_mute in marker_pos[pos]:
           if the_mute not in marker: continue
           node = marker[the_mute]
           if args.scan or call == mute_alleles[the_mute][1]:
               if node in probes:
                   probes[node].append(the_mute)
               else:
                   probes[node]=[the_mute]
           elif (not args.scan):
              if node in not_m:
                  not_m[node].append(the_mute)
              else:
                  not_m[node]=[the_mute]
    for the_dict in [probes, not_m]:
        for k in the_dict.keys():
            the_dict[k]=list(set(the_dict[k]))
    return probes, not_m




nstyle = {}
nstyle[0]=NodeStyle()
nstyle[0]["shape"]="square"
nstyle[0]["size"]=4
nstyle[0]["fgcolor"]="blue"


def getLeavesScore(depth,TP,FN,node,probes,not_match):
    if node.name in probes:
        TP_curr = TP + probes[node.name]
    else:
        TP_curr = TP
    if node.name in not_match:
        FN_curr = FN + not_match[node.name]
    else:
        FN_curr = FN
    TP_set = set(TP_curr)
    FN_set = set(FN_curr)
    if node.is_leaf() or depth>= args.depth or node.name in terminals:
        return [[node.name,TP_set, FN_set]]
    res =[]
    for child in node.children:
        descs  = getLeavesScore(depth+1,TP_curr,FN_curr,child,probes,not_match)
        for this_descendant in descs:
            if len(TP_curr) < len(this_descendant[1]): 
                #only bother about descendants who have more info than us!
                res.append(this_descendant)
    if len(res)==0:
        return [[node.name,TP_set, FN_set]]
    return res

        
    
convertfn = (lambda x:x) if args.score_node_only else math.sqrt

def getSimpleStyle(score):
     curr_style = NodeStyle()
     curr_style["shape"] = "sphere"
     curr_style["size"] = score
     curr_style["fgcolor"] = "darkred"
     return curr_style

def recReweigh(curr,tr,probes):
    score = 0
    for child in curr.children:
        score=score+recReweigh(child,tr,probes)
    if args.score_node_only: score=0
    score=score+len(probes.get(curr.name,[]))
    if score not in nstyle:
        curr_style = NodeStyle()
        curr_style["shape"] = "sphere"
        curr_style["size"] = convertfn(5*score)+10
        curr_style["fgcolor"] = "darkred"
        nstyle[score]=curr_style
    else:
        curr_style = nstyle[score]
    curr.set_style(curr_style)
    return score


def reweight(curr_node,children,tr,probes):
      node = tr.get_tree_root()
      recReweigh(node, tr,probes)

      
def showTable(sample,probes):
    hgs = sorted(probes.keys())
    f=open(join(args.out_dir,sample),"w")
    for hg in hgs:
        the_probes = sorted(probes[hg])
        f.write("%s\t%s\n"%(hg," ".join(the_probes)))
    f.close()

def prune(tree, prune_set):
   def r_prune(node):
       keep = node.name in prune_set
       rem_list = []
       for child in node.children:
           child_got, child_list = r_prune(child)
           keep = keep or child_got
           rem_list = rem_list + child_list
       if keep:
           return (True, rem_list)
       else:
           return (False, [node]+rem_list)
   ok, del_nodes = r_prune(tree)
   return del_nodes

def produceClassification(out,base,tr,probes,not_match):
    out.write(base+TAB)
    result = getLeavesScore(0,[],[],tr,probes,not_match)    
    leaf_names, TPs, FNs = zip(*result)
    all_pos = set([]).union(*TPs)
    all_neg = set([]).union(*FNs)
    leaf_attrs= []
    best_leaf=""
    best_score=0
    for i, lname in enumerate(leaf_names):
        FP  = len(all_pos - TPs[i])
        TP  = len(TPs[i])
        F1 = 2*TP/(2*TP+len(FNs[i])+FP)
        if F1<args.cut: continue
        if F1>best_score:
            best_leaf=lname
            best_score=F1
        res = "%s:%d/%4.2f"%(lname,TP,F1)
        if args.show_probes:
            res=res+ "(%s)"%(",".join(TPs[i]))
        if args.show_bad_probes:
            res=res+ "(-%s)"%(",".join(FNs[i]))
        leaf_attrs.append((F1,res))
    leaf_attrs.sort(reverse=True)
    leaf_data=list(map(lambda x:x[1],leaf_attrs))
    if args.best: leaf_data=[leaf_data[0]]
    out.write(",".join(leaf_data))
    out.write("\n")
    return best_leaf
    

def processSample(out,sample_name,sample_data):
   probes, not_match  = classify(sample_data)
   reweight("Root",children,tr,probes)
   m = re.search(".*/(.*)\..*",sample_name)
   if m:
       base = m.group(1)
   else:
       base = splitext(sample_name)[0]
   if args.table:
       showTable(base,probes)
   prune_set = probes.keys()
   if args.pdf and len(prune_set)>0:
       if args.prune:  # for this part option of pruning
           rem_list = prune(tr,prune_set)
           for r in rem_list:
             r.detach()
       drawTree(base,tr)
   if args.result:
       if len(prune_set)==0:
           out.write("%s%sNIL\n"%(base,TAB))
           return
       if not args.prune: 
           # for this we only want a pruned tree  so prune if not done yet
           rem_list = prune(tr,prune_set)
           for r in rem_list:
               r.detach()
       best_leaf = produceClassification(out,base,tr,probes,not_match)
       if len(best_leaf)>=1:
          for gfn in groups:
              this_g = gfn(sample_name)
              if this_g not in overall: 
                  overall[this_g]=defaultdict(int)
              overall[this_g][best_leaf]=overall[this_g].get(best_leaf,0)+1
              overall[this_g][count]=overall[this_g].get(count,0)+1


def presentOverall(tree):
    want = overall[group_all].keys()
    for name in want:
        if name==count:continue
        these = tree.search_nodes(name=name)
        style  = getSimpleStyle(overall[group_all][name])
        these[0].set_style(style)
    del_nodes = prune(tree,want)
    for node in del_nodes: node.detach()
        

def getPie(depth,node,group):
    """Finds the wedges of a pie chart.
       
       returns a pair 
         first element is a list of wedges [(name, weight), (name, weight), ... ]
         second element   integer -- unallocated weights
    """
    N=overall[group][count]
    weight = overall[group].get(node.name,0)
    if node.is_leaf() or depth>args.depth or node.name in terminals:
        if weight/N*100 >= args.pie:
            return ([(node.name,weight)],0)
        else:
            return  ([], weight)
    wedges  = []
    for child in node.children:
        (w,u) = getPie(depth+1,child,group)
        wedges = wedges + w
        weight = weight+u
    if weight/N*100 >= args.pie: 
        # weight of node + children above cut-off
        wedges.append((node.name, weight))
        weight = 0  # all weight allocated
    return (wedges, weight)


def drawWedges(outbase,group,cdict,wedges):
    slices = wedges[0] #  wedges[0] wedges big enough to show
    if wedges[1]*overall[group][count]/100>2: #wedges[1] size of "other"
        slices.append(("Other",wedges[1]))
    labels, weights = zip(*slices)
    colours = [cdict[label] for label in labels]
    fig1, ax1 = plt.subplots()
    ax1.pie(weights, labels=labels, autopct='%1.0f%%',startangle=90,\
            colors = colours,\
            textprops={'fontsize': args.points})
    ax1.axis('equal') 
    if group==group_all: 
        g_fname=""
    else:
        g_fname=rdict.get(group,group)+"-"
    plt.savefig(join(args.out_dir,"%s-%spie.pdf")%(outbase,g_fname))
    g=open(join(args.out_dir,"%s-%ssumm.csv"%(outbase,g_fname)),"w")
    for i, label in enumerate(labels):
        g.write("%s\t%d\n"%(label,weights[i]))
    g.write("Total\t%d\n"%(overall[group][count]))
    g.close()

    

def getGroups(fname, id_col, phe_col):
    gdf = pd.read_csv(fname,delim_whitespace=True,dtype={phe_col:str})
    err_msg = ""
    for col in [id_col, phe_col]:
        if  col not in gdf.columns:
            err_msg=err_msg+\
                     "The column <%s> cannot be found in the file <%s>"%(col,fname)
    if err_msg: sys.exit(err_msg)
    gdf.set_index(id_col,inplace=True)
    return gdf


def getFixedColours(wedges):
    ''' need to keep colours consistent as not all haplogroups appear in all
        e.g., to be there must appear at a frequency above the cut-off so could
        appear in a sub-group but not the overall list '''
    ccount = {}
    for (slices, rest) in wedges:
        for (label, weight) in slices:
            ccount[label]=ccount.get(label,0)+1
    leaf_names  = [lf for lf in ccount.keys()]
    leaf_names.sort(key=lambda leaf: ccount[leaf],reverse=True)
    colour_hex = ['#1f77b4', '#ff1f8e', '#2ca02c', '#d62728', '#9467bd', '#8c564b',\
                  '#e377c2', '#7f7f7f', '#bcbd22', '#17becf', '#4c24b4',\
                  '#ef8f8e','#bcfde2','#1070cf','#20af3c' ]
    cdict = dict(zip(leaf_names,colour_hex))
    if "Other" not in cdict:
        cdict['Other']=colour_hex[len(cdict)]
    return cdict

rdict={}
if args.labels:
    rules = args.labels.split(",")
    rdict[group_all]="ALL"
    for r in rules:
        old,new=r.split("=")
        rdict[old]=new


terminals = args.terminals.split(",")        

groups = [lambda x: group_all]
if args.group:
    gdf=getGroups(*args.group)
    def getfn(sample):
        g_n=gdf.loc[sample][args.group[2]]
        if str(g_n) in ["NA","nan"]:
            g_n="none"
        return g_n
    groups.append(getfn)
                       


overall={} 
marker, parent = readTree(open(args.tree))
marker_pos, mute_alleles = getMarkers(open(args.mute))
children, tree_string=cvtNewick(parent)

orig = Tree(tree_string,format=8)

out=open(join(args.out_dir,args.out),"w")
outbase =  splitext(args.out)[0]

if args.format=="plink":
    if len(args.sample)!=1:
        sys.exit("When using PLINK fomat there can only be some sample -- the base of a plink data set")
    (bim,fam,bed) =read_plink(args.sample[0],verbose=False)
    N=len(fam) if args.num==0 else args.num
    for p in range(N):
        tr = orig.copy(method='deepcopy')
        if fam.fid[p]==fam.iid[p]:
            the_id = fam.fid[p]
        else:
            the_id = fam.fid[p]+"_"+fam.iid[p]
        print(p,the_id)
        data = getInd(p,bed,bim)
        processSample(out,the_id,data)
else:
    N=len(args.sample) if args.num==0 else args.num
    for sample in args.sample[:N]:
        tr = orig.copy(method='deepcopy')
        processSample(out,sample,open(sample))
out.close()


group_names = [g for g in overall.keys() if overall[g][count]>0]
if args.pie:
    all_wedges = []
    for group in group_names:
        all_wedges.append(getPie(0,orig,group))
    # now we have all the wedges, we can fix the colours
    cdict = getFixedColours(all_wedges)
    # now draw theedges
    for i,group in enumerate(group_names):
        drawWedges(outbase,group,cdict, all_wedges[i])


#if not args.overall: sys.exit(0)
#presentOverall(orig)
#fn=outbase+"-overall"
#drawTree(fn,orig)

lnames = set()
# zero out any values for which the leaf does not appear
for g in group_names:
    del overall[g][count]
    lnames = lnames | set(overall[g].keys())
for g in group_names:
    for lf in lnames:
        if lf not in overall[g]: overall[g][lf]=0


classification = pd.DataFrame.from_dict(overall,dtype=int)
classification.rename(rdict,inplace=True,axis=1)
classification.to_csv(join(args.out_dir,outbase)+"-table.csv",sep="\t")
for col in classification.columns:
  classification[col]=classification[col]*100/sum(classification[col])
  classification=classification.round(1)

classification.to_csv(join(args.out_dir,outbase)+"-norm-table.csv",sep="\t")
    

