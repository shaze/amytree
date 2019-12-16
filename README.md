# ymthaplotools




## plotY.py

This code does Y-haplogroup assigmment and has some utilities for analysing and viewing the results. 

I started with AMY-TREE, which was very useful but I ran into some trouble and couldn't easly change the code.

But this code takes as input the tree and mutation files that AMY-TREE uses.

I tested the code against what AMY-TREE produced and it was broadly consistent. I checked by hand discrepant samples and preferred my results.

Later, someone suggested SNAPPY and I then ran SNAPPY against the same data. The results are very similar. The differences are used by slightly different databases used (I used the AMY-TREE data with some small amendments)


### Installing

Uses the ete3 package which has a number of requirements
* On my CentOS system I had to have the packags libxkbcommon-x11 and libxkbcommon-x11-devel installed
* To run you need to have an active X11 session. If you don't want that, e.g. when running remotely, I have found success in using a  virtual X11 session
   * Install xorg-x11-server-Xvfb
   * run like this `xvfb-run python3 bin/plotY.py .... `

### Running

This code takes a data set and reference files and classifies samples according to their Y-haplogroup, showing alternatives

```
usage: plotY.py [-h] [--format FORMAT] [--depth DEPTH] [--scan] [--table]
                [--pdf] [--result] [--best] [--points POINTS]
                [--pie min-wedge] [--prune] [--show-probes]
                [--show-bad-probes] [--ignore] [--overall] [--cut CUT]
                [--num NUM] [--out OUT] [--out-dir OUT_DIR]
                [--score-node-only]
                tree mute sample [sample ...]

positional arguments:
  tree
  mute               list of mutations
  sample             sample

optional arguments:
  -h, --help         show this help message and exit
  --format FORMAT
  --depth DEPTH
  --scan
  --table
  --pdf              produce PDF of the tree
  --result           show text file of result
  --best             use internal node instead of leaves if better
  --points POINTS    font size of labels
  --group FNAME COL  file with group info and column name
  --pie min-wedge    Produce pie chart
  --prune            prune tree for PDF
  --show-probes      show relevant SNPs in output file
  --show-bad-probes  show SNPs inconsistent with call in output file
  --terminal-nodes   list of terminal node
  --labels            comma separated list of relabellings
  --ignore
  --overall          Show overall tree
  --cut CUT
  --num NUM          number of samples
  --out OUT
  --out-dir OUT_DIR  output directory
  --score-node-only
```



The compulsory arguments are
* `tree`: A tree in the form that AMY-TREE expects
* `mute`: A file with mutations in the form that AMY-tree expects
* `sample`: A list of 1 or more sample files.

The program takes each sample tries to classify it.  There are different output options available but in essence it will return any haplogroups for which there is evidence that that the sample belongs to. Note that there could be multiple classifications  where there is some ambiguity in the tree, or errors in the data. For each classification the F1 score is computed. Each SNP in the data can be classified with respect to each haplogroup as a
* a true positive (TP) - evidence for that haplogroup
* a true negative (TN) - evidence against _other_ haplogroups 
* a false negative (FN) - a call on the path from the root to the haplogroup node called which is inconsistent with 
  the haplogroup call. Typically this means that we expect to see a mutation at this point but don't
* a false positive (FP) -- a call that is evidence for some other haplogroup

There are many formulas available to classifying the quality

See https://en.wikipedia.org/wiki/Sensitivity_and_specificity for a discussion. In this application even a very poor classifier is likely to have a very high TN value so I think it best to use a formula that combines TP, FP and FN and the F1 score is one such formula


The optional arguments are explained below. If not given the default values are used
* `best` By default, any leaf nodes in the tree that have an F1 score above the cut-off are returned. If `simplify` is chosen, then we try return the single best node. The heuristic used is at each step at each node in the tree to return the child of the node, or the node itself with the best sensitivity and if the sensitivity is the same the best 
* `cut` Only show haplogroups with F1-score above the cut-off (default is 0)
* `format` The format of the input: 
     * `plink` -- the default. The data can be found in a PLINK data set (bed,bim,fam). In this case, the `sample`
        parameter should consist of exactly one string which is the base of the PLINK data set. All individuals
        in the PLINK data set are analysed. No checking of the PLINK file is made so it is the user's responsibility
        to ensure it is only Y-data and only contains males.
     * `amy` The file is in the format that AMY-tree expects. In this case there should be at least 1 but could be many        (e.g. specified with globbing) data files. Each person in their own file
     * `bim`  The input is a PLINK bim file. This really only intended to be used for the "scan" option below.
* `--depth`  Default very big. How deep in the tree should samples be classifed -- see `terminal-nodes` for more detail.
* `--group`  Takes three arguments: (1) the name of a file with group information. The columns of the file should be tab or space separated, and the first line should have headings. (2) The label of one of the columns which contains the  ID of the individuals; and (3) the label of a column used to categorise the samples. If this sample is given then the overall charts are categorised by group fond in the column. See also the `mingroup` and `labels` options which modify this.
* `--ignore` Ignore the SNPs that AMY-tree would ignore -- default is NOT (so different behaviour to AMY-tree)
* ``--labels` (only useful if the `--group` option given) a comma separated list of relabelling rules of the form old1=new1,old2=new2,old3=new3. The goal of this is to make more human friendly output. For example, if  your group file has groups 1,2,3 you might want this displayed as Zulu,Xhosa,Venda. You then say `--labels 1=Zulu,2=Xhosa,3=Venda`
* `--num` takes an integer N and only analyses the first N samples. Probably only useful for debugging the code
* `--out-dir` directory into which all output should be placed -- default is current working directory
* `--out`  the name of the output file -- default is "output.txt"
* `--overall` if chosen, produces a PDF file summarising the haplogroup assignments for all the sample. The nodes in red are covered by samples, size of node adjusted by number of samples.
* `--pdf` for every sample, create a PDF picture of the tree, nodes that are covered are put in red, size adjusted by coverage
* `--pie n`: Produce a pie chart (PDF format) of the distribution of haplotypes of all the samples. The parameter is a floating point number that describes the minimum percentage that a slice in the pie must have in order to be shown. If a slice has less than the given percentage, then the percentage is is allocated to the parent haplogroup. If  there is a slice for E1b with 15 and E1 for 12 this means, that 15 per cent of the sample have haplotype E1b and 12 percent of the sample have a haplogroup in E1 _other_ than E1b. At the root of tree an "Other" slice is shown if needed provided the "Other" unallocated component is less than 2 per cent. The output PDF file is put in the directory specified by `out_dir` and has a name with prefix given by `out` and the suffix `-pie.pdf`.
* `--points n`: The font size of the labels in the pie chart
* `--prune` If producing PDF only show the path of the haplogroup
* `--result` put a one line summary in the output file for each sample. Each line contains the ID of the sample, followed by a comma separated list of haplogroup assignments. After each assignment is the number of TPs and the F1 score.
* `--scan` this is used to show the potential nodes in the Y-tree that are covered by the BIM file or the AMY data rather than doing classifying
* `--score-node-only` for the PDF option, the score of the node is increased by the score of descendant nodes. When this option is chosen only the score of the actual node is used
* `--show-bad-probes` show the FNs for each tree (thes are shown with a "-" in front) 
* `--show-probes` show the TPs for each tree
* `--table` produces for each sample given a text file showing all mutations that match the tree
* `--terminal-nodes` a comma separated list of nodes to be treated as terminal or leaf nodes. For some work, rather than going as deep as possible `E1ba1ba1a2` you might prefer to stop at `E1b`. There are two ways of doing this, using `depth` which is an integer parameter specifying how deep in the tree you should traverse. However, this is a crude measure so better but more work is to manually specify if any nodes should be treated as terminal

###

Example run

xvfb-run python3 bin/plotY.py --best  --prune --pdf --pie 9 --result --depth 9  /opt/exp_soft/bioinf/AMY-tree/african_tree.txt /opt/exp_soft/bioinf/AMY-tree/african_mutation.txt sa-males --out-dir sa-males --out sa-males.txt


## pimt.py

This program takes the output of Haplogrep and produces a variery of pie charts and tables

```
usage: pimt.py [-h] [--phe PHE PHE PHE] [--tree TREE] [--cut-off CUT_OFF]
               [--labels str] [--normalise] [--terminal-nodes str]
               haplos out

positional arguments:
  haplos
  out

optional arguments:
  -h, --help         show this help message and exit
  --phe PHE PHE PHE
  --tree TREE
  --cut-off CUT_OFF
  --labels str       relabelling rules
  --terminal-nodes str    Terminal nodes (comma-separated)
```

The compulsory arguments are
* `haplos`: the name of an output file produced by Haplogrep
* `out`: base name of all output files

The optional arguments are
* `--phe` which takes three parameters: the name  of phenotype file, the label of a column in that file, and a comma separated list of the elements in that column which are of interest. If this option is given then only those elements in the _haplos_ file that match the specified values are used. In addition, a categorised analysis is done (1) with all the elements that match any of the conditions, and then (2) for each condition separately.
*  `--labels` (see plotY above)
*  `--terminal-nodes` (see plotY above)


## mtparse.py

Takes two arguments
* XML MT tree (or I suppose a Y tree) in the format provided by Haplogrep
* A list of SNP positions on a chip

Produces a list of haplogroups covered by the chip -- for each haplogroup, which SNPs on the chip can be found