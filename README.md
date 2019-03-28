# ymthaplotools




## plotY.py

This code does Y-haplogroup assigmment and has some utilities for analysing and viewing the results. 

I started with AMY-TREE, which was very useful but I ran into some trouble and couldn't easly change the code.
It takes as input the tree and mutation files that AMY-TREE uses.

I tested the code against what AMY-TREE produced and it was broadly consistent. I checked by hand discrepant samples and preferred my results.

Later, someone suggested SNAPPY and I then ran SNAPPY against the same data. The results are very similar. The differences are used by slightly different databases used (I used the AMY-TREE data with some small amendments)


### Installing

Uses the ete3 package which has a number of requirements that I found difficult to install.
* On my CentOS system I had to have the packags libxkbcommon-x11 and libxkbcommon-x11-devel insatlled


### Running

This code takes a data set and reference files and classifies samples according to their Y-haplogroup, showing alternatives

```
usage: plotY.py [-h] [--format FORMAT] [--scan] [--table] [--pdf] [--result]
                [--prune] [--show-probes] [--show-bad-probes] [--ignore]
                [--cut CUT] [--out OUT] [--score-node-only]
                tree mute sample [sample ...]

positional arguments:
  tree
  mute               list of mutations
  sample             sample

optional arguments:
  -h, --help         show this help message and exit
  --format FORMAT
  --scan
  --table
  --pdf              produce PDF of the tree
  --result           show text file of result
  --prune            prune tree for PDF
  --show-probes      show relevant SNPs in output file
  --show-bad-probes  show SNPs inconsistent with call in output file
  --ignore
  --cut CUT
  --out OUT
  --out-dir DIR   
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
* `format` The format of the input: 
     * `plink` -- the default. The data can be found in a PLINK data set (bed,bim,fam). In this case, the `sample`
        parameter should consist of exactly one string which is the base of the PLINK data set. All individuals
        in the PLINK data set are analysed. No checking of the PLINK file is made so it is the user's responsibility
        to ensure it is only Y-data and only contains males.
     * `amy` The file is in the format that AMY-tree expects. In this case there should be at least 1 but could be many        (e.g. specified with globbing) data files. Each person in their own file
     * `bim`  The input is a PLINK bim file. This really only intended to be used for the "scan" option below.
* `out`  the name of the output file -- default is "output.txt"
* `out-dir` directory into which all output should be placed -- default is current working directory
* `result` put a one line summary in the output file for each sample. Each line contains the ID of the sample, followed by a comma separated list of haplogroup assignments. After each assignment is the number of TPs and the F1 score.
* `show-probes` show the TPs for each tree
* `show-bad-probes` show the FNs for each tree (thes are shown with a "-" in front) 
* `pdf` for every sample, create a PDF picture of the tree, nodes that are covered are put in red, size adjusted by coverage
* `prune` If producing PDF only show the path of the haplogroup
* `score-node-only` for the PDF option, the score of the node is increased by the score of descendant nodes. When this option is chosen only the score of the actual node is used
* `cut` Only show haplogroups with F1-score above the cut-off (default is 0)
* `ignore
* `scan` this is used to show the potential nodes in the Y-tree that are covered by the BIM file or the AMY data rather than doing classifying
* `table` produces for each sample given a text file showing all mutations that match the tree
* `overall` if chosen, produces a PDF file summarising the haplogroup assignments for all the sample. The nodes in red are covered by samples, size of node adjusted by number of samples.
* `num` takes an integer N and only analyses the first N samples. Probably only useful for debugging the code


## mtparse.py

Takes two arguments
* XML MT tree (or I suppose a Y tree) in the format provided by Haplogrep
* A list of SNP positions on a chip

Produces a list of haplogroups covered by the chip -- for each haplogroup, which SNPs on the chip can be found