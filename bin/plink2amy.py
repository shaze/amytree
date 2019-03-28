#!/usr/bin/env /usr/bin/python36

from pandas_plink import read_plink
import sys


(bim,fam,bed) =read_plink(sys.argv[1],verbose=False)

bed = bed.compute()



def getInd(p):
  gt = bed[:,p]
  def cvt(s):
        pref = "chrY\t%d\t%s\t%s\t%s\n"
        name = "x" if len(bim.snp[s])==0 else bim.snp[s]
        if gt[s]==0:
            return pref%(bim.pos[s],bim.a1[s],bim.a0[s],name)
        elif gt[s]==1:
            print(bim.snp[s],gt[s],"hz")
            return pref%(bim.pos[s],".",".",name)
        elif gt[s]==2:
            return pref%(bim.pos[s],bim.a1[s],".",name)
        else:
            return pref%(bim.pos[s],".",".",name)            
  res =  list(map(cvt,range(len(bim))))
  return res

for p in range(len(fam)):
    if fam.fid[p]==fam.iid[p]:
        the_id = fam.fid[p]
    else:
        the_id = fam.fid[p]+"_"+fam.iid[p]
    g=open("%s.amy"%the_id,"w")
    g.writelines(getInd(p))
    g.close()
