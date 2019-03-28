

params.problems = 0
mut_ch  =   file(params.mutes)
ref_ch  = file(params.refy)
qc_ch   = file(params.qc)
tree_ch = file(params.tree)


empties = [false,"False","false", "FALSE",0,"","0"]

if (empties.contains(params.problems))
    problems_ch = Channel.from("main.nf") 
else
    problems_ch = file(params.problems)

plink_ch = Channel
   .fromFilePairs("${params.plink}.{bed,bim,fam}",size:3, flat : true){ file -> file.baseName}





fam = new File("${params.plink}.fam").readLines().collect ([]) { it -> d   = it.split(' '); return d[0]+" "+d[1] }



process getMales {
  input:
     set val(base), file(bed), file(bim), file(fam) from  plink_ch
     file (problems) from problems_ch
  output:
     set file("males.bed"), file("males.bim"), file("males.fam") into males_ch
  script:
     if (empties.contains(params.problems))
        prob = "" 
     else
        prob = "--exclude $problems"
     """
     plink --bfile $base $prob --maf 0.075 --geno 0.1 --keep-allele-order --filter-males --chr 24 --make-bed --out males
     """
}


process convertAmy {
  input:
  set file(bed), file(bim), file(fam) from males_ch
  output:
     file("*amy") into amy_ch mode 'flatten'
  publishDir 'yhaplo/amy'
  script:
     base = bed.baseName
     """
     plink2amy.py $base
     """
}


process getStatusReference {
  input:
   file(ref) from ref_ch
   file(mutation) from mut_ch
  output:
   file(status) into status_ch
  publishDir  'yhaplo/'
  script:
    status="status.txt"
    """
     getStatusReference.pl $ref $mutation hg19 $status
    """
}




data_ch = status_ch.combine(amy_ch)



process amy {
   memory  '1G'
  input:
   set file(status), file(amy) from data_ch
   file (tree) from tree_ch
   file (ref)  from ref_ch
   file (mutation) from mut_ch
   file (qc)   from qc_ch
  publishDir "yhaplo/trees"
  output:
  set val(base), file("*analysis.txt") into amy_result_ch
  script:
   base=amy.baseName
   """
    AMY-tree_v2.0.pl $amy ./ $tree $mutation $ref $status  $qc hg19
   """
}

process grabIndivResult {
  input:
    set val(base), file(res) from amy_result_ch
  output:
    file ("${base}.rsy") into overall_result_ch
  script:
    """
    grab_result.py $base $res ${base}.rsy
    """
}


process summariseResult {
  input:
    file(results) from overall_result_ch.toSortedList()
  publishDir "yhaplo/"
  output:
  file(out)
  script:
   out = params.out
   """
   cat *rsy > $out
   """
}



