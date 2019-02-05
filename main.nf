

params.amy   = "/opt/exp_soft/bioinf/AMY-tree"
params.mutes = "MutationConversion_v2.2.txt"
params.tree  = "UpdatedTree_v2.2.txt"
params.qc    = "qualityControl_v2.0.txt"
params.ref_dir = "/dataB/aux/37"

amy   = params.amy
mut1_ch = Channel.create()
mut_ch =   Channel.fromPath("$amy/${params.mutes}").tap(mut1_ch)




ref1_ch = Channel.fromPath("${params.ref_dir}/yref.fa")
ref2_ch = Channel.fromPath("${params.ref_dir}/yref.fa")
qc_ch   = Channel.fromPath("${amy}/${params.qc}")

plink_ch = Channel
   .fromFilePairs("${params.plink}.{bed,bim,fam}",size:3, flat : true){ file -> file.baseName}


tree_ch = Channel.fromPath("$amy/${params.tree}")



fam = new File("${params.plink}.fam").readLines().collect ([]) { it -> d   = it.split(' '); return d[0]+" "+d[1] }



process getMales {
  input:
     set val(base), file(bed), file(bim), file(fam) from  plink_ch
  output:
     file("male.ped") into ped_ch
     file("male.map") into map_ch
  script:
     """
     plink --bfile $base --keep-allele-order --filter-males --chr 24 --recode  --out male
     """
}


process getIndPlinkfile {
  input:
  set val(ped_str), file(map) from ped_ch.splitText().combine(map_ch)
  output:
     file("${out}.vcf") into vcf_ch
  script:
     genotype = ped_str.trim()
     p = ped_str.split(' ')
     out = p[0]+"_"+p[1]
     """
     echo $genotype > ped
     plink --keep-allele-order --ped ped --map $map  --recode vcf --out $out
     """
}


process why {
   module 'perl526'
  input:
     file(vcf) from vcf_ch
  output:
     file(why) into why_ch 
  publishDir "yhaplo/$base"
  script:
     base = vcf.baseName
     why  = "${base}_AMY-tree.txt"
     """
     perl $amy/WHY_v1.0.pl $vcf vcf
     """
}



process getStatusReference {
   module 'perl526'
  input:
   file(ref) from ref1_ch
   file(mutation) from mut_ch
  output:
   file(status) into status_ch
  publishDir  'yhaplo/'
  script:
    status="status.txt"
    """
    perl $amy/getStatusReference.pl $ref $mutation hg19 $status
    """
}




data_ch = status_ch.merge(tree_ch,ref2_ch,qc_ch,mut1_ch).combine(why_ch)



process amy {
   module 'perl526'
   time '92h'
   memory  '2G'
  input:
   set file(status), file(tree), file(ref), file(qc), file(mutation), file(why) from data_ch
  publishDir "yhaplo/"
  output:
   file("$ind/*")
  script:
   ind = why.baseName.replace("_AMY-tree","")
   """
   mkdir -p $ind
   perl $amy/AMY-tree_v2.0.pl $why $ind/ $tree $mutation $ref $status  $qc hg19
   """
}






