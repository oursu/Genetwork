from optparse import OptionParser
import os
import subprocess
import re
import gzip

def main():
    parser=OptionParser()
    
    parser.add_option('--gtf',dest='gtf',default='/mnt/data/annotations/by_organism/human/hg19.GRCh37/GENCODE_ann/gencodeAnno/v19/gencode.v19.annotation.gtf.gz')
    parser.add_option('--bed',dest='bed',default='')
    opts,args=parser.parse_args()

    outname=opts.bed+os.path.basename(opts.gtf)+'.bed'
    out=open(outname,'w')

    for line in gzip.open(opts.gtf,'r').readlines():
        if line[0]=='#':
            continue
        items=line.strip().split('\t')
        element=items[2]
        if element!='gene':
            continue
        chromo,start,end,strand=items[0],items[3],items[4],items[6]
        chromo=re.sub('chr','',chromo)
        geneID=re.sub('"','',re.sub('gene_id "','',items[8].split(';')[0]))
        geneName=re.sub('"','',re.sub('gene_name "','',items[8].split(';')[4]))
        if strand=='+':
            tss=max(0,int(start)-1)
        if strand=='-':
            tss=max(0,int(end)-1)

        if len(geneName.strip())==0:
            geneName='NotAvailableGeneSymbol'+geneID
        out.write(chromo+'\t'+str(tss)+'\t'+str(int(tss)+1)+'\t'+geneName+';'+geneID+'\n')

    out.close()
    os.system('zcat -f '+outname+' | gzip > '+outname+'.gz')
    os.system('rm '+outname)

main()
