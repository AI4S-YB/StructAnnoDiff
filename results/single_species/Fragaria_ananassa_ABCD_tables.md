# F. ananassa A/B/C/D analysis tables

## Figure A table: global quantities

figure,layer,metric,metric_type,before,after,delta,pct_change
A,Global annotation quantity,Genes,count,"108,087","108,447",360,0.33
A,Global annotation quantity,mRNAs,count,"108,087","108,447",360,0.33
A,Global annotation quantity,Exons,count,"582,088","558,292","-23,796",-4.09
A,Global annotation quantity,CDS features,count,"108,087","108,447",360,0.33
A,Global annotation quantity,5' UTR features,count,"47,356","63,020","15,664",33.08
A,Global annotation quantity,3' UTR features,count,"51,473","64,952","13,479",26.19
A,Global annotation quantity,Single-exon genes,count,"17,715","17,526",-189,-1.07
A,Global annotation quantity,Mean exons per mRNA,mean,5.40,5.10,-0.30,-5.56
A,Global annotation quantity,Mean CDS length (bp),mean,"1,144","1,087",-57,-4.98
A,Global annotation quantity,Mean gene length (bp),mean,"3,158","2,992",-166,-5.26
A,Global annotation quantity,Median gene length (bp),median,"2,275","2,359",84,3.69


## Figure B table: locus fate

figure,layer,category,before_genes,before_pct,after_genes,after_pct
B,Locus fate,Syntenic 1:1,"96,298",89.09,"96,298",88.80
B,Locus fate,Split,"3,964",3.67,"8,401",7.75
B,Locus fate,Merge,671,0.62,327,0.30
B,Locus fate,Complex,358,0.33,396,0.37
B,Locus fate,Unresolved weak-overlap,0,0,0,0
B,Locus fate,Strict unmatched,"6,796",6.29,"3,025",2.79


## Figure C table: 1:1 structural attributes

figure,layer,attribute,count,pct_of_syntenic
C,Confirmed 1:1 structural attributes,Exact 1:1,"50,854",52.81
C,Confirmed 1:1 structural attributes,Gene boundary changed,"43,122",44.78
C,Confirmed 1:1 structural attributes,UTR added,"13,299",13.81
C,Confirmed 1:1 structural attributes,UTR lost,796,0.83
C,Confirmed 1:1 structural attributes,UTR exon gained,"7,125",7.40
C,Confirmed 1:1 structural attributes,UTR exon removed,"3,959",4.11
C,Confirmed 1:1 structural attributes,UTR refined,"34,352",35.67
C,Confirmed 1:1 structural attributes,Coding exon gain,"3,854",4.00
C,Confirmed 1:1 structural attributes,Coding exon loss,"10,676",11.09
C,Confirmed 1:1 structural attributes,Exon boundary refined,"26,628",27.65
C,Confirmed 1:1 structural attributes,CDS changed,"21,085",21.90
C,Confirmed 1:1 structural attributes,CDS boundary refined,222,0.23
C,Confirmed 1:1 structural attributes,Isoform changed,"9,023",9.37


## Figure D table: paired change magnitude

figure,layer,metric,n_pairs,mean_before,mean_after,mean_delta,median_delta,p05_delta,p95_delta,mean_pct_delta,decreased,unchanged,increased,changed,changed_pct
D,Confirmed 1:1 change magnitude,Gene span length,"96,298","3,010.00","2,939.86",-70.14,0,-901,749,6.05,"21,081","51,292","23,925","45,006",46.74
D,Confirmed 1:1 change magnitude,Model span length,"96,298","3,010.00","2,933.24",-76.77,0,-920,742.15,5.91,"21,239","51,294","23,765","45,004",46.73
D,Confirmed 1:1 change magnitude,CDS length,"96,298","1,131.32","1,238.75",107.43,0,-219,"1,032",10.02,"12,278","71,663","12,357","24,635",25.58
D,Confirmed 1:1 change magnitude,Exon count,"96,298",5.23,6.46,1.23,0,-1,7,14.72,"8,233","75,129","12,936","21,169",21.98
D,Confirmed 1:1 change magnitude,mRNA count,"96,298",1,1.16,0.16,0,0,1,15.96,0,"87,275","9,023","9,023",9.37

