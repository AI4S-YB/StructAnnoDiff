# Cucumber A/B/C/D analysis tables

## Figure A table: global quantities

figure,layer,metric,metric_type,before,after,delta,pct_change
A,Global annotation quantity,Genes,count,"24,317","24,145",-172,-0.71
A,Global annotation quantity,mRNAs,count,"24,317","24,145",-172,-0.71
A,Global annotation quantity,Exons,count,"127,080","127,343",263,0.21
A,Global annotation quantity,CDS features,count,"24,317","24,142",-175,-0.72
A,Global annotation quantity,5' UTR features,count,"15,483","20,636","5,153",33.28
A,Global annotation quantity,3' UTR features,count,"15,449","20,636","5,187",33.57
A,Global annotation quantity,Single-exon genes,count,"4,893","5,299",406,8.30
A,Global annotation quantity,Mean exons per mRNA,mean,5.20,5.30,0.10,1.92
A,Global annotation quantity,Mean CDS length (bp),mean,"1,128","1,130",2,0.18
A,Global annotation quantity,Mean gene length (bp),mean,"4,001","4,124",123,3.07
A,Global annotation quantity,Median gene length (bp),median,"2,883","3,028",145,5.03


## Figure B table: locus fate

figure,layer,category,before_genes,before_pct,after_genes,after_pct
B,Locus fate,Syntenic 1:1,"19,965",82.10,"19,965",82.69
B,Locus fate,Split,4,0.02,8,0.03
B,Locus fate,Merge,2,0.01,1,0.00
B,Locus fate,Complex,0,0,0,0
B,Locus fate,Unresolved weak-overlap,"2,300",9.46,"2,092",8.66
B,Locus fate,Strict deleted / novel,"2,046",8.41,"2,079",8.61


## Figure C table: 1:1 structural attributes

figure,layer,attribute,count,pct_of_syntenic
C,Confirmed 1:1 structural attributes,Exact 1:1,"15,508",77.68
C,Confirmed 1:1 structural attributes,Gene boundary changed,"4,285",21.46
C,Confirmed 1:1 structural attributes,UTR added,"1,940",9.72
C,Confirmed 1:1 structural attributes,UTR lost,2,0.01
C,Confirmed 1:1 structural attributes,UTR exon gained,"1,206",6.04
C,Confirmed 1:1 structural attributes,UTR exon removed,770,3.86
C,Confirmed 1:1 structural attributes,UTR refined,"1,778",8.91
C,Confirmed 1:1 structural attributes,Coding exon gain,320,1.60
C,Confirmed 1:1 structural attributes,Coding exon loss,"1,076",5.39
C,Confirmed 1:1 structural attributes,Exon boundary refined,762,3.82
C,Confirmed 1:1 structural attributes,CDS changed,"1,792",8.98
C,Confirmed 1:1 structural attributes,CDS boundary refined,36,0.18
C,Confirmed 1:1 structural attributes,Isoform changed,0,0


## Figure D table: paired change magnitude

figure,layer,metric,n_pairs,mean_before,mean_after,mean_delta,median_delta,p05_delta,p95_delta,mean_pct_delta,decreased,unchanged,increased,changed,changed_pct
D,Confirmed 1:1 change magnitude,Gene span length,"19,965","4,537.74","4,490.16",-47.58,0,-614.60,642.80,1.28,"2,574","15,582","1,809","4,383",21.95
D,Confirmed 1:1 change magnitude,Model span length,"19,965","4,460.69","4,489.89",29.20,0,-207.80,730.80,2.28,"1,257","16,599","2,109","3,366",16.86
D,Confirmed 1:1 change magnitude,CDS length,"19,965","1,262.12","1,251.26",-10.85,0,0,0,-0.52,985,"18,209",771,"1,756",8.80
D,Confirmed 1:1 change magnitude,Exon count,"19,965",5.86,5.83,-0.03,0,-1,1,-0.32,"1,534","17,368","1,063","2,597",13.01
D,Confirmed 1:1 change magnitude,mRNA count,"19,965",1,1,0,0,0,0,0,0,"19,965",0,0,0

