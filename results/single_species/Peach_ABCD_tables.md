# Peach A/B/C/D analysis tables

## Figure A table: global quantities

figure,layer,metric,metric_type,before,after,delta,pct_change
A,Global annotation quantity,Genes,count,"30,181","31,751","1,570",5.20
A,Global annotation quantity,mRNAs,count,"30,181","31,751","1,570",5.20
A,Global annotation quantity,Exons,count,"152,350","156,559","4,209",2.76
A,Global annotation quantity,CDS features,count,"30,181","31,749","1,568",5.20
A,Global annotation quantity,5' UTR features,count,"13,535","16,309","2,774",20.50
A,Global annotation quantity,3' UTR features,count,"13,967","16,826","2,859",20.47
A,Global annotation quantity,Single-exon genes,count,"4,212","5,702","1,490",35.38
A,Global annotation quantity,Mean exons per mRNA,mean,5,4.90,-0.10,-2
A,Global annotation quantity,Mean CDS length (bp),mean,"1,211","1,220",9,0.74
A,Global annotation quantity,Mean gene length (bp),mean,"2,921","2,882",-39,-1.34
A,Global annotation quantity,Median gene length (bp),median,"2,178","2,224",46,2.11


## Figure B table: locus fate

figure,layer,category,before_genes,before_pct,after_genes,after_pct
B,Locus fate,Syntenic 1:1,"29,305",97.10,"29,305",92.30
B,Locus fate,Split,678,2.25,"1,454",4.58
B,Locus fate,Merge,66,0.22,33,0.10
B,Locus fate,Complex,6,0.02,6,0.02
B,Locus fate,Unresolved weak-overlap,0,0,0,0
B,Locus fate,Strict unmatched,126,0.42,952,3.00


## Figure C table: 1:1 structural attributes

figure,layer,attribute,count,pct_of_syntenic
C,Confirmed 1:1 structural attributes,Exact 1:1,"26,561",90.64
C,Confirmed 1:1 structural attributes,Gene boundary changed,"2,029",6.92
C,Confirmed 1:1 structural attributes,UTR added,"1,612",5.50
C,Confirmed 1:1 structural attributes,UTR lost,51,0.17
C,Confirmed 1:1 structural attributes,UTR exon gained,183,0.62
C,Confirmed 1:1 structural attributes,UTR exon removed,472,1.61
C,Confirmed 1:1 structural attributes,UTR refined,709,2.42
C,Confirmed 1:1 structural attributes,Coding exon gain,493,1.68
C,Confirmed 1:1 structural attributes,Coding exon loss,"1,304",4.45
C,Confirmed 1:1 structural attributes,Exon boundary refined,457,1.56
C,Confirmed 1:1 structural attributes,CDS changed,"2,041",6.96
C,Confirmed 1:1 structural attributes,CDS boundary refined,19,0.06
C,Confirmed 1:1 structural attributes,Isoform changed,541,1.85


## Figure D table: paired change magnitude

figure,layer,metric,n_pairs,mean_before,mean_after,mean_delta,median_delta,p05_delta,p95_delta,mean_pct_delta,decreased,unchanged,increased,changed,changed_pct
D,Confirmed 1:1 change magnitude,Gene span length,"29,305","2,881.29","2,871.95",-9.34,0,0,64.80,4.38,507,"27,270","1,528","2,035",6.94
D,Confirmed 1:1 change magnitude,Model span length,"29,305","2,881.29","2,871.95",-9.34,0,0,64.80,4.38,507,"27,270","1,528","2,035",6.94
D,Confirmed 1:1 change magnitude,CDS length,"29,305","1,436.98","1,415.47",-21.51,0,0,30,0.98,774,"26,925","1,606","2,380",8.12
D,Confirmed 1:1 change magnitude,Exon count,"29,305",6.19,5.90,-0.29,0,-1,0,-2.35,"1,754","27,318",233,"1,987",6.78
D,Confirmed 1:1 change magnitude,mRNA count,"29,305",1.16,1.14,-0.02,0,0,0,-0.86,541,"28,764",0,541,1.85

