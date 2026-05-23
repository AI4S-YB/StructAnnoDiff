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
B,Locus fate,Syntenic 1:1,"28,986",96.04,"28,986",91.29
B,Locus fate,Split,322,1.07,679,2.14
B,Locus fate,Merge,0,0,0,0
B,Locus fate,Complex,0,0,0,0
B,Locus fate,Unresolved weak-overlap,784,2.60,"1,294",4.08
B,Locus fate,Strict deleted / novel,89,0.29,792,2.49


## Figure C table: 1:1 structural attributes

figure,layer,attribute,count,pct_of_syntenic
C,Confirmed 1:1 structural attributes,Exact 1:1,"26,562",91.64
C,Confirmed 1:1 structural attributes,Gene boundary changed,"1,695",5.85
C,Confirmed 1:1 structural attributes,UTR added,"1,247",4.30
C,Confirmed 1:1 structural attributes,UTR lost,42,0.14
C,Confirmed 1:1 structural attributes,UTR exon gained,110,0.38
C,Confirmed 1:1 structural attributes,UTR exon removed,556,1.92
C,Confirmed 1:1 structural attributes,UTR refined,802,2.77
C,Confirmed 1:1 structural attributes,Coding exon gain,438,1.51
C,Confirmed 1:1 structural attributes,Coding exon loss,"1,112",3.84
C,Confirmed 1:1 structural attributes,Exon boundary refined,387,1.34
C,Confirmed 1:1 structural attributes,CDS changed,"1,982",6.84
C,Confirmed 1:1 structural attributes,CDS boundary refined,259,0.89
C,Confirmed 1:1 structural attributes,Isoform changed,567,1.96


## Figure D table: paired change magnitude

figure,layer,metric,n_pairs,mean_before,mean_after,mean_delta,median_delta,p05_delta,p95_delta,mean_pct_delta,decreased,unchanged,increased,changed,changed_pct
D,Confirmed 1:1 change magnitude,Gene span length,"28,986","2,916.14","2,902.58",-13.57,0,0,0,1.03,510,"27,285","1,191","1,701",5.87
D,Confirmed 1:1 change magnitude,Model span length,"28,986","2,916.14","2,902.58",-13.57,0,0,0,1.03,510,"27,285","1,191","1,701",5.87
D,Confirmed 1:1 change magnitude,CDS length,"28,986","1,457.60","1,426.98",-30.62,0,0,0,0.34,752,"26,904","1,330","2,082",7.18
D,Confirmed 1:1 change magnitude,Exon count,"28,986",6.32,5.99,-0.34,0,-1,0,-2.66,"1,627","27,248",111,"1,738",6.00
D,Confirmed 1:1 change magnitude,mRNA count,"28,986",1.16,1.14,-0.02,0,0,0,-0.92,567,"28,419",0,567,1.96

