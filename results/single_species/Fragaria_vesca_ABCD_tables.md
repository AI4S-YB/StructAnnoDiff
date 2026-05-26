# F. vesca A/B/C/D analysis tables

## Figure A table: global quantities

figure,layer,metric,metric_type,before,after,delta,pct_change
A,Global annotation quantity,Genes,count,"28,588","34,006","5,418",18.95
A,Global annotation quantity,mRNAs,count,"28,588","34,007","5,419",18.96
A,Global annotation quantity,Exons,count,"156,498","164,975","8,477",5.42
A,Global annotation quantity,CDS features,count,"28,588","34,006","5,418",18.95
A,Global annotation quantity,5' UTR features,count,"16,001","19,239","3,238",20.24
A,Global annotation quantity,3' UTR features,count,"17,263","19,711","2,448",14.18
A,Global annotation quantity,Single-exon genes,count,"5,004","5,481",477,9.53
A,Global annotation quantity,Mean exons per mRNA,mean,5.50,4.90,-0.60,-10.91
A,Global annotation quantity,Mean CDS length (bp),mean,"1,178","1,155",-23,-1.95
A,Global annotation quantity,Mean gene length (bp),mean,"3,214","2,875",-339,-10.55
A,Global annotation quantity,Median gene length (bp),median,"2,406","2,253",-153,-6.36


## Figure B table: locus fate

figure,layer,category,before_genes,before_pct,after_genes,after_pct
B,Locus fate,Syntenic 1:1,"22,873",80.01,"22,873",67.26
B,Locus fate,Split,"1,441",5.04,"3,128",9.20
B,Locus fate,Merge,547,1.91,267,0.79
B,Locus fate,Complex,48,0.17,52,0.15
B,Locus fate,Unresolved weak-overlap,0,0,0,0
B,Locus fate,Strict unmatched,"3,679",12.87,"7,686",22.60


## Figure C table: 1:1 structural attributes

figure,layer,attribute,count,pct_of_syntenic
C,Confirmed 1:1 structural attributes,Exact 1:1,"2,835",12.39
C,Confirmed 1:1 structural attributes,Gene boundary changed,"19,018",83.15
C,Confirmed 1:1 structural attributes,UTR added,"4,436",19.39
C,Confirmed 1:1 structural attributes,UTR lost,"2,064",9.02
C,Confirmed 1:1 structural attributes,UTR exon gained,"1,238",5.41
C,Confirmed 1:1 structural attributes,UTR exon removed,"1,793",7.84
C,Confirmed 1:1 structural attributes,UTR refined,"13,050",57.05
C,Confirmed 1:1 structural attributes,Coding exon gain,"1,875",8.20
C,Confirmed 1:1 structural attributes,Coding exon loss,"2,508",10.96
C,Confirmed 1:1 structural attributes,Exon boundary refined,"13,802",60.34
C,Confirmed 1:1 structural attributes,CDS changed,"7,569",33.09
C,Confirmed 1:1 structural attributes,CDS boundary refined,51,0.22
C,Confirmed 1:1 structural attributes,Isoform changed,"8,896",38.89


## Figure D table: paired change magnitude

figure,layer,metric,n_pairs,mean_before,mean_after,mean_delta,median_delta,p05_delta,p95_delta,mean_pct_delta,decreased,unchanged,increased,changed,changed_pct
D,Confirmed 1:1 change magnitude,Gene span length,"22,873","3,175.40","3,254.54",79.14,66,"-1,388","1,366.40",20.30,"5,190","3,314","14,369","19,559",85.51
D,Confirmed 1:1 change magnitude,Model span length,"22,873","3,175.40","3,254.15",78.75,66,"-1,388","1,364.20",20.27,"5,191","3,314","14,368","19,559",85.51
D,Confirmed 1:1 change magnitude,CDS length,"22,873","1,231.52","3,112.10","1,880.58",0,-114,"9,126",115.84,"2,328","9,318","11,227","13,555",59.26
D,Confirmed 1:1 change magnitude,Exon count,"22,873",5.53,15.59,10.06,0,-1,54,109.01,"2,192","10,482","10,199","12,391",54.17
D,Confirmed 1:1 change magnitude,mRNA count,"22,873",1,2.18,1.18,0,0,6,117.89,0,"13,977","8,896","8,896",38.89

