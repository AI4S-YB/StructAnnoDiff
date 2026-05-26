# Pineapple A/B/C/D analysis tables

## Figure A table: global quantities

figure,layer,metric,metric_type,before,after,delta,pct_change
A,Global annotation quantity,Genes,count,"26,162","26,657",495,1.89
A,Global annotation quantity,mRNAs,count,"26,162","26,657",495,1.89
A,Global annotation quantity,Exons,count,"131,045","137,990","6,945",5.30
A,Global annotation quantity,CDS features,count,"26,162","26,582",420,1.61
A,Global annotation quantity,5' UTR features,count,"11,845","16,668","4,823",40.72
A,Global annotation quantity,3' UTR features,count,"12,371","17,182","4,811",38.89
A,Global annotation quantity,Single-exon genes,count,"5,603","5,183",-420,-7.50
A,Global annotation quantity,Mean exons per mRNA,mean,5,5.20,0.20,4
A,Global annotation quantity,Mean CDS length (bp),mean,"1,153","1,146",-7,-0.61
A,Global annotation quantity,Mean gene length (bp),mean,"4,244","4,392",148,3.49
A,Global annotation quantity,Median gene length (bp),median,"2,894","2,986",92,3.18


## Figure B table: locus fate

figure,layer,category,before_genes,before_pct,after_genes,after_pct
B,Locus fate,Syntenic 1:1,"24,639",94.18,"24,639",92.43
B,Locus fate,Split,16,0.06,32,0.12
B,Locus fate,Merge,2,0.01,1,0.00
B,Locus fate,Complex,4,0.02,4,0.02
B,Locus fate,Unresolved weak-overlap,"1,103",4.22,"1,051",3.94
B,Locus fate,Strict unmatched,398,1.52,930,3.49


## Figure C table: 1:1 structural attributes

figure,layer,attribute,count,pct_of_syntenic
C,Confirmed 1:1 structural attributes,Exact 1:1,"15,053",61.09
C,Confirmed 1:1 structural attributes,Gene boundary changed,"6,995",28.39
C,Confirmed 1:1 structural attributes,UTR added,"3,850",15.63
C,Confirmed 1:1 structural attributes,UTR lost,353,1.43
C,Confirmed 1:1 structural attributes,UTR exon gained,"3,865",15.69
C,Confirmed 1:1 structural attributes,UTR exon removed,0,0
C,Confirmed 1:1 structural attributes,UTR refined,"5,507",22.35
C,Confirmed 1:1 structural attributes,Coding exon gain,180,0.73
C,Confirmed 1:1 structural attributes,Coding exon loss,806,3.27
C,Confirmed 1:1 structural attributes,Exon boundary refined,"5,120",20.78
C,Confirmed 1:1 structural attributes,CDS changed,"1,912",7.76
C,Confirmed 1:1 structural attributes,CDS boundary refined,718,2.91
C,Confirmed 1:1 structural attributes,Isoform changed,0,0


## Figure D table: paired change magnitude

figure,layer,metric,n_pairs,mean_before,mean_after,mean_delta,median_delta,p05_delta,p95_delta,mean_pct_delta,decreased,unchanged,increased,changed,changed_pct
D,Confirmed 1:1 change magnitude,Gene span length,"24,639","4,304.30","4,284.21",-20.09,0,-413.10,431.10,1.09,"3,476","17,526","3,637","7,113",28.87
D,Confirmed 1:1 change magnitude,Model span length,"24,639","4,285.30","4,282.04",-3.26,0,-384.10,432,1.16,"3,533","17,457","3,649","7,182",29.15
D,Confirmed 1:1 change magnitude,CDS length,"24,639","1,175.22","1,167.01",-8.20,0,0,0,-0.13,882,"23,445",312,"1,194",4.85
D,Confirmed 1:1 change magnitude,Exon count,"24,639",5.09,5.23,0.14,0,0,1,5.23,650,"20,366","3,623","4,273",17.34
D,Confirmed 1:1 change magnitude,mRNA count,"24,639",1,1,0,0,0,0,0,0,"24,639",0,0,0

