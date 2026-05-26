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
B,Locus fate,Syntenic 1:1,"20,523",84.40,"20,523",85.00
B,Locus fate,Split,274,1.13,559,2.32
B,Locus fate,Merge,"1,208",4.97,591,2.45
B,Locus fate,Complex,24,0.10,24,0.10
B,Locus fate,Unresolved weak-overlap,0,0,0,0
B,Locus fate,Strict unmatched,"2,288",9.41,"2,448",10.14


## Figure C table: 1:1 structural attributes

figure,layer,attribute,count,pct_of_syntenic
C,Confirmed 1:1 structural attributes,Exact 1:1,"16,381",79.82
C,Confirmed 1:1 structural attributes,Gene boundary changed,"3,931",19.15
C,Confirmed 1:1 structural attributes,UTR added,"2,381",11.60
C,Confirmed 1:1 structural attributes,UTR lost,4,0.02
C,Confirmed 1:1 structural attributes,UTR exon gained,"1,101",5.36
C,Confirmed 1:1 structural attributes,UTR exon removed,808,3.94
C,Confirmed 1:1 structural attributes,UTR refined,"1,934",9.42
C,Confirmed 1:1 structural attributes,Coding exon gain,360,1.75
C,Confirmed 1:1 structural attributes,Coding exon loss,"1,390",6.77
C,Confirmed 1:1 structural attributes,Exon boundary refined,"1,102",5.37
C,Confirmed 1:1 structural attributes,CDS changed,"2,314",11.28
C,Confirmed 1:1 structural attributes,CDS boundary refined,49,0.24
C,Confirmed 1:1 structural attributes,Isoform changed,0,0


## Figure D table: paired change magnitude

figure,layer,metric,n_pairs,mean_before,mean_after,mean_delta,median_delta,p05_delta,p95_delta,mean_pct_delta,decreased,unchanged,increased,changed,changed_pct
D,Confirmed 1:1 change magnitude,Gene span length,"20,523","4,410.18","4,291.24",-118.94,0,"-1,013",693.90,10.07,"2,940","15,545","2,038","4,978",24.26
D,Confirmed 1:1 change magnitude,Model span length,"20,523","4,332.07","4,290.98",-41.09,0,-647.70,762.80,11.76,"1,620","16,560","2,343","3,963",19.31
D,Confirmed 1:1 change magnitude,CDS length,"20,523","1,213.77","1,206.57",-7.20,0,-39,0.90,0.08,"1,238","18,258","1,027","2,265",11.04
D,Confirmed 1:1 change magnitude,Exon count,"20,523",5.61,5.57,-0.05,0,-1,0,-0.08,"1,836","17,689",998,"2,834",13.81
D,Confirmed 1:1 change magnitude,mRNA count,"20,523",1,1,0,0,0,0,0,0,"20,523",0,0,0

