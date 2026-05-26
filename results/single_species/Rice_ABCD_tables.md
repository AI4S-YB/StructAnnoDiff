# Rice A/B/C/D analysis tables

## Figure A table: global quantities

figure,layer,metric,metric_type,before,after,delta,pct_change
A,Global annotation quantity,Genes,count,"39,406","40,326",920,2.33
A,Global annotation quantity,mRNAs,count,"39,406","40,331",925,2.35
A,Global annotation quantity,Exons,count,"180,582","182,652","2,070",1.15
A,Global annotation quantity,CDS features,count,"39,406","40,140",734,1.86
A,Global annotation quantity,5' UTR features,count,"16,156","26,844","10,688",66.15
A,Global annotation quantity,3' UTR features,count,"18,235","27,963","9,728",53.35
A,Global annotation quantity,Single-exon genes,count,"8,888","9,362",474,5.33
A,Global annotation quantity,Mean exons per mRNA,mean,4.60,4.50,-0.10,-2.17
A,Global annotation quantity,Mean CDS length (bp),mean,"1,093","1,064",-29,-2.65
A,Global annotation quantity,Mean gene length (bp),mean,"3,132","3,104",-28,-0.89
A,Global annotation quantity,Median gene length (bp),median,"2,046","2,357",311,15.20


## Figure B table: locus fate

figure,layer,category,before_genes,before_pct,after_genes,after_pct
B,Locus fate,Syntenic 1:1,"32,478",82.42,"32,478",80.54
B,Locus fate,Split,746,1.89,"1,629",4.04
B,Locus fate,Merge,"1,019",2.59,498,1.23
B,Locus fate,Complex,987,2.50,"1,062",2.63
B,Locus fate,Unresolved weak-overlap,0,0,0,0
B,Locus fate,Strict unmatched,"4,176",10.60,"4,658",11.55


## Figure C table: 1:1 structural attributes

figure,layer,attribute,count,pct_of_syntenic
C,Confirmed 1:1 structural attributes,Exact 1:1,"17,415",53.62
C,Confirmed 1:1 structural attributes,Gene boundary changed,"11,144",34.31
C,Confirmed 1:1 structural attributes,UTR added,"5,872",18.08
C,Confirmed 1:1 structural attributes,UTR lost,64,0.20
C,Confirmed 1:1 structural attributes,UTR exon gained,"2,263",6.97
C,Confirmed 1:1 structural attributes,UTR exon removed,880,2.71
C,Confirmed 1:1 structural attributes,UTR refined,"6,541",20.14
C,Confirmed 1:1 structural attributes,Coding exon gain,785,2.42
C,Confirmed 1:1 structural attributes,Coding exon loss,"3,413",10.51
C,Confirmed 1:1 structural attributes,Exon boundary refined,"6,538",20.13
C,Confirmed 1:1 structural attributes,CDS changed,"6,213",19.13
C,Confirmed 1:1 structural attributes,CDS boundary refined,34,0.10
C,Confirmed 1:1 structural attributes,Isoform changed,"6,863",21.13


## Figure D table: paired change magnitude

figure,layer,metric,n_pairs,mean_before,mean_after,mean_delta,median_delta,p05_delta,p95_delta,mean_pct_delta,decreased,unchanged,increased,changed,changed_pct
D,Confirmed 1:1 change magnitude,Gene span length,"32,478","3,075.03","3,023.58",-51.45,0,-884.15,576.15,11.39,"6,378","21,082","5,018","11,396",35.09
D,Confirmed 1:1 change magnitude,Model span length,"32,478","3,074.93","3,023.58",-51.35,0,-883.15,576.15,11.39,"6,371","21,087","5,020","11,391",35.07
D,Confirmed 1:1 change magnitude,CDS length,"32,478","1,954.58","1,842.84",-111.73,0,"-1,287.45","1,053",4.75,"5,605","22,434","4,439","10,044",30.93
D,Confirmed 1:1 change magnitude,Exon count,"32,478",9.77,9.29,-0.48,0,-5,6,8.32,"4,304","23,769","4,405","8,709",26.82
D,Confirmed 1:1 change magnitude,mRNA count,"32,478",1.60,1.61,0.01,0,-1,1,6.70,"2,005","27,186","3,287","5,292",16.29

