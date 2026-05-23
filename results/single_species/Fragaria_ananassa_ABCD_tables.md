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
B,Locus fate,Syntenic 1:1,"93,658",86.65,"93,658",86.36
B,Locus fate,Split,747,0.69,"1,496",1.38
B,Locus fate,Merge,0,0,0,0
B,Locus fate,Complex,0,0,0,0
B,Locus fate,Unresolved weak-overlap,"7,322",6.77,"10,825",9.98
B,Locus fate,Strict deleted / novel,"6,360",5.88,"2,468",2.28


## Figure C table: 1:1 structural attributes

figure,layer,attribute,count,pct_of_syntenic
C,Confirmed 1:1 structural attributes,Exact 1:1,"50,847",54.29
C,Confirmed 1:1 structural attributes,Gene boundary changed,"40,143",42.86
C,Confirmed 1:1 structural attributes,UTR added,"10,402",11.11
C,Confirmed 1:1 structural attributes,UTR lost,882,0.94
C,Confirmed 1:1 structural attributes,UTR exon gained,"6,239",6.66
C,Confirmed 1:1 structural attributes,UTR exon removed,"4,215",4.50
C,Confirmed 1:1 structural attributes,UTR refined,"34,093",36.40
C,Confirmed 1:1 structural attributes,Coding exon gain,"3,881",4.14
C,Confirmed 1:1 structural attributes,Coding exon loss,"9,697",10.35
C,Confirmed 1:1 structural attributes,Exon boundary refined,"25,216",26.92
C,Confirmed 1:1 structural attributes,CDS changed,"20,402",21.78
C,Confirmed 1:1 structural attributes,CDS boundary refined,723,0.77
C,Confirmed 1:1 structural attributes,Isoform changed,"8,812",9.41


## Figure D table: paired change magnitude

figure,layer,metric,n_pairs,mean_before,mean_after,mean_delta,median_delta,p05_delta,p95_delta,mean_pct_delta,decreased,unchanged,increased,changed,changed_pct
D,Confirmed 1:1 change magnitude,Gene span length,"93,658","3,067.84","3,009.90",-57.94,0,-758.15,505,1.15,"20,438","51,618","21,602","42,040",44.89
D,Confirmed 1:1 change magnitude,Model span length,"93,658","3,067.84","3,005.11",-62.72,0,-775,499,1.04,"20,584","51,620","21,454","42,038",44.88
D,Confirmed 1:1 change magnitude,CDS length,"93,658","1,142.35","1,269.53",127.19,0,-171,"1,086",10.92,"11,264","70,458","11,936","23,200",24.77
D,Confirmed 1:1 change magnitude,Exon count,"93,658",5.42,6.78,1.36,0,-1,8,14.08,"7,386","74,138","12,134","19,520",20.84
D,Confirmed 1:1 change magnitude,mRNA count,"93,658",1,1.16,0.16,0,0,1,16.45,0,"84,846","8,812","8,812",9.41

