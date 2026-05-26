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
B,Locus fate,Syntenic 1:1,"32,163",81.62,"32,163",79.76
B,Locus fate,Split,13,0.03,26,0.06
B,Locus fate,Merge,120,0.30,60,0.15
B,Locus fate,Complex,22,0.06,22,0.05
B,Locus fate,Unresolved weak-overlap,"3,135",7.96,"3,673",9.11
B,Locus fate,Strict unmatched,"3,953",10.03,"4,382",10.87


## Figure C table: 1:1 structural attributes

figure,layer,attribute,count,pct_of_syntenic
C,Confirmed 1:1 structural attributes,Exact 1:1,"17,711",55.07
C,Confirmed 1:1 structural attributes,Gene boundary changed,"10,394",32.32
C,Confirmed 1:1 structural attributes,UTR added,"5,156",16.03
C,Confirmed 1:1 structural attributes,UTR lost,65,0.20
C,Confirmed 1:1 structural attributes,UTR exon gained,"2,002",6.22
C,Confirmed 1:1 structural attributes,UTR exon removed,885,2.75
C,Confirmed 1:1 structural attributes,UTR refined,"6,565",20.41
C,Confirmed 1:1 structural attributes,Coding exon gain,702,2.18
C,Confirmed 1:1 structural attributes,Coding exon loss,"3,209",9.98
C,Confirmed 1:1 structural attributes,Exon boundary refined,"6,339",19.71
C,Confirmed 1:1 structural attributes,CDS changed,"5,846",18.18
C,Confirmed 1:1 structural attributes,CDS boundary refined,56,0.17
C,Confirmed 1:1 structural attributes,Isoform changed,"6,927",21.54


## Figure D table: paired change magnitude

figure,layer,metric,n_pairs,mean_before,mean_after,mean_delta,median_delta,p05_delta,p95_delta,mean_pct_delta,decreased,unchanged,increased,changed,changed_pct
D,Confirmed 1:1 change magnitude,Gene span length,"32,163","3,190.20","3,097.93",-92.27,0,-813.90,390,0.16,"6,301","21,514","4,348","10,649",33.11
D,Confirmed 1:1 change magnitude,Model span length,"32,163","3,190.07","3,097.93",-92.15,0,-812,390,0.16,"6,294","21,519","4,350","10,644",33.09
D,Confirmed 1:1 change magnitude,CDS length,"32,163","2,079.99","1,923.07",-156.92,0,"-1,428","1,050",2.71,"5,491","22,520","4,152","9,643",29.98
D,Confirmed 1:1 change magnitude,Exon count,"32,163",10.60,9.75,-0.84,0,-6,6,4.59,"4,208","23,870","4,085","8,293",25.78
D,Confirmed 1:1 change magnitude,mRNA count,"32,163",1.64,1.63,-0.01,0,-1,1,5.98,"2,128","26,849","3,186","5,314",16.52

