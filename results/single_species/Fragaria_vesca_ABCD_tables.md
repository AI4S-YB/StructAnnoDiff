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
B,Locus fate,Syntenic 1:1,"21,425",74.94,"21,425",63.00
B,Locus fate,Split,34,0.12,68,0.20
B,Locus fate,Merge,0,0,0,0
B,Locus fate,Complex,0,0,0,0
B,Locus fate,Unresolved weak-overlap,"3,691",12.91,"5,312",15.62
B,Locus fate,Strict unmatched,"3,438",12.03,"7,201",21.18


## Figure C table: 1:1 structural attributes

figure,layer,attribute,count,pct_of_syntenic
C,Confirmed 1:1 structural attributes,Exact 1:1,"2,918",13.62
C,Confirmed 1:1 structural attributes,Gene boundary changed,"17,473",81.55
C,Confirmed 1:1 structural attributes,UTR added,"3,617",16.88
C,Confirmed 1:1 structural attributes,UTR lost,"1,750",8.17
C,Confirmed 1:1 structural attributes,UTR exon gained,"1,062",4.96
C,Confirmed 1:1 structural attributes,UTR exon removed,"1,692",7.90
C,Confirmed 1:1 structural attributes,UTR refined,"13,088",61.09
C,Confirmed 1:1 structural attributes,Coding exon gain,"1,516",7.08
C,Confirmed 1:1 structural attributes,Coding exon loss,"2,382",11.12
C,Confirmed 1:1 structural attributes,Exon boundary refined,"13,037",60.85
C,Confirmed 1:1 structural attributes,CDS changed,"6,979",32.57
C,Confirmed 1:1 structural attributes,CDS boundary refined,179,0.84
C,Confirmed 1:1 structural attributes,Isoform changed,"8,797",41.06


## Figure D table: paired change magnitude

figure,layer,metric,n_pairs,mean_before,mean_after,mean_delta,median_delta,p05_delta,p95_delta,mean_pct_delta,decreased,unchanged,increased,changed,changed_pct
D,Confirmed 1:1 change magnitude,Gene span length,"21,425","3,370.27","3,408.83",38.56,53,"-1,094",918.80,5.63,"4,873","3,411","13,141","18,014",84.08
D,Confirmed 1:1 change magnitude,Model span length,"21,425","3,370.27","3,408.41",38.14,53,"-1,094",917.60,5.60,"4,874","3,411","13,140","18,014",84.08
D,Confirmed 1:1 change magnitude,CDS length,"21,425","1,314.21","3,342.29","2,028.08",0,-99,"9,760.80",119.64,"2,053","8,832","10,540","12,593",58.78
D,Confirmed 1:1 change magnitude,Exon count,"21,425",5.93,16.90,10.97,0,-1,58,111.12,"1,791","10,013","9,621","11,412",53.26
D,Confirmed 1:1 change magnitude,mRNA count,"21,425",1,2.26,1.26,0,0,6,126.22,0,"12,628","8,797","8,797",41.06

