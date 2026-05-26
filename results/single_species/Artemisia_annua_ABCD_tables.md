# A. annua A/B/C/D analysis tables

## Figure A table: global quantities

figure,layer,metric,metric_type,before,after,delta,pct_change
A,Global annotation quantity,Genes,count,"54,347","39,322","-15,025",-27.65
A,Global annotation quantity,mRNAs,count,"54,347","39,322","-15,025",-27.65
A,Global annotation quantity,Exons,count,"226,160","203,145","-23,015",-10.18
A,Global annotation quantity,CDS features,count,"54,347","39,322","-15,025",-27.65
A,Global annotation quantity,5' UTR features,count,"45,872","34,818","-11,054",-24.10
A,Global annotation quantity,3' UTR features,count,"46,516","35,459","-11,057",-23.77
A,Global annotation quantity,Single-exon genes,count,"17,013","7,932","-9,081",-53.38
A,Global annotation quantity,Mean exons per mRNA,mean,4.20,5.20,1,23.81
A,Global annotation quantity,Mean CDS length (bp),mean,975,"1,253",278,28.51
A,Global annotation quantity,Mean gene length (bp),mean,"3,485","4,240",755,21.66
A,Global annotation quantity,Median gene length (bp),median,"2,676","3,284",608,22.72


## Figure B table: locus fate

figure,layer,category,before_genes,before_pct,after_genes,after_pct
B,Locus fate,Syntenic 1:1,"34,444",63.38,"34,444",87.59
B,Locus fate,Split,28,0.05,56,0.14
B,Locus fate,Merge,"2,805",5.16,"1,325",3.37
B,Locus fate,Complex,81,0.15,76,0.19
B,Locus fate,Unresolved weak-overlap,0,0,0,0
B,Locus fate,Strict unmatched,"16,989",31.26,"3,421",8.70


## Figure C table: 1:1 structural attributes

figure,layer,attribute,count,pct_of_syntenic
C,Confirmed 1:1 structural attributes,Exact 1:1,"18,932",54.96
C,Confirmed 1:1 structural attributes,Gene boundary changed,"12,225",35.49
C,Confirmed 1:1 structural attributes,UTR added,"1,940",5.63
C,Confirmed 1:1 structural attributes,UTR lost,401,1.16
C,Confirmed 1:1 structural attributes,UTR exon gained,"1,824",5.30
C,Confirmed 1:1 structural attributes,UTR exon removed,"1,054",3.06
C,Confirmed 1:1 structural attributes,UTR refined,"12,487",36.25
C,Confirmed 1:1 structural attributes,Coding exon gain,"2,146",6.23
C,Confirmed 1:1 structural attributes,Coding exon loss,467,1.36
C,Confirmed 1:1 structural attributes,Exon boundary refined,"9,412",27.33
C,Confirmed 1:1 structural attributes,CDS changed,"3,836",11.14
C,Confirmed 1:1 structural attributes,CDS boundary refined,182,0.53
C,Confirmed 1:1 structural attributes,Isoform changed,"4,249",12.34


## Figure D table: paired change magnitude

figure,layer,metric,n_pairs,mean_before,mean_after,mean_delta,median_delta,p05_delta,p95_delta,mean_pct_delta,decreased,unchanged,increased,changed,changed_pct
D,Confirmed 1:1 change magnitude,Gene span length,"34,444","3,993.50","4,231.96",238.46,0,-608,"1,585.70",18.14,"4,065","20,604","9,775","13,840",40.18
D,Confirmed 1:1 change magnitude,Model span length,"34,444","3,993.01","4,231.96",238.95,0,-603,"1,586.85",18.19,"4,032","20,626","9,786","13,818",40.12
D,Confirmed 1:1 change magnitude,CDS length,"34,444","1,255.16","1,510.40",255.24,0,0,"1,755",18.78,"1,106","26,820","6,518","7,624",22.13
D,Confirmed 1:1 change magnitude,Exon count,"34,444",5.10,6.45,1.36,0,0,9,22.63,754,"26,610","7,080","7,834",22.74
D,Confirmed 1:1 change magnitude,mRNA count,"34,444",1,1.15,0.15,0,0,1,15.20,0,"30,195","4,249","4,249",12.34

