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
B,Locus fate,Syntenic 1:1,"24,870",95.06,"24,870",93.41
B,Locus fate,Split,205,0.78,421,1.58
B,Locus fate,Merge,589,2.25,288,1.08
B,Locus fate,Complex,47,0.18,53,0.20
B,Locus fate,Unresolved weak-overlap,0,0,0,0
B,Locus fate,Strict unmatched,451,1.72,993,3.73


## Figure C table: 1:1 structural attributes

figure,layer,attribute,count,pct_of_syntenic
C,Confirmed 1:1 structural attributes,Exact 1:1,"8,759",35.22
C,Confirmed 1:1 structural attributes,Gene boundary changed,"16,041",64.50
C,Confirmed 1:1 structural attributes,UTR added,"16,069",64.61
C,Confirmed 1:1 structural attributes,UTR lost,0,0
C,Confirmed 1:1 structural attributes,UTR exon gained,"3,930",15.80
C,Confirmed 1:1 structural attributes,UTR exon removed,0,0
C,Confirmed 1:1 structural attributes,UTR refined,0,0
C,Confirmed 1:1 structural attributes,Coding exon gain,232,0.93
C,Confirmed 1:1 structural attributes,Coding exon loss,843,3.39
C,Confirmed 1:1 structural attributes,Exon boundary refined,"11,597",46.63
C,Confirmed 1:1 structural attributes,CDS changed,"1,324",5.32
C,Confirmed 1:1 structural attributes,CDS boundary refined,15,0.06
C,Confirmed 1:1 structural attributes,Isoform changed,0,0


## Figure D table: paired change magnitude

figure,layer,metric,n_pairs,mean_before,mean_after,mean_delta,median_delta,p05_delta,p95_delta,mean_pct_delta,decreased,unchanged,increased,changed,changed_pct
D,Confirmed 1:1 change magnitude,Gene span length,"24,870","4,273.21","4,263.85",-9.36,0,-551,460,6.81,"3,660","17,491","3,719","7,379",29.67
D,Confirmed 1:1 change magnitude,Model span length,"24,870","3,737.49","4,259.57",522.08,277,0,"1,929.55",28.97,453,"8,790","15,627","16,080",64.66
D,Confirmed 1:1 change magnitude,CDS length,"24,870","1,158.17","1,153.62",-4.55,0,0,0,-0.00,934,"23,561",375,"1,309",5.26
D,Confirmed 1:1 change magnitude,Exon count,"24,870",5.00,5.16,0.16,0,0,1,5.92,686,"20,460","3,724","4,410",17.73
D,Confirmed 1:1 change magnitude,mRNA count,"24,870",1,1,0,0,0,0,0,0,"24,870",0,0,0

