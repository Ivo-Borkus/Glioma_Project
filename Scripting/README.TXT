Questions:
Can we highlight composition specific differences between oligo and astrocytomas?
    This might be at a deeper level than what we are seeing now. --> BAM subtypes
Using this, can we figure out what is the population or cell state driving the disease subtypes? 

Can we identify differences between grades in comparing astrocytoma and GBM?
Are there differences between glioblastoma recurrent and primary tumours.



To do:
Get rid of the bad quality populations. --> Redo processing
    Perform 'naive' filtering of bad quality cells

Signature scoring (pro/anti-inf and resident like signature)
compute the overall signatures only on myeloid cells 
    compare the three diseases 
Stratify further into microglia, macrophages and monocytes to compare (overall) signatures.
Pseudobulk myeloid cells per patient --> Score them on pro/anti-inf and resident-infiltrating signatures.
Test significance per model.

Differential abundance testing (MiloR)
First: Split lineage (Myeloid and T cells)
--> Compare abundances between models of interest: 
1. Astrocytoma ~ oligodendroglioma
2. GBM recurrent ~ GBM primary 
3. GBM ~ Astrocytoma (grade 1) + (Grade 2)... --> Can we see something with a trend?

Confirm findings with scCODA



