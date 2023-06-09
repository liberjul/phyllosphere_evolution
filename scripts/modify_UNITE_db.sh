REF=../data/reference_data/sh_general_release_dynamic_29.11.2022.fasta
REF_MOD=../data/reference_data/sh_general_release_dynamic_29.11.2022_mod.fasta

sed 's/o__Microbotryomycetes_ord_Incertae_sedis;f__Chrysozymaceae;g__Sampaiozyma/o__Microbotryales;f__Sampaiozymaceae;g__Sampaiozyma/' $REF > $REF_MOD
sed -i 's/o__Microbotryomycetes_ord_Incertae_sedis;f__Chrysozymaceae;g__Yunzhangia/o__Microbotryales;f__Microbotryaceae;g__Yunzhangia/' $REF_MOD
sed -i 's/o__Leucosporidiales;f__Leucosporidiaceae/o__Microbotryales;f__Leucosporidiaceae/' $REF_MOD
sed -i 's/f__Chrysozymaceae;g__Slooffia/f__Microbotryomycetes_fam_Incertae_sedis;g__Slooffia/' $REF_MOD
sed -i 's/o__Microbotryomycetes_ord_Incertae_sedis;f__Microbotryomycetes_fam_Incertae_sedis;g__Pseudoleucosporidium;/o__Sporidiobolales;f__Curvibasidiaceae;g__Pseudoleucosporidium;/' $REF_MOD
sed -i 's/o__Microbotryomycetes_ord_Incertae_sedis;f__Microbotryomycetes_fam_Incertae_sedis;g__Curvibasidium/o__Sporidiobolales;f__Curvibasidiaceae;g__Curvibasidium/' $REF_MOD
sed -i 's/o__Kriegeriales;f__Camptobasidiaceae/o__Microbotryomycetes_ord_Incertae_sedis;f__Camptobasidiaceae/' $REF_MOD
sed -i 's/f__Chrysozymaceae;g__Spencerozyma/f__Microbotryomycetes_fam_Incertae_sedis;g__Spencerozyma/' $REF_MOD
sed -i 's/f__Chrysozymaceae;g__Trigonosporomyces/f__Microbotryomycetes_fam_Incertae_sedis;g__Trigonosporomyces/' $REF_MOD
sed -i 's/f__Chrysozymaceae;g__Pseudohyphozyma/f__Microbotryomycetes_fam_Incertae_sedis;g__Pseudohyphozyma/' $REF_MOD
sed -i 's/f__Chrysozymaceae;g__Vonarxula/f__Microbotryomycetes_fam_Incertae_sedis;g__Vonarxula/' $REF_MOD
sed -i 's/f__Chrysozymaceae;g__Udeniozyma/f__Microbotryomycetes_fam_Incertae_sedis;g__Udeniozyma/' $REF_MOD
sed -i 's/f__Chrysozymaceae;g__Oberwinklerozyma/f__Microbotryomycetes_fam_Incertae_sedis;g__Oberwinklerozyma/' $REF_MOD

echo ">Aimania_erigeronia|OQ388269|SH0000000.01FU|refs|k__Fungi;p__Basidiomycota;c__Microbotryomycetes;o__Microbotryales;f__Sampaiozymaceae;g__Aimania;s__Aimania_erigeronia" >> $REF_MOD
echo "CATTAGTGAATATAGCGTGTCTTCGGAGCGCGACTCTCACTTTACACACTGTGCATTCTTTCAGCAGATGGACAAGACTTTCGGGTCGAGTCCCTGCGCTCTCATTAAACACGAGTTAATGTATGTGAAAATTACAAAACAAAGAAAAACTTTCAACAACGGATCTCTTGGCTCTCGCATCGATGAAGAACGCAGCGAATTGCGATAAGTAATGTGAATTGCAGAATTCAGTGAATCATCGAATCTTTGAACGCATCTTGCGCTCTCTGGTATTCCGGAGAGCATGTCTGTTTGAGTGTCATGAACTCTTCAACCCACCAGTTTCTAGTTAAACTGCGTGGCGCTTGGATCCTGAGCGCTTGCTCCGCTTCTTCCAGAGCTCGTTCGTAATACATTAGCTTTCTTGACTCGAATCGGATTGACTCGGCGTAATAGACTATTCGCTGAGGACAAGCTTCGCGCTTGGCCGGCTCGATTGTCAGATGAAGGCTTCTAGTCTCGTCAATTTTTAAGATTAGACCTCAAATCAGATAGGACTACCCGCTGAACTTAAGCATATCAATAAGCGGAGGAAGAGAAACTAACAAGGATTCCCCTAGTAACGGCGAGTGAAGCGGGAAAAGCTCAAATTTGAAATCTGGCACTTTCAGTGTCCGAGTTGTAATCTCGAGAAACGTTTTCCGCGCTGGACTGCATAGAAGTCTGTTGGAATACAGCGTCATAGTGGTGAGAACCCCGTACGTGATGCAGATGCCCAGTGCTTTGTGATACGTTTTCGAAGAGTCGAGTTGTTTGGGAATGCAGCTCAAAATGGGAGGTAAATTCCTTCTAAAGCTAAATATTGGCGAGAGACCGATAGCGAACAAGTACCGTGAGGGAAAGATGAAAAGCACTTTGGAAAGAGAGTTAACAGTACGTGAAATTGTTGGAAGGGAAACGCTTGAAGTCAAACTTGCTTGCCGGGCAACCGGTTTGCAGGCCAGCATCAGTTTCCAGGGGC" >> $REF_MOD

echo ">Aimania_cardamina|OQ388277.1|SH0000001.01FU|refs|k__Fungi;p__Basidiomycota;c__Microbotryomycetes;o__Microbotryales;f__Sampaiozymaceae;g__Aimania;s__Aimania_cardamina" >> $REF_MOD
echo "CATTAGTGAATTCAGCGAATCTCCGGATCGCGACCTCTCACTTTTCACACTGTGCACCTAGAATTCTTGGCGAGCCGCTAGGCGAACCCAAGATTTTTATTTAAACACGAGTTGATGTATGTCAATAAAATATAAAACAAAACAAAACTTTCAACAACGGATCTCTTGGCTCTCGCATCGATGAAGAACGCAGCGAAATGCGATAAGTAATGTGAATTGCAGAATTCAGTGAATCATCGAATCTTTGAACGCATCTTGCGCTCTCTGGTATTCCGGAGAGCATGTCTGTTTGAGTGTCATGAACTCTTCAACCCACCGGTTTCTAGTTAAACTGCCGTGGTGTTTGGATCTTGAGCGTCTGCTTTTCTTTTTGGAAAGCTCGTTCGTAATGCATTAGCTTTCTCGACTTTGATCGGATTGACTCGGCGTAATAGATTATTCGCTGAGGACAAGCTTCGCGCTTGGCCGACCATATGTCTACCGAAGGCTTCTAATCAGTCCTTTTGGACATACCTTTTAAGATTAGACCTCAAATCAGATAGGACTACCCGCTGAACTTAAGCATATCAATAAGCGGAGGAAAAGAAACTAACAAGGATTCCCCTAGTAACGGCGAGTGAAGCGGGAAAAGCTCAAATTTGAAATCTGGTACTTTCAGTACCCGAGTTGTAATCTCGAGAAACGTTTTCCGCGCCGGACTGCATAGAAGTCTGTTGGAATACAGCGTCATAGTGGTGAGAACCCCGTACGTGATGCAGATGCCCGGTGCTTTGTGATACGTTTTCGAAGAGTCGAGTTGTTTGGGAATGCAGCTCAAAATGGGAGGTAAATTCCTTCTAAAGCTAAATATTGGCGAGAGACCGATAGCGAACAAGTACCGTGAGGGAAAGATGAAAAGCACTTTGGAAAGAGAGTTAACAGTACGTGAAATTGTTGGAAGGGAAACGCTTGAAGTCAAACTTGCTTGCCGGGCAACCGGTTTGCAGGCCAGCATCAGTTTTCGGGGGTTGAAAAG" >> $REF_MOD
