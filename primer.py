from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
import primer3
from Bio.SeqUtils import GC
from orffinder import orffinder

seq = 'cgcaaatgggcggtaggcgtgtacggtgggaggtctatataagcagagctggtttagtgaaccgtcagatccgctagcgattacgccaagctcgaaattaaccctcactaaagggaacaaaagctggagctccaccgcggtggcggccgcatggagcaaacagtgcttgtaccaccaggacctgacagcttcaacttcttcaccagagaatctcttgcggctattgaaagacgcattgcagaagaaaaggcaaagaatcccaaaccagacaaaaaagatgacgacgaaaatggcccaaagccaaatagtgacttggaagctggaaagaaccttccatttatttatggagacattcctccagagatggtgtcagagcccctggaggacctggacccctactatatcaataagaaaacttttatagtattgaataaagggaaggccatcttccggttcagtgccacctctgccctgtacattttaactcccttcaatcctcttaggaaaatagctattaagattttggtacattcattattcagcatgctaattatgtgcactattttgacaaactgtgtgtttatgacaatgagtaaccctcctgattggacaaagaatgtagaatacaccttcacaggaatatatacttttgaatcacttataaaaattattgcaaggggattctgtttagaagattttactttccttcgggatccatggaactggctcgatttcactgtcattacatttgcgtacgtcacagagtttgtggacctgggcaatgtctcggcattgagaacattcagagttctccgagcattgaagacgatttcagtcattccaggcctgaaaaccattgtgggagccctgatccagtctgtgaagaagctctcagatgtaatgatcctgactgtgttctgtctgagcgtatttgctctaattgggctgcagctgttcatgggcaacctgaggaataaatgtatacaatggcctcccaccaatgcttccttggaggaacatagtatagaaaagaatataactgtgaactacaacggtacacttataaatgaaactgtctttgagtttgactggaagtcatatattcaagattcaagatatcattatttcctggagggttttttagatgcactactatgtggaaatagctctgatgcaggccaatgtccagagggatatatgtgtgtgaaagctggtagaaatcccaattatggctacacaagctttgataccttcagttgggcttttttgtccttgtttcgactaatgactcaggacttctgggaaaatctttatcaactgacattacgtgctgctgggaaaacgtacatgatattcttcgtattggtcattttcttgggctcattctacctaataaatttgatcctggctgtggtggccatggcctacgaggaacagaatcaggccaccttggaagaagcagaacagaaagaggccgaatttcagcagatgattgaacagcttaaaaagcaacaggaggcagctcagcaggcagcaacggcaactgcctcagaacattccagagagcccagtgcagcaggcaggctctcagacagctcatctgaagcctctaagttgagttccaagagtgctaaggaaagaagaaatcggaggaagaaaagaaaacagaaagagcagtctggtggggaagagaaagatgaggatgaattccaaaaatctgaatctgaggacagcatcaggaggaaaggttttcgcttctccattgaagggaaccgattgacatatgaaaagaggtactcctccccacaccagtctttgttgagcatccgtggctccctattttcaccaaggcgaaatagcagaacaagccttttcagctttagagggcgagcaaaggatgtgggatctgagaacgacttcgcagatgatgagcacagcacctttgaggataacgagagccgtagagattccttgtttgtgccccgacgacacggagagagacgcaacagcaacctgagtcagaccagtaggtcatcccggatgctggcagtgtttccagcgaatgggaagatgcacagcactgtggattgcaatggtgtggtttccttggttggtggaccttcagttcctacatcgcctgttggacagcttctgccagaggtgataatagataagccagctactgatgacaatggaacaaccactgaaactgaaatgagaaagagaaggtcaagttctttccacgtttccatggactttctagaagatccttcccaaaggcaacgagcaatgagtatagccagcattctaacaaatacagtagaagaacttgaagaatccaggcagaaatgcccaccctgttggtataaattttccaacatattcttaatctgggactgttctccatattggttaaaagtgaaacatgttgtcaacctggttgtgatggacccatttgttgacctggccatcaccatctgtattgtcttaaatactcttttcatggccatggagcactatccaatgacggaccatttcaataatgtgcttacagtaggaaacttggttttcactgggatctttacagcagaaatgtttctgaaaattattgccatggatccttactattatttccaagaaggctggaatatctttgacggttttattgtgacgcttagcctggtagaacttggactcgccaatgtggaaggattatctgttctccgttcatttcgattgctgcgagttttcaagttggcaaaatcttggccaacgttaaatatgctaataaagatcatcggcaattccgtgggggctctgggaaatttaaccctcgtcttggccatcatcgtcttcatttttgccgtggtcggcatgcagctctttggtaaaagctacaaagattgtgtctgcaagatcgccagtgattgtcaactcccacgctggcacatgaatgacttcttccactccttcctgattgtgttccgcgtgctgtgtggggagtggatagagaccatgtgggactgtatggaggttgctggtcaagccatgtgccttactgtcttcatgatggtcatggtgattggaaacctagtggtaagtatcaaggttacaagacaggtttaaggagaccaatagaaactgggcttgtcgagacagagaagactcttgcgtttctgataggcacctattggtcttactgacatccactttgcctttctctccacaggtcctgaatctctttctggccttgcttctgagctcatttagtgcagacaaccttgcagccactgatgatgataatgaaatgaataatctccaaattgctgtggataggatgcacaaaggagtagcttatgtgaaaagaaaaatatatgaatttattcaacagtccttcattaggaaacaaaagattttagatgaaattaaaccacttgatgatctaaacaacaagaaagacagttgtatgtccaatcatacagcagaaattgggaaagatcttgactatcttaaagatgtaaatggaactacaagtggtataggaactggcagcagtgttgaaaaatacattattgatgaaagtgattacatgtcattcataaacaaccccagtcttactgtgactgtaccaattgctgtaggagaatctgactttgaaaatttaaacacggaagactttagtagtgaatcggatctggaagaaagcaaagagaaactgaatgaaagcagtagctcatcagaaggtagcactgtggacatcggcgcacctgtagaagaacagcccgtagtggaacctgaagaaactcttgaaccagaagcttgtttcactgaaggctgtgtacaaagattcaagtgttgtcaaatcaatgtggaagaaggcagaggaaaacaatggtggaacctgagaaggacgtgtttccgaatagttgaacataactggtttgagaccttcattgttttcatgattctccttagtagtggtgctctggcatttgaagatatatatattgatcagcgaaagacgattaagacgatgttggaatatgctgacaaggttttcacttacattttcattctggaaatgcttctaaaatgggtggcatatggctatcaaacatatttcaccaatgcctggtgttggctggacttcttaattgttgatgtttcattggtcagtttaacagcaaatgccttgggttactcagaacttggagccatcaaatctctcaggacactaagagctctgagacctctaagagccttatctcgatttgaagggatgagggtggttgtgaatgcccttttaggagcaattccatccatcatgaatgtgcttctggtttgtcttatattctggctaattttcagcatcatgggcgtaaatttgtttgctggcaaattctaccactgtattaacaccacaactggtgacaggttcgacatcgaagacgtgaataatcatactgattgcctaaaactaatagaaagaaatgagactgctcgatggaaaaatgtgaaagtaaactttgataatgtaggatttgggtatctctctttgcttcaagttgccacattcaaaggatggatggatataatgtatgcagcagttgattccagaaatgtggaactccagcctaagtatgaagaaagtctgtacatgtatctttactttgttattttcatcatctttgggtccttcttcaccttgaacctgtttattggtgtcatcatagataatttcaaccagcagaaaaagaagtttggaggtcaagacatctttatgacagaagaacagaagaaatactataatgcaatgaaaaaattaggatcgaaaaaaccgcaaaagcctatacctcgaccaggaaacaaatttcaaggaatggtctttgacttcgtaaccagacaagtttttgacataagcatcatgattctcatctgtcttaacatggtcacaatgatggtggaaacagatgaccagagtgaatatgtgactaccattttgtcacgcatcaatctggtgttcattgtgctatttactggagagtgtgtactgaaactcatctctctacgccattattattttaccattggatggaatatttttgattttgtggttgtcattctctccattgtaggtatgtttcttgccgagctgatagaaaagtatttcgtgtcccctaccctgttccgagtgatccgtcttgctaggattggccgaatcctacgtctgatcaaaggagcaaaggggatccgcacgctgctctttgctttgatgatgtcccttcctgcgttgtttaacatcggcctcctactcttcctagtcatgttcatctacgccatctttgggatgtccaactttgcctatgttaagagggaagttgggatcgatgacatgttcaactttgagacctttggcaacagcatgatctgcctattccaaattacaacctctgctggctgggatggattgctagcacccattctcaacagtaagccacccgactgtgaccctaataaagttaaccctggaagctcagttaagggagactgtgggaacccatctgttggaattttcttttttgtcagttacatcatcatatccttcctggttgtggtgaacatgtacatcgcggtcatcctggagaacttcagtgttgctactgaagaaagtgcagagcctctgagtgaggatgactttgagatgttctatgaggtttgggagaagtttgatcccgatgcaactcagttcatggaatttgaaaaattatctcagtttgcagctgcgcttgaaccgcctctcaatctgccacaaccaaacaaactccagctcattgccatggatttgcccatggtgagtggtgaccggatccactgtcttgatatcttatttgcttttacaaagcgggttctaggagagagtggagagatggatgctctacgaatacagatggaagagcgattcatggcttccaatccttccaaggtctcctatcagccaatcactactactttaaaacgaaaacaagaggaagtatctgctgtcattattcagcgtgcttacagacgccaccttttaaagcgaactgtaaaacaagcttcctttacgtacaataaaaacaaaatcaaaggtggggctaatcttcttataaaagaagacatgataattgacagaataaatgaaaactctattacagaaaaaactgatctgaccatgtccactgcagcttgtccaccttcctatgaccgggtgacaaagccaattgtggaaaaacatgagcaagaaggcaaagatgaaaaagccaaagggaaataacgatcgcctgcaggagtcaagggcgaattcgtttatcagggaaagtaaaaagtaaggatccgcccctctccctcccccccccctaacgttactggccgaagccgcttggaataaggccggtgtgcgtttgtctatatgttattttccaccatattgccgtcttttggcaatgtgagggcccggaaacctggccctgtcttcttgacgagcattcctaggggtctttcccctctcgccaaaggaatgc'
a = Seq(seq)
print(len(seq))
print(orffinder.getORFs(a, minimum_length=75, remove_nested=True))

size = 26
max_GC = 60
min_GC = 40
max_Tm = 60
min_Tm = 50
location_start = int(input('specify location interval for primer design, start:'))
location_end = int(input('end:'))
seq=seq[location_start:location_end]

forward_primers = []
reverse_primers = []
good_GC_forward_primers = []
good_GC_reverse_primers = []
good_GC_Tm_forward_primers = []
good_GC_Tm_reverse_primers = []
good_GC_Tm_No_Hairpin_forward_primers = []
good_GC_Tm_No_Hairpin_reverse_primers = []
ideal_forward_primers = []
ideal_reverse_primers = []

for i in range(len(seq) - (size - 1)):
    a = seq[i:i + size]
    if a.find('AAAA') == -1 and a.find('TTTT') == -1 and a.find('CCCC') == -1 and a.find('GGGG') == -1:
        primer = Seq(a)
        forward_primers.append(primer)

for k in forward_primers:
    reverse_primers.append(k.reverse_complement())

for m in forward_primers:
    if GC(m) > min_GC and GC(m) < max_GC:
        good_GC_forward_primers.append(m)
for l in good_GC_forward_primers:
    if mt.Tm_NN(l) > min_Tm and mt.Tm_NN(l) < max_Tm:
        good_GC_Tm_forward_primers.append(l)

for t in good_GC_Tm_forward_primers:
    a = primer3.bindings.calc_hairpin(str(t), mv_conc=50.0, dv_conc=1.5, dntp_conc=0.6, dna_conc=50.0, temp_c=37.0,
                                     max_loop=30, output_structure=False)
    if a.structure_found == False:
        good_GC_Tm_No_Hairpin_forward_primers.append(t)

for m in reverse_primers:
    if GC(m) > min_GC and GC(m) < max_GC:
        good_GC_reverse_primers.append(m)
for l in good_GC_reverse_primers:
    if mt.Tm_NN(l) > min_Tm and mt.Tm_NN(l) < max_Tm:
        good_GC_Tm_reverse_primers.append(l)

for t in good_GC_Tm_reverse_primers:
    a = primer3.bindings.calc_hairpin(str(t), mv_conc=50.0, dv_conc=1.5, dntp_conc=0.6, dna_conc=50.0, temp_c=37.0,
                                     max_loop=30, output_structure=False)
    if a.structure_found == False:
        good_GC_Tm_No_Hairpin_reverse_primers.append(t)

for l in good_GC_Tm_No_Hairpin_forward_primers:
    a = primer3.bindings.calc_homodimer(str(l), mv_conc=50.0, dv_conc=1.5, dntp_conc=0.6, dna_conc=50.0, temp_c=37.0,
                                       max_loop=30, output_structure=False)
    if a.tm < 0:
        ideal_forward_primers.append(l)

for l in good_GC_Tm_No_Hairpin_reverse_primers:
    a = primer3.bindings.calc_homodimer(str(l), mv_conc=50.0, dv_conc=1.5, dntp_conc=0.6, dna_conc=50.0, temp_c=37.0,
                                       max_loop=30, output_structure=False)
    if a.tm < 0:
        ideal_reverse_primers.append(l)

location_forward = []
location_reverse = []
F = []
R = []
bp = []

for i in ideal_forward_primers:
    for k in ideal_reverse_primers:
        if seq.find(str(k.reverse_complement())) - seq.find(str(i)) > len(str(i)):
            F.append(i)
            R.append(k)
            location_forward.append(seq.find(str(i)))
            location_reverse.append(seq.find(str(k.reverse_complement())))
            bp.append(seq.find(str(k.reverse_complement())) - seq.find(str(i)))

for i in range(len(F)):

    k = primer3.bindings.calc_heterodimer(str(F[i]), str(R[i]), mv_conc=50.0, dv_conc=1.5, dntp_conc=0.6, dna_conc=50.0,
                                         temp_c=37, max_loop=30, output_structure=False)

    if k.tm < 0 and abs((mt.Tm_NN(F[i]) - mt.Tm_NN(R[i])) < 2) and abs(
            GC(F[i]) - GC(R[i])) < 2 and 200<bp[i]<300:
        print(('FORWARD', F[i], 'LOCATION', location_start + location_forward[i], 'Tm', mt.Tm_NN(F[i]), 'GC',
               GC(F[i]),
               'REVERSE', R[i], 'LOCATION', location_start + location_reverse[i], 'Tm', mt.Tm_NN(R[i]), 'GC', GC(R[i]),
               'Amplicon Size', bp[i]))
        with open('deneme1.txt', 'a') as f:
            a = ('FORWARD', F[i], 'LOCATION', location_start + location_forward[i], 'Tm', mt.Tm_NN(F[i]), 'GC',
                 GC(F[i]),
                 'REVERSE', R[i], 'LOCATION', location_start + location_reverse[i], 'Tm', mt.Tm_NN(R[i]), 'GC', GC(R[i]),
                 'Amplicon Size', bp[i])
            f.writelines(str(a))
            f.writelines('\n')
