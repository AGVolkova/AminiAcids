Alphabet = {57: 'G', 71: 'A', 87: 'S', 97: 'P', 99: 'V', 101: 'T', 103: 'C', 113: 'L', 114: 'N', 115: 'D', 128: 'Q', 129: 'E', 131: 'M', 137: 'H', 147: 'F', 156: 'R', 163: 'Y', 186: 'W'}
AminoAcids_Mass={'G':57, 'A':71, 'S':87, 'P':97, 'V':99, 'T':101, 'C':103, 'I':113, 'L':113, 'N':114, 'D':115, 'K':128, 'Q':128, 'E':129, 'M':131, 'H':137, 'F':147, 'R':156, 'Y':163, 'W':186}
All_AminoAcids=['G', 'A', 'S', 'P', 'V', 'T', 'C', 'I', 'L', 'N', 'D', 'K', 'Q', 'E', 'M', 'H', 'F', 'R', 'Y', 'W']
def Expand(Peptides):
    if len(Peptides)==1:
        Peptides= [[int(x)] for x in Alphabet]
        #Peptides.append([0])
    else:
        l=[]
        for i in Peptides:
            for j in Alphabet:
                a = i.copy()
                a.append(int(j))
                l.append(a)
        Peptides=l
        #Peptides.append([0])
    return Peptides
#print(Expand([[57], [71], [87], [97], [99], [101], [103], [113], [114], [115], [128], [129], [131], [137], [147], [156], [163], [186]]))


def MassesToPeptode(Masses):
    Masses=[i for i in Masses if i !=0]
    if Masses and len(Masses):
        P=Alphabet[Masses[0]]
        for i in range(1, len(Masses)):
            P = P + Alphabet[Masses[i]]
    else:
        return None
    return P

def CyclicSpectrum(Peptide, Alphabet, AminoAcidMass):
    PrefixMass=[0]
    for i in range(1, len(Peptide)+1):
        for j in Alphabet:
            if Peptide[i-1]==j:
                PrefixMass.append(PrefixMass[i-1]+AminoAcidMass[j])
    peptidemass=max(PrefixMass)
    CyclicSpectrum=[0]
    for i in range(0, len(Peptide)):
        for j in range(i+1, len(Peptide)+1):
            CyclicSpectrum.append(PrefixMass[j]-PrefixMass[i])
            if i >0 and j< (len(Peptide)):
                CyclicSpectrum.append(peptidemass-(PrefixMass[j]-PrefixMass[i]))
    Sorted_list=sorted(CyclicSpectrum)
    return Sorted_list
#print('yes', CyclicSpectrum(MassesToPeptode([113, 128, 186]), All_AminoAcids, AminoAcids_Mass))

def LinearSpectrum(Peptide, Alphabet, AminoAcidMass):
    PrefixMass=[0]
    for i in range(1, len(Peptide)+1):
        for j in Alphabet:
            if Peptide[i-1]==j:
                PrefixMass.append(PrefixMass[i-1]+AminoAcidMass[j])
    LinearSpectrum=[0]
    for i in range(0, len(Peptide)):
        for j in range(i+1, len(Peptide)+1):
            LinearSpectrum.append(PrefixMass[j]-PrefixMass[i])
    Sorted_list=sorted(LinearSpectrum)
    return Sorted_list
#print(LinearSpectrum('W',All_AminoAcids, AminoAcids_Mass))

def Consistency(Maly, Bolshoyy):
    Bolshoy=Bolshoyy.copy()
    for i in Maly:
        if i in Bolshoy:
            Bolshoy.remove(i)
        else: return 'No'
    return 'Yes'
#print(Consistency(LinearSpectrum('GKLENW',All_AminoAcids, AminoAcids_Mass), [0, 57, 113, 114, 128, 129, 185, 186, 241, 242, 243, 298, 300, 356, 370, 427, 429, 484, 541, 542, 670, 727]))
#print(Consistency([128], [0, 113, 128, 186, 241, 299, 314, 427]))


def Consistency_check(M, BB):
    B=BB.copy()
    for i in M:
        if i in B:
            B.remove(i)
        else:
            print(i)
            return i

    return 'Yes'
#print(Consistency(LinearSpectrum(MassesToPeptode(i),All_AminoAcids, AminoAcids_Mass), Spectrum))
#print(Consistency_check(LinearSpectrum(MassesToPeptode([128]),All_AminoAcids, AminoAcids_Mass), [0, 113, 128, 186, 241, 299, 314, 427]))
def CyclopeptideSequencing(S):
    kkk=S.split(' ')
    Spectrum=[int(x) for x in kkk]
    CandidatePeptides=['']
    FinalPeptides=[]
    while len(CandidatePeptides)!=0:
        CandidatePeptides=Expand(CandidatePeptides)
        for i in CandidatePeptides:
            if sum(i)==max(Spectrum):
                if CyclicSpectrum(MassesToPeptode(i), All_AminoAcids, AminoAcids_Mass)==Spectrum and i not in FinalPeptides:
                    FinalPeptides.append(i)
                CandidatePeptides=[x for x in CandidatePeptides if x !=i]
            elif Consistency(LinearSpectrum(MassesToPeptode(i),All_AminoAcids, AminoAcids_Mass), Spectrum)=='No':
                CandidatePeptides=[x for x in CandidatePeptides if x !=i]
    for i in FinalPeptides:
        stroka=[str(k) for k in i]
        print('-'.join(stroka))
    return FinalPeptides,

print(CyclopeptideSequencing('0 71 71 71 103 113 114 115 129 142 156 174 184 185 186 232 242 245 256 270 271 299 303 313 341 342 345 359 374 385 412 416 416 428 455 456 456 487 487 488 515 527 527 531 558 569 584 598 601 602 630 640 644 672 673 687 698 701 711 757 758 759 769 787 801 814 828 829 830 840 872 872 872 943'))
