Alphabet = {57: 'G', 71: 'A', 87: 'S', 97: 'P', 99: 'V', 101: 'T', 103: 'C', 113: 'L', 114: 'N', 115: 'D', 128: 'Q', 129: 'E', 131: 'M', 137: 'H', 147: 'F', 156: 'R', 163: 'Y', 186: 'W'}
AminoAcids_Mass={'G':57, 'A':71, 'S':87, 'P':97, 'V':99, 'T':101, 'C':103, 'I':113, 'L':113, 'N':114, 'D':115, 'K':128, 'Q':128, 'E':129, 'M':131, 'H':137, 'F':147, 'R':156, 'Y':163, 'W':186}
All_AminoAcids=['G', 'A', 'S', 'P', 'V', 'T', 'C', 'I', 'L', 'N', 'D', 'K', 'Q', 'E', 'M', 'H', 'F', 'R', 'Y', 'W']
def Expand(Peptides):
    if len(Peptides)==1:
        Peptides= [[x] for x in Alphabet]
        #Peptides.append([0])
    else:
        l=[]
        for i in Peptides:
            for j in Alphabet:
                a = i.copy()
                a.append(j)
                l.append(a)
        Peptides=l
        #Peptides.append([0])
    return Peptides
#print(Expand([[57], [71], [87], [97], [99], [101], [103], [113], [114], [115], [128], [129], [131], [137], [147], [156], [163], [186]]))

def Consistency(Maly, Bolshoy):
    for i in Maly:
        if i in Bolshoy:
            Bolshoy=[x for x in Bolshoy if x!=i]
        else: return 'No'
    return 'Yes'

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

def CyclopeptideSequencing(S):
    Spectrum=S.split(' ')
    Spectrum=[int(x) for x in Spectrum]
    CandidatePeptides=['']
    FinalPeptides=[]
    while len(CandidatePeptides)!=0:
        CandidatePeptides=Expand(CandidatePeptides)
        for i in CandidatePeptides:
            if sum(i)==max(Spectrum):
                if CyclicSpectrum(MassesToPeptode(i), All_AminoAcids, AminoAcids_Mass)==Spectrum and i not in FinalPeptides:
                    FinalPeptides.append(i)
                CandidatePeptides=[x for x in CandidatePeptides if x !=i]
            elif Consistency(i, Spectrum)=='No':
                CandidatePeptides = [x for x in CandidatePeptides if x != i]
    for i in FinalPeptides:
        stroka=[str(k) for k in i]
        print('-'.join(stroka))
    return FinalPeptides

print(CyclopeptideSequencing('0 71 99 101 103 114 128 129 147 156 170 200 213 229 242 248 250 259 284 285 299 341 343 351 356 376 388 406 412 413 442 455 459 479 490 507 513 535 541 558 569 589 593 606 635 636 642 660 672 692 697 705 707 749 763 764 789 798 800 806 819 835 848 878 892 901 919 920 934 945 947 949 977 1048'))
# for i in FinalPeptides:
#     print(' '.join(i))