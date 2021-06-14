Alphabet = {57: 'G', 71: 'A', 87: 'S', 97: 'P', 99: 'V', 101: 'T', 103: 'C', 113:'I/L', 114: 'N', 115: 'D', 128: 'K/Q', 129: 'E', 131: 'M', 137: 'H', 147: 'F', 156: 'R', 163: 'Y', 186: 'W'}
AminoAcids_Mass={'G':57, 'A':71, 'S':87, 'P':97, 'V':99, 'T':101, 'C':103, 'I':113, 'L':113, 'N':114, 'D':115, 'K':128, 'Q':128, 'E':129, 'M':131, 'H':137, 'F':147, 'R':156, 'Y':163, 'W':186}

def Expand(Peptides):
    if len(Peptides)==1:
        Peptides= [[x] for x in Alphabet]
        Peptides.append([0])
    else:
        l=[]
        for i in Peptides:
            for j in Alphabet:
                a = i.copy()
                a.append(j)
                l.append(a)
        Peptides=l
        Peptides.append([0])
    return Peptides
print(Expand([[57], [71], [87], [97], [99], [101], [103], [113], [114], [115], [128], [129], [131], [137], [147], [156], [163], [186]]))

def Consistency(Maly, Bolshoy):
    for i in Maly:
        if i in Bolshoy:
            Bolshoy=[x for x in Bolshoy if x!=i]
        else: return 'No'
    return 'Yes'

def CyclicSpectrum(Peptide, Alphabet, AminoAcidMass):
    PrefixMass=[0]
    for i in range(1, len(Peptide)+1):
        for j in Alphabet:
            if Peptide[i-1]==j:
                PrefixMass.append(PrefixMass[i-1]+AminoAcidMass[j])
    peptidemass=max(PrefixMass)
    print(PrefixMass)
    CyclicSpectrum=[0]
    for i in range(0, len(Peptide)):
        for j in range(i+1, len(Peptide)+1):
            CyclicSpectrum.append(PrefixMass[j]-PrefixMass[i])
            if i >0 and j< (len(Peptide)):
                CyclicSpectrum.append(peptidemass-(PrefixMass[j]-PrefixMass[i]))
                print('Yes')
    Sorted_list=sorted(CyclicSpectrum)
    return Sorted_list

def CyclopeptideSequencing(S):
    Spectrum=S.split(' ')
    Spectrum=[int(x) for x in Spectrum]
    CandidatePeptides=['']
    FinalPeptides=[]
    while len(CandidatePeptides)!=0:
        CandidatePeptides=Expand(CandidatePeptides)
        print(CandidatePeptides)
        for i in CandidatePeptides:
            if sum(i)==max(Spectrum):
                if CyclicSpectrum(i, Alphabet, AminoAcids_Mass)==Spectrum and i not in FinalPeptides:
                    FinalPeptides.append(i)
                CandidatePeptides=[x for x in CandidatePeptides if x !=i]
            elif Consistency(i, Spectrum)=='No':
                CandidatePeptides = [x for x in CandidatePeptides if x != i]
    return FinalPeptides

print(CyclopeptideSequencing('0 113 128 186 241 299 314 4270 113 128 186 241 299 314 427'))