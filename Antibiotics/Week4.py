Alphabet = {57: 'G', 71: 'A', 87: 'S', 97: 'P', 99: 'V', 101: 'T', 103: 'C', 113: 'L', 114: 'N', 115: 'D', 128: 'Q', 129: 'E', 131: 'M', 137: 'H', 147: 'F', 156: 'R', 163: 'Y', 186: 'W'}
AminoAcids_Mass={'G':57, 'A':71, 'S':87, 'P':97, 'V':99, 'T':101, 'C':103, 'I':113, 'L':113, 'N':114, 'D':115, 'K':128, 'Q':128, 'E':129, 'M':131, 'H':137, 'F':147, 'R':156, 'Y':163, 'W':186}
All_AminoAcids=['G', 'A', 'S', 'P', 'V', 'T', 'C', 'I', 'L', 'N', 'D', 'K', 'Q', 'E', 'M', 'H', 'F', 'R', 'Y', 'W']

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
#print(CyclicSpectrum('NQEL', All_AminoAcids, AminoAcids_Mass))

def Cyclopeptide_Scoring(Theoretical, Exper):
    # Exper=Experimental.split(' ')
    # Exp=[int(x) for x in Exper]
    # print(Exp)
    Exp=Exper.copy()
    Score=0
    if isinstance(Theoretical, list)==True:
        for i in Theoretical:
            if i in Exp:
                Score=Score+1
                Exp.remove(i)
    else:
        if Theoretical in Exp:
            Score = Score + 1
            Exp.remove(Theoretical)
    #print(Score)
    # if Score==83:
    #     print(Theoretical)
    #print(Score)
    return Score
#print(Cyclopeptide_Scoring([113], [0, 71, 113, 129, 147, 200, 218, 260, 313, 331, 347, 389, 460]))
def Trim(Leaderboard, Spectrum, N):
    dict={}
    for i in Leaderboard:
        #print(i)
        l=[str(k) for k in i]
        dict[' '.join(l)]=Cyclopeptide_Scoring(i, Spectrum)
    TL=[k for k in dict.values()]
    sorted_values=sorted(dict.values())
    sorted_dict={}
    sorted_keys=sorted(dict, key=dict.get, reverse=True)
    #print(sorted_keys)
    #print(dict)
    return sorted_keys[0:N]
#print(Trim([[57], [347, 460], [113], [147, 218], [71,71,71], [89, 129], [218, 269, 71, 147], [113, 113, 113, 113, 113], [50], [23, 331], [347, 113], [113], [71]], [0, 71, 113, 129, 147, 200, 218, 260, 313, 331, 347, 389, 460], 10))

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

def MassesToPeptode(Masses):
    Masses=[i for i in Masses if i !=0]
    if Masses and len(Masses):
        P=Alphabet[Masses[0]]
        for i in range(1, len(Masses)):
            P = P + Alphabet[Masses[i]]
    else:
        return None
    return P

def Consistency(Maly, Bolshoyy):
    Bolshoy=Bolshoyy.copy()
    for i in Maly:
        if i in Bolshoy:
            Bolshoy.remove(i)
        else: return 'No'
    return 'Yes'

def LeaderboardCyclopeptideSequencing(S, N):
    kkk = S.split(' ')
    Spectrum = [int(x) for x in kkk]
    print(max(Spectrum))
    Leaderboard=['']
    LeaderPeptide=''
    while len(Leaderboard)!=0:
        Leaderboard=Expand(Leaderboard)
        LeadB=[]
        for j in Leaderboard:
            LeadB.append([int(k) for k in j])
        Leaderboard=LeadB
        #print(Leaderboard)
        for i in Leaderboard:
            if sum(i)==max(Spectrum):
                if Cyclopeptide_Scoring(CyclicSpectrum(MassesToPeptode(i), All_AminoAcids, AminoAcids_Mass), Spectrum)>Cyclopeptide_Scoring(LeaderPeptide, Spectrum):
                    LeaderPeptide=i
                    #Leaderboard=[k for k in Leaderboard if k!=i]
            elif sum(i)>max(Spectrum):
                #print('Yes1', sum(i))
                Leaderboard=[k for k in Leaderboard if k!=i]
            #elif Consistency(LinearSpectrum(MassesToPeptode(i),All_AminoAcids, AminoAcids_Mass), Spectrum)=='No':
             #   Leaderboard=[k for k in Leaderboard if k!=i]

        Trimmed=Trim(Leaderboard, Spectrum, N)
        Leaderboard=[k.split( ) for k in Trimmed]
    print(Leaderboard)
    print(len(LeaderPeptide))
    return '-'.join([str(x) for x in LeaderPeptide])
#print(LeaderboardCyclopeptideSequencing('0 97 99 113 114 115 128 128 147 147 163 186 227 241 242 244 244 256 260 261 262 283 291 309 330 333 340 347 385 388 389 390 390 405 435 447 485 487 503 504 518 544 552 575 577 584 599 608 631 632 650 651 653 672 690 691 717 738 745 770 779 804 818 819 827 835 837 875 892 892 917 932 932 933 934 965 982 989 1039 1060 1062 1078 1080 1081 1095 1136 1159 1175 1175 1194 1194 1208 1209 1223 1322', 1000))
def Convolution(S):
    kkk = S.split(' ')
    Spectrum = sorted([int(x) for x in kkk])
    Spec=Spectrum[:-1]
    print(Spectrum, Spec)
    res=[]
    k=0
    for i in Spectrum[::-1]:
        print(i)
        for j in Spec[::-1]:
            print(i, j)
            res.append(i-j)

        Spectrum=Spectrum[:-1]
        Spec=Spec[:-1]
    otbor=list(filter(lambda x:x!=0, res))
    return ' '.join([str(x) for x in otbor])
print(Convolution('0 57 71 87 99 101 113 114 128 137 156 185 188 214 224 227 227 236 242 293 298 301 323 325 328 341 355 364 380 399 415 424 438 451 454 456 478 481 486 537 543 552 552 555 565 591 594 623 642 651 665 666 678 680 692 708 722 779'))