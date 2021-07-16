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
        print(Leaderboard)
        for i in Leaderboard:
            if sum(i)==max(Spectrum):
                if Cyclopeptide_Scoring(CyclicSpectrum(MassesToPeptode(i), All_AminoAcids, AminoAcids_Mass), Spectrum)>Cyclopeptide_Scoring(LeaderPeptide, Spectrum):
                    LeaderPeptide=i
                    Leaderboard=[k for k in Leaderboard if k!=i]
            elif sum(i)>max(Spectrum):
                print('Yes1', sum(i))
                Leaderboard=[k for k in Leaderboard if k!=i]
            elif Consistency(LinearSpectrum(MassesToPeptode(i),All_AminoAcids, AminoAcids_Mass), Spectrum)=='No':
                Leaderboard=[k for k in Leaderboard if k!=i]

        Trimmed=Trim(Leaderboard, Spectrum, N)
        Leaderboard=[k.split( ) for k in Trimmed]
    return '-'.join([str(x) for x in LeaderPeptide])
print(LeaderboardCyclopeptideSequencing('0 71 87 97 97 97 99 99 101 114 114 128 128 129 131 131 147 156 156 170 185 186 186 186 196 196 211 215 215 226 242 244 253 257 257 262 275 284 285 293 312 313 314 314 317 317 340 342 343 352 371 371 372 382 382 385 389 401 412 440 441 443 445 448 448 449 460 468 471 473 479 481 486 496 499 499 500 532 557 557 559 568 571 574 576 578 578 585 593 596 597 597 604 627 631 634 646 654 656 656 663 675 682 685 688 692 694 702 705 724 725 743 753 753 755 762 762 774 783 785 789 790 793 813 816 819 822 833 838 842 849 850 852 853 871 890 891 909 914 916 918 936 939 939 944 947 948 949 950 960 967 970 1005 1005 1008 1019 1019 1036 1037 1038 1045 1046 1057 1064 1064 1067 1070 1075 1091 1095 1102 1104 1106 1133 1133 1135 1135 1137 1156 1165 1167 1175 1178 1188 1193 1201 1205 1220 1222 1223 1224 1230 1231 1232 1234 1249 1253 1261 1266 1276 1279 1287 1289 1298 1317 1319 1319 1321 1321 1348 1350 1352 1359 1363 1379 1384 1387 1390 1390 1397 1408 1409 1416 1417 1418 1435 1435 1446 1449 1449 1484 1487 1494 1504 1505 1506 1507 1510 1515 1515 1518 1536 1538 1540 1545 1563 1564 1583 1601 1602 1604 1605 1612 1616 1621 1632 1635 1638 1641 1661 1664 1665 1669 1671 1680 1692 1692 1699 1701 1701 1711 1729 1730 1749 1752 1760 1762 1766 1769 1772 1779 1798 1798 1800 1808 1820 1823 1827 1850 1857 1857 1858 1861 1869 1876 1876 1878 1880 1883 1886 1895 1897 1897 1922 1954 1955 1955 1958 1968 1973 1975 1981 1983 1986 1994 2005 2006 2006 2009 2011 2013 2014 2042 2053 2065 2069 2072 2072 2082 2083 2083 2102 2111 2112 2114 2137 2137 2140 2140 2141 2142 2161 2169 2170 2179 2192 2197 2197 2201 2210 2212 2228 2239 2239 2243 2258 2258 2268 2268 2268 2269 2284 2298 2298 2307 2323 2323 2325 2326 2326 2340 2340 2353 2355 2355 2357 2357 2357 2367 2383 2454', 197))