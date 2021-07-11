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

def Cyclopeptide_Scoring(Theoretical, Experimental):
    Exper=Experimental.split(' ')
    Exp=[int(x) for x in Exper]
    print(Exp)
    Score=0
    for i in Theoretical:
        if i in Exp:
            Score=Score+1
            Exp.remove(i)
    return Score
print(Cyclopeptide_Scoring(CyclicSpectrum('MPHRNDKMKWAWRMVDHEICGFFMLPMSHEMHIGCFGNRD', All_AminoAcids, AminoAcids_Mass), '0 57 57 57 71 97 97 97 99 99 103 103 103 103 113 113 113 114 114 115 115 128 128 129 129 131 131 131 131 131 137 137 137 137 147 147 147 154 156 156 160 170 171 186 186 204 204 210 212 214 216 226 227 227 228 229 230 234 234 234 234 240 242 243 244 250 250 251 257 257 259 260 266 267 268 270 274 274 278 287 294 301 307 307 307 311 314 317 327 331 337 339 339 341 342 342 345 347 347 351 355 357 364 365 366 369 371 371 373 381 385 386 388 390 391 396 397 397 413 413 414 420 421 425 430 430 436 438 440 442 443 448 450 454 454 456 468 470 473 474 477 478 478 482 484 485 487 493 494 494 500 500 501 502 503 504 510 533 534 538 541 543 544 551 553 557 561 561 567 571 571 572 577 579 581 584 585 587 591 593 595 598 599 600 601 603 605 608 612 614 616 624 625 631 631 634 637 640 641 643 647 650 656 658 664 670 680 684 687 688 690 692 692 698 704 707 714 716 718 719 721 727 727 727 728 728 728 730 730 731 734 737 738 745 747 750 755 758 761 761 768 770 772 781 784 787 795 795 798 805 807 811 817 821 826 827 829 829 834 838 840 841 847 847 850 851 854 855 858 858 859 859 865 868 871 874 875 875 875 881 884 898 902 907 908 908 913 913 918 926 932 932 934 937 938 944 954 954 954 955 957 957 958 962 965 978 978 981 984 984 984 985 987 988 989 990 994 994 1001 1011 1012 1015 1026 1027 1029 1031 1033 1035 1037 1039 1041 1041 1044 1045 1046 1056 1057 1061 1069 1069 1072 1078 1085 1085 1086 1097 1097 1102 1104 1114 1115 1115 1118 1118 1118 1123 1125 1125 1132 1132 1134 1141 1141 1141 1142 1144 1146 1149 1154 1161 1164 1166 1168 1169 1170 1171 1180 1182 1183 1184 1188 1188 1189 1194 1198 1200 1212 1215 1215 1217 1228 1228 1232 1239 1245 1245 1247 1249 1249 1251 1262 1263 1264 1265 1268 1271 1272 1272 1279 1281 1283 1285 1288 1289 1291 1295 1296 1298 1299 1299 1301 1301 1314 1320 1327 1329 1331 1343 1345 1346 1348 1348 1352 1352 1353 1359 1368 1375 1375 1376 1378 1380 1380 1384 1384 1388 1392 1396 1397 1398 1398 1404 1408 1408 1409 1411 1411 1413 1414 1419 1426 1432 1432 1437 1446 1451 1457 1459 1461 1462 1466 1466 1468 1474 1475 1479 1479 1483 1483 1487 1489 1491 1505 1508 1510 1511 1515 1515 1516 1522 1525 1528 1529 1531 1532 1535 1537 1539 1539 1545 1550 1554 1563 1565 1575 1579 1579 1579 1581 1582 1584 1588 1592 1596 1602 1603 1605 1607 1610 1612 1612 1618 1618 1618 1619 1620 1620 1625 1632 1634 1635 1636 1638 1638 1640 1645 1664 1665 1668 1676 1676 1676 1682 1685 1685 1687 1695 1702 1706 1709 1712 1713 1715 1716 1718 1718 1719 1721 1726 1731 1731 1733 1734 1735 1739 1749 1749 1749 1749 1753 1754 1762 1765 1766 1766 1773 1783 1784 1788 1788 1789 1789 1792 1802 1805 1805 1808 1811 1816 1821 1823 1829 1834 1846 1846 1846 1848 1849 1850 1850 1852 1852 1853 1859 1862 1862 1863 1867 1869 1876 1880 1880 1886 1886 1886 1891 1892 1899 1899 1902 1902 1907 1915 1918 1919 1936 1939 1943 1944 1945 1949 1952 1952 1958 1959 1964 1965 1965 1975 1975 1977 1978 1979 1979 1983 1983 1989 1990 1993 1996 2000 2000 2004 2006 2009 2009 2015 2017 2018 2021 2023 2023 2028 2030 2040 2054 2055 2057 2058 2062 2062 2072 2073 2075 2076 2078 2080 2086 2086 2089 2096 2096 2101 2103 2103 2106 2110 2112 2115 2115 2116 2120 2122 2122 2125 2126 2127 2130 2133 2135 2141 2153 2154 2156 2158 2168 2174 2187 2189 2192 2193 2193 2199 2201 2202 2204 2205 2206 2209 2209 2213 2217 2218 2218 2219 2223 2225 2229 2230 2234 2238 2238 2243 2243 2246 2247 2256 2259 2259 2263 2266 2266 2267 2282 2286 2288 2289 2295 2305 2305 2310 2312 2315 2316 2321 2322 2330 2330 2330 2332 2333 2334 2335 2336 2337 2337 2340 2342 2345 2346 2349 2349 2359 2360 2365 2369 2379 2379 2386 2390 2390 2392 2392 2396 2403 2403 2413 2417 2422 2423 2433 2433 2436 2437 2440 2442 2445 2445 2446 2447 2448 2449 2450 2452 2452 2452 2460 2461 2466 2467 2470 2472 2477 2477 2487 2493 2494 2496 2500 2515 2516 2516 2519 2523 2523 2526 2535 2536 2539 2539 2544 2544 2548 2552 2553 2557 2559 2563 2564 2564 2565 2569 2573 2573 2576 2577 2578 2580 2581 2583 2589 2589 2590 2593 2595 2608 2614 2624 2626 2628 2629 2641 2647 2649 2652 2655 2656 2657 2660 2660 2662 2666 2667 2667 2670 2672 2676 2679 2679 2681 2686 2686 2693 2696 2696 2702 2704 2706 2707 2709 2710 2720 2720 2724 2725 2727 2728 2742 2752 2754 2759 2759 2761 2764 2765 2767 2773 2773 2776 2778 2782 2782 2786 2789 2792 2793 2799 2799 2803 2803 2804 2805 2807 2807 2817 2817 2818 2823 2824 2830 2830 2833 2837 2838 2839 2843 2846 2863 2864 2867 2875 2880 2880 2883 2883 2890 2891 2896 2896 2896 2902 2902 2906 2913 2915 2919 2920 2920 2923 2929 2930 2930 2932 2932 2933 2934 2936 2936 2936 2948 2953 2959 2961 2966 2971 2974 2977 2977 2980 2990 2993 2993 2994 2994 2998 2999 3009 3016 3016 3017 3020 3028 3029 3033 3033 3033 3033 3043 3047 3048 3049 3051 3051 3056 3061 3063 3064 3064 3066 3067 3069 3070 3073 3076 3080 3087 3095 3097 3097 3100 3106 3106 3106 3114 3117 3118 3137 3142 3144 3144 3146 3147 3148 3150 3157 3162 3162 3163 3164 3164 3164 3170 3170 3172 3175 3177 3179 3180 3186 3190 3194 3198 3200 3201 3203 3203 3203 3207 3217 3219 3228 3232 3237 3243 3243 3245 3247 3250 3251 3253 3254 3257 3260 3266 3267 3267 3271 3272 3274 3277 3291 3293 3295 3299 3299 3303 3303 3307 3308 3314 3316 3316 3320 3321 3323 3325 3331 3336 3345 3350 3350 3356 3363 3368 3369 3371 3371 3373 3374 3374 3378 3384 3384 3385 3386 3390 3394 3398 3398 3402 3402 3404 3406 3407 3407 3414 3423 3429 3430 3430 3434 3434 3436 3437 3439 3451 3453 3455 3462 3468 3481 3481 3483 3483 3484 3486 3487 3491 3493 3494 3497 3499 3501 3503 3510 3510 3511 3514 3517 3518 3519 3520 3531 3533 3533 3535 3537 3537 3543 3550 3554 3554 3565 3567 3567 3570 3582 3584 3588 3593 3594 3594 3598 3599 3600 3602 3611 3612 3613 3614 3616 3618 3621 3628 3633 3636 3638 3640 3641 3641 3641 3648 3650 3650 3657 3657 3659 3664 3664 3664 3667 3667 3668 3678 3680 3685 3685 3696 3697 3697 3704 3710 3713 3713 3721 3725 3726 3736 3737 3738 3741 3741 3743 3745 3747 3749 3751 3753 3755 3756 3767 3770 3771 3781 3788 3788 3792 3793 3794 3795 3797 3798 3798 3798 3801 3804 3804 3817 3820 3824 3825 3825 3827 3828 3828 3828 3838 3844 3845 3848 3850 3850 3856 3864 3869 3869 3874 3874 3875 3880 3884 3898 3901 3907 3907 3907 3908 3911 3914 3917 3923 3923 3924 3924 3927 3928 3931 3932 3935 3935 3941 3942 3944 3948 3953 3953 3955 3956 3961 3965 3971 3975 3977 3984 3987 3987 3995 3998 4001 4010 4012 4014 4021 4021 4024 4027 4032 4035 4037 4044 4045 4048 4051 4052 4052 4054 4054 4054 4055 4055 4055 4061 4063 4064 4066 4068 4075 4078 4084 4090 4090 4092 4094 4095 4098 4102 4112 4118 4124 4126 4132 4135 4139 4141 4142 4145 4148 4151 4151 4157 4158 4166 4168 4170 4174 4177 4179 4181 4182 4183 4184 4187 4189 4191 4195 4197 4198 4201 4203 4205 4210 4211 4211 4215 4221 4221 4225 4229 4231 4238 4239 4241 4244 4248 4249 4272 4278 4279 4280 4281 4282 4282 4288 4288 4289 4295 4297 4298 4300 4304 4304 4305 4308 4309 4312 4314 4326 4328 4328 4332 4334 4339 4340 4342 4344 4346 4352 4352 4357 4361 4362 4368 4369 4369 4385 4385 4386 4391 4392 4394 4396 4397 4401 4409 4411 4411 4413 4416 4417 4418 4425 4427 4431 4435 4435 4437 4440 4440 4441 4443 4443 4445 4451 4455 4465 4468 4471 4475 4475 4475 4481 4488 4495 4504 4508 4508 4512 4514 4515 4516 4522 4523 4525 4525 4531 4532 4532 4538 4539 4540 4542 4548 4548 4548 4548 4552 4553 4554 4555 4555 4556 4566 4568 4570 4572 4578 4578 4596 4596 4611 4612 4622 4626 4626 4628 4635 4635 4635 4645 4645 4645 4645 4651 4651 4651 4651 4651 4653 4653 4654 4654 4667 4667 4668 4668 4669 4669 4669 4679 4679 4679 4679 4683 4683 4685 4685 4685 4711 4725 4725 4725 4782'))

def Trim(Leaderboard, Spectrum, N):
    dict={}
    for i in Leaderboard:
        dict[i]=Cyclopeptide_Scoring(i, Spectrum)
    trimmed_list=sorted(dict.values)[0:N]

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

def LeaderboardCyclopeptideSequencing(S, N):
    kkk = S.split(' ')
    Spectrum = [int(x) for x in kkk]
    Leaderboard=set()
    LeaderPeptide=[]
    while len(Leaderboard)!=0:
        Leaderboard=Expand(Leaderboard)
        for i in Leaderboard:
            if sum(i)==max(Spectrum):
                if Cyclopeptide_Scoring(i, Spectrum)>Score(LeaderPeptide, Spectrum):
                    LeaderPeptide=i
            elif sum (i)>max(Spectrum):
                Leaderboard.remove(i)
        Leaderboard=Trim(Leaderboard, Spectrum, N)
        return LeaderPeptide
print(LeaderboardCyclopeptideSequencing('0, 71, 113, 129, 147, 200, 218, 260, 313, 331, 347, 389, 460,', 10))