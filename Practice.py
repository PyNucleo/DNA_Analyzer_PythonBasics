Valid_DNA=["A", "T", "G", "C"]
DNA_Sequence=input()
DNA_Sequence=DNA_Sequence.upper()
def valid_base(DNA_Sequence):
        for i in DNA_Sequence:
            if i not in Valid_DNA:
                return False
        return True
def not_stop(DNA_Sequence):
        return DNA_Sequence!="STOP"
def total_length(DNA):
    Length=len(DNA)
    return Length
def base_number(DNA):
    Base_Counts={"A":0,"T":0,"G":0,"C":0}
    for Base in DNA:
        Base_Counts[Base]+=1
    return Base_Counts
def gc_content_per(DNA):
    S = base_number(DNA)
    C_Content=S["C"]
    G_Content=S["G"]
    GC_Hundo=(((C_Content+G_Content)/total_length(DNA))*100)
    return GC_Hundo
def reverse_complement(DNA):
    Empty=[]
    for i in range(0, len(DNA)):
        Empty.append(DNA[i])
    Emp=''
    for j in range(len(Empty)-1,-1,-1):
        Emp+=Empty[j]
    E=''
    for comp in range(0, len(Emp)):
        if Emp[comp]=="A":
            E+="T"
        elif Emp[comp]=="T":
            E+="A"
        elif Emp[comp]=="G":
            E+="C"
        elif Emp[comp]=="C":
            E+="G"
    return E
def mrna_template(DNA_Sequence):
    empty_mRNA=[]
    DICT_CONVERT={"A":"U",
                  "G":"C",
                  "T":"A",
                  "C":"G"}
    for base in DNA_Sequence: 
        if base!="U":
            empty_mRNA.append(DICT_CONVERT[base])
        else:
            empty_mRNA.append("U")
    mRNA= ''.join(empty_mRNA)
    return mRNA
def translate_template(DNA_Sequence):
    AA_Seq=''
    count=0
    mRNA=list(mrna_template(DNA_Sequence))
    CODON_DICT={}
    Phe_Key=["UUU","UUC"]
    Leu_Key=["UUA","UUG","CUU","CUC","CUA","CUG"]
    Ile_Key=["AUU","AUC","AUA"]
    Val_Key=["GUU","GUC","GUA","GUG"]
    Ser_Key=["UCU","UCC","UCA","UCG", "AGU","AGC"]
    Pro_Key=["CCU","CCC","CCA","CCG"]
    Thr_Key=["ACU","ACC","ACA","ACG"]
    Ala_Key=["GCU","GCC","GCA","GCG"]
    Tyr_Key=["UAU","UAC"]
    His_Key=["CAU","CAC"]
    Gln_Key=["CAA","CAG"]
    Asn_Key=["AAU","AAC"]
    Lys_Key=["AAA","AAG"]
    Asp_Key=["GAU","GAC"]
    Glu_Key=["GAA","GAG"]
    Cys_Key=["UGU","UGC"]
    Trp_Key=["UGG"]
    Arg_Key=["CGU","CGC","CGA","CGG", "AGA","AGG"]
    Gly_Key=["GGU","GGC","GGA","GGG"]
    STOP_Key=["UAA","UAG","UGA"]
    Met_Key=["AUG"]
    CODON_DICT.update({Phe:"Phe" for Phe in Phe_Key})
    CODON_DICT.update({Ile:"Ile" for Ile in Ile_Key})
    CODON_DICT.update({Val:"Val" for Val in Val_Key})
    CODON_DICT.update({Ser:"Ser" for Ser in Ser_Key})         
    CODON_DICT.update({Pro:"Pro" for Pro in Pro_Key})            
    CODON_DICT.update({Leu:"Leu" for Leu in Leu_Key})
    CODON_DICT.update({Thr:"Thr" for Thr in Thr_Key})         
    CODON_DICT.update({Ala:"Ala" for Ala in Ala_Key})
    CODON_DICT.update({Tyr:"Tyr" for Tyr in Tyr_Key})
    CODON_DICT.update({His:"His" for His in His_Key})
    CODON_DICT.update({Gln:"Gln" for Gln in Gln_Key})
    CODON_DICT.update({Asn:"Asn" for Asn in Asn_Key})
    CODON_DICT.update({Lys:"Lys" for Lys in Lys_Key})            
    CODON_DICT.update({Asp:"Asp" for Asp in Asp_Key})
    CODON_DICT.update({Glu:"Glu" for Glu in Glu_Key})
    CODON_DICT.update({Cys:"Cys" for Cys in Cys_Key})
    CODON_DICT.update({Trp:"Trp" for Trp in Trp_Key})
    CODON_DICT.update({Arg:"Arg" for Arg in Arg_Key})
    CODON_DICT.update({Gly:"Gly" for Gly in Gly_Key})
    CODON_DICT.update({STOP:"STOP" for STOP in STOP_Key})
    CODON_DICT.update({Met:"Met" for Met in Met_Key})
    while len(mRNA)>=3:
        c_out = mRNA[0:3]
        s=''.join(c_out)
        del mRNA[0:3]
        AA_Seq += CODON_DICT[s] + ' '
        c_out=''
        count+=3
        if CODON_DICT[s]=="STOP":
            return AA_Seq
    return AA_Seq
def main_functional(DNA_Sequence):
    if valid_base(DNA_Sequence) and not_stop(DNA_Sequence):
        print("Sequence accepted!")
        print(f'Length: {total_length(DNA_Sequence)} bases')
        print(f'Base counts: {base_number(DNA_Sequence)}')
        print(f'GC content: {gc_content_per(DNA_Sequence):.2f}%')
        print(f'Reverse complement: {reverse_complement(DNA_Sequence)}') 
        print(f'Transcription of template strand: {mrna_template(DNA_Sequence)}')
        print(f'Respective amino acid chain: {translate_template(DNA_Sequence)}')
main_functional(DNA_Sequence)