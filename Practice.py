Valid_DNA=["A", "T", "G", "C"]
DNA_Sequence=input()
if DNA_Sequence!="STOP":
    print("Sequence accepted!")
    DNA_Sequence=DNA_Sequence.upper()
    def Test(DNA):
        x=True
        for i in DNA_Sequence:
            if i not in Valid_DNA:
                return False
        return x
    if Test(DNA_Sequence):
        def Total_Length(DNA):
            Length=len(DNA)
            return Length
        def Base_Number(DNA):
            Base_Counts={"A":0,"T":0,"G":0,"C":0}
            for Base in DNA:
                Base_Counts[Base]+=1
            return Base_Counts
        def GC_Content_Per(DNA):
            S = Base_Number(DNA)
            C_Content=S["C"]
            G_Content=S["G"]
            GC_Hundo=(((C_Content+G_Content)/Total_Length(DNA))*100)
            return GC_Hundo
        def Reverse_Complement(DNA):
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
        def main(DNA):
            print(f'Length: {Total_Length(DNA)} bases')
            print(f'Base counts: {Base_Number(DNA):}')
            print(f'GC content: {GC_Content_Per(DNA):.2f}%')
            print(f'Reverse complement: {Reverse_Complement(DNA)}')
        main(DNA_Sequence)
else:
    print("Program terminated.")
