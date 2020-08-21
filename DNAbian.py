import numpy as np

#DNA-Encoding RULE #1 A = 00, T=11, G=01, C=10
dna1 = {}
dna1["00"] = "A"
dna1["01"] = "G"
dna1["10"] = "C"
dna1["11"] = "T"
dna1["A"] = [0,0]
dna1["G"] = [0,1]
dna1["C"] = [1,0]
dna1["T"] = [1,1]

#DNA xor
dna1["AA"]=dna1["TT"]=dna1["GG"]=dna1["CC"]="A"
dna1["AG"]=dna1["GA"]=dna1["TC"]=dna1["CT"]="G"
dna1["AC"]=dna1["CA"]=dna1["GT"]=dna1["TG"]="C"
dna1["AT"]=dna1["TA"]=dna1["CG"]=dna1["GC"]="T"

dna2 = {}
dna2["00"] = "A"
dna2["01"] = "C"
dna2["10"] = "G"
dna2["11"] = "T"
dna2["A"] = [0,0]
dna2["C"] = [0,1]
dna2["G"] = [1,0]
dna2["T"] = [1,1]

dna3 = {}
dna3["00"] = "G"
dna3["01"] = "A"
dna3["10"] = "T"
dna3["11"] = "C"
dna3["G"] = [0,0]
dna3["A"] = [0,1]
dna3["T"] = [1,0]
dna3["C"] = [1,1]

dna4 = {}
dna4["00"] = "C"
dna4["01"] = "A"
dna4["10"] = "T"
dna4["11"] = "G"
dna4["C"] = [0,0]
dna4["A"] = [0,1]
dna4["T"] = [1,0]
dna4["G"] = [1,1]

dna5 = {}
dna5["00"] = "G"
dna5["01"] = "T"
dna5["10"] = "A"
dna5["11"] = "C"
dna5["G"] = [0,0]
dna5["T"] = [0,1]
dna5["A"] = [1,0]
dna5["C"] = [1,1]

dna6 = {}
dna6["00"] = "C"
dna6["01"] = "T"
dna6["10"] = "A"
dna6["11"] = "G"
dna6["C"] = [0,0]
dna6["T"] = [0,1]
dna6["A"] = [1,0]
dna6["G"] = [1,1]

dna7 = {}
dna7["00"] = "T"
dna7["01"] = "G"
dna7["10"] = "C"
dna7["11"] = "A"
dna7["T"] = [0,0]
dna7["G"] = [0,1]
dna7["C"] = [1,0]
dna7["A"] = [1,1]

dna8 = {}
dna8["00"] = "T"
dna8["01"] = "C"
dna8["10"] = "G"
dna8["11"] = "A"
dna8["T"] = [0,0]
dna8["C"] = [0,1]
dna8["G"] = [1,0]
dna8["A"] = [1,1]

def dna_encode(b, g, r):
    b = np.unpackbits(b, axis=1)
    g = np.unpackbits(g, axis=1)
    r = np.unpackbits(r, axis=1)
    m, n = b.shape
    r_enc = np.chararray((m, int(n / 2)))
    g_enc = np.chararray((m, int(n / 2)))
    b_enc = np.chararray((m, int(n / 2)))

    for color, enc in zip((b, g, r), (b_enc, g_enc, r_enc)):
        idx = 0
        for j in range(0, m):
            for i in range(0, n, 2):
                enc[j, idx] = dna1["{0}{1}".format(color[j, i], color[j, i + 1])]
                idx += 1
                if (i == n - 2):
                    idx = 0
                    break

    b_enc = b_enc.astype(str)
    g_enc = g_enc.astype(str)
    r_enc = r_enc.astype(str)
    return b_enc, g_enc, r_enc