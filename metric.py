import numpy as np
from Bio import SeqIO

def Kmer(sequence, K):
    m=4**K
    na_vect=[0]*(3*m)
    pos_sum=[0]*m
    squa_sum=[0]*m
    n=len(sequence)-(K-1)
    index_map = {  'a':0, 'A':0, 'c':1, 'C':1, 'g':2, 'G':2, 't':3, 'T':3  }
    for i in range(0, n):
        flag=1
        for l in range(0,K):
            if sequence[i+l] not in index_map.keys():
                flag=0
        if flag == 0:
            continue
        tem=index_map[sequence[i]]
        for l in range(1,K):
            tem=4*tem+index_map[sequence[i+l]]
        na_vect[tem] += 1
        pos_sum[tem] += i+1
    for k in range(0,m):
        if na_vect[k] != 0:
            na_vect[k+m] = pos_sum[k] / na_vect[k]
        else:
            na_vect[k+m]=0
    for i in range(0, n):
        flag=1
        for l in range(0,K):
            if sequence[i+l] not in index_map.keys():
                flag=0
        if flag == 0:
            continue
        tem=index_map[sequence[i]]
        for l in range(1,K):
            tem=4*tem+index_map[sequence[i+l]]
        squa_sum[tem] += ( i + 1 - na_vect[tem+m] ) ** 2
    for k in range(0,m):
        if na_vect[k] != 0:
            na_vect[k+2*m] = squa_sum[k] / (n * na_vect[k])
        else:
            na_vect[k+2*m]=0
    return na_vect

def distance(X):
  n,m = X.shape
  G = np.dot(X,X.T)
  H = np.tile(np.diag(G), (n,1))
  return np.sqrt(H + H.T - 2*G)

Fasta_Filename="SARSCoV2.fasta"

seq=[]
for read in SeqIO.parse(Fasta_Filename, "fasta"):
    seq.append((read.seq).__str__())
N=len(seq)

for K in range(1,10):
    nv=np.zeros((N,3*(4**K)))
    for j in range(N):
        nv[j,:]=Kmer(seq[j],K)
    for i in range(0,3):
        nvO=nv[:,i*(4**K):(i+1)*(4**K)]
        disO=distance(nv)
        np.save("order_dis\\"+str(K)+"mer_order"+str(i)+"_distance.npy", disO)

dis=[]
weight=[0.03549459364589514, 0.039082434930630405, 0.30412906453093275, 0.030801860439204753, 0.013913759681552054, 0.0847541874627664, 0.2147908821824287, 0.009746970476623207, 0.037468677500208866, 0.4921144274939822, 0.0017718520847580079, 0.0025613627212840574, 0.7629788655456013, 0.0012773638058877917, 0.005374227624313233, 1.0, 0.0017544541737680417, 0.001553408700637301, 0.8418668476692482, 0.0010680411101253774, 0.00010663595705133611, 0.713477302936144, 0.0010440696319391177, 0.00092914523021887, 0.1843582232771541, 0.000813690120208786, 5.0762169732243335e-06]
for K in range(1,10):
    for i in range(0,3):
        dis.append(np.load("order_dis\\"+str(K)+"mer_order"+str(i)+"_distance.npy", allow_pickle=True))
disfinal=weight[0]*dis[0]
for i in range(1,len(weight)):
    disfinal+=weight[i]*dis[i]
np.save("Optimal_distance.npy", disfinal)
    






