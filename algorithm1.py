import numpy as np

dis = np.load('Optimal_distance.npy', allow_pickle=True)
n,_=dis.shape
dis=dis+np.diag([float("inf")] * n)
Csv_Filename="family.npy"
family = np.load(Csv_Filename, allow_pickle=True)
Uniquefamily= np.unique(family)
faindex={}
for fa in Uniquefamily:
    faindex[fa]=[]
for i in range(len(family)):
    faindex[family[i]].append(i)

MaxL=0
for key in faindex.keys():
    if MaxL<len(faindex[key]):
        MaxL=len(faindex[key])
kkk=int(MaxL/1000)+1
print(kkk)

def mindis(i, fa):
    vec = dis[i, faindex[fa]]
    return np.sort(vec)[kkk-1]  

def distribution(fa1, fa2): 
    vec=[]
    for i in faindex[fa2]:
        vec.append(mindis(i,fa1))
    return np.array(vec)

def distribution_in():
    vec=[]
    for key in faindex.keys():
        vec.extend(distribution(key,key))
    return vec

def distribution_out():
    vec=[]
    for key1 in faindex.keys():
        for key2 in faindex.keys():
            if key1 != key2:
                vec.extend(distribution(key1,key2))
    return vec

def nearest_neighbor_classification(dist_matrix, family):
    ALL=len(family)
    ACC=0
    for i in range(ALL):
        nearest_neighbor_index = np.argmin(dist_matrix[i,:])
        nearest_neighbor_family = family[nearest_neighbor_index]
        if nearest_neighbor_family==family[i]:
            ACC+=1
    return float(ACC)/float(ALL)

ACC1=0
ACC2=0
ALL=0
vec=distribution_in()
bound=np.percentile(vec,99)
for name1 in faindex.keys():
    acc1=0
    acc2=0
    all=len(faindex[name1])
    for i in faindex[name1]:
        flag1=0
        flag2=0
        for name2 in faindex.keys():
            if name2==name1:
                continue
            elif mindis(i,name2)<bound:
                flag1=1
        for name2 in faindex.keys():
            if mindis(i,name2)<bound:
                flag2=1
        if flag1==0:
            acc1+=1
        if flag2==1:
            acc2+=1
    ACC1+=acc1
    ACC2+=acc2
    ALL+=all
print("{:.2f}".format(100*(1-(ACC1/ALL)))+"\%","{:.2f}".format(100*(1-(ACC2/ALL)))+"\%")
