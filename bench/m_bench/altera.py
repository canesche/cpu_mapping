def f1():
    f = open('./bench/m_bench/txt')
    h = open('./bench/m_bench/txxt')
    st = f.readlines()
    ht = h.readlines()
    for i in range(len(ht)):
        ht[i] = ht[i].replace(" ","")
        ht[i] = ht[i].replace("\n","")

    for j in range(len(ht)):
        if ht[j] == "" : continue
        for i in range(len(st)):
            if ht[j]+" " in st[i]:
                st[i] = st[i].replace(ht[j], str(j))
    f.close()
    f = open('./bench/m_bench/txt','w')
    for s in st:
        f.write(s)
    h.close()
    h = open('./bench/m_bench/txxt','w')
    counter = 0
    for s in ht:
        h.write("%d [label=%s op=]\n"%(counter,s))
        counter+=1
    a = 1

def f2():
    f = open('./bench/m_bench/txt')
    st = f.readlines()
    for j in range(len(st)):
        st[j] = st[j].replace('\n','')
        st[j] = st[j].replace(' ','')
        st[j] = st[j].split("->")
    
    for i in range(len(st)):
        for j in range(i+1,len(st)):
            if int(st[i][1]) > int(st[j][1]):
                st[i], st[j] = st[j] , st[i]
    f.close()
    
    wstr = ""
    b = None
    for s in st:
        if b is not None:
            if b == s[1]:
                wstr+=("%s -> %s [port=1 weight=0]\n"%(s[0],s[1]))    
            else: 
                wstr+=("%s -> %s [port=0 weight=0]\n"%(s[0],s[1]))    
            b = s[1]    
        else:
            wstr+=("%s -> %s [port=0 weight=0]\n"%(s[0],s[1]))
            b = s[1]

    f = open('./bench/m_bench/txt','w')
    f.write(wstr)
    f.close()
   
    a = 1

f2()