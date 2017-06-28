import numpy as np
import matplotlib.pyplot as plt
print("Rodando")
data='/home/francisco/myfile6.txt'
f=open(data)
linha=[]
grosso=[]
linha = f.readlines()
L=len(linha)
print("numero de linhas", L)
col=0
conta=0
tabela=[]
for i in range(len(linha)):

    if linha[i]=='split\n':
        col=col+1
        conta=len(linha)-1-i
    else:
        grosso.append(linha[i].split(';'))

print ( "valor de conta ", conta)

for i in range(len(grosso)):
    grosso[i][1]=grosso[i][1].replace('\n','')
#for i in range(len(linha)-len(linha)/conta):
x=np.zeros((col,conta))
y=np.zeros((col,conta))
tudogrosso=np.zeros((len(grosso),2))

for i in range(len(grosso)):
    tudogrosso[i][0]=float(grosso[i][0])
    tudogrosso[i][1]=float(grosso[i][1])
m=0
n=0
for i in range(col):
    for j in range(conta):

        x[i][j]=tudogrosso[j+i*conta][0]
        y[i][j]=tudogrosso[j+i*conta][1]
print(col)


plt.ylim([0,3])
plt.xlim([-1,1])
plt.ylabel('Densidade (Rho)')
plt.xlabel('x')
plt.title('t='+str(1499*0.01)+'s')
plt.plot(x[1499],y[1499],'o')
#plt.hist(x[29])
plt.savefig('t=3000')
plt.clf()
for i in range(col):
    #plt.hist(x[i])
    plt.ylim([0,4])
    plt.xlim([-2,2])
    plt.ylabel('Densidade (Rho)')
    plt.xlabel('x')
    plt.title('t='+str(i*0.01)+'s')
    plt.plot(x[i],y[i],'o')
    if (i<10):
        plt.savefig('imagem000'+str(i))
    else:
        if (i<100):
            plt.savefig('imagem00'+str(i))
        else:
            if (i<1000):
                plt.savefig('imagem0'+str(i))
            else:
                if(i<10000):
                    plt.savefig('imagem'+str(i))

    plt.clf()



f.close()
