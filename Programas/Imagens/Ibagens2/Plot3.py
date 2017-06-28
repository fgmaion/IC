import numpy as np
import matplotlib.pyplot as plt

print ("Running")
data = "/home/francisco/Documents/Fisica/Mecflu_IC/Programas/Imagens/Ibagens2/myfile10.txt"
f = open(data)

linha = []
linha  = f.readlines()
linha2 = []

col  = 0
conta  = 810

for i in range(len(linha)):
    if (linha[i] == "split\n"):
        col = col + 1
    else:
        linha2.append(linha[i].split(";"))
for i in range(len(linha2)):
    linha2[i][1] = linha2[i][1].replace("\n","")
linha3 = np.zeros((len(linha2),2))
for i in range(len(linha2)):
    linha3[i][0] = float(linha2[i][0])
    linha3[i][1] = float(linha2[i][1])
print(linha3[0][0],"\n")
print("colunas,", col,"\nconta", conta,"\n")
x = np.zeros((col,conta))
y = np.zeros((col,conta))
print((conta-1)+(col-1)*conta)
for i in range(col):
    for j in range(conta):
        x[i][j] = linha3[j+i*conta][0]
        y[i][j] = linha3[j+i*conta][1]

for i in range(col):
    #plt.hist(x[i])
    plt.plot(x[i],y[i],"o")
    plt.ylim([-1,1])
    plt.xlim([0,1.8])
    plt.ylabel('Posicao(y)')
    plt.xlabel('Posicao(x)')
    plt.title('t='+str(i*0.01)+'s')
    plt.axes().set_aspect("equal")
    if(i<10):
        plt.savefig("imagem000"+str(i))
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
