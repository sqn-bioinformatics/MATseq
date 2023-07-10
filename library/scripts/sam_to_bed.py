#loads sam and removes the header
inputsam=open(snakemake.input[0]).readlines()
inputsam=inputsam[93337:]

#merges Fw and Rv reads
outputfile=[]
for i in range(int(len(inputsam)/2)):
    a=inputsam[i*2].split('\t')
    b=inputsam[i*2+1].split('\t')
    if a[0]!=b[0]:
        print('error reads pair not matching')
    if int(a[8])>0:
        temp=a[2]+'\t'+a[3]+'\t'+str(int(a[3])+int(a[8]))+'\t'+a[0]+'\t+\n'
    elif int(a[8])<0:
        temp=a[2]+'\t'+b[3]+'\t'+str(int(b[3])+int(b[8]))+'\t'+a[0]+'\t-\n'
    else:
        print('error selecting strand')
    outputfile.append(temp)
open(snakemake.output[0],'w').writelines(outputfile)
print(snakemake.output[0],' done!')