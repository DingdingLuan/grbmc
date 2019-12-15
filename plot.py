import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib

eisoindex=pd.read_csv('result/const/Eiso_index.txt')
eisoindexbpl=pd.read_csv('result/bpl/Eiso_index.txt')
eisoindexgauss=pd.read_csv('result/gauss/Eiso_index.txt')
fluenceindex=pd.read_csv('result/const/fluence_index.txt')
fluenceindexbpl=pd.read_csv('result/bpl/fluence_index.txt')
fluenceindexgauss=pd.read_csv('result/gauss/fluence_index.txt')
#fluenceindex=pd.read_csv('result/const/fluence_index.txt')
#zindex=pd.read_csv('result/const/z_index.txt')
zindex=pd.read_csv('result/const/z_cdf_index.txt')
zindexbpl=pd.read_csv('result/bpl/z_cdf_index.txt')
zindexgauss=pd.read_csv('result/gauss/z_cdf_index.txt')
lpindex=pd.read_csv('result/const/Lp_index.txt')
lpindexbpl=pd.read_csv('result/bpl/Lp_index.txt')
lpindexgauss=pd.read_csv('result/gauss/Lp_index.txt')
peaindex=pd.read_csv('result/const/pea_index.txt')
peaindexbpl=pd.read_csv('result/bpl/pea_index.txt')
peaindexgauss=pd.read_csv('result/gauss/pea_index.txt')
t90index=pd.read_csv('result/const/t90_index.txt')
t90indexbpl=pd.read_csv('result/bpl/t90_index.txt')
t90indexgauss=pd.read_csv('result/gauss/t90_index.txt')

eisoconst=pd.read_csv('result/const/output_Eiso.txt')
fluenceconst=pd.read_csv('result/const/output_fluence.txt')
zconst=pd.read_csv('result/const/output_z.txt')
lpconst=pd.read_csv('result/const/output_Lp.txt')
peaconst=pd.read_csv('result/const/output_pea.txt')
t90const=pd.read_csv('result/const/output_t90.txt')
z_cdf_const=pd.read_csv('result/const/output_z_cdf.txt')

eisogauss=pd.read_csv('result/gauss/output_Eiso.txt')
fluencegauss=pd.read_csv('result/gauss/output_fluence.txt')
zgauss=pd.read_csv('result/gauss/output_z.txt')
lpgauss=pd.read_csv('result/gauss/output_Lp.txt')
peagauss=pd.read_csv('result/gauss/output_pea.txt')
t90gauss=pd.read_csv('result/gauss/output_t90.txt')
z_cdf_gauss=pd.read_csv('result/gauss/output_z_cdf.txt')


eisobpl=pd.read_csv('result/bpl/output_Eiso.txt')
fluencebpl=pd.read_csv('result/bpl/output_fluence.txt')
zbpl=pd.read_csv('result/bpl/output_z.txt')
lpbpl=pd.read_csv('result/bpl/output_Lp.txt')
peabpl=pd.read_csv('result/bpl/output_pea.txt')
t90bpl=pd.read_csv('result/bpl/output_t90.txt')
z_cdf_bpl=pd.read_csv('result/bpl/output_z_cdf.txt')



eisoobs=pd.read_csv('obsresult/eiso.txt',usecols=[1])
eisoobsindex=pd.read_csv('obsresult/eiso.txt',usecols=[0])
f=pd.read_csv('obsresult/fluence_pdf.txt',delimiter=",")
lpobs=pd.read_csv('obsresult/lp.txt',usecols=[1])
pobs=pd.read_csv('obsresult/peakflux.txt',usecols=[1])
pobsindex=pd.read_csv('obsresult/peakflux.txt',usecols=[0])
t90obs=pd.read_csv('obsresult/t90.txt',usecols=[1])
t90obsindex=pd.read_csv('obsresult/t90.txt',usecols=[0])


lposbindex=pd.read_csv('obsresult/lp.txt',usecols=[0])
zobsindex=pd.read_csv('obsresult/z.txt',usecols=[0])
pobsindex=pd.read_csv('obsresult/peakflux.txt',usecols=[0])

zobs_cdf=pd.read_csv('obsresult/z_cdf.txt',usecols=[1])
zobsindex_cdf=pd.read_csv('obsresult/z_cdf.txt',usecols=[0])

szobs_cdf=pd.read_csv('obsresult/sz.txt',usecols=[1])
szobsindex_cdf=pd.read_csv('obsresult/sz.txt',usecols=[0])

# zobss=pd.read_csv('obsresult/z.txt',delimiter=",")
# zobs=zobss['z']
# zobsindex=zobss[' index']
# zbar=zobss['bin']

fluenceobsindex=-f[' index']
fluenceobs=f['s']
fluencebar=f['bin']

# zratioindex=pd.read_csv('result/const/z_ratioindex.txt')
# zconstratio=pd.read_csv('result/const/z_ratio.txt')
# zgaussratio=pd.read_csv('result/gauss/z_ratio.txt')
# zplratio=pd.read_csv('result/bpl/z_ratio.txt')




matplotlib.rcParams['xtick.direction']='in'
matplotlib.rcParams['ytick.direction']='in'

# plt.figure(figsize=(9,9))
# plt.title('z_ratio')
# plt.plot(zratioindex,zconstratio,c='r',label='const')
# plt.plot(zratioindex,zgaussratio,c='g',label='gauss')
# plt.plot(eisoindex,eisobpl,c='b',label='bpl')
# plt.plot(zratioindex,zplratio,c='black',label='obs')
# plt.legend(loc='best')
# plt.savefig('result/plot/z_ratio.eps')

#plot Eiso:
plt.figure(figsize=(10,9))
plt.title('Eiso compared')
plt.step(eisoobsindex,eisoobs,linewidth=1.3,linestyle='-',c='black',label='obs')
plt.step(eisoindex,eisoconst,linewidth=1.1,linestyle='--',c='r',label='const')
plt.step(eisoindexgauss,eisogauss,linewidth=1.1,linestyle='--',c='g',label='gauss')
plt.step(eisoindexbpl,eisobpl,linewidth=1.1,linestyle='--',c='b',label='bpl')
plt.xlim(48,56)
plt.ylim(0,0.8)
plt.legend(loc='best')
plt.savefig('result/plot/eiso_step.jpg')
# plt.show()

# plt.figure(figsize=(9,9))
# plt.title('Eiso compared')
# plt.plot(eisoindex,eisoconst,c='r',label='const')
# plt.plot(eisoindex,eisogauss,c='g',label='gauss')
# # plt.plot(eisoindex,eisobpl,c='b',label='bpl')
# plt.plot(eisoindex,eisoobs,c='black',label='obs')
# plt.legend(loc='best')
# plt.savefig('result/plot/eiso_line.eps')

# #plot t90:
plt.figure(figsize=(9,9))
plt.title('t90 compared')
plt.step(t90obsindex,t90obs,linewidth=1.3,linestyle='-',c='black',label='obs')
plt.step(t90index,t90const,linewidth=2.7,linestyle='--',c='r',label='const')
plt.step(t90indexgauss,t90gauss,linewidth=2.3,linestyle='--',c='g',label='gauss')
plt.step(t90indexbpl,t90bpl,linewidth=0.5,linestyle='--',c='b',label='bpl')
plt.legend(loc='best')
plt.savefig('result/plot/t90_step.jpg')

# fluencebar=[]
# for i in range(20):
# 	fluencebar=np.append(fluencebar,0.125)

#plot fluence

a=[]
b=[]
c=[]
for i in range(len(fluenceobsindex)):
	if fluenceobs[i]!=0:
		a=np.append(a,fluenceobsindex[i])
		b=np.append(b,fluenceobs[i])
		c=np.append(c,fluencebar[i])

fluenceobsindex=a
fluenceobs=b
fluencebar=c

plt.figure(figsize=(13 ,9))
plt.title('fluence compared')
plt.errorbar(fluenceobsindex,fluenceobs,xerr=fluencebar,fmt='s',elinewidth=1.3,capsize=3,color='black',label='obs')
plt.plot(fluenceindex,fluenceconst,linewidth=1.1,linestyle='--',c='r',label='const')
plt.plot(fluenceindexgauss,fluencegauss,linewidth=1.1,linestyle='-',c='g',label='gauss')
plt.plot(fluenceindexbpl,fluencebpl,linewidth=1.1,linestyle='--',c='b',label='bpl')
plt.xlim((-8,-3))
plt.ylim((0,0.8))
# plt.errorbar(fluenceobsindex,fluenceobs,fmt='-o',elinewidth=0.8,capsize=3,color='black',label='obs')
# plt.step(fluenceobsindex,fluenceobs,linewidth=2.5,linestyle='-',c='black',label='obs')
plt.legend(loc='best')
plt.savefig('result/plot/fluence_step.jpg')
#plt.show()


# print(fluenceindex,fluenceobs)
# plt.figure(figsize=(9,9))
# plt.title('fluence compared')
# plt.plot(fluenceindex,fluenceconst,c='r',label='const')
# plt.plot(fluenceindex,fluencegauss,c='g',label='gauss')
# # plt.plot(fluenceindex,fluencebpl,c='b',label='bpl')
# plt.plot(fluenceindex,fluenceobs,c='black',label='obs')
# plt.legend(loc='best')
# plt.savefig('result/plot/fluence_line.eps')

#####################################################
# a=[]
# b=[]
# c=[]
# for i in range(12):
# 	if zobs[i]!=0:
# 		a=np.append(a,zobsindex[i])
# 		b=np.append(b,zobs[i])
# 		c=np.append(c,zbar[i])

# zobsindex=a
# zobs=b
# zbar=c

#plot z
# plt.figure(figsize=(13,9))
# plt.title('z compared')
# plt.plot(zindex,zconst,linewidth=1.1,linestyle='--',c='r',label='const')
# plt.plot(zindex,zgauss,linewidth=1.1,linestyle='--',c='g',label='gauss')
# plt.plot(zindex,zbpl,linewidth=1.1,linestyle='--',c='b',label='bpl')
# plt.errorbar(zobsindex,zobs,xerr=zbar,fmt='s',elinewidth=1.3,capsize=3,color='black',label='obs')
# plt.legend(loc='best')
# plt.savefig('result/plot/z.jpg')

#plot z_cdf
plt.figure(figsize=(13,9))
plt.title('z compared')
plt.xlim(0,8)
plt.ylim(0,1)
#plt.step(szobsindex_cdf,1-szobs_cdf,linewidth=20.3,linestyle='-',c='yellow',label='s_obs')
plt.step(zobsindex_cdf,1-zobs_cdf,linewidth=1.3,linestyle='-',c='black',label='obs')
plt.plot(zindex,1-z_cdf_const,linewidth=1.1,linestyle='--',c='r',label='const')
plt.plot(zindexgauss,1-z_cdf_gauss,linewidth=1.1,linestyle='--',c='g',label='gauss')
plt.plot(zindexbpl,1-z_cdf_bpl,linewidth=1.1,linestyle='--',c='b',label='bpl')
plt.legend(loc='best')
plt.savefig('result/plot/z_cdf.jpg')

# plt.figure(figsize=(9,9))
# plt.title('z compared')
# plt.plot(zindex,zobs,c='black',label='obs')
# plt.plot(zindex,zconst,c='r',label='const')
# plt.plot(zindex,zgauss,c='g',label='gauss')
# # plt.plot(zindex,zbpl,c='b',label='bpl')
# plt.legend(loc='best')
# plt.savefig('result/plot/z_line.eps')
###########################################################################
#plot Lp
#plot Lp oohohhh yes!
plt.figure(figsize=(9,9))
plt.title('Lp compared')
plt.step(lposbindex,lpobs,linewidth=1.3,linestyle='-',c='black',label='obs')
plt.plot(lpindex,lpconst,linewidth=1.1,linestyle='-',c='r',label='const')
plt.plot(lpindexgauss,lpgauss,linewidth=1.1,linestyle='--',c='g',label='gauss')
plt.plot(lpindexbpl,lpbpl,linewidth=1.1,linestyle='--',c='b',label='bpl')
plt.legend(loc='best')
plt.savefig('result/plot/lp_step.jpg')

# plt.figure(figsize=(9,9))
# plt.title('Lp compared')
# plt.plot(lpindex,lpconst,c='r',label='const')
# plt.plot(lpindex,lpgauss,c='g',label='gauss')
# # plt.plot(lpindex,lpbpl,c='b',label='bpl')
# plt.plot(lpindex,lpobs,c='black',label='obs')
# plt.legend(loc='best')
# plt.savefig('result/plot/lp_line.eps')

# # #plot peak_flux
# plt.figure(figsize=(9,9))
# plt.title('peak_flux compared')
# plt.plot(peaindex,peaconst,c='r',label='const')
# plt.plot(peaindex,peagauss,c='g',label='gauss')
# # plt.plot(peaindex,peabpl,c='b',label='bpl')
# plt.plot(peaindex,pobs,c='black',label='obs')
# plt.legend(loc='best')
# plt.savefig('result/plot/peak_flux_line.eps')

# plt.figure(figsize=(13,9))
# plt.title('peak_flux compared')
# plt.plot(peaindex,peaconst,linewidth=1.1,linestyle='-',c='r',label='const')
# plt.plot(peaindex,peagauss,linewidth=1.1,linestyle='--',c='g',label='gauss')
# plt.plot(peaindex,peabpl,linewidth=1.1,linestyle='--',c='b',label='bpl')
# plt.step(pobsindex,pobs,linewidth=1.3,linestyle='-',c='black',label='obs')
# plt.legend(loc='best')
# plt.savefig('result/plot/peak_flux_step.jpg')
