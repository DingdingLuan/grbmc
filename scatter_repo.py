import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd

# Read the data of obs:
obs_z=pd.read_csv('obsresult/z_repo.txt')
obs_lp=pd.read_csv('obsresult/lp_repo.txt')
obs_pea=pd.read_csv('obsresult/pea_repo.txt')
obs_eiso=pd.read_csv('obsresult/eiso_repo.txt')
obs_t90=pd.read_csv('obsresult/t90_repo.txt')
obs_fluence=pd.read_csv('obsresult/fluence_zrepo.txt')
obs_lp=np.log10(obs_lp)
obs_eiso=np.log10(obs_eiso)
obs_pea=np.log10(obs_pea)
obs_fluence=np.log10(obs_fluence)

#Read the data of const model:
const_z=pd.read_csv('result/const/z_repo.txt')
const_lp=pd.read_csv('result/const/lp_repo.txt')
const_pea=pd.read_csv('result/const/pea_repo.txt')
const_eiso=pd.read_csv('result/const/eiso_repo.txt')
const_t90=pd.read_csv('result/const/t90_repo.txt')
const_fluence=pd.read_csv('result/const/fluence_repo.txt')
const_lp=np.log10(const_lp)
const_eiso=np.log10(const_eiso)
const_pea=np.log10(const_pea)
const_fluence=np.log10(const_fluence)
const_parameters=pd.read_csv('result/const/output_parameters.txt')
const_parameters=np.array(const_parameters)
const_chi=const_parameters[0]
const_ec=const_parameters[1]
const_sigma=const_parameters[2]
const_gammac=const_parameters[3]


#Read the data of gauss model:
gauss_z=pd.read_csv('result/gauss/z_repo.txt')
gauss_lp=pd.read_csv('result/gauss/lp_repo.txt')
gauss_pea=pd.read_csv('result/gauss/pea_repo.txt')
gauss_eiso=pd.read_csv('result/gauss/eiso_repo.txt')
gauss_t90=pd.read_csv('result/gauss/t90_repo.txt')
gauss_fluence=pd.read_csv('result/gauss/fluence_repo.txt')
gauss_lp=np.log10(gauss_lp)
gauss_eiso=np.log10(gauss_eiso)
gauss_pea=np.log10(gauss_pea)
gauss_fluence=np.log10(gauss_fluence)
gauss_parameters=pd.read_csv('result/gauss/output_parameters.txt')
gauss_parameters=np.array(gauss_parameters)
gauss_chi=gauss_parameters[0]
gauss_ec=gauss_parameters[1]
gauss_sigma=gauss_parameters[2]
gauss_gammac=gauss_parameters[3]

#Read the data of bpl model:
bpl_z=pd.read_csv('result/bpl/z_repo.txt')
bpl_lp=pd.read_csv('result/bpl/lp_repo.txt')
bpl_pea=pd.read_csv('result/bpl/pea_repo.txt')
bpl_eiso=pd.read_csv('result/bpl/eiso_repo.txt')
bpl_t90=pd.read_csv('result/bpl/t90_repo.txt')
bpl_fluence=pd.read_csv('result/bpl/fluence_repo.txt')
bpl_lp=np.log10(bpl_lp)
bpl_eiso=np.log10(bpl_eiso)
bpl_pea=np.log10(bpl_pea)
bpl_fluence=np.log10(bpl_fluence)
bpl_parameters=pd.read_csv('result/bpl/output_parameters.txt')
bpl_parameters=np.array(bpl_parameters)
bpl_chi=bpl_parameters[0]
bpl_ec=bpl_parameters[1]
bpl_sigma=bpl_parameters[2]
bpl_gammac=bpl_parameters[3]
bpl_b_slope=bpl_parameters[4]



#Plot lp
plt.figure(figsize=(13,9))
plt.title('z-lp')
plt.xlabel('Redshift')
plt.ylabel('lp')
plt.scatter(const_z,const_lp,s=0.07,alpha=1,marker='o',c='r',label='const')
plt.scatter(obs_z,obs_lp,s=30,alpha=1,marker='^',c='black',label='obs')
plt.xlim(0,10)
plt.text(8,49,'$log(E_{\gamma\_c})=%2.2f$'%const_ec,fontsize='14',ha='left',wrap=True)
plt.text(8,48.5,'$\sigma_e=%.4f$'%const_sigma,fontsize='14',ha='left',wrap=True)
plt.text(8,48,'$\Gamma_c==%3.1f$'%const_gammac,fontsize='14',ha='left',wrap=True)
plt.text(8,47.5,'$\chi^2==%.3f$'%const_chi,fontsize='14',ha='left',wrap=True)
plt.legend(loc='best')
plt.savefig('result/plot/z_lp_const.jpg')

plt.figure(figsize=(13,9))
plt.title('z-lp')
plt.xlabel('Redshift')
plt.ylabel('lp')
plt.scatter(gauss_z,gauss_lp,s=0.07,alpha=1,marker='*',c='g',label='gauss')
plt.scatter(obs_z,obs_lp,s=30,alpha=1,marker='^',c='black',label='obs')
plt.text(8,49,'$log(E_{\gamma\_c})=%2.2f$'%gauss_ec,fontsize='14',ha='left',wrap=True)
plt.text(8,48.5,'$\sigma_e=%.4f$'%gauss_sigma,fontsize='14',ha='left',wrap=True)
plt.text(8,48,'$\Gamma_c==%3.1f$'%gauss_gammac,fontsize='14',ha='left',wrap=True)
plt.text(8,47.5,'$\chi^2==%.3f$'%gauss_chi,fontsize='14',ha='left',wrap=True)
plt.xlim(0,10)
plt.legend(loc='best')
plt.savefig('result/plot/z_lp_gauss.jpg')

plt.figure(figsize=(13,9)) 
plt.title('z-lp') 
plt.xlabel('Redshift')
plt.ylabel('lp')
plt.scatter(bpl_z,bpl_lp,s=0.07,alpha=1,marker='o',c='blue',label='bpl')
plt.scatter(obs_z,obs_lp,s=30,alpha=1,marker='^',c='black',label='obs')
plt.text(8,49,'$log(E_{\gamma\_c})=%2.2f$'%bpl_ec,fontsize='14',ha='left',wrap=True)
plt.text(8,48.5,'$\sigma_e=%.4f$'%bpl_sigma,fontsize='14',ha='left',wrap=True)
plt.text(8,48,'$\Gamma_c==%3.1f$'%bpl_gammac,fontsize='14',ha='left',wrap=True)
plt.text(8,47.5,'$b_{slope}==%1.f$'%bpl_b_slope,fontsize='14',ha='left',wrap=True)
plt.text(8,47,'$\chi^2==%.3f$'%bpl_chi,fontsize='14',ha='left',wrap=True)
plt.xlim(0,10)
plt.legend(loc='best')
plt.savefig('result/plot/z_lp_bpl.jpg')






#Plot pea
plt.figure(figsize=(13,9))
plt.title('z-pea')
plt.xlabel('Redshift')
plt.ylabel('pea')
plt.scatter(const_z,const_pea,s=0.07,alpha=1,marker='o',c='r',label='const')
plt.scatter(obs_z,obs_pea,s=30,alpha=1,marker='^',c='black',label='obs')
plt.xlim(0,10)
plt.text(8,6,'$log(E_{\gamma\_c})=%2.2f$'%const_ec,fontsize='14',ha='left',wrap=True)
plt.text(8,5.5,'$\sigma_e=%.4f$'%const_sigma,fontsize='14',ha='left',wrap=True)
plt.text(8,5,'$\Gamma_c==%3.1f$'%const_gammac,fontsize='14',ha='left',wrap=True)
plt.text(8,4.5,'$\chi^2==%.3f$'%const_chi,fontsize='14',ha='left',wrap=True)
plt.legend(loc='best')
plt.savefig('result/plot/z_pea_const.jpg')

plt.figure(figsize=(13,9))
plt.title('z-pea')
plt.xlabel('Redshift')
plt.ylabel('pea')
plt.scatter(gauss_z,gauss_pea,s=0.07,alpha=1,marker='*',c='g',label='gauss')
plt.scatter(obs_z,obs_pea,s=30,alpha=1,marker='^',c='black',label='obs')
plt.xlim(0,10)
plt.text(8,6,'$log(E_{\gamma\_c})=%2.2f$'%gauss_ec,fontsize='14',ha='left',wrap=True)
plt.text(8,5.5,'$\sigma_e=%.4f$'%gauss_sigma,fontsize='14',ha='left',wrap=True)
plt.text(8,5,'$\Gamma_c==%3.1f$'%gauss_gammac,fontsize='14',ha='left',wrap=True)
plt.text(8,4.5,'$\chi^2==%.3f$'%gauss_chi,fontsize='14',ha='left',wrap=True)
plt.legend(loc='best')
plt.savefig('result/plot/z_pea_gauss.jpg')

plt.figure(figsize=(13,9))
plt.title('z-pea')
plt.xlabel('Redshift')
plt.ylabel('pea')
plt.scatter(bpl_z,bpl_pea,s=0.07,alpha=1,marker='o',c='blue',label='bpl')
plt.scatter(obs_z,obs_pea,s=30,alpha=1,marker='^',c='black',label='obs')
plt.xlim(0,10)
plt.text(8,6,'$log(E_{\gamma\_c})=%2.2f$'%bpl_ec,fontsize='14',ha='left',wrap=True)
plt.text(8,5.5,'$\sigma_e=%.4f$'%bpl_sigma,fontsize='14',ha='left',wrap=True)
plt.text(8,5,'$\Gamma_c==%3.1f$'%bpl_gammac,fontsize='14',ha='left',wrap=True)
plt.text(8,4,'$\chi^2==%.3f$'%bpl_chi,fontsize='14',ha='left',wrap=True)
plt.text(8,4.5,'$b_{slope}==%1.f$'%bpl_b_slope,fontsize='14',ha='left',wrap=True)
plt.legend(loc='best')
plt.savefig('result/plot/z_pea_bpl.jpg')




#Plot eiso
plt.figure(figsize=(13,9))
plt.title('z-eiso')
plt.xlabel('Redshift')
plt.ylabel('eiso')
plt.scatter(const_z,const_eiso,s=0.07,alpha=1,marker='o',c='r',label='const')
plt.scatter(obs_z,obs_eiso,s=30,alpha=1,marker='^',c='black',label='obs')
plt.xlim(0,10)
plt.text(8,50,'$log(E_{\gamma\_c})=%2.2f$'%const_ec,fontsize='14',ha='left',wrap=True)
plt.text(8,49.5,'$\sigma_e=%.4f$'%const_sigma,fontsize='14',ha='left',wrap=True)
plt.text(8,49,'$\Gamma_c==%3.1f$'%const_gammac,fontsize='14',ha='left',wrap=True)
plt.text(8,48.5,'$\chi^2==%.3f$'%const_chi,fontsize='14',ha='left',wrap=True)
plt.legend(loc='best')
plt.savefig('result/plot/z_eiso_const.jpg')

plt.figure(figsize=(13,9))
plt.title('z-eiso')
plt.xlabel('Redshift')
plt.ylabel('eiso')
plt.scatter(gauss_z,gauss_eiso,s=0.07,alpha=1,marker='*',c='g',label='gauss')
plt.scatter(obs_z,obs_eiso,s=30,alpha=1,marker='^',c='black',label='obs')
plt.xlim(0,10)
plt.text(8,50,'$log(E_{\gamma\_c})=%2.2f$'%gauss_ec,fontsize='14',ha='left',wrap=True)
plt.text(8,49.5,'$\sigma_e=%.4f$'%gauss_sigma,fontsize='14',ha='left',wrap=True)
plt.text(8,49,'$\Gamma_c==%3.1f$'%gauss_gammac,fontsize='14',ha='left',wrap=True)
plt.text(8,48.5,'$\chi^2==%.3f$'%gauss_chi,fontsize='14',ha='left',wrap=True)
plt.legend(loc='best')
plt.savefig('result/plot/z_eiso_gauss.jpg')

plt.figure(figsize=(13,9))
plt.title('z-eiso')
plt.xlabel('Redshift')
plt.ylabel('eiso')
plt.scatter(bpl_z,bpl_eiso,s=0.07,alpha=1,marker='o',c='blue',label='bpl')
plt.scatter(obs_z,obs_eiso,s=30,alpha=1,marker='^',c='black',label='obs')
plt.xlim(0,10)
plt.text(8,50,'$log(E_{\gamma\_c})=%2.2f$'%bpl_ec,fontsize='14',ha='left',wrap=True)
plt.text(8,49.5,'$\sigma_e=%.4f$'%bpl_sigma,fontsize='14',ha='left',wrap=True)
plt.text(8,49,'$\Gamma_c==%3.1f$'%bpl_gammac,fontsize='14',ha='left',wrap=True)
plt.text(8,48,'$\chi^2==%.3f$'%bpl_chi,fontsize='14',ha='left',wrap=True)
plt.text(8,48.5,'$b_{slope}==%1.f$'%bpl_b_slope,fontsize='14',ha='left',wrap=True)
plt.legend(loc='best')
plt.savefig('result/plot/z_eiso_bpl.jpg')




#Plot fluence
plt.figure(figsize=(13,9))
plt.title('z-fluence')
plt.xlabel('Redshift')
plt.ylabel('fluence')
plt.scatter(const_z,const_fluence,s=0.07,alpha=1,marker='o',c='r',label='const')
plt.scatter(obs_z,obs_fluence,s=30,alpha=1,marker='^',c='black',label='obs')
plt.xlim(0,10)
plt.ylim(-8,0)
plt.text(8,-1.5,'$log(E_{\gamma\_c})=%2.2f$'%const_ec,fontsize='14',ha='left',wrap=True)
plt.text(8,-2,'$\sigma_e=%.4f$'%const_sigma,fontsize='14',ha='left',wrap=True)
plt.text(8,-2.5,'$\Gamma_c==%3.1f$'%const_gammac,fontsize='14',ha='left',wrap=True)
plt.text(8,-3,'$\chi^2==%.3f$'%const_chi,fontsize='14',ha='left',wrap=True)
plt.legend(loc='best')
plt.savefig('result/plot/z_fluence_const.jpg')

plt.figure(figsize=(13,9))
plt.title('z-fluence')
plt.xlabel('Redshift')
plt.ylabel('fluence')
plt.scatter(gauss_z,gauss_fluence,s=0.07,alpha=1,marker='*',c='g',label='gauss')
plt.scatter(obs_z,obs_fluence,s=30,alpha=1,marker='^',c='black',label='obs')
plt.xlim(0,10)
plt.ylim(-8,0)
plt.text(8,-1.5,'$log(E_{\gamma\_c})=%2.2f$'%gauss_ec,fontsize='14',ha='left',wrap=True)
plt.text(8,-2,'$\sigma_e=%.4f$'%gauss_sigma,fontsize='14',ha='left',wrap=True)
plt.text(8,-2.5,'$\Gamma_c==%3.1f$'%gauss_gammac,fontsize='14',ha='left',wrap=True)
plt.text(8,-3,'$\chi^2==%.3f$'%gauss_chi,fontsize='14',ha='left',wrap=True)
plt.legend(loc='best')
plt.savefig('result/plot/z_fluence_gauss.jpg')

plt.figure(figsize=(13,9))
plt.title('z-fluence')
plt.xlabel('Redshift')
plt.ylabel('fluence')
plt.scatter(bpl_z,bpl_fluence,s=0.07,alpha=1,marker='o',c='blue',label='bpl')
plt.scatter(obs_z,obs_fluence,s=30,alpha=1,marker='^',c='black',label='obs')
plt.xlim(0,10)
plt.ylim(-8,0)
plt.text(8,-1.5,'$log(E_{\gamma\_c})=%2.2f$'%bpl_ec,fontsize='14',ha='left',wrap=True)
plt.text(8,-2,'$\sigma_e=%.4f$'%bpl_sigma,fontsize='14',ha='left',wrap=True)
plt.text(8,-2.5,'$\Gamma_c==%3.1f$'%bpl_gammac,fontsize='14',ha='left',wrap=True)
plt.text(8,-3.5,'$\chi^2==%.3f$'%bpl_chi,fontsize='14',ha='left',wrap=True)
plt.text(8,-3,'$b_{slope}==%1.f$'%bpl_b_slope,fontsize='14',ha='left',wrap=True)
plt.legend(loc='best')
plt.savefig('result/plot/z_fluence_bpl.jpg')
