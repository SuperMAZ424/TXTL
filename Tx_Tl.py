import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from math import exp
import os


Kt_IPTG = 0.92                                                           #min-1
kdr1 = 3e-7                                                              #nM-2.min-1
k_dr1 = 12                                                               #min-1
k2R=	50	                                                             #nM−1·min−1	LacI dimerization rate constant‡
k_2R=	10**-3	                                                         #min−1	LacI dimer dissociation rate constant‡
ksMR=	0.23	                                                         #nM·min−1	lacI transcription rate†
nMR =	0.462	                                                         #min−1	lacI mRNA degradation constant†
nR=	0.2	                                                                 #min−1	repressor monomer degradation constant†
ksR=	15	                                                             #min−1	LacI monomer translation rate constant†
nR2=	0.2	                                                             #min−1	LacI  dimer degradation constant†
kdr2=	3e-7	                                                         #nM−2·min−1 association rate constant for 2nd derepression mechanism
k_dr2=	4.8e3	                                                         #nM−1·min−1 dissociation rate constant for 2nd derepression mechanism
kr=	960	                                                                 #nM−1·min−1 association rate constant for repression
k_r=	2.4	                                                             #min−1	dissociation rate constant for repression§
nI2R2=	0.2	                                                             #min−1	repressor-inducer degradation constant†
ktx = 2400/int(900)#input('enter number of basepairs of the part DNA: '))     #min-1 transcription rate
delta_mRNA = 0.462                                                       #min-1 mRNA degradation rate
ktl = 480/int(274)#input('protein sequence length: '))                        #min-1 translation rate
delta_prot= 0.2                                                         #min-1 protein degradation rate
O0 = str(500) #input('enter plasmid copy number, if multiple separate by space: ') #plasmid copy-number
O0 = O0.split(" ")                                                            #make list of the plasmid copy number
iptgint = str(10e5) #input('Enter the IPTG concentration in nM, if multiple separate by space: ')         #initial IPTG concentration(s)
iptgint = iptgint.split(" ")                                                  #if multiple iptg concentrations are input it will be converted to list, they must be separated by space
path_of_the_directory= "D:\Igem\igem meetings\\2022_10_7\clt_\CLT" #input('enter the thermodynamics parameters files copied from mfold directory: ')
D_gibbs_of_hyb = str(-12.5) #input('Enter the values of energy of hybridization: ')
path_of_the_directory += '\\'
name_list = ''
energy_list = ''
sec_energy = ''
base_number = ''

"""**define all constants, enter all inputs, post input processing, applying the defined functions, and print the outputs**"""

Rt = 57000                                                                   #constant input('Enter the total number of ribosomes in a cell (Rt): ')
S = 20                                                                       #constant input('Enter the number of ribosomes per polysome (S): ')
#Tmrna = 5700                                                                 #constant input('Enter the total numbers of mRNAs in the cell (TmRNA): ')
Gas_constant = 1.99*(10**(-3))                                               #kCal/mol.K constant
timef = 1200 #input('enter the final time for simulation in minutes: ')            #the final time for simulation in minutes
times = 1200 #input('the number of time steps: ')                                  #number of time steps


for filename in os.listdir(path_of_the_directory):
    if filename.endswith(".txt"):
        z = path_of_the_directory + filename
        file = open(z, 'r')
        file_list = file.readlines()
# file= open('C:\\Users\\Mostafa Zaki\\Desktop\\2022_18_9\\str1.txt', 'r')

        file.close()
    counter = 1

    for n in file_list:
        if counter == 1:
            if (len(name_list) <1):
                name_list += n.replace('\n','')
            else:
                name_list += ',' + n.replace('\n','')
        if counter == 3:
            if (len(energy_list) <1) or (energy_list[-1] == '/'):
                energy_list += n[6:].replace('\n','')
            else:
                energy_list += ',' + n[6:].replace('\n','')

        if n[0:5] == "Helix":
            if (len(sec_energy) <1) or (sec_energy[-1] == '/'):
                sec_energy += n[6:12].replace('\t','')
            else:
                sec_energy += ',' + n[6:12].replace('\t', '')
            if (len(base_number) < 1) or (base_number[-1] == '/'):
                m= n[12:16].replace('\t', "")
                m= m.replace('b', '')
                m= m.replace('a','')
                base_number += m
            else:
                m = ',' + n[12:16].replace('\t', "")
                m = m.replace('b', '')
                m = m.replace('a', '')
                base_number += m

        counter += 1
    name_list += '/'
    sec_energy += '/'
    base_number += '/'



print(name_list)
print(energy_list)
print(sec_energy)
print(type(sec_energy))
base_number = base_number.replace(' ', '')
print(base_number)
print(type(base_number))




"""define functions for calculating required parameters"""

def P_Si():
  global temperature, D_Gibbs_list, Gas_constant
  summation = 0
  for i in D_Gibbs_list:
    m = exp(-(float(i)/(float(Gas_constant)*float(temperature))))
    summation += m
  p_si_list = []
  for n in D_Gibbs_list:
    probability = (exp(-(float(n)/(float(Gas_constant)*float(temperature)))))/ float(summation)
    p_si_list.append(probability)
  return(p_si_list)

"""calculate the theta and probability of unpairing """

def theta_calculator():
  global Gas_constant, temperature, DG_SS_list, L_list
  x = 0
  Pi_list = []
  for x in range(len(DG_SS_list)):

      DG_list2 = DG_SS_list[x].split(',')
      L_list2 = L_list[x].split(',')

      #print(DG_list2)
      #print(L_list2)
      if len(DG_list2) == len(L_list2):
          s = 0
          Pi_list2 = []
          for s in range(len(L_list2)):
              #print('s is', s)
              gibbs_h_list = [float(DG_list2[s])] * int(L_list2[s])
              #print('gibbs_h_list is ', gibbs_h_list)
              theta_list = []
              for i in gibbs_h_list:
                  theta = ((1 + exp(-(float(i) / (float(Gas_constant) * float(temperature))))) ** (-1)) ** (
                              1 / float(L_list2[s]))
                  theta_list.append(theta)
              #print(theta_list)
              probability_of_unpairing = 1
              for b in theta_list:
                  probability_of_unpairing = float(probability_of_unpairing) * float(b)
              Pi_list2.append(probability_of_unpairing)
              s += 1
          if len(Pi_list2) > 1:
              b = 1
              for o in Pi_list2:
                  b = b * o
              Pi_list.append(b)
          else:
              Pi_list.append(Pi_list2[0])
          #print(Pi_list)
  return(Pi_list)

"""calculation of probability of exposure"""

def p_ex(p_si_list, probability_of_unpairing):
  if len(p_si_list) == len(probability_of_unpairing):
    pexp = 0
    q = 0
    for q in range(len(p_si_list)):
      pexp += float(p_si_list[q]) * float(probability_of_unpairing[q])
      return(pexp)
  else:
    return('two lists are not equal')

"""calculate Kr (hybridization constant)"""

def Kr():
  global DGh_list
  Kr_list = []
  for u in DGh_list:
    Kru = exp((-(float(u)))/(float(Gas_constant)*float(temperature)))
    Kr_list.append(Kru)
  return(Kr_list)

"""calculate alpha 

α = 1 + KR . Pex . (Rt/S) + Kr . Pex . Tmrna
"""

def alpha(pex, Kr_list, Tmrna):
  global Rt, S
  alpha_list = []
  for n in Kr_list:
    alpha = 1 + (n * pex * (Rt/S)) + (n * pex * Tmrna)
    alpha_list.append(alpha)
  return(alpha_list)

"""calculate Pc"""

def Pc_calc(alpha, pex, Kr_list, Tmrna):
  global Rt, S
  h = 0
  Pc_list = []
  if Tmrna == 0:
      Pc = 0
      Pc_list.append(Pc)
  else:
    for h in range(len(Kr_list)):
      nal = float(alpha[h])
      kr = float(Kr_list[h])
      Pc = ((nal)-(((nal**2)-(4 * (kr)**(2)) * ((pex)**(2)) * (Rt/S) * (Tmrna)))**(0.5)) / (2 * (kr) * (pex) * (Tmrna))
      Pc_list.append(Pc)
  return(Pc_list)



#############################################################################################################################################

temperature = '298.5' #float(input('Enter the working temperature: '))
#D_Gibbs_of_folding = '-3.30,-3.00,-2.80,-2.60,-2.50,-2.30' #input('Enter all energies seperated by comma (,): ')

#DG_sec_struc = '-5.70/-5.70,-3.80/-5.70,-4.20/-5.70,-4.50,-1.90,-2.20/-5.70,-1.30/-4.30,-2.20' #'''input('enter the energy of each stack (if more than 1 stack in the same secundary structure separate (,)) (for different secondary structure separate by (/)) : )'''
#L = '3/3,5/3,5/3,3,3,3/3,2/5,3' #input('enter the number of nucleutides in stack form separated by comma (,): ')

#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################

D_Gibbs_list = energy_list.split(',')
DGh_list = D_gibbs_of_hyb.split(',')
DG_SS_list = sec_energy[:-1].split('/')
L_list = base_number[:-1].split('/')

#min_DGh = min(DGh_list)
#print('minimum energy of hybridization is: ', min_DGh)
'''
print('DGh_list is ',DGh_list)
print('D_Gibbs_list is ',D_Gibbs_list)
print('DG_SS_list is ', DG_SS_list)
print('L_list is ', L_list)
'''

#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################

A = P_Si()
#print('the probability of folds is: ', A)
B = theta_calculator()
#print('the probability of unpairing is: ',B)
C = p_ex(A,B)
#print('the probability of exposure is: ', C)
D = Kr()
#print('list of Kr constants is: ', D)
#E = alpha( C, D)
#print('alpha is',E)
#F = Pc_calc(E, C, D)
#print('Pc_calc is: ', F)


def central_dogma(x, t):
    global C, D
    iptgex = x[0]
    Iptgin = x[1]
    laci = x[2]
    laci2_Iin = x[3]
    laci2 = x[4]
    Mlaci = x[5]
    R2O = x[6]          ##########################
    O = x[7]         ###########################
    mRNA = x[8]
    prot = x[9]

    E = alpha(C, D, Tmrna=mRNA)
    F = Pc_calc(E, C, D, Tmrna=mRNA)

    rxn1_f = Kt_IPTG*iptgex   #diffusion of iotg into the cell
    rxn1_b = Kt_IPTG*Iptgin   #diffusion of iotg out of the cell
    rxn2_f = kdr1*laci2*(Iptgin**2) #formation of laci2-iptg complex
    rxn2_b = k_dr1*laci2_Iin   #dissociation of laci2-iptg complex
    rxn3_f = k2R*(laci**2)         #formation of laci dimer
    rxn3_b = k_2R*laci2       #dissociation of laci dimer
    rxn4 = ksMR               #formation of laci mRNA
    rxn5 = nMR*Mlaci          #degredation of laci mRNA
    rxn7 = ksR * Mlaci        #laci translation
    rxn6 = nR*laci            #laci monomer degradation
    rxn8 = nR2*laci2          #laci2 degradation
    rxn9_f = kdr2*R2O*(Iptgin**2)   #second derepresion forward
    rxn9_b = k_dr2*laci2_Iin*O      #second derepresion backward
    rxn10_f = kr*laci2*O            #association of repression
    rxn10_b = k_r*R2O               #dissociation of repression
    rxn11 = nI2R2*laci2_Iin         #degredation of laci dimer-iptg complex
    rxn12_l = 0.001*ktx*R2O         #leakage of lacI2-lacO
    rxn12 = ktx*O                   #mRNA transcription
    rxn13 = nMR*mRNA         #mRNA degradation
    rxn14 = float(ktl)*float(F[0])*mRNA              #translation
    rxn15 = delta_prot*prot         #protein degradation

    diptgexdt = -rxn1_f + rxn1_b                                              #2 rxns
    diptgindt = +rxn1_f - rxn2_f - rxn1_b + rxn2_b -rxn9_f + rxn9_b           #6 rxns
    dlaci2_Iindt = rxn2_f - rxn2_b + rxn9_f - rxn9_b - rxn11                  #5 rxns
    dlacidt =   rxn3_b - rxn3_f - rxn6 + rxn7                                 #4 rxns
    dlaci2dt = rxn3_f - rxn3_b - rxn2_f + rxn2_b - rxn8 + rxn10_b - rxn10_f   #7 rxns
    dMlacidt = rxn4 - rxn5                                                    #2 rxns
    dR2Odt = rxn9_b - rxn9_f +rxn10_f - rxn10_b                               #4 rxns
    dOdt = -rxn9_b + rxn9_f - rxn10_f + rxn10_b                               #4 rxns
    dmRNAdt = rxn12 -rxn13 + rxn12_l                                          #3 rxns
    dprotdt = rxn14 - rxn15                                                   #2 rxns

    return [diptgexdt, diptgindt, dlacidt, dlaci2_Iindt, dlaci2dt, dMlacidt, dR2Odt, dOdt, dmRNAdt, dprotdt]



t = np.linspace(0., float(timef), int(times))
counter = 1
print(O0, iptgint)
for l in O0:
    for r in iptgint:
        x0 = [r, 0, 0, 0, 0, 0,0,l, 0,0]
        x = odeint(central_dogma, x0, t)

        iptgex = x[:,0]
        Iptgin = x[:,1]
        laci = x[:,2]
        laci2_Iin = x[:,3]
        laci2 = x[:,4]
        Mlaci = x[:,5]
        R2O = x[:,6]
        O = x[:,7]
        mRNA = x[:,8]
        prot = x[:,9]
        params = {
            'axes.labelsize':10,
            'font.size':15,
            'legend.fontsize':10,
            'xtick.labelsize':8,
            'ytick.labelsize':8,
            'figure.figsize': [8,8],
        }
        plt.rcParams.update(params)
        #plt.title("Variation of concentration with time")
        #plt.xlabel("Time (min)")
        #plt.ylabel("Concentration (nM)")

        #figure, axis = plt.subplots(3, 3)

        #axis[0, 0].plot(t, iptgex, label= 'iptg external')
        #axis[0, 0].set_title('iptg external')
        #axis[0, 0].grid()
        #axis[0, 0].legend()
        #axis[0, 1].plot(t, Iptgin, label= 'iptg internal')
        #axis[0, 1].set_title('iptg internal')
        #axis[0, 1].grid()
        #axis[0, 1].legend()
        #axis[0, 2].plot(t, laci, label= 'laci')
        #axis[0, 2].set_title('laci')
        #axis[0, 2].grid()
        #axis[0, 2].legend()
        #axis[1, 0].plot(t, laci2, label= 'laci2')
        #axis[1, 0].set_title('laci2')
        #axis[1, 0].grid()
        #axis[1, 0].legend()
        #axis[1, 1].plot(t, laci2_Iin, label= 'laci-iptg')
        #axis[1, 1].set_title('laci-iptg')
        #axis[1, 1].grid()
        #axis[1, 1].legend()
        #axis[1, 2].plot(t, Mlaci, label= 'LacI mRNA')
        #axis[1, 2].set_title('LacI mRNA')
        #axis[1, 2].grid()
        #axis[1, 2].legend()
        #axis[2, 0].plot(t, R2O, label= 'R2-O')
        #axis[2, 0].set_title('R2-O')
        #axis[2, 0].grid()
        #axis[2, 0].legend()
        #axis[2, 1].plot(t, O, label= 'lacO')
        #axis[2, 1].set_title('lacO')
        #axis[2, 1].grid()
        #axis[2, 1].legend()
        #axis[2, 2].plot(t, mRNA, label= ('target mRNA' + str(counter)))
        #axis[2, 2].set_title('target mRNA' + str(counter))
        #axis[2, 2].grid()
        #axis[2, 2].legend()
        ##plt.savefig("laco.svg")
        #plt.show()
        #counter += 1
        #o = 1
        #te = 'TE'
        #plt.plot(t, prot, label = 'protein')
        #plt.plot(t, mRNA, label = 'mRNA')
        #plt.ylabel('Protein concentration (nM)')
        #plt.xlabel('Time (min)')
        #plt.grid()
        #plt.legend()
        #plt.show()

        plt.plot(mRNA, prot, label = 'protein-mRNA')
        plt.ylabel('Protein concentration (nM)')
        plt.xlabel('mRNA concentration (nM)')
        plt.grid()
        plt.legend()
        plt.show()
