
#-----------------------------------

def PID_RT(SP, PV, Man, MVMan, MVFF, KC, Ti, Td, alpha, Ts, MVMin, MVMax, MV, MVP, MVI, MVD, E, ManFF=False, PVInit=0, method='EBD'):
    """
    The function "PID_RT" needs to be included in a "for or while loop".
    SP: (or SetPoint) vector
    PV: (or Process Value) vector
    Man: (or Manual controller mode) vector (True or False)
    MVMan:(or Manual value for MV) vector
    MVFF: (or Feedforward) vector
    KC: controller gain
    Ti: integral time constant [s]
    Td: derivative time constant [s]
    alpha: TFD = alpha*Td where TFD is the derivative filter time constant [s]
    Ts: sampling period [s]
    MVMin: minimum value for MV (used for saturation and anti wind-up) :MVMax: maximum value for MV (used for saturation and anti wind-up)
    MV: MV (or Manipulated Value) vector
    MVP: MVP (or Propotional part of MV) vector
    MVI: MVI (or Integral part of MV) vector
    MVD: MVD (or Derivative part of MV) vector
    E: E (or control Error) vector
    ManFF: Activated FF in manual mode (optional: default boolean value is False)
    PVInit: Init value for PV (optional: default value is 0): used if PID_RT is ran first in the squence and no value of PV available yet.
    method: discretisation method (optional: default value is 'EBD')
    EBD-EBD: EBD for integral action and EBD for derivative action 
    EBD-TRAP: EBD for integral action and TRAP for derivative action 
    TRAP-EBD: TRAP for integral action and EBD for derivative action 
    TRAP-TRAP: TRAP for integral action and TRAP for derivative action
    The function "PID_RT" appends new values to the vectors "MV", "MVP", "MVI", and "MVD".
    The appended values are based on the PID algorithm, the controller mode, and feedforward.
    Note that saturation of "MV" within the limits [MVMin MVMax] is implemented with anti wind-up.
"""
    #Initialiwation
    TFD = Td * alpha
    #-----------------ERROR----------------
    if len(PV) ==0 :
        print("E = ", E, "  \len(E) = ", len(E), "Rest = ", SP[-1] - PVInit)
        E.append(SP[-1] - PVInit)
    else:
        E.append(SP[-1] - PV[-1])

    #----------------Proportionnal_ Action-----------
    if len(MVP)==0 :
        MVP.append((KC*E[-1]))
    else:
        if method == 'EBD':
            MVP.append(KC*E[-1])
        elif method == 'TRAP':  
            MVP.append(KC*E[-1])        
        else: 
            MVP.append(KC*E[-1])
    
    #-----------------Integral Action---------
    if len(MVI) == 0 :      
        MVI.append((KC*Ts/Ti)*E[-1])    # MV[k] is MV[-1] and MV[k-1] is MV[-2]
    else:
        if method == 'EBD':
            MVI.append(MVI[-1] + ((KC*Ts)/Ti)*(E[-1]))
        elif method == 'TRAP':
            MVI.append(MVI[-1] + (0.5*KC*Ts/Ti)*(E[-1] + E[-2]))
        else:      #default :EBD
            MVI.append(MVI[-1] + (KC*Ts/Ti)*E[-1])

    #------------------Derivative Action-------------
    if len(MVD)==0:
        MVD.append(((KC*Td)/(TFD+Ts))*(E[-1]))
    else:
        if method == 'EBD':
            MVD.append((TFD/(TFD+Ts))*MVD[-1] + ((KC*Td)/(TFD+Ts))*(E[-1] - E[-2])) 
        elif method =='TRAP':
            MVD.append((((TFD-(Ts/2))/(TFD+(Ts/2)))*MVD[-1] + ((KC*Td)/(TFD+(Ts/2)))*(E[-1] - E[-2])))     
        else:  # EBD
            MVD.append((TFD/(TFD+Ts))*MVD[-1] + ((KC*Td)/(TFD+Ts))*(E[-1] - E[-2]))



    #------------Mode manuel + anti-wind up (integrator help)----------
    if Man[-1] == True:
        if ManFF:
            MVI[-1] = MVMan[-1] - MVP[-1] - MVD[-1] 
        else:
            MVI[-1] = MVMan[-1] - MVP[-1] - MVD[-1] - MVFF[-1]

    # -----------------Saturation----------------------
    if (MVP[-1] + MVI[-1] + MVD[-1] + MVFF[-1]) > MVMax:
        MVI[-1] = MVMax - MVP[-1] - MVD[-1] - MVFF[-1]
    elif (MVP[-1] + MVI[-1] + MVD[-1] + MVFF[-1]) < MVMin:
        MVI[-1] = MVMin - MVP[-1] - MVD[-1] - MVFF[-1]

    # Ajout sur MV
    MV.append(MVP[-1] + MVI[-1] + MVD[-1] + MVFF[-1])



#----------------------------------------------------------------------------------------------   
def IMCTuning(Kp, Tlag1, Tlag2, theta, gamma, model="SOPDT"):
    """
    Computes the optimised PID controller settings for FOPDT/SOPDT processes.
    :Kp: process gain
    :Tlag1: processfirst lag time constant.
    :Tlag2: (SOPDT only) process second lag time constant.
    :theta: process delay.
    :gamma: loop response time as a ratio of T1
    :model:
        FOPDT_PI: First Order Plus Dead Time for P-I control (IMC tuning case G)
        FOPDT_PID: First Order Plus Dead Time for P-I-D control (IMC tuning case H)
        SOPDT :Second Order Plus Dead Time for P-I-D control (IMC tuning case I)
        
    return : PID controller parameters Kc, Ti and Td
    """
    T_CLP = gamma * Tlag1
    
    if model=="FOPDT_PI":
        Kc = (Tlag1/(T_CLP+theta))/Kp
        Ti = Tlag1
        Td = 0
    elif model=="FOPDT_PID":
        Kc= ((Tlag1 + theta/2)/(T_CLP + theta/2))/Kp
        Ti = Tlag1 + theta/2
        Td = (Tlag1*theta)/(2*Tlag1+theta)
    elif model=="SOPDT": 
        Kc = ((Tlag1 + Tlag2)/(T_CLP + theta))/Kp
        Ti = (Tlag1 +Tlag2)
        Td = ((Tlag1*Tlag2))/(Tlag1+Tlag2)
    
    return (Kc, Ti, Td)




#-----------------------------------
def Margin(Ps,C,omega,Show=True):
    """
    Calcule la marge de gain et la marge de phase. Elles nous permettent d'analyser la robustesse du PID.
    :Ps : Process
    :C: Fonction de transfert du Contrôleur 
    :omega : vecteur de la fréquence 
    :show : autorise l'affichage graphique
    """

    # Initialisation des paramètres
    s = 1j*omega
    Kc = C.parameters['Kc']
    Ti = C.parameters['Ti']
    Td = C.parameters['Td']
    Tfd = C.parameters['Tfd']

    # Calcul du Controller 
    Cs = Kc*(1 + 1/(Ti*s)+ (Td*s)/(Tfd*s +1))
    
    # Loop gain L(s) = P(s)C(s)
    Ls =Cs*Ps 

    # Plot de L(s)
    if Show ==True:
        fig, (ax_freq, ax_time) = plt.subplots(2,1)
        fig.set_figheight(12)
        fig.set_figwidth(22)

        # Amplitude
        ax_freq.semilogx(omega,20*np.log10(np.abs(Ls)),label='L(s)')
        gain_min = np.min(20*np.log10(np.abs(Ls)/5))
        gain_max = np.max(20*np.log10(np.abs(Ls)*5))
        ax_freq.set_xlim([np.min(omega), np.max(omega)])
        ax_freq.set_ylim([gain_min, gain_max])
        ax_freq.set_ylabel('Amplitude |P| [db]')
        ax_freq.set_title('Bode plot of P')
        ax_freq.legend(loc='best')
            
        # Phase
        ax_time.semilogx(omega, (180/np.pi)*np.unwrap(np.angle(Ls)),label='L(s)')
        ax_time.set_xlim([np.min(omega), np.max(omega)])
        ph_min = np.min((180/np.pi)*np.unwrap(np.angle(Ps))) - 10
        ph_max = np.max((180/np.pi)*np.unwrap(np.angle(Ps))) + 10
        ax_time.set_ylim([np.max([ph_min, -200]), ph_max])
        ax_time.set_ylabel(r'Phase $\angle P$ [°]')
        ax_time.legend(loc='best')
        ax_freq.axhline(y=0,color='green')
        ax_time.axhline(y=-180,color='cyan')

    # Crossover frequency
    i = 0
    for value in Ls:   # slide 69
        i+=1
        dB = 20*np.log10(np.abs(value))
        if dB < 0.05 and dB > -0.05:
            OmegaC =  omega[i-1]
            PhaseC = np.angle(value,deg=True)
            break        

    # Ultimate Frequency
    n = 0
    for value in Ls:
        n+=1
        deg = np.angle(value,deg=True)
        if deg < -179.5 and deg > -180.5:
            OmegaU = omega[n-1]
            u_freq = 20*np.log10(np.abs(value))
            break
    
    # Affichage graphique
    if Show ==True:
        ax_freq.plot([OmegaU,OmegaU],[0,u_freq]) 
        ax_time.plot([OmegaC,OmegaC],[PhaseC,-180])
    print('Gain margin :',-u_freq,'dB at the ultimate frequency :',OmegaU,'rad/s')
    print('Phase margin : ',PhaseC +180,'° at the crossover frequency :',OmegaC,'rad/s')
    
    # # Save picture
    # nameFile = 'Plots/Gain_phase_margin'

    # if not os.path.exists('Plots'):
    #     os.makedirs('Plots')

#-----------------------------------        
# class Controller:  // pour plus tard
    
    # def __init__(self, parameters):
        
    #     self.parameters = parameters
    #     self.parameters['Kp'] = parameters['Kp'] if 'Kp' in parameters else 1.0
    #     self.parameters['theta'] = parameters['theta'] if 'theta' in parameters else 0.0
    #     self.parameters['Tlead1'] = parameters['Tlead1'] if 'Tlead1' in parameters else 0.0
    #     self.parameters['Tlead2'] = parameters['Tlead2'] if 'Tlead2' in parameters else 0.0
    #     self.parameters['Tlag1'] = parameters['Tlag1'] if 'Tlag1' in parameters else 0.0
    #     self.parameters['Tlag2'] = parameters['Tlag2'] if 'Tlag2' in parameters else 0.0
    #     self.parameters['nInt'] = parameters['nInt'] if 'nInt' in parameters else 0        

#-----------------------------------        


