import numpy as np
import scipy.signal as signal

class PID_IMC:
    def __init__(self, Kp, Ti, Td, Ts):
        """
        Implémente un contrôleur PID basé sur la méthode IMC.
        
        Arguments :
        - Kp : Gain proportionnel
        - Ti : Temps intégral (Ti = Kp / Ki)
        - Td : Temps dérivé
        - Ts : Temps d'échantillonnage
        """
        self.Kp = Kp
        self.Ti = Ti
        self.Td = Td
        self.Ts = Ts
        
        self.prev_error = 0
        self.integral = 0

    def compute(self, setpoint, measurement):
        """
        Calcule la sortie du PID optimisé par IMC.
        
        Arguments :
        - setpoint : Consigne (SP)
        - measurement : Valeur mesurée (PV)
        
        Retourne :
        - MV : Commande calculée
        """
        error = setpoint - measurement
        self.integral += error * self.Ts
        derivative = (error - self.prev_error) / self.Ts

        MV = self.Kp * (error + (1 / self.Ti) * self.integral + self.Td * derivative)

        self.prev_error = error

        return MV

def IMC_PID(Tg, Tu, Kp, lambda_tuning):
    """
    Calcule les paramètres PID optimaux avec la méthode IMC.
    
    Arguments :
    - Tg : Constante de temps du processus
    - Tu : Temps mort du processus
    - Kp : Gain du processus
    - lambda_tuning : Paramètre d'ajustement IMC (plus grand = plus stable, plus petit = plus réactif)
    
    Retourne :
    - (Kp_PID, Ti, Td) : Paramètres optimaux du PID
    """
    Kp_PID = (Tg / (Kp * (lambda_tuning + Tu)))
    Ti = Tg
    Td = Tu / 2

    return Kp_PID, Ti, Td

def compute_stability_margins(num, den):
    """
    Calcule les marges de stabilité (Gain Margin, Phase Margin).
    
    Arguments :
    - num : Numérateur de la fonction de transfert en boucle ouverte
    - den : Dénominateur de la fonction de transfert en boucle ouverte
    
    Retourne :
    - Gain Margin (GM), Phase Margin (PM)
    """
    system = signal.TransferFunction(num, den)
    w, mag, phase = signal.bode(system)

    GM = 1 / (10**(np.min(mag) / 20))  # Gain Margin
    PM = 180 + np.min(phase)  # Phase Margin
    
    return GM, PM
