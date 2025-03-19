import numpy as np
from scipy.signal import lfilter

def lead_lag(MV, K, T_lead, T_lag, Ts, method='EBD'):
    """
    Implémente une fonction Lead-Lag.
    
    Arguments :
    - MV : Entrée du système (Manipulated Variable)
    - K : Gain statique
    - T_lead : Constante de temps du numérateur (lead)
    - T_lag : Constante de temps du dénominateur (lag)
    - Ts : Temps d'échantillonnage
    - method : 'EBD' (Euler Backward Difference) par défaut

    Retourne :
    - PV : Sortie du système après Lead-Lag
    """

    # Définition des coefficients de la fonction de transfert Lead-Lag
    num = [K * T_lead, K]  # Numérateur : K(T_lead*s + 1)
    den = [T_lag, 1]  # Dénominateur : (T_lag*s + 1)

    # Normalisation pour éviter des problèmes numériques
    num = np.array(num) / den[0]
    den = np.array(den) / den[0]

    # Simulation avec lfilter (équivalent à un filtre passe-bas Lead-Lag)
    PV = lfilter(num, den, MV)

    return PV


