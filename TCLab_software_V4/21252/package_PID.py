import numpy as np

class PIDController:
    def __init__(self, Kp, Ki, Kd, Ts, FF_gain=0):
        """
        Implémente un contrôleur PID discret avec Feedforward.

        Arguments :
        - Kp : Gain proportionnel
        - Ki : Gain intégral
        - Kd : Gain dérivé
        - Ts : Temps d'échantillonnage
        - FF_gain : Gain du feedforward (0 par défaut)

        """
        self.Kp = Kp
        self.Ki = Ki * Ts
        self.Kd = Kd / Ts
        self.FF_gain = FF_gain

        self.prev_error = 0
        self.integral = 0
        self.prev_DV = 0  # Pour la structure feedforward

    def compute(self, setpoint, measurement, disturbance=0):
        """
        Calcule la sortie du PID avec Feedforward.

        Arguments :
        - setpoint : Consigne (SP)
        - measurement : Mesure du processus (PV)
        - disturbance : Signal de perturbation pour le feedforward (DV)

        Retourne :
        - MV : Signal de commande (Manipulated Variable)
        """
        error = setpoint - measurement  # Erreur
        self.integral += error  # Terme intégral
        derivative = error - self.prev_error  # Terme dérivé

        # PID classique
        MV = self.Kp * error + self.Ki * self.integral + self.Kd * derivative

        # Ajout du Feedforward
        MV += self.FF_gain * (disturbance - self.prev_DV)

        # Mise à jour des valeurs précédentes
        self.prev_error = error
        self.prev_DV = disturbance

        return MV
