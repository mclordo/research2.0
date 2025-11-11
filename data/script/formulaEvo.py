from datetime import datetime
import math
from scipy.interpolate import UnivariateSpline
import numpy as np

# Daten von https://kraege.de/wp-content/uploads/2024/11/www_Katalog_25_dt.pdf Seite 33 Frigo
KWIdeal = np.array([18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36])
WaterNeedIdeal = np.array([1,1.8,2.3,2.2,1.9,1.7,1.7,1.75,1.9,2,1.95,1.85,1.85,2,2.15,2.1,1.9,1.8,1.75])
KW = np.array([18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39])
WaterNeed = np.array([0.5,0.6,0.7,0.8,2.6,1,0.9,1,1.1,1.7,1.95,1.3,1.5,1.8,2.15,2.1,1.9,1.6,1.5,1.4,1.3,1.2])

# Werte für die Ausgleichskurve berechnen
KcSplineFit = UnivariateSpline(KW, WaterNeed, s=0, k=2)
KW_fit = np.linspace(KW.min(), KW.max(), 200)
Kc_fit = KcSplineFit(KW_fit)

KcSplineFitIdeal = UnivariateSpline(KWIdeal, WaterNeedIdeal, s=0)
KW_fitIdeal = np.linspace(KWIdeal.min(), KWIdeal.max(), 200)
Kc_fit_ideal = KcSplineFitIdeal(KW_fitIdeal)

def krKorean(T_avg):
    return (0.00185 * (T_avg)**2) - (0.0433 * (T_avg)) + 0.4023

def krSpanish(T_avg):
    return 0.0008 * (T_avg)**2 - 0.0279 * (T_avg) + 0.4017

def rso_fao56(dayOfYear, lat_deg, elev_m, kwh=False):
    """
    Berechnet die Clear-Sky-Kurzwelleneinstrahlung R_so nach FAO-56.

    Parameter
    ---------
    dayOfYear : int
        DOY (1..366)
    lat_deg : float
        Geografischer Breitengrad in Grad (+N, -S).
    elev_m : float
        Höhe über Meeresspiegel in Metern.
    kwh : bool, optional
        Wenn True, Rückgabe in kWh/m²/Tag statt MJ/m²/Tag.

    Rückgabe
    --------
    float
        R_so als Tageswert (MJ/m²/Tag oder kWh/m²/Tag).
    """
    # 2) Konstanten und Winkel
    G_SC = 0.0820  # MJ m^-2 min^-1
    phi = math.radians(lat_deg)
    pi = math.pi

    # 3) Erdbahn und Geometrie (FAO-56)
    dr = 1.0 + 0.033 * math.cos(2.0 * pi * dayOfYear / 365.0)                          # relative Entfernung Erde-Sonne
    delta = 0.409 * math.sin(2.0 * pi * dayOfYear / 365.0 - 1.39)                      # Solardeklination
    cos_omega_s = -math.tan(phi) * math.tan(delta)                                   # Sonnenuntergangswinkel
    cos_omega_s = max(-1.0, min(1.0, cos_omega_s))                                   # numerische Klammer
    omega_s = math.acos(cos_omega_s)

    # 4) Extraterrestrische Strahlung Ra (MJ/m²/Tag)
    Ra = (24.0 * 60.0 / pi) * G_SC * dr * (
        omega_s * math.sin(phi) * math.sin(delta)
        + math.cos(phi) * math.cos(delta) * math.sin(omega_s)
    )

    # 5) Clear-Sky-Strahlung Rso (MJ/m²/Tag)
    Rso = (0.75 + 2.0e-5 * elev_m) * Ra

    # Optional: Ausgabe in kWh/m²/Tag
    if kwh:
        return Rso / 3.6
    return Rso

def calculate_evaporation_solar(date,Tdiff_daily,Kr,radiation):
    """
    Berechnet die Verdunstung basierend auf Solarstrahlung und Temperatur.
    Parameter:
    date : datetime
        Datum für die Berechnung.
    Tdiff_daily : float
        Tägliche Temperaturdifferenz in °C.
    radiation : float
        Solarstrahlung in MJ/m²/Tag.
    Kr : float
        Korrekturfaktor basierend auf der durchschnittlichen Tagestemperatur.

    Rückgabe:
    float array
        Verdunstung in mm/m²*Tag.
    """

    # Anpassung des Wasserverbrauches 4m² gemessen auf 1m² berechnung
    AreaCompensationFactor = 4  

    # constant parameters
    tau = 0.75
    hightOverSea = 596  # in m, can be set to the height of the location    
        
    # datum = datetime(2025, 7, 24)  # Beispiel: 24. Juli 2025
    DOY = date.timetuple().tm_yday

    KW_fraction = date.isocalendar().week + (date.isocalendar().weekday - 1) / 7 
    Kc = KcSplineFit(KW_fraction)  # Example value, can be adjusted based on the crop stage

    # Berechnung Rs
    Rs = Kr * (Tdiff_daily)** 0.5 * radiation
    RsmaxApprox = (0.75 + 2*10**-5 * hightOverSea) * radiation
    Rsmax = rso_fao56(DOY, 47.47337, hightOverSea, kwh=False)  # Example for Berlin

    # Berrechung Rs_green
    Rs_green = Rs * tau

    # ETo_green Annäherung
    if DOY <= 220:
        ETo_green = (0.288 + 0.0019* DOY) * Rs_green 
    elif DOY > 220:
        ETo_green = (1.3398 - 0.00288* DOY) * Rs_green 

    # Wasserverbrauch
    calculate_evaporation_solar = ETo_green * Kc * AreaCompensationFactor
    #correction factor for measuring 2m² ad calculating 1m²

    return calculate_evaporation_solar