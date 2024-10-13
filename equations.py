import numpy as np
import matplotlib.pyplot as plt
from dummy_data_generator import DummyDataGenerator
from scipy.integrate import solve_ivp



class ModelEquations:

    def __init__(self, patient_data):
        self.patient_data = patient_data
        self.length_6mp_dosages_list = len(self.patient_data["6MP_Daily_Dose_mg"].tolist())
        self.length_mtx_dosages_list = len(self.patient_data["MTX_Weekly_Dose_mg"].tolist())
        self.solution_6mp_model = None
        self.solution_mtx_model = None
        self.x3 = None
        self.x7 = None

    def get_6mp_model_parameters(self):
        k_a = 4.8 # emilme oranı, 6-MP absorption rate from GI tract
        alpha = (10**12) / 152177 # 10^^12/152177
        F = 0.45 #  Bioavailability factor
        # D6mp =  # 6mp ilaç alım dozajı, veriden gelecek.
        T_dur = 1/24 # ilaç emilimi için geçen süre (day)
        k_e = 5 #  6-MP elimination rate from plasma
        k_pt = 29.8 # 6-MP to 6-TGN conversion rate
        e_rel = 0.5 # TPMT enzyme activity constant
        K_t = 4.04 * (10**5) # M-M constant for 6-TGN
        k_pm = 655.8 # 6-MP to MeMP Conversion Rate
        K_m =  3.28 * 10**5# makalede değerini bulamadım. 
        V_pt = 1 # 6-TGN elimination rate from RBCs
        k_te = 0.0714 # 6-TGN elimination rate from RBCs
        return k_a,alpha,F,T_dur,k_e, k_pt, e_rel, K_t, k_pm, K_m, V_pt, k_te
    
    def get_D6mp(self, t):
        mp_dosages_list = self.patient_data["6MP_Daily_Dose_mg"].tolist()
        return mp_dosages_list[int(t)]
        # return 50

    def D6mp_model(self, t, y):
        k_a,alpha,F,T_dur,k_e, k_pt, e_rel, K_t, k_pm, K_m, V_pt, k_te = self.get_6mp_model_parameters()
        x1, x2, x3 = y
        D6mp = self.get_D6mp(t)
        dx1_dt = (-k_a * x1) + (alpha * F * D6mp / T_dur) # GI Bölgesindeki 6MP Yoğunluğu
        dx2_dt = (k_a * x1) - (k_e * x2) - ( (k_pt * (1 - e_rel) * x2) / (K_t + x2) ) - ( (k_pm * e_rel * x2) / (K_m + x2) ) # Plazma bölgesindeki 6MP Yoğunluğu
        dx3_dt = ( (V_pt * k_pt * (1 - e_rel) * x2) / (K_t + x2) ) - (k_te * x3) # Eritrositlerdeki 6-TGN Yoğunluğu
        return [dx1_dt, dx2_dt, dx3_dt]

    def apply_6mp_model(self):
        y0 = [0, 0, 0] # Başlangıç koşulları
        t_span = list(range(self.length_6mp_dosages_list))
        t_eval = np.linspace(t_span[0], t_span[1], 500)
        self.solution_6mp_model = solve_ivp(self.D6mp_model, t_span, y0, t_eval=t_eval)

    def get_mtx_model_parameters(self):
        k_a = 26.64 # GI'da emilim oranı.
        k_e = 5.76 # Plazmadan MTX emilim oranı.
        k_p = 9.6
        V = 11.606
        V_ml = 2.3895 * 10**4
        K_ml = 2.898
        K_eff = 179.76
        V_m_fpgs = 7.0119 * 10**3
        K_m_fpgs = 35.262
        K_ggh = 4.992
        F = 0.45 #  Bioavailability factor
        T_dur = 1/24
        # BSA = 152# Body surface area, veriden gelecek
        # D_mtx = _  # veriden gelecek
        beta = (10**6) / 454440
        return k_a, k_e, k_p, V, V_ml, K_ml, K_eff, V_m_fpgs, K_m_fpgs, K_ggh, F, T_dur , beta
    
    def get_mtx(self,t):
        mtx_dosages_list = self.patient_data["MTX_Weekly_Dose_mg"].tolist()
        return mtx_dosages_list[int(t)]
    
    def mtx_model(self,t,y, BSA):
        k_a, k_e, k_p, V, V_ml, K_ml, K_eff, V_m_fpgs, K_m_fpgs, K_ggh, F, T_dur , beta = self.get_mtx_model_parameters()
        x4,x5,x6,x7 = y
        D_mtx = self.get_mtx(t)
        dx4_dt = (-k_a * x4) + ( (beta * F * D_mtx) / (T_dur * BSA) )
        dx5_dt = (k_a * x4) - (k_e * x5)
        dx6_dt = ( ( (V_ml * x5) / V ) / (K_ml + ( x5 / V ) ) ) + ( (k_p * x5) / V ) - ( K_eff * x6 ) - ( (V_m_fpgs * x6) / (K_m_fpgs + x6)) + (K_ggh * x7)
        dx7_dt = ( (V_m_fpgs * x6) / (K_m_fpgs + x6)) - (K_ggh * x7)
        return [dx4_dt, dx5_dt, dx6_dt, dx7_dt]
    
    def apply_mtx_model(self):
        BSA = self.patient_data["BSA"].values[0]
        y0 = [0, 0, 0, 0] # Başlangıç koşulları
        t_span = list(range(self.length_mtx_dosages_list))
        t_eval = np.linspace(t_span[0], t_span[1], 500)
        self.solution_mtx_model = solve_ivp(self.mtx_model, t_span, y0, t_eval=t_eval, args= (BSA,))
    
    # def get_ALL_model_parameters(self):
    #     k_prol = 0###
    #     k_tr = 1###
    #     k_circ = 0.5346
    #     gamma = 1###
    #     Base = 1.0
    #     slope = 1 ###
    #     P_6MP = 1.0
    #     P_MTX = 1.0
    #     E_drug = np.log((slope * (P_6MP * self.x3 + P_MTX * self.x7) / 1000) + 1)
        
    #     return k_prol, k_tr, k_circ, gamma, Base, E_drug

    # def ALL_model(self,t,y):
    #     x8, x9, x10, x11, x12 = y
    #     k_prol, k_tr, k_circ, gamma, Base, E_drug = self.get_ALL_model_parameters()
    #     dx8_dt = (k_prol * x8) * (1 - E_drug) * ((Base / x12)**gamma) - (k_tr * x8)
    #     dx9_dt = k_tr * (x8 - x9)
    #     dx10_dt = k_tr * (x9 - x10)
    #     dx11_dt = k_tr * (x10 - x11)
    #     dx12_dt = (k_tr * x11) - (k_circ * x12)
    #     return [dx8_dt, dx9_dt, dx10_dt, dx11_dt, dx12_dt]
    
    # def apply_ALL_model(self):
    #     x8_9_10_11_ic = (Base * k_circ) / k_tr
    #     x12_ic = Base
    #     y0 = [x8_9_10_11_ic, x8_9_10_11_ic, x8_9_10_11_ic, x8_9_10_11_ic, x12_ic]
    #     t_span = list(range(len(mtx_dosages_list)))
    #     t_eval = np.linspace(t_span[0], t_span[1], 500)
    #     solution_ALL_model = solve_ivp(self.ALL_model, t_span, y0, t_eval=t_eval)


    def plot_6_mp_model_solutions(self):
        t = self.solution_6mp_model.t
        x1, x2, x3 = self.solution_6mp_model.y
        # Sonuçları çiz
        plt.figure()
        plt.plot(t, x1, label='x1')
        plt.plot(t, x2, label='x2')
        plt.plot(t, x3, label='x3')
        plt.xlabel('Zaman')
        plt.ylabel('Değer')
        plt.legend()
        plt.title('6MP Model Denklemleri Çözümü')
        plt.show()
    
    def plot_mtx_model_solutions(self):
        t = self.solution_mtx_model.t
        x4, x5, x6, x7 = self.solution_mtx_model.y
        plt.figure()
        plt.plot(t, x4, label='x4')
        plt.plot(t, x5, label='x5')
        plt.plot(t, x6, label='x6')
        plt.plot(t, x7, label='x7')
        plt.xlabel('Zaman')
        plt.ylabel('Değer')
        plt.legend()
        plt.title('MTX Model Denklemleri Çözümü')
        plt.show()