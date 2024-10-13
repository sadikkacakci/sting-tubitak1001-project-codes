import numpy as np
import pandas as pd


class DummyDataGenerator:

    def __init__(self,treatment_process,number_of_patients):
        self.treatment_process = treatment_process
        self.number_of_patients = number_of_patients
        self.features = {
            'Age': {'median': 4.75, 'range': (1.1, 17.1)},
            'Weight': {'median': 22, 'range': (10, 90)},
            'Height': {'median': 112.45, 'range': (80, 182.7)},
            'Body surface area': {'median': 0.82, 'range': (0.47, 1.98)},
            '6MP daily dose': {'median': 40, 'range': (5, 150)},
            'MTX weekly dose': {'median': 15, 'range': (1.25, 60)},
            'ANC': {'median': 1.8, 'range': (0.0, 19.9)}
            }
        self.dummy_data = None


    def calculate_bsa(self, weight, height):
        """
        Mosteller'ın vücut yüzey alanı (BSA) formülüne göre BSA hesaplar.
        
        Parametreler:
        weight (float): Kilogram cinsinden ağırlık
        height (float): Santimetre cinsinden boy
        
        Dönüş:
        float: Metrekare cinsinden vücut yüzey alanı (BSA)
        """
        return round(np.sqrt((weight * height) / 3600), 2)

    def get_dosages(self):
        def generate_positive_normal(median, range_min, range_max):
            std_dev = (range_max - range_min) / 6
            value = np.random.normal(median, std_dev)
            while value <= 0:
                value = np.random.normal(median, std_dev)
            return value
        total_day = 0
        daily_6mp_dose_list = []
        weekly_mtx_dose_list = []
        loop_control = True
        while loop_control:
            # daily_6mp_dose = np.random.normal(features['6MP daily dose']['median'], (features['6MP daily dose']['range'][1] - features['6MP daily dose']['range'][0]) / 6)
            # weekly_mtx_dose = np.random.normal(features['MTX weekly dose']['median'], (features['MTX weekly dose']['range'][1] - features['MTX weekly dose']['range'][0]) / 6)
            daily_6mp_dose = generate_positive_normal(self.features['6MP daily dose']['median'], self.features['6MP daily dose']['range'][0], self.features['6MP daily dose']['range'][1])
            weekly_mtx_dose = generate_positive_normal(self.features['MTX weekly dose']['median'], self.features['MTX weekly dose']['range'][0], self.features['MTX weekly dose']['range'][1])
            period = np.random.randint(4, 15)
            if total_day + 14 > self.treatment_process:
                period = self.treatment_process - total_day
                loop_control = False
            for _ in range(period):
                daily_6mp_dose_list.append(daily_6mp_dose)
                weekly_mtx_dose_list.append(weekly_mtx_dose)
            total_day += period
        return daily_6mp_dose_list, weekly_mtx_dose_list

    def create_data(self):
        # Zaman serisi uzunluğu (3 ay = 90 gün)

        # DataFrame için boş liste oluştur
        data_list = []

        # Her bir hasta için veri oluştur
        for patient_id in range(1000, 1000 + self.number_of_patients):  # 10 hasta için örnek oluşturuyoruz
            # Rastgele hasta özellikleri oluştur
            age = np.random.normal(self.features['Age']['median'], (self.features['Age']['range'][1] - self.features['Age']['range'][0]) / 6)
            weight = np.random.normal(self.features['Weight']['median'], (self.features['Weight']['range'][1] - self.features['Weight']['range'][0]) / 6)
            height = np.random.normal(self.features['Height']['median'], (self.features['Height']['range'][1] - self.features['Height']['range'][0]) / 6)
            bsa = self.calculate_bsa(weight, height)
            # daily_6mp_dose = np.random.normal(features['6MP daily dose']['median'], (features['6MP daily dose']['range'][1] - features['6MP daily dose']['range'][0]) / 6)
            # weekly_mtx_dose = np.random.normal(features['MTX weekly dose']['median'], (features['MTX weekly dose']['range'][1] - features['MTX weekly dose']['range'][0]) / 6)
            anc_initial = np.random.normal(self.features['ANC']['median'], (self.features['ANC']['range'][1] - self.features['ANC']['range'][0]) / 6)

            # ANC değerini simüle etmek için günlük rastgele sapma (örneğin ±0.2 G/L)
            anc_std_dev = 0.2
            anc_values = np.random.normal(anc_initial, anc_std_dev, self.treatment_process).clip(0, 19.9)
            
            # Zaman serisi verisini oluştur
            date_range = pd.date_range(start='2024-01-01', periods=self.treatment_process, freq='D')
            # daily_6mp = [daily_6mp_dose] * days
            # weekly_mtx = [weekly_mtx_dose] * days
            daily_6mp, weekly_mtx = self.get_dosages()
            weights = [weight] * self.treatment_process
            heights = [height] * self.treatment_process

            # Age fixing
            if age < 0 :
                age = 1
            # Dataframe'e ekle
            data = {
                'Patient_ID': [patient_id] * self.treatment_process,
                'Age' : int(age),
                'Date': date_range,
                'Weight_kg': weights,
                'Height_cm': heights,
                '6MP_Daily_Dose_mg': daily_6mp,
                'MTX_Weekly_Dose_mg': weekly_mtx,
                'ANC_G_L': anc_values,
                'BSA': bsa
            }
            data_list.append(pd.DataFrame(data))
        self.dummy_data = pd.concat(data_list, ignore_index=True)

    def get_dummy_data(self):
        self.create_data()
        return self.dummy_data