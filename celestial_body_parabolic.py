
        # Parabolic Anomaly
        if self.type == "parabolic":
            self.parabolic_anomaly = np.sqrt(self.parameter) * np.tan(self.true_anomaly / 2.0) # D
