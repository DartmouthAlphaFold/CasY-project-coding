import numpy as np
# import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from sklearn.linear_model import LinearRegression

# Generate some sample data
np.random.seed(0)
num_samples = 1000
num_features = 5
X = np.random.rand(num_samples, num_features)
true_coeffs = np.array([2, -1, 0.5, 10, 7])
noise = np.random.normal(0, 1, num_samples)
y = np.dot(X, true_coeffs) + noise

# Define the linear function
def linear_function(X, *coeffs):
    return np.dot(X, coeffs)

# Define the weights
weights = np.random.uniform(0.5, 2, num_samples)  # Example weights

# Perform weighted linear regression using curve_fit
initial_guess = np.zeros(num_features)
popt, pcov = curve_fit(linear_function, X, y, p0=initial_guess, sigma=1/np.sqrt(weights))

# Get the estimated coefficients
coeffs = popt

# Print the results
print("Estimated Coefficients:", coeffs)



# Perform weighted linear regression using sklearn
model = LinearRegression(fit_intercept=True)
model.fit(X, y, sample_weight=weights)

# Get the estimated parameters
coeffs = model.coef_

print("Estimated Coefficients:", coeffs)