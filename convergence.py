import pandas as pd
import matplotlib.pyplot as plt

# Load the convergence data from the CSV file
data = pd.read_csv('convergence 5.csv', header=None, names=['Iteration', 'Function Value'])

# Plotting the convergence without log scale
plt.figure(figsize=(10, 6))
plt.plot(data['Iteration'], data['Function Value'], marker='o', linestyle='-', color='b')
plt.title('Convergence Plot for Trid Function: Function Value vs. Iteration Number')
plt.xlabel('Iteration Number')
plt.ylabel('Function Value')
plt.grid()
plt.show()