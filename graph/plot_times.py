import sys
import matplotlib.pyplot as plt
import pandas as pd

type_of_value = sys.argv[2]
# Read the CSV file
data = pd.read_csv(sys.argv[1])
# Sort the data by value in descending order
data.sort_values(by='seconds', ascending=False, inplace=True)
# Create the horizontal bar plot
plt.barh(data['rule'], data['seconds'])

# Set labels and title
plt.ylabel('Rule')
plt.xlabel(type_of_value)
plt.title('Difference between rules')
plt.gca().set_xticklabels(plt.gca().get_xticklabels(), rotation=45, ha='right')
# Display the plot
plt.show()
