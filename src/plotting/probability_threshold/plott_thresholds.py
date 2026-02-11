import numpy as np
import matplotlib.pyplot as plt

# Load ML scores
positive_scores = np.loadtxt('positive_ML_scores.txt')

# Thresholds to filter positive control
thresholds = [0.1, 0.2, 0.3, 0.4, 0.5]
filters = ['Unfiltered'] + [f'>{t}' for t in thresholds]

# Custom colors (add a color for unfiltered)
colors = ["#d62728", "#241fb4", "#1757cf", "#67aabd", "#4b768c", "#2ca08d"]

# Prepare filtered data for each threshold, include unfiltered first
filtered_data = [positive_scores] + [positive_scores[positive_scores > t] for t in thresholds]

plt.figure(figsize=(10, 6))

# Create boxplot
bp = plt.boxplot(filtered_data, patch_artist=True, labels=filters)

# Apply custom colors to boxes
for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.8)

# Style whiskers, caps, medians
for whisker in bp['whiskers']:
    whisker.set_color('black')
    whisker.set_linewidth(1.2)
for cap in bp['caps']:
    cap.set_color('black')
    cap.set_linewidth(1.2)
for median in bp['medians']:
    median.set_color('yellow')
    median.set_linewidth(2)

plt.xlabel('Minimum ML Filter')
plt.ylabel('ML Score (probability)')
plt.title('Positive Control ML Score Distributions')
plt.grid(axis='y', alpha=0.3)

plt.tight_layout()
plt.savefig('Positive_control_boxplots_with_unfiltered.png', dpi=300)
plt.show()
