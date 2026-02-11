import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

# -----------------------------
# Load ML scores
# -----------------------------
positive_scores = np.loadtxt('positive_ML_scores.txt')
negative_scores = np.loadtxt('negative_ML_scores.txt')

# -----------------------------
# Thresholds to filter positive control
# -----------------------------
thresholds = [0.1, 0.2, 0.3, 0.4, 0.5]

for thresh in thresholds:
    # Filter positive scores
    filtered_pos = positive_scores[positive_scores > thresh]
    
    # Fit normal distributions
    mu_pos, std_pos = filtered_pos.mean(), filtered_pos.std()
    mu_neg, std_neg = negative_scores.mean(), negative_scores.std()
    
    # -----------------------------
    # Create figure
    # -----------------------------
    plt.figure(figsize=(10, 6))
    
    # Positive control histogram
    plt.hist(filtered_pos, bins=50, density=True, alpha=0.6, color='blue',
             edgecolor='black', label=f'Positive (BrdU+), n={len(filtered_pos):,}')
    
    # Negative control histogram (unfiltered)
    plt.hist(negative_scores, bins=50, density=True, alpha=0.6, color='orange',
             edgecolor='black', label=f'Negative (BrdU-), n={len(negative_scores):,}')
    
    # Normal distribution curves
    x_pos = np.linspace(filtered_pos.min(), filtered_pos.max(), 200)
    plt.plot(x_pos, norm.pdf(x_pos, mu_pos, std_pos), 'b-', linewidth=2,
             label=f'Positive Normal fit (μ={mu_pos:.3f}, σ={std_pos:.3f})')
    
    x_neg = np.linspace(0, 1, 200)
    plt.plot(x_neg, norm.pdf(x_neg, mu_neg, std_neg), 'r-', linewidth=2,
             label=f'Negative Normal fit (μ={mu_neg:.3f}, σ={std_neg:.3f})')
    
    # Labels, title, legend
    plt.xlabel('ML Score (probability)')
    plt.ylabel('Density')
    plt.title(f'ML Score Distributions: Positive vs Negative Control\nFiltered Positive > {thresh}')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Save figure
    plt.tight_layout()
    filename = f'ML_distributions_filtered_positive_{thresh:.1f}.png'
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.show()
    
    # -----------------------------
    # Print statistics
    # -----------------------------
    print("="*60)
    print(f"POSITIVE CONTROL (BrdU+) Filtered > {thresh}")
    print("="*60)
    print(f"Mean ML score: {mu_pos:.3f}")
    print(f"Std deviation: {std_pos:.3f}")
    print(f"Median: {np.median(filtered_pos):.3f}")
    print(f"% above 0.3: {100*np.mean(filtered_pos >= 0.3):.1f}%")
    print(f"% above 0.4: {100*np.mean(filtered_pos >= 0.4):.1f}%")
    print(f"% above 0.5: {100*np.mean(filtered_pos >= 0.5):.1f}%")
    
    print("\nNEGATIVE CONTROL (BrdU-)")
    print("="*60)
    print(f"Mean ML score: {mu_neg:.3f}")
    print(f"Std deviation: {std_neg:.3f}")
    print(f"Median: {np.median(negative_scores):.3f}")
    print(f"% above 0.3: {100*np.mean(negative_scores >= 0.3):.1f}%")
    print(f"% above 0.4: {100*np.mean(negative_scores >= 0.4):.1f}%")
    print(f"% above 0.5: {100*np.mean(negative_scores >= 0.5):.1f}%")
