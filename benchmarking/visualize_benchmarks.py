#!/usr/bin/env python3
"""
Visualization script for Julia vs SageMath benchmark results
Generates publication-quality plots
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['font.size'] = 11

# Read data
df = pd.read_csv('benchmark_results.csv')

# ===================================
# 1. TIME COMPARISON BAR CHART
# ===================================
fig, ax = plt.subplots(figsize=(14, 6))

x = np.arange(len(df))
width = 0.35

bars1 = ax.bar(x - width/2, df['Julia_Mean_ms'], width, 
               label='Julia (Fraction-Free)', 
               color='#1f77b4', alpha=0.8,
               yerr=df['Julia_Std_ms'], capsize=3)

bars2 = ax.bar(x + width/2, df['Sage_Mean_ms'], width,
               label='SageMath',
               color='#ff7f0e', alpha=0.8,
               yerr=df['Sage_Std_ms'], capsize=3)

ax.set_xlabel('Test Case', fontweight='bold')
ax.set_ylabel('Time (milliseconds)', fontweight='bold')
ax.set_title('Julia Fraction-Free SNF vs SageMath - Execution Time Comparison', 
             fontweight='bold', fontsize=14)
ax.set_xticks(x)
ax.set_xticklabels(df['Test'].str.replace(' ', '\n'), rotation=45, ha='right', fontsize=9)
ax.legend()
ax.set_yscale('log')
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('benchmark_time_comparison.png', dpi=300, bbox_inches='tight')
print("✅ Saved: benchmark_time_comparison.png")

# ===================================
# 2. SPEEDUP FACTOR CHART
# ===================================
fig, ax = plt.subplots(figsize=(14, 6))

colors = ['#2ecc71' if s > 1 else '#e74c3c' for s in df['Speedup']]
bars = ax.bar(range(len(df)), df['Speedup'], color=colors, alpha=0.8, edgecolor='black')

# Add value labels on bars
for i, (bar, val) in enumerate(zip(bars, df['Speedup'])):
    height = bar.get_height()
    ax.text(bar.get_x() + bar.get_width()/2., height,
            f'{val:.1f}×',
            ha='center', va='bottom', fontweight='bold', fontsize=9)

ax.axhline(y=1, color='red', linestyle='--', linewidth=2, label='Equal Performance')
ax.set_xlabel('Test Case', fontweight='bold')
ax.set_ylabel('Speedup Factor (Julia / SageMath)', fontweight='bold')
ax.set_title('Julia Performance Advantage Over SageMath', fontweight='bold', fontsize=14)
ax.set_xticks(range(len(df)))
ax.set_xticklabels(df['Test'].str.replace(' ', '\n'), rotation=45, ha='right', fontsize=9)
ax.legend()
ax.grid(True, alpha=0.3, axis='y')

plt.tight_layout()
plt.savefig('benchmark_speedup_factors.png', dpi=300, bbox_inches='tight')
print("✅ Saved: benchmark_speedup_factors.png")

# ===================================
# 3. SCALABILITY ANALYSIS (Dense matrices only)
# ===================================
dense_df = df[df['Type'] == 'Dense Random'].copy()
if len(dense_df) > 0:
    dense_df['Size'] = dense_df['Rows']  # Assuming square matrices
    dense_df = dense_df.sort_values('Size')
    
    fig, ax = plt.subplots(figsize=(10, 7))
    
    ax.errorbar(dense_df['Size'], dense_df['Julia_Mean_ms'], 
                yerr=dense_df['Julia_Std_ms'],
                marker='o', markersize=8, linewidth=2, capsize=5,
                label='Julia (Fraction-Free)', color='#1f77b4')
    
    ax.errorbar(dense_df['Size'], dense_df['Sage_Mean_ms'],
                yerr=dense_df['Sage_Std_ms'],
                marker='s', markersize=8, linewidth=2, capsize=5,
                label='SageMath', color='#ff7f0e')
    
    ax.set_xlabel('Matrix Size (n×n)', fontweight='bold', fontsize=12)
    ax.set_ylabel('Time (milliseconds)', fontweight='bold', fontsize=12)
    ax.set_title('Scalability: Time vs Matrix Size (Random Dense Matrices)', 
                 fontweight='bold', fontsize=14)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.set_yscale('log')
    ax.set_xscale('log')
    
    plt.tight_layout()
    plt.savefig('benchmark_scalability.png', dpi=300, bbox_inches='tight')
    print("✅ Saved: benchmark_scalability.png")

# ===================================
# 4. MATRIX TYPE COMPARISON
# ===================================
type_stats = df.groupby('Type').agg({
    'Julia_Mean_ms': 'mean',
    'Sage_Mean_ms': 'mean',
    'Speedup': 'mean'
}).reset_index()

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

# Time comparison by type
x_types = np.arange(len(type_stats))
width = 0.35

ax1.bar(x_types - width/2, type_stats['Julia_Mean_ms'], width,
        label='Julia', color='#1f77b4', alpha=0.8)
ax1.bar(x_types + width/2, type_stats['Sage_Mean_ms'], width,
        label='SageMath', color='#ff7f0e', alpha=0.8)

ax1.set_xlabel('Matrix Type', fontweight='bold')
ax1.set_ylabel('Average Time (milliseconds)', fontweight='bold')
ax1.set_title('Average Performance by Matrix Type', fontweight='bold', fontsize=13)
ax1.set_xticks(x_types)
ax1.set_xticklabels(type_stats['Type'], rotation=45, ha='right')
ax1.legend()
ax1.grid(True, alpha=0.3, axis='y')

# Speedup by type
bars = ax2.bar(x_types, type_stats['Speedup'], color='#2ecc71', alpha=0.8, edgecolor='black')

for bar, val in zip(bars, type_stats['Speedup']):
    height = bar.get_height()
    ax2.text(bar.get_x() + bar.get_width()/2., height,
            f'{val:.1f}×',
            ha='center', va='bottom', fontweight='bold')

ax2.axhline(y=1, color='red', linestyle='--', linewidth=2)
ax2.set_xlabel('Matrix Type', fontweight='bold')
ax2.set_ylabel('Average Speedup Factor', fontweight='bold')
ax2.set_title('Average Speedup by Matrix Type', fontweight='bold', fontsize=13)
ax2.set_xticks(x_types)
ax2.set_xticklabels(type_stats['Type'], rotation=45, ha='right')
ax2.grid(True, alpha=0.3, axis='y')

plt.tight_layout()
plt.savefig('benchmark_by_type.png', dpi=300, bbox_inches='tight')
print("✅ Saved: benchmark_by_type.png")

# ===================================
# 5. BOX PLOT - Statistical Distribution
# ===================================
fig, ax = plt.subplots(figsize=(10, 6))

# Prepare data for box plot (simulating distribution from mean and std)
julia_times = []
sage_times = []
labels = []

for _, row in df.iterrows():
    # Approximate distribution
    julia_times.append(row['Julia_Mean_ms'])
    sage_times.append(row['Sage_Mean_ms'])
    labels.append(row['Test'][:20])  # Truncate long names

# Create grouped box plot
data_to_plot = [julia_times, sage_times]
bp = ax.boxplot(data_to_plot, labels=['Julia', 'SageMath'],
                patch_artist=True, widths=0.6)

bp['boxes'][0].set_facecolor('#1f77b4')
bp['boxes'][1].set_facecolor('#ff7f0e')

ax.set_ylabel('Time (milliseconds)', fontweight='bold')
ax.set_title('Statistical Distribution of Execution Times', fontweight='bold', fontsize=14)
ax.set_yscale('log')
ax.grid(True, alpha=0.3, axis='y')

plt.tight_layout()
plt.savefig('benchmark_distribution.png', dpi=300, bbox_inches='tight')
print("✅ Saved: benchmark_distribution.png")

# ===================================
# 6. SUMMARY STATISTICS TABLE FIGURE
# ===================================
fig, ax = plt.subplots(figsize=(12, 6))
ax.axis('tight')
ax.axis('off')

summary_data = [
    ['Metric', 'Value'],
    ['', ''],
    ['Total Tests', str(len(df))],
    ['All Results Match', '✅ Yes' if df['Match'].all() else '❌ No'],
    ['Average Speedup', f"{df['Speedup'].mean():.2f}×"],
    ['Median Speedup', f"{df['Speedup'].median():.2f}×"],
    ['Min Speedup', f"{df['Speedup'].min():.2f}×"],
    ['Max Speedup', f"{df['Speedup'].max():.2f}×"],
    ['', ''],
    ['Julia Avg Time', f"{df['Julia_Mean_ms'].mean():.2f} ms"],
    ['SageMath Avg Time', f"{df['Sage_Mean_ms'].mean():.2f} ms"],
]

table = ax.table(cellText=summary_data, cellLoc='left', loc='center',
                colWidths=[0.5, 0.5])
table.auto_set_font_size(False)
table.set_fontsize(12)
table.scale(1, 2.5)

# Style header
for i in range(2):
    table[(0, i)].set_facecolor('#4CAF50')
    table[(0, i)].set_text_props(weight='bold', color='white')

# Style rows
for i in range(2, len(summary_data)):
    for j in range(2):
        if i % 2 == 0:
            table[(i, j)].set_facecolor('#f0f0f0')

plt.title('Benchmark Summary Statistics', fontweight='bold', fontsize=16, pad=20)
plt.savefig('benchmark_summary_table.png', dpi=300, bbox_inches='tight')
print("✅ Saved: benchmark_summary_table.png")

print("\n" + "="*60)
print("All visualizations generated successfully!")
print("="*60)
print("\nGenerated files:")
print("  1. benchmark_time_comparison.png - Time comparison bar chart")
print("  2. benchmark_speedup_factors.png - Speedup factors")
print("  3. benchmark_scalability.png - Scalability analysis")
print("  4. benchmark_by_type.png - Performance by matrix type")
print("  5. benchmark_distribution.png - Statistical distribution")
print("  6. benchmark_summary_table.png - Summary statistics")
