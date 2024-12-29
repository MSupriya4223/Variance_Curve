import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from sklearn.preprocessing import MinMaxScaler

# Mention the Dataset name
x = 'GSE18842.xlsx'

# Load the Dataset using pandas
df = pd.read_excel(x)
# print(df.head(5))

# Define the keys to specify 'Healthy' and 'Disease' samples
Healthy = 'Control'
Disease = 'Tumor'

# Preprocessing
df = df.drop(['!Sample_title'], axis=1)
df = df.dropna(subset=['Gene Symbol'])

# Handle the Duplicate rows based on Gene name
duplicates = df[df.duplicated(subset='Gene Symbol', keep=False)]
avg_duplicates = duplicates.groupby('Gene Symbol').mean().reset_index()
filtered_df = df[~df.duplicated(subset='Gene Symbol', keep=False)]
df = pd.concat([filtered_df, avg_duplicates])

# Normalize the matrix
scaler = MinMaxScaler()
normalized = scaler.fit_transform(df.iloc[:, 1:])
normalized_df = pd.DataFrame(normalized, columns=df.columns[1:]).reset_index(drop=True)
normalized_df['Gene Symbol'] = df['Gene Symbol'].reset_index(drop=True)
normalized_df.set_index('Gene Symbol', inplace=True)

# Transpose the matrix
transposed_df = normalized_df.transpose()

# Calculate variance
variance_df = transposed_df.var()

# Sorting top variance
top_variance = variance_df.nlargest(40)

# print("Top Genes are: ", top_variance.index.tolist())

# Set threshold
threshold = 0.06

# Prepare the plot
plt.figure(figsize=(12, 8))

# Plot all points and segments
for i in range(len(top_variance) - 1):
    y1 = top_variance.iloc[i]
    y2 = top_variance.iloc[i + 1]

    if y1 > threshold and y2 > threshold:
        plt.plot([i, i + 1], [y1, y2], 'b-')  # Blue line segment above threshold
        plt.plot(i, y1, 'bo', markersize=8)  # Blue point above threshold
    elif y1 <= threshold and y2 <= threshold:
        plt.plot([i, i + 1], [y1, y2], 'gray')  # Gray line segment below threshold
        plt.plot(i, y1, 'o', color='gray', markersize=8)  # Gray point below threshold
    else:
        # Intersection case
        plt.plot([i, i + 0.5], [y1, threshold], 'b-')  # Blue part to threshold
        plt.plot([i + 0.5, i + 1], [threshold, y2], 'gray')  # Gray part from threshold
        plt.plot(i, y1, 'bo', markersize=8)  # Point on original side

# Plot the last point
if top_variance.iloc[-1] > threshold:
    plt.plot(len(top_variance) - 1, top_variance.iloc[-1], 'bo', markersize=8)  # Blue point
else:
    plt.plot(len(top_variance) - 1, top_variance.iloc[-1], 'o', color='gray', markersize=8)  # Gray point

# Add horizontal threshold line
plt.axhline(y=threshold, color='r', linestyle='-', label='Threshold')

# Create custom legend
blue_line = mlines.Line2D([], [], color='blue', marker='o', linestyle='-', markersize=10, label='Above Threshold')
gray_line = mlines.Line2D([], [], color='gray', marker='o', linestyle='-', markersize=10, label='Below Threshold')
threshold_line = mlines.Line2D([],[], color='r', linestyle='-', label='Threshold Line')
                               
# Set x-axis labels and title
plt.xticks(ticks=range(40), labels=top_variance.index, rotation=90, fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel('Genes', fontsize=24)
plt.ylabel('Variance', fontsize=24)
plt.title('Line Graph of the Top 40 Variance')

# Add grid and legend
plt.grid(axis='both', linestyle='--', alpha=0.5)
plt.legend(fontsize=22)
plt.legend(handles=[blue_line, gray_line, threshold_line], fontsize=22)
plt.tight_layout()

# # Save the figure as a PNG file
# plt.savefig('GSE18842__Variance_Line_Diagram 02.png', dpi=500)  # Adjust dpi for higher resolution if needed

# Show plot
plt.show()

# Print the genes above threshold
variance_significant_genes = top_variance[top_variance > threshold].index.tolist()

print("Genes above threshold are: ", variance_significant_genes)
print(f"\nThere are {len(variance_significant_genes)} Genes.")
