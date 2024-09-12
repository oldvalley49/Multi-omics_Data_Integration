import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

df = pd.read_csv("/Users/wunfs1010/Downloads/jilab_code/comparison/var/coembed/pancreas.csv", index_col=0)
print(df)

df = df.reset_index().rename(columns={'index': 'index'})

# Extract the algorithm and the number of features (var_num) from the row names
df['algorithm'] = df['index'].apply(lambda x: x.split('_v')[0])
df['var_num'] = df['index'].apply(lambda x: int(x.split('_v')[1]))

# Create the plot using sns.lineplot to plot FOSCTTM vs var_num for each algorithm
plt.figure(figsize=(10, 6))
sns.lineplot(data=df, x='var_num', y='foscttm', hue='algorithm', marker='o')

# Customize the plot
plt.title('Pancreas: Alignment vs Number of Features')
plt.xlabel('Number of Features (var_num)')
plt.ylabel('FOSCTTM Score')
plt.legend(title='Algorithm', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.grid(True)

# Show the plot
plt.tight_layout()
plt.show()
