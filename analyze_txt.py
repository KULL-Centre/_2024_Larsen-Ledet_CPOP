import json
import pandas as pd
import argparse
from functions import *

parser = argparse.ArgumentParser(description='Reading the sequence data json')
parser.add_argument('-m', '--main_path', type=str,
                    help='Path to the folder with all txt files (see README)')
parser.add_argument('-ij', '--input_json', type=str,
                    help='Path to the json with the data required for score calculation (see README)')
parser.add_argument('-t', '--threshold', type=int,
                    help='Read count threshold for the control condition, default = 0', default = 0)
parser.add_argument('-p', '--printing', type=bool,
                    default=True,
                    help='True or False for printing additional stats while processing')
parser.add_argument('-n', '--norm_type',
                    choices=['syn_sum', 'wt'],
                    default='sum_syn',
                    help='Normalisation type: to wild type counts (wt), or to sum of synonimous variants at all positions in the tile (sum_syn, default)')

args = parser.parse_args()
json_path = args.input_json
main_path = args.main_path
threshold = args.threshold
printing = args.printing
norm_type = args.norm_type

# Open and read the JSON file
with open(json_path, 'r') as file:
    sequence_data_dict = json.load(file)
    
tile_shared_averaged, df_full, tile_list, wt_df, raw_var_call, filtered_out, pre_rescale_shared_averaged, tile_mod_dict = analysis(main_path, 
    sequence_data_dict, normalisation_scaling_type = norm_type, print_yn = printing, threshold = threshold)

df_full.to_csv('extra_data_df.csv', index=False)
tile_shared_averaged.to_csv('all_tiles_data.csv', index=False)

# Save count files per tile for Sven
list_counts = list()

for tile in raw_var_call.keys():
    df_list = list()
    for cond in raw_var_call[tile].keys():
        df = raw_var_call[tile][cond].copy()
        df = df.loc[df['mut_type']!='wt']
        modifier = tile_mod_dict[tile]
        df['pos'] = df['pos'] + modifier
        count_cols = [c for c in df.columns if c.startswith('count')]
        new_cols = [f'{c}_{cond}' for c in count_cols]
        df[new_cols] = df[count_cols]
        df.drop(columns=count_cols, inplace=True)
        df_list.append(df)
    df_tile = functools.reduce(lambda x, y: x.merge(y, on=['mut_type', 'pos', 'sequence'], how='left'), df_list)
    df_tile.to_csv(f'per_tile_variant_counts_{tile}.csv', index=False)
    
    
#Heatmaps for all tiles all conditions
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colorbar import ColorbarBase

# Designed mutation classes
mut_list = [mut for mut in sequence_data_dict['values'] if mut != 'wt']
for mut in mut_list:
    for d in ['score', 'err']:
        df = tile_shared_averaged.loc[tile_shared_averaged['mut_type']==mut].copy()
        df = df[['pos', 'aa']+[col for col in df.columns if col.startswith(d)]]

        cmax = df[[col for col in df.columns if col.startswith(d)]].max().max()
        cmin = df[[col for col in df.columns if col.startswith(d)]].min().min()

        # Create a custom colormap
        if d == 'score':
            colors = [(0.66, 0.33, 0.33), (1, 1, 1)]  # red to white
        else:
            colors = [(1, 1, 1), (1, 0.5, 0)]  # white to orange
        
        cmap_name = 'custom_cmap'
        cm = LinearSegmentedColormap.from_list(cmap_name, colors)

        plt.figure(figsize=(80, 8))  # Set the figure size

        plt.title(f'{mut} {d}', fontsize=18)  # Set the title of the entire figure

        #colorbar_ax = plt.gcf().add_axes([0.95, 0.15, 0.02, 0.7])  # Adjust the position of the colorbar


        #pivot_df = df.pivot_table(values=condition, index='mut_type', columns='pos')
        pivot_df = df.set_index('pos').drop(columns=['aa']).T

        # Reorder mut_type according to the specified order
        # pivot_df = pivot_df.reindex(mut_type_order)

        # Create a mask for NaN values
        mask = pivot_df.isna()

        # Update the colormap to show NaN values as gray
        cmap = cm.copy()
        cmap.set_bad('gray')

        sns.heatmap(pivot_df, annot=False, cmap=cmap, 
                    vmin=cmin, vmax=cmax,
                    fmt='.2f', cbar=True,
                    cbar_kws={"orientation": "vertical", "pad": 0.005, "fraction": 0.005, "shrink":0.5},
                    linewidths=0.5, linecolor='black', mask=mask, square=True)  # Set square to True

        # Add wild type amino acids to x-axis labels with larger font size
        x_labels = [f'{col}\n{df[df["pos"] == col]["aa"].values[0]}' for col in pivot_df.columns]
        plt.xticks(np.arange(0.5, len(pivot_df.columns) + 0.5), x_labels, fontsize=10, rotation ='horizontal')
        plt.yticks(fontsize=14)

        plt.xlabel('Position', fontsize=16)
        plt.ylabel('Conditions', fontsize=16)

        # Set linewidths and adjust spines to add borders
        plt.gca().spines['top'].set_visible(True)
        plt.gca().spines['bottom'].set_visible(True)
        plt.gca().spines['left'].set_visible(True)
        plt.gca().spines['right'].set_visible(True)
        plt.gca().spines['right'].set_linewidth(1.5)
        plt.gca().spines['top'].set_linewidth(1.5)
        plt.gca().spines['left'].set_linewidth(1.5)
        plt.gca().spines['bottom'].set_linewidth(1.5)


        # Adjust the layout and save the figure
        #plt.tight_layout(rect=[0, 0.03, 0.92, 0.95])  # Adjust the layout to fit the title and colorbar
        plt.savefig(f'heatmap_{mut}_{d}.pdf', bbox_inches='tight')