import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_ind
import os

# read all .xlsx files in directory into pd dfs 
data_dir = './data/Shrock2020_all-data/'
main_file_name = 'Shrock_2020.xlsx'
main_df = pd.read_excel(os.path.join(data_dir, main_file_name), sheet_name=None)
main_df
