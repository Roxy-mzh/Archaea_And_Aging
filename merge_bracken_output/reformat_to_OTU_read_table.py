
import pandas as pd
import numpy as np
bracken_merged = pd.read_csv ('bracken_merged.csv')
bracken_wide_table = pd.pivot(bracken_merged, index=['name'], columns = 'sample_id',values ='new_est_reads') #Reshape from long to wide
bracken_wide_table.to_csv ('wide_bracken_merged.csv')
