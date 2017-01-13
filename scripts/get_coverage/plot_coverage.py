import os
import math
import pandas as pd
import sys
import numpy as np

src=sys.argv[1]
dest=sys.argv[2]

for subdir, dirs, files in os.walk(src):
    for file in files:
        df = pd.DataFrame()
        try:
            df = pd.read_table(src+file, sep='\t', names=['id', 'position', 'depth'])
        except:
            print(src+file)
            continue
        if df.empty:
            continue
        log_scale = []
        for i in df['depth']:
            if i != 0:
                log_scale.append(math.log(i))
            else:
                log_scale.append(0)
        df['log(depth)'] = log_scale
        depth = np.mean(df['depth'])
        per_cov = (len(df[df['depth']>=1])/len(df['depth'])) * 100
        print('|'+file+' | '+str(per_cov)+'|' +str(depth)+ '|')
        # plt.clf()
        # df['log(depth)'].plot()
        # plt.xlabel('Position', fontsize=16)
        # plt.ylabel('ln(depth)', fontsize=16)
        # plt.savefig(dest+file.replace(".tsv", ".png"))
