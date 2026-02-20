import pandas as pd
import numpy as np

file = 'peri_gluc.csv'

df_dynamic = pd.read_csv(file)

# normalization:
# 1) convert insulin secretion from µU/mL to pmol/L (1 µU/mL = 6 pmol/L)
# 2) divide the insulin secretion to value per islet (65 islets per sample)
# 3) substract the average baseline (3mM)

df_norm = df_dynamic.loc[:,[col for col in df_dynamic.columns if col.startswith('time_') or col.startswith('record_')]].groupby('record_id').mean().mul(6).div(65)

baselines = list(df_norm.loc[:,['time_5', 'time_10', 'time_15']].mean(axis=1))

for i, idx in enumerate(df_norm.index):
    df_norm.loc[idx,:] = df_norm.loc[idx,:].sub(baselines[i]).clip(lower=0)


df_norm.to_csv(file.replace('.csv', '.norm_TdW.csv'))

df_summary = pd.DataFrame(index=list(df_norm.index), columns=['peak_15', 'peak_6', 'peak_KCl', 'auc_15', 'auc_6', 'auc_KCl'])

x = [float(col.strip('time_')) for col in  df_norm.columns]

def sub_section(y, t1, t2, t=x):
    t = np.array(t)
    dy = y[(t>=t1) & (t<= t2)]
    dt = t[(t>=t1) & (t<= t2)]
    return dt, dy

for donor in df_norm.index:
    t = np.array(x)
    y = np.array(df_norm.loc[donor])

    #15
    dt, dy = sub_section(y, 20, 65)

    df_summary.loc[donor,'peak_15'] = np.max(dy)
    df_summary.loc[donor,'auc_15'] = np.trapz(dy, dt)

    #6
    dt, dy = sub_section(y, 90, 140)

    df_summary.loc[donor,'peak_6'] = np.max(dy)
    df_summary.loc[donor,'auc_6'] = np.trapz(dy, dt)

    #KCl
    dt, dy = sub_section(y, 160, 180)

    df_summary.loc[donor,'peak_KCl'] = np.max(dy)
    df_summary.loc[donor,'auc_KCl'] = np.trapz(dy, dt)

df_summary.to_csv(file.replace('.csv', '.norm_TdW.summary.csv'))