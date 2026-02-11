# my-codex-test
chatgpt AI assistance


## MATLAB usage (CSV auto process)

Use `qprbs13_cei_tx_postprocess_from_csv` to read waveform CSV and run the full CEI/QPRBS13 flow automatically:

```matlab
cfg = struct('M',64,'T_Np',29,'T_Dp',4,'T_Dw',4,'T_Nw',20, ...
             'save_mat',true,'output_mat','result.mat');
out = qprbs13_cei_tx_postprocess_from_csv('waveform.csv', 'symbols.csv', cfg);
```

- `waveform.csv`: at least 2 numeric columns (time and voltage). If column names are present, the loader prefers `t/time` and `v_tx/v/voltage`.
- `symbols.csv`: one numeric symbol column, exactly 8191 entries in `{-1,-1/3,1/3,1}`.

