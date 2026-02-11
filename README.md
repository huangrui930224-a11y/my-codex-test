# my-codex-test
chatgpt AI assistance

## MATLAB usage (only data.csv with t and v)

If you only have a waveform CSV (time + voltage), run:

```matlab
cfg = struct('M',64,'T_Np',29,'T_Dp',4,'T_Dw',4,'T_Nw',20, ...
             'save_mat',true,'output_mat','result.mat');
out = qprbs13_cei_tx_postprocess_from_csv('data.csv', [], cfg);
```

- `data.csv`: at least 2 numeric columns for `t` and `v` (column names can be `t/time` and `v_tx/v/voltage`, or any first two numeric columns).
- `t` does **not** need to start from 0.
- Cropped captured waveform is supported as long as it contains at least `20*N` UIs (`N=8191`).
- When `sym_src=[]`, symbols are automatically inferred from waveform levels.

## MATLAB usage (with symbols)

```matlab
out = qprbs13_cei_tx_postprocess_from_csv('waveform.csv', 'symbols.csv', cfg);
```
