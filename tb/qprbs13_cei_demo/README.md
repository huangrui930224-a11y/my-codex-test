# qprbs13_cei_gen 使用示例（Spectre/Verilog-A）

下面给出一个最小化思路，说明如何实例化 `qprbs13_cei_gen`、如何给时钟、以及如何观察 `A/B/pam4`。

## 1) 实例化

```verilog
`include "disciplines.vams"
`include "constants.vams"
`include "va/blocks/qprbs13_cei_gen.va"

module qprbs13_cei_tb;
  electrical clk, rst_n, A, B, pam4;

  // 时钟与复位源（示意）
  // 可换成 analogLib/vpulse 等任何可在 Spectre 中使用的激励
  parameter real vclk_lo = 0.0;
  parameter real vclk_hi = 0.75;
  parameter real t_ui    = 20p;

  qprbs13_cei_gen #(
    .v0(0.0), .v1(0.25), .v2(0.5), .v3(0.75),
    .vbit_low(0.0), .vbit_high(0.75),
    .tr(2p), .tf(2p),
    .vth(0.375),
    .seed(1),
    .gray_enable(1),
    .pam4_enable(1)
  ) dut (
    .clk(clk),
    .rst_n(rst_n),
    .A(A),
    .B(B),
    .pam4(pam4)
  );

  analog begin
    // 50% 占空比时钟
    V(clk) <+ transition(((($abstime % t_ui) < (t_ui/2.0)) ? vclk_hi : vclk_lo), 0, 1p, 1p);

    // 上电后一段时间拉低复位，再释放
    V(rst_n) <+ transition(($abstime < 5*t_ui) ? 0.0 : 0.75, 0, 1p, 1p);
  end
endmodule
```

## 2) 观察要点

- `A` / `B`：每个 `clk` 上升沿更新一次，电平为 `vbit_low` 或 `vbit_high`。
- `pam4`：当 `pam4_enable=1` 时，每个 `clk` 上升沿更新一次，取值为 `v0/v1/v2/v3` 之一。
- 若 `gray_enable=1`，则 `A/B` 先进行 Gray 映射（`A=A_bin`, `B=A_bin^B_bin`）后再输出与做 PAM4 映射。

## 3) 建议瞬态仿真设置

- `tran` 仿真时间至少覆盖几百个 UI（如 `stop=200*t_ui`）。
- 保存 `V(clk), V(rst_n), V(A), V(B), V(pam4)`。
- 用游标检查在每个时钟上升沿后，`A/B/pam4` 是否按预期跳变。
