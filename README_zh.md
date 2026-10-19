# APLE-LM 开源复现代码

[English](README.md) · [MIT 许可证](LICENSE) · [论文引用](CITATION.cff)

本包覆盖目前选定的博士论文正文图 M1–M8。采用已审查的定位数值实现，
保留低网格、低 RF 链等压力场景中的真实失败结果。软件作者署名为 **Boyu Teng**；
论文五位作者列在引用文件中。原数值源文件的作者注释完整保留。

本次 v0.2 将主图按期刊顺序重新编号，阵列分组与距离图已拆开。
[新旧编号对照](docs/FIGURES.md)列明对应关系。脚本可使用 `grid`、`snr`、`cdf`、
`aperture`、`rf`、`range`、`coverage`、`iterations` 等主题名称，避免依赖图号。

## 直接查看和重绘已有结果

将本目录设为 MATLAB 当前目录：

```matlab
setup_aple;
plot_results('reference','Language','zh');
```

不需要重新仿真，也不需要历史 MAT 文件。`reference/figures/` 已提供八张中文正文 PDF；
上述命令由 `reference/summary.csv` 与 `reference/cdf.csv` 重绘并输出到
`runs/reference_plots/`。默认导出 PDF 和 PNG；需要可编辑 MATLAB 图时增加
`'SaveFig',true`。图中不放总标题，实验名称和条件见 [图表说明](docs/FIGURES.md)。

## 从零运行实验

已在 MATLAB R2023b 验证；兼容目标为 R2021a 及以上，其他版本尚未验证。
串行计算只需要 MATLAB 及其 JVM。并行选项需要 Parallel Computing Toolbox。
中文图需要宋体等中文字体；没有中文字体时使用 `'Language','en'`。

```matlab
% 首次使用建议：默认阵列、1 个种子、6 个信噪比、4 种方法。
out = run_experiments('M2','smoke');
plot_results(out,'Figures',{'M2','M3'},'Language','zh');

% 所有实验的少量代表配置，主要用于检查流程。
run_experiments('all','smoke');

% 完整复现：所有配置，每配置 1000 个种子。
out = run_experiments('all','full');
plot_results(out,'Language','zh');

% 单独运行 RF 链实验；可选并行，默认每批保存一次。
out = run_experiments('M5','full','UseParallel',true,'Workers',2);
plot_results(out,'Figures',{'M5'},'Language','zh');
```

完整实验共有 181,000 个“试验种子 × 配置”作业，其中热图占 143,000 个，
需要较长计算时间。`smoke` 只减少种子和配置数量，保留算法迭代预算，不作为性能结论。
重复同一命令即可断点续跑；完成的作业会跳过。更改种子、数值源码或 MATLAB 版本时
使用新的 `OutputDir`，避免混用结果。批大小和 worker 数量可以调整。

## 文件和存储策略

- `src/`：46 个数值依赖文件，与当前审查后的实现一致。
- `experiments/`：统一实验配置、按种子生成输入、断点保存和统计。
- `plotting/`：M1–M8 中英文绘图。
- `reference/`：轻量 CSV 参考结果、元数据和八张 PDF；不包含 MAT。
- `docs/`：实验协议、技术改动、指标定义和验证范围。
- `runs/`：本机新生成的结果，Git 忽略，发布包不包含。

每个实验族只保存一份紧凑检查点和一份恢复备份，写入期间临时增加一份待提交文件。
仅持久化逐次标量误差、理论界、完成状态和少量运行信息。不保存观测、信道矩阵、
字典、Hessian、逐次定位轨迹、源码快照或每次试验的单独文件。
输入和固定字典在内存中按种子重建；续跑最多重算未提交的一批。
详见 [存储与恢复](docs/STORAGE.md)。

## 结果口径

M1、M4、M5 使用单次试验 RMS ≤ 3 m 的条件 RMSE，未将大误差数值裁剪成 3 m；
M1/M5 同时给出全样本分母下的超限概率，不绘制区间。
完整 RMSE、保留样本数等仍在汇总 CSV 中。M3 的 CDF 保留所有样本。
M4 保留五种分组和三条 BCRB；M5 左图保留三条算法线与两条 BCRB；
M8 上下两个子图显示前 10 次真实迭代，统计中保留全部 20 次；M7 使用原有 13×11 扫描点。

为复现已完成结果，本包保留历史 `2/Ns` 用户分离参数和 OR 接受判据。
Ns=9、20、25 与 Ns=15 的门限不同，不能表述为统一固定门限的公平性补实验。
低 RF 链失效机制没有在本次开源整理中修改；相关限制见
[实验协议](docs/REPRODUCIBILITY.md) 和 [技术调整](docs/NUMERICAL_NOTES.md)。

```matlab
setup_aple;
results = runtests('tests/test_release.m');
assert(all([results.Passed]));
```

[验证说明](docs/VALIDATION.md) 区分本次实际执行的抽样回归与此前完整的 1000 次参考结果。
许可证为 MIT，版权行使用 `Copyright (c) 2026 Boyu Teng`。
