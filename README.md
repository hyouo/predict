# predict — 跨细胞背景扰动预测

本仓库是后续研究的版本化工作位置：保存算法、统计假设、数据来源、实验协议、负结果及复现入口。当前接入基线为 **v0.7 核心代码**，不是一次新的算法领先声明。

## 当前证据

现有数值研究主要来自辅助磷蛋白面板：5 个背景、每背景 72 个保留条件、18 个读数。它是少样本目标校准，不是仅凭未处理基线的单细胞 RNA 零样本预测。关键历史汇总位于 `reports/v0.7/`；解释和限制见 [STATUS.md](STATUS.md)。

原研究包的 15 个源代码/测试/入口文件按字节接入，原 37 项测试保留；新增 8 项仓库工具测试。原始归档 SHA256 和文件映射在 `docs/history/v0.7_import.json`。完整历史数组及其他旧研究流水线仍在原归档中，不声称已全部迁移或全部重跑。

## 从干净环境运行

```bash
git clone https://github.com/hyouo/predict.git
cd predict
# 合并前使用 research/bootstrap-v0.7 分支；合并后直接使用 main。
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
python -m tools.assets --asset auxiliary
python -m unittest discover -s tests -v
python -m tools.smoke --output runs/my-first-check
```

`smoke` 用真实辅助数据做一次固定划分的集成检查，记录代码/数据/预测哈希、环境、划分和 Git 提交，并比较分块与非分块实现。它不是独立生物验证；输出目录存在时拒绝覆盖。

已有本地原始数据时：

```bash
python -m tools.assets --asset auxiliary --source /path/to/hepatocyte_signaling_raw.csv
```

固定资产的 SHA256 不匹配、下载失败或数据缺失都会报错，不会偷偷使用模拟数据。

## RNA 数据入口

```bash
python -m tools.assets --asset op3
# sci-Plex 文件更大，仅显式运行：
python -m tools.assets --asset sciplex
```

资产版本、来源和大小上限在 `assets/manifest.json`。HDF5 结构读取成功只说明文件可读，**不代表原始计数、归一化、供体/孔/重复、对照或数据划分语义已经核验，更不代表模型已训练**。

GitHub Actions 的 CI 用于测试和一次小型集成检查；独立的 OP3 工作流尝试获取固定公开文件、核验哈希和输出结构审计。运行是否成功以 Actions 日志为准。大矩阵只作短期 artifact，不进入 Git；没有定时训练、付费 GPU 或自主后台科研进程。

## 研究和协作

先读 [研究状态](STATUS.md)、[路线图](docs/ROADMAP.md) 和 [工作约定](AGENTS.md)。新增研究使用独立分支和 PR，保留负结果及预算变化，不覆写既有实验。

原始入口 `run_covariance.py`、`run_sentinel_probe.py` 为未改写的历史程序，可能覆写其输出目录；只在独立工作树/全新目录运行，正式迭代优先使用有运行记录的入口。`python -m tools.verify_import` 仅核验初始化基线；以后算法修改造成差异是正常现象，不应修改历史哈希来掩盖变化。
