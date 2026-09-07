# Predict — v0.1.0 可用研究基线

跨细胞背景的化学扰动定量响应预测。现在提供可安装的 `perturb-predict` 命令、可保存/独立加载的模型、严格的数据合同，以及真实 OP3 的可重复验收。

**可用范围：单源、少样本目标校准、已在源背景测过的查询药物。不是全新药物/完全未见细胞背景零样本，不是临床产品，不宣称 SOTA。**

## 安装与运行

```bash
python -m pip install -r requirements-release.txt
python -m pip install .
perturb-predict --version
perturb-predict fetch-op3 --out data/rna/op3_standardized_processed.h5ad
perturb-predict benchmark-op3 --data data/rna/op3_standardized_processed.h5ad \
  --out runs/usable-001 --accept-conditional --accept-public-reuse
```

`--out` 必须不存在。已有固定哈希 OP3 文件可直接复用。网络失败明确报错，不使用合成数据替代。发行验收使用 Python 3.13.5、NumPy 2.3.5、h5py 3.15.1；无需 GPU。没有发布到 PyPI。

## 文档

- [快速开始与分步CLI](docs/releases/QUICKSTART_ZH.md)
- [自有数据格式与Python API](docs/releases/CUSTOM_DATA_ZH.md)
- [预先提交的发行验收标准](docs/releases/v0.1.0-acceptance-plan.md)
- [验收结果与局限](docs/releases/VALIDATION_ZH.md)
- [实际评分JSON](reports/releases/v0.1.0/scores.json)

模型采用全谱的单源部分汇聚斜率与截距收缩；训练、预测、评分三个操作互相分离。推理仅需模型和不含目标标签的查询包。完整药物留一选参数，按基因ID对齐，拒绝缺失/重复ID和非有限值；读取模型禁用pickle；输出附SHA256、参数和执行状态。

旧研究代码保留在 `src/` 和 `tools/`，其协议与历史分数不变。v0.8不同校准预算的两条本地研究轨迹不混合；本次明确提炼最新提供的11/10校准药物轨迹。产品版本v0.1.0与历史研究版本不是同一编号序列。

## 验收与科研边界

本地51项发行测试及50项原仓库测试通过；安装后的wheel完成两个目标的真实OP3回放和独立推理；与旧v0.8对应估计器数值一致。当前GitHub执行结果以Actions为准。

本轮公开OP3平均MSE 0.198067，对照源复制0.219748；这是重复使用公开开发数据的工程验收，不是独立生物学验证。上游基因筛选、同孔源/目标共享、供体映射缺失等问题未被包装解决。下一科学关卡仍是独立研究/实验条件下的验证，以及目标基线调控信息是否带来增量。

不配置定时科研循环、付费GPU或自动发布科学结论。
