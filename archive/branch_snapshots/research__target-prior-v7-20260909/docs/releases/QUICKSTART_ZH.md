# Predict v0.1.0：第一个可用研究版本

这是可以安装、保存模型并独立推理的**单源、少样本化学扰动响应预测基线**，不是新药物零样本、全新细胞背景、临床产品或 SOTA 声明。研究版本 v0.8 与产品版本 v0.1.0 是不同版本序列。

## 适用场景

已知一批药物在源背景中的定量表达效应；目标背景有至少4种不同校准药物的表达效应；希望预测已在源背景测量、但尚未在目标背景测量的其他药物。相同基因、同种属、同剂量、同时间、同效应定义必须对齐。原始计数需先经过明确的归一化和对照匹配，不能当效应矩阵直接输入。

固定模型是 v0.8 的 T-only 部分汇聚斜率＋合并对照截距收缩。没有增加谱截断，也没有增加目标测试标签。逐基因斜率向共享斜率收缩，使用完整药物留一选择25组惩罚/截距候选；没有单细胞数充当独立药物样本的交叉验证。

## 安装

```bash
python -m venv .venv
# Linux/macOS:
source .venv/bin/activate
# Windows PowerShell 使用 .venv\Scripts\Activate.ps1
python -m pip install -r requirements-release.txt
python -m pip install .
perturb-predict --version
```

本地发行验收使用 Python 3.13.5、NumPy 2.3.5、h5py 3.15.1。声明的 Python 最低版本为3.11，其他组合未经本轮全面验收。包不需要 GPU；wheel 可用 `pip install perturb_predict-0.1.0-py3-none-any.whl` 安装。没有发布到 PyPI，不要把包名直接当成已存在的 PyPI 项目。

## 一条命令运行真实 OP3 验收

先取得固定文件（约23.7 MB；已有同哈希文件可以直接复用，不会重新下载）：

```bash
perturb-predict fetch-op3 --out data/rna/op3_standardized_processed.h5ad
perturb-predict benchmark-op3 \
  --data data/rna/op3_standardized_processed.h5ad \
  --out runs/usable-001 \
  --accept-conditional --accept-public-reuse
```

也可直接使用已提供的数据包里的 `data/rna/op3_standardized_processed.h5ad`，避免网络下载。不要重复使用已经存在的 `--out` 目录。

输出含两个目标细胞模型、每个药物的5288基因效应预测、宽表CSV、与源复制/零效应比较的评分、逐药物误差、输入读取记录、代码/输入/模型/预测 SHA256、完整拟合参数和 COMPLETE/FAILED 状态。

两组预测全部生成后才开启公开评分。私有测试表达始终拒绝读取。公开数据此前已经用于方法开发，本命令是工程验收/复现，不是新的独立验证。上游基因筛选、供体映射缺失、同孔源/目标共享和参考误差仍然限制生物学推断。

## 分步训练与无目标标签推理

```bash
perturb-predict prepare-op3 --data data/rna/op3_standardized_processed.h5ad \
  --out runs/prepared-001 --accept-conditional
perturb-predict train --data runs/prepared-001/B_cells/train.npz --out runs/model-001
perturb-predict predict --model runs/model-001/model.npz \
  --data runs/prepared-001/B_cells/query.npz --out runs/prediction-001
perturb-predict score-op3 --data data/rna/op3_standardized_processed.h5ad \
  --prediction runs/prediction-001/predictions.npz --out runs/score-001 --accept-public-reuse
```

单独推理只需要 `model.npz` 和无目标标签的 `query.npz`，无需原始HDF5，也无需训练目录。`predictions.npz` 的列顺序是模型基因顺序；输入基因同集合不同顺序时自动对齐。缺基因、多基因、重复ID、NaN或剂量/时间/细胞背景不一致均拒绝，不会静默补零。评分时同时需要同目录的 `freeze.json`；它用于检测预测文件变更，不是加密签名或防篡改第三方认证。

## 使用自己的数据

见 `CUSTOM_DATA_ZH.md`。Python API 可直接导入：

```python
from perturb_predict import TransferModel
from perturb_predict.io import write_bundle
```

`source`和`target`均为“不同药物×基因”的效应矩阵，而非原始计数；训练两矩阵必须逐行同药物、逐列同基因。每个药物先按你自己的实验合同聚合；不能把多孔简单命名成不同药物以增加样本数。模型只接受一个源背景，不自动猜测哪个源最好。

## 结果怎样解释

MSE 衡量规定效应尺度下的平方误差；centered_mse＋bias_mse＝mse；retrieval_top1 衡量药物响应匹配，不是基因正确率。`relative_source_distance` 只是源响应到校准中心的描述性距离；>3的提醒是固定启发式，不是置信概率。没有提供经生物重复校准的可信区间。

所有结果仅供研究。没有远端自动研究循环、付费GPU或自动发布科学结论。
