# 快速开始

从仓库的 `main` 使用 `perturb-predict`。本页依次介绍安装、离线演示、自有数据入口和公开 OP3 工程重放。无需在历史研究分支之间选择；当前软件版本可用 `perturb-predict --version` 查看。

模型的适用条件是：已知一个源背景中一批药物的定量表达效应，目标背景中有至少 4 种不同校准药物的表达效应，查询药物也已在源背景中测量。基因、物种、剂量、时间与效应定义需要一致。原始计数必须先经过明确的归一化和对照匹配，不能直接作为效应矩阵输入。

## 安装

需要 Python 3.11 或更高版本。软件不依赖 GPU；普通安装只需要 `pyproject.toml` 声明的依赖。

```bash
git clone https://github.com/hyouo/predict.git
cd predict
python -m venv .venv
```

Linux / macOS：

```bash
source .venv/bin/activate
```

Windows PowerShell：

```powershell
.venv\Scripts\Activate.ps1
```

然后在仓库根目录安装当前代码：

```bash
python -m pip install .
perturb-predict --version
perturb-predict --help
```

如果 Windows 的执行策略不允许激活脚本，可以直接使用 `.venv\Scripts\python.exe -m pip install .`，之后以 `.venv\Scripts\perturb-predict.exe` 调用命令，无需改变系统执行策略。

本项目提供从源码和 wheel 安装的流程，不声明已经发布到 PyPI。不要把 `pip install perturb-predict` 当作本文验证的安装来源。

[Package CI](../../.github/workflows/package.yml) 配置了 Linux 上的 Python 3.11–3.14 与 Windows 上的 Python 3.13 检查；通过情况以 [Actions 的具体运行](https://github.com/hyouo/predict/actions/workflows/package.yml) 为准。历史 OP3 验收另用 Python 3.13.5 与 `requirements-release.txt` 中的固定依赖；普通安装不需要先安装历史研究的完整依赖集。

## 离线操作演示

安装后执行：

```bash
perturb-predict demo --out runs/demo-001
```

该命令不访问网络，使用固定种子构造 6 种校准药物、3 种查询药物与 8 个基因的纯合成数组，再依次运行训练、预测和评分。所有名称、效应与数据合同均标记为合成；生物学样本数为 0。**它用于理解软件流程和检查安装，不能用于估计真实生物学性能。**

主要输出如下：

| 文件 | 用途 |
| --- | --- |
| `protocol.json` | 合成配方、种子、划分、信息预算与运行环境 |
| `inputs/train.npz` | 含源端和目标端校准效应的训练包 |
| `inputs/query.npz` | 仅含源端查询效应的查询包 |
| `inputs/truth.npz` | 单独保存的合成查询真值，只传给评分步骤 |
| `model/model.npz`、`model/fit.json` | 可独立加载的模型和拟合记录 |
| `prediction/predictions.npz`、`prediction/predictions.csv`、`prediction/freeze.json` | 预测、可查看的表格及冻结记录 |
| `evaluation/metrics.json`、`evaluation/per_compound.csv` | 模型、源响应复制与零效应对照的合成评分 |
| `summary.json` | 汇总、明确的合成证据标记与产物 SHA256 |
| `COMPLETE.json` | 当前运行完成标记，各阶段也有独立标记 |

同一记录的数值环境下，生成规则与数值结果可复现；时间戳等运行元信息会变化，不应要求两个运行目录逐字节相同。

所有运行输出目录必须不存在。重新运行时改用 `runs/demo-002` 等新路径，不要提前创建该目录。运行阶段失败时会尽可能保留 `FAILED.json`；发生进程终止或存储故障时可能来不及写入失败记录，因此目录存在本身不表示运行完成。

## 分步练习与自有数据

离线演示产生的输入可以直接用于分步练习；每个命令使用一个新输出目录：

```bash
perturb-predict train --data runs/demo-001/inputs/train.npz --out runs/step-model-001
perturb-predict predict --model runs/step-model-001/model.npz --data runs/demo-001/inputs/query.npz --out runs/step-prediction-001
perturb-predict evaluate --prediction runs/step-prediction-001/predictions.npz --truth runs/demo-001/inputs/truth.npz --out runs/step-score-001
```

独立推理只需要 `model.npz` 和不含目标查询标签的 `query.npz`。评分还需要预测文件同目录中的 `freeze.json` 与完成标记；移动预测时应保留整个预测目录。`freeze.json` 用于检查产物变更，不是加密签名或第三方认证。

接入实验数据请阅读 [自有数据合同与 Python API](CUSTOM_DATA_ZH.md)。`source` 和 `target` 都是“不同药物 × 基因”的效应矩阵；训练矩阵必须逐行对应同一药物、逐列对应同一基因。不能把重复孔改名为不同药物以增加校准样本数。模型只接受一个源背景，不自动挑选最优源。

常用 API：

```python
from perturb_predict import TransferModel
from perturb_predict.io import write_bundle
```

输入基因集合相同但顺序不同时会按 ID 对齐。缺失、多余或重复 ID、NaN、无穷值，以及数据合同不一致均拒绝，不能用零填补缺失观测来绕过检查。

## 公开 OP3 工程重放

这一步需要实测公开数据，独立于离线演示。它重放项目已经用于开发的 OP3 数据，不产生新的独立生物学验证。原有模型与验收记录见 [VALIDATION_ZH.md](VALIDATION_ZH.md)。

若目标是重放原固定环境，请另建 Python **3.13.5** 环境，并在仓库根目录安装：

```bash
python -m pip install -r requirements-release.txt
python -m pip install .
```

为减少数值线程差异，Linux / macOS 设置 `export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1`；Windows PowerShell 的对应设置见 [贡献指南](../../CONTRIBUTING.md#开发环境与软件检查)。

取得固定文件（约 23.7 MB），然后运行重放：

```bash
perturb-predict fetch-op3 --out data/rna/op3_standardized_processed.h5ad
perturb-predict benchmark-op3 --data data/rna/op3_standardized_processed.h5ad --out runs/op3-replay-001 --accept-conditional --accept-public-reuse
```

已有同哈希文件时可直接复用。`fetch-op3` 会核对固定 SHA256；下载失败或缓存不匹配时明确报错，不使用合成数据替代。离线环境可直接通过 `--data` 指向已经取得的同一固定文件。使用数据前请核查 [第三方说明](../../THIRD_PARTY_NOTICES.md) 与 [数据署名](../../data/ATTRIBUTION.md)。

两个开关分别确认上游已筛选基因面板等条件，以及公开标签已在开发中使用。运行会先冻结两个目标的预测，再开启公开评分；拒绝读取 private 测试表达。输出包含两个目标的模型、预测、源响应复制/零效应对照、逐药物误差、访问记录与 SHA256。

也可分步操作：

```bash
perturb-predict prepare-op3 --data data/rna/op3_standardized_processed.h5ad --out runs/op3-prepared-001 --accept-conditional
perturb-predict train --data runs/op3-prepared-001/B_cells/train.npz --out runs/op3-model-001
perturb-predict predict --model runs/op3-model-001/model.npz --data runs/op3-prepared-001/B_cells/query.npz --out runs/op3-prediction-001
perturb-predict score-op3 --data data/rna/op3_standardized_processed.h5ad --prediction runs/op3-prediction-001/predictions.npz --out runs/op3-score-001 --accept-public-reuse
```

正式 OP3 评分要求完整的查询集合，不能只提交容易预测的药物。上游基因筛选、供体映射缺失、同孔源/目标共享及参考误差仍限制解释；软件包装没有解决这些科学问题。

## 常见问题

| 现象 | 处理 |
| --- | --- |
| 找不到 `perturb-predict` 命令 | 确认激活了安装所用的虚拟环境，也可运行 `python -m perturb_predict --help`。 |
| 输出目录已存在 | 指定新的 `--out` 路径，保留原运行用于排查。 |
| 文件哈希不匹配 | 核对资产来源与文件完整性；不要修改代码里的固定哈希绕过验证。 |
| 数据合同或 ID 不一致 | 按 [自有数据合同](CUSTOM_DATA_ZH.md) 检查效应尺度、背景、物种、剂量、时间与 ID 集合。 |
| 环境安装失败或运行报错 | 提交最小复现、安装命令、版本与完整错误；入口见 [贡献指南](../../CONTRIBUTING.md)。 |

## 结果解释

`mse` 是规定效应尺度下的平方误差；`centered_mse + bias_mse = mse`。`retrieval_top1` 衡量药物响应匹配，不是基因正确率。`relative_source_distance` 是源响应到校准中心的描述性距离，大于 3 的提醒为固定启发式，不是置信概率。

模型需要源端实测查询响应与目标校准，不能用于声称全新药物或全新背景零样本能力。当前不提供经生物学重复校准的置信区间。研究路线与尚未满足的独立证据关卡见 [当前状态](../../STATUS.md)。
