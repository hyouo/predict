# 自有数据合同与API

一套模型只对应一个目标背景、一个源背景和一个明确的效应尺度。至少4个不同校准药物；数量少时不保证性能。查询药物需要有源端实测效应，不可用SMILES或药物名称代替数值矩阵。

## 数据格式

每个NPZ由普通数值/Unicode数组和一个标量 `metadata_json` 组成。schema_version=1。程序使用 allow_pickle=False，不支持Python对象序列化。文件及解压数组总量限制512 MiB；需要更大规模时应另行修改和测试，不能假设无限扩展。

|kind|必须且仅允许的数组|
|---|---|
|train|source[N,G]、target[N,G]、ids[N]、genes[G]|
|query|source[Q,G]、ids[Q]、genes[G]；不能有target|
|truth|target[Q,G]、ids[Q]、genes[G]|

contract必须显式包含 effect_space、source_context、target_context、organism、dose、time、gene_namespace，值均为非空字符串，推理与模型必须一致。

```python
import numpy as np
from perturb_predict.io import write_bundle

# 以下数组应由真实数据预处理得到，而不是自动生成/模拟替代。
x = np.load('source_calibration_effects.npy', allow_pickle=False)
y = np.load('target_calibration_effects.npy', allow_pickle=False)
q = np.load('source_query_effects.npy', allow_pickle=False)
ids = np.load('calibration_compound_ids.npy', allow_pickle=False)
query_ids = np.load('query_compound_ids.npy', allow_pickle=False)
genes = np.load('gene_ids.npy', allow_pickle=False)
contract = {
  'effect_space': 'log2_CPM_plus1_minus_matched_control_mean',
  'source_context': 'source_A', 'target_context': 'target_B',
  'organism': 'human', 'dose': '1 uM', 'time': '24 h',
  'gene_namespace': 'Ensembl_no_version',
}
write_bundle('custom/train.npz', 'train', ids, genes, contract, source=x, target=y)
write_bundle('custom/query.npz', 'query', query_ids, genes, contract, source=q)
```

```bash
perturb-predict train --data custom/train.npz --out runs/custom-model
perturb-predict predict --model runs/custom-model/model.npz --data custom/query.npz --out runs/custom-prediction
```

后来取得目标查询实验结果后，另建truth文件，以 `perturb-predict evaluate --prediction runs/custom-prediction/predictions.npz --truth custom/truth.npz --out runs/custom-score` 评分。评分按完整药物与基因集合对齐，禁止少报几个难例。不同目标背景应分别训练，不要把两套校准混在同一个ID空间。

## 外部预处理的责任

通用API能够验证形状、标识、契约、训练查询不重叠、有限数值；它**不能验证你提供的效应是否来自独立实验、是否正确归一化、是否已经利用测试答案筛基因**。供体/板/孔划分、对照匹配、剂量、采样时间、基因选择和组成变化必须由上游分析记录。缺失效应不能用零填补；零只有在它确实是观测值/定义值时才合法。

OP3适配器只接受固定哈希文件。其他H5AD应先转换为上述显式效应包，不能通过改扩展名绕过数据审计。
