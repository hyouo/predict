# Predict

可安装、可审计的化学扰动响应预测研究基线。给定药物在一个源背景中的实测表达效应，以及目标背景中至少 4 种校准药物的实测效应，Predict 可以训练模型、保存模型，并预测其他已在源背景测过的药物在目标背景中的效应。

**从 `main` 开始使用和贡献。** `perturb-predict` 是当前维护的软件入口；历史研究材料用于追溯，不需要切换研究分支才能使用。软件版本号仅用于追溯安装包，历史研究编号不代表其他维护版本。

## 快速开始

需要 Python 3.11 或更高版本，无需 GPU。先克隆仓库并创建独立环境：

```bash
git clone https://github.com/hyouo/predict.git
cd predict
python -m venv .venv
```

Linux / macOS 使用 `source .venv/bin/activate`；Windows PowerShell 使用 `.venv\Scripts\Activate.ps1`。然后安装当前检出的代码并运行离线演示：

```bash
python -m pip install .
perturb-predict --version
perturb-predict demo --out runs/demo-001
```

`demo` 使用确定性的纯合成数据，不下载实验数据；输出模型、预测、评分和 `runs/demo-001/summary.json`，用于确认安装和操作流程。**合成演示分数不构成生物学证据。** 每次运行须使用不存在的输出目录，例如下次改为 `runs/demo-002`。

安装依赖由 [pyproject.toml](pyproject.toml) 声明。[软件兼容性工作流](.github/workflows/package.yml) 覆盖 Python 3.11–3.14 与 Windows 检查，实际运行状态见 [GitHub Actions](https://github.com/hyouo/predict/actions)。固定 Python 3.13.5 的历史复现环境另行保留。本文提供从源码安装的方法，不声明已发布到 PyPI。

## 能做什么

- 分开训练、预测和评分；独立推理只需模型与不含目标查询标签的输入包。
- 按药物和基因 ID 对齐，校验效应尺度、物种、背景、剂量和时间；拒绝缺失或重复 ID、非有限值及训练/查询重叠。
- 使用不含 pickle 的数值模型文件，记录输入与产物 SHA256、拟合参数和运行完成或失败状态。
- 提供固定哈希的公开 OP3 数据适配器，以及显式的自有数据合同。

该模型需要源端实测响应和目标端少样本校准，不能据此声称支持全新药物或完全未测目标背景的零样本预测。它是研究软件，不提供临床决策或已校准的预测置信区间。缺失观测不能补作生物学零值。

## 下一步

| 目的 | 入口 |
| --- | --- |
| 安装、离线演示与分步命令 | [快速开始](docs/releases/QUICKSTART_ZH.md) |
| 接入自己的数据、调用 Python API | [自有数据合同](docs/releases/CUSTOM_DATA_ZH.md) |
| 重放公开 OP3 工程验收 | [快速开始中的 OP3 流程](docs/releases/QUICKSTART_ZH.md#公开-op3-工程重放) |
| 理解已记录结果及其局限 | [历史验收记录](docs/releases/VALIDATION_ZH.md) |
| 查看本次开源工程验收范围 | [开源就绪验收](docs/releases/OPEN_SOURCE_READINESS_ZH.md) |
| 报告问题或贡献代码 | [贡献指南](CONTRIBUTING.md) · [Issues](https://github.com/hyouo/predict/issues) |
| 报告安全问题 | [安全策略](SECURITY.md) |

公开 OP3 曾用于方法开发，重放它只验证工程可复现性。独立制备或外部研究中的有效性仍是待解决的科学关卡，见 [当前状态](STATUS.md) 和 [研究路线](research/ACTIVE_ROUTE_ZH.md)。历史材料见 [历史归档](archive/README.md)；归档的存在不代表其中所有实现或结论均通过当前软件验收。

## 许可与数据来源

项目原创软件及对应文档采用 [MIT 许可](LICENSE)。第三方数据、第三方材料及历史归档不由该许可重新授权，具体范围见 [第三方与历史材料说明](THIRD_PARTY_NOTICES.md) 和 [数据署名](data/ATTRIBUTION.md)。使用或再分发实验数据前，应核对其原始条款。
