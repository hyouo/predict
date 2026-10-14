# 贡献指南

欢迎提交可以重现的问题、文档修正和范围明确的代码改进。请以 `main` 为基准；当前维护的安装包位于 `perturb_predict/`。历史研究归档供追溯使用，不作为新的软件入口。

## 开始前

先阅读 [README](README.md)、[AGENTS.md](AGENTS.md)、[当前状态](STATUS.md)、[路线图](docs/ROADMAP.md) 和相关的开放 [Issues](https://github.com/hyouo/predict/issues)。算法、数据或科研结论的变更还应阅读 [研究路线](research/ACTIVE_ROUTE_ZH.md) 与 [整合政策](research/INTEGRATION_POLICY.md)。

普通修复可以直接提交 Pull Request。较大的 API、数据合同或默认算法变更，请先用 issue 描述问题、预期行为和验证方式。不要把多个无关改动或新的科学主张塞进同一个修复。

## 开发环境与软件检查

使用 Python 3.11 或更高版本，在仓库根目录执行：

```bash
python -m venv .venv
```

Linux / macOS 激活环境：

```bash
source .venv/bin/activate
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
```

Windows PowerShell 激活环境并设置线程数：

```powershell
.venv\Scripts\Activate.ps1
$env:OPENBLAS_NUM_THREADS = "1"
$env:OMP_NUM_THREADS = "1"
$env:MKL_NUM_THREADS = "1"
```

安装可编辑包，运行发行软件测试和纯合成离线演示：

```bash
python -m pip install -e .
python -m unittest discover -s tests_release -v
perturb-predict demo --out runs/contribution-demo-001
```

测试使用 Python 标准库 `unittest`。演示不会下载数据，测试和演示的通过只说明实现行为符合检查。不要将合成结果用于支持生物学性能主张。

涉及打包的改动还需构建 wheel，并在仓库目录之外、独立虚拟环境中安装后运行 `perturb-predict --version`、`perturb-predict demo` 和发行测试，避免源码路径掩盖安装问题。可在仓库根目录构建：

```bash
python -m pip wheel --no-deps . --wheel-dir dist
```

再将生成的 wheel 安装到独立环境；完整的自动检查步骤见 [package.yml](.github/workflows/package.yml)。Python 版本矩阵和操作系统结果以本次 PR 的 Actions 为准。

## 历史基线和研究相关检查

涉及 `src/`、`tools/` 或旧基线的改动，应另用 Python **3.13.5** 环境，保留上述单线程设置，并在仓库根目录执行：

```bash
python -m pip install -r requirements.txt
python -m tools.assets --asset auxiliary
python -m unittest discover -s tests -v
python -m tools.smoke --output runs/contribution-smoke-001
```

辅助资产会从固定公开来源下载并核对 SHA256。有已取得的原始文件时，可用 `python -m tools.assets --asset auxiliary --source /path/to/original.csv` 导入；仍必须通过同一哈希校验。下载失败须报告，不能用合成数据替代实测基线。每次集成检查使用新的输出目录。

涉及桥接设计模块的改动，另运行：

```bash
python -m unittest discover -s studies/bridge_design_v14/tests -v
```

涉及 OP3 读取、训练、预测或评分的改动，需要评估 [固定 OP3 工程重放](docs/releases/QUICKSTART_ZH.md#公开-op3-工程重放) 和 [既有发行验收协议](docs/releases/v0.1.0-acceptance-plan.md) 中的相关检查。网络或资源限制导致未运行时，明确记录未运行的检查与原因，不据此声称通过。

## 数据与科学变更

- 训练和调参只使用允许的校准信息；目标查询标签不能决定预处理、调参、探针或模型选择。
- 明确区分未见配对、未见背景、未见药物及少样本校准，记录所有目标处理、对照、探针、重复和失败的预算。
- 不把缺失标签补成零，不把基因数、掩码种子或共享对照产生的观测当作独立生物学重复。
- 新科学运行应先提交协议，再进行测试评价；使用新目录，记录代码版本、输入哈希、划分、种子、参数、指标及失败状态。
- 保留同信息的简单对照和失败结果，分别报告合成检查、公开数据回放与独立证据。工程提速或重写不能视为算法创新。

不得覆盖 `reports/v0.7/`、`docs/history/v0.7_import.json` 或修改 `archive/branch_snapshots/` 中的原始树。旧入口可能覆盖结果，运行前参照 [AGENTS.md](AGENTS.md) 使用隔离工作区。不要提交原始大矩阵、私人数据、账户凭据或本地环境文件。

## 提交 Pull Request

1. 从最新 `main` 创建临时工作分支；保留其他贡献者的修改。
2. 遵循现有代码风格，修复根因，并为实际变化的行为补充必要的回归检查。
3. 运行相关检查，检查 `git diff`，更新对应文档；提交说明写明问题及行为变化。
4. 向 `main` 提交 PR，填写动机、改动、实际执行的命令与结果、未执行项及风险。涉及研究时说明证据级别和信息预算。

提交贡献前请确认有权分享这些内容。项目原创软件与对应文档使用 [MIT 许可](LICENSE)；第三方文件必须保留原有许可、署名和来源，范围见 [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md)。不得用本项目的 MIT 许可替换第三方数据或历史材料的原有条款。
