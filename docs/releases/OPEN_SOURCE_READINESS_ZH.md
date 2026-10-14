# 开源工程就绪验收

本轮目标是让使用者从 `main` 找到一个清楚、可安装、可贡献的软件入口，完成离线操作演示，并能区分软件检查与科学证据。当前包版本仍为 `0.1.0`，本文件不宣布新发行、PyPI 发布、默认算法更换或新的科学验证。

## 交付范围

- README 与快速开始聚焦 `main` 和 `perturb-predict`；历史研究保留为可追溯材料。
- 从 `pyproject.toml` 安装当前软件，保留固定 Python 3.13.5 的历史重放环境。
- `demo --out <新目录>` 运行确定性的纯合成流程，输出模型、预测、独立评分及汇总，并显式标注生物学样本数为 0。
- 提供贡献指南、安全报告说明、问题与 PR 模板、原创软件的 MIT 许可及第三方材料排除说明。
- CI 检查构建产物、安装后的命令与测试、Python 3.11–3.14 和 Windows；远端是否成功必须由实际运行记录支持。

## 验收标准与结果记录

本地验收于 2026-10-14 使用 Python 3.13.5、NumPy 2.3.5、h5py 3.15.1；原基线另使用 requirements.txt 的固定依赖。OPENBLAS/OMP/MKL 均为单线程。以下结果来自本轮实际执行，旧的 [历史验收记录](VALIDATION_ZH.md) 不作为本次运行证据。

| 项目 | 验收标准 | 本次结果 |
| --- | --- | --- |
| 源码安装 | 从当前源码安装，`perturb-predict --version` 与包版本一致，`pip check` 通过 | `uv pip install --no-deps -e .`；版本与元数据一致；`uv pip check` 通过。 |
| 分发包 | 构建 sdist 和 wheel，元数据检查通过；安装包包含所需代码、许可和说明，不携带研究数据或历史归档 | `python -m build --no-isolation`（sdist 再构建 wheel）；`python -m twine check --strict dist/*` 通过。检查 wheel 含 CLI、demo 和许可，包内无数据/归档。 |
| 独立安装 | 在仓库外导入已安装包，执行发行测试和离线 demo，避免源码路径掩盖缺文件问题 | 独立虚拟环境安装 wheel，仓库外导入 site-packages；发行测试 78/78 通过，`perturb-predict demo` 完成。 |
| 合成演示 | 无网络也可完成；训练/查询分离；先冻结预测再评分；汇总明确标注合成来源和非独立证据 | 5 项 demo 回归通过；固定种子、完整流水线、查询标签隔离、拒绝覆盖及失败标记均已检查。 |
| 运行产物 | 拒绝覆盖已存在输出；保留完成或可记录的失败状态；产物哈希与实际文件一致 | 继承的 5 项 JSON 原子写入测试、新增 17 项归档/状态/下载边界及 4 项通用下载竞争测试通过。 |
| 回归检查 | 发行测试、相关原基线测试和桥接设计检查按变更范围执行；失败或跳过逐项说明 | `unittest discover -s tests_release -v` 78/78；`-s tests -v` 54/54；`-s studies/bridge_design_v14/tests -v` 18/18；`tools.smoke` 完成。 |
| 兼容性 CI | Linux Python 3.11、3.12、3.13、3.14，以及 Windows Python 3.13 的安装后测试有真实运行结果 | 工作流已配置；本地只验证 Linux Python 3.13.5，远端矩阵结果等待本次 PR 的 Actions。 |
| 文档与模板 | 文档中的命令、路径和相对链接与当前代码一致；issue 模板 YAML 可解析 | 47 个相对文件链接、4 个 issue YAML 与本地标题锚点检查通过；`git diff --check` 通过。 |
| 许可与归档 | MIT 文本与官方标准一致；第三方 notice 未被改写；原归档树和科研状态不被此次工程改动覆盖 | 标准 MIT 文本已核对，原 data/CELLNOPT_LICENSE.txt 未改；36 个旧分支 SHA 与 29 份精确快照已记录，合入前以 tools.verify_history 核验祖先关系。 |

本地运行无法替代未执行的操作系统或 Python 版本检查。网络受限、依赖不可用或 CI 未执行时，保留明确的阻塞记录，不以旧成功记录或配置存在来推断通过。

## 固定数据回放

从仓库外运行已安装 wheel 的 `tests_release/acceptance.py`，两个目标的两次完整重放与独立模型推理逐元素一致；公开查询子集被拒绝且没有读取真值。既有参考 MSE `0.1980665327273829` 未改变。这是已暴露公开数据的工程回归，不是新生物学结果。本地记录目录为 `runs/installed-op3-20261014/`，辅助集成记录为 `runs/open-source-smoke-20261014/`，二者不随源码包分发。

## 科学边界

离线演示检查操作和实现，公开 OP3 重放检查既有工程流程；二者都不能证明独立制备上的生物学效果。当前唯一科研主线及其新增证据要求保留在 [研究路线](../../research/ACTIVE_ROUTE_ZH.md)。本轮不关闭独立证据关卡，不修改历史分数，不将其他研究分支的结果合并为默认模型的证据。

## 许可审阅范围

MIT 许可用于项目原创软件及对应文档；第三方数据、第三方材料与历史归档排除在本次授权范围之外，详见 [THIRD_PARTY_NOTICES.md](../../THIRD_PARTY_NOTICES.md)。新增许可证的合并审阅需确认项目有权授予这些原创材料的许可；本轮不声称已完成所有历史文件或上游数据的权利审计。

标准文本核对来源为 [SPDX MIT](https://spdx.org/licenses/MIT.html) 和 [OSI MIT](https://opensource.org/license/mit)。许可范围说明与软件的研究用途说明分别记录；标准 MIT 正文未添加用途限制。
