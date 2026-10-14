# 许可范围、第三方材料与历史归档

根目录 [LICENSE](LICENSE) 中的 MIT 许可适用于本项目原创软件及对应的原创文档，在相应权利人有权授权的范围内生效。它不替换文件已有的第三方许可，不为下载的数据、引用的论文、第三方代码或历史归档授予新的权利。权利归属不明确的材料不能仅凭它出现在仓库中就视为 MIT 授权。

## 范围

| 材料 | 处理方式 |
| --- | --- |
| 当前项目原创代码、测试、工具和对应原创文档 | 适用根目录 MIT 许可，已有更具体许可声明的文件除外。 |
| `data/` 下的数据、下载到本地的数据与上游数据内容 | 不由根目录 MIT 许可重新授权；遵循原数据来源条款，见 [数据署名](data/ATTRIBUTION.md) 和 [固定资产清单](assets/manifest.json)。 |
| `archive/`、`docs/history/` 与 `reports/` 中的既有历史材料 | 不纳入本次 MIT 授权范围；保留既有来源与声明。未带明确许可的历史文件不因归档而获得新的再分发许可。 |
| 引用或复制的第三方代码、图表、论文文字及其他第三方内容 | 继续由原权利人与原条款管理；使用和再分发时单独核对。 |
| Python 及运行、构建依赖 | 使用各上游项目自己的许可；本项目的许可不覆盖它们。 |

历史研究的证据范围另见 [整合政策](research/INTEGRATION_POLICY.md)。本说明不改变归档文件、研究结果或数据的原有许可状态，也不认定所有历史文件已经完成逐项权利审计。

## 辅助信号数据

辅助磷蛋白面板来自 CellNOpt 项目收录的 HepG2CancerRes2011 数据。原有说明保留在 [data/ATTRIBUTION.md](data/ATTRIBUTION.md)，原软件许可原文保留在 [data/CELLNOPT_LICENSE.txt](data/CELLNOPT_LICENSE.txt)。该原文未因新增 MIT 许可而修改。

仓库既有审计未发现该上游目录中独立的数据专用许可；CellNOpt 软件许可不能自动解释为期刊文章或所有数据使用方式的许可。使用者应核查原始来源，保留 Julio Saez-Rodriguez 及合作者的论文署名与 DOI `10.1158/0008-5472.CAN-10-4453`。

## OP3、sci-Plex 与其他研究数据

OP3、sci-Plex 下载入口、固定文件哈希与来源记录见 [assets/manifest.json](assets/manifest.json) 和 [data/ATTRIBUTION.md](data/ATTRIBUTION.md)。应保留 Chem-PerturBridge 分发来源、Open Problems 2023（GEO GSE279945）与 Srivatsan 等 sci-Plex 研究的署名，核查对应数据的原始使用和再分发条款。上游数据仓库的许可元信息不能替代对各原始数据集条款的核查。

历史研究中涉及 L1000、GSE101406 或其他来源时，沿用相应记录中的来源、引用与限制；本文件不新增对这些数据的许可声明。`demo` 生成的明确标注为合成的示例数组不应描述为这些数据集的子集或替代实验观测。

## 维护要求

新增第三方材料时，在同一 PR 中记录来源、版本或提交、原许可和修改情况，并保留原声明。缺少许可或许可范围有疑问时，不将其重新标为 MIT；应先澄清再决定是否分发。

根目录 MIT 文本已按 [SPDX 的 MIT 条目](https://spdx.org/licenses/MIT.html) 与 [Open Source Initiative 的 MIT 文本](https://opensource.org/license/mit) 核对（2026-10-14）。本说明只界定项目授权范围，不在标准 MIT 文本中添加使用领域限制。
