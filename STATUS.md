# 项目状态

更新：2026-10-14。`main` 是唯一维护版本，正常使用从 [README](README.md) 开始。

## 当前软件

- 安装包：`perturb-predict`，Python API：`perturb_predict.TransferModel`。
- 任务：单源、少样本目标校准的化学扰动响应预测；查询药物必须在源背景有实测响应。
- 使用流程：离线 `demo` → 自有数据 `train` / `predict` / `evaluate`；固定 OP3 回放保留作工程验收。
- 本轮：补齐许可、贡献指南、包元数据、安装检查与离线示例；加强文件持久化、损坏输入和失败记录处理。默认数值算法未更换。
- 验收范围及实际运行记录：[OPEN_SOURCE_READINESS_ZH.md](docs/releases/OPEN_SOURCE_READINESS_ZH.md)。这不是新性能发行或独立生物学验证。

## 历史与研究

旧分支通过 [归档清单](archive/branch-heads-20261014.json) 和主线提交历史保留；研究编号不作为另一套软件版本。以前的整合状态原文保存在 [历史状态](docs/history/STATUS_before_single_main.md)。

研究仍受 [ACTIVE_ROUTE_ZH.md](research/ACTIVE_ROUTE_ZH.md) 与 [INTEGRATION_POLICY.md](research/INTEGRATION_POLICY.md) 约束。独立制备桥接与新证据仍未验证；原始科学问题不因工程完善而解决。#17 继续作为独立证据关卡。没有自动科研循环、付费计算、作者邮件或性能发布。
