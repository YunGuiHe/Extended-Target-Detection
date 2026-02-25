# Extended Target Detection (扩展目标稳定检测)

本项目旨在实现并复现“面向扩展目标检测的概率鲁棒雷达波形设计（Probabilistically Robust Radar Waveform Design for Extended Target Detection）”一文中的核心算法与实验。

## 工程内容说明

1. **核心论文与相关文档**
   - `Probabilistically_Robust_Radar_Waveform_Design_for_Extended_Target_Detection.pdf`: 论文原文，本工程复现的基础理论。
   - `论文解读.md`: 针对该论文方法、公式推导以及算法流程的具体梳理和解读。

2. **算法代码库 (`code/`)**
   - **`algorithms/`**：包含基础的波形设计算法以及论文提出的不同条件下的凸优化与求解方法。
   - **`experiments/`**：包含论文中用来评估收敛性、输出信杂噪比(SINR)、合成损失以及功率谱密度(PSD)对比的一系列验证实验代码。
   - **`utils/`**：包含用来计算协方差矩阵、波形综合分析的通用工具函数。
   - **执行主入口**：可通过运行 `code/` 根目录下的各类 `.m` 脚本（如 `main_run.m`，或单独的具体实验脚本）来进行效果仿真与波形生成。
   - **结果分析**：对应的评估指标或仿真截图会被生成并保存在本目录。

## 环境
建议在 MATLAB 中运行该库中的代码。涉及 SDP 求解的部分可能需要安装如 CVX 等相应的优化计算工具箱。
