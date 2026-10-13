# 铜表面重构与次表层氢：项目介绍

> 本文的文件状态说明基于原始交付材料的静态检查。后续已修复 GCMC 启动入口并验证交付势模型；当前运行状态、资源实测及与原文的对比见 [`复现尝试与资源检查.md`](复现尝试与资源检查.md)。

## 1. 项目定位与研究问题

本目录对应论文 **“Subsurface hydrogen as a hidden driver of copper surface reconstruction in CO₂ electroreduction”**（Zhang 等，*National Science Review* **13**, nwag128, 2026；[DOI](https://doi.org/10.1093/nsr/nwag128)；本地论文见 [`nwag128.pdf`](nwag128.pdf)）。研究对象是 CO₂ 电还原反应（CO₂RR）过程中铜催化剂的表面结构演化。

核心问题是：在约 **−1 V vs. RHE** 的典型工作电位下，表面吸附的 CO（`*CO`）覆盖度有限，而氢（`*H`）较多；除了已知的 `*CO`／表面氢（`*Hsur`）协同作用，**进入铜次表层的氢（`*Hsub`）能否独立驱动铜原子迁移与表面粗化？** 论文进一步比较 Cu(100)、Cu(111)、Cu(211)、Cu(711) 的晶面差异，并讨论这种结构变化对 C₂⁺ 产物倾向和催化剂稳定性的影响，以及 Zn、Al、Ga 合金化能否抑制次表层氢。

这里的“重构”主要指 Cu(100) 样表面在高氢覆盖度下形成**局部、无序且动态变化的三重空位密排区域**；论文借用“氢诱导六角重构”描述其几何趋势，不表示获得了完整、严格周期性的六角晶相。

## 2. 计算路线：模型、算法与输入输出

| 环节 | 功能及主要输入 | 主要输出／在本项目中的作用 |
| --- | --- | --- |
| **PH-SA 位点采样** | 从 Cu 晶面、团簇和吸附体系的原子坐标及晶胞出发，识别表面与体相可用位点；论文使用研究组的 PH-SA 方法。 | 用于构建 `*CO`、`*Hsur`、`*Hsub` 的初始构型；本目录的 [`adsorption_sites.json`](GCMC/GCMC/adsorption_sites.json) 是针对随附 Cu(100) `POSCAR` 的预计算位点，含 H 位点坐标及 CO 位点坐标和法向量。**未附 PH-SA 的生成程序**。 |
| **VASP 密度泛函理论（DFT）／AIMD** | 输入 `POSCAR`、`INCAR`，实际计算还需赝势 `POTCAR`、k 点设置等；电子结构方法为 revPBE、PAW、400 eV 截断、DFT-D3。短时 NVT AIMD 生成候选构型，DFT 给出标注。 | 原子结构、总能和受力；论文初始数据 **58,740** 个结构。关键反应能垒也由 DFT 计算。目录中的 [`DFT calculation file`](DFT%20calculation%20file/DFT%20calculation%20file) 主要保存结构与部分 `INCAR`，并非完整 VASP 作业或能量输出。 |
| **DP-GEN 主动学习与 DeePMD 深度势（MLIP）** | 以 Cu–C–O–H 的结构、晶胞、DFT 总能和逐原子受力为训练数据；每轮训练四个势模型，以模型间最大受力偏差筛选需重新做 DFT 标注的构型。论文阈值为 **0.10–0.35 eV/Å**；候选构型由 LAMMPS + DeePMD 分子动力学探索。 | 最终 **98,933** 个独立结构；论文报告势模型能量误差 **5.42 meV/atom**、受力误差 **93.5 meV/Å**。目录提供训练／验证／测试数据及一个冻结势文件 [`model_4.pb`](DP-GEN_file/DP-GEN_file/model_4.pb)，未提供完整 DP-GEN 迭代配置、四个模型及训练日志。 |
| **DPMD（深度势分子动力学）** | 深度势、晶面或约 2.5 nm Cu 纳米团簇和高覆盖 `*Hsur` 初态。 | 在比直接 AIMD 更大的模型上观察结构演化；论文用团簇模拟确认 (100) 面的类似重构。本目录有 [`wulff_cluster/POSCAR`](DFT%20calculation%20file/DFT%20calculation%20file/wulff_cluster/POSCAR) 等结构，但无完整 DPMD 轨迹或运行配置。 |
| **巨正则蒙特卡洛（GCMC）+ 势模型局部弛豫** | Cu 晶面 `POSCAR`、位点 JSON、`*CO` 和 `*H` 的化学势、温度、深度势；尝试 CO/H 插入、删除和位点迁移，每次候选结构用 ASE BFGS 与深度势弛豫，再按能量与化学势变化做 Metropolis 接受判断。 | 吸附物覆盖度、末态和可供局部环境统计的轨迹；论文在四个晶面各做 **15,000 步**，并做延长及独立轨迹检查。本目录仅有 [`gcmc.py`](GCMC/GCMC/gcmc.py)、Cu(100) 初始 [`POSCAR`](GCMC/GCMC/POSCAR) 和位点 JSON；其他晶面现成的终态结构位于 DFT 目录。 |
| **CI-NEB、振动检验、GC-DFT、LOBSTER/COHP** | `is`（初态）、`ts`（过渡态）、`fs`（终态）构型和相应 DFT 参数；GC-DFT 调整电子数／电势；LOBSTER 对 DFT 电子结构做 Cu–Cu 键的 COHP/ICOHP 分析。 | 氢向次表层迁移、Cu 加合原子形成及 `*CO → *CHO` 质子化的势垒；确认过渡态只有一个虚频，评估电势敏感性及 Cu–Cu 键强度。目录保存大量 `is/ts/fs` 的 `POSCAR`/`CONTCAR`，但缺少 NEB 图像、频率、LOBSTER 输出及完整能量表。 |
| **局部环境、GCN 与动力学解释** | 从 GCMC 构型统计覆盖度、Cu 配位、广义配位数（GCN）、原子位移；将代表位点送入 DFT 势垒计算。 | 联系“氢进入次表层 → Cu 键弱化与低配位位点增多 → Cu 加合原子更易形成、C–C 偶联可能更有利”。论文图 1 提及动力学建模，但**本目录没有独立的微观动力学求解脚本或产物选择性原始计算表**；这里的选择性判断主要由结构和反应势垒趋势支持。 |

论文方法部分说明：GCMC 大体系主要依靠 MLIP 计算能量；氢迁移、铜加合原子形成及 CO 质子化等**关键势垒仍由 DFT** 计算。电化学条件在 GCMC 中通过吸附物化学势间接表示；提供的 Python 脚本没有显式水、离子或电极电势求解器。

## 3. 论文主要结果与机制分析

1. **晶面决定氢能否进入次表层。** 在 −1 V vs. RHE 的论文 GCMC 结果中，Cu(100) 的 `*CO` 约 **0.29 ML**、`*Hsub` 约 **0.18 ML**；具有 (100) 台面的 Cu(711) 分别约 **0.29 ML** 和 **0.26 ML**。Cu(211) 的 `*CO` 约 **0.21 ML**、`*Hsub` 仅 **0.02 ML**；Cu(111) 基本未见 `*CO` 吸附或氢渗入。这里的 ML 是相对于表面位点数的覆盖度，不是文件中的原子比例。
2. **先有局部重构，再容易渗氢。** 高 `*Hsur` 覆盖使 Cu(100) 部分四重空位区域向三重空位密排区域重排；剩余四重区域被拉开。论文图 3 统计的剩余四重位平均面积从理想表面的 **6.53 Å²** 增至 **6.96 Å²**，局部 Cu–Cu 距离可达 **2.80 Å**（理想值约 **2.53 Å**）。理想 Cu(100) 上氢进入八面体／四面体次表层位的势垒分别是 **1.06／0.72 eV**；在重构的局部四面体通道中降至 **0.11 eV**。这支持论文强调的几何开口机制。−0.2 V vs. RHE 的对比模拟以 `*CO` 为主（约 0.70 ML），未出现同类明显重构。
3. **`*Hsub` 可独立促进铜加合原子形成。** 无邻近 `*CO` 时，仅邻近 2–3 个 `*Hsur` 的 Cu 加合原子形成势垒仍约 **1.45 eV**；在 `2*Hsur + 1*Hsub` 和 `2*Hsur + 2*Hsub` 局部环境中分别降至 **0.97 eV**、**0.61 eV**。若 Cu 同时邻近一个 `*CO` 和一个 H，H 位于表面时势垒约 **1.0 eV**，位于次表层时约 **0.73 eV**。LOBSTER 的 ICOHP 分析将这一趋势与 Cu–Cu 键弱化联系起来。
4. **活性与稳定性相互牵连。** 重构增加低 GCN Cu 位点；论文的 DFT 势垒趋势表明这些位点既可能有利于 C–C 偶联，又较不易发生 `*CO → *CHO` 质子化，因此可能提高 C₂⁺ 倾向。但同一重构也促进表面 Cu 迁移，可能损害长期稳定性。这里是**基于原子模拟和势垒的机制推断**，不能直接读作已测得的 C₂⁺ 法拉第效率或寿命提升。
5. **合金化提供抑制途径。** 对 Cu(100) 掺入 Zn、Al、Ga 的 DFT 计算显示 H 吸附变弱、向次表层迁移势垒上升；1 ML 表面 H 和 50% 表面掺杂的 AIMD 模型中，表面 Cu 保留更高配位，六角样重构受到抑制。论文未系统优化掺杂浓度，也未以本目录文件给出可直接复算的器件级性能数据。

论文还指出 pH、电位范围、界面水、阳离子电场及其他吸附物可能改变上述趋势；`*Hsub` 只是表面演化的一个关键驱动因素。四晶面的 15,000 步 GCMC 及额外 5,000 步检查属于**采样步数**，不等同于真实反应时间。

## 4. 目录中的数据及可用格式

| 路径／规模（静态清点） | 物理意义与程序可读形式 |
| --- | --- |
| [`DP-GEN_training_data/training_data`](DP-GEN_training_data/training_data)：**720** 个按化学计量命名的体系、**79,147** 帧，约 **244 MB** | DeePMD 训练集。每个 `Cu…C…O…H…` 子目录中的 `type_map.raw` 固定为 `Cu C O H`；`type.raw` 为逐原子物种索引；`set.000/coord.npy` 是逐帧笛卡尔坐标（形状 `Nframe × 3Natom`），`box.npy` 为晶胞（`Nframe × 9`），`energy.npy` 为逐帧总能，`force.npy` 为逐帧逐原子力（`Nframe × 3Natom`）。坐标／晶胞通常按 Å、能量按 eV、力按 eV/Å 理解，并应与 DFT 标注设置保持一致。 |
| [`DP-GEN_validation_data/validation_data`](DP-GEN_validation_data/validation_data)：**699** 个体系、**9,893** 帧，约 **31 MB** | 同格式验证集，供训练期模型选择和误差监测；不是新的化学体系种类定义。 |
| [`DP-GEN_test_data/test_data`](DP-GEN_test_data/test_data)：**702** 个体系、**9,893** 帧，约 **31 MB** | 同格式测试集，供训练后能量、受力泛化误差评估。三个分区合计恰为论文所述 **98,933 帧**；单帧原子数约 **24–139**，各计量组可有不同帧数。`.raw` 是文本，`.npy` 为 NumPy 二进制数组；上述数据可由 DeePMD 标准数据集读取器或 NumPy 读取。 |
| [`model_4.pb`](DP-GEN_file/DP-GEN_file/model_4.pb)：约 **16.6 MiB** | Cu–C–O–H 冻结 Deep Potential 模型，供 DeePMD 计算器从原子类型、坐标、晶胞预测总能和原子力；这是**一个**交付模型，不能单凭它复现 DP-GEN 四模型偏差筛选。 [`DP-GEN_file/INCAR`](DP-GEN_file/DP-GEN_file/INCAR) 是 AIMD 标注示例参数。 |
| [`GCMC/GCMC`](GCMC/GCMC)：`POSCAR`、JSON、Python 脚本、说明 | `POSCAR` 是 **432 个 Cu 原子**的 3 层 12×12 Cu(100) 表面模型；JSON 保存 **576 个 H 位点、576 个 CO 位点及对应法向量**，坐标为三维数值列表。脚本用 ASE 读取 VASP 结构，用 JSON 选位点，以 DeepMD + ASE BFGS 求能量，写出 `gcmc_main.traj`、逐步 `trajectories/*.traj`、`optimized_structures/*.vasp` 和当前 `CONTCAR`。位点文件依赖该结构的晶胞与坐标，不能直接套用于其他晶面。 |
| [`DFT calculation file/DFT calculation file`](DFT%20calculation%20file/DFT%20calculation%20file)：**208** 个文件、约 **1.4 MB** | VASP `POSCAR`／`CONTCAR` 为晶胞、元素、原子数、坐标及部分约束的文本结构；`INCAR` 为计算参数。子目录按科学问题划分：Cu adatom formation（103 文件）、Hsub incorporation（52）、CO protonation to CHO（25）、alloy（15）、GC-DFT calculations（7）、GCMC simulation results（4）及 wulff_cluster（2）。`is/ts/fs` 分别表示反应路径初态／过渡态／终态。GCMC 结果结构不是完整覆盖度时间序列。 |

数据集计量目录名（例如 `Cu24C4O4H8`）表示该组原子组成，**不是该组样本帧数**。`C` 和 `O` 数目相同时通常对应完整 CO 分子；H 数目表示氢原子总数，单靠目录名无法区分 `*Hsur` 与 `*Hsub`，须结合坐标、表面层位置和局部环境判定。对于 DFT 结构，`CONTCAR` 是已保存的优化或动力学构型，不能仅凭该文件还原总能、收敛性或势垒。

## 5. 目前可进行模拟的指令集合

以下命令在 **Linux 项目根目录**执行，各任务独立选择。新计算仍须先估算资源并取得许可；命令示例不代表已获运行许可。已完成结果及与原文的差异见 [复现尝试与资源检查.md](复现尝试与资源检查.md)。目前可运行的是单模型验证、PH-SA 位点生成、近似气相参考和 Cu(100) GCMC；DFT/NEB、合金 AIMD 和完整 DP-GEN 的输入及环境仍待补齐。

### 5.1 环境、输入与通用运行方式

已有 `.venv` 时直接使用；新环境使用 Python 3.9。当前机器可按下列路径创建，其他机器将解释器路径替换为自己的 Python 3.9：

```bash
/opt/miniconda3/envs/ace_gcn_39/bin/python -m venv .venv
.venv/bin/python -m pip install -r reproduction/requirements.txt
.venv/bin/python -m pip check
```

Git 仓库不含训练/验证/测试数据集。测试命令需先从 [论文 Zenodo 存档](https://zenodo.org/records/18375029) 获取并解压测试集，确认存在 `DP-GEN_test_data/test_data/*/set.000/coord.npy`。GCMC 使用随附 `DP-GEN_file/DP-GEN_file/model_4.pb`、`GCMC/GCMC/POSCAR` 和 `adsorption_sites.json`；PH-SA 使用 `reproduction/PH-SA_upstream/` 的固定版本源码。

每次操作创建一个新目录，以下变量和函数在同一个 Bash 终端中有效：

```bash
mkdir -p reproduction/results
export CU_RUN_DIR=$(mktemp -d reproduction/results/manual_XXXXXX)
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
export NUMEXPR_NUM_THREADS=4 TF_INTRA_OP_PARALLELISM_THREADS=4
export TF_INTER_OP_PARALLELISM_THREADS=1 CUDA_VISIBLE_DEVICES=-1

# 参数：日志名称、超时秒数、地址空间上限 GiB、Python 脚本及其参数。
run_cpu() {
  local cu_name="$1" cu_seconds="$2" cu_memory="$3"
  shift 3
  if [[ -e "$CU_RUN_DIR/$cu_name.log" || -e "$CU_RUN_DIR/$cu_name.resources.txt" ]]; then
    echo "日志已存在，请更换名称或新建 CU_RUN_DIR。" >&2
    return 1
  fi
  /usr/bin/time -v -o "$CU_RUN_DIR/$cu_name.resources.txt" \
    timeout --signal=TERM --kill-after=30s "$cu_seconds" \
    prlimit --as="$((cu_memory * 1024 * 1024 * 1024))" \
    taskset -c 0-3 .venv/bin/python "$@" \
    > "$CU_RUN_DIR/$cu_name.log" 2>&1
}
```

`taskset` 要求 CPU 0–3 在当前任务允许范围内；受调度器分配时按实际 CPU 编号调整。`prlimit` 限制地址空间，包含虚拟映射，比单纯 RSS 预算更保守。输出目录总量需另行检查。函数前台执行，资源记录中的 `Exit status: 0` 表示成功，`124` 通常表示超时。需断开终端时可使用 `tmux`；函数不能直接作为 `nohup` 的可执行程序。

### 5.2 交付模型测试集验证

100 帧用于兼容性检查；9,893 帧用于复算论文报告误差。预算 4 CPU、8 GiB、无 GPU；分别预留 100 MiB/20 MiB 输出，超时 30/90 分钟。全测试此前实测约 13 分 35 秒，共享环境耗时会变化。

```bash
run_cpu test100 1800 8 reproduction/evaluate_model.py \
  --model DP-GEN_file/DP-GEN_file/model_4.pb \
  --data DP-GEN_test_data/test_data --frames 100 --threads 4 \
  --output "$CU_RUN_DIR/test100.json"
```

```bash
run_cpu test_all 5400 8 reproduction/evaluate_model.py \
  --model DP-GEN_file/DP-GEN_file/model_4.pb \
  --data DP-GEN_test_data/test_data --frames 9893 --threads 4 \
  --output "$CU_RUN_DIR/test_all.json"
```

### 5.3 Cu(100) PH-SA 表面位点重建

使用 `k=1.0`、H 半径和 CO 的 C 锚点，仅生成表面位点；已验证各 576 点及 CO 法向量与交付输入匹配。参数来自逆向分析，尚不能唯一证明作者实际调用方式。预算 4 CPU、16 GiB、无 GPU、100 MiB 输出，30 分钟超时；此前实测 50.14 秒。候选不会自动替换 GCMC 输入。

```bash
run_cpu phsa_surface 1800 16 reproduction/compare_phsa_sites.py \
  --k 1.0 --co-anchor C --surface-only \
  --output-dir "$CU_RUN_DIR/phsa_surface"
```

输出 `comparison.json`、`raw_sites.json`、`candidate_sites.json`。省略 `--surface-only` 会增加内部位点计算，需重新评估预算。

### 5.4 MLIP 气相诊断与键长弛豫

这是近似参考方案：`delta_H=delta_CO=0` 表示把气相减吸附态的热校正差设为零，不能当作作者数值或 DFT 验证。诊断共 34 个两原子单点；弛豫每种最多 50 次搜索加末态求值。每项预算 4 CPU、8 GiB、无 GPU、10 MiB 输出、15 分钟超时。

```bash
run_cpu gas_probe 900 8 reproduction/probe_gas_references.py \
  --model DP-GEN_file/DP-GEN_file/model_4.pb \
  --delta-h 0 --delta-co 0 --u-rhe -1 \
  --output "$CU_RUN_DIR/gas_probe.json"
```

仅在诊断成功、检查其筛查标记后执行弛豫；不能在前一步失败时继续：

```bash
run_cpu gas_refined 900 8 reproduction/refine_gas_references.py \
  --probe "$CU_RUN_DIR/gas_probe.json" \
  --output "$CU_RUN_DIR/gas_refined.json"
```

### 5.5 Cu(100) 近似 GCMC

以下沿用已弛豫的 MLIP 气相能、零热校正及 −1 V RHE，取 `mu_H=−2.4292938222817155 eV`、`mu_CO=−14.28483274990397 eV`；**不是作者精确参数**。继续使用交付位点和随附接受规则，不自动修正补充公式中的体积/热波长因子。

10 次提案流程测试（需 5.4 的 `gas_refined.json`，入口检查两种气体是否收敛）：预算 4 CPU、16 GiB、无 GPU、100 MiB 输出、30 分钟超时，每候选/初态最多 20 步 BFGS。

```bash
run_cpu gcmc_pilot 1800 16 reproduction/run_gcmc_pilot.py \
  --reference "$CU_RUN_DIR/gas_refined.json" \
  --output-dir "$CU_RUN_DIR/gcmc_pilot"
```

200 次提案独立短链：预算 4 CPU、16 GiB、无 GPU、2 GiB 输出，预计 60–150 分钟、180 分钟超时，每候选/初态最多 400 步 BFGS。两条已完成链实测约 65–84 分钟，均未达到采样收敛。

```bash
run_cpu gcmc200 10800 16 GCMC/GCMC/gcmc.py GCMC/GCMC/POSCAR \
  --dp_model DP-GEN_file/DP-GEN_file/model_4.pb \
  --sites-file GCMC/GCMC/adsorption_sites.json \
  --mu-h -2.4292938222817155 --mu-co -14.28483274990397 \
  --steps 200 --temperature 300 --opt-steps 400 --seed 20261010 \
  --output-dir "$CU_RUN_DIR/gcmc200"
```

独立种子对照可将 `--seed` 改为 `20261013`，同时更换日志名称及输出目录。新计算分别申请预算，不将相同接受率解释为覆盖度收敛。

### 5.6 检查点续跑与一致性验证

连续 4 次对比分段 2+2 次（总共 8 次提案、两次初态优化）：预算 4 CPU、16 GiB、无 GPU、100 MiB 输出、20 分钟超时，实测约 6 分 38 秒。

```bash
run_cpu resume_check 1200 16 reproduction/verify_gcmc_resume.py \
  --output-root "$CU_RUN_DIR/resume_check"
```

已完成的 `gcmc200` 可在新的预算许可下续至总计 400 次。下面新增最多 200 次，示例超时 180 分钟、16 GiB 地址空间保护，输出预留 2 GiB；实际预算需结合当前覆盖度确认。`--steps 400` 是**总目标**，不是再做 400 次；输入、化学势、优化设置及种子须与原运行一致。只加载可信本地 `checkpoint.pkl`。跨终端时先将 `CU_RUN_DIR` 设为原运行目录。

```bash
run_cpu gcmc_resume400 10800 16 GCMC/GCMC/gcmc.py GCMC/GCMC/POSCAR \
  --dp_model DP-GEN_file/DP-GEN_file/model_4.pb \
  --sites-file GCMC/GCMC/adsorption_sites.json \
  --mu-h -2.4292938222817155 --mu-co -14.28483274990397 \
  --steps 400 --temperature 300 --opt-steps 400 --seed 20261010 \
  --resume --output-dir "$CU_RUN_DIR/gcmc200"
```

新脚本每个完成提案保存检查点；旧的 50 提案目录没有检查点，不能精确续跑。任意崩溃时刻可能造成日志/轨迹与检查点不一致，须核对后恢复，不能直接删除记录绕过检查。

### 5.7 查看与分析已有结果

以下只读取已有记录和结构，不产生新提案或势能计算：

```bash
tail -n 20 "$CU_RUN_DIR/gcmc200.log"
cat "$CU_RUN_DIR/gcmc200.resources.txt"
du -sh "$CU_RUN_DIR/gcmc200"

.venv/bin/python reproduction/summarize_gcmc.py "$CU_RUN_DIR/gcmc200" \
  --output "$CU_RUN_DIR/gcmc_summary.json"
.venv/bin/python reproduction/analyze_gcmc_final.py "$CU_RUN_DIR/gcmc200" \
  --output "$CU_RUN_DIR/gcmc_geometry.json"
```

对比项目内已完成的两条独立轨迹并生成图：

```bash
MPLCONFIGDIR=/tmp/cu_matplotlib_cache .venv/bin/python \
  reproduction/compare_gcmc_seeds.py \
  --first reproduction/results/gcmc_baseline200_mu0 \
  --second reproduction/results/gcmc_seed20261013_mu0 \
  --output-dir "$CU_RUN_DIR/seed_comparison"
```

结构分析的 Cu–Cu 截断 3.0 Å、GCN 分母 12、H 深度±0.2 Å 模糊带均为本次预设。总 H 数量不能当作 Hsub；这些统计不替代作者判据、采样收敛及 DFT 势垒验证。
