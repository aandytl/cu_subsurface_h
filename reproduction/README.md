# 复现准备

环境独立于系统和共享 conda：

```bash
source .venv/bin/activate
python -m pip install -r reproduction/requirements.txt
```

每次科学计算先报告资源预算并取得用户许可。以下命令尚未执行。

## 第一阶段：100 帧测试集验证

预算：4 CPU 线程、8 GiB 内存、无 GPU、输出 <100 MiB，墙钟时间上限 30 分钟。`timeout` 控制时间；8 GiB 是规划值，实际硬限制须由调度器/运行环境设置。不要把共享环境全部资源视作独占分配。

```bash
timeout --signal=TERM --kill-after=30s 1800s .venv/bin/python reproduction/evaluate_model.py \
  --model DP-GEN_file/DP-GEN_file/model_4.pb \
  --data DP-GEN_test_data/test_data --frames 100 --threads 4 \
  --output reproduction/results/test_100.json
```

指标：能量误差按每帧原子数归一后统计 RMSE；受力按所有 Cartesian 分量统计 RMSE。完整样本选择写入结果。作者误差指标的具体测试范围和权重还需核对，100 帧误差不能直接宣称复现论文误差。

阶段 A 实际结果及首次 ABI 失败记录已写入根目录《复现尝试与资源检查.md》。此 PyPI DeePMD wheel 要求 CPU TensorFlow 2.16.1；以 `requirements.txt` 和 `environment.lock.txt` 为准。

## 第二阶段：全测试集验证（须另行许可）

基于阶段 A 估算约 52 分钟，规划 30–90 分钟；4 CPU、8 GiB 内存预算，无 GPU，结果 <20 MiB。

```bash
/usr/bin/time -v -o reproduction/results/test_all_resources.txt \
  timeout --signal=TERM --kill-after=30s 5400s taskset -c 0-3 \
  .venv/bin/python reproduction/evaluate_model.py \
  --model DP-GEN_file/DP-GEN_file/model_4.pb \
  --data DP-GEN_test_data/test_data --frames 9893 --threads 4 \
  --output reproduction/results/test_all.json \
  > reproduction/results/test_all.log 2>&1
```

## Cu(100) GCMC

正文支持默认 15,000 次提案；代码保留 300 K、H 上限 1 ML、BFGS 最多 400 步。绝对化学势仍须补齐气相参考能与热力学校正，不提供猜测值。入口必须传 `--mu-co`、`--mu-h`、`--dp_model` 和全新 `--output-dir`；`--seed` 默认 20261010。先经许可运行短程基准，再估算完整作业。

代码仍是作者随附接受率实现；其与补充第 5 页的差异及其他尚待核对项详见根目录《复现尝试与资源检查.md》。其他晶面位点、主动学习完整配置和 DFT 作业输入仍需补全。

## 气相参考能近似诊断（须单独许可）

34 个两原子单点，用于筛查势模型在气相构型上的外推风险。两个 delta=0 是显式的“气相减吸附态热校正之差为零”近似，非作者数值；不自动启动 GCMC。即使筛查通过，仍缺气相 DFT/振动校正验证。

```bash
/usr/bin/time -v -o reproduction/results/gas_probe_resources.txt \
  timeout --signal=TERM --kill-after=30s 900s taskset -c 0-3 \
  .venv/bin/python reproduction/probe_gas_references.py \
  --model DP-GEN_file/DP-GEN_file/model_4.pb \
  --delta-h 0 --delta-co 0 --u-rhe -1 \
  --output reproduction/results/gas_probe.json \
  > reproduction/results/gas_probe.log 2>&1
```

## 已完成的弛豫和短程 GCMC（本次已获独立许可）

拟预算：顺序执行，4 CPU、16 GiB 内存、无 GPU、输出预留 100 MiB。气相弛豫 15 分钟超时，GCMC 30 分钟超时。未知热力学校正仍设为零；气相数值收敛不等于 DFT 验证。短程 GCMC 仅 10 次提案，每次/初态最多 20 步 BFGS，不作为论文覆盖度或重构复现。

```bash
/usr/bin/time -v -o reproduction/results/gas_refined_resources.txt \
  timeout --signal=TERM --kill-after=30s 900s taskset -c 0-3 \
  .venv/bin/python reproduction/refine_gas_references.py \
  --probe reproduction/results/gas_probe.json \
  --output reproduction/results/gas_refined.json \
  > reproduction/results/gas_refined.log 2>&1
```

仅在获许可且 gas_refined.json 的 both_converged=true 后运行；入口也会检查此条件：

```bash
/usr/bin/time -v -o reproduction/results/gcmc_pilot_resources.txt \
  timeout --signal=TERM --kill-after=30s 1800s taskset -c 0-3 \
  .venv/bin/python reproduction/run_gcmc_pilot.py \
  --reference reproduction/results/gas_refined.json \
  --output-dir reproduction/results/gcmc_pilot_mu0 \
  > reproduction/results/gcmc_pilot.log 2>&1
```

## 50 提案基线（已独立许可并完成）

预算与近似限制见根目录运行记录。使用独立输出目录；仅在取得本阶段许可后执行以下命令。

```bash
OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 NUMEXPR_NUM_THREADS=4 TF_INTRA_OP_PARALLELISM_THREADS=4 TF_INTER_OP_PARALLELISM_THREADS=1 CUDA_VISIBLE_DEVICES=-1 \
/usr/bin/time -v -o reproduction/results/gcmc_baseline50_resources.txt \
  timeout --signal=TERM --kill-after=30s 5400s taskset -c 0-3 \
  .venv/bin/python GCMC/GCMC/gcmc.py GCMC/GCMC/POSCAR \
  --dp_model DP-GEN_file/DP-GEN_file/model_4.pb \
  --sites-file GCMC/GCMC/adsorption_sites.json \
  --mu-h -2.4292938222817155 --mu-co -14.28483274990397 \
  --steps 50 --temperature 300 --opt-steps 400 --seed 20261010 \
  --output-dir reproduction/results/gcmc_baseline50_mu0 \
  > reproduction/results/gcmc_baseline50.log 2>&1
```
