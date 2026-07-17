# rho-prime UPC 四 pion 事件产生程序：学生交接说明

这份代码用于计算重离子 UPC 中

`rho'(1450) -> rho0 sigma -> pi+ pi- pi+ pi-`

的产生截面网格，并从网格中抽样生成四 pion 事件。你拿到 zip 后，先不要改代码，按下面顺序跑通最小例子。

## 1. 目录里有什么

- `include/RhoPrime/`：C++ 头文件。
- `src/`：主要 C++ 源码。
- `CMakeLists.txt`：CMake 编译配置。
- `UPC_Probabilities_*.root`：UPC/Glauber 概率表。生成截面网格时会自动读取。
- `calc_UPC_EMD_Glauber_unified.C`：重新生成 UPC 概率表的 ROOT 宏。一般不用重跑。
- `analyze_bose.cxx`、`analyze_nonbose.cxx`：分析生成的事件 ROOT 文件。
- `compare_bose_modulations.cxx`：比较 Bose / non-Bose 的角分布调制。
- `review/TECHNICAL_REVIEW_NOTES.md`：以前的技术检查记录，比较长，遇到高级问题再看。
- `review/scripts/compare_event_samples.C`：比较加权样本和 unit-weight 样本的 ROOT 宏。
- `reference_plots/`：旧的参考图，只用于看大概趋势。
- `note.pdf`：已经写好的物理说明，建议先读这个理解模型和物理动机。
- `note.tex`：`note.pdf` 的 LaTeX 源文件，以后需要修改 note 时用。

没有打包 `build/` 目录，因为里面有几十 GB 的历史输出文件。你需要在自己的机器上重新编译。

## 2. 需要先安装什么

需要 Linux/WSL 环境，并且能使用：

- CMake，建议版本 `>= 3.16`
- C++17 编译器，例如 `g++`
- ROOT，已测试版本为 ROOT `6.30`

检查命令：

```bash
root-config --version
cmake --version
g++ --version
```

如果 `root-config` 找不到，说明 ROOT 环境没有加载。常见做法是先执行类似：

```bash
source /path/to/root/bin/thisroot.sh
```

具体路径取决于你自己的 ROOT 安装位置。

如果在 lxplus/CVMFS 上使用本目录里的 `build_cmssw_1511` 可执行文件，运行前必须加载同一个 CMSSW runtime：

```bash
cd /eos/cms/store/group/phys_heavyions/jianjie/CMSSW_15_1_1/src
cmsenv
cd ../../MC_UPC_RhoPrime4Pi_new
```

也可以在 `MC_UPC_RhoPrime4Pi_new` 目录里直接使用包装脚本：

```bash
./run_cmssw_1511.sh ./build_cmssw_1511/generate_bose --help
```

否则运行 `build_cmssw_1511/generate_bose` 时可能会报 `libtbb.so.12: cannot open shared object file` 或 `GLIBCXX_3.4.30 not found`。

## 3. 第一次编译

进入解压后的目录，例如：

```bash
cd rhoprime_student_package_20260708
```

创建 build 目录并编译：

```bash
cmake -S . -B build
cmake --build build -j4
```

编译成功后，应该看到这些可执行程序：

```bash
ls build/generate_grid_opt build/generate_bose build/analyze_bose build/analyze_nonbose build/compare_bose_modulations
```

## 4. 最小跑通例子

第一步：生成一个很小的 PbPb 5.36 TeV、NoTag 触发截面网格。

```bash
./build/generate_grid_opt PbPb 5360 NoTag build/grid_PbPb5360_NoTag_smoke.root 40 100 1 4 1 20 400 1
```

参数含义按顺序是：

- `PbPb`：碰撞系统，也可用 `AuAu`
- `5360`：sqrt(s_NN)，单位 GeV
- `NoTag`：触发标签，也可用 `XnXn`
- `build/grid_PbPb5360_NoTag_smoke.root`：输出网格文件
- `40`：质量 M 的 bin 数，测试时用小值，正式计算要加大
- `100`：横动量网格每边 bin 数，测试时用小值
- `1`：小 b 区域积分步数，测试时用小值
- `4`：大 b 区域步长，单位 fm
- `1`：快速度范围为 `-1 < y < 1`
- `20`：y 的 bin 数
- `400`：全局横向空间盒子大小，单位 fm
- `1`：要求必须找到匹配的 UPC 概率 ROOT 文件；找不到就直接报错，避免误用近似结果

成功时终端会打印一行 `RESULT ... output=...`。

第二步：从已有生产网格生成 5,000,000 个 Bose 对称事件。测试时可以把事件数改成 `1000`。

```bash
./run_cmssw_1511.sh ./build_cmssw_1511/generate_bose grids/grid_PbPb5360_NoTag_prod.root build/events_bose_5000000.root 5000000 12345 500 1
```

参数含义：

- `grids/grid_PbPb5360_NoTag_prod.root`：输入网格
- `build/events_bose_5000000.root`：输出事件 ROOT 文件
- `5000000`：事件数；测试时可改成 `1000`
- `12345`：随机数种子
- `500`：每个质量点估计衰变权重归一化的 trial 数；测试用 500，正式可用 1000 或更高
- `1`：打开 Bose symmetrization；如果填 `0` 就是 non-Bose 对照样本

成功时会看到：

```text
Generated Bose-symmetrized rho-prime events: 5000000/5000000
Output ROOT: build/events_bose_5000000.root
```

第三步：分析事件文件，输出直方图。

```bash
./run_cmssw_1511.sh ./build_cmssw_1511/analyze_bose build/events_bose_5000000.root build/analyze_bose_5000000.root 5000000
```

输出文件 `build/analyze_bose_5000000.root` 里包含质量、快速度、pT、pair mass、角分布等直方图。

## 5. 生成 non-Bose 对照样本

non-Bose 样本只需要把 `generate_bose` 的第 6 个可选参数设为 `0`：

```bash
./run_cmssw_1511.sh ./build_cmssw_1511/generate_bose grids/grid_PbPb5360_NoTag_prod.root result/events_nonbose_5000000.root 5000000 12345 500 0
./run_cmssw_1511.sh ./build_cmssw_1511/analyze_nonbose result/events_nonbose_5000000.root build_cmssw_1511/analyze_nonbose_5000000.root 5000000
```

```bash
./run_cmssw_1511.sh ./build_cmssw_1511/generate_bose grids/grid_PbPb5360_NoTag_prod.root result/events_bose_5000000.root 5000000 12345 500 1
./run_cmssw_1511.sh ./build_cmssw_1511/analyze_nonbose result/events_bose_5000000.root build_cmssw_1511/analyze_bose_5000000.root 5000000
```

比较 Bose 和 non-Bose：

```bash
./run_cmssw_1511.sh ./build_cmssw_1511/compare_bose_modulations build/events_bose_5000000.root build/events_nonbose_5000000.root
```

## 6. 事件生成模式怎么选

默认模式是加权事件，最稳，建议先用这个：

```bash
./build/generate_bose input_grid.root output_weighted.root 20000 6161 1000 1
```

如果需要 `event_weight=1` 的 unit-weight 事件，可以用后面的参数打开：

```bash
./build/generate_bose input_grid.root output_unit.root 20000 6161 1000 1 1 50000 1.0 3 1
```

这里额外参数含义是：

- `1`：打开 unit-weight 输出
- `50000`：候选池或 trial 数
- `1.0`：安全因子，mode 2/3 中主要保留参数位置
- `3`：unit-weight 模式。`2` 是有放回重采样，快但可能重复事件；`3` 是 reservoir，无重复候选，更适合一般使用
- `1`：ROOT 输出压缩等级

初学者建议：

- 做物理检查：先用默认加权事件。
- 老师明确要求 unit-weight：用 mode `3`。
- 每次改变 unit-weight 设置后，用 `review/scripts/compare_event_samples.C` 和加权样本做比较。

比较命令示例：

```bash
root -l -b -q 'review/scripts/compare_event_samples.C("build/events_weighted.root","build/events_unit.root")'
```

## 7. 正式计算时怎么改参数

上面的 smoke test 参数很小，只是为了检查程序能跑。正式网格要增加 bin 数，例如：

```bash
./build/generate_grid_opt AuAu 200 XnXn build/grid_AuAu200_XnXn_prod.root 30 100 38 4 1.0 20 40 1
```

这会比 smoke test 慢很多，但结果更细。

常用系统和概率文件：

- `AuAu 200` 会读取 `UPC_Probabilities_AuAu_200GeV_Glauber_Final.root`
- `PbPb 5020` 会读取 `UPC_Probabilities_PbPb_5020GeV_Glauber_Final.root`
- `PbPb 5360` 会读取 `UPC_Probabilities_PbPb_5360GeV_Glauber_Final.root`

如果你算别的能量，代码可能没有对应概率表。正式结果不要随便使用 fallback 近似。建议最后一个参数用 `1`，让程序在找不到概率表时直接停止。

## 8. 输出 ROOT 文件里主要有什么

网格文件，例如 `grid_AuAu200_XnXn_smoke.root`：

- `h_dy`：d sigma / dy
- `h_dM`：d sigma / dM
- `h2_dM_dy`：M-y 二维截面
- `h_pxpy_...` 或类似名字：横动量/线偏振相关网格
- 若干 `TParameter`：记录系统、能量、触发、是否读取概率表等元数据

事件文件，例如 `events_bose_1000.root`：

- `Events`：TTree
- `event_weight`：事件权重
- `M_rhoprime`、`y_rhoprime`、`pt_rhoprime`：母粒子运动学
- `m_rho`、`m_sigma`：中间态质量
- `pi1`、`pi2`、`pi3`、`pi4`：四个 pion 的 TLorentzVector
- `pi_charge`：pion 电荷
- `pi_parent`：来自 rho 还是 sigma 的 truth 标记

## 9. 常见错误

### `root-config: command not found`

ROOT 没有安装或环境没加载。先找到 ROOT 的 `thisroot.sh` 并 `source`。

### CMake 找不到 ROOT

同样是 ROOT 环境问题。确认：

```bash
root-config --version
echo $ROOTSYS
```

### `libtbb.so.12: cannot open shared object file`

这是因为正在运行 `build_cmssw_1511` 里的程序，但当前 shell 没有加载 CMSSW_15_1_1 runtime。先执行：

```bash
cd /eos/cms/store/group/phys_heavyions/jianjie/CMSSW_15_1_1/src
cmsenv
cd ../../MC_UPC_RhoPrime4Pi_new
```

或者直接用 `./run_cmssw_1511.sh ./build_cmssw_1511/generate_bose ...`。

### `cannot load required UPC probability ROOT`

你选择的系统/能量没有对应概率表，或者运行目录不对。先确认当前目录下有：

```bash
ls UPC_Probabilities_*_Glauber_Final.root
```

如果要算新能量，需要先用 `calc_UPC_EMD_Glauber_unified.C` 生成概率表，或者让老师确认是否允许使用近似 fallback。

### 程序很慢

先确认你不是把 smoke 参数直接改成了很大的生产参数。网格计算和大统计事件生成都会变慢。先用小参数跑通，再逐步加大。

### 输出文件很大

事件数越多 ROOT 文件越大。测试阶段用 `1000` 或 `10000` 个事件即可。正式生产再加大。

## 10. 建议学生交付什么

每次跑一组结果，请至少保存：

1. 你执行过的完整命令。
2. 终端输出里 `RESULT ...` 那一行。
3. 生成的 grid ROOT 文件名。
4. 生成的 event ROOT 文件名。
5. 分析 ROOT 文件名。
6. 简短说明：系统、能量、触发、事件数、是否 Bose、是否 unit-weight。

不要只发一张图。没有命令和 ROOT 文件，别人很难检查你的结果。
