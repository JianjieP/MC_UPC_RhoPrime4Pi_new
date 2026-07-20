# MC UPC RhoPrime4Pi Tarball

这个目录提供一个面向 CMS/CRAB 生产的 rho-prime UPC LHE 生成器打包方案。核心目标是把已经在 CMSSW ROOT 环境中编译好的 `rhoprime_lhe`、生产用 ROOT grid、配置 card 和运行入口一起封装成一个可搬运的 tarball。解包后只需要调用 `runcmsgrid.sh`，即可生成标准输出文件 `cmsgrid_final.lhe`。

适用的默认生产通道是：

| 项目 | 默认值 |
| --- | --- |
| 碰撞系统 | `PbPb` |
| 能量 | `sqrt(sNN) = 5360 GeV` |
| neutron class | `NoTag` |
| 过程 | `gamma_A_to_rhoprime_A` |
| 衰变链 | `rho' -> rho0 sigma -> pi+ pi- pi+ pi-` |
| 默认 card | `cards/production/rhoprime_PbPb5360_NoTag.card` |
| 默认 grid | `grids/grid_PbPb5360_NoTag_prod.root` |

## Tarball 里有什么

`packaging/make_tarball.sh` 会把临时 staging 目录中的内容打包到项目根目录。当前打包布局如下：

```text
.
|-- runcmsgrid.sh
|-- README_RUNTIME.md
|-- SHA256SUMS
|-- bin/
|   `-- rhoprime_lhe
|-- cards/
|   `-- rhoprime_PbPb5360_NoTag.card
|-- data/
|   `-- grid_PbPb5360_NoTag_prod.root
`-- metadata/
    |-- build_info.txt
    `-- git_commit.txt
```

各文件用途：

| 路径 | 用途 |
| --- | --- |
| `runcmsgrid.sh` | tarball 解包后的标准运行入口，接收 `NEVENTS RANDOM_SEED NCPU`。 |
| `bin/rhoprime_lhe` | 实际 LHE 生成器可执行文件。 |
| `data/grid_PbPb5360_NoTag_prod.root` | 生产用二维/多维 grid 输入。 |
| `cards/rhoprime_PbPb5360_NoTag.card` | 生产 card，保存物理设置、unweighting 设置和数值参数。 |
| `metadata/git_commit.txt` | 打包时项目 Git commit。 |
| `metadata/build_info.txt` | 打包时的 `SCRAM_ARCH`、`CMSSW_VERSION`、编译器、ROOT 版本和 binary RPATH/RUNPATH。 |
| `SHA256SUMS` | tarball 内关键文件的 SHA256 校验。 |
| `README_RUNTIME.md` | 解包后最小运行说明。 |

## 一键构建和打包

推荐在 CMSSW 环境中完成编译，避免 binary 链到系统 ROOT。项目已经提供了构建记录：

```bash
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /eos/cms/store/group/phys_heavyions/jianjie/CMSSW_15_1_1/src
eval "$(scram runtime -sh)"

cd /eos/cms/store/group/phys_heavyions/jianjie/MC_UPC_RhoPrime4Pi_new

ROOT_BASE=$(readlink -f "$(dirname "$(root-config --libdir)")")

cmake -S . -B build_cmssw_1511 \
  -DROOT_DIR="${ROOT_BASE}/cmake" \
  -DCMAKE_PREFIX_PATH="${ROOT_BASE}:${CMAKE_PREFIX_PATH}" \
  -DCMAKE_BUILD_TYPE=Release

cmake --build build_cmssw_1511 --target rhoprime_lhe -j 4
```

确认 `build_cmssw_1511/rhoprime_lhe` 存在后执行：

```bash
./packaging/make_tarball.sh
```

默认会生成：

```text
rhoprime_PbPb5360_NoTag_v0p1.tgz
```

如果 binary 在其他 build 目录，可以用 `BUILD_DIR` 指定：

```bash
BUILD_DIR=/path/to/build ./packaging/make_tarball.sh
```

打包脚本会做几件保护性检查：

| 检查 | 目的 |
| --- | --- |
| `bin/rhoprime_lhe` 必须存在且可执行 | 防止空包或错包。 |
| `readelf -d` 不应出现 `/usr/lib64/root` | 防止把系统 ROOT 链接的 binary 放进生产包。 |
| 写入 `metadata/build_info.txt` | 记录可复现性信息。 |
| 写入 `SHA256SUMS` | 解包后可校验关键文件。 |

## 解包和运行

在任意工作目录解包：

```bash
mkdir -p run_area
tar -xzf rhoprime_PbPb5360_NoTag_v0p1.tgz -C run_area
cd run_area
```

运行格式：

```bash
./runcmsgrid.sh NEVENTS RANDOM_SEED NCPU
```

例子：

```bash
./runcmsgrid.sh 1000 12345 4
```

成功后会生成：

```text
cmsgrid_final.lhe
```

`runcmsgrid.sh` 会自动完成：

| 步骤 | 说明 |
| --- | --- |
| 参数检查 | `NEVENTS`、`RANDOM_SEED`、`NCPU` 必须是整数；实际生产中 `NEVENTS` 和 `NCPU` 应大于 0。 |
| 输入检查 | 确认 binary、grid、card 都存在。 |
| 临时输出 | 先写 `cmsgrid_final.lhe.tmp`，检查通过后再改名。 |
| 事件数检查 | 统计 `<event>` 数量，必须等于请求的 `NEVENTS`。 |
| XML 检查 | 如果系统有 `xmllint`，会做 stream validation。 |
| SHA256 打印 | 输出最终 LHE 的 SHA256，方便记录和追踪。 |

## 可调参数

参数分为三层：打包参数、运行参数、card/生成器参数。

### 1. 打包参数

| 参数 | 设置方式 | 默认值 | 说明 |
| --- | --- | --- | --- |
| `BUILD_DIR` | 环境变量 | `${PROJECT_ROOT}/build_cmssw_1511` | 指定包含 `rhoprime_lhe` 的 CMake build 目录。 |
| `SCRAM_ARCH` | CMSSW 环境变量 | 未设置时写为 `unset`/`unknown_arch` | 记录到 metadata；用于确认生产架构。 |
| `CMSSW_VERSION` | CMSSW 环境变量 | 未设置时写为 `unset`/`unknown_cmssw` | 记录到 metadata；用于确认 CMSSW 版本。 |
| `VERSION` | 环境变量 | `v0p1` | tarball 名称版本号。 |
| `PACKAGE_TAG` | 环境变量 | `PbPb5360_NoTag` | tarball 名称中的生产标签，例如 `PbPb5360_NoTag_Bose`。 |
| `GRID_SOURCE` | 环境变量 | `${PROJECT_ROOT}/grids/grid_PbPb5360_NoTag_prod.root` | 要放入 tarball 的 grid 源文件；打包后仍复制为 runtime 固定名 `data/grid_PbPb5360_NoTag_prod.root`。 |
| `CARD_SOURCE` | 环境变量 | `${PROJECT_ROOT}/cards/production/rhoprime_PbPb5360_NoTag.card` | 要放入 tarball 的 card 源文件；打包后仍复制为 runtime 固定名 `cards/rhoprime_PbPb5360_NoTag.card`。 |

### 2. 运行参数

`runcmsgrid.sh` 的用户接口固定为：

```bash
./runcmsgrid.sh NEVENTS RANDOM_SEED NCPU
```

| 参数 | 类型 | 例子 | 说明 |
| --- | --- | --- | --- |
| `NEVENTS` | 整数 | `1000` | 要生成的 LHE event 数。 |
| `RANDOM_SEED` | 整数 | `12345` | 随机数种子；生产 job 之间应不同。 |
| `NCPU` | 整数 | `4` | 传给生成器的线程数，即 `--threads`。 |

内部固定传给生成器的路径：

| 生成器参数 | 当前值 |
| --- | --- |
| `--grid` | `${HERE}/data/grid_PbPb5360_NoTag_prod.root` |
| `--config` | `${HERE}/cards/rhoprime_PbPb5360_NoTag.card` |
| `--output` | `${PWD}/cmsgrid_final.lhe.tmp` |

### 3. Card 参数

生产 card 在：

```text
cards/production/rhoprime_PbPb5360_NoTag.card
```

打包后位于：

```text
cards/rhoprime_PbPb5360_NoTag.card
```

当前 card 支持以下具体字段：

| 参数 | 当前值 | 说明 |
| --- | --- | --- |
| `beam_system` | `PbPb` | 碰撞系统记录字段。 |
| `sqrt_snn_gev` | `5360.0` | 核子对质心能量；生成器会用其一半作为默认 beam energy。 |
| `neutron_class` | `NoTag` | neutron tagging 类别。 |
| `process` | `gamma_A_to_rhoprime_A` | 过程记录字段。 |
| `rhoprime_decay` | `rho0_sigma` | rho-prime 衰变模式记录字段。 |
| `rho_decay` | `pi+_pi-` | rho0 衰变模式记录字段。 |
| `sigma_decay` | `pi+_pi-` | sigma 衰变模式记录字段。 |
| `rapidity_min` | `-2.4` | 生成相空间记录字段。 |
| `rapidity_max` | `2.4` | 生成相空间记录字段。 |
| `bose_symmetrization` | `false` | 是否启用 Bose symmetrization；设为 `true` 时会把相同电荷 pion 的交换振幅加入权重计算。 |
| `rho_mass_gev` | `0.77526` | rho0 质量记录字段。 |
| `rho_width_gev` | `0.1491` | rho0 宽度记录字段。 |
| `sigma_mass_gev` | `0.500` | sigma 质量记录字段。 |
| `sigma_width_gev` | `0.400` | sigma 宽度记录字段。 |
| `event_mode` | `weighted` | `weighted` 保留 event weight；`unweighted` 输出 unit-weight LHE。 |
| `unweighting_mode` | `1` | unit-weight 生成模式：`1=accept-reject`，`2=resample`，`3=reservoir`。 |
| `unweighting_trials` | `200000` | unweighting 试验次数或 pool size。 |
| `unweighting_safety_factor` | `1.25` | mode 1 envelope 安全因子。 |
| `grid_file` | `data/grid_PbPb5360_NoTag_prod.root` | 直接运行 `rhoprime_lhe` 时的 grid fallback；`runcmsgrid.sh` 会显式传 `--grid`。 |
| `decay_norm_trials_per_mass` | `20000` | 每个质量 bin 的 decay normalization trials。 |
| `output_compression_level` | `1` | 中间 ROOT 输出压缩级别。 |
| `process_id` | `81` | 写入 LHE 的 process id。 |

注意：其中部分字段是记录/文档字段，只有生成器源码中解析的字段会改变运行行为。当前会直接影响生成器行为的字段包括：

```text
bose_symmetrization
event_mode
unweighting_mode
unweighting_trials
unweighting_safety_factor
grid_file
decay_norm_trials_per_mass
output_compression_level
process_id
sqrt_snn_gev
beam_energy_gev
```

## 制作 Bose Tarball

`bose_symmetrization` 会直接影响生成器行为。代码在读取 card 时会把 `bose_symmetrization` 或 `bose_symmetrize` 写入 `cfg.bose_symmetrize`；当它为 `true` 时，事件权重会包含四种 pion 配对振幅，而不是只使用默认的 `(pi1, pi2) -> rho`、`(pi3, pi4) -> sigma` 配对。

项目提供了一个 Bose card：

```text
cards/production/rhoprime_PbPb5360_NoTag_Bose.card
```

制作 Bose tarball：

```bash
CARD_SOURCE=/eos/cms/store/group/phys_heavyions/jianjie/MC_UPC_RhoPrime4Pi_new/cards/production/rhoprime_PbPb5360_NoTag_Bose.card \
PACKAGE_TAG=PbPb5360_NoTag_Bose \
./packaging/make_tarball.sh
```

输出文件：

```text
rhoprime_PbPb5360_NoTag_Bose_v0p1.tgz
```

解包后运行方式不变：

```bash
tar -xzf rhoprime_PbPb5360_NoTag_Bose_v0p1.tgz -C run_area
cd run_area
./runcmsgrid.sh 1000 12345 4
```

如果只是想临时测试，也可以直接编辑 tarball 内或项目内的 card，把：

```text
bose_symmetrization = false
```

改为：

```text
bose_symmetrization = true
```

然后重新打包或直接运行生成器。

### 4. 直接运行生成器时的命令行参数

如果不通过 `runcmsgrid.sh`，也可以直接运行：

```bash
./bin/rhoprime_lhe \
  --grid data/grid_PbPb5360_NoTag_prod.root \
  --config cards/rhoprime_PbPb5360_NoTag.card \
  --events 1000 \
  --seed 12345 \
  --threads 4 \
  --output cmsgrid_final.lhe
```

可用命令行覆盖项：

| 参数 | 说明 |
| --- | --- |
| `--grid PATH` | 输入 ROOT grid。 |
| `--config PATH` | 输入 card。 |
| `--events N` | event 数，必须大于 0。 |
| `--seed SEED` | 随机数种子。 |
| `--threads N` | 线程数，必须大于 0。 |
| `--output PATH` | 输出 LHE 文件。 |
| `--beam-energy GEV` | 每束 beam energy，默认 `2680`。 |
| `--weighted` | 写 weighted events。 |
| `--unweighted` | 写 unit-weight events。 |
| `--unweighting-mode MODE` | `1=accept-reject`，`2=resample`，`3=reservoir`。 |
| `--unweighting-trials N` | unweighting trials 或 pool size。 |
| `--unweighting-safety FACTOR` | mode 1 envelope 安全因子。 |
| `--decay-norm-trials N` | 每个 mass bin 的 decay normalization trials。 |

## 切换生产通道

如果要从 `PbPb5360_NoTag` 切到其他通道，例如 `XnXn`，需要同步替换三类文件/字段：

| 位置 | 需要修改 |
| --- | --- |
| `cards/production/*.card` | `neutron_class`、`grid_file`、必要的物理参数。 |
| 打包命令 | 用 `PACKAGE_TAG`、`GRID_SOURCE`、`CARD_SOURCE` 指定 tarball 标签、grid 源文件和 card 源文件。 |
| `packaging/runcmsgrid.sh` | 通常不需要改；打包脚本会把选中的 grid/card 复制到 runtime 固定路径。只有想改变 tarball 内部文件名时才需要同步修改。 |

替换后重新编译或确认 binary 可用，再重新运行 `./packaging/make_tarball.sh`。

## 常见问题

### `Generator executable not found`

说明 `BUILD_DIR/rhoprime_lhe` 不存在或不可执行。先重新构建：

```bash
cmake --build build_cmssw_1511 --target rhoprime_lhe -j 4
```

或指定正确的 build 目录：

```bash
BUILD_DIR=/path/to/build ./packaging/make_tarball.sh
```

### `Refusing to package ... linked to system ROOT`

说明 binary 链到了 `/usr/lib64/root`，生产环境不推荐。请进入 CMSSW `cmsenv`/`scram runtime` 环境后重新 CMake configure 和 build。

### `Requested N events, but found M`

`runcmsgrid.sh` 在最终输出前会检查 `<event>` 数量。这个报错表示生成器实际写出的事件数和请求值不一致，应保留日志并检查 seed、grid、card 和生成器输出。

### `xmllint is unavailable`

这只是警告。没有 `xmllint` 时会跳过 XML stream validation，但仍会检查输出非空和 event 数量。
