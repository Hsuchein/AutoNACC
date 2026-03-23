# Amber 非标准残基（非天然氨基酸）模拟完整流程

> 综合整理自多个教程来源：Jerkwin教程、Amber官方Tutorial 5 (GFP/CRO)、
> Amber Tutorial 4 (plastocyanin)、antechamber官方教程(pro4.html)、
> BioExcel AmberTools教程、PyRED/R.E.D. Server文档等。

---

## 一、总体工作流程概览

```
PDB预处理 → 提取非标准残基 → 加帽基(ACE/NME)
    → Gaussian结构优化 → ESP静电势计算
    → RESP电荷拟合(antechamber/resp)
    → prepgen生成残基模板(prepin)
    → parmchk/parmchk2检查缺失参数(frcmod)
    → tleap加载参数、构建体系
    → 溶剂化 + 加反离子
    → 能量最小化 → 升温 → 平衡 → 成品MD → 轨迹分析
```

---

## 二、PDB文件预处理

### 2.1 基本清理

```bash
# 使用pdb4amber自动处理
pdb4amber -i input.pdb -o output.pdb --dry --reduce
```

- `--dry`：删除结晶水
- `--reduce`：添加氢原子并优化位置

### 2.2 残基名修改

根据蛋白质实际状态修改残基名：

| 原始残基 | 修改为 | 含义 |
|---------|--------|------|
| CYS (形成二硫键) | CYX | 参与二硫键的半胱氨酸 |
| CYS (配位金属) | CYM | 去质子化半胱氨酸 |
| HIS (Nδ质子化) | HID | δ-氮质子化组氨酸 |
| HIS (Nε质子化) | HIE | ε-氮质子化组氨酸 |
| HIS (双质子化) | HIP | 双质子化组氨酸(带正电) |
| MSE (硒代甲硫氨酸) | MET | 将SE原子改为SD |

### 2.3 手动编辑要点

- 删除所有已有氢原子（让后续工具重新加氢）
- 删除距蛋白质 >3Å 的游离水分子
- 确认非标准残基的残基名唯一（如SIN, SEP, PPH, CRO, OLS等）
- 添加TER记录分隔不连续的链

---

## 三、非标准残基结构准备

### 3.1 提取非标准残基

从蛋白PDB中提取非标准残基的坐标，用分子编辑器（GaussView、Chimera、MOE、MarvinSketch等）添加氢原子。

### 3.2 加帽基（Capping Groups）

**为什么需要帽基？** RESP电荷拟合需要在蛋白质主链环境中进行，因此将非标准残基制作成二肽形式：

```
ACE - [非标准残基] - NME
```

- **ACE (Acetyl)**：CH₃CO-，N端帽基
- **NME (N-Methylamide)**：-NHCH₃，C端帽基

帽基的作用是模拟蛋白质主链环境，在RESP拟合时帽基原子的电荷被约束为零，拟合完成后去除帽基原子。

### 3.3 准备多构象（推荐）

为提高RESP电荷质量，建议准备2个以上构象：
- **α-螺旋构象**：φ ≈ -60°, ψ ≈ -40°
- **β-折叠构象**：φ ≈ -120°, ψ ≈ -140°

可用TINKER、Amber的distance geometry例程或手动调整二面角生成。

---

## 四、Gaussian量子化学计算

### 4.1 方法一：antechamber自动生成Gaussian输入文件

```bash
antechamber -fi mol2 -fo gcrt -i residue_capped.mol2 -o residue.com -nc <净电荷>
# 或
antechamber -fi mol2 -fo gzmat -i residue_capped.mol2 -o residue.gau
```

- `-fi mol2`：输入mol2格式
- `-fo gcrt`：输出Gaussian卡氏坐标输入文件
- `-fo gzmat`：输出Gaussian Z-matrix输入文件

### 4.2 方法二：手动编写Gaussian输入文件

#### 4.2.1 结构优化 + ESP静电势一步完成

```
%chk=residue.chk
%mem=4GB
%nproc=8
# Opt HF/6-31G* SCF=Tight Pop=MK IOp(6/33=2)

Title: residue optimization and ESP

-1 1
 C     x.xxxxx   y.yyyyy   z.zzzzz
 N     ...
 ...

```

**关键词解释：**

| 关键词 | 含义 |
|--------|------|
| `Opt` | 进行结构优化（能量最小化） |
| `HF/6-31G*` | 使用HF方法和6-31G*基组（**AMBER标准RESP电荷推荐**，适用于ff94/ff99/ff14SB/ff19SB/GAFF） |
| `SCF=Tight` | 使用更严格的SCF收敛标准 |
| `Pop=MK` | 输出Merz-Kollman静电势拟合电荷 |
| `IOp(6/33=2)` | 输出拟合ESP的网格点坐标和精确静电势值（**antechamber读取RESP所必需**） |

**净电荷和自旋多重度行：** `-1 1` 表示净电荷-1、单重态

#### 4.2.2 分两步进行（优化 → ESP）

**第一步：结构优化**
```
%chk=residue.chk
%mem=4GB
%nproc=8
#P B3LYP/6-31G* Opt

Title

0 1
[坐标]

```

**第二步：ESP计算（使用优化后的结构）**
```
%chk=residue.chk
%mem=4GB
%nproc=8
#P HF/6-31G* SCF=Tight Geom=AllCheck Guess=Read Pop=MK IOp(6/33=2,6/41=10,6/42=17)

```

注意：使用`--Link1--`可以在同一个输入文件中串联两步。

#### 4.2.3 额外IOp选项说明

| IOp | 含义 |
|-----|------|
| `IOp(6/33=2)` | 输出ESP网格点和静电势（**必需**） |
| `IOp(6/41=10)` | 每个原子使用10层同心球面点 |
| `IOp(6/42=6)` | 每层点密度≈600点/原子（默认） |
| `IOp(6/42=17)` | 高密度≈2500点/原子（大分子用10以减少计算量） |

#### 4.2.4 部分优化（冻结主链二面角）

对于作为蛋白链一部分的残基，优化时应**冻结φ/ψ主链二面角**：
- 在Z-matrix格式中将对应二面角标记为`F`（frozen）
- 或在笛卡尔坐标格式中使用`Opt=ModRedundant`并添加冻结约束

### 4.3 运行Gaussian

```bash
g16 < residue.com > residue.out
# 或
nohup g16 < residue.com > residue.log &
```

### 4.4 检查计算是否正常完成

在输出文件末尾查找：
```
Normal termination of Gaussian
```

---

## 五、RESP电荷拟合

### 5.1 方法一：antechamber一步完成（最简便）

```bash
antechamber -i residue.out -fi gout -o residue.mol2 -fo mol2 \
  -c resp -s 2 -rn RES -at amber -nc <净电荷>
```

**或输出prepin格式：**
```bash
antechamber -i residue.out -fi gout -o residue.prepin -fo prepi \
  -c resp -s 2 -rn RES -at amber -nc <净电荷>
```

**参数详解：**

| 参数 | 含义 |
|------|------|
| `-i residue.out` | Gaussian输出文件 |
| `-fi gout` | 输入格式为Gaussian输出 |
| `-o residue.mol2` | 输出文件 |
| `-fo mol2` | 输出格式为mol2（也可用prepi） |
| `-c resp` | 使用RESP电荷拟合方法 |
| `-c bcc` | 使用AM1-BCC经验电荷方法（**快速替代方案**） |
| `-s 2` | 详细状态输出级别（2=详细） |
| `-rn RES` | 设定残基名（3个字母） |
| `-at amber` | 使用AMBER原子类型（推荐用于修饰氨基酸） |
| `-at gaff` | 使用GAFF原子类型（推荐用于小分子配体） |
| `-at gaff2` | 使用GAFF2原子类型（更新版GAFF） |
| `-nc -1` | 指定体系净电荷 |
| `-pf y` | 清理中间文件 |

### 5.2 方法二：分步RESP拟合（多构象）

适用于需要多构象RESP拟合的高精度需求：

**步骤1：提取ESP数据**
```bash
espgen -i residue_conf1.out -o conf1.esp
espgen -i residue_conf2.out -o conf2.esp
cat conf1.esp conf2.esp > all.esp
```

**步骤2：生成RESP输入文件**
```bash
respgen -i residue_conf1.mol2 -o resp-step1.respin -f resp1
respgen -i residue_conf1.mol2 -o resp-step2.respin -f resp2
```
> 注意：多构象时respgen生成的文件需要**手动修改**以包含两个构象的信息。

**步骤3：运行两阶段RESP拟合**

**第一阶段**（初始拟合）：
```bash
resp -O -i resp-step1.respin -o resp-step1.respout \
  -e all.esp -t qout_stage1
```

**第二阶段**（约束拟合，固定除甲基/亚甲基碳以外的所有原子）：
```bash
resp -O -i resp-step2.respin -o resp-step2.respout \
  -e all.esp -q qout_stage1 -t qout_stage2
```

**RESP两阶段说明：**

| 阶段 | qwt权重 | 约束对象 |
|------|---------|---------|
| Stage 1 | 0.0005 | 化学等价原子约束相等 |
| Stage 2 | 0.001 | 固定所有原子，仅放开甲基/亚甲基碳；同一碳上的H约束相等 |

**步骤4：将电荷写回mol2**
```bash
antechamber -fi mol2 -fo mol2 -i residue_conf1.mol2 -o residue_resp.mol2 \
  -c rc -cf qout_stage2
```

### 5.3 方法三：使用AM1-BCC快速电荷（简化替代）

```bash
antechamber -i residue.pdb -fi pdb -o residue.mol2 -fo mol2 \
  -c bcc -nc 0 -rn RES -at gaff2
```

AM1-BCC是一种快速的半经验方法，适合初步测试，但**对于发表级别的工作，推荐使用RESP电荷**。

### 5.4 resp.qin文件格式（固定帽基电荷）

当使用分步RESP拟合时，需要准备`resp.qin`文件来固定ACE/NME帽基原子的电荷：

```
# 格式：每行8个值，每个值占10个字符宽度
# 非零值 = 固定电荷（ACE/NME原子）
# 零值 = 待拟合电荷（残基原子）
-0.5696  0.0000  0.3119  0.0000 -0.5577  0.0000  0.0000  0.0000
 0.0000 -0.5157  0.0000  0.3432  0.0000 -0.6123 -0.7580  0.0000
```

在`resp.in`中，对应原子的约束标志：
- `-1`：电荷固定为qin中的值
- `0`：自由拟合
- 正整数n：约束等于第n个原子的电荷

### 5.5 常见RESP错误及修复

| 错误信息 | 原因 | 修复方法 |
|---------|------|---------|
| `Bad value during floating point read` | resp.in中分子数为整数`1` | 改为`1.0000` |
| `End-of-file during read` | resp.in末尾缺少空行 | 在文件末尾添加空行 |
| Fortran读取错误 | resp.in中有注释行 | 删除所有注释行 |
| `At line 403` | AmberTools版本过旧 | 升级到AmberTools 18+ |

### 5.6 关于帽基电荷的处理

在RESP拟合中：
- ACE和NME帽基原子的电荷被**约束为零**（intra-molecular charge constraint = 0）
- 拟合完成后移除帽基原子
- prepgen工具会自动进行**电荷再分配**，确保去除帽基后残基总电荷正确

### 5.5 使用R.E.D. Server（在线工具）

R.E.D. Server提供自动化的RESP电荷拟合服务：
1. 上传PDB文件（含多构象，用MODEL分隔）
2. 设置分子电荷和自旋多重度
3. 选择电荷模型（RESP-A1推荐）
4. 定义帽基约束
5. 服务器自动运行QM计算和RESP拟合
6. 下载mol2/lib文件

---

## 六、力场参数准备

### 6.1 使用prepgen生成残基模板

当非标准残基需要**嵌入蛋白质链**中（不是独立配体）时，需要用prepgen定义残基的连接方式。

**准备主链定义文件（mainchain.mc）：**
```
HEAD_NAME N
TAIL_NAME C
MAIN_CHAIN CA
MAIN_CHAIN C
MAIN_CHAIN N
OMIT_NAME H2
OMIT_NAME HN11
OMIT_NAME OXT
OMIT_NAME HXT
OMIT_NAME CH3_of_ACE
OMIT_NAME CH3_of_NME
PRE_HEAD_TYPE C
POST_TAIL_TYPE N
CHARGE 0.0
```

**各字段解释：**

| 字段 | 含义 |
|------|------|
| `HEAD_NAME` | 连接前一残基的原子（通常为N） |
| `TAIL_NAME` | 连接后一残基的原子（通常为C） |
| `MAIN_CHAIN` | 主链原子名列表 |
| `OMIT_NAME` | 帽基原子（拟合时包含，最终残基中去除） |
| `PRE_HEAD_TYPE` | 前一残基连接原子的类型（通常为"C"） |
| `POST_TAIL_TYPE` | 后一残基连接原子的类型（通常为"N"） |
| `CHARGE` | 残基的目标净电荷 |

**运行prepgen：**
```bash
prepgen -i residue.ac -o residue.prepin -m mainchain.mc -rn RES -rf residue.res
```

- `-i`：输入.ac文件（antechamber coordinate格式）
- `-o`：输出prepin文件
- `-m`：主链定义文件
- `-rn`：残基名（3字母）
- `-rf`：输出残基信息文件

prepgen会自动：
1. 移除OMIT原子
2. 重新分配电荷使总电荷匹配CHARGE值
3. 生成NEWPDB.PDB用于可视化验证

### 6.2 从CIF文件直接处理（PDB化学组分字典）

对于已有CIF定义的残基：
```bash
antechamber -fi ccif -i CRO.cif -bk CRO -fo ac -o cro.ac -c bcc -at amber
```

### 6.3 使用parmchk/parmchk2检查缺失参数

```bash
# 基本用法
parmchk2 -i residue.prepin -f prepi -o residue.frcmod

# 检查与ff14SB的兼容性
parmchk2 -i residue.prepin -f prepi -o residue_ff14SB.frcmod \
  -a Y -p $AMBERHOME/dat/leap/parm/parm10.dat

# 使用GAFF2
parmchk2 -i residue.mol2 -f mol2 -o residue.frcmod -s gaff2
```

**参数说明：**

| 参数 | 含义 |
|------|------|
| `-i` | 输入文件（prepin或mol2） |
| `-f` | 输入格式（prepi/mol2） |
| `-o` | 输出frcmod文件 |
| `-s gaff2` | 使用GAFF2力场数据库 |
| `-a Y` | 输出所有参数（包括完美匹配的） |
| `-p parm10.dat` | 指定参考力场参数文件 |

### 6.4 frcmod文件结构

```
remark goes here
MASS
CU 65.36

BOND
NB-CU  70.000   2.05000
CU-S   70.000   2.10000

ANGLE
CU-NB-CV  50.000  126.700
NB-CU-NB  10.000  110.000

DIHE
X -NB-CU-X  1  0.000  180.000  3.000

IMPROPER
...

NONBON
CU  2.20  0.200
```

**重要：** 检查frcmod文件中是否有`ATTN: needs revision`标记，这些参数是通过类比估计的，可能需要手动修正。

### 6.5 双frcmod策略（ff14SB + GAFF混合）

当残基中同时有标准蛋白质原子类型和新原子类型时：

```bash
# 生成ff14SB参数（可能有缺失）
parmchk2 -i cro.prepin -f prepi -o frcmod1.cro -a Y \
  -p $AMBERHOME/dat/leap/parm/parm10.dat
grep -v "ATTN" frcmod1.cro > frcmod1_clean.cro

# 生成GAFF参数（覆盖缺失项）
parmchk2 -i cro.prepin -f prepi -o frcmod2.cro
```

**在tleap中加载顺序很重要：先加载GAFF参数，再加载ff14SB参数**，后者会覆盖前者中的重复项：
```
loadAmberParams frcmod2.cro   # GAFF先加载
loadAmberParams frcmod1.cro   # ff14SB后加载（覆盖）
```

---

## 七、在tleap/xleap中构建体系

### 7.1 创建残基库文件（.lib/.off）

**方法A：从prepin创建**
```
# 在xleap中
loadamberprep residue.prepin
edit RES                    # 可视化检查
charge RES                  # 检查总电荷
set RES head RES.1.N        # 设置头原子
set RES tail RES.1.C        # 设置尾原子
set RES.1 restype protein   # 设置残基类型
saveoff RES residue.lib     # 保存库文件
```

**方法B：从mol2创建**
```
RES = loadmol2 residue.mol2
set RES restype protein
set RES name "RES"
set RES head RES.1.N
set RES tail RES.1.C
saveoff RES residue.lib
```

**head/tail说明：**
- `head`：连接到前一个残基的原子（通常N）
- `tail`：连接到后一个残基的原子（通常C）
- 如果残基是链末端，设为`null`：`set RES tail null`

**特殊情况—独立配体（非链内残基）：**
```
set RES head null
set RES tail null
```

### 7.2 完整的tleap构建脚本

```bash
# tleap.in
source leaprc.protein.ff14SB      # 蛋白力场
source leaprc.gaff                 # GAFF力场（如果使用gaff原子类型）
source leaprc.water.tip3p          # TIP3P水模型
set default PBRadii mbondi3        # PB半径（用于隐式溶剂）

# 加载非标准残基
loadamberprep residue.prepin       # 或 loadoff residue.lib
loadamberparams residue.frcmod

# 如有多个非标准残基
loadoff sin.lib
loadoff sep.lib
loadoff pph.lib
loadamberparams sin.frcmod
loadamberparams pph.frcmod

# 加载蛋白PDB
x = loadpdb protein.pdb

# 手动添加特殊键（二硫键、金属配位键等）
bond x.29.SG x.45.SG              # 二硫键
bond x.122.SG x.187.SG
bond x.181.OG x.228.P             # 磷酸化位点

# 检查体系
check x
charge x

# 溶剂化
solvateoct x TIP3PBOX 12.0        # 截角八面体水盒子，缓冲距离12Å
# 或
solvatebox x TIP3PBOX 15.0 iso    # 立方水盒子，缓冲距离15Å

# 添加反离子中和体系
addions x Na+ 0                   # 0表示自动中和
addions x Cl- 0

# 添加额外盐浓度（可选，如0.15M NaCl）
# N_ions = 0.0187 × [浓度M] × N_water
addionsRand x Na+ 24 Cl- 24

# 保存拓扑和坐标文件
saveamberparm x system.prmtop system.inpcrd
savepdb x system_solvated.pdb

quit
```

**运行tleap：**
```bash
tleap -s -f tleap.in
```

### 7.3 水模型选择

| 水模型 | leaprc命令 | 溶剂盒名 | 特点 |
|--------|------------|----------|------|
| TIP3P | `source leaprc.water.tip3p` | TIP3PBOX | 最常用，计算快 |
| OPC | `source leaprc.water.opc` | OPCBOX | 更准确，推荐搭配ff19SB |
| SPC/E | `source leaprc.water.spce` | SPCBOX | 适用于某些特殊情况 |
| TIP4P-Ew | `source leaprc.water.tip4pew` | TIP4PEWBOX | 四点水模型 |

### 7.4 solvatebox vs solvateoct

- `solvatebox`：创建长方体/立方体水盒子
- `solvateoct`：创建截角八面体水盒子（**推荐**，原子数减少约29%）

---

## 八、能量最小化

### 8.1 阶段一：约束蛋白质，仅优化溶剂

**输入文件 min1.in：**
```
Minimization Stage 1: Restrain Protein
 &cntrl
  imin=1,           ! 运行最小化
  maxcyc=500,       ! 总步数500
  ncyc=250,         ! 前250步最陡下降法，后250步共轭梯度法
  ntb=1,            ! 周期性边界，恒定体积
  ntr=1,            ! 使用位置约束
  cut=10.0,         ! 非键截断距离10Å
  ntpr=50,          ! 每50步输出能量
  restraint_wt=100.0, ! 约束力常数 100 kcal/mol/Å²
  restraintmask='@CA,C,N,O', ! 约束主链原子
  ! 或 restraintmask=':1-228', ! 约束残基1-228
 /
```

**运行命令：**
```bash
sander -O -i min1.in -o min1.out \
  -c system.inpcrd -p system.prmtop \
  -r min1.rst -ref system.inpcrd
# 或使用pmemd（GPU加速）
pmemd.cuda -O -i min1.in -o min1.out \
  -c system.inpcrd -p system.prmtop \
  -r min1.rst -ref system.inpcrd
```

### 8.2 阶段二：全系统最小化

**输入文件 min2.in：**
```
Minimization Stage 2: Full System
 &cntrl
  imin=1,
  maxcyc=2500,      ! 增加到2500步
  ncyc=1000,        ! 前1000步最陡下降
  ntb=1,
  ntr=0,            ! 无约束
  cut=10.0,
  ntpr=100,
 /
```

```bash
sander -O -i min2.in -o min2.out \
  -c min1.rst -p system.prmtop \
  -r min2.rst
```

---

## 九、分子动力学模拟

### 9.1 升温（Heating）：0K → 300K

**输入文件 heat.in：**
```
Heating: 0 to 300K over 10ps (NVT)
 &cntrl
  imin=0,           ! 运行MD（非最小化）
  irest=0,          ! 不从重启文件读速度（新模拟）
  ntx=1,            ! 仅读坐标
  ntb=1,            ! 周期性边界，恒定体积(NVT)
  ntp=0,            ! 无压力控制
  tempi=0.0,        ! 初始温度0K
  temp0=300.0,      ! 目标温度300K
  ntt=3,            ! Langevin恒温器
  gamma_ln=1.0,     ! Langevin碰撞频率 1.0 ps⁻¹
  ig=-1,            ! 随机数种子（推荐-1自动生成）
  nstlim=5000,      ! 总步数（5000×0.002ps = 10ps）
  dt=0.002,         ! 时间步长2fs
  ntc=2,            ! SHAKE约束含氢键
  ntf=2,            ! 不计算含氢键的力（与ntc=2配合）
  cut=10.0,         ! 非键截断距离
  ntr=1,            ! 使用位置约束
  restraint_wt=10.0, ! 约束力常数 10 kcal/mol/Å²
  restraintmask='@CA,C,N,O',
  ntpr=500,         ! 每500步输出能量
  ntwx=500,         ! 每500步输出轨迹
  ntwr=5000,        ! 每5000步写重启文件
  ioutfm=1,         ! NetCDF格式轨迹
  nmropt=1,         ! 开启NMR约束（用于线性升温）
 /
 &wt
  TYPE='TEMP0',
  ISTEP1=1,
  ISTEP2=5000,
  VALUE1=0.0,
  VALUE2=300.0,
 /
 &wt TYPE='END' /
```

```bash
pmemd.cuda -O -i heat.in -o heat.out \
  -c min2.rst -p system.prmtop \
  -r heat.rst -ref min2.rst \
  -x heat.nc
```

### 9.2 平衡（Equilibration）：NPT 100ps

**输入文件 equil.in：**
```
Equilibration: NPT at 300K, 100ps
 &cntrl
  imin=0,
  irest=1,          ! 重启模拟（读取前一步的速度）
  ntx=5,            ! 读坐标和速度
  ntb=2,            ! 周期性边界，恒定压力(NPT)
  ntp=1,            ! 各向同性压力缩放
  pres0=1.0,        ! 目标压力 1 atm
  taup=2.0,         ! 压力弛豫时间 2ps
  tempi=300.0,
  temp0=300.0,
  ntt=3,
  gamma_ln=1.0,
  ig=-1,
  nstlim=50000,     ! 50000×0.002 = 100ps
  dt=0.002,
  ntc=2,
  ntf=2,
  cut=10.0,
  ntr=0,            ! 无位置约束
  ntpr=1000,
  ntwx=1000,
  ntwr=50000,
  ioutfm=1,
 /
```

```bash
pmemd.cuda -O -i equil.in -o equil.out \
  -c heat.rst -p system.prmtop \
  -r equil.rst -x equil.nc
```

### 9.3 成品MD（Production）

**输入文件 prod.in：**
```
Production MD: NPT at 300K
 &cntrl
  imin=0,
  irest=1,
  ntx=5,
  ntb=2,
  ntp=1,
  pres0=1.0,
  taup=2.0,
  temp0=300.0,
  ntt=3,
  gamma_ln=1.0,
  ig=-1,
  nstlim=500000,    ! 500000×0.002 = 1ns（根据需要调整）
  dt=0.002,
  ntc=2,
  ntf=2,
  cut=10.0,
  ntr=0,
  ntpr=5000,
  ntwx=5000,
  ntwr=500000,
  ioutfm=1,
 /
```

```bash
pmemd.cuda -O -i prod.in -o prod.out \
  -c equil.rst -p system.prmtop \
  -r prod.rst -x prod.nc
```

---

## 十、关键参数对照表

| 参数 | 值 | 含义 |
|------|-----|------|
| `imin` | 0/1 | 0=MD, 1=最小化 |
| `ntb` | 0/1/2 | 0=无周期边界, 1=NVT, 2=NPT |
| `ntt` | 0/1/3 | 0=无温控, 1=Berendsen, 3=Langevin(**推荐**) |
| `ntp` | 0/1 | 0=无压控, 1=各向同性压力缩放 |
| `ntc` | 1/2/3 | 1=无SHAKE, 2=约束含H键(**最常用**), 3=约束所有键 |
| `ntf` | 1/2 | 1=计算所有力, 2=忽略含H键力(**与ntc=2配合**) |
| `dt` | 0.001/0.002 | 时间步长ps（SHAKE开启时可用0.002） |
| `cut` | 8-12 | 非键截断距离Å（显式溶剂通常8-10） |
| `igb` | 0/1/5/8 | 隐式溶剂模型（0=不用, 5=OBC, 8=GBn2） |
| `ioutfm` | 0/1 | 0=ASCII轨迹, 1=NetCDF(**推荐**) |

---

## 十一、轨迹分析

### 11.1 RMSD分析

```
# cpptraj输入文件 rmsd.in
parm system.prmtop
trajin prod.nc
reference min2.rst
rms reference out backbone_rmsd.dat @CA,C,N
run
```

```bash
cpptraj -i rmsd.in
```

### 11.2 平衡检验

使用`process_mdout.perl`提取能量信息：
```bash
process_mdout.perl equil.out prod.out
```

监测指标：
- **势能/动能/总能量**：应趋于稳定
- **温度**：应在300K附近波动
- **压力**：应在1 atm附近波动
- **体积**：应平滑变化
- **密度**：应稳定在约 1.0 g/cm³
- **RMSD**：应缓慢增加后趋于稳定（无剧烈跳跃）

---

## 十二、常见问题与注意事项

### 12.1 原子类型选择

- 修饰氨基酸（链内残基）：使用 `-at amber` → AMBER原子类型（大写，如CT, N, C）
- 独立小分子配体：使用 `-at gaff` 或 `-at gaff2` → GAFF原子类型（小写，如c3, n, ca）
- GAFF小写与AMBER大写不冲突，可以混合使用

### 12.2 电荷方法选择

| 方法 | 精度 | 速度 | 推荐场景 |
|------|------|------|---------|
| RESP (HF/6-31G*) | 最高 | 需Gaussian | 发表级工作 |
| RESP (B3LYP/cc-pVTZ + IEFPCM) | 高 | 需Gaussian | ff03力场 |
| AM1-BCC | 中等 | 快速(无需QM) | 快速测试、大量分子筛选 |

### 12.3 力场搭配

| 蛋白力场 | RESP基组 | 推荐水模型 |
|---------|---------|-----------|
| ff14SB | HF/6-31G* | TIP3P |
| ff19SB | HF/6-31G* | OPC |
| ff03 | B3LYP/cc-pVTZ (IEFPCM ε=4) | TIP3P |

### 12.4 手动校正主链电荷

对于嵌入蛋白链的非标准残基，建议将**主链原子**(N, H, CA, HA, C, O)的电荷替换为ff14SB中标准残基的值（从`$AMBERHOME/dat/leap/lib/amino12.lib`获取），仅保留侧链原子的RESP电荷。

### 12.5 antechamber的局限性

antechamber**只能处理完整分子，不能处理分子片段**（如金属配位位点）。对于片段类残基：
- 使用xleap手动绘制键、分配原子类型和电荷
- 或使用R.E.D. Server进行自动化处理
- R.E.D. Server注意：**不能处理脯氨酸及其衍生物**（因为胺氮连接两个烷基）

### 12.6 预构建的非天然氨基酸力场库

- **Forcefield_NCAA**：包含147种非天然氨基酸的AMBER ff03兼容参数（含β-氨基酸和N-甲基化氨基酸）
- **Forcefield_PTM**：32种翻译后修饰的参数，使用ACE-XXX-NME二肽模型 + 两阶段RESP + GAFF参数借用
- **AMBER Parameter Database** (http://amber.manchester.ac.uk/)：社区贡献的各种残基参数

### 12.7 位置约束 vs BELLY

**永远使用位置约束(ntr=1)**，而非BELLY方法。BELLY已被证明会导致不稳定性和虚假行为。

### 12.6 Langevin vs Berendsen恒温器

- **Langevin (ntt=3)**：在维持和均衡温度方面明显更好，推荐用于加热和平衡
- 但Langevin会影响短时间尺度动力学，对于需要精确动力学的成品MD，可考虑使用Berendsen (ntt=1)

---

## 十三、完整工作流程命令汇总

```bash
# === 1. PDB预处理 ===
pdb4amber -i raw.pdb -o clean.pdb --dry --reduce

# === 2. 提取非标准残基并加帽基（手动操作或脚本） ===
# 使用分子编辑器(GaussView, Chimera)创建 ACE-RES-NME.mol2

# === 3. Gaussian优化 + ESP ===
antechamber -fi mol2 -fo gcrt -i ACE-RES-NME.mol2 -o residue.com -nc <charge>
# 手动编辑residue.com添加 Pop=MK IOp(6/33=2)
g16 < residue.com > residue.out

# === 4. RESP电荷 + 原子类型 ===
antechamber -i residue.out -fi gout -o residue.ac -fo ac \
  -c resp -s 2 -rn RES -at amber -nc <charge>

# === 5. 生成残基模板 ===
# 编写mainchain.mc文件
prepgen -i residue.ac -o residue.prepin -m mainchain.mc -rn RES

# === 6. 检查缺失参数 ===
parmchk2 -i residue.prepin -f prepi -o residue.frcmod

# === 7. 构建体系 ===
tleap -s -f tleap.in
# tleap.in内容见第七节

# === 8. 最小化 ===
pmemd.cuda -O -i min1.in -o min1.out -c system.inpcrd -p system.prmtop -r min1.rst -ref system.inpcrd
pmemd.cuda -O -i min2.in -o min2.out -c min1.rst -p system.prmtop -r min2.rst

# === 9. 升温 ===
pmemd.cuda -O -i heat.in -o heat.out -c min2.rst -p system.prmtop -r heat.rst -ref min2.rst -x heat.nc

# === 10. 平衡 ===
pmemd.cuda -O -i equil.in -o equil.out -c heat.rst -p system.prmtop -r equil.rst -x equil.nc

# === 11. 成品MD ===
pmemd.cuda -O -i prod.in -o prod.out -c equil.rst -p system.prmtop -r prod.rst -x prod.nc

# === 12. 分析 ===
cpptraj -i analysis.in
```

---

## 十四、参考教程来源

1. Jerkwin Amber教程 - 涉及非标准残基的模拟 (SIN/SEP/PPH)
2. Amber Official Tutorial 5 - GFP CRO荧光团参数化
3. Amber Official Tutorial 4 - Plastocyanin铜蛋白
4. Amber antechamber Tutorial (pro4.html) - 4-羟基脯氨酸
5. BioExcel AmberTools Tutorial - GAFF2小分子参数化
6. PyRED/R.E.D. Server Documentation - 自动化RESP电荷拟合
7. Robin Betz Blog - 异肽键参数化
8. AMBER-hub - 自定义核苷酸创建
9. Carlos Ramos Blog - Amber自定义残基参数化
10. FF_PTM - 翻译后修饰力场参数（Princeton）
