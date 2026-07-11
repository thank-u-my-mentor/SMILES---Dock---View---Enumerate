
  
> 目的：这不是长篇教程，而是一份“每次开跑前照着查”的操作代码本。核心参考是 Amber 官方 MCPB.py heme 教程：<https://ambermd.org/tutorials/advanced/tutorial20/mcpbpy_heme.php>。官方例子是 heme；我们的体系是 non-heme Fe + His/Glu + ACT/UNK/UNL radical 中间体，所以流程骨架相同，但化学判断必须更严格。

# 0. 总流程
人工/文献给定活性中心构象
- -> 清理 PDB 和 ligand 名称
- -> 明确 Fe 配位 donor、His 质子化、ligand 电荷/自由基/自旋
- -> 准备 mol2/frcmod
- -> MCPB.py -s 1 生成 small/large model
- -> GAMESS/Gaussian 做 small_opt、small_fc/Hessian、large_mk
- -> MCPB.py -s 2/3/4 生成 Fe bonded force field
- -> tleap/Amber -> GROMACS
- > EM/NVT/NPT/production
- > PBC 修复 + RDC/距离/Fe-donor 审计
# 1. Official 硬规则
Amber 官方教程强调，最容易**出错**的是第一步：`PDB` 和 `mol2` 准备。
必须检查：
- PDB 使用 Amber 可识别命名：**HID/HIE/HIP**，而不是模糊 **HIS**。
	- 每个 residue 内 atom name 必须唯一，尤其 ligand。
	- 金属原子 residue name 和 atom name 都用大写，例如 FE/FE。
	- 金属原子要单独作为一个 residue，不要藏在 ligand residue 里。
	- 非标准残基，包括金属、配位水、底物、radical ligand，都要有 mol2。
	- mol2 文件名、residue name、PDB residue name 大小写要一致。

官方会用 `reduce` 给普通 ligand 加氢，但也明确说要人工检查氢是否正确。我们的 radical 中间体不能盲用 `reduce`，因为它会倾向于把缺氢自由基“补成普通饱和分子”。

# 2. 我们体系和官方 heme 例子的差异

Fe 在 heme 里，配位关系相对清楚。
HEM 的 charge 可以按质子化状态定义，例如 -4。
Gaussian 或 GAMESS 做 QM。
MCPB.py 用 **Hessian/Seminario** 生成**金属键和角参数**。
## 我们的 Fe-radical 中间体：
- Fe 是 non-heme iron。
- Fe 通常要明确六配位：2 His + Glu/Asp + ACT/UNK/UNL/HOH 等。
- UNK/UNL 可能是 radical：少一个 H，且 multiplicity 不是普通 singlet。
- His 的 Fe-donor N 不能带 H；另一侧 imidazole N 通常要带 H。
- charge 只解决总电荷，不解决 radical；radical 还需要正确 atom count 和 MULT。
# 3. His 质子化：我们已经踩过的大坑
## Fe 如果配位到 His 的 `NE2`，常规中性 His 应该是：
**resname = HID**
- ND1-H 存在
- NE2-H 不存在
- NE2 -> Fe 配位

## Fe 如果配位到 `ND1`，常规中性 His 应该是：
**resname = HIE**
NE2-H 存在
ND1-H 不存在
ND1 -> Fe 配位

***不应该默认接受***：
```text

两个 N 都无 H：更像 imidazolate/负电模型，必须有明确理由。
两个 N 都有 H：HIP，带正电，一般不是 Fe donor 默认状态。

```

# 19. 收尾整理与 TER/OXT 审计规则

## 19.1 为什么 MCPB.py 可能在连续主链中间写入很多 TER

`TER` 在 PDB 语义里表示一条 polymer chain 结束。`tleap` 会严格相信它：如果一个普通氨基酸后面出现 `TER`，即使下一个残基的 N 原子就在正常肽键距离内，`tleap` 也会把前一个残基当作 C-terminal residue，并可能自动补 `OXT`。

MCPB.py 生成 `TJ2unlNH_mcpbpy.pdb` 时会重写金属中心附近的 residue name、custom atom type、标准/非标准残基块，并把水、配体、金属模型拆成多个 segment。若输入 PDB 的 chain ID、residue number、insertion code、原始 TER 或金属位点重命名不规整，MCPB.py 输出 PDB 可能按内部“残基块”机械地插入 TER。这不是可靠的化学判断。

固定规则：

```text
不要直接相信 MCPB.py 输出 PDB 的 TER。
在 tleap 前必须审计 TER：
如果 TER 前一残基 C 与后一残基 N 距离 <= 1.75 A，且同属蛋白连续主链，应删除该 TER。
水、离子、配体、真实 chain break 的 TER 保留。
```

本项目实测：

```text
原 TJ2unlNH_mcpbpy.pdb: 178 个 TER
错误后果: TRP71 被 tleap 当成 C 端并补 OXT；TRP71 OXT 与 GLY72 N 距离约 0.14 A
结果: GROMACS EM 出现 infinite force
修复: cleanTER 版本删除 peptide-like TER，保留真实分段 TER
验证: tleap Errors = 0；EM/NVT/NPT 通过；Fe-donor 键和主链连接在 prmtop 中存在
```

固定脚本：

```bash
python3 /mnt/e/Codex/skills/mcpb-md/scripts/clean_mcpb_pdb_ter_by_geometry.py \
  --in TJ2unlNH_mcpbpy.pdb \
  --out TJ2unlNH_mcpbpy.cleanTER.pdb

python3 /mnt/e/Codex/skills/mcpb-md/scripts/check_mcpb_topology_bonds.py \
  TJ2unlNH_solv.prmtop TJ2unlNH_solv.inpcrd
```

## 19.2 ANTECHAMBER_*、ATOMTYPE.INF、sqm.* 是什么

这些通常是 AmberTools `antechamber` / `sqm` / `parmchk2` 的中间文件：

```text
ANTECHAMBER_*: antechamber 的 atom typing、bond typing、临时 mol2/ac/mol 输出
ATOMTYPE.INF: atom type assignment 的辅助信息
sqm.in/out/pdb: AM1-BCC 或 semiempirical charge 计算输入/输出/结构
```

它们不是最终 MD 必须直接读取的文件。最终 MD 更关键的是：

```text
*.mol2
*.frcmod
MCPB.py 生成的 *_mcpbpy.frcmod
tleap 生成的 prmtop/inpcrd
GROMACS top/gro/tpr/cpt/xtc/edr/log
```

整理原则：

```text
如果最终 mol2/frcmod/prmtop/inpcrd 已经稳定，ANTECHAMBER_* 与 sqm.* 可归档。
不要在没有保存 mol2/frcmod 和命令记录时直接删除它们。
radical/open-shell ligand 不应盲目信任 sqm/AM1-BCC；最终 charge 应以 QM/RESP 或明确记录的修正方案为准。
```

## 19.3 本项目最小保留集

当前整理目录：

```text
E:\TJ\260706_mcpb_m2_200ns\project_organized_20260708
```

保留逻辑：

```text
00_inputs:
  mcpb.in, mcpb_original.pdb, ligand/residue mol2/frcmod, small/standard/large PDB/fingerprint

01_qm_success_M2:
  成功 M2 OPT input/log
  成功 M2 Hessian input/log

02_mcpb_s2_s4_forcefield:
  *_mcpbpy.frcmod, cleanTER PDB/report, Amber prmtop/inpcrd, tleap logs

03_charge_model_smallM2:
  small-M2 MK/RESP input/log/ESP, RESP files, charge mapping reports

04_md_current_scripts:
  mdp/top/itp/run scripts and Fe-distance checkpoint files

98_failed_attempts_index:
  少量失败/混合自旋参考，不保存全部试错垃圾
```

正在运行的 production 目录不要移动：

```text
E:\TJ\260706_mcpb_m2_200ns\md_200ns
```

## 19.4 skill 命名建议

不要把 PLE 重命名成 MCPB 工具。PLE 更适合保留为 pocket/ligand exploration 相关 skill。

金属中心参数化应使用单独 skill：

```text
mcpb-md
```

理由：

```text
MCPB/GAMESS/Amber/GROMACS 是重计算、文件审计和参数化工作流；
PLE 是 docking/pocket/ligand exploration 语义；
混名会降低以后独立运行时的可维护性。
```

PLE 里现在有审计/修复脚本：
```bash

python /mnt/e/Codex/skills/pocket-ligand-explorer/scripts/fix_fe_histidine_protonation.py \
  --pdb input_complex.pdb \
  --out input_complex_fe_his_fixed.pdb \
  --cutoff 2.8

python /mnt/e/Codex/skills/pocket-ligand-explorer/scripts/audit_mcpb_histidine_protonation.py \
  --pdb input_complex_fe_his_fixed.pdb \
  --cutoff 2.8

```
# 4. ligand / radical 的 mol2 原则
普通 ligand 可以这样做：
```bash
antechamber -fi pdb -fo mol2 -i ACT.pdb -o ACT.mol2 -c bcc -pf y -nc -1
parmchk2 -i ACT.mol2 -f mol2 -o ACT.frcmod -s gaff2
```

这里 `-nc` 是 net charge，例如：
```text
ACT acetate:      -nc -1
普通中性 ligand:  -nc 0
```
***但 radical ligand 要额外注意***：
```text

如果苄位 C radical，就必须真的删掉那个 C 上的一个 H。
如果 N 要和 Fe 配位，那个 Fe-bound N 通常不该再多带一个错误的 H。
charge 不是 radical 的全部定义；MULT/电子数同样关键。
antechamber/sqm 对 radical anion 可能失败，不能当成全自动真理。

```


实际建议：
```text

1. 先用 PDB/SDF/mol2 明确 atom names、bond order、缺失的 radical H。
2. mol2 提供键级、atom type、初始 charge。
3. MCPB/QM 层用 total charge + MULT 描述整个 Fe 小模型电子态。
4. 真正严肃 charge 最好走 QM/RESP，而不是完全依赖 AM1-BCC。

```
# 5. Fe 六配位 donor 清单必须手工确认
不要让 MCPB “猜化学”。先明确 donor 表，例如：
```text
HID A161 NE2
HID A241 NE2
GLU A320 OE1
UNK A501 N1
UNK A501 O2
ACT A502 O2 或 O1
可选：HOH A500 O

```
Fe-donor 距离通常先放在约 1.95=2.3 A。太短如 1.6-1.7 A，或太长如 3 A 以上，都要回到 PyMOL/化学模型检查。
# 6. MCPB    step 1：生成模型，不是最终参数
`MCPB.py -s 1` 做的事：
```text
读取 mcpb.in 和 mcpb_original.pdb。
生成 small model、standard model、large model。
生成 Gaussian/GAMESS 输入。
生成 fingerprint，记录 Fe-link 信息和新 atom type。
```
## 示例命令：
```bash
source /mnt/l/WSL/conda_envs/AmberTools25/amber.sh
cd /mnt/e/ZFY/mcpb_5for_tleap_bidentate_M7_HIDfixed_20260624_v1/mcpb
MCPB.py -i mcpb.in -s 1

```
## step 1 后立刻检查：
```text

small model atom/electron count 是否与 charge/MULT 匹配。
Fe 六配位是否还在。
HID/HIE 是否有正确 N-H。
UNK/ACT atom names 是否和 mol2 一致。

```
## 电子奇偶规则：
```text
total_electrons = sum(atomic_numbers) - total_charge
偶数电子 -> MULT 通常必须是奇数：1, 3, 5, 7...
奇数电子 -> MULT 通常必须是偶数：2, 4, 6...

```
Fe(III) high-spin + organic radical 常见候选：
```text
MULT=5：反铁磁耦合候选
MULT=7：铁磁耦合候选
```
# 7. GAMESS 三个 QM 任务分别解决什么
## 7.1 fixed-geometry SCF / ENERGY
只问：当前坐标、电荷、自旋下，电子结构能不能自洽？
1. **不移动**原子。
2. **不生成**优化后结构。
3. **成功标志是 DENSITY CONVERGED**，且没有 SCF IS UNCONVERGED。
```bash

tail -f gamess_logs/JOB.log

```
## 7.2 small_opt / OPTIMIZE
问：电子结构自洽后，原子坐标能不能放松到局部低能构象？
```text
会移动 small model 原子坐标。
可能轻微改变 Fe-donor 距离。
如果 Fe 配位断了，不能继续 Hessian。

```
成功标志：
```text
......END OF GEOMETRY SEARCH......
EXECUTION OF GAMESS TERMINATED NORMALLY
梯度足够小
Fe-donor 没断
```
## 7.3 small_fc / Hessian
问：在优化后构象附近，原子位移会带来多大能量惩罚？
```text
这一步给 MCPB/Seminario 提供 Fe-bond/angle 力常数。
它非常慢，因为要算二阶导数/力常数。
它不是 production MD，也不是采样。
```
## 7.4 large_mk
问：更大模型的静电势如何，用于 RESP/MK charge fitting。
```text
通常不做完整大模型优化，主要为了电荷。
MCPB.py step 3 会用它生成 metal-site residue mol2。
```
# 8. GAMESS 官方对应命令骨架
官方教程说 GAMESS 可以替代 Gaussian。典型流程是：
```bash
rungms GROUP_small_opt 00 8 > gamess_logs/GROUP_small_opt.log 2>&1
# 把 small_opt 的优化后坐标填入 small_fc.inp
rungms GROUP_small_fc 00 8 > gamess_logs/GROUP_small_fc.log 2>&1
rungms GROUP_large_mk 00 8 > gamess_logs/GROUP_large_mk.log 2>&1
```

PLE 里通常用脚本封装，避免 terminal 引号错误：

```bash
cd /mnt/e/ZFY/mcpb_5for_tleap_bidentate_M7_HIDfixed_20260624_v1/mcpb
bash run_ZFY5B_HID_M7_FE_scf_energy_soscf_finegrid_tmux.sh
tail -f gamess_logs/ZFY5B_HID_M7_FE_scf_energy_soscf_finegrid.log

```

我建议以后少手敲长命令，多运行生成好的 `run_*.sh`。这比在 shell 里写一大串引号可靠很多。

# 9. MCPB step 2/3/4
**QM 都完成后：**
```bash
source /mnt/l/WSL/conda_envs/AmberTools25/amber.sh
cd /path/to/mcpb

MCPB.py -i mcpb.in -s 2 --logf GROUP_small_fc.log
MCPB.py -i mcpb.in -s 3 --logf GROUP_large_mk.log
MCPB.py -i mcpb.in -s 4

```
**含义**
```text
-s 2：从 Hessian 生成 Fe bond/angle 参数，写 mcpbpy.frcmod。
-s 3：从 large_mk 做 charge fitting，生成 custom residue mol2。
-s 4：生成 mcpbpy.pdb 和 tleap 输入。
```
# 10. tleap 前必须审计 custom residue
`tleap` 不够聪明。*它可能接受一个语法正确但化学断裂的 mol2* 
我们踩过的坑：
```text
HD1.mol2 缺 CA-CB 键。
tleap 不报错，因为它认为用户模板就是这样。
GROMACS 也照跑，最后 His 侧链漂走。
```
必须检查：
- **HIS**-like custom residue:
- **N-CA, CA-C, C-O, CA-CB, CB-CG, ring bonds** 都要在。
- **Fe-bound His** 如果是中性 **HID/HIE**，要有正确的一个 **imidazole N-H**。
- **GLU/ASP-like** custom residue:
- **CA-CB, CB-CG, carboxylate bonds** 要在。
- **Fe-donor bond** 要在 **topology** 里。
## PLE 审计命令：
```bash
python /mnt/e/Codex/skills/pocket-ligand-explorer/scripts/audit_mcpb_custom_residue_bonds.py \
  --mcpb-dir /path/to/03_mcpb \
  --pdb /path/to/03_mcpb/GROUP_mcpbpy.pdb
  
python /mnt/e/Codex/skills/pocket-ligand-explorer/scripts/audit_mcpb_histidine_protonation.py \
  --mol2 /path/to/03_mcpb/HD1.mol2 /path/to/03_mcpb/HD2.mol2 \
  --mol2-donor NE2
```
# 11. Amber/GROMACS 前处理
```bash
tleap -s -f GROUP_tleap.in > tleap.out
acpype -p GROUP_solv.prmtop -x GROUP_solv.inpcrd
```
***EM  ->  NVT  ->  NPT  ->  production***

每一步后检查：

EM 不爆炸。
NVT/NPT 温度压力正常。
Fe-donor 距离没有突然飞走。
custom His 的 CA-CB 仍约 1.5 A。
***轨迹结束后自动生成 nojump/centered/compact 版本***，不直接看 raw xtc。

# 12. 官方流程里没有彻底解决 radical 怎么办
Amber 官方 heme 教程没有给出“radical ligand 一键处理”的答案。它给的是框架：
- 非标准残基自己准备 mol2/frcmod。
- charge 根据实际质子化/氧化态指定。
- QM 计算的 charge/spin 是复杂问题，需要按体系判断。
- 若 metal coordination 在优化后断裂，要换方法/基组/策略。
*所以我们的 radical 规则是*：
```text
radical 必须体现在结构上：少一个 H。
radical 必须体现在电子态上：正确 charge + MULT。
radical ligand 不应被 reduce/OpenBabel 自动补氢。
antechamber 只能做格式/初始 GAFF 检查，不能证明 radical 化学正确。
真正严肃的 radical charge 最好来自 QM/RESP 或文献参数。
```
# 13. 当前经验总结
1. **Hessian/Fe force constants** 可以是有用。
2. 但不能盲目把 **small model** 优化坐标 **patch** 回完整蛋白。
3. **His CA-CB** 被拉坏会导致 MCPB mol2 缺键。
4. Fe-bound His 不能丢失 **imidazole N-H**。
# 14. 每次开跑前的最短 checklist
- [ ] PDB 里 Fe 是 FE/FE，单独 residue。
- [ ] Chain 选对，没有 A/B 混淆。
- [ ] 老 ligand/阻碍物删掉；真正参与反应的 ligand 加回。
- [ ] H++/人工质子化已经转成 Amber 名称。
- [ ] Fe-bound His 的 HID/HIE 和 N-H 位置正确。
- [ ] ligand residue 内 atom name 唯一，和 mol2 致。
- [ ] radical H 确实删除，没有被自动补回。
- [ ] ACT/UNK/UNL charge 合理，MULT 和电子奇偶匹配。
- [ ] MCPB step1 后 small/large model electron count 重新算过。
- [ ] small_opt 后 Fe-donor 没断。
- [ ] Hessian 后再跑 MCPB -s 2/3/4。
- [ ] tleap 前审计 custom mol2 键和 His N-H。
- [ ] GROMACS EM/NVT/NPT 后审计 Fe-donor 和 CA-CB。
- [ ] production 后先做 PBC 修复，再做 PyMOL/RDC。
# 15. 失败时怎么判断下一步

SCF IS UNCONVERGED：电子态没站稳；别跑 Hessian。

FAILURE TO LOCATE STATIONARY POINT：优化没有找到稳定点；检查梯度和 Fe-donor。

Fe-donor 断裂：换化学模型/约束/方法，不要继续参数化。

custom mol2 缺 CA-CB 或 N-H：回到 MCPB/tleap 前修复，不要 production。

tleap 没报错不代表化学正确。

PyMOL 看 raw xtc 很怪：先做 PBC nojump/centered，不要直接下结论。

# 16. open-shell Fe large_mk 急救流程：不要只看 MCPB 是否接受

这一节是为了在没有 AI 帮忙时，自己判断 large_mk 是否值得继续、是否能进入 MCPB.py `-s 3`。

## 16.1 先定死一个原则

```text
small_opt / Hessian 用 M2，最终 large_mk/RESP 也必须用 M2。
M2 Hessian + M6 large_mk 只能当预览，不适合作为最终 MD 参数。
```

原因是：

```text
Hessian/Seminario 决定 Fe 周围 bonded force constants。
large_mk/RESP 决定 metal-site residue charges。
二者描述的应该是同一个电子态。
```

旧 TJ/M7 目录里 large_mk 曾经被 MCPB.py 接受，但那不等于 QM 干净成功。必须自己检查 log：

```bash
grep -E "SCF IS UNCONVERGED|DENSITY CONVERGED|DIIS CONVERGED|FINAL .*ENERGY|S-SQUARED|ELECTROSTATIC POTENTIAL|TERMINATED" GROUP_large_mk.log
```

如果出现：

```text
SCF IS UNCONVERGED
FINAL ... ENERGY IS 0.0000000000
S-SQUARED 离目标自旋很远
```

即使 MCPB.py `-s 3` 后说 charges and atom numbers match，也只能说明文件格式被解析了，不能说明电子态合理。

## 16.2 M2 large_mk 的推荐分阶段路线

不要一上来就让大模型从 Hückel 初猜直接跑 RO/U-B3LYP ESP。对 Fe radical/open-shell 大模型，优先分阶段：

```text
阶段 A：M2-only ROHF/HF checkpoint
目的：得到一个真正来自 M2 的可读回 $VEC。
不做 DFT，不做 ESP，不拟合电荷。

阶段 B：M2-only ROHF/HF refine
目的：从阶段 A 的 $VEC 继续把 SCF 收稳。

阶段 C：M2 RO-B3LYP large_mk
目的：从阶段 B 的 M2 $VEC 读入，打开 DFT + ELPOT/PDC，生成 ESP。
```

当前 TJ M2 的经验：

```text
M6 large_mk 可以收敛，但不能作为最终 M2 电荷。
M2 from M6 $VEC 虽然能正常结束，但 S-SQUARED = 6.886，说明被高自旋/污染态拖走，拒绝作为最终结果。
M2-only ROHF/HF 从 Hückel 直接跑到 70 步附近会震荡，所以不能无限等。
```

## 16.3 什么时候停止一条 SCF

不是看“程序还在跑”，而是看趋势。

可以继续看：

```text
DENSITY CHANGE 总体下降到 1e-2 以下，并继续下降。
DIIS ERROR 总体下降。
能量不再几十到几百 hartree 来回跳。
```

应该换策略：

```text
40-70 步后 DENSITY CHANGE 还在 0.03-0.3 之间震荡。
DIIS ERROR 卡在 1e-2 到 1e-1，不再下降。
能量反复大跳。
S-SQUARED 远离目标自旋。
```

重要：探路 job 不要手动 `pkill` 当作常规流程。更好的做法是在输入里设定探路步数，例如 `MAXIT=70`，让 GAMESS 正常结束。正常结束更可能写出可复用的 `$VEC`。

## 16.4 判断 `$VEC` 是否真的可用

只有 `.dat` 里有完整 `$VEC`，才可以 `GUESS=MOREAD`。

```bash
ls -lh /home/qin/softwares/gamess/restart/GROUP.dat
grep -n '\$VEC' /home/qin/softwares/gamess/restart/GROUP.dat | head
```

经验判断：

```text
十几 KB 的 .dat 通常只是输入/空壳，不是完整 large model $VEC。
几十 MB 的 .dat 才像完整 large model orbital vector。
没有完整 $VEC 时，不要假装能 MOREAD，只能重新跑 checkpoint。
```

如果一条 SCF 在某个中间步明显最好，不要手动掐进程来“保存它”。正确做法是：

```text
1. 从 log 里找出最好的一步，例如 energy 平稳、density change 或 orbital residual 最低的 NSCF step。
2. 复制同一个输入，把 MAXIT 改成这个步数，例如 MAXIT=33。
3. 重新跑，让 GAMESS 自己正常结束。
4. 检查新 .dat 是否有完整 $VEC。
5. 用这个 $VEC 开 MOREAD refine。
```

这比人工 `pkill` 可靠得多，因为 SIGTERM 经常只留下 10 KB 级别的空壳 `.dat`。

ROHF checkpoint 可以用来生成 M2 UHF 初猜，但不能直接给 UHF `MOREAD`。UHF 通常需要 alpha/beta 两套轨道；如果直接把 ROHF `$VEC` 给 UHF，可能报：

```text
PREMATURE END OF ORBITAL INPUT ENCOUNTERED
```

可行的急救办法是把 ROHF `$VEC` 轨道块复制两遍，放在同一个 `$VEC ... $END` 里，作为 UHF 的 alpha/beta 初猜。这样做的用途是：

```text
ROHF/HF 负责产生 M2-only、S^2=0.750 的干净初猜。
duplicated $VEC 让 UHF 能从这个 M2 初猜出发，而不是借 M6 高自旋历史。
U-B3LYP large_mk 再从这个 duplicated M2 $VEC 开始生成 ESP。
```

这个方法仍然必须通过 CHECK 和最终 S-SQUARED 审计；它不是绕过电子态判断的捷径。

## 16.5 GAMESS 输入硬规则

多重态不是 MCPB.py 自动猜的，必须在输入和 `mcpb.in` 里一致。

```text
mcpb.in:
smmodel_spin 2
lgmodel_spin 2

GAMESS:
SCFTYP=ROHF 或 UHF
ICHARG=总电荷
MULT=2
```

如果 `MULT>1` 却忘了写 `SCFTYP=UHF/ROHF`，GAMESS 会默认 RHF，然后报：

```text
SCFTYP=RHF MUST HAVE MULT=1
```

如果遇到 `SYMDIA` 或 `EINVIT` 对角化错误，先加：

```text
$SYSTEM KDIAG=3 $END
```

如果 Hessian 或位移点 SCF 容易坏，默认更保守：

```text
$SCF DIRSCF=.T. DIIS=.T. SOSCF=.F. DAMP=.T. SHIFT=.T. MAXIT=200 $END
```

如果 OPT/large_mk 要减少数值噪音，优先从一开始固定 fine grid，而不是中途换网格。

## 16.6 当前项目可直接照抄的监控命令

```bash
cd /mnt/e/TJ/260706_mcpb_m2_200ns/succeed_MCPBpy_project

tail -f logs/JOB.log

grep -E "DENSITY CONVERGED|DIIS CONVERGED|SCF IS UNCONVERGED|FINAL .*ENERGY|S-SQUARED|ELECTROSTATIC POTENTIAL|TERMINATED|^\s+[0-9]+\s+[0-9]+\s+-[0-9]" logs/JOB.log | tail -80

ps -eo pid,ppid,pcpu,pmem,etime,cmd | grep -E "gamess|ddikick|rungms|TJ2unlNH" | grep -v grep
```

如果电脑明显卡顿，先查是不是同时有多个 GAMESS：

```bash
tmux ls
ps -eo pid,ppid,pcpu,pmem,etime,cmd | grep gamess.00.x | grep -v grep
```

不要同时开两条大 GAMESS。确认要停时再用：

```bash
pkill -TERM -f JOB_STEM
```

## 16.7 从 large_mk 到 MCPB.py step 3 的最低验收线

最终可用的 M2 large_mk 至少应满足：

```text
输入是 MULT=2。
SCFTYP 是 ROHF 或 UHF，不是默认 RHF。
没有 SCF IS UNCONVERGED。
有 DENSITY CONVERGED，或 DIIS CONVERGED 且 density change 很小。
S-SQUARED 接近 M2 合理范围；doublet 理想值约 0.75。
log 里有 ELECTROSTATIC POTENTIAL / PDC 点。
GAMESS terminated normally。
```

然后再跑：

```bash
source /mnt/l/WSL/conda_envs/AmberTools25/amber.sh
MCPB.py -i mcpb.in -s 3 --logf logs/GROUP_large_mk.log
MCPB.py -i mcpb.in -s 4
```

`--logf` 不是官方唯一写法，但在我们这种有很多尝试文件的目录里很重要。它明确告诉 MCPB.py 解析哪一个 log，避免自动匹配到旧的、失败的、或自旋不一致的文件。

---

一句话记忆：MCPB.py 不是“自动理解金属酶”的黑盒；它是把你定义好的金属配位化学，用 QM 的 Hessian/电荷转换成 Amber 能跑的 bonded force field。你给错 His、H、charge、MULT 或 donor，它会认真地把错误参数化。

# 17. M2 large_mk 现场攻略：70 步上限与保守收尾

这一节专门写给“我不在的时候你自己判断”的场景。核心思想不是死背某个脚本，而是把 GAMESS log 当成仪表盘：看趋势，保存好窗口，别让坏路径无限烧 CPU。

## 17.1 先认清当前 TJ M2 的状态

当前我们要的最终参数是：

```text
small OPT/Hessian: M2，已成功。
large_mk/RESP: 也必须是 M2。
M6 large_mk: 可以当方法参考，不能作为最终 M2 电荷。
```

现在最重要的经验是：

```text
M2-only ROHF/HF checkpoint33 是一个比 checkpoint70 更好的救援起点。
原因不是 33 这个数字神奇，而是 checkpoint33 的 residual/density 状态更干净。
checkpoint70 虽然也是正常 MAXIT 结束并写出了 $VEC，但后半段已经在震荡，步数更大不等于波函数更好。
```

对 SCF 来说，“最后一步”不一定最好；要看 density change、DIIS error、能量是否平滑，以及有没有突然几十 Hartree 的大跳。

## 17.2 每次监控只看这几类行

```bash
cd /mnt/e/TJ/260706_mcpb_m2_200ns/succeed_MCPBpy_project

log=logs/JOB.log

grep -E "DENSITY CONVERGED|DIIS CONVERGED|SCF IS UNCONVERGED|FINAL .*ENERGY|S-SQUARED|ELECTROSTATIC POTENTIAL|TERMINATED|NEXT ENERGY RISES|^\s+[0-9]+\s+[0-9]+\s+-[0-9]" "$log" | tail -120
```

SCF iteration 表里最关键的是：

```text
ITER: 当前 SCF 迭代步。
TOTAL ENERGY: 总能。
E CHANGE: 能量变化。
DENSITY CHANGE: 密度矩阵变化，最终要很小。
DIIS ERROR: 轨道/残差错误，应该总体下降。
VIR. SHIFT / DAMPING: 如果突然巨大，说明程序在强行稳定电子态。
```

粗略判断：

```text
好趋势：
energy 小幅下降或小幅摆动；
density change 从 1、0.1、0.01 往更小走；
DIIS error 进入 1e-3 甚至更低；
没有反复 NEXT ENERGY RISES。

坏趋势：
energy 一步跳 10-30 Hartree；
density change 回到 0.5-1；
DIIS error 卡在 0.03、0.04、0.2 这种平台；
level shift/damping 突然巨大，例如 64、256、4096；
出现多次 NEXT ENERGY RISES 后没有恢复。
```

## 17.3 70 步规则

对于这种 open-shell Fe large_mk，不要让一个坏 SCF 无限跑。当前规则：

```text
最多看到 ITER=70。
如果 70 步前 DENSITY CONVERGED/DIIS CONVERGED，就继续后处理。
如果到 70 还在震荡，就停掉这条路径。
```

停止当前 job：

```bash
stem=TJ2unlNH_large_mk_M2_UB3LYP_from_checkpoint33_dupvec_fullgrid
pkill -TERM -f "$stem"
```

注意：

```text
手动 pkill 后，这条 job 的 .dat 常常只有 14 KB 左右，通常没有完整 $VEC。
不要拿这种被杀掉的 .dat 做 MOREAD。
要用正常 MAXIT 结束过、文件大小约 17 MB 或更大的 checkpoint .dat。
```

## 17.4 如果看到“好窗口”，不要硬等到坏掉

当前 U-B3LYP duplicated-M2-vector 路径里，第 20-21 步一度很好：

```text
ITER 20: density change 约 0.0075，DIIS error 约 0.0032。
ITER 21: 也还可接受。
之后 level shift/DIIS 重置把体系推走。
```

这种情况的正确思路是：

```text
1. 记录最好窗口的 ITER，例如 20 或 21。
2. 复制同一个输入，把 MAXIT 改成这个 ITER。
3. 让 GAMESS 正常跑到 MAXIT 后退出，得到完整 $VEC。
4. 从这个 $VEC 再开 refine，而不是从已经崩掉的第 70/100 步继续。
```

这比“越跑越久总会好”可靠得多。

## 17.5 checkpoint33 -> ROHF/HF 保守收尾

如果 U-B3LYP 到 70 步仍然震荡，就回到 M2-only ROHF/HF checkpoint33。推荐的保守收尾输入骨架是：

```text
$SYSTEM MEMDDI=400 MWORDS=250 KDIAG=3 $END
$CONTRL
 SCFTYP=ROHF RUNTYP=ENERGY EXETYP=RUN
 ICHARG=1 MULT=2 MAXIT=70
$END
$SCF
 DIRSCF=.T. DIIS=.T. SOSCF=.F.
 DAMP=.T. SHIFT=.T. ETHRSH=10.0 MAXDII=3
$END
$BASIS GBASIS=N31 NGAUSS=6 NDFUNC=1 $END
$GUESS GUESS=MOREAD NORB=1048 $END
```

这里的含义：

```text
ROHF/HF: 先不做 DFT 和 ESP，只求一个更干净的 M2 波函数。
MAXIT=70: 正常退出并尽量写完整 $VEC，不无限占资源。
MAXDII=3: 缩小 DIIS 记忆，避免旧错误方向拖着体系大跳。
SOSCF=.F.: 这个体系里 SOSCF 曾把坏 vector 推得更糟，先不用它。
DAMP/SHIFT=.T.: 宁可慢，也尽量减少乱跳。
KDIAG=3: 避开我们见过的对角化问题。
```

如果这条仍然出现巨大能量跳变，可以再试一条更慢的 damping-only 探路：

```text
DIIS=.F. SOSCF=.F. DAMP=.T. SHIFT=.T. MAXIT=40
```

它不一定最终收敛，但可以判断体系是否能沿着同一个电子态慢慢走。

## 17.6 ROHF $VEC 转 UHF/U-B3LYP 的规则

ROHF 的 `$VEC` 通常只有一套轨道；UHF 需要 alpha/beta 两套。如果直接把 ROHF `$VEC` 给 UHF，可能报：

```text
PREMATURE END OF ORBITAL INPUT ENCOUNTERED
```

急救办法：

```text
把 ROHF $VEC 轨道块复制两遍，放进同一个 $VEC ... $END。
第一遍作为 alpha 初猜，第二遍作为 beta 初猜。
然后先跑 EXETYP=CHECK。
CHECK 通过后再正式跑。
```

这一步的目的不是制造物理结果，而是给 U-B3LYP 一个 M2-only 的干净起点，避免从 M6 高自旋历史继承错误电子态。

## 17.7 独立操作的三条路线

遇到 large_mk 不收敛时，按这个顺序选，不要同时开多条大 job。

```text
路线 A：趋势还在变好
继续观察到 40-70 步，但设置上限。

路线 B：中间出现好窗口，后来崩掉
用 MAXIT=N 重跑到好窗口，正常保存 $VEC，再从这个 $VEC refine。

路线 C：U-B3LYP 一直找不到 M2 态
回到 M2-only ROHF/HF checkpoint33，先做保守 refine，得到更稳的 M2 $VEC。
```

最终进入 MCPB.py step 3 前，large_mk 必须满足：

```text
MULT=2。
SCFTYP=UHF/ROHF，不是 RHF。
没有 SCF IS UNCONVERGED。
有 ESP/PDC 输出。
S-SQUARED 没有跑到高自旋态。
GAMESS 正常结束。
```

如果只是 MCPB.py 接受 log，但 log 里有 `SCF IS UNCONVERGED`、能量为 0、或 S-SQUARED 明显不对，不算真正成功。

## 17.8 资源安全

常用检查：

```bash
tmux ls
ps -eo pid,ppid,pcpu,pmem,etime,cmd | grep -E "gamess.00.x|ddikick|rungms|TJ2unlNH" | grep -v grep
```

原则：

```text
大 GAMESS 一次只开一条。
large_mk 默认 8 核，不用 16/17 核硬冲。
不要依赖被 pkill 的 .dat。
真正有价值的是正常 MAXIT 结束后写出的完整 $VEC。
```

## 17.9 2026-07-07 的新增教训：UHF/B3LYP checkpoint20 不是最终救命药

我们从 M2 ROHF/HF checkpoint33 duplicated `$VEC` 出发，重跑了一个 U-B3LYP `MAXIT=20` checkpoint，目的是保存原始长跑里第 20 步附近的好窗口。

结果：

```text
checkpoint20 正常结束，restart .dat 约 35 MB，$VEC 完整。
但它是 MAXIT=20 强制结束，不是 SCF 收敛。
最终 S-SQUARED = 2.616，明显高于 doublet 理想值 0.75。
从这个 checkpoint20 继续 U-B3LYP refine，前 7 步迅速出现大能量跳变和 density 爆炸。
```

结论：

```text
这个 checkpoint20 可以作为诊断证据，说明 UHF/B3LYP 在 large model 上有严重自旋/占据污染；
它不能作为最终 large_mk 电荷来源；
也不应继续沿 UHF/B3LYP checkpoint20 路线硬拧。
```

下一步更合理的方向：

```text
优先考虑 restricted-open-shell B3LYP，也就是 SCFTYP=ROHF DFTTYP=B3LYP，并从 M2 ROHF checkpoint33 的 $VEC 读入。
如果 RO-B3LYP 仍不能稳定，应考虑缩小/重定义 large model，或把电荷拟合模型分块处理。
不要把 M6 large_mk 电荷当作最终 M2 参数。
```

## 17.10 2026-07-07 的 RO-B3LYP 尝试结论

我们按更保守路线试了三条：

```text
1. RO-B3LYP + checkpoint33 ROHF $VEC + short DIIS
   第 3-5 步出现 10-950 Hartree 级大跳，失败。

2. RO-B3LYP + checkpoint33 ROHF $VEC + damping-only
   不再爆炸，但卡在 density 约 0.2-0.7、DIIS error 约 0.04-0.14 的平台。
   说明关 DIIS 能避免失控，但不能自然收敛。

3. 保存 damping-only 的 checkpoint12，再 short-DIIS refine
   checkpoint12 正常写出 18 MB $VEC，S-SQUARED=0.750；
   但 refine 第 1 步 density=144.77，立即失败。
```

结论：

```text
当前 full large model 的 M2 large_mk 不是简单 SCF 参数能救的问题。
UHF/B3LYP 有自旋污染；
RO-B3LYP 可以保持 S^2=0.750，但仍找不到稳定收敛路径。
下一步优先考虑缩小/重定义 large model，或接受 M6 large_mk 作为临时预览参数并明确标注 mixed-spin provisional。
```

## 17.11 cut_off 2.5 小 large model 测试流程

目的：

```text
原成功项目不动；
新建一个 M2 / cut_off 2.5 的 MCPB.py -s1 测试项目；
先检查 Fe donor 和 large model 原子数；
只有当 large model 真的变小且 donor 完整时，才值得继续跑新的 M2 large_mk。
```

可复现建项命令：

```bash
src=/mnt/e/TJ/260706_mcpb_m2_200ns/succeed_MCPBpy_project
dst=/mnt/e/TJ/260706_mcpb_m2_200ns/succeed_MCPBpy_project_cutoff25

mkdir -p "$dst"

cp -p "$src"/mcpb_original.pdb "$dst"/
cp -p "$src"/FE.mol2 "$src"/ACT.mol2 "$src"/UNL.mol2 "$src"/HOH.mol2 "$dst"/
cp -p "$src"/ACT.frcmod "$src"/UNL.frcmod "$dst"/

cp -p "$src"/mcpb_M2_expected.in "$dst"/mcpb.in
sed -i 's/^cut_off .*/cut_off 2.5/' "$dst"/mcpb.in
```

确认 `mcpb.in` 里必须是：

```text
cut_off 2.5
smmodel_chg 1
smmodel_spin 2
lgmodel_chg 1
lgmodel_spin 2
```

运行 step 1：

```bash
source /mnt/l/WSL/conda_envs/AmberTools25/amber.sh
cd /mnt/e/TJ/260706_mcpb_m2_200ns/succeed_MCPBpy_project_cutoff25
MCPB.py -i mcpb.in -s 1 > mcpb_s1_cutoff25.log 2>&1
tail -80 mcpb_s1_cutoff25.log
```

检查 donor：

```bash
grep -E "Selected Metal|is in 2.5|large model contains|Totally there are" mcpb_s1_cutoff25.log
```

2026-07-07 实际结果：

```text
Fe donor 完整：
UNL@N1
HID187@NE2
HID270@NE2
GLU349@OE1
ACT501@O2
HOH875@O

small model: 77 atoms / 303 electrons
standard model: 82 atoms
large model: 109 atoms / 473 electrons
large residues:
[1, 186, 187, 188, 269, 270, 271, 348, 349, 350, 431, 501, 875]
```

判断：

```text
cut_off 2.5 没有实际缩小当前 large model；
它仍然生成 109 atom large_mk，和当前已知困难模型一样大；
因此不建议直接继续跑新的 M2 large_mk，因为大概率只是重复消耗 CPU。
```

下一步如果还要走“缩小 large model”路线：

```text
1. 继续只跑 MCPB.py -s1 测试更小 cut_off，例如 2.3 或 2.2；
2. 每次都先检查 donor 是否完整；
3. 如果 donor 开始丢失，就不能再靠 cut_off 缩小；
4. 如果 cut_off 无法缩小但 donor 又必须保留，就需要手工重定义 large model/cap，而不是继续改 SCF。
```

## 17.12 cut_off 2.4 / 2.3 / 2.2 的实际结果

重要概念：

```text
mcpb.in 里的 cut_off 不是只针对 large model。
它先用于识别金属附近 donor / metal-site residues；
small / standard / large model 都会受这个识别结果影响。
large model 还会额外包含相邻残基并转成 ACE/NME cap。
因此 cut_off 变小不一定立刻缩小 large model；只要 donor 集合没变，large model 可能完全不变。
```

2026-07-07 继续测试：

```bash
bash /mnt/e/Codex/setup_mcpb_cutoff_project.sh 2.4
bash /mnt/e/Codex/setup_mcpb_cutoff_project.sh 2.3
bash /mnt/e/Codex/setup_mcpb_cutoff_project.sh 2.2

source /mnt/l/WSL/conda_envs/AmberTools25/amber.sh
cd /mnt/e/TJ/260706_mcpb_m2_200ns/succeed_MCPBpy_project_cutoff24
MCPB.py -i mcpb.in -s 1 > mcpb_s1_cutoff24.log 2>&1

cd /mnt/e/TJ/260706_mcpb_m2_200ns/succeed_MCPBpy_project_cutoff23
MCPB.py -i mcpb.in -s 1 > mcpb_s1_cutoff23.log 2>&1

cd /mnt/e/TJ/260706_mcpb_m2_200ns/succeed_MCPBpy_project_cutoff22
MCPB.py -i mcpb.in -s 1 > mcpb_s1_cutoff22.log 2>&1
```

结果：

```text
cut_off 2.4:
donor 完整；small model 77 atoms；large model 109 atoms。

cut_off 2.3:
donor 完整；small model 77 atoms；large model 109 atoms。

cut_off 2.2:
UNL@N1 和 HID187@NE2 丢失；
small model 39 atoms；large model 60 atoms。
这个模型不再是正确六配位 Fe active site，不能用于最终参数化。
```

关键 Fe-donor 距离：

```text
GLU349@OE1   2.092 A
ACT501@O2    2.102 A
HOH875@O     2.130 A
HID270@NE2   2.172 A
UNL1@N1      2.223 A
HID187@NE2   2.226 A
```

判断：

```text
cut_off-only 缩小路线卡在临界点：
2.3 仍然不缩小；
2.2 已经丢关键 donor。

因此不能靠继续降低 cut_off 得到“正确但更小”的 large model。
如果还要缩小，只能手工重定义 large model/cap，或者接受当前 109-atom large model 的 M2 收敛困难。
```

## 17.13 手动缩小 large model：不要直接删 ACE/NME cap 当最终路线

背景：
```text
cut_off-only 路线已经测试过：
2.4 / 2.3 不缩小，large model 仍是 109 atoms；
2.2 会丢 UNL@N1 和 HID187@NE2，active site 不再正确。

所以如果 M2 large_mk 继续因为 large model 太大/太难收敛而失败，
下一条现实路线不是继续调 cut_off，而是手动重定义 large model/cap。
```

关键修正：
```text
不要直接删 large_mk 里的 ACE/NME cap 当最终路线。

原因是 MCPB.py 的 small Hessian 和 large_mk 不是两个互不相干的文件：
small Hessian 负责 Fe 周围 bonded force constants；
large_mk 负责 RESP/ESP 电荷；
二者最后都要通过 fingerprint / atom mapping 回到同一个 MCPB 参数体系。

ACE/NME cap 即使不一定作为最终 MD 里的真实残基保留，
也会参与 QM 模型边界、ESP 场和电荷约束。
只从 large_mk 中手动删掉它们，可能造成：
1. 电荷拟合边界变形；
2. atom mapping 或 charge redistribution 不一致；
3. MCPB.py -s3/-s4 表面跑通，但得到的局部电荷不可信；
4. 最终 MD 能运行，却把错误藏进 Fe 附近 electrostatics。
```

已经筛出的远端 cap 信息只能作为诊断参考，不作为默认删除方案：
```text
最远 ACE group: 269-ACE，一共 5 atoms
large_mk 原子序号：50, 51, 52, 53, 54
标签：269-ACE-HH31, 269-ACE-CH3, 269-ACE-C, 269-ACE-O, 269-ACE-HH32
Fe 距离约：8.45-9.52 A

第二远 ACE group: 348-ACE，一共 5 atoms
large_mk 原子序号：76, 77, 78, 79, 80
标签：348-ACE-HH31, 348-ACE-CH3, 348-ACE-C, 348-ACE-O, 348-ACE-HH32
Fe 距离约：7.44-9.03 A

第三远 ACE group: 186-ACE，一共 5 atoms
large_mk 原子序号：24, 25, 26, 27, 28
标签：186-ACE-HH31, 186-ACE-CH3, 186-ACE-C, 186-ACE-O, 186-ACE-HH32
Fe 距离约：7.68-8.66 A
```

不推荐的删法：
```text
不要单独删除最远 H，例如 350-NME-HH32 或 269-ACE-HH31。
这不是合理化学模型，只是制造断键/电子数问题。

不要只在 large_mk.inp 里删 ACE/NME，而不重新生成 MCPB.py 的相关 fingerprint/PDB/mapping。
这会让 Hessian-derived bonded parameters 和 RESP charge model 不再严格自洽。
```

操作建议：
```text
1. 原成功项目不要动。
2. 如果要缩小，必须新建独立测试目录。
3. 不要只编辑 GAMESS large_mk.inp；要让 PDB/fingerprint/large_mk/mapping 自洽。
4. 先跑 MCPB.py -s1 或等价的模型生成步骤，再跑 EXETYP=CHECK。
5. 再跑短 SCF probe。
6. 最后必须验证 MCPB.py -s3 和 -s4 输出：
   - 是否所有最终 MD 原子都有 charge/type；
   - 金属中心总电荷是否合理；
   - donor residue charges 是否没有异常跳变；
   - tleap/gromacs topology 是否没有缺参数；
   - Fe-donor 几何在 EM 后没有崩。
```

结论：
```text
手动删 ACE/NME 可以作为“为什么 SCF 难”的诊断思路，
但不应作为最终 M2 参数化默认路线。

最终可接受路线必须是：
M2 Hessian bonded parameters 与 M2 large_mk RESP charges 来自同一个自洽 MCPB 模型定义。

如果 109-atom large_mk M2 仍不能收敛，下一步优先考虑：
1. 找 MCPB.py 是否能使用 standard/smaller charge model；
2. 重新定义 large model/cap 并完整再生 mapping；
3. 作为临时预览才使用 M6 charge，并明确标注 mixed-spin provisional。
```

## 17.14 关于“删远端主链原子”的判断

用户提出的合理担心：
```text
large model 比 small model 多出的很多原子只是远端主链/边界原子；
它们离 Fe 比较远，而且使 M2 large_mk 很难收敛；
能不能删 5-10 个这种原子？
```

结论：
```text
不能裸删 donor residue 的主链重原子，例如 N、CA、C、O。
这些原子虽然离 Fe 远，但它们决定配位侧链的价态边界和 ESP 环境。
裸删会造成断键、悬空价态、错误电子数，甚至让 M2 更难收敛。

如果要缩小，正确操作不是 delete，而是 truncate + recap：
按 small model 的方式把远端主链替换成 H/CH3 等封端原子，
并重新生成 GAMESS input、fingerprint、charge/multiplicity 和 MCPB mapping。
```

当前 large model 中最“诱人但不能裸删”的远端主链原子：
```text
270-HID-HA   7.95 A from Fe
187-HID-HA   7.83 A
270-HID-N    7.47 A
270-HID-H    7.26 A
349-GLU-O    7.13 A
187-HID-N    6.99 A
270-HID-CA   6.97 A
349-GLU-C    6.88 A
187-HID-CA   6.82 A
```

这些原子可以作为“重封端边界”的参考，
但不作为单独删除清单。

更稳的缩小路线：
```text
1. 优先检查 MCPB.py 已生成的 standard model：
   TJ2unlNH_standard.fingerprint 约 82 atom lines；
   它保留 UNL/HID187/HID270/GLU349/ACT/HOH 和真实 donor residues，
   去掉了 large model 的 ACE/NME 扩边界。

2. 如果要用 standard/small-like 模型做 RESP，
   不能只改 large_mk.inp；
   必须保证 PDB/fingerprint/GAMESS log/MCPB.py -s3 解析对象一致。

3. 如果构造 reduced RESP model：
   - 保留 Fe、全部 donor 原子、UNL、ACT、HOH；
   - 对 HID/GLU 只在远端主链边界 truncate；
   - 所有断开的价态必须补 H/CH3；
   - 重新计算 total charge 与 MULT=2 的电子数奇偶性；
   - 先 EXETYP=CHECK，再短 SCF probe，再完整 MK。
```

## 17.15 visual_check 里的 ACE/NME 不是最终残基，但也是 QM cap

容易误解的地方：
```text
TJ2unlNH_large_visual_check.pdb 主要用于可视化检查；
文件头已经注明：
CONECT FROM LIGAND MOL2 AND METAL-DONOR DISTANCE ONLY

也就是说，它没有把普通蛋白主链肽键 CONECT 都写出来。
因此在 PyMOL/VMD 里，ACE/NME 可能看起来像孤立的、无关的邻近残基。
```

实际几何关系：
```text
ACE186:C  -- HID187:N    1.33 A
HID187:C  -- NME188:N    1.33 A
ACE269:C  -- HID270:N    1.33 A
HID270:C  -- NME271:N    1.33 A
ACE348:C  -- GLU349:N    1.32 A
GLU349:C  -- NME350:N    1.32 A
```

判断：
```text
这些距离是标准酰胺/肽键距离。
所以 ACE/NME 不是最终 MD 里真正要保留的普通残基，
但它们是 large_mk QM 模型里用于封端 His/Glu donor residue 的 cap。

它们可以作为 reduced charge model 的缩小对象，
但不能因为 visual_check 里没有 CONECT 就直接裸删。
裸删后必须重新检查：
1. 总电荷；
2. MULT=2 的电子数奇偶性；
3. GAMESS EXETYP=CHECK；
4. MCPB.py -s3/-s4 是否还能自洽映射；
5. 最终 donor residue charge 是否异常。
```

实用优先级：
```text
如果只是为了少 5-10 个原子做 M2 large_mk 测试，
优先测试“去掉一组或两组 ACE cap”的 diagnostic model，
因为 ACE group 总核电荷为 22，删除后保持 doublet 电子数奇偶性更容易。

不优先单独删 NME group：
NME group 总核电荷为 15；
单删一个 NME 会改变电子数奇偶性，和 MULT=2 冲突，除非同时调整总电荷或配对删除另一个奇数片段。
```

## 17.16 ACE269-only trimmed large_mk 实测结果

2026-07-07 测试目录：
```text
/mnt/e/TJ/260706_mcpb_m2_200ns/succeed_MCPBpy_project_manualtrim_ACE269
```

模型改动：
```text
从 109-atom large model 中删除 269-ACE：
large_mk 原子序号 50-54
269-ACE-HH31, 269-ACE-CH3, 269-ACE-C, 269-ACE-O, 269-ACE-HH32

trim 后：
104 atoms
451 electrons
ICHARG=1
MULT=2
```

通过的检查：
```text
1. TJ2unlNH_large.fingerprint 从 109 行变为 104 行。
2. TJ2unlNH_large.pdb 中 ACE A 269 已删除。
3. GAMESS $DATA 为 104 atoms。
4. RO-B3LYP/MK EXETYP=CHECK 通过。
5. ROHF/HF EXETYP=CHECK 通过。
```

实际 SCF 结果：
```text
直接 RO-B3LYP/MK from HUCKEL：
第 1-9 步 energy / density 大幅震荡；
出现和原 109-atom large model 类似的坏趋势；
主动停止，避免烧到 MAXIT=200。

ROHF/HF preconv, MAXIT=70：
没有 density converged；
日志出现 SCF IS UNCONVERGED, TOO MANY ITERATIONS；
中间有较好点，例如 29、40-45 附近，但后续继续回弹；
最终第 70 步也不是好波函数。

ROHF/HF checkpoint29：
能正常写出 .dat/$VEC；
但 MAXIT=29 重跑轨迹不完全等同于 70 步 job 的第 29 行；
该 checkpoint 只能作为诊断初猜，不能视为严格收敛 HF。

RO-B3LYP/MK from HFcheckpoint29:
MOREAD CHECK 通过；
正式计算读入 $VEC 成功；
但第 9-11 步 density 又反弹到 7-10；
主动停止。
```

结论：
```text
只删除 269-ACE 的 104-atom M2 large_mk 没有解决核心 SCF 问题。
这说明当前困难不是单纯来自最远 5 个 cap 原子；
更像是 M2 large charge model 的 open-shell 占据/电子态问题。

不要把 ACE269-only log 用于 MCPB.py -s3。
它没有给出干净的 M2 B3LYP ESP。
```

下一步选择：
```text
如果坚持 M2 final charges：
1. 试更激进但仍保持偶电子删除的 269-ACE + 348-ACE，减少 10 atoms；
2. 或构造 MCPB.py standard/small-like charge model，而不是继续 small tweaks；
3. 或尝试更强 occupation control / 更低级别预收敛，但成本会继续上升。

如果目标是尽快进入 GROMACS 做结构预览：
使用已经成功的 M2 Hessian bonded parameters + M6 large_mk RESP charge，
明确标注 mixed-spin provisional，只作为结构/流程预览，不作为最终发表级参数。
```

## 17.17 ACE269 + ACE348 双 ACE trimmed large_mk 实测结果

2026-07-07 测试目录：
```text
/mnt/e/TJ/260706_mcpb_m2_200ns/succeed_MCPBpy_project_manualtrim_ACE269_ACE348
```

模型改动：
```text
从原始 109-atom large model 中删除两组 ACE cap：

269-ACE:
large_mk 原子序号 50-54
269-ACE-HH31, 269-ACE-CH3, 269-ACE-C, 269-ACE-O, 269-ACE-HH32

348-ACE:
large_mk 原子序号 76-80
348-ACE-HH31, 348-ACE-CH3, 348-ACE-C, 348-ACE-O, 348-ACE-HH32

trim 后：
99 atoms
429 electrons
ICHARG=1
MULT=2
```

通过的检查：
```text
1. TJ2unlNH_large.fingerprint 从 109 行变为 99 行。
2. TJ2unlNH_large.pdb 中 ACE A 269 和 ACE A 348 已删除。
3. GAMESS $DATA 为 99 atoms。
4. 269-ACE + 348-ACE 总核电荷为 44，删除后仍保持 M2 doublet 所需的奇电子数。
5. RO-B3LYP/MK EXETYP=CHECK 通过。
```

直接 RO-B3LYP/MK from HUCKEL 实测：
```text
第 1-10 步已经明显失控：

ITER 1   density change 5.68
ITER 2   density change 13.05
ITER 3   density change 17.36
ITER 7   density change 104.62
ITER 8   energy jump +2149.89 Hartree, DIIS error 1.97
ITER 10  density change 75.66

主动 Ctrl-C 停止，避免继续浪费 CPU。
```

结论：
```text
删除 269-ACE + 348-ACE 后，模型尺寸从 109 atoms 降到 99 atoms，
但是 M2 direct RO-B3LYP/MK 仍然不能稳定收敛，甚至比 ACE269-only 更早失控。

这进一步说明当前问题不是单纯由 5-10 个远端 ACE cap 原子造成；
核心困难仍是 M2 large charge model 的 open-shell 占据/SCF 稳定性。

不要把该 log 用于 MCPB.py -s3。
```

下一步判断：
```text
继续裸删 ACE cap 的收益已经很低。
如果还要追 M2 final charges，应该转向：
1. chemically recapped reduced model，而不是裸删 cap；
2. MCPB.py standard/small-like charge model；
3. 更明确的占据控制/片段初猜；
4. 或者承认 M2 large RESP 太难，先用 mixed-spin provisional 进入 GROMACS 做结构流程预览。
```

## 17.18 small_model M2 MK/RESP 救援路线实测

2026-07-07 测试目录：
```text
/mnt/e/TJ/260706_mcpb_m2_200ns/succeed_MCPBpy_project_small_mk_M2
```

核心问题：
```text
如果 M2 small model 的 OPT/Hessian 都已经干净成功，为什么不直接用 small_model 做 MK/RESP 电荷？
```

实测结论：
```text
可以算，而且算得非常顺。

输入来源：
TJ2unlNH_M2_hessian_from_punchhess_opt12_diis_shift.inp

计算：
TJ2unlNH_M2_small_mk_from_opt12_UHF_B3LYP_fullgrid

结果：
DENSITY CONVERGED
FINAL U-B3LYP ENERGY = -3001.8042149377
SCF iterations = 4
S-SQUARED = 0.778
MK fitting points = 2776
GAMESS terminated normally
```

判断：
```text
这证明“M2 本身不能算”不是事实。
真正困难来自 large_model 的开壳层占据、边界原子、模型尺寸和 ESP 拟合模型定义。
M2 small_model 的波函数与成功 OPT/Hessian 完全一致，是一个干净的 M2 电荷来源候选。
```

生成 small-MK 输入的脚本：
```text
/mnt/e/Codex/make_small_mk_from_hessian.sh
```

运行示例：
```bash
cd /mnt/e/TJ/260706_mcpb_m2_200ns/succeed_MCPBpy_project_small_mk_M2

/home/qin/softwares/gamess/rungms \
  TJ2unlNH_M2_small_mk_from_opt12_UHF_B3LYP_fullgrid 00 8 \
  > gamess_logs/TJ2unlNH_M2_small_mk_from_opt12_UHF_B3LYP_fullgrid_8core.log 2>&1
```

检查命令：
```bash
grep -aE "DENSITY CONVERGED|SCF IS UNCONVERGED|FINAL U-B3LYP ENERGY|S-SQUARED|NUMBER OF POINTS SELECTED|EXECUTION OF GAMESS" \
  gamess_logs/TJ2unlNH_M2_small_mk_from_opt12_UHF_B3LYP_fullgrid_8core.log
```

## 17.19 small_model RESP 不能直接等同于 MCPB.py -s3

MCPB.py 默认 `-s3` 的逻辑：
```text
1. 用 large.pdb 生成 RESP 输入。
2. 从 large_mk.log 提取 ESP。
3. 用 large.fingerprint 把 resp2.chg 一一对应到 large model。
4. 再筛出 standard.fingerprint 里的 final atoms，生成最终 mol2。
```

因此：
```text
small_model log 不能直接硬塞给普通 MCPB.py -s3。
原因不是 QM 失败，而是 atom number/order/fingerprint 不匹配。
```

当前模型的原子数关系：
```text
small model    77 atoms
standard model 82 atoms
large model   109 atoms

standard 与 small 重叠：65 atoms
standard 中 small 缺失：17 atoms
small 中额外 cap：12 atoms
```

small 缺失的 standard atoms：
```text
187-HID-N, 187-HID-CA, 187-HID-C, 187-HID-O, 187-HID-H, 187-HID-HA
270-HID-N, 270-HID-CA, 270-HID-C, 270-HID-O, 270-HID-H, 270-HID-HA
349-GLU-N, 349-GLU-CA, 349-GLU-C, 349-GLU-O, 349-GLU-HA
```

small 里额外的 cap atoms：
```text
187-HID-H1, 187-HID-CH3, 187-HID-H2, 187-HID-H3
270-HID-H1, 270-HID-CH3, 270-HID-H2, 270-HID-H3
349-GLU-H1, 349-GLU-CH3, 349-GLU-H2, 349-GLU-H3
```

RESP 映射审计脚本：
```text
/mnt/e/Codex/small_m2_resp_map.py
```

运行示例：
```bash
source /mnt/l/WSL/conda_envs/AmberTools25/amber.sh

python /mnt/e/Codex/small_m2_resp_map.py \
  --small-pdb /mnt/e/TJ/260706_mcpb_m2_200ns/succeed_MCPBpy_project_cutoff23/TJ2unlNH_small.pdb \
  --standard-fp /mnt/e/TJ/260706_mcpb_m2_200ns/succeed_MCPBpy_project_cutoff23/TJ2unlNH_standard.fingerprint \
  --small-mk-log /mnt/e/TJ/260706_mcpb_m2_200ns/succeed_MCPBpy_project_small_mk_M2/gamess_logs/TJ2unlNH_M2_small_mk_from_opt12_UHF_B3LYP_fullgrid_8core.log \
  --outdir /mnt/e/TJ/260706_mcpb_m2_200ns/succeed_MCPBpy_project_small_mk_M2/small_resp_map \
  --mol2 \
  /mnt/e/TJ/260706_mcpb_m2_200ns/succeed_MCPBpy_project_cutoff23/FE.mol2 \
  /mnt/e/TJ/260706_mcpb_m2_200ns/succeed_MCPBpy_project_cutoff23/ACT.mol2 \
  /mnt/e/TJ/260706_mcpb_m2_200ns/succeed_MCPBpy_project_cutoff23/UNL.mol2 \
  /mnt/e/TJ/260706_mcpb_m2_200ns/succeed_MCPBpy_project_cutoff23/HOH.mol2
```

审计结果：
```text
RESP small total charge = 0.999998
overlap RESP charge sum = 0.960724
discarded cap charge sum = 0.039274
ff19SB missing-backbone charge sum = -0.426400

如果 overlap 用 small-M2 RESP，而缺失主链直接用 ff19SB：
standard total before correction = 0.534324
需要补 +0.465676 才能回到总电荷 +1。
```

这说明：
```text
small-M2 RESP 是干净的量子结果；
但从 capped small model 映射到 final standard model 时，必须定义电荷闭合策略。
不能把 0.534324 总电荷的 mol2 直接拿去做最终 MD。
```

两个可选校正：
```text
方案 A：cap-group correction
把每个 donor residue 的 small cap group 总电荷转移给 standard 缺失主链原子组。
优点：严格保留 small capped model 的分组总电荷。
缺点：GLU349 缺失主链平均每原子需要 +0.077735 e，偏大。

方案 B：global-overlap correction
缺失主链保持 ff19SB；把 +0.465676 均匀分摊到 65 个 small-RESP 重叠原子。
每个重叠原子只加 +0.007164 e。
这更接近 MCPB.py 默认“保留 backbone，调整 metal-site RESP atoms”的思想。
```

生成的审计文件：
```text
small_M2_RESP_to_standard_charge_map.tsv
small_M2_RESP_to_standard_charge_map_capgroup_corrected.tsv
small_M2_RESP_to_standard_charge_map_global_overlap_corrected.tsv
small_M2_RESP_mapping_report.txt
```

当前建议：
```text
最终发表级别不要使用 M6 large_mk charge。
如果 M2 large_mk 持续不收敛，可以把 small-M2 RESP + global-overlap correction 作为救援版 M2 charge route。
但在进入最终 MD 前，必须：
1. 用 corrected charge map 生成/替换 mol2 charges；
2. 检查每个 mol2 residue 总电荷；
3. 跑 tleap/parmed 转 GROMACS；
4. 做真空/水中短 minimization，确认 Fe 配位键和电荷没有异常；
5. 在方法部分明确说明 charge model 是 spin-consistent M2 capped-small-model RESP with backbone fixed/corrected normalization。
```

## 17.20 永久提醒：mcpb.in 的 spin 不能靠记忆

`mcpb.in` 中这两行必须每次人工检查：
```text
smmodel_spin
lgmodel_spin
```

本项目曾经出现：
```text
smmodel_spin 6
lgmodel_spin 6
```

这会导致 MCPB.py `-s1` 生成的 small_fc/large_mk 输入继承 M6。
即使后来手动改了 Hessian 为 M2，如果 large_mk 仍来自 M6，就会产生：
```text
M2 Hessian bonded parameters + M6 RESP charges
```

该组合只可用于临时流程预览，不可作为最终发表级参数。

## 17.21 ChgModB/C/D 约束测试

2026-07-07 测试目录：
```text
/mnt/e/TJ/260706_mcpb_m2_200ns/succeed_MCPBpy_project_chgmod_tests_20260707_v2
```

注意：
```text
这次测试使用的是已经正常结束的 M6 large_mk log：
TJ2unlNH_large_mk_uhf_direct_fullgrid_8core.log

所以它不是最终 M2 电荷结果。
它只用于理解 MCPB.py 官方 ChgModB/C/D 对 RESP/backbone 约束的行为。
```

运行脚本：
```text
/mnt/e/Codex/run_chgmod_tests.sh
```

三种模式都成功：
```text
MCPB.py -s 3b  OK
MCPB.py -s 3c  OK
MCPB.py -s 3d  OK

large model:    109 atoms
standard model: 82 atoms
最终 mol2 总电荷均约为 +1
```

各 residue 电荷总和：
```text
3b:
UL1  +0.409649
HD1  +0.329778
HD2  +0.403256
GU1  -0.435577
FE1  +0.128024
AT1  +0.005690
HH1  +0.159179
TOTAL +0.999999

3c:
UL1  +0.402030
HD1  +0.341347
HD2  +0.380885
GU1  -0.435210
FE1  +0.155409
AT1  +0.000166
HH1  +0.155373
TOTAL +1.000000

3d:
UL1  +0.278250
HD1  +0.399730
HD2  +0.322469
GU1  -0.516036
FE1  +0.418729
AT1  -0.060343
HH1  +0.157196
TOTAL +0.999995
```

关键观察：
```text
ChgModB/3b 和 ChgModC/3c 会让 GLU349 的 CB/CG 出现极端 RESP 电荷：

3b:
GU1-CB  +1.766215
GU1-CG  -2.256107

3c:
GU1-CB  +1.495556
GU1-CG  -2.049262

这类正负巨大抵消的电荷不适合直接作为最终参数。
```

ChgModD/3d 的表现：
```text
3d 额外约束 CB，因此避免了 GLU349 CB/CG 被拉到 ±2 e 附近。
3d 中最主要的金属核心电荷：
FE1-FE  +0.418729
UNL-N1  -0.545783
GLU-OE1 -0.595456
HOH-O   -0.657000
```

判断：
```text
如果后续使用官方 large_mk log 做 RESP，当前体系不宜优先用默认 3b；
应优先检查 3d，或者至少明确比较 3b/3c/3d 的逐原子电荷异常。

如果后续构建 reduced-M2 charge model，也应借鉴 3d 思想：
backbone 和 CB 尽量保持 ff19SB/固定，
把 RESP 自由度主要留给 metal donor 原子、Fe、UNL、ACT、HOH 和真正需要重拟合的侧链末端。
```

## 17.22 small-M2 RESP + ChgModD-like 候选参数

目标：
```text
large_model M2 GAMESS SCF 长期不收敛；
但 small_model M2 MK/RESP 已经干净收敛。

因此构建一个更像 ChgModD 的救援电荷模型：
1. backbone 固定；
2. donor residue 的 CB 固定；
3. small-M2 RESP 电荷用于 Fe、donor 原子、UNL、ACT、HOH、侧链末端；
4. 总电荷差额只在自由 RESP 原子上归一化。
```

生成脚本：
```text
/mnt/e/Codex/make_smallM2_chgmod3like_mol2.py
```

候选 mol2 目录：
```text
/mnt/e/TJ/260706_mcpb_m2_200ns/succeed_MCPBpy_project_smallM2_chgmod3like_candidate
```

运行命令：
```bash
python3 /mnt/e/Codex/make_smallM2_chgmod3like_mol2.py \
  --standard-fp /mnt/e/TJ/260706_mcpb_m2_200ns/succeed_MCPBpy_project_cutoff23/TJ2unlNH_standard.fingerprint \
  --small-map /mnt/e/TJ/260706_mcpb_m2_200ns/succeed_MCPBpy_project_small_mk_M2/small_resp_map/small_M2_RESP_to_standard_charge_map.tsv \
  --template-dir /mnt/e/TJ/260706_mcpb_m2_200ns/succeed_MCPBpy_project_chgmod_tests_20260707_v2/3d \
  --outdir /mnt/e/TJ/260706_mcpb_m2_200ns/succeed_MCPBpy_project_smallM2_chgmod3like_candidate
```

电荷归一化结果：
```text
standard atoms: 82
fixed backbone/CB atoms: 20
adjustable small-M2 RESP atoms: 62
total before normalization: 0.523102
normalization correction: +0.476898
per adjustable atom shift: +0.00769190
total after normalization: +1.000000
```

固定的 backbone/CB 原子：
```text
187-HID-N, CA, C, O, CB, H, HA
270-HID-N, CA, C, O, CB, H, HA
349-GLU-N, CA, C, O, CB, HA
```

候选 mol2 的 residue 总电荷：
```text
UL1  +0.779233
HD1  +0.520368
HD2  +0.352779
GU1  -0.699617
FE1  +0.136432
AT1  -0.300960
HH1  +0.211771
TOTAL +1.000006
```

极端电荷检查：
```text
最大正电荷约 +0.666
最大负电荷约 -0.586

没有 ChgModB/C 中 GLU-CB/CG 出现的 +1.5 到 +1.8 / -2.0 到 -2.25 e 异常抵消。
```

MCPB.py -s4 测试目录：
```text
/mnt/e/TJ/260706_mcpb_m2_200ns/succeed_MCPBpy_project_smallM2_chgmod3like_s4
```

`MCPB.py -s4` 结果：
```text
成功生成：
TJ2unlNH_tleap.in
TJ2unlNH_mcpbpy.pdb
```

注意：
```text
自动生成的 TJ2unlNH_tleap.in 使用原始 PDB residue number，例如 mol.431.FE。
tleap 读入后内部 residue index 不等于原始编号，因此会出现 bond 命令解析失败。
这不是电荷问题。

需要使用原成功项目中的 internal-index 修正版：
TJ2unlNH_tleap_internalidx.in
```

使用 internal-index tleap 输入后：
```text
tleap 成功，Errors = 0

生成：
TJ2unlNH_dry.prmtop
TJ2unlNH_dry.inpcrd
TJ2unlNH_solv.prmtop
TJ2unlNH_solv.inpcrd
TJ2unlNH_dry.pdb
TJ2unlNH_solv.pdb
```

当前判断：
```text
这是目前比 M6 mixed-spin RESP 更合理的 M2-consistent 候选路线。
它仍然不是“官方默认 large_mk M2 RESP”，方法部分必须如实说明：
small-model M2 RESP + ChgModD-like fixed backbone/CB normalization。

进入正式 MD 前，下一步应：
1. 转换为 GROMACS top/gro；
2. 做 dry/solv 结构检查；
3. 跑能量最小化；
4. 检查 Fe-donor 距离、体系总电荷、是否有断裂或爆炸。
```
