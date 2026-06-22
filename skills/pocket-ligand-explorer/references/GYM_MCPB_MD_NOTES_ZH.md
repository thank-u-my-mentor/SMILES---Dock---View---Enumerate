# GYM 4X8B Fe(III)-MCPB 到 40 ns MD 工作记录

## 当前体系

- 原始 PDB：`/mnt/e/GYM/4X8B.pdb`
- docking receptor：删除 `GOL` 和 `HOH`，保留 chain A 的 Fe
- 选定 pose：`S000001_mode1`
- ligand SMILES：`C=C(CCC(ONC(OC)=O)=O)C1=CC=CC=C1`
- H++ pkout：`/mnt/e/GYM/4X8B_drop_GOL_HOH.pkout`
- Amber protonation map：`/mnt/e/GYM/protonation_map_chainA_unique.csv`

## Fe 配位核心

H++/坐标审计识别到 Fe 周围 2.8 A 内的三个配位 His：

- `HID A51 NE2`，Fe-N 距离约 2.136 A
- `HID A134 NE2`，Fe-N 距离约 2.068 A
- `HID A138 NE2`，Fe-N 距离约 2.117 A

这次 MCPB.py 使用 Fe(III) 高自旋假设：

- `software_version gms`
- small model spin multiplicity：6
- large model spin multiplicity：6
- ligand 暂时作为非配位底物，不纳入 Fe bonded model

## 关键输出

- receptor：
  `/mnt/e/GYM/receptor/receptor_chainA_noGOL_noHOH_keepFe_hppAmber.pdb`
- pose1 ligand：
  `/mnt/e/GYM/receptor/ligand_S000001_mode1.pdb`
- MCPB 目录：
  `/mnt/e/GYM/metal_model/S000001_mode1_FeIII_hppAmber`
- GAMESS 运行脚本：
  `/mnt/e/GYM/metal_model/S000001_mode1_FeIII_hppAmber/run_gamess_mcpb_qm.sh`
- PyMOL pose1：
  `/mnt/e/GYM/pymol/S000001_mode1_palegreen_surface.pml`

## 下一步运行顺序

进入 MCPB 目录：

```bash
cd /mnt/e/GYM/metal_model/S000001_mode1_FeIII_hppAmber
```

运行 GAMESS QM：

```bash
bash run_gamess_mcpb_qm.sh
```

这个脚本顺序运行：

1. `FE_A508_6707_small_opt.inp`
2. `FE_A508_6707_small_fc.inp`
3. `FE_A508_6707_large_mk.inp`

完成后检查：

```bash
grep -i "terminated normally" gamess_logs/*.log
```

如果三个任务都正常结束，再继续 MCPB.py 后处理：

```bash
MCPB.py -i mcpb.in -s 2
MCPB.py -i mcpb.in -s 3
MCPB.py -i mcpb.in -s 4
```

成功后会得到 Fe bonded model 的 Amber 参数文件。之后再把这些参数接入
`ple md-run --time 40`，进行 40 ns 常规 MD。

## 注意

不要再使用 harmonic Fe-His 弹簧作为正式模型。harmonic restraint 只能作为
negative/positive control 或临时稳定手段。GYM 的正式路线应以 MCPB.py + GAMESS
生成的 bonded Fe 参数为准。
