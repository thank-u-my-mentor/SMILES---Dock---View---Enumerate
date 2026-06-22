# H++ 输入 PDB 准备规则

这份规则用于把晶体结构、对接后的底物、金属离子整理成 H++ 更容易接受的 PDB 文件。目标不是让结构变成最终 MD 文件，而是先让 H++ 正确判断蛋白质残基的质子化状态。

## 基本原则

1. 每次只处理一个chain。
   - 如果只研究 chain A，就删除 chain B。
   - 如果保留多条链，必须确认每条链都有正确的 `TER`。

2. 先删除占据目标口袋的旧晶体 ligand。
   - 例如 `GOL`、`HHH`、原底物、抑制剂、缓冲剂分子。
   - 如果目标是预测新底物进入这个位点，旧 ligand 必须删掉。

3. 金属离子可以保留。
   - Fe、Zn、Mg、Ca 等金属可以保留在 H++ 输入里。
   - 金属附近的 His/Asp/Glu 质子化状态通常要人工复核。
   - 对 Fe 酶尤其要记录 Fe 周围 2.5-3.0 A 内的配位原子。

4. 对接后的目标 ligand 可以保留。
   - H++ 已经证明可以接受“蛋白 + 金属 + docked ligand”的输入。
   - ligand 的 residue name 建议统一，例如 `LIG`。
   - ligand 前后要有 `TER`，不要让它和蛋白链连在一起。

5. 水分子一般先删除。
   - H++ 判断蛋白质 pKa 时通常不需要保留所有晶体水。
   - 对接阶段也建议删除水，因为水可能挡住搜索空间或引入冲突。
   - 只有当某个水是明确催化水、桥连金属水、或文献强支持的保守水时，才单独保留并记录理由。

## H++ 容易报错的常见原因

1. 蛋白链中间缺 residue。
   - H++ 会报告 sequence discontinuity。
   - 如果缺口远离活性中心，可以把该链拆成两个 terminal segment，但要人工接受这个近似。
   - 如果缺口靠近活性中心，建议用 Rosetta/Modeller 补 loop 后再送 H++。

2. 多条链之间没有 `TER`。
   - H++ 可能误以为 chain A 和 chain B 是一条连续链。
   - 每条蛋白链结束后都要有 `TER`。

3. 非标准 residue name 混在蛋白主链里。
   - 标准氨基酸必须用标准三字母名。
   - H++ 输出里的 `HID/HIE/HIP` 后续要转成 Amber 可接受的 residue name。

4. alternate conformation 没处理。
   - 如果同一个 residue 有 A/B 两套构象，建议只保留 occupancy 较高或人工选择的一套。
   - 不要把同一个原子位置的 A/B 构象都交给 H++。

5. ligand 或金属被误连到蛋白链。
   - HETATM 前后加 `TER`。
   - ligand residue name 使用 `LIG`。
   - 金属单独作为 HETATM 保留，不要改成 ATOM。

## 推荐工作流

1. 从原始 PDB 复制一份工作文件，不修改原始文件。
2. 删除目标口袋旧 ligand，例如 `GOL`。
3. 删除普通水 `HOH`。
4. 如只研究单体，删除其他 chain。
5. 保留催化金属 Fe。
6. 加入你选定的 docked pose ligand，并命名为 `LIG`。
7. 确认蛋白、金属、ligand 之间有合理 `TER`。
8. 检查是否有明显断链、重复 residue、altloc。
9. 送 H++。
10. 将 H++ 的 `.pkout` 保存到项目目录，
11. #然后用 `ple hpp-map` 转成 Amber protonation map。此步可交给秦天翼让他独立完成

## 后续 MD 提醒

H++ 只解决蛋白质质子化趋势，不会自动给 ligand 或 Fe 配位中心生成严肃力场参数。进入 MD 前还需要：

- ligand 用 antechamber/GAFF2 或更严肃方法生成参数；
- Fe 中心用 MCPB.py、文献参数，或临时 harmonic restraint 模型处理；
- 检查 His/Asp/Glu 的最终 Amber residue name；
- 检查 ligand 是否有错误键长、错误手性或异常构象。
