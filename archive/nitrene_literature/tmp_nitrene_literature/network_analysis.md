# 引用网络拓扑分析

> **PageRank (PR Score)**：引用网络中的『稳态影响力』。可以理解为：如果一只随机游走的学术兔子在网络中跳来跳去，它最终停留在某篇文献上的概率。PR 越高，说明该文献被越多重要文献引用，是领域的枢纽/综述/里程碑。
> 
> **Betweenness Centrality (BC Score)**：桥梁价值。如果移除该文献，网络中多少最短路径会断裂？BC 高的文献通常是跨子领域的『连接器』，把光化学和酶工程、或合成和方法学连在一起。

## 二、网络核心（PageRank Top 20）

| Rank | DOI | Title | Year | Venue | Citations | PR Score | Seed? |
|------|-----|-------|------|-------|-----------|----------|-------|
| 1 | [10.1021/ja509308v...](https://doi.org/10.1021/ja509308v) | Enzyme-Controlled Nitrogen-Atom Transfer Enables Regiod... | 2014 |  | 0 | 0.0805 | ✅ |
| 2 | [10.1038/s41570-021-0...](https://doi.org/10.1038/s41570-021-00291-4) | Nitrene transfer catalysts for enantioselective C-N bon... | 2021 |  | 0 | 0.0666 | ✅ |
| 3 | [10.1021/acscentsci.5...](https://doi.org/10.1021/acscentsci.5b00056) | Enantioselective Enzyme-Catalyzed Aziridination Enabled... | 2015 |  | 0 | 0.0581 | ✅ |
| 4 | [10.1021/jacs.9b11608...](https://doi.org/10.1021/jacs.9b11608) | Nitrene Transfer Catalyzed by a Non-Heme Iron Enzyme an... | 2019 |  | 0 | 0.0513 | ✅ |
| 5 | [10.1002/cbic.2014022...](https://doi.org/10.1002/cbic.201402286) | Non‐natural Olefin Cyclopropanation Catalyzed by Divers... | 2014 | ChemBioChem | 43 | 0.0289 |  |
| 6 | [10.1002/cbic.2014022...](https://doi.org/10.1002/cbic.201402241) | Structural, Functional, and Spectroscopic Characterizat... | 2014 | ChemBioChem | 46 | 0.0289 |  |
| 7 | [10.1021/ja503593n...](https://doi.org/10.1021/ja503593n) | Enantioselective Imidation of Sulfides via Enzyme-Catal... | 2014 | Journal of the Ameri | 122 | 0.0289 |  |
| 8 | [10.1039/D3SC04661C...](https://doi.org/10.1039/D3SC04661C) | Catalytic, asymmetric carbon-nitrogen bond formation us... | 2023 |  | 0 | 0.0257 | ✅ |
| 9 | [10.1016/j.bmc.2014.0...](https://doi.org/10.1016/j.bmc.2014.05.015) | Intramolecular C(sp3)—H amination of arylsulfonyl azide... | 2014 | Bioorganic & Medicin | 111 | 0.0218 |  |
| 10 | [10.1002/cbic.2014000...](https://doi.org/10.1002/cbic.201400060) | Enhancing the Efficiency and Regioselectivity of P450 O... | 2014 | ChemBioChem | 73 | 0.0218 |  |
| 11 | [10.1038/nchembio.143...](https://doi.org/10.1038/nchembio.1438) | Direct nitration and azidation of aliphatic carbons by ... | 2014 | Nature Chemical Biol | 137 | 0.0218 |  |
| 12 | [10.1021/cs400893n...](https://doi.org/10.1021/cs400893n) | P450-Catalyzed Intramolecular sp3 C–H Amination with Ar... | 2014 | ACS Catalysis | 164 | 0.0218 |  |
| 13 | [10.1002/anie.2019025...](https://doi.org/10.1002/anie.201902576) | Nickel-Catalyzed Asymmetric Hydrogenation of N-Sulfonyl... | 2019 | Angewandte Chemie | 114 | 0.0208 |  |
| 14 | [10.1021/jacs.8b13906...](https://doi.org/10.1021/jacs.8b13906) | Repurposing Nonheme Iron Hydroxylases To Enable Catalyt... | 2019 | Journal of the Ameri | 19 | 0.0208 |  |
| 15 | [10.1021/acs.biochem....](https://doi.org/10.1021/acs.biochem.8b00730) | Thermodynamics of Iron(II) and Substrate Binding to the... | 2018 | Biochemistry | 15 | 0.0208 |  |
| 16 | [10.1146/annurev-bioc...](https://doi.org/10.1146/annurev-biochem-061516-044724) | 2-Oxoglutarate-Dependent Oxygenases. | 2018 | Annual Review of Bio | 408 | 0.0208 |  |
| 17 | [10.1021/jacs.8b02769...](https://doi.org/10.1021/jacs.8b02769) | Unprecedented Cyclization Catalyzed by a Cytochrome P45... | 2018 | Journal of the Ameri | 88 | 0.0208 |  |
| 18 | [10.1002/anie.2014094...](https://doi.org/10.1002/anie.201409470) | Expanding the Enzyme Universe: Accessing Non-Natural Re... | 2015 | Angewandte Chemie | 454 | 0.0191 |  |
| 19 | [10.1021/cs5018612...](https://doi.org/10.1021/cs5018612) | Enzymatic C(sp3)-H Amination: P450-Catalyzed Conversion... | 2015 | ACS Catalysis | 141 | 0.0191 |  |
| 20 | [10.1002/anie.2014028...](https://doi.org/10.1002/anie.201402809) | Improved cyclopropanation activity of histidine-ligated... | 2014 | Angewandte Chemie | 162 | 0.0191 |  |

## 三、桥梁文献（Betweenness Centrality Top 15）

| Rank | DOI | Title | Year | Venue | Citations | BC Score |
|------|-----|-------|------|-------|-----------|----------|
| 1 | [10.1021/jacs.6c04805...](https://doi.org/10.1021/jacs.6c04805) | Photoactive Iminobismuthanes for Catalytic C–H Aminatio... | 2026 | Journal of the Ameri | 0 | 0.0000 |
| 2 | [10.1021/jacs.6c05174...](https://doi.org/10.1021/jacs.6c05174) | In Crystallo Synthesis of a Triplet Silver Nitrene | 2026 | Journal of the Ameri | 0 | 0.0000 |
| 3 | [10.1126/science.aee3...](https://doi.org/10.1126/science.aee3321) | Chiral S(VI) platform unifies selective C–H amination o... | 2026 | Science | 1 | 0.0000 |
| 4 | [10.1002/chem.71018...](https://doi.org/10.1002/chem.71018) | Enantioselective Intramolecular Benzylic C–H Amination ... | 2026 | Chemistry – A Europe | 0 | 0.0000 |
| 5 | [10.1002/chem.70977...](https://doi.org/10.1002/chem.70977) | Nitrene Transfer‐Triggered Novel
 Oxa
 −Pictet−Spengler... | 2026 | Chemistry – A Europe | 0 | 0.0000 |
| 6 | [10.1039/d6sc01440b...](https://doi.org/10.1039/d6sc01440b) | Highly stereoselective synthesis of allylic β-lactams v... | 2026 | Chemical Science | 0 | 0.0000 |
| 7 | [10.1016/j.checat.202...](https://doi.org/10.1016/j.checat.2025.101331) | Enantioselective palladium-catalyzed α-arylation of pri... | 2025 | Chem Catalysis | 2 | 0.0000 |
| 8 | [10.1126/science.adh8...](https://doi.org/10.1126/science.adh8753) | Mechanistic snapshots of rhodium-catalyzed acylnitrene ... | 2023 | Science | 61 | 0.0000 |
| 9 | [10.1021/jacs.3c05258...](https://doi.org/10.1021/jacs.3c05258) | Regio- and Enantioselective Catalytic δ-C-H Amidation o... | 2023 | Journal of the Ameri | 28 | 0.0000 |
| 10 | [10.1021/jacs.3c03587...](https://doi.org/10.1021/jacs.3c03587) | Aspartyl β-Turn-Based Dirhodium(II) Metallopeptides for... | 2023 | Journal of the Ameri | 17 | 0.0000 |
| 11 | [10.1021/acs.orglett....](https://doi.org/10.1021/acs.orglett.3c00940) | Catalytic Enantioselective Amination of Enol Silyl Ethe... | 2023 | Organic Letters | 11 | 0.0000 |
| 12 | [10.1021/jacs.3c00693...](https://doi.org/10.1021/jacs.3c00693) | Substrate-Directed Enantioselective Aziridination of Al... | 2023 | Journal of the Ameri | 30 | 0.0000 |
| 13 | [10.1002/anie.2022185...](https://doi.org/10.1002/anie.202218577) | Chiral Iron Porphyrins Catalyze Enantioselective Intram... | 2023 | Angewandte Chemie | 27 | 0.0000 |
| 14 | [10.1021/jacs.2c08285...](https://doi.org/10.1021/jacs.2c08285) | Enzymatic Nitrogen Insertion into Unactivated C–H Bonds | 2022 | Journal of the Ameri | 56 | 0.0000 |
| 15 | [10.1021/jacs.2c07337...](https://doi.org/10.1021/jacs.2c07337) | Rhodium(II)-Catalyzed Enantioselective Intermolecular A... | 2022 | Journal of the Ameri | 26 | 0.0000 |

## 四、新兴聚类（2022 年以来）

> 未发现足够大的新兴聚类（聚类大小 < 3）