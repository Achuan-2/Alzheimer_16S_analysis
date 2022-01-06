# 阿尔茨海默症患者肠道菌群16S-测序数据分析

Video: https://www.bilibili.com/video/BV1Sm4y1S7Mi
## Methods
![流程图_2022-01-06](https://cdn.jsdelivr.net/gh/Achuan-2/PicBed@pic/assets/README/流程图_2022-01-06.svg)
## Result
171 个样本共得到 1,584 个 ASVs。其中 AD 患者组和 C 对照组共有的 ASVs 有 1126 个，AD 患者组独有的 ASVs 有 330 个，C 对照组独有的 ASVs 有 128 个。
![Group_Venn_plot_2022-01-06](https://cdn.jsdelivr.net/gh/Achuan-2/PicBed@pic/assets/README/Group_Venn_plot_2022-01-06.svg)

对两组样本的数据进行 alpha 多样性分析（图 1A），在群落丰富度指标（包括 Observed species, ACE, 和 Chao1 指数）以及群落多样性指标（包括 Shannon 和 Simpson 指数）上，AD 患者组肠道菌群的物种丰富度以及多样性总体上低于健康对照组，存在显著差异（Wilcoxon 检验，p 值均小于 0.05)。

对两组样本的数据进行 beta 多样性分析，分别计算 Jaccard, Bray–Curtis, Unweighted UniFrac 和 Weighted UniFrac 距离，使用 PCoA 在二维平面投影（图 1B），可知无论是否考虑物种丰度变化和物种进化关系，尽管组内个体存在差异，AD 组和 C 组样本组成群落分布明显分开，组间微生物群落整体结构均存在显著性差异（PERMANOVA 检验，p 值均小于 0.05）。

通过 alpha 多样性和 beta 多样性分析可以得知，阿尔兹海默症会影响到人肠道菌群的整体结构。
![alpha+beta_2022-01-06](https://cdn.jsdelivr.net/gh/Achuan-2/PicBed@pic/assets/README/alpha+beta_2022-01-06.svg)
> 图 1 （A ）两组间 Observed species、 ACE、 Chao1 以及 Shannon、Simpson 指数得分的箱式图：红色代表 AD 患者组，绿色代表健康对照组；箱子的上下限，分别是数据的上四分位数和下四分位数；箱子中间的一条线，是数据的中位数，代表了组内数据的平均水平；每个指数都使用 Wilcoxon 检验分析组间显著性差异，'\*\*' 代表 p 值 < 0.01， '\*' 代表 p 值 < 0.05  （B）使用 PCoA 降维可视化不同距离算法下的样本距离：四张图从左到右，从上到下分别是：Jaccard, Bray–Curtis, Unweighted UniFrac 和 Weighted UniFrac 距离的 PCoA 降维距离分布图。横坐标（PCoA1）表示第一主成分，百分比则表示第一主成分对样品差异的贡献值；纵坐标（PCoA2）表示第二主成分，百分比表示第二主成分对样品差异的贡献值；每一个点代表一个样本，相同颜色的点来自同一个分组，红色代表 AD 患者组，绿色代表健康对照组；虚线椭圆代表的是 90% 置信度下，对该组的总体区间估计；每张图的左上方标注了使用 PERMANOVA 检验得到的 p 值。


### 阿尔兹海默症患者肠道微生物组成的相对失调

在门水平上，发现 AD 患者肠道微生物中 Verrucomicrobiota（疣微菌门） 和 Actinobacteriota （放线菌门）明显增多；而 Firmicutes （厚壁菌门）显著减少。而 Proteobacteria （变形菌门）和 Bacteroidetes（拟杆菌门）没有显著变化，这与一些研究不一致

在科水平上，AD 患者肠道微生物中 Bifidobacteriaceae（双歧杆菌科）、Akkermansiaceae（艾克曼菌科）、Enterococcaceae（肠球菌科）、Corynebacteriaceae（棒状杆菌科）、Eubacteriaceae（优杆菌科） 在内的细菌家族显著增加；而 Ruminococcaceae（瘤胃菌科）、Lachnospiraceae（毛螺菌科）、Peptostreptococcaceae（消化链球菌科）、Monoglobaceae、Butyricicoccaceae 显著下降。在 AD 患者减少的科水平中，Ruminococcaceae 和 Lachnospiraceae 生产不同类型的短链脂肪酸（SCFA）。在 SCFA 中，丁酸盐因其对维持健康的有益作用而在研究中受到特别关注。丁酸盐可影响胃肠生理学、肝脏代谢的外周免疫和血脑屏障的完整性，从而间接促进大脑功能，它可以驱动小胶质细胞的成熟，并且是维持成熟小胶质细胞所必需的。

在 Genus 水平上，Bifidobacterium（双歧杆菌属）、Akkermansia （艾克曼菌属）、[Ruminococcus]_gnavus_group（瘤胃球菌属）、Erysipelatoclostridium、Enterococcus（肠球菌属）、Hungatella 显著增加；Faecalibacterium（粪杆菌属）、Agathobacter、Roseburia（罗斯氏菌属）、Romboutsia、Lachnospiraceae_NK4A136_group（毛螺菌属）、Clostridia_UCG-014、Coprococcus（粪球菌属）、Fusicatenibacter、Butyricicoccus、Lachnospira（毛螺菌属）、[Eubacterium]_ruminantium_group 显著下降。

在 AD 患者增加的属水平菌种中，Enterococcus 可在原代大鼠皮层神经元中产生早期阿尔茨海默样神经原纤维表位 ，在 AD 病因中可作为有害细菌。有趣的是 Bifidobacterium 和 Akkermansia 在传统上被认为是有益菌，却在 AD 患者中显著增加。Bifidobacteriaceae 和 Enterococcaceae 主要是产生乳酸，在 AD 患者中增加。Bifidobacteriaceae 主要是一种乳酸产生菌，对人类非常有益，并已被用作乳制品中的食品添加剂。Kobayashi 等人进行的一项初步临床研究发现，使用 Bifidobacterium breve A1 的口服补充剂可以通过抑制炎症和免疫反应基因的基因表达，改善认知功能，维持老年人的生活质量。Akkermansia 是一个专门降解粘蛋白的属，可利用粘蛋白衍生的糖，如岩藻糖，通过丙二醇途径生产丙酸盐。先前的研究表明，Akkermansia muciniphila（典型菌株）与预防肥胖、促进伤口愈合、增强抗肿瘤反应和诱导免疫反应有关 )。令人惊讶的是，我们的数据表明 Bifidobacterium 和 Akkermansia 是 AD 患者肠道富集微生物群中最丰富的属之一。这表明两者可能在 AD 的发病和发展中起关键作用。


在减少的属水平菌种，Faecalibacterium （典型菌株 F.prausnitzii）是厚壁菌门的主要成员，被认为是健康肠道最重要的细菌指标之一，可调节肠道上皮水平的炎症，有研究也发现老年帕金森病患者 Faecalibacterium 比例下降，Bifidobacterium 增多，先前的研究发现，Faecalibacterium 具有抗炎特性，因为其能够产生丁酸盐并诱导 耐受性细胞因子，这些变化都可能导致促炎性肠道环境，从而导致健康状况下降的老年人出现慢性低度炎症。Roseburia 属是属于产生 SCFA（尤其是丁酸）的共生细菌，可影响免疫维持、结肠运动和抗炎特性。Coprococcus 粪球菌属是大肠中数量较少的细菌，从果糖中产生丁酸，从乳酸中产生丙酸（通过丙烯酸途径）。Faecalibacterium 和 Coprococcus，与多个生活质量分数（quality-of-life ，QoL）呈正相关。Butyricicoccus 也是产生丁酸的菌属，发现丁酸菌与临床指标 MMSE、WAIS 和 Barthel 以及抗炎细胞因子 IFN-γ 呈正相关，与促炎细胞因子 TNF-α 呈负相关，Zhang 等人 的研究表明，与年龄匹配的对照组相比，AD 小鼠模型中的丁酸菌数量明显减少，与我们的数据相符。


所有这些与 AD 相关的粪便属水平降低，与富含 AD 的属相互作用，或许导致 SCFA 的变化，可能参与了 AD 的发病和发展。


![difference_analy_2022-01-06](https://cdn.jsdelivr.net/gh/Achuan-2/PicBed@pic/assets/README/difference_analy_2022-01-06.svg)
> 图 2 A） AD 患者与健康对照者肠道细菌物种的 LEfse 分析。 左图是 LEfse 绘制的**进化分支图**，由内至外辐射的圆圈代表了由门至属的分类级别，在不同分类级别上的每一个小圆圈代表该水平下的一个分类，小圆圈直径大小与相对丰度大小呈正比；着色原则：无显著差异的物种统一着色为黄色，差异物种 Biomarker 跟随组进行着色，红色节点表示在 AD 组别中起到重要作用的微生物类群，绿色节点表示在 C 组别中起到重要作用的微生物类群。图中英文字母表示的物种名称在右侧图例中进行展示。右图是 LDA 值分布柱状图，颜色代表对应分组，柱状图的长度代表差异物种的贡献度大小（即为 LDA Score），图中展示了 LDA Score 大于 3 的条件下不同组间丰度有显著差异的物种，即组内丰度显著高于其它组的 Biomarker。B）AD 患者与健康对照者肠道细菌的物种差异分析：使用 分别对门水平、科水平、属水平进行 Wilcoxon 检验，只展示 Bonferroni 矫正后的 p 值 <0.05 的物种。


### 基于肠道微生物标志物能有效鉴别阿尔兹海默症患者

对 171 个样本的 ASV 丰度数据随机划分为训练集（n=85）和测试集（n=86），先使用全部属水平的丰度数据构建随机森林分类器，通过 MeanDecreaseAccuracy 指标排序，选取了分值最高的 15 个 属水平的丰度数据，重新构建随机森林分类器（命名为 MC），OOB 从 21.18% 降到 17.65%。热图展示了这些标志物在所有样本的分布情况。该随机森林分类器对于训练数据集的预测达到了 100% 的精确度，对于测序数据集的预测则达到 0.78% 的精确度，绘制 ROC 曲线(接收者操作特征曲线，receiver operating characteristic curve），其 AUC（Area Under Curve） 面积为 0.85。

利用 Wilcoxon 秩和检验挑选的组间有显著性差异的属水平特征，提取总丰度前 15 的特征，构建另一个随机森林分类器（命名为 WC），对于测序数据集的预测则达到 0.79% 的精确度，绘制 ROC 曲线，其 AUC 面积为 0.88。WC 分类器的效果好于 MC。

通过 WC 分类器的良好效果，可知利用肠道微生物标志物来区别 AD 患者和健康人是具有很大潜力的。

![RandomForest_2022-01-06](https://cdn.jsdelivr.net/gh/Achuan-2/PicBed@pic/assets/README/RandomForest_2022-01-06.svg)
> 图3 随机森林分类器：A） WC 分类器针对测试集进行预测的混淆矩阵，横坐标代表预测的结果，纵坐标代表实际的标签，方块颜色代表对应区间的数目；B）MC分类器和WC分类器的ROC曲线图，黄色代表MC分类器，绿色代表WC分类器，图例旁标注有分类器各自的AUC值；C）MC分类器挑选的特征在所有样本的热力图，横坐标代表的是样本，红色代表AD组，绿色代表C对照组，纵坐标是通过 MeanDecreaseAccuracy 指标排序挑选的分值最高的 15 个 属水平

### 阿尔兹海默症患者的肠道微生物功能失调

在二级 KEGG 通路上，比较了 45 个通路，并确定了 7 个在 AD 患者组和对照组之间具有明显差异丰度的 KEGG 类别，发现 AD 患者肠道菌群的 Glycan biosynthesis and metabolism（聚糖生物合成与代谢）、Carbohydrate metabolism（碳水化合物代谢），Chemical structure transformation maps （化学结构转换图）显著增加；而在 Cell motility（细胞运动）、“Folding,sorting and degration”（折叠分类讲解）、Energy metabolism（能量代谢)，lipid metabolism(脂质代谢）显著减少（图4)

在三级 KEGG 通路上，比较了 163 个通路，并确定了 14 个在 AD 患者组和对照组之间具有明显差异丰度的 KEGG 类别，发现 AD 患者肠道菌群中的 Lipoic acid metabolism（硫辛酸代谢增加）、Other glycan degradation（其他聚糖降解）、Biosynthesis of terpenoids and steriods（萜类化合物和类固醇的生物合成）、Glycosaminoglycan degradation（糖胺聚糖降解）、Glycosphingolipid biosynthesis - globo and isoglobo series（糖鞘脂的生物合成 - globo 和 isoglobo 系列） 和 Glycosphingolipid biosynthesis - ganglio series （糖鞘脂的生物合成 - ganglio 系列）共 6 条途经显著增加；而 Fatty acid biosynthesis（脂肪酸合成）、Peptidoglycan biosynthesis（肽聚糖生物合成）、Biotin metabolism （生物素代谢）、Porphyrin and chlorophyll metabolism（卟啉和叶绿素代谢）、Bacterial chemotaxis（细菌趋化性）、Flagellar assembly（鞭毛组件）、Sulfur relay system（硫磺中继系统）和 Thiamine metabolism（硫胺素代谢）共 8 条途径在 AD 患者肠道微生物群中显著减少（图5)。

从而推知，肠道微生物群的功能失调可能参与 AD 的发病和发展。
![KEGG2_2022-01-06](https://cdn.jsdelivr.net/gh/Achuan-2/PicBed@pic/assets/README/KEGG2_2022-01-06.svg)
> 图4：AD患者和健康对照组在KEGG Level2 水平上的代谢通路差异

![KEGG3_2022-01-06](https://cdn.jsdelivr.net/gh/Achuan-2/PicBed@pic/assets/README/KEGG3_2022-01-06.svg)

> 图5：AD患者和健康对照组在KEGG Level2 水平上的代谢通路差异


## Discussion

近年来，多组学技术表明，肠道微生物群在促进人类健康方面起着至关重要的作用，因此它们常常被称为“被遗忘的器官”。越来越多的证据表明，肠道微生物群通过各种复杂机制构成维持肠道内环境稳定的关键因素。肠道微生物群不仅被认为是每种胃肠道疾病的促因，而且其影响的分析也被扩展到其他器官，如中枢神经系统（CNS）。大量研究已证实中枢神经系统（central nervous system，CNS）与胃肠道的双向活动，在肠道运动、吸收、内分泌、免疫功能，维持消化道内环境稳定等起重要作用。由于肠道微生物在脑-肠轴之间扮演着十分重要的角色， 我们可以将中枢神经系统、自主神经系统、肠神经系统、消化道以及种类繁多的肠道菌群视为一个整体，即菌-脑-肠轴  。在这个系统中，信号转导发挥调节作用，这个作用是双向的。探索肠道微生物群在神经退行性疾病中的作用和机制是一个新兴的研究领域。肠-脑轴微生物群-肠-脑轴信号已经揭开了精神病学的一个新纪元，有望为精神疾病的诊断和治疗提供新的靶点，并破译其病因。

本文通过分析，发现 AD 患者的肠道微生物群在门、科和属的水平上的分布与健康对照组有着显著不同，物种丰富度降低，Alpha 多样性指数降低，Beta 多样性指数改变，与前人的研究相符合。Alpha 多样性和 Beta  多样性指数都提供了 AD 患者中肠道菌群组成和多样性改变的有力证据。观察到的 AD 患者肠道菌群整体结构失调也表明，与之相关的肠道微生物群的组成也发生了显著变化。通过组间物种组成差异分析，观察到在 AD 患者中 Faecalibacterium、Butyricicoccus、Coprococcus 等产丁酸的菌属减少，Bifidobacterium、Enterococcaceae 等乳酸产生菌显著增加，这貌似说明 AD 患者的肠道菌群倾向于产生乳酸分泌物。通过对肠道菌群的功能预测分析也表明，改变的肠道菌群与患者功能和代谢活动的改变有关，因为 SCFA 等代谢物的变化，可能参与了 AD 的发病和发展。

本文通过分析，发现 AD 患者的肠道微生物群在门、科和属的水平上的分布与健康对照组有着显著不同，物种丰富度降低，Alpha 多样性指数降低，Beta 多样性指数改变。Alpha 多样性和 Beta  多样性指数都提供了 AD 患者中肠道菌群组成和多样性改变的有力证据。观察到的 AD 患者肠道菌群整体结构失调也表明，与之相关的肠道微生物群的组成也发生了显著变化。通过组间物种组成差异分析，观察到在 AD 患者中 Faecalibacterium、Butyricicoccus、Coprococcus 等产丁酸的菌属减少，Bifidobacterium、Enterococcaceae 等乳酸产生菌显著增加，这貌似说明 AD 患者的肠道菌群倾向于产生乳酸分泌物。通过对肠道菌群的功能预测分析也表明，改变的肠道菌群与患者功能和代谢活动的改变有关，因为 SCFA 等代谢物的变化，可能参与了 AD 的发病和发展。

本文使用 AD 患者和健康对照组显著差异的 15 个菌属构建了随机森林分类器，显示了很好的预测性能。虽然没有对其他数据集进行进一步的测试，但也体现了使用肠道微生物作为生物标志物辅助 AD 患者早期非侵入性诊断的潜力。但要确定具体的微生物标志物，还需要更多的样本数据，针对具体地区进行探索，以保证模型的鲁棒性。


## Code Availability

see [qime2_pipeline.md](pipeline/qime2_pipeline.md)

## Reference

see [reference.md](pipeline/reference.md)
