# Gate 5D参数锁定

日期：2026-07-25

## 分析范围

- 仅深化chr16:86,736,383–87,736,383（GRCh37）的ALPS–WMH共享区域。
- 不增加新的疾病、MR或共定位性状。
- 候选变异定义：Gate 5C SuSiE 95%可信集并集，覆盖WMH、mALPS、left ALPS、right ALPS。

## GTEx eQTL

- 数据库：GTEx v8。
- 数据来源：GTEx Portal API v2。
- 组织：13个GTEx脑组织。
- 查询对象：全部SuSiE可信集并集变异。
- 仅使用GTEx已定义为显著的single-tissue cis-eQTL记录。
- NES表示GTEx替代等位基因对标准化表达的效应；不直接等同于ALPS升高/降低方向，除非完成等位基因方向协调。

## MAGMA

- 软件：MAGMA v1.10。
- LD参考：1000 Genomes EUR，N=489。
- GWAS：mALPS、left ALPS、right ALPS。
- 基因坐标：NCBI37.3。
- SNP到基因窗口：基因体上下游各10 kb。
- 基因模型：SNP-wise mean。
- 基因显著性：Bonferroni 0.05/各表型实际检验基因数；同时报告BH-FDR。
- 通路：MAGMA competitive gene-set analysis。
- 基因集：GO Biological Process 2023、Reactome 2022、MSigDB Hallmark 2020。
- 通路显著性：各表型全部预设基因集内BH-FDR<0.05；Bonferroni作为严格标准。

## Gene prioritization

证据维度：

1. chr16 SuSiE可信集或高PIP变异证据；
2. GTEx脑组织显著eQTL证据；
3. MAGMA基因关联证据；
4. 位置映射（基因体±10 kb）；
5. 三个ALPS表型证据一致性。

优先级仅表示统计和功能证据汇总，不证明该基因为因果基因。
