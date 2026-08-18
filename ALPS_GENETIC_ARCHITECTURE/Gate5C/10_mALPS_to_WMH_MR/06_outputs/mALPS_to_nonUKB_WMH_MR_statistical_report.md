# mALPS → 非UKB WMH探索性MR统计报告

## 1. Gate与执行范围

本阶段仅评估并执行mALPS→WMH，不运行left ALPS或right ALPS MR。

- Exposure：mALPS，European ancestry，样本量中位数约30,952；
- Outcome：Traylor et al. 2016 WMH，European ancestry，N=3,670；
- 样本重叠：未发现UK Biobank或其他已知重叠；
- 初始mALPS独立IV：12；
- Traylor直接覆盖：10；
- 保守删除缺少结局EAF的回文SNP后：8；
- 1000G EUR代理补救：加入2个r²=1.0非回文代理；
- 最终IV：10，恰好达到预设MR启动标准。

最终工具强度：

- minimum F=29.84；
- median F=37.32；
- 所有IV的F>10。

Gate判定：**GO — exploratory MR completed**。

## 2. 主要MR结果

| 方法 | IV | Beta | SE | 95% CI | P |
|---|---:|---:|---:|---:|---:|
| IVW | 10 | -0.4697 | 0.1550 | -0.7734至-0.1660 | 0.00243 |
| Weighted median | 10 | -0.4845 | 0.2013 | -0.8790至-0.0900 | 0.0161 |
| MR-Egger slope | 10 | -0.8170 | 0.8248 | -2.4335至0.7996 | 0.351 |
| MR-RAPS | 10 | -0.4610 | 0.1635 | -0.7814至-0.1405 | 0.00481 |

主分析IVW、weighted median和MR-RAPS方向一致且均达到P<0.05。MR-Egger斜率同方向，但由于只有10个IV，估计不精确。

效应方向表示：

> 遗传预测的较高mALPS与较低WMH burden相关。

Beta使用原始mALPS暴露尺度，WMH为Traylor研究中的标准化表型；不能直接转换为临床ALPS阈值或WMH体积变化。

## 3. 多效性与异质性

### Cochran Q

- IVW：Q=6.180，df=9，P=0.722；
- MR-Egger：Q=5.997，df=8，P=0.648。

未发现明显SNP间异质性。

### MR-Egger截距

- intercept=0.0195；
- SE=0.0454；
- P=0.680。

未发现明显方向性水平多效性。

### MR-PRESSO

- Global test：RSSobs=7.377；
- P=0.762；
- outlier=0。

未发现MR-PRESSO异常SNP，因此无outlier-corrected estimate。

## 4. 稳定性分析

### Leave-one-out

逐一删除任何一个SNP后：

- IVW beta均保持负向；
- beta范围约-0.536至-0.389；
- 所有leave-one-out P<0.016。

结果不由单个SNP驱动。

### 代理SNP敏感性

仅使用8个直接匹配、非代理IV：

- IVW beta=-0.4751；
- 95% CI -0.8231至-0.1271；
- P=0.00746。

说明主结果不依赖两个代理SNP。

### 去除chr16工具变量

删除rs4843555后：

- IVW beta=-0.4659；
- 95% CI -0.7859至-0.1459；
- P=0.00432。

结果不由既有chr16 ALPS–WMH共享区域单独驱动。

### 去除结局关联最强SNP

删除rs10817103后：

- IVW beta=-0.3893；
- 95% CI -0.7058至-0.0728；
- P=0.0159。

方向和统计学显著性保持。

### Steiger近似检验

- 总R² exposure=0.01188；
- 总R² outcome=0.00428；
- 9/10个SNP支持mALPS→WMH方向；
- 汇总方向判定为正确。

该检验将两项性状视为标准化定量表型，并以暴露EAF近似结局EAF，因此仅作辅助证据。

## 5. 统计学结论

本次预设的独立非UKB WMH探索分析获得一致的负向MR估计：

> Genetic liability to higher mean DTI-ALPS index was associated with lower WMH burden in the independent Traylor et al. 2016 cohort.

证据特点：

- IVW、weighted median和MR-RAPS一致；
- 无明显异质性；
- Egger截距和MR-PRESSO未提示明显水平多效性；
- leave-one-out、去代理、去chr16及去最强结局SNP后方向稳定；
- 不受已知UK Biobank样本重叠影响。

## 6. 解释边界

该结果应定位为**支持性、探索性方向证据**，不能作为确定性因果结论，原因包括：

1. 最终IV数量恰好为10；
2. WMH结局样本仅N=3,670，统计精度有限；
3. WMH数据来自缺血性卒中患者队列，外推至一般人群需谨慎；
4. 两个原始IV由r²=1.0代理替代；
5. 尚未在更大CHARGE非UKB完整数据中重复；
6. 本结果不能替代chr16正式非UKB共定位。

因此，当前论文中宜使用“supports a potential directional relationship”或“exploratory MR evidence”，不应使用“proves that reduced ALPS causes WMH”。

## 7. Gate更新

| 补充目标 | 当前状态 |
|---|---|
| 非UKB WMH正式chr16共定位 | 等待CHARGE完整signed区域统计；尚未完成 |
| mALPS→WMH MR可行性 | PASS |
| mALPS→非UKB WMH探索MR | 已完成；结果阳性且敏感性一致 |
| left ALPS→WMH MR | 不执行，IV<10 |
| right ALPS→WMH MR | 不执行，结局协调后IV<10 |

