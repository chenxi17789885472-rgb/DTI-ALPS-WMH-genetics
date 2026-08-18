# mALPS → WMH探索性MR参数锁定

- Exposure：mALPS GWAS，GRCh37，European ancestry。
- Outcome：Traylor et al. 2016 WMH GWAS，European ancestry，N=3,670，不含UK Biobank。
- 原始工具：P<5×10^-8，1000G EUR clumping，r²<0.001，10 Mb，共12个。
- 正式启动门槛：协调后可用独立IV不少于10个。
- 缺失或因回文而方向不确定的IV：允许使用1000G EUR中r²≥0.8的非回文代理。
- 代理选择：优先r²最高，其次结局文件存在、非回文、等位基因可协调。
- 主方法：IVW。
- 辅助方法：weighted median；IV数允许时报告MR-Egger、Cochran Q和leave-one-out。
- 若代理后仍少于10个：停止MR，不报告因果估计。
- 本分析为独立小样本WMH结局的探索性可行性分析，不替代正式CHARGE非UKB chr16共定位复制。
