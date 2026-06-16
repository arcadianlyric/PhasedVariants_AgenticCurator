## PhasedVariants AgenticCurator

基于多智能体 LLM 系统和 RAG 技术，从分相全基因组测序数据中自动完成基因/变异解读。

### 背景

基因组分相（phasing）确定哪些遗传变异共同位于同一条染色体（单倍型）上，这对于解析复合杂合性和顺式调控效应至关重要。Complete Genomics 的 [cWGS pipeline](https://github.com/Complete-Genomics/DNBSEQ_Complete_WGS/tree/dev) 使用 DeepVariant 和 HapCUT2 生成高保真度的分相 VCF 文件。

将分相变异与基因功能、致病机制和临床意义关联起来，传统上依赖领域专家的人工解读。这一过程耗时长、难以规模化，且不同解读者之间缺乏一致性。

本项目使用多智能体 LLM 工作流（结合检索增强生成 RAG、知识图谱和迭代质量控制）自动化完成这一解读过程。

**关键词：** 多智能体系统, RAG, FAISS, LLM, 单倍型分相, 基因/变异解读, 知识图谱, PubMed, GeneCards, arXiv, Tavily, 幻觉抑制, LangChain, LangGraph, LangSmith

---

### 系统架构

```mermaid
flowchart TD
    subgraph Input
        VCF[Phased VCF<br/>DeepVariant + HapCUT2]
        KG[(PrimeKG<br/>Knowledge Graph)]
    end

    subgraph Step1[Step 1: Variant Discovery]
        VEP[VEP Annotation]
        HAP[Haplotype Analysis<br/>HIGH impact on both copies]
        NET[Gene Network<br/>Disease / Pathway / Phenotype]
    end

    subgraph Step2[Step 2: Literature Retrieval]
        PM[PubMed]
        GC[GeneCards]
        AX[arXiv]
        TV[Tavily Web Search]
    end

    subgraph Step3[Step 3: LangGraph Multi-Agent Curation]
        SG[StateGraph Orchestrator]
        OA[Generate Node<br/>Output Agent / DeepSeek]
        RA[Review Node<br/>Review Agent / Grok]
        RF[Refine Node<br/>Revise with feedback]
        LS[(LangSmith Trace<br/>optional)]
        SG --> OA
        OA -->|analysis| RA
        RA -->|score >= 7.0| RPT
        RA -->|needs revision| RF
        RF -->|updated analysis| RA
        SG -. trace .-> LS
    end

    subgraph Output
        RPT[Curated Report<br/>per gene, scored 0-10]
    end

    VCF --> VEP --> HAP --> NET
    KG --> NET
    NET -->|gene_list.json| Step2
    PM & GC & AX --> FAISS[FAISS Vector Store]
    TV --> DirectCtx[Direct Context]
    FAISS & DirectCtx --> OA
    KG --> OA
```

系统分三步运行：

1. **变异发现** -- 解析分相 VCF，识别两个单倍型上均具有 VEP HIGH 影响的变异，通过 PrimeKG 构建基因网络。
2. **文献检索** -- 多源渐进式搜索（基因+疾病+变异 -> 基因+疾病 -> 仅基因），覆盖 PubMed、GeneCards、arXiv 和 Tavily。
3. **多智能体解读** -- 输出智能体生成基于文献的分析；审查智能体独立评估质量（5 个维度，0-10 分制）并标记幻觉。循环迭代直至达到质量阈值。

---

### 结果

#### 第一步：变异发现与基因网络

```bash
cd src
vcf_file=../data/HG002_exon.vep.vcf.gz
kg_path=../db/kg.csv
ref_fai=../db/GCA_000001405.15_GRCh38_no_alt_analysis_set.ercc.fa.fai
python explore_phased_vcf.py --vcf_file $vcf_file --kg_path $kg_path --ref_fai $ref_fai
```

- **输入：** 带 VEP 注释的分相 VCF（格式：`Uploaded_variation|Location|Allele|Gene|Feature|...`）
- **输出：** `results/network_graph.html`（交互式基因-疾病-通路网络），`results/gene_associations.json`
- 用户从网络中选择感兴趣的基因，写入 `gene_list.json` 供后续解读使用。

![network_graph](images/network_graph.jpg)

#### 第二步：文献检索

```bash
cd src
python literature_retrieval.py
```

- **输入：** `gene_list.json`
- **输出：** `results/{gene}_comprehensive_literature.json`, `results/{gene}_pubmed_response.txt`
- 数据源：PubMed（生物医学摘要）、GeneCards（基因功能/别名/通路）、arXiv（计算生物学预印本）、Tavily（实时网络搜索 + AI 综合答案）。
- 渐进式搜索策略，每条结果包含 `query_used` 字段以保证可追溯性。

#### 第三步：多智能体解读（推荐）

```bash
cd src
python llm_agentic.py            # v2 多智能体（默认，推荐）
python llm_agentic.py --legacy   # v1 规划+执行+反思
```

- **输入：** `gene_list.json`、第二步的文献文件、知识图谱
- **输出：** `results/{gene}_multiagent_analysis.json`, `results/{gene}_multiagent_report.md`

用于对比的单智能体方案：

```bash
python llm_queryAlone.py     # 仅 LLM（无上下文，基线）
python llm_augmented.py      # LLM + PubMed 文本（简单增强）
python llm_rag.py            # LLM + FAISS RAG（单智能体）
```

#### 性能指标

| 指标 | 数值 |
|------|------|
| 多智能体解读耗时 | 每个基因约 100 秒 |
| 多源文献检索耗时 | 每个基因约 12 秒 |
| 质量评分范围 | 6-9/10 |
| 幻觉降低率（加入 Tavily 后） | 约 40% |

---

### 材料与方法

#### 数据源

| 数据源 | 用途 | 说明 |
|--------|------|------|
| [PrimeKG](https://zitniklab.hms.harvard.edu/projects/PrimeKG/) | 知识图谱：基因-疾病-通路-表型关联 | 下载 `kg.csv` 到 `./db` |
| [PubMed](https://pubmed.ncbi.nlm.nih.gov/) | 同行评审的生物医学摘要 | E-utilities API |
| [GeneCards](https://www.genecards.org/) | 基因功能、别名、通路、疾病关联 | HTML 解析 |
| [arXiv](https://arxiv.org/) | 计算生物学预印本 | Atom API |
| [Tavily](https://tavily.com/) | 实时网络搜索 + AI 综合答案 | REST API（需 API key） |
| [ClinVar](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/) | 变异级别的临床显著性 | 表格下载 |

#### 工具与算法

| 组件 | 工具 | 选择理由 |
|------|------|----------|
| 分相 | [HapCUT2](https://github.com/vibansal/HapCUT2) | 短读长单倍型组装的主流工具 |
| 变异检测 | [DeepVariant](https://github.com/google/deepvariant) | 深度学习 caller，PrecisionFDA 基准测试中准确度最高 |
| 变异注释 | [VEP](https://www.ensembl.org/vep) | Ensembl 标准工具，支持后果预测和 IMPACT 分类 |
| 向量库 | [FAISS](https://github.com/facebookresearch/faiss) | 快速近似最近邻搜索，可处理百万级向量 |
| 文本嵌入 | [sentence-transformers/all-MiniLM-L6-v2](https://www.sbert.net/) | 速度与质量的平衡，适合生物医学文本 |
| LLM（输出智能体） | [DeepSeek](https://deepseek.com/) | 最初（2025 年）因性价比选择；备选：[MiniMax](https://platform.minimaxi.com/)、[Kimi](https://platform.moonshot.cn/) |
| LLM（审查智能体） | [Grok](https://x.ai/) | xAI 模型，搜索锚定能力强；不同 LLM 减少审查中的共同盲区 |
| LLM 框架 | [LangChain](https://github.com/langchain-ai/langchain) | 文本分块、FAISS 集成、检索器抽象 |
| 智能体编排 | [LangGraph](https://github.com/langchain-ai/langgraph) + 可选 [LangSmith](https://smith.langchain.com/) | 有状态 Output -> Review -> Refine 图编排；LangSmith 是可选项，仅用于 trace 观测 |

#### 多智能体架构（v2）

v2 工作流将生成和评估分离为两个独立智能体：

- **输出智能体（Output Agent / DeepSeek）** -- 接收 RAG 上下文（FAISS: PubMed + GeneCards + arXiv）、Tavily 网络上下文和知识图谱关联。生成结构化的基因分析。后续迭代中，接收审查智能体的反馈并据此修改。
- **审查智能体（Review Agent / Grok）** -- 使用不同的 LLM，独立地在 5 个维度（完整性、准确性、证据支持、清晰度、临床实用性）上评分。通过交叉引用文献上下文标记潜在幻觉。返回结构化 JSON 反馈。
- **LangGraph 编排器（Orchestrator）** -- 将输出、审查和修订定义为 `StateGraph` 节点。条件边在质量达到阈值（7.0/10）时提前停止，否则继续进入修订，直到达到最大迭代次数。

生成（DeepSeek）和审查（Grok）使用不同 LLM 是有意为之：消除同一模型评估自身输出时的共同盲区。Grok 的搜索锚定能力也有助于事实核查。这种设计确保审查者不受生成过程的偏见影响（不同模型、独立系统提示词、独立关注点）。

LangSmith 是可选项。不设置 `LANGSMITH_API_KEY` 时 workflow 仍可正常运行；只有需要云端 tracing/debugging 时才需要设置。启用后，LangGraph 的生成、审查和修订节点 trace 会记录到 `LANGCHAIN_PROJECT`（默认 `phasedvariants-agentic-curator`）。

#### 方案对比

| 方案 | 上下文 | 智能体数 | 质量控制 | 幻觉抑制 |
|------|--------|----------|----------|----------|
| `llm_queryAlone` | 无 | 1 | 无 | 无 |
| `llm_augmented` | PubMed 文本 | 1 | 无 | 文献锚定 |
| `llm_rag` | FAISS (PubMed+GeneCards+arXiv) + Tavily | 1 | 无 | RAG + 网络事实 |
| `llm_agentic --legacy` | 同 RAG | 7（规划/执行/反思） | 自我反思循环 | RAG + 反思 |
| `llm_agentic`（v2，默认） | 同 RAG | 2（DeepSeek + Grok） | 独立审查智能体（不同 LLM） | RAG + 跨模型审查 + 幻觉标记 |

---

### 环境配置

1. 创建 conda 环境：
   ```bash
   conda env create -f environment.yml
   conda activate dev
   ```

2. 在 `.env` 中设置 API key（参考 `.env.example`）：
   ```
   DEEPSEEK_API_KEY=your-key       # 输出智能体
   GROK_API_KEY=your-key            # 审查智能体（https://console.x.ai/）
   TAVILY_API_KEY=your-key          # 网络搜索
   PUBMED_EMAIL=your.email@example.com
   LANGSMITH_API_KEY=your-key       # 可选：仅用于 LangGraph 云端 tracing/debugging
   LANGCHAIN_PROJECT=phasedvariants-agentic-curator  # 可选 LangSmith 项目名
   ```

3. 下载 PrimeKG 知识图谱：
   ```bash
   # 从 https://zitniklab.hms.harvard.edu/projects/PrimeKG/ 下载 kg.csv
   # 放置到 ./db/kg.csv
   ```

4. 准备分相 VCF 输入（来自 [cWGS pipeline](https://github.com/Complete-Genomics/DNBSEQ_Complete_WGS/tree/dev) 的输出）。

---

### 项目结构

```
.
├── src/
│   ├── explore_phased_vcf.py       # 第一步：变异发现
│   ├── vep.py                      # VEP 注释解析
│   ├── primeKG.py                  # 知识图谱查询
│   ├── visualize_gene_pathway_disease_phenotype.py
│   ├── literature_retrieval.py     # 第二步：多源文献检索
│   ├── llm_agentic.py             # 第三步：主入口（默认 v2）
│   ├── multi_agent_workflow.py     # v2 多智能体（输出智能体 + 审查智能体）
│   ├── agentic_framework.py        # v1 旧版（规划+执行+反思）
│   ├── planning_agent.py           # v1 规划智能体
│   ├── reflection_agent.py         # v1 反思智能体
│   ├── llm_rag.py                  # 单智能体 RAG 方案
│   ├── llm_augmented.py            # 简单文献增强
│   ├── llm_queryAlone.py           # 仅 LLM 基线
│   └── config.py                   # API key 管理
├── data/                           # 分相 VCF 输入
├── db/                             # PrimeKG kg.csv, 参考基因组索引
├── results/                        # 所有输出（JSON, 报告, FAISS 索引）
├── images/                         # 文档图片
├── gene_list.json                  # 感兴趣的基因（用户创建）
├── environment.yml                 # Conda 环境
└── .env                            # API keys（不提交到版本库）
```

---

### 讨论与后续工作

**当前局限：**
- 变异解读（`_variant_curator_agent`）目前为占位实现，完整 ClinVar 集成正在进行中。
- 调控元件分析（启动子、增强子）尚未实现。
- 质量评分依赖 LLM 的自评能力，可能存在不一致性。
- GeneCards 的网页抓取可能因网站改版而失效，直接 API 集成更为稳健。

**计划改进：**
1. 完整的基于 ClinVar 的变异解读，包含结构化证据等级（致病性/可能致病/VUS）。
2. 调控元件注释（启动子、增强子、UTR 变异）。
3. 系统评估：将多智能体输出与人工解读的金标准基因报告进行对比。
4. 增加 LangGraph 持久化 checkpoint，支持长任务断点续跑。
5. 大基因 panel（100+ 基因）的批处理优化。

---

### 参考文献

1. [cWGS](https://github.com/Complete-Genomics/DNBSEQ_Complete_WGS/tree/dev) - Complete WGS pipeline (DeepVariant + HapCUT2)
2. [DeepVariant](https://github.com/google/deepvariant) - 深度学习变异检测
3. [HapCUT2](https://github.com/vibansal/HapCUT2) - 单倍型分相
4. [VEP](https://www.ensembl.org/vep) - 变异效应预测
5. [PrimeKG](https://zitniklab.hms.harvard.edu/projects/PrimeKG/) - 精准医学知识图谱
6. [PubMed](https://pubmed.ncbi.nlm.nih.gov/) - 生物医学文献
7. [GeneCards](https://www.genecards.org/) - 基因数据库
8. [arXiv](https://arxiv.org/) - 预印本
9. [Tavily](https://tavily.com/) - 网络搜索 API
10. [ClinVar](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/) - 变异临床显著性
11. [FAISS](https://github.com/facebookresearch/faiss) - 向量相似性搜索
12. [LangChain](https://github.com/langchain-ai/langchain) - LLM 应用框架
13. [DeepSeek](https://deepseek.com/) - 输出智能体 LLM（2025 年最初选择；备选：MiniMax/Kimi）
14. [Grok](https://x.ai/) - 审查智能体 LLM（xAI，搜索锚定能力强）
15. [Sentence Transformers](https://www.sbert.net/) - 文本嵌入
16. [stLFR](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6499310/) - 共条形码 NGS 读段
17. [liftOver](http://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/) - 基因组坐标转换
