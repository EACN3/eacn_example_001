# 机器学习 (Machine Learning)

## 在本问题中的工作范畴

从机器学习的专业视角理解和解决单细胞批次整合中未知稀有亚群保留问题。

## 核心职责

从现有机器学习研究中发掘新的看待问题的视角，给出合适的解决方案，并提供参考论文或链接。

## 职责边界

工作在算法与模型的思想层面，通过文献研究形成对问题的新理解。

## MCP 工具包

### 文献检索
- **paper-search-mcp** — 20+学术源一站检索（arXiv、PubMed、Semantic Scholar、Google Scholar、Crossref、OpenAlex 等） — `pip install paper-search-mcp`
- **PaperMCP** — arXiv、HuggingFace、Google Scholar、OpenReview、DBLP、PapersWithCode 统一检索 — github.com/ScienceAIHub/PaperMCP
- **research_hub_mcp** — Rust 实现，CrossRef（1.3亿+论文）、Semantic Scholar、arXiv、PMC、OpenReview、SSRN 等，智能路由 — github.com/Ladvien/research_hub_mcp
- **Academix** — OpenAlex（10万+免费调用/天）、DBLP、Semantic Scholar、arXiv、CrossRef 统一检索 — github.com/xingyulu23/Academix
- **arxiv-mcp-server** — arXiv 论文搜索、下载、阅读 — github.com/blazickjp/arxiv-mcp-server
- **semanticscholar-MCP-Server** — Semantic Scholar 论文、作者、引用网络 — github.com/JackKuo666/semanticscholar-MCP-Server
- **Google-Scholar-MCP-Server** — Google Scholar 检索 — github.com/JackKuo666/Google-Scholar-MCP-Server

### ML 专业资源
- **openreview-mcp-server** — ICML/ICLR/NeurIPS 会议论文检索 — github.com/anyakors/openreview-mcp-server
- **PapersWithCode MCP** — 论文+代码+数据集+benchmark — pulsemcp.com/servers/paperswithcode
- **Hugging Face MCP Server（官方）** — 模型、数据集、Spaces、论文 — huggingface.co/docs/hub/en/hf-mcp-server
- **wandb-mcp-server（官方）** — W&B 实验追踪与分析 — github.com/wandb/wandb-mcp-server
- **mlflow-mcp** — MLflow 实验追踪与模型注册 — github.com/kkruglik/mlflow-mcp

### 向量数据库（语义检索）
- **Qdrant MCP（官方）** — 语义记忆层 — `uvx mcp-server-qdrant`
- **Chroma MCP（官方）** — ChromaDB 向量检索 — github.com/chroma-core/chroma-mcp
- **Pinecone MCP** — Pinecone 向量数据库 — pulsemcp.com/servers/sirmews-pinecone

### 代码原型验证
- **jupyter-mcp-server** — Notebook 交互执行 — `pip install jupyter-mcp-server`
- **mcp-run-python** — Pyodide/WebAssembly 沙箱 Python 执行（无需 Docker） — github.com/pydantic/mcp-run-python
- **E2B MCP Server** — 云端沙箱 Jupyter 执行 — github.com/e2b-dev/mcp-server

### 基础工具
- **@modelcontextprotocol/server-memory** — 持久化知识图谱记忆 — 官方
- **@modelcontextprotocol/server-sequential-thinking** — 结构化多步推理 — 官方

### Skills（直接装在 CC 上的技能）
- **K-Dense Scientific Skills** — 170+ 科研技能，含 ML 相关库使用指导 — github.com/K-Dense-AI/claude-scientific-skills
- **Context7** — 拉取最新版本库文档注入上下文，帮助写出正确的代码 — `claude mcp add context7 -- npx -y @upstash/context7-mcp`

### 代码质量指导
- **MCP Code Checker** — pylint + pytest + mypy 一键运行并返回分析报告 — github.com/MarcusJellinghaus/mcp-code-checker
- **Ruff + Mypy Quality Skill** — Ruff lint + Mypy 类型检查 — mcpmarket.com
