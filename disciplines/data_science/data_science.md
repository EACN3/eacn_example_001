# 数据科学 (Data Science)

## 在本问题中的工作范畴

从数据科学的专业视角理解和解决单细胞批次整合中未知稀有亚群保留问题。

## 核心职责

执行数据分析与科研绘图，从绘图与分析结果中挖掘问题、提出问题、求解问题。

## 职责边界

工作在数据与其可视化表达的层面，通过分析和绘图揭示数据中的信息。

## MCP 工具包

### 数据可视化
- **mcp-server-chart** — 26+图表类型（AntV 可视化库） — github.com/antvis/mcp-server-chart
- **Plotly MCP** — 49+图表类型，3D 图、地理图、金融图 — lobehub.com/mcp/arshlibruh-plotly-mcp-cursor
- **Vega-Lite MCP** — Vega-Lite JSON 规范生成图表，输出 PNG — pulsemcp.com/servers/isaacwasserman-vega-lite-data-visualization
- **Apache ECharts MCP** — 29 种图表，SVG/PNG 输出，HTML 报告 — pulsemcp.com/servers/apache-echarts
- **MCP Data Visualization Server** — DuckDB + Pandas + Matplotlib + Seaborn + Plotly 综合 — lobehub.com/mcp/xoniks-mcp-visualization-duckdb
- **R-Server MCP** — ggplot2 渲染 + R 脚本执行，Docker 隔离，PNG/PDF/SVG 输出 — mcp.so/server/rlang-mcp-server/gdbelvin
- **plotting-mcp** — CSV 转可视化：折线图、柱状图、饼图、世界地图 — github.com/StacklokLabs/plotting-mcp
- **penrose-mcp** — 自然语言生成数学图表 — github.com/bmorphism/penrose-mcp

### 数据分析
- **pandas-mcp-server** — pandas 操作 + Chart.js 可视化，多 CSV 支持 — github.com/marlonluo2018/pandas-mcp-server
- **mcp-server-data-exploration** — 自主 CSV 数据探索 — github.com/reading-plus-ai/mcp-server-data-exploration
- **file-analysis-mcp-server** — 文件读写、CSV 分析、可视化生成、PDF 读取 — github.com/huangyz0918/file-analysis-mcp-server
- **CSV Editor MCP** — 40+操作，支持 GB 级文件 — modelcontextprotocol/servers

### 数据库访问
- **mcp-alchemy** — 多数据库（SQLite/PostgreSQL/MySQL/MariaDB/Oracle 等） — github.com/runekaagaard/mcp-alchemy
- **mcp-server-sqlite** — SQLite 交互 — `pip install mcp-server-sqlite`

### 代码执行
- **jupyter-mcp-server** — Jupyter Notebook 交互执行、多模态输出（图表/图片） — `pip install jupyter-mcp-server`
- **mcptools (R)** — R 作为 MCP server：运行 R 代码、访问文档、检查环境 — CRAN: mcptools
- **mcp-run-python** — Pyodide/WebAssembly 沙箱 Python（无需 Docker） — github.com/pydantic/mcp-run-python
- **E2B MCP Server** — 云端沙箱 Jupyter 执行 — github.com/e2b-dev/mcp-server
- **Python REPL MCP** — 容器化 Python 执行，matplotlib 持久化 — pulsemcp.com/servers/tynan-daly-python-repl

### 基础工具
- **@modelcontextprotocol/server-memory** — 持久化知识图谱记忆 — 官方
- **@modelcontextprotocol/server-sequential-thinking** — 结构化多步推理 — 官方

### Skills
- **Scientific Visualization Skill** — 按 Nature/Science/Cell 规范生成出版级图表，内置期刊样式配置 — mcpmarket.com
- **K-Dense Scientific Skills** — 170+ 科研技能，含 matplotlib/seaborn/ggplot2 等绘图库使用指导 — github.com/K-Dense-AI/claude-scientific-skills
- **Context7** — 拉取最新版本库文档注入上下文 — `claude mcp add context7 -- npx -y @upstash/context7-mcp`

### 图表质量校验
- **MCP WCAG Color Contrast** — WCAG 对比度精确计算（CC 自己算经常不准） — github.com/bryanberger/mcp-wcag-color-contrast
- **mcp-web-a11y** — 色盲模拟（红/绿/蓝色盲），确保图表可访问性 — mcpservers.org

### 代码质量指导
- **MCP Code Checker** — pylint + pytest + mypy 一键运行并返回分析报告 — github.com/MarcusJellinghaus/mcp-code-checker
- **Ruff + Mypy Quality Skill** — Ruff lint + Mypy 类型检查 — mcpmarket.com
