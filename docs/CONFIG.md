# 配置系统说明 / Configuration Guide

## 分层与优先级 / Layers and precedence

配置分三层,优先级从高到低:

```
命令行参数  >  阶段配置(configs/0X_*.yaml)  >  公共配置(configs/00_project.yaml)  >  代码默认值
CLI args    >  stage config                  >  project config                     >  code defaults
```

- **阶段配置**:每个步骤一个文件,只放该步骤的参数;
- **公共配置** `00_project.yaml`:仅存放被多个步骤引用的路径(EGFAnalysisPy 目录、DSurfTomo 可执行文件等)。目前只有 05 与 06b 的入口读取它;01–04 的脚本完全由各自阶段配置驱动;
- 所有 `*.example.yaml` 是模板:复制去掉 `.example` 后缀(或任意命名)再修改,模板本身不要改。

```bash
cp configs/02_correlation.example.yaml configs/02_correlation.yaml
python scripts/02_correlation.py --config configs/02_correlation.yaml
```

## 配置文件一览 / File inventory

| 配置文件 | 服务的入口 | 说明 |
|----------|-----------|------|
| `00_project.example.yaml` | 05、06b | 项目级公共路径 |
| `01_clean.example.yaml` | `scripts/01_clean.py` | 键名与 CLI 选项一一对应 |
| `02_correlation.example.yaml` | `scripts/02_correlation.py` | NoisePy 互相关与叠加参数 |
| `03_stack_qc.example.yaml` | `scripts/03_stack_qc.py` | moveout 可视化参数 |
| `03b_stack_panels.example.yaml` | `scripts/03b_stack_panels.py` | 多面板图参数 |
| `04_export_egf.example.yaml` | `scripts/04_export_egf.py` | 键名与 CLI 选项一一对应 |
| `05_dispersion.example.yaml` | `scripts/05_dispersion.py` | 频散分析全部参数(替代手改 config.py) |
| `06a_line_inversion/line.example.yaml` | 06a 两个入口 | 一条测线一份配置(invert 与 stitch 共用) |
| `06b_area_inversion/dsurftomo.example.yaml` | 06b 三个入口 | 一个反演工程一份配置 |

## 键名约定 / Key naming

01–04 的 YAML 键名与命令行长选项**同拼写**(`input-root`、`stack-dir`……,写成下划线也可),显式命令行参数永远覆盖配置文件;未知键会警告并忽略。05、06a、06b 使用分节结构(`paths:` / `run:` / `prepare:` 等),见各 example 文件内的逐项注释。

## 运行时配置机制 / Runtime-config mechanism

底层程序只认自有格式时,入口脚本负责转换,用户只维护根目录 YAML:

- **05 步**:`05_dispersion.py` 把 YAML 渲染为 `<result_path>/runtime_config.py`(继承 `EGFAnalysisPy/config.py` 并覆盖参数、重算派生量,文件头记录来源配置与生成时间),再交给 `EGFAnalysisPy/run_parallel.py --config` 执行。仓库内模板永不被修改;结果目录里的运行时配置就是这次运行的完整参数存档;
- **06b-1**:`01_prepare.py` 把 `prepare:` 段翻译成 `DSurfTomo/scripts/PrepareDSurfTomoProject.py` 的命令行参数(生成的 `DSurfTomo.in` 等文件都写入工程目录);
- **06b-2/3**:可执行文件路径取 `paths.dsurftomo_bin`(阶段配置),为空则回退到 `00_project.yaml` 的 `dsurftomo_bin` —— 这就是优先级规则的一个实例。

## 覆盖保护 / Overwrite protection

各入口默认**不覆盖**已有结果:01 需 `--overwrite` 才重写输出文件;05 发现 `qc_result/` 非空、06b-1 发现工程目录非空时会报错停止,需显式 `--overwrite`。

## 路径书写建议 / Path conventions

- 示例配置里的 `/path/to/...` 是占位符,必须改成你的真实路径;
- 相对路径相对于**运行命令的目录**(建议始终从项目根目录运行);
- 不要把真实数据路径提交回仓库 —— 复制出的个人配置(不带 `.example`)建议加进 `.gitignore` 或放在仓库外。
