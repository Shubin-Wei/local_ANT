# 06b 面状台阵速度结构反演(DSurfTomo)/ Areal-Array Inversion

## 适用场景

台站覆盖一个二维区域(面状台阵),射线路径在区域内交叉。使用 **DSurfTomo**(Fang et al., 2015)将面波频散数据直接反演为三维横波速度模型,无需先构建相/群速度图。

> DSurfTomo 只属于本分支;线状台阵请使用 [06a 线状反演](line_inversion.md)。

## 前置条件

1. **编译 DSurfTomo**(见 [INSTALL.md](INSTALL.md) §4):

   ```bash
   cd DSurfTomo && python configure && cd ..
   # configure 会进入 src/ 执行 make 并把可执行文件复制到 bin/DSurfTomo
   ```

2. **GMT 6+**(绘图,INSTALL §5);
3. 第 05 步的 `qc_result/` 频散结果。

在 `configs/00_project.yaml` 里设置 `dsurftomo_bin: DSurfTomo/bin/DSurfTomo`(或在阶段配置 `paths.dsurftomo_bin` 中指定,后者优先)。

## 编号入口(日常使用)

一个反演工程一份配置:

```bash
cp configs/06b_area_inversion/dsurftomo.example.yaml configs/06b_area_inversion/my_area.yaml
# 编辑 my_area.yaml:input_dir、project_dir、period_step、网格与深度等
conda activate noisepy

# 1. 准备工程:qc_result → surfdata.dat + MOD + DSurfTomo.in(写入 project_dir)
python scripts/06b_area_inversion/01_prepare.py --config configs/06b_area_inversion/my_area.yaml

# 2. 运行反演(在 project_dir 内执行,控制台输出同时存 dsurftomo_run.log)
python scripts/06b_area_inversion/02_invert.py --config configs/06b_area_inversion/my_area.yaml

# 3. GMT 绘图(默认深度切片;图件写入 project_dir/figs)
python scripts/06b_area_inversion/03_plot.py --config configs/06b_area_inversion/my_area.yaml
```

`prepare:` 段的键与 `DSurfTomo/scripts/PrepareDSurfTomoProject.py` 的命令行参数一一对应(入口自动转换),要点:

| 键 | 含义 |
|----|------|
| `curve_type` | `both` / `phase` / `group` |
| `period_step` | 周期网格步长(s);不必与 05 步 DeltaT 完全一致,只要生成的周期点能在 05 步输出曲线的周期采样中找到(取 DeltaT 或其整数倍最稳妥) |
| `grid_spacing` 或 `grid_count` | 水平网格,二选一:`[dvxd,dvzd]`(度)或 `[nx,ny]` |
| `depths` | 深度节点(km,逗号分隔字符串) |
| `minvel` / `velgrad` / `t1` / `t2` / `model_vel` | 初始模型估计(留 null 用自动值) |
| `weight` / `damp` / `sublayers` / `nsrcmax` | 反演平滑、阻尼等 |

## 内部脚本(调试用)

`DSurfTomo/scripts/` 下的原有脚本是 06b 编号入口调用的底层实现,全部保留;日常使用编号入口即可,需要单步调试时也可手动运行:

```bash
cd DSurfTomo/scripts

# 一体式(交互式问答,推荐手动方式用这个)
python PrepareDSurfTomoProject.py

# 或分步:
python extractnew.py                    # qc_result → surfdata.dat(+ surfdata.json)
python GenerateIniMOD.py                # → MOD 初始模型
python GenerateDSurfTomoInputFile.py    # → DSurfTomo.in
python GenerateTrueMOD.py               # 可选 → MOD.true(棋盘测试)

# 运行反演(在工程目录内)
cd /path/to/project && /path/to/local_ANT/DSurfTomo/bin/DSurfTomo
```

`extractnew.py` 读取 `qc_result/` 下的 `CDisp.*.txt`、`GDisp.*.txt`,台站坐标取自文件头,无需额外 station 文件;周期步长取 05 步 `DeltaT` 或其整数倍。

## 输出

反演在工程目录产生(已实测确认):

- `DSurfTomo.inMeasure.dat` — 最终模型(经度 纬度 深度 速度 四列,绘图脚本默认读取);
- `DSurfTomo.inMeasure.dat.iterNNN` — 每次迭代的模型;
- `residualFirst.dat` / `residualLast.dat` — 首/末次迭代残差;
- `DSurfTomo.in.log` — 程序日志;`dsurftomo_run.log` — 编号入口保存的完整控制台输出。

## 绘图

全部绘图统一由 `03_plot.py` 驱动:把要跑的脚本写进配置的 `plot.scripts` 列表(按序执行,`--project-dir`/`--output-dir` 自动传入;需要额外参数的脚本用 `name:`+`args:` 形式):

```yaml
plot:
  scripts:
    - plotdepth_subplot.gmt              # 深度切片(自动传 --depths)
    - name: plotcross_subplot.gmt        # 剖面
      args: ["--profiles", "lon1/lat1/lon2/lat2/label"]
    - plotray_subplot.gmt                # 射线路径
    - plotchecker_subplot.gmt            # 棋盘测试
```

四个 GMT 脚本均接受统一的 `--project-dir`/`--output-dir` 接口;各自的更多可选参数见其 `--help`。频散数据统计图可另行运行 `python DSurfTomo/scripts/plot_dispersion_period_counts.py --surfdata <工程目录>/surfdata.dat`。

## 常见注意事项

- 工程目录(`project_dir`)放在仓库外,不要把反演结果提交进 Git;
- `01_prepare.py` 发现工程目录非空会拒绝执行,需 `--overwrite`;
- 棋盘测试:用内部脚本 `GenerateTrueMOD.py` 生成 `MOD.true` 后按 DSurfTomo 文档设置正演开关。
