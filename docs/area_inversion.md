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

## 方式一:编号入口(推荐)

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
| `period_step` | 周期网格步长(s),**必须与 05 步 DeltaT 一致** |
| `grid_spacing` 或 `grid_count` | 水平网格,二选一:`[dvxd,dvzd]`(度)或 `[nx,ny]` |
| `depths` | 深度节点(km,逗号分隔字符串) |
| `minvel` / `velgrad` / `t1` / `t2` / `model_vel` | 初始模型估计(留 null 用自动值) |
| `weight` / `damp` / `sublayers` / `nsrcmax` | 反演平滑、阻尼等 |

## 方式二:交互式手动流程

`DSurfTomo/scripts/` 下的原有脚本全部保留,可手动分步执行:

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

`extractnew.py` 读取 `qc_result/` 下的 `CDisp.*.txt`、`GDisp.*.txt`,台站坐标取自文件头,无需额外 station 文件;周期步长须与 05 步 `DeltaT` 一致。

## 输出

反演在工程目录产生(已实测确认):

- `DSurfTomo.inMeasure.dat` — 最终模型(经度 纬度 深度 速度 四列,绘图脚本默认读取);
- `DSurfTomo.inMeasure.dat.iterNNN` — 每次迭代的模型;
- `residualFirst.dat` / `residualLast.dat` — 首/末次迭代残差;
- `DSurfTomo.in.log` — 程序日志;`dsurftomo_run.log` — 编号入口保存的完整控制台输出。

## 绘图

`03_plot.py` 默认调用 `plotdepth_subplot.gmt`(深度切片,自动从 `DSurfTomo.in` 读网格参数,连续色标)。其余 GMT 脚本可加入配置的 `plot.scripts` 列表或手动运行:

```bash
cd DSurfTomo/scripts
bash plotcross_subplot.gmt   --project-dir /path/to/project --profiles "lon1/lat1/lon2/lat2/label"   # 剖面
bash plotray_subplot.gmt     --project-dir /path/to/project    # 射线路径
bash plotchecker_subplot.gmt --project-dir /path/to/project    # 棋盘测试
python plot_dispersion_period_counts.py --surfdata /path/to/surfdata.dat   # 频散数据统计
```

## 常见注意事项

- 工程目录(`project_dir`)放在仓库外,不要把反演结果提交进 Git;
- `01_prepare.py` 发现工程目录非空会拒绝执行,需 `--overwrite`;
- 棋盘测试:方式二里用 `GenerateTrueMOD.py` 生成 `MOD.true` 后按 DSurfTomo 文档设置正演开关。
