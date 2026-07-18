# TODO

## 功能 / Features

- [ ] 在线状台阵速度结构反演流程中增加敏感核计算模块。
      Add a sensitivity-kernel computation module to the linear-array inversion workflow.

## 测试 / Tests

- [ ] 用 TestData 小数据集完整跑通 01–05 步,记录各步耗时与输出规模。
      Run steps 01–05 end-to-end on the small TestData set; record runtime and output sizes.
- [ ] 编译 DSurfTomo 并用示例工程验证 06b 三个入口(01_prepare / 02_invert / 03_plot)端到端可用。
      Compile DSurfTomo and verify the three 06b entries end-to-end on a sample project.
- [ ] 安装 evodcinv 后用少量台站对(`--max-files`)冒烟验证 06a 两个入口。
      After installing evodcinv, smoke-test both 06a entries with a few pairs (`--max-files`).

## 文档 / Documentation

- [ ] 补充 docs/FAQ.md 与 docs/FAQ_en.md:随实际使用持续收集安装、配置、运行与结果分析问题。
      Keep extending docs/FAQ.md and docs/FAQ_en.md as real-world issues come up.

## 待确认 / To be confirmed

- [ ] evodcinv 的版本号:原项目环境中已无法考证(正演内核 disba 0.6.1 确认可用);确认后在 docs/INSTALL.md §2e 与 environments/noisepy.yml 中锁定。
      The evodcinv version could not be recovered from the original environments (disba 0.6.1 confirmed working); pin it in docs/INSTALL.md §2e and environments/noisepy.yml once known.
- [ ] `configs/05_dispersion.example.yaml` 中 `picker.ref_T` / `ref_T2` 高级选项的推荐用法,需原作者补充说明。
      Recommended usage of the advanced `picker.ref_T` / `ref_T2` options needs input from the original author.
