# Development rules

[中文](engineering.zh.md) | English

1. Make the smallest sufficient change; leave unrelated code alone.
2. Keep one purpose per commit; separate features, fixes and refactoring.
3. Give each module a clear responsibility and each function one explicit task.
4. Reduce deep nesting and duplication. At about 500 lines, review responsibilities; line count is a reminder, not a limit.
5. Use descriptive names, the existing style and automated formatting.
6. Explain reasons, constraints and unusual handling in comments; keep them current.
7. Preserve error context and tracebacks; never silently swallow failures or log secrets.
8. Validate external input, parameterize queries where applicable and keep credentials out of source control.
9. Run relevant checks; add meaningful regression tests for defects.
10. Review the diff before committing; remove debugging, temporary artifacts and unrelated changes; explain the purpose.
11. Update documentation when interfaces, configuration or startup changes.
12. Solve current needs with proportionate design; avoid speculative abstractions and dependencies.

## Local checks

Install `pip install -e ".[dev]"`, then run `ruff check lant tests tools examples` and `python -m pytest tests --cov=lant`. New infrastructure is formatted with ruff; CI lists its scoped format check. Existing scientific modules are not bulk reformatted in functional commits. Build both distributions with `python -m build`, then verify a wheel outside the checkout.
