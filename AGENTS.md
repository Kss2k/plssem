# Repository Guidelines

## Project Structure & Module Organization

- `R/`: Package source (core modeling code, S4 classes/methods).
- `man/`: Generated Rd docs (from roxygen in `R/`).
- `tests/testthat/`: Unit/integration tests (`test_*.R`) and fixtures in `tests/testthat/test_data/`.
- `vignettes/`: Long-form documentation (built with `knitr`/`rmarkdown`).
- `docs/`: Generated pkgdown site (deployed to GitHub Pages). Avoid hand-editing.
- `data/` and `inst/`: Packaged datasets and installed files.

## Build, Test, and Development Commands

Run these from the repository root in an R session:

- `devtools::install()`: Install the package locally (useful before testing).
- `devtools::test(stop_on_failure = TRUE)`: Run `testthat` tests.
- `devtools::check(error_on = "warning")`: Full R CMD check (mirrors CI in `.github/workflows/checks.yml`).
- `devtools::document()`: Re-generate `NAMESPACE` and `man/` from roxygen comments.
- `pkgdown::build_site()`: Rebuild the website into `docs/`.

## Coding Style & Naming Conventions

- R style: 2-space indentation, use `<-`, and prefer explicit namespaces (e.g., `stats::median`).
- Naming: prefer `snake_case` for functions/variables (e.g., `pls_predict()`), `UpperCamelCase` for S4 classes.
- Public APIs must have roxygen2 docs in `R/` and corresponding examples kept lightweight.

## Testing Guidelines

- Framework: `testthat` (tests live in `tests/testthat/`).
- Keep tests deterministic and reasonably fast; when bootstrapping, use small `boot.R` values unless performance is the subject.
- Naming: `tests/testthat/test_<topic>.R` (e.g., `test_multilevel_random_slope.R`).

## Commit & Pull Request Guidelines

- Commits: follow existing history—short, imperative subjects (e.g., “Fix reindexing logic …”, “Update DESCRIPTION”).
- PRs: include a clear description, link relevant issues, and note user-facing changes.
- Before opening a PR: run `devtools::check()` locally and update docs via `devtools::document()` when APIs change.

## Configuration Tips

- CI runs on macOS/Windows/Linux; avoid OS-specific paths and rely on base R/temp directories.
- Do not add large binaries or generated artifacts beyond `docs/` (pkgdown output).
