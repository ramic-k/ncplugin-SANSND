# Repository Guidelines

Contributor notes for the `SANSND` NCrystal plugin. Keep edits focused, documented, and reproducible.

## Project Structure & Module Organization
- `src/`: C++ plugin sources and headers (`NCPhysicsModel.*`, `NCPluginFactory.*`, boilerplate registration).
- `data/`: Example `.ncmat` and auxiliary files consumed by the plugin when `@CUSTOM_SANSND` sections request external inputs.
- `testcode/`: Developer-only validation helpers (C++ helpers in `src/`, ctypes bindings in `python/`, and quick-run scripts in `scripts/`); `utils/bootstrap.sh` wires the build helper.
- Build/meta files: `CMakeLists.txt` (module build+install), `pyproject.toml` (scikit-build core), `ncplugin_name.txt` (authoritative plugin name).

## Build, Test, and Development Commands
- Prereqs: `ncrystal-config` on `PATH`, CMake ≥3.20, Python ≥3.6, NumPy for Python helpers.
- Recommended workflow:
  ```bash
  . testcode/utils/bootstrap.sh      # defines ncpluginbuild and checks prereqs
  ncpluginbuild                      # configure + build + install into test cache
  ncpluginbuild --force              # clean reconfigure if cache looks stale
  ncpluginbuild -DCMAKE_BUILD_TYPE=Release  # pass CMake opts through
  ```
- Raw CMake (if you prefer manual control):
  ```bash
  cmake -S . -B build -DCMAKE_INSTALL_PREFIX=./local -DNCPLUGIN_INSTALLSYMLINKS=ON
  cmake --build build
  cmake --install build
  ```
- Environment after `ncpluginbuild` is updated so plugin library, test scripts, and data are found via `PATH`, `PYTHONPATH`, and `NCRYSTAL_DATA_PATH`.

## Coding Style & Naming Conventions
- C++ code mirrors NCrystal core: 2-space indentation, braces on the same line, early guard/throw macros (`NCRYSTAL_THROW2`, `nc_assert_always`), and `NCPLUGIN_MSG` for diagnostics.
- Classes/types use PascalCase (`PhysicsModel`, `PluginFactory`); functions/methods use camelCase; constants/macros remain uppercase.
- Keep the plugin name consistent across `ncplugin_name.txt`, symbols (`NCPLUGIN_NAME_CSTR`), and `pyproject.toml`.
- Data files live under `data/` and are referenced relative to that root; avoid spaces/backup suffixes.

## Testing Guidelines
- No formal test runner; rely on the dev helpers installed by `ncpluginbuild`.
- Quick smoke checks:
  - `plotscat` or `plotxs` (after `ncpluginbuild`) to visualize sampled angles or cross sections.
  - Python hook sanity check:
    ```bash
    python - <<'PY'
    from testcode.python import cbindings
    print('hooks:', sorted(cbindings.hooks))
    PY
    ```
- When adding models/parameters, supply a minimal `.ncmat` example in `data/` and exercise it with the scripts to confirm the `@CUSTOM_SANSND` parsing path.

## Commit & Pull Request Guidelines
- Use short, present-tense commit subjects (mirroring history: “cleaning files in testcode/”, “Basic test of plugin including CI”) and group related changes.
- Before opening a PR: rerun `ncpluginbuild`, confirm scripts still execute, and note any required NCrystal version bumps.
- PR descriptions should mention motivation, key implementation notes, affected scripts/data, and include plots/log snippets if behaviour changed. Link issues or discussions when available.

## Security & Configuration Tips
- Never commit generated binaries or local cache from `testcode/utils/cache/`.
- Treat `ncrystal-config` output as the single source for include/lib paths; avoid hardcoding absolute locations.
- When distributing new data, cite the source and keep file names predictable (lowercase, dash/underscore separated).
