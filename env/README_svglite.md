# svglite: separate R libraries (not in the conda envs)

Every figure is saved as `.png` and `.svg` (`save_plot()` / `png_svg()` in
`code/analysis/common.R`). For the SVGs to keep **editable text**, they are
written with `svglite`. svglite is not in either pinned conda env, and
installing into `tick_micro` risks a solver-driven upgrade (see
`tick_micro.yml`). So it is built from CRAN source into a separate R library
per env:

| env           | R     | library                    | contents                              |
|---------------|-------|----------------------------|---------------------------------------|
| `tick_micro`  | 4.3.3 | `env/r_lib_svglite/`        | svglite 2.1.3 (uses env's systemfonts 1.2.3) |
| `tick_report` | 4.3.3 | `env/r_lib_svglite_report/` | systemfonts 1.2.3 + svglite 2.1.3     |

`common.R` picks the library that matches the running R (`R.home()`), and
appends it **last** to `.libPaths()`, so it can never shadow an env package.
If the library is missing, SVGs fall back to Cairo `grDevices::svg()`, which
draws text as outlined glyphs.

The pinned version is svglite 2.1.3: 2.2.x needs systemfonts >= 1.3.0, newer
than tick_micro's.

## Rebuild

R inside the envs cannot reach CRAN directly, so download with curl. Build in a
directory that allows execution (`configure` scripts must run).

```bash
cd <project root>
mkdir -p env/.build_src env/.build_tmp
curl -s -o env/.build_src/svglite_2.1.3.tar.gz \
  https://cloud.r-project.org/src/contrib/Archive/svglite/svglite_2.1.3.tar.gz
curl -s -o env/.build_src/systemfonts_1.2.3.tar.gz \
  https://cloud.r-project.org/src/contrib/Archive/systemfonts/systemfonts_1.2.3.tar.gz

# tick_micro (the compiler is in the miniforge base bin)
E=/data/hb-austick/work/miniforge3/envs/tick_micro
mkdir -p env/r_lib_svglite
PATH=/data/hb-austick/work/miniforge3/bin:$E/bin:$PATH TMPDIR=$PWD/env/.build_tmp \
  $E/bin/R CMD INSTALL --library=env/r_lib_svglite env/.build_src/svglite_2.1.3.tar.gz

# tick_report (no pkg-config path for fontconfig/freetype: pass them explicitly)
R=/data/hb-austick/work/miniforge3/envs/tick_report
mkdir -p env/r_lib_svglite_report
/data/hb-austick/work/miniforge3/bin/conda run -n tick_report --cwd $PWD bash -c "
  export TMPDIR=$PWD/env/.build_tmp INCLUDE_DIR='$R/include/freetype2 -I$R/include' LIB_DIR=$R/lib
  R CMD INSTALL --library=env/r_lib_svglite_report env/.build_src/systemfonts_1.2.3.tar.gz &&
  R_LIBS=env/r_lib_svglite_report R CMD INSTALL --library=env/r_lib_svglite_report env/.build_src/svglite_2.1.3.tar.gz"

rm -rf env/.build_src env/.build_tmp
```
