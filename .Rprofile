# source("renv/activate.R")
if (interactive() && Sys.getenv("RSTUDIO") == "") {
  source(file.path(Sys.getenv(
    if (.Platform$OS.type == "windows") "USERPROFILE" else "HOME"
  ), ".vscode-R", "init.R"))
  if ("httpgd" %in% .packages(all.available = TRUE)) {
    options(vsc.plot = FALSE, vsc.use_httpgd = FALSE)
    options(device = function(...) {
      httpgd::hgd(silent = TRUE)
      .vsc.browser(httpgd::hgd_url(history = FALSE), viewer = "Beside")
    })
  }
}

if (Sys.getenv("RSTUDIO")=="1"){
  Sys.setenv("CONDA_PREFIX"=NaN)
  Sys.setenv("RETICULATE_PYTHON_FALLBACK"=NaN)
}
# try(bigparallelr::set_blas_ncores(1), silent = TRUE)
