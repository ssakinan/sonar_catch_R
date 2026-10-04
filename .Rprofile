# Project startup: the first time this project opens on a computer, open the
# analysis scripts in the editor. RStudio remembers open tabs after that.
local({
  flag <- file.path(".Rproj.user", "scripts_opened")
  if (interactive() && !file.exists(flag)) {
    setHook("rstudio.sessionInit", function(newSession) {
      if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
        for (f in sort(list.files(pattern = "^0[0-9]_.*\\.R$"))) rstudioapi::navigateToFile(f)
        dir.create(".Rproj.user", showWarnings = FALSE)
        file.create(flag)
      }
    }, action = "append")
  }
})
