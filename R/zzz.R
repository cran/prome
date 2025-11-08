.onAttach <- function(libname, pkgname) {
  ver <- utils::packageVersion(pkgname)
  packageStartupMessage(
    paste0(pkgname, " version ", ver, " successfully loaded.")
  )
}
