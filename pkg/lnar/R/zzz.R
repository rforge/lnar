.onLoad <- function(libname, pkgname)
{
# library.dynam("lnar", pkgname, libname)
}
onUnload <- function(libpath)
{
  library.dynam.unload("lnar",libpath)
}

