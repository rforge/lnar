compmod <- function(cout,name="derivs")
{
setCMethod(name,signature(neq="integer",t="numeric",y="numeric",
           fout="numeric",vthetas="numeric"),cout$ccode,
           convention=".C",language="C",cxxargs="-fpic",PACKAGE="lnar")
}
