# this is a stub package that points to its newer version deposited
# in http://www.bioconductor.org

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(paste("\n  The package qp has been replaced by the package qpgraph",
                              "that forms part of the Bioconductor project and can be",
                              "downloaded from http://www.bioconductor.org as follows:\n",
                              "source(\"http://www.bioconductor.org/biocLite.R\")",
                              "biocLite(\"qpgraph\")",
                              "\n  Please do 'help(qp)' to find out more information",
                              "about transiting from qp to qpgraph", sep="\n  "))
}
