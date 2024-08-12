# # got this in the check() so i removed them from Imports: in the DESCRIPTION
# ❯ checking dependencies in R code ... NOTE
# Namespaces in Imports field not imported from:
#   ‘grDevices’ ‘graphics’ ‘methods’
# All declared Imports should be used.
#
#
# # what do all these mean:
# # <function>: no visible binding for global variable for '<variable>'
# # <function>: no visible global function definition for '<variable>
#
# Consider adding
# importFrom("grDevices", "dev.off", "pdf", "rgb")
# importFrom("graphics", "abline", "axis", "mtext", "par", "points",
#            "rect", "segments", "title")
# importFrom("methods", "new")
# importFrom("utils", "read.csv", "write.csv")
# to your NAMESPACE file (and ensure that your DESCRIPTION Imports field
#                         contains 'methods').


# TODOs
# 1
# if you want to store data in some raw, non-R-specific form and make it available to the user, put it in inst/extdata/.
# For example, readr and readxl each use this mechanism to provide a collection of delimited files and Excel workbooks, respectively.
# See Section 7.3.
# might want to move cytobands and rgdObject from sysdata to extdata
