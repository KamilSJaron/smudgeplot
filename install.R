library(devtools)

document()
test()
install.packages(".", repos = NULL, type="source")
