library(testthat)
library(miRspongeR)

miR2Target <- system.file("extdata", "miR2Target.csv", package="miRspongeR")
miRTarget <- read.csv(miR2Target, header=TRUE, sep=",")
ExpDatacsv <- system.file("extdata", "ExpData.csv", package="miRspongeR")
ExpData <- read.csv(ExpDatacsv, header=FALSE, sep=",")

# miRHomology method
miRHomologyceRInt <- spongeMethod(miRTarget, method = "miRHomology")

# pc method
pcceRInt <- spongeMethod(miRTarget, ExpData, method = "pc")

# sppc method
sppcceRInt <- spongeMethod(miRTarget, ExpData, senscorcutoff = 0.1, method = "sppc")


test_that("Test spongeMethod", {
    expect_equal(spongeMethod(miRTarget, method = "miRHomology"), miRHomologyceRInt)
    expect_equal(spongeMethod(miRTarget, ExpData, method = "pc"), pcceRInt)
    expect_equal(spongeMethod(miRTarget, ExpData, senscorcutoff = 0.1, method = "sppc"), sppcceRInt)    
})
