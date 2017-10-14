Dataset1 <- function() {

    Dataset <<- read.csv("~/Documents/HEAT_LILLY_Heat.Stat.14_2016.csv")

    return(list(Dataset))
}
foo <- Dataset1()
 print(Dataset)


