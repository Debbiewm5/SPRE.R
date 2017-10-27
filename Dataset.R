Dataset <- function() {
    Dataset <<- read.csv("HEAT_LILLY_Heat.Stat.14_2016.csv")

    return(list(Dataset))
}
foo <- Dataset()
print(Dataset)


