a <- read.csv("AZ vs Mus(inac.) NV Low.csv")
a1 <- a$X
a2 <- a$logFC
b <- read.csv("EB vs Mus(inac.) NV Low.csv")
b1 <- b$X
b2 <- b$logFC
c <- read.csv("Int vs Mus(inac.) NV Low.csv")
c1 <- c$X
c2 <- c$logFC
d <- read.csv("Isp.1 vs Mus(inac.) NV Low.csv")
d1 <- d$X
d2 <- d$logFC
e <- read.csv("Mus vs Mus(inac.) NV Low.csv")
e1 <- e$X.
e2 <- e$logFC
f <- read.csv("N2 vs Mus(inac.) NV Low.csv")
f1 <- f$X
f2 <- f$logFC
g <- read.csv("Neu vs Mus(inac.) NV Low.csv")
g1 <- g$X
g2 <- g$logFC

housekeeping1 <- intersect(a1,b1)
housekeeping2 <- intersect(c1,d1)
housekeeping3 <- intersect(e1,f1)
housekeeping4 <- intersect(housekeeping1, housekeeping2)
housekeeping5 <- intersect(housekeeping3, housekeeping4)
housekeeping6 <- intersect(g1, housekeeping5)


y <- read.csv("Combined Low.csv")
y1 <- y[!duplicated(y$X), ]

write.table(y1, file="Combined Low wo dup.csv", sep = ",", col.names = NA, qmethod = "double")



y1 <- y[duplicated(df), ]

y1[y1$X %in% housekeeping6, ]



aa <- c(rep("A", 3), rep("B", 3), rep("C",2))
bb <- c(1,1,2,4,1,1,2,2)
df <-data.frame(aa,bb)
