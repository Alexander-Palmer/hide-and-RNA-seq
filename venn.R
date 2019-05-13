z1 <- read.csv('AZ vs N2 Top.csv')
z2 <- read.csv('EB vs N2 Top.csv')
z3 <- read.csv('Int vs N2 Top.csv')
z4 <- read.csv('Isp.1 vs N2 Top.csv')
z5 <- read.csv('Mus vs N2 Top.csv')
z6 <- read.csv('Neu vs N2 Top.csv')
z <- zz5

zz1 <- union(z1$X, z3$X)
zz2 <- union(z4$X, z5$X)
zz3 <- union(zz1, zz2)
zz4 <- union(z6$X, zz3)
zz5 <- union(y6$X, zz4)

y1 <- read.csv('AZ vs Mus(inac.) Top.csv')
y2 <- read.csv('EB vs Mus(inac.) Top.csv')
y3 <- read.csv('Int vs Mus(inac.) Top.csv')
y4 <- read.csv('Isp.1 vs Mus(inac.) Top.csv')
y5 <- read.csv('Mus vs Mus(inac.) Top.csv')
y6 <- read.csv('N2 vs Mus(inac.) Top.csv')
y7 <- read.csv('Neu vs Mus(inac.) Top.csv')
y <- yy6

yy1 <-union(y1$X, y2$X)
yy2 <-union(y3$X, y4$X)
yy3 <-union(y5$X, y6$X)
yy4 <-union(yy1, yy2)
yy5 <-union(yy3, yy4)
yy6 <-union(yy5, y7$X)


library(VennDiagram)

venn.plot <- venn.diagram(
  list("N2" = z, "Mus(Inac.)" = y), fill = c("cornflowerblue", "red"), cex = 1.1, 
  scaled=TRUE,
  "All sig miRNA between Mus.inac and N2.tiff"
)

