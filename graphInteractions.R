require(igraph)

distTable <- read.table("distanceAcinetobacter.txt")
distMatrix <- as.matrix(distTable)
distances <- as.dist(distMatrix)
coordinates <- cmdscale(distances)


interactTable <- read.table("interactionsAcinetobacter.txt")
interactMatrix <- as.matrix(interactTable)
antagMatrix <- t(interactMatrix)
interactGraph <- graph.adjacency(antagMatrix)
iG <- interactGraph

v1 <- '5711'
v2 <- c('334', '4052')
v3 <- c('4026', '5075', '87')
v4 <- c('51', '4490', '74', '52')

e1 <- 977
e2 <- c(47,241)
e3 <- c(224,830,1975)
e4 <- c(851,524,1509,904)

edges = c(e1, e2, e3, e4)

V(iG)$color <- 'yellow'
E(iG)$color <- NA
E(iG)$width <- 0.001

V(iG)$color <- ifelse(V(iG)$name == v1, 'red', 
                                 ifelse(V(iG)$name %in% v2, 'blue',
                                        ifelse(V(iG)$name %in% v3, 'green',
                                               ifelse(V(iG)$name %in% v4, 'purple',
                                                      'yellow'))))

E(iG)[e1]$color <- rgb(1,0,0) 
#E(iG)[e2]$color <- 'blue'
#E(iG)[e3]$color <- 'green'
#E(iG)[e4]$color <- 'purple'


E(iG)[e1]$width <- 2
#E(iG)[edges]$width <- 2

pdf(file = 'antagonismAcinetobacterUnary.pdf', bg = 'white')
plot.igraph(iG,
            layout = coordinates,
            vertex.size = 4, 
            vertex.color = 'yellow', 
            vertex.frame.color = 'gray', 
            vertex.shape = 'circle', 
            vertex.label.cex = 0.2,
            vertex.label.color = 'gray',
            edge.arrow.size = 2, 
            edge.arrow.width = 0.4,
            edge.curved = TRUE, 
            rescale = TRUE,
            margin = 0.5
)
dev.off()