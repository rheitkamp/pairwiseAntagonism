require(igraph)

distTable <- read.table("distanceVibrio.txt")
distMatrix <- as.matrix(distTable)
distances <- as.dist(distMatrix)
coordinates <- cmdscale(distances)

interactTable <- read.table("interactionsVibrio.txt")
interactMatrix <- as.matrix(interactTable)
antagMatrix <- t(interactMatrix)
interactGraph <- graph.adjacency(antagMatrix)

svg(file = 'antagonismVibrio.svg', width=15, length=13)
plot.igraph(interactGraph,
            layout = coordinates, 
            vertex.size = 4, 
            vertex.color = 'cyan', 
            vertex.frame.color = 'gray', 
            vertex.shape = 'circle', 
            vertex.label.cex = 0.2,
            vertex.label.color = 'gray',
            edge.color = 'gray', 
            edge.width = 0.01, 
            edge.arrow.size = 0.1, 
            edge.arrow.width = 0.4,
            edge.curved = TRUE, 
            rescale = TRUE,
            margin = 0.5
)
dev.off()