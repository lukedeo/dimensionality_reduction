isomap <- function(X, k = 6, d = 2, mode = "classical", verbose = FALSE, weighted = FALSE)
{
	if(verbose)
	{
		cat("Forming Neighborhoods...")
	}
	A <- fastPdist(X, X)
	if(!weighted)
	{
		A <- neighbor_graph(A, k, TRUE)
	}
	else
	{
		A <- weighted_neighbor_graph(A, k)
	}
	if(verbose)
	{
		cat("Done.\nForming graph construction...")
	}
	if(!weighted)
	{
		G <- graph.adjacency(A, mode = "undirected")
	}

	else
	{
		G <- graph.adjacency(A, mode = "undirected", weighted=TRUE)
	}
	
	if(verbose)
	{
		cat("Done.\nFinding geodesic distances...")
	}
	A <- shortest.paths(G, mode="all")
	# A <- shortest.paths(G, mode="all")

	if(verbose)
	{
		cat("Done.\nPeforming CMDS...")
	}
	A <- cmds(A, d)
	if(verbose)
	{
		cat("Done.")
	}
	A$description <- paste("ISOMAP, k = ", k, sep = "")
	return(A)
}