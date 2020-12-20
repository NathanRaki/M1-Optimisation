Prim = function(X, A) {
  visited = c(sample(X,1)) # Initialisation de la liste des sommets visités par un sommet pris au hasard
  mst = c() # Initialisation de notre Minimum Spanning Tree
  edges = which(A!=0, arr.ind=T) # Récupération des arêtes à partir de notre matrice d'adjacence
  
  while(length(visited) != length(X)) {
    possible = list() # Liste des arêtes possibles
    for (node in visited) {
      neighbours = edges[which(edges[,'row']==node), 'col'] # On récupère la liste des voisins d'un noeud
      neighbours = neighbours[which(!(neighbours %in% visited))] # On prend uniquement ceux qui ne sont pas visités
      for (neighbour in neighbours) {
        possible[[length(possible)+1]] = c(node, neighbour) # On les ajoute à la liste des arêtes possibles
      }
    }
    minval = Inf # On itialise un minimum à l'infini
    cursor = c() # Variable utilisée pour contenir notre arête minimale
    for (edge in possible) { # Pour chaque arête possible
      # Si sa valeur est inférieure au minimum stocké, on met le curseur dessus et on change le minimum
      if (A[edge[1], edge[2]] < minval) { 
        minval = A[edge[1], edge[2]]
        cursor = edge
      }
    }
    visited = append(visited, cursor[2]) # On ajoute notre nouveau noeud visité
    mst = append(mst, paste(cursor[1],'-',cursor[2], sep="")) # On ajoute l'arête possible minimale à notre arbre
  }
  return(list(visited, mst))
}

Ford_Bellman = function(vertices, adjacency, source)
{
  "
  INPUT:
    vertices : Liste des sommets du graphe
    adjacency : Matrice d'adjacence du graphe
    source : Sommet à partir duquel on calcule les plus courts chemins
  "
  dist = c()
  end = FALSE
  edges = which(adjacency!=0, arr.ind=T)
  
  print(edges)
  
  for (i in vertices) {
    dist[i] = Inf
  }
  dist[source] = 0
  
  while(!end)
  {
    end = TRUE
    for (i in vertices[-source])
    {
      pred = unname(edges[which(edges[,'col']==i),'row'])
      if (length(pred) > 0)
      {
        for (j in pred)
        {
          weight = dist[j] + adjacency[j,i]
          if (weight < dist[i])
          {
            dist[i] = weight
            end = FALSE
          }
        }
      }
    }
  }
  return(dist)
}

Ford_Fulkerson = function(vertices, adjacency, src, dst)
{
  maxFlow = 0 # Maximum flow (return value)
  edges = which(adjacency!=0, arr.ind=T) # Matrix of edges (row->col)
  
  flot = matrix(0, nrow=length(vertices), ncol=length(vertices)) # Flow matrix initialized to 0
  
  search = TRUE # Flag used for the 'GoTo Line 1'
  while(search) {
    
    mark = list() # Initializing a marking list
    mark[[src]] = c(-1, Inf, 1) # Initialize source mark (-1 for no predecessor, 1 for positive)
    S = c(src) # Used to store our visited vertices
    S_cache = c() # Vector used for a check at the end of the while(TRUE) part
    
    # First part of the pseudo-code, going from our source vertex to our destination
    stop=FALSE # Flag used for the 'GoTo Line 14'
    while(TRUE) # Used to repeat, we'll exit this loop with break statements
    {
      reloadS=FALSE # Flag used to come back here if the S vector has been changed and need to be considered by the for loop
      for (j in setdiff(X,S)) # For each node not in S
      { 
        iTab = intersect(union(edges[which(edges[,'col']==j), 'row'], edges[which(edges[,'row']==j), 'col']), S) # Store all of j's neighbours that are in S
        for (i in iTab) # For each vertex in iTab
        { 
          cij = adjacency[i,j] # cij is the edge capacity
          aij = flot[i,j] # aij is our current flow for the ij edge
          aji = flot[j,i] # aji is our current flow for the ji edge
          if ((cij - aij > 0) | (aji > 0)) # If there is something to do about the edge flow
          { 
            if (cij - aij > 0) # If positive
            { 
              aj = min( c( mark[[i]][2], cij-aij ) ) # j flow equals the minimum between i flow and cij-aij
              mark[[j]] = c(i, aj, 1) # We update j mark (i the used predecessor, aj the flow, 1 for positive)
            }
            else if (aji > 0) # If negative
            { 
              aj = min( c( mark[[i]][2], aji ) ) # j flow equals the minimum between i flow and aji
              mark[[j]] = c(i, aj, -1) # We update j mark (i the used predecessor, aj the flow, -1 for negative)
            }
            S = append(S, j) # We add j to the visited vertices
            reloadS = TRUE # Now we need to exit the for loop so it considers that we changed S
            if (j == dst) # If j is out destination
            {
              maxFlow = maxFlow + mark[[dst]][2] # We increment the maxFlow by j's flow
              stop=TRUE # Now we need to exit the while(TRUE)
            }
          }
          if(stop | reloadS) {break} # Flags IF statement
        }
        if(stop | reloadS) {break} # Flags IF statement
      }
      if(identical(S, S_cache)) {break} else {S_cache=S} # If S hasn't changed between 2 iterations, it means we couldn't reach the destination, so we exit
      if(stop) {break} # Flags IF statement
    }
    
    if (dst %in% S) # If our destination is in S
    {
      while(j != src)
      {
        if(mark[[j]][3] == 1) # If j mark is positive
          flot[mark[[j]][1],j] = flot[mark[[j]][1],j] + mark[[dst]][2] # Update our flow matrix
        
        else if (mark[[j]][3] == -1) # If j mark is negative
          flot[j,mark[[j]][1]] = flot[j,mark[[j]][1]] - mark[[dst]][2] # Update our flow matrix
        
        j = mark[[j]][1] # Change j to its predecessor
      }
    }
    else {search = FALSE} # If he's not, we end the algorithm
  }
  return(maxFlow)
}


