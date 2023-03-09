
class ArticulationPoint:
    
    def APUtil(self, G, u): 
  
        #Count of children in current node  
        children = 0
    
        # Mark the current node as visited and print it 
        self.visited[u] = True
    
        # Initialize discovery time and low value 
        self.disc[u] = self.time
        self.low[u] = self.time
        self.time += 1
        max_low = -1
        #Recur for all the vertices adjacent to this vertex 
        for v in G.iterNeighbors(u): 
            # If v is not visited yet, then make it a child of u 
            # in DFS tree and recur for it 
            if self.visited[v] == False: 
                self.parent[v] = u 
                children += 1
                self.APUtil(G, v) 
            
                # Check if the subtree rooted with v has a connection to 
                # one of the ancestors of u
                max_low = max(self.low[v], max_low)
            
                self.low[u] = min(self.low[u], self.low[v])


            elif v != self.parent[u]:
                # Update low value of u for parent function calls     
                self.low[u] = min(self.low[u], self.disc[v]) 
            
        # u is an articulation point in following cases 
        # (1) u is root of DFS tree and has two or more chilren. 
        if self.parent[u] == -1 and children > 1: 
            self.ap.append(u)
        #(2) If u is not root and low value of one of its child is more 
        # than discovery value of u. 
        if self.parent[u] != -1 and max_low >= self.disc[u]:
            self.ap.append(u)

  
    #The function to do DFS traversal. It uses recursive APUtil() 
    def AP(self, G): 
   
        # Mark all the vertices as not visited  
        # and Initialize parent and visited,  
        # and ap(articulation point) arrays
        n = max([u for u in G.iterNodes()]) + 1
        self.visited = [False] * n
        self.disc = [float("Inf")] * n
        self.low = [float("Inf")] * n
        self.parent = [-1] * n
        self.ap = []
        #To store articulation points 
        
        # Call the recursive helper function 
        # to find articulation points 
        # in DFS tree rooted with vertex 'i'
        self.time = 0
        for u in G.iterNodes():
            if self.visited[u] == False: 
                self.APUtil(G, u) 

        return self.ap

    def APfrom(self, G, u): 
   
        # Mark all the vertices as not visited  
        # and Initialize parent and visited,  
        # and ap(articulation point) arrays
        n = max([u for u in G.iterNodes()]) + 1
        self.visited = [False] * n
        self.disc = [float("Inf")] * n
        self.low = [float("Inf")] * n
        self.parent = [-1] * n
        self.ap = []
        #To store articulation points 
        
        # Call the recursive helper function 
        # to find articulation points 
        # in DFS tree rooted with vertex 'i'
        self.time = 0
        self.APUtil(G, u) 

        return self.ap
