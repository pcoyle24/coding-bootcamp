def getNumberOfWays(N, Coins): 
  
    # Create the ways array to 1 plus the amount 
    # to stop overflow 
    ways = [0] * (N + 1); 
  
    # Set the first way to 1 because its 0 and 
    # there is 1 way to make 0 with 0 coins 
    ways[0] = 1; 
  
    # Go through all of the coins 
    for i in range(len(Coins)): 
  
        # Make a comparison to each index value 
        # of ways with the coin value. 
        for j in range(len(ways)): 
            if (Coins[i] <= j): 
                # Update the ways array 
                ways[j] += ways[(int)(j - Coins[i])];
                # print(ways) 
  
    # return the value at the Nth position 
    # of the ways array. 
    # print(ways) 
    return ways[N]; 
  

N = 7
S = [1, 2, 3]

numways = getNumberOfWays(N, S)
print(numways)

  


