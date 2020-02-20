/*********************************************
 * OPL 12.10.0.0 Model
 * Creation Date: 29 Ara 2019
 *********************************************/
 	using CPLEX;
 	int NHole = ...;
 	range Hole = 1..NHole;
 	int NBlock = ...;
 	range Block = 1..NBlock;
 	
 	
 	//arrays for x and y coordinates of holes and blocks
 	int xPoint[Hole] = ...;
 	int yPoint[Hole] = ...;
 	int xBlock[Block] = ...;
 	int yBlock[Block] = ...;
	int Width[Block] = ...;
	int Height[Block] = ...; 	
	
	
	int NxPoint = 50;
 	range xCoordinate = 1..NxPoint;
 	int NyPoint = 50;
 	range yCoordinate = 1..NyPoint;
 	
 	//This range is for U in our model. 
 	range uRange = 2..NHole;
	
	//array to check whether there exists a block at a point (x, y)
	int ContainsBlock[xCoordinate][yCoordinate];
	
	//This array holds whether there exist a between 2 holes or not. This array is same as Yij in our model. 
	int ContainsAWay[Hole][Hole];
	
	
	//Decision variables. U is same as in the model. Path is Xij in our model.
	dvar int u[2..NHole];
	dvar int Path[Hole][Hole] in 0..1;
	
	

	execute INIT_PARAMETERS{
	  
	  //This block checks all points(x = 1..50 and y = 1..50) to see whether they are part of a block.
	  for(var i in xCoordinate){  
	    for(var j in yCoordinate){
	      for(var b in Block){
	        if(i >= xBlock[b] && i <= xBlock[b] + Width[b]  && j >= yBlock[b] && j <= yBlock[b] + Height[b]){
	         	 ContainsBlock[i][j] = 1;
	         	 break;
	        }
	        ContainsBlock[i][j] = 0;
	      }
	    }
	  }
	  
	  
	 //This block checks whether a way exists between hole i and j. 
	 for(var i in Hole){
	   for(var j in Hole){
	     var smallX, largeX, smallY, largeY;
	     
	     //This calculation is necessary to initiliaze arrays
	     if(xPoint[i] < xPoint[j]){
	       smallX = xPoint[i];
	       largeX = xPoint[j];
	     }
	     else{
	       smallX = xPoint[j];
	       largeX = xPoint[i];
	     }
	     if(yPoint[i] < yPoint[j]){
	       smallY = yPoint[i];
	       largeY = yPoint[j];
	     }
	     else{
	       smallY = yPoint[j];
	       largeY = yPoint[i];
	     } 
   		
   		//This variable holds the number of ways going from a starting point to other points. 
   		var Ways = new Array(largeX - smallX + 1);
   		
   		for(var m = 0; m <= largeX - smallX; m++){
   			Ways[m] = new Array(largeY- smallY + 1);
   			for(var x = 0; x <= largeY - smallY; x++){
   			  Ways[m][x] = 0;
   			}
   		}   		
   		
        // If the starting point is inside a block, then there is no path from this point to other point. 
        if (ContainsBlock[xPoint[i]][yPoint[i]] == 1) {
            ContainsAWay[i][j] = 0;
            break;
        }

        // Number of ways of reaching the starting point = 1.
        else
        	Ways[0][0] = 1;


		var x1Larger = xPoint[i] >= xPoint[j];
		var y1Larger = yPoint[i] >= yPoint[j];
		
        // Filling the values for the first column
        for (var k = 1; k <= (largeX-smallX); k++) {
          //If x point of starting point is larger than end point's, we should go left, otherwise, we should go right.
          if((x1Larger && ContainsBlock[(xPoint[i] - k)][(yPoint[i])] == 0) || (!x1Larger && ContainsBlock[xPoint[i] + k][yPoint[i]] == 0))
           Ways[k][0] = Ways[k-1][0];
        }

        // Filling the values for the first row
        for (var l = 1; l <= (largeY - smallY); l++) {
          //If y point of starting point is larger than end point's, we should go down, otherwise, we should go up.
          if((y1Larger && ContainsBlock[xPoint[i]][yPoint[i] - l] == 0) || (!y1Larger && ContainsBlock[xPoint[i]][yPoint[i] + l] == 0))
           Ways[0][l] = Ways[0][l-1];
        }

        // Starting from cell(1,1) fill up the values
        // No. of ways of reaching cell[i][j] = cell[i - 1][j] + cell[i][j - 1]
        for (var m = 1; m <= (largeX - smallX); m++) {
            for (var n = 1; n <= (largeY - smallY); n++) {
              if ((x1Larger && y1Larger && ContainsBlock[(xPoint[i] - m)][(yPoint[i]) - n] == 0) || (!x1Larger && !y1Larger && ContainsBlock[xPoint[i] + m][yPoint[i] + n] == 0) ||
                 (x1Larger && !y1Larger && ContainsBlock[xPoint[i] - m][yPoint[i] + n] == 0) || (!x1Larger && y1Larger && ContainsBlock[xPoint[i] + m][yPoint[i] -n] == 0)){
                    Ways[m][n] = Ways[m - 1][n] + Ways[m][n - 1];
                } 
            }
        } 
        
        //If there exists at least one way to move between two holes, then it contains a way. 
        if(Ways[largeX - smallX][largeY - smallY] > 0)
        	ContainsAWay[i][j] = 1;
        else 
        	ContainsAWay[i][j] = 0;
        	
       }
       
     }               	
        	
   } 
   
   	//We should minimize Manhattan distance we move while we visiting all holes. Path[i][j] = 1 if we choose to move between holes i and j. 
   	dexpr float cost = sum(i in Hole, j in Hole)((abs(xPoint[i] - xPoint[j]) + abs(yPoint[i] - yPoint[j])) * Path[i][j]); 	    		 
	
	minimize 
		cost;
			
			
	subject to{	
		
		//Each hole should have an edge coming inside. However, it should not be more than one to prevent subroutines in the graph. 
		ctOneEdgeComing:
			forall(i in Hole)
				(sum(j in Hole) (Path[i][j])) == 1;
				
		//Each hole should have an edge going outside. However, it should not be more than one to prevent subroutines in the graph. 
		ctOneEdgeGoing:	
			forall(j in Hole)
		   		(sum(i in Hole) (Path[i][j])) == 1;
	
		
		//If there is no way between to holes, then there should not be an edge between two holes. 
		ctContainsAWay:
			forall(i in Hole, j in Hole)
				Path[i][j] <= ContainsAWay[i][j];	
				
				
		//Miller's method to prevent subroutine.
		ctNoSubroutine:
			forall(i in uRange, j in uRange){
				u[i] - u[j] + 50*(Path[i][j]) <= 49;	  
			}
			  			  
		//This constraint is a part of Miller's method. 	  
		ctURange:
			forall(i in uRange)
			  1 <= u[i] <= 49;
			  //2 <= u[i] <= 50;
		
		
	/* These constraints are redundant. However, we did not want to delete them. 		
		ct1:
		     forall(i in Hole)
			 	Path[i][i] == 0;		  
			  
		ct2:
			forall(i in Hole, j in Hole)	
				Path[i][j] == 0 && Path[j][i] == 0 || !(Path[i][j] == Path[j][i]);
				
		ct3:
			forall(j in Hole)
			  	Graph[j] <= sum(i in Hole)(Path[i][j] == 1 && Graph[i] == 1);
			  	
	*/
	} 

	//This is written to be able to display path using numbers of points. 
	execute display{
		var i = 1;
	  	writeln("Path is: ")
	  	write( i );
	  	for(var j in Hole){
  	  		if(Path[i][j] == 1){
  				i = j;
  		    	break;
       		}	  		     
  		}
  		while(i != 1){
  		  for(var j in Hole){
  		    if(Path[i][j] == 1){
  		    	write(" -> " + j);
  		     	i = j;
  		     	break;
       		}  		     
  		  }
  		}   
    } 
