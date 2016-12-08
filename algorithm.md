#algorithm!

#airfoil section
1. Input the airfoil database, insert the coordinate to x and y direction.
2. Make the index goes clockwise (from TE to LE to TE in a clockwise manner).
3. Compute the length of each element.
4. Compute the angle beta on each node.
5. Compute control point coordinate of each element (middle of element). > dont wanna do this i guess!
6. Compute nx and ny on nodes, and on control points. nx and ny of nodes is the average of the control points neighbooring nx and ny.
7. Compute N1 and N2 on each node. Compute \eta first.

#solver
8. Looping process from i = 0; i < N
  1. Compute the right hand side of equation
    1. Compute G,i,N,- x dpsi,inf/dn^(N-)
    2. Compute G,i,N,+ x dpsi,inf/dn^(N+)
    3. Compute sigma,j=0,to,N-1 (H,i,j x psi,inf^j)
    4. Add 1, 2, and 3, to the right hand side coefficient.
  2. Compute the left hand side of equation
    1. Compute from,j,0,to,N-2 (G,i,j) -> (matrix index 0 to N - 2)
    2. Compute sigma,j,0,to,N-1 (H,i,j) -> (last index of matrix (N-1))
    3. Thus we will have the A matrix composed of N-1 x N-1 elements. 
    4. Thus we have Ax = y problem, which both A and y are known because they were already computed.
    5. Compute each element of vector x by using gauss method.
9. Finished. We know the from,0,to,N-2 (dpsi/dn) + psi,total^TE from each of the element of x vector.
10. Compute the psi on each node. (integrate from the derivated psi wrt n)
11. Compute velocity u and v
## stop first, after we reach point 11, continue! 
