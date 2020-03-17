function Build_Matrix_From_HSS(hss)
#Takes in an HSS matrix (structure) and returns the dense matrix A
#(nxn double). This code should be use for error checking/debugging only,
#in practice you should not construct the full matrix A since the amount of
#memory required to do so will be large. Uses feild labled Bu instead of Br
#
#INPUT:             hss    (structure) contains matrices (Us, Vs, Bs Rs and
#                           Ws) that compose the hss representation of the
#                           input matrix, as well as partition dimensions
#                           and if branch is a leaf.
#
#OUTPUT:            A       dense matrix A

    if isa(hss, Leaf) #if leaf node
        U = hss.U;
        V = hss.V;
        D = hss.D;

        D, U, V
    else
        #left recursive call
        Dl, Ul, Vl = Build_Matrix_From_HSS(hss.treeL);

        #right recursive call
        Dr, Ur, Vr = Build_Matrix_From_HSS(hss.treeR);

        #construct upper right and lower left block of current level
        urb = Ul*hss.Bu*Vr';
        llb = Ur*hss.Bl*Vl';

        #create D, U and V to pass to parent
        D = [Dl urb; llb Dr];
        U = [Ul*hss.Rl; Ur*hss.Rr];
        V = [Vl*hss.Wl; Vr*hss.Wr];

        D, U, V
    end
end
