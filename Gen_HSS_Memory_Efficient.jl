function Gen_Tree_Complete(minPart::Int64,N::Int64,d::Int64)
#function(minPart,N)
#This code is for testing purposes only.  Generates a tree
#structure that is not associated with any function.  Function takes in the
#size of the desired matrix and generates a tree by splitting this
#repeatedly in half on the left and right for a 'complete' tree.
#INPUT:             minPart     (Int64) minimum number of partitions
#                               type of tree that will be generated
#                   N           (Uint) number of points in desired
#                               function
#OUTPUT:            tree        (Union(Node_Spine,Leaf_Spine) contains split dimensions and depth at
#                               each level

    if N/2 < minPart # leaf
        depth = 0;
        leafL = Leaf_Spine(N,N,depth,d);

        d = d + N;
        leafL, depth, d
    else # node
        treeL, depth_l, d_l = Gen_Tree_Complete(minPart,convert(Int64,N/2),d);
        treeR, depth_r, d_r = Gen_Tree_Complete(minPart,convert(Int64,N/2),d_l);

        depth = 1 + max(depth_l,depth_r);
        tree = Node_Spine(treeL,treeR,N,N,depth,d);
        tree, depth, d_r
    end
end

function Gen_Tree_Mtree(minPart::Int64,N::Int64)
#function [tree] = Gen_TestTree(tree,N)
#This code is for testing purposes only. Artificually generates a tree
#structure that is not associated with any function.  Function takes in the
#number of desired nodes and generates a Worst case memory tree
# must have the line 'using Polynomials' in the main routine
#
#INPUT:
#                   N           (Int64) number of points in desired
#                               function
#                   minPart     (Int64) minimum number of partitions
#OUTPUT:            tree        (Tree) contains split dimensions at
#                               each level, as well as whether or not the
#                               node is a leaf


    # number of nodes
    n = 2*floor(N/minPart) -1; # n = N_N = N_L -1, and N_L = N/minPart.

    #depth of tree
    nodeRoots = roots(Poly([1-n,1,1]));
    d_exact  = nodeRoots[nodeRoots.>0]; #choose the positive root.

    #find maximum depth, d_max, of the tree we will generate
    d_max = convert(Int64,ceil(d_exact[1]));

    #Call the function initially on the root node of the tree.
    #   .<-
    #  / \
    # /\ /\
    #/\ /\/\
    tree = Gen_MTree_Right(N,minPart,d_max);
    tree
end

function Gen_MTree_Right(m,minPart,depth)
#Descend into right branch of the root of the current subtree (pictured below)
#   .
#  / \<-
# /\ /\
#/\ /\/\
#INPUT:
#                   m           (Int64) partition dimension of current
#                               node
#                   minPart     (Int64) minimum number of partitions
#                   depth       (Int64) depth of current node
#                   pB          (Bool) a flag that denotes whether the
#                               current node is on the prime branch
#OUTPUT:            tree        (Tree) contains split dimensions at
#                               each level, as well as whether or not the
#                               node is a leaf
    if !(m < 2*minPart)
    #Case 1: If we do have enough rows to split into 2 blocks of minimum
    #partition size split and recurse on both children.  If we do not, then do
    #not partition further; Return two leaf nodes.

        #max number of blocks of minimum partition size we can have on this left subtree
        numBlocks = convert(Int64,floor(m/minPart));
        #remaining depth of the left subtree.
        r_depth = numBlocks-1;

        #if we cannot generate a full left subtree - only generate children
        #corresponding to the the number of partitions we have left.
        if m <(depth +1)*minPart
            mL = r_depth*minPart;
        else
            mL = depth*minPart;
        end

        #Left Subtree Call
        treeL = Gen_MTree_Left(mL,minPart,depth);

        #Right Subtree Call
        #do we have enough rows to partition into two blocks of minumum partition size?
        if  m >= 2*minPart && m < 3*minPart
            #Case A: do we only have enough to split into two blocks and no
            #more?(blocks must be of at least minimum partition size and no
            #more than 2 times the minimum partition size)
            #Ex: minPart = 60
            #      /\170
            #   110  60
            #

            mL = m- minPart; #extra rows tacked onto left child
            leafL = Leaf_Spine(mL,mL,-1,-1);

            mR = minPart;
            leafR = Leaf_Spine(mR,mR,-1,-1);

            tree  = Node_Spine(leafL,leafR,m,m,-1,-1);
            tree
        elseif m<(depth+1)*minPart
            #Case B: do we have enough to split into more than two blocks, but cannot
            #generate a full left going subtree?
            #Ex:(N = 4096) Depth of full tree = 12. MinPart = 60.
            #(Though here I am only showing partition  dimensions for 3 levels.
            #Dots indicate a part of the tree not shown. )
            #   . . /\
            #  .     /\204
            #  . 144/\ 60
            #  .  60  84
            #  /\

            mR = m - (numBlocks-1)*minPart;

            treeR = Leaf_Spine(mR,mR,-1,-1);
            tree = Node_Spine(treeL,treeR,m,m,-1,-1);
            tree
        else #(depth+1)*minPart < m
            #Case C: We have enough rows to generate a full leftgoing subtree of depth d_max

            mR = m- depth*minPart;
            treeR = Gen_MTree_Right(mR,minPart,depth);

            tree = Node_Spine(treeL,treeR,m,m,-1,-1);
            tree
        end

    else #Case 2: we did not have enough rows to partition - label node as a leaf and return
        leaf = Leaf_Spine(m,m,-1,-1);
        leaf

    end

end

function Gen_MTree_Left(m,minPart,depth)
#Generate left branch of the root of the current subtree (pictured below)
#     .
#    / \
# ->/\ /\
#  /\ /\/\
#INPUT:
#                   m           (Int64) partition dimension of current
#                               node
#                   minPart     (Int64) minimum number of partitions
#                   depth       (Int64) depth of current node
#                   pB          (Bool) a flag that denotes whether the
#                               current node is on the prime branch
#OUTPUT:            tree        (Tree) contains split dimensions at
#                               each level, as well as whether or not the
#                               node is a leaf

    if m < 2*minPart  # if Leaf

        #pB = false;
        leafL = Leaf_Spine(m,m,-1,-1);

        leafL
    else #node
        depth -= 1;
        mL = m - minPart;
        treeL = Gen_MTree_Left(mL,minPart,depth);

        treeR = Leaf_Spine(minPart,minPart,-1,-1);

        tree = Node_Spine(treeL,treeR,m,m,-1,-1);
        tree
    end

end

function label_tree(tree,d)
#This function labels each node with the starting index (row/col value)
# of its corresponding diagonal block. Can also label the depth if the
# commented lines are uncommented
#
#INPUT:             tree    (Tree) contains partition dimensions for
#                           each subdivision of the input matrix, for each
#                           of which there are a left and right child
#                           structure
#                   d       (Int64) row and col index which corresponds to
#                           the first element in the current diagonal block
#                           at every node.
#OUTPUT:            tree    (Tree) same as input,  that contains valid
#                           diagonal (row/col) index at each node
#                           (and depth if commented lines are uncommented.
#
# Author: Kristen Lessel - Sept 2014

    if !isa(tree,Leaf_Spine) # if not a leaf node

        tree.d = d;
        tree.treeL,d = label_tree(tree.treeL,d);
        depth_l = tree.treeL.depth;

        tree.treeR,d = label_tree(tree.treeR,d);
        depth_r = tree.treeR.depth;

        tree.depth = 1 + max(depth_l,depth_r);
        tree, d
    else
        tree.depth = 0;
        tree.d = d;
        d = d + tree.m;
        tree, d
    end
end

function Gen_HSS_Memory_Efficient(tree,r)
#This is a memory efficient, two pass algorithm that takes in a Tree
#which is a partition tree for a given matrix which is described by the
#function, f, below and returns its corresponding Hss representation.
#Computation of U,R,V,W is done via deepest first post-ordering.  B
#matrices are computed by in descending order (starting at root and finishing
#at the child).  Relavant matrices are then multiplied and added
#from child to root in order to compute each B matrix.  (Bottom Up Routine
# to compute B).
#INPUT:             r       (Int64) largest allowable rank of hankel blocks -
#                           this determines the amount of compression, and
#                           should be compatible with your input tree.
#                           (if the maximum allowable rank is p, then the
#                           tree partitions should be no smaller than 3p)
#                           Corresponding singular values below this rank
#                           will be dropped.
#                   tree    (Tree) contains partition dimensions for
#                           each subdivision of the input matrix, for each
#                           of which there are a left and right child
#                           structure
#OUTPUT:            hss    (Hss) contains matrices (Us, Vs, Bs Rs and
#                           Ws) that compose the hss representation of the
#                           input matrix, in addition to the partition
#                           dimensions
#                   peakMem (Array{Int64}(1)) first entry is a memory counter
#                           and second entry is the maximum memory
# Author: Kristen Lessel - June 2015
    N = tree.m;
    pm_count = 0;
    pm_max = 0;
    peakMem = [pm_count,pm_max];
    hss, dummy1, dummy2, peakMem = HSS_Basis_TranslationOp(tree,N,r,peakMem);
    row_idx = 1;
    col_idx = 1;
    hss, peakMem = HSS_Expansion_Coeffs_Diag(hss,row_idx,col_idx,N,peakMem);
    hss, peakMem
end

function HSS_Basis_TranslationOp(tree::Tree,N::Int64,r::Int64,peakMem::Array{Int64,1})
#Computes U's, V's, R's, W's and D's of HSS structure, and stores these
#heirarchically
#INPUT:             N       (Int64) grid size
#                   rank    (Int64) largest allowable rank of hankel blocks -
#                           this determines the amount of compression, and
#                           should be compatible with your input tree.
#                           (if the maximum allowable rank is p, then the
#                           tree partitions should be no smaller than 3p)
#                           Corresponding singular values below this rank
#                           will be dropped.
#                   tree    (Union(Node_Spine,Leaf_Spine))
#                           contains partition dimensions for each
#                           subdivision of the input matrix, for each
#                           of which there are a left and right child
#                           structure
#OUTPUT:            hss     (Union(Node,Leaf)) contains matrices
#                           (U, V, R and W) that compose the
#                           hss representation of the input matrix, in
#                           addition to the partition dimensions
#                   rH      (Array{Float64,2}) row hankel block with 'current' U
#                           removed
#                   cH      (Array{Float64,2}) column hankel block with 'current' V'
#                           removed
  if isa(tree,Leaf_Spine)
      m = tree.m;
      d = tree.d;

      #generate indices for currrent diagonal block
      m_idx = collect(d:(d+m)-1);
      #generate indices for current off diagonal hankel block
      butm_idx = [1:m_idx[1]-1; m_idx[end]+1:N];

      #function call to generate lowest level row and column hankel blocks
      rH2 = f(m_idx,butm_idx,N);
      cH2= f(butm_idx,m_idx,N);
      #keep track of peak memory
      count = 2*length(m_idx)*length(butm_idx);
      peakMem[1] = peakMem[1]+count;
      peakMem[2] = max(peakMem[1],peakMem[2]);

      #take svds of upper row/left column and lower row/right column hankel blocks
      Ur, Rr, Vr = svd_ranktol(rH2,r);
      Uc, Rc, Vc = svd_ranktol(cH2,r);
      #keep track of peak memory
      count = length(Ur)+ length(Vc)+length(Rr)+length(Vr)+length(Uc)+length(Rc);
      peakMem[1] = peakMem[1]+count;
      peakMem[2] = max(peakMem[1],peakMem[2]);

      #keep track of peak memory
      count = length(rH2)+length(cH2);
      peakMem[1] = peakMem[1]-count;
      peakMem[2] = max(peakMem[1],peakMem[2]);

      rH2 = Array{Float64}(undef,0);
      cH2 = Array{Float64}(undef,0);

      #D = f(m_idx,m_idx,N)
      leaf = Leaf(m,m,d,tree.depth,Ur,Vc,f(m_idx,m_idx,N));

      rH = diagm(Rr)*Vr';
      cH = Uc*diagm(Rc);

      count = length(Rr)+length(Rc);
      peakMem[1] = peakMem[1]-count;
      peakMem[2] = max(peakMem[1],peakMem[2]);

      leaf, rH, cH, peakMem

  else #If not a Leaf

      #Descend into the tree, depth first
      if tree.treeL.depth >= tree.treeR.depth #left node deeper than right
          #left recursive call
          hssL, upperRow, leftCol, peakMem = HSS_Basis_TranslationOp(tree.treeL,N,r,peakMem);

          #right recursive call
          hssR, lowerRow, rightCol, peakMem = HSS_Basis_TranslationOp(tree.treeR,N,r,peakMem);
      else #right node deeper than left
          #right recursive call
          hssR, lowerRow, rightCol, peakMem = HSS_Basis_TranslationOp(tree.treeR,N,r,peakMem);

          #left recursive call
          hssL, upperRow, leftCol, peakMem = HSS_Basis_TranslationOp(tree.treeL,N,r,peakMem);
      end

      #starting index (row/col value) of its current diagonal block.
      d = tree.d;

      #indexes of rol/col hankel blocks excluding corresponding part of diagonal block
      uR_idx = [1:(d-1); (d+hssR.m):(N-hssL.m)];
      lR_idx = [1:(d-1); (d+hssL.m):(N-hssR.m)];

      rH_top = upperRow[:, uR_idx];
      rH_bottom = lowerRow[:, lR_idx];
      cH_left = leftCol[uR_idx,:];
      cH_right = rightCol[lR_idx,:];
      #keep track of peak memory
      count = length(rH_top)+ length(rH_bottom)+length(cH_left)+length(cH_right);
      peakMem[1] = peakMem[1]+count;
      peakMem[2] = max(peakMem[1],peakMem[2]);

      #merge remainder of Hankel blocks from children
      rH2 = [rH_top; rH_bottom];
      cH2 = [cH_left cH_right];

      count = length(rH_top)+ length(rH_bottom)+length(cH_left)+length(cH_right);
      peakMem[1] = peakMem[1]-count;
      peakMem[2] = max(peakMem[1],peakMem[2]);

      #take svd of remaining portions of blocks to determine Rs and Ws
      Ur, Sr, Vr = svd_ranktol(rH2,r);
      Uc, Sc, Vc = svd_ranktol(cH2,r);
      #keep track of peak memory
      count = length(Ur) +length(Sr) +length(Vr) + length(Uc) +length(Sc) +length(Vc)
      peakMem[1] = peakMem[1]+count;
      peakMem[2] = max(peakMem[1],peakMem[2]);

      #keep track of peak memory %just added 10/8/15
      count = length(rH2)+length(cH2);
      peakMem[1] = peakMem[1]-count;
      peakMem[2] = max(peakMem[1],peakMem[2]);

      rH2 = Array{Float64}(undef,0);
      cH2 = Array{Float64}(undef,0);

      #partition Us to get R's, partition V's to get Ws
      Rl = Ur[1:size(upperRow,1),:];
      Rr = Ur[size(upperRow,1)+1:end,:];
      Wl = Vc[1:size(leftCol,2),:];
      Wr = Vc[size(leftCol,2)+1:end,:];
      #keep track of peak memory
      count = length(Rl)+length(Rr)+length(Wl)+length(Wr);
      peakMem[1] = peakMem[1]+count;
      peakMem[2] = max(peakMem[1],peakMem[2]);

      count = length(Ur) +length(Vc);
      peakMem[1] = peakMem[1]-count;
      peakMem[2] = max(peakMem[1],peakMem[2]);


      hss = Node(hssL,hssR,tree.m,tree.m,tree.depth,d,Array{Float64,2}(undef,0,0),Array{Float64,2}(undef,0,0),Rl,Rr,Wl,Wr);
      #return relavant portions of hankel blocks
      rH = diagm(Sr)*Vr';
      cH = Uc*diagm(Sc);

      #keep track of peak memory
      count = 2*r;
      peakMem[1] = peakMem[1]-count;
      peakMem[2] = max(peakMem[1],peakMem[2]);


      hss, rH, cH, peakMem
  end
end


function HSS_Expansion_Coeffs_Diag(hss::Hss,row_idx::Int64,col_idx::Int64,N::Int64,peakMem::Array{Int64,1})
#Recursively computes Expansion Coefficients (B's) of HSS structure for the
#upper left and lower right block of the current node
#   -------
#  \ x \   \
#  \-------\
#  \   \ x \
#   -------
#INPUT:             hss     (Union(Node,Leaf)) branch of current node,
#                   row_idx (Int64) row index of current node
#                   col_idx (Int64) col index of current node
#                   N       (Int64) grid size
#OUTPUT:            hss     (Union(Node,Leaf)) contains matrices (U, V, B R and
#                           W) that compose the HSS representation of the
#                           input tree for a corresponding matrix, as well
#                           as corresponding partition dimensions and
#                           depth for each node

    if isa(hss.treeL,Leaf) && isa(hss.treeR,Leaf) #both nodes are leaves
        Au = Array{Float64}(undef,2);
        Bu = Array{Float64}(undef,2);
        Al = Array{Float64}(undef,2);
        Bl = Array{Float64}(undef,2);

        col_idx = col_idx + hss.treeL.m;
        m_idx = collect(row_idx:(row_idx+hss.treeL.m-1));
        butm_idx = collect(col_idx:(col_idx+hss.treeR.m-1));

        Au = f(m_idx,butm_idx,N);
        Bu = hss.treeL.U'*Au*hss.treeR.V;
        Au =  Array{Float64}(undef,0);
        #keep track of peak memory
        count = length(m_idx)*length(butm_idx) +  length(Bu);
        peakMem[1] = peakMem[1]+count;
        peakMem[2] = max(peakMem[1],peakMem[2]);

        count = length(m_idx)*length(butm_idx);
        peakMem[1] = peakMem[1]-count;
        peakMem[2] = max(peakMem[1],peakMem[2]);


        Al = f(butm_idx,m_idx,N);
        Bl = hss.treeR.U'*Al*hss.treeL.V;
        Al =  Array{Float64}(undef,0);
        #keep track of peak memory
        count = length(m_idx)*length(butm_idx) +  length(Bl);
        peakMem[1] = peakMem[1]+count;
        peakMem[2] = max(peakMem[1],peakMem[2]);

        count = length(m_idx)*length(butm_idx);
        peakMem[1] = peakMem[1]-count;
        peakMem[2] = max(peakMem[1],peakMem[2]);

        hss.Bu = Bu;
        hss.Bl = Bl;
        hss, peakMem

    elseif !isa(hss.treeL,Leaf) && !isa(hss.treeR, Leaf) #neither node is a leaf

        #Calculate B's for off-diagonal blocks
        col_idx = col_idx + hss.treeL.m;
        Bu, Bl, peakMem = HSS_Expansion_Coeffs_OffDiag(hss.treeL,hss.treeR,row_idx,col_idx,N,peakMem);
        col_idx = col_idx - hss.treeL.m;

        hss.Bl = Bl;
        hss.Bu = Bu;

        #Calculate B's for upper left diagonal block
        hssL, peakMem = HSS_Expansion_Coeffs_Diag(hss.treeL,row_idx,col_idx,N,peakMem);

        #Calculate B's for lower right diagonal block
        row_idx = row_idx + hss.treeL.m;
        col_idx = col_idx + hss.treeL.m;
        hssR, peakMem = HSS_Expansion_Coeffs_Diag(hss.treeR,row_idx,col_idx,N,peakMem);

        hss.treeL = hssL;
        hss.treeR = hssR;
        hss, peakMem

    elseif isa(hss.treeL,Leaf) && !isa(hss.treeR, Leaf) #left node is a leaf, right node is not
        bu_l = Array{Float64}(undef,2);
        bu_r = Array{Float64}(undef,2);
        bl_l = Array{Float64}(undef,2);
        bl_r = Array{Float64}(undef,2);

        #Calculate B's for off-diagonal blocks

        #left block of upper right corner and upper block of lower left corner
        #(these are computed in the same pass)
        #    ---------------
        #   \       \   \   \
        #   \       \ x \   \
        #   \       \   \   \
        #   \---------------\
        #   \   x   \   \   \
        #   \-------\---\---\
        #   \       \   \   \
        #    ---------------
        #           (9)

        col_idx = col_idx + hss.treeL.m;
        Bu_l, Bl_l, peakMem = HSS_Expansion_Coeffs_OffDiag(hss.treeL,hss.treeR.treeL,row_idx,col_idx,N,peakMem);
        bu_l = Bu_l*hss.treeR.Wl;
        bl_l = hss.treeR.Rl'*Bl_l;
        #keep track of peak memory
        count = length(bu_l)+length(bl_l);
        peakMem[1] = peakMem[1]+count;
        peakMem[2] = max(peakMem[1],peakMem[2]);

        count = length(Bu_l) + length(Bl_l); #gives number of elements
        peakMem[1] = peakMem[1]-count;
        peakMem[2] = max(peakMem[1],peakMem[2]);

        Bu_l =  Array{Float64}(undef,0);
        Bl_l =  Array{Float64}(undef,0);

        #right block of upper right corner and lower block of lower left corner
        #(these are computed in the same pass)
        #    ---------------
        #   \       \   \   \
        #   \       \   \ x \
        #   \       \   \   \
        #   \---------------\
        #   \       \   \   \
        #   \-------\---\---\
        #   \   x   \   \   \
        #    ---------------
        #           (10)

        col_idx = col_idx + hss.treeR.treeL.m;
        Bu_r, Bl_r, peakMem = HSS_Expansion_Coeffs_OffDiag(hss.treeL,hss.treeR.treeR,row_idx,col_idx,N,peakMem);
        col_idx = col_idx - hss.treeL.m - hss.treeR.treeL.m;
        bu_r = Bu_r*hss.treeR.Wr;
        bl_r = hss.treeR.Rr'*Bl_r;
        #keep track of peak memory
        count = length(bu_r)+length(bl_r);
        peakMem[1] = peakMem[1]+count;
        peakMem[2] = max(peakMem[1],peakMem[2]);

        count = length(Bu_r)+length(Bl_r);
        peakMem[1] = peakMem[1]-count;
        peakMem[2] = max(peakMem[1],peakMem[2]);

        Bu_r = Array{Float64}(undef,0);
        Bl_r = Array{Float64}(undef,0);

        hss.Bl = bl_l + bl_r;
        hss.Bu = bu_l + bu_r;

        count = length(bl_l)+length(bl_r)+length(bu_l)+length(bu_r);
        peakMem[1] = peakMem[1]-count;
        peakMem[2] = max(peakMem[1],peakMem[2]);

        bl_l = Array{Float64}(undef,0);
        bl_r = Array{Float64}(undef,0);
        bu_l = Array{Float64}(undef,0);
        bu_r = Array{Float64}(undef,0);

        #Calculate B's for lower right diagonal block upper right
        #block of lower right corner and lower left block of lower
        #right corner
        #    ---------------
        #   \       \   \   \
        #   \       \   \   \
        #   \       \   \   \
        #   \---------------\
        #   \       \   \ x \
        #   \-------\---\---\
        #   \       \ x \   \
        #    ---------------
        #           (11)

        row_idx = row_idx + hss.treeL.m;
        col_idx = col_idx + hss.treeL.m;
        hssR, peakMem = HSS_Expansion_Coeffs_Diag(hss.treeR,row_idx,col_idx,N,peakMem);

        hss.treeR = hssR;
        hss, peakMem

    elseif !isa(hss.treeL,Leaf) && isa(hss.treeR, Leaf) #right node is a leaf, left node is not
        bu_l = Array{Float64}(undef,2);
        bu_r = Array{Float64}(undef,2);
        bl_l = Array{Float64}(undef,2);
        bl_r = Array{Float64}(undef,2);

        #Calculate B's for off-diagonal blocks

        #upper block of upper right corner and left block of lower left corner
        #(these are computed in the same pass)
        #    ---------------
        #   \   \   \   x   \
        #   \---\---\-------\
        #   \   \   \       \
        #   \---------------\
        #   \   \   \       \
        #   \ x \   \       \
        #   \   \   \       \
        #    ---------------
        #           (12)

        col_idx = col_idx + hss.treeL.m;
        Bu_l, Bl_l, peakMem = HSS_Expansion_Coeffs_OffDiag(hss.treeL.treeL,hss.treeR,row_idx,col_idx,N,peakMem);
        bu_l = hss.treeL.Rl'*Bu_l;
        bl_l = Bl_l*hss.treeL.Wl;
        #keep track of peak memory
        count = length(bu_l)+length(bl_l);
        peakMem[1] = peakMem[1]+count;
        peakMem[2] = max(peakMem[1],peakMem[2]);

        count = length(Bu_l)+length(Bl_l);
        peakMem[1] = peakMem[1]-count;
        peakMem[2] = max(peakMem[1],peakMem[2]);

        Bu_l = Array{Float64}(undef,0);
        Bl_l = Array{Float64}(undef,0);

        #lower block of upper right corner and right block of lower left corner
        #(these are computed in the same pass)
        #    ---------------
        #   \   \   \       \
        #   \---\---\-------\
        #   \   \   \   x   \
        #   \---------------\
        #   \   \   \       \
        #   \   \ x \       \
        #   \   \   \       \
        #    ---------------
        #           (13)

        row_idx  = row_idx + hss.treeL.treeL.m;
        Bu_r, Bl_r, peakMem = HSS_Expansion_Coeffs_OffDiag(hss.treeL.treeR,hss.treeR,row_idx,col_idx,N,peakMem);
        row_idx = row_idx - hss.treeL.treeL.m;
        col_idx = col_idx - hss.treeL.m;
        bu_r = hss.treeL.Rr'*Bu_r;
        bl_r = Bl_r*hss.treeL.Wr;
        #keep track of peak memory
        count = length(bu_r)+length(bl_r);
        peakMem[1] = peakMem[1]+count;
        peakMem[2] = max(peakMem[1],peakMem[2]);

        count = length(Bu_r)+length(Bl_r);
        peakMem[1] = peakMem[1]-count;
        peakMem[2] = max(peakMem[1],peakMem[2]);

        Bu_r = Array{Float64}(undef,0);
        Bl_r = Array{Float64}(undef,0);

        hss.Bl = bl_l + bl_r;
        hss.Bu = bu_l + bu_r;

        #keep track of peak memory
        count = length(bl_l)+length(bl_r)+length(bu_l)+length(bu_r);
        peakMem[1] = peakMem[1]-count;
        peakMem[2] = max(peakMem[1],peakMem[2]);

        bl_l = Array{Float64}(undef,0);
        bl_r = Array{Float64}(undef,0);
        bu_l = Array{Float64}(undef,0);
        bu_r = Array{Float64}(undef,0);

        #Calculate B's for upper right diagonal block

        #lower left block of upper left corner and upper right block of upper
        # left corner
        #    ---------------
        #   \   \ x \       \
        #   \---\---\-------\
        #   \ x \   \       \
        #   \---------------\
        #   \   \   \       \
        #   \   \   \       \
        #   \   \   \       \
        #    ---------------
        #           (14)

        hssL, peakMem = HSS_Expansion_Coeffs_Diag(hss.treeL,row_idx,col_idx,N,peakMem);

        hss.treeL = hssL;
        hss, peakMem
    end
end

function HSS_Expansion_Coeffs_OffDiag(treeL::Hss,treeR::Hss,row_idx::Int64,col_idx::Int64,N::Int64,peakMem::Array{Int64,1})
#Computes Expansion Coefficients (B's) of HSS structure for the upper right
#and lower left block of the current node
#   -------
#  \   \ x \
#  \-------\
#  \ x \   \
#   -------
#INPUT:             treeL   (Union(Node,Leaf)) left branch of current node
#                   treeR   (Union(Node,Leaf)) right branch of current node
#                   row_idx (Int64) row index of current node
#                   col_idx (Int64) col index of current node
#                   N       (Int64) grid size
#OUTPUT:            Bu      (Array{Float62,2}) Expansion coefficient matrix B at the
#                           current level which corresponds to the upper
#                           right block
#                   Bl      (Array{Float62,2}) Expansion coefficient matrix B at the
#                           current level which corresponds to the lower
#                           left block

    Au = Array{Float64,2}(undef);
    Bu = Array{Float64,2}(undef);
    Al = Array{Float64,2}(undef);
    Bl = Array{Float64,2}(undef);

    bu_l = Array{Float64,2}(undef);
    bu_r = Array{Float64,2}(undef);
    bl_l = Array{Float64,2}(undef);
    bl_r = Array{Float64,2}(undef);
    if isa(treeL, Leaf) && isa(treeR, Leaf) #If node is a leaf
        m_idx = collect((row_idx):(row_idx+treeL.m-1));
        butm_idx = collect(col_idx:(col_idx+treeR.m-1));

        Au = f(m_idx,butm_idx,N);
        Bu = treeL.U'*Au*treeR.V;
        #keep track of peak memory
        count = length(Au)+length(Bu);
        peakMem[1] = peakMem[1]+count;
        peakMem[2] = max(peakMem[1],peakMem[2]);

        count = length(Au);
        peakMem[1] = peakMem[1]-count;
        peakMem[2] = max(peakMem[1],peakMem[2]);

        Au = Array{Float64}(undef,0);

        Al = f(butm_idx,m_idx,N);
        Bl = treeR.U'*Al*treeL.V;
        #keep track of peak memory
        count = length(Al)+length(Bl);
        peakMem[1] = peakMem[1]+count;
        peakMem[2] = max(peakMem[1],peakMem[2]);

        count = length(Al);
        peakMem[1] = peakMem[1]-count;
        peakMem[2] = max(peakMem[1],peakMem[2]);

        Al = Array{Float64}(undef,0);

        Bu, Bl, peakMem

    elseif !isa(treeL, Leaf) && !isa(treeR,Leaf)  #If neither node is a leaf
        #COMPUTE EXPANSION COEFFICIENTS

        #upper left corner of both upper right and lower left blocks
        #(these are computed at the same time. One block will be of
        #dimension MxN, the other will always be of dimension NxM)
        #
        #    ---------------
        #   \       \ x \<--\-- MxN
        #   \        -------\
        #   \       \   \   \
        #   \---------------\
        #   \ x \<--\-------\-- NxM
        #   \-------\       \
        #   \   \   \       \
        #    ---------------
        #           (1)

        Bu_ll, Bl_ll, peakMem = HSS_Expansion_Coeffs_OffDiag(treeL.treeL,treeR.treeL,row_idx,col_idx,N,peakMem);
        bu_ll = treeL.Rl'*Bu_ll*treeR.Wl;
        bl_ll = treeR.Rl'*Bl_ll*treeL.Wl;
        #keep track of peak memory
        count = length(bu_ll)+length(bl_ll);
        peakMem[1] = peakMem[1]+count;
        peakMem[2] = max(peakMem[1],peakMem[2]);

        count = length(Bu_ll)+length(Bl_ll);
        peakMem[1] = peakMem[1]-count;
        peakMem[2] = max(peakMem[1],peakMem[2]);

        Bu_ll = Array{Float64}(undef,0);
        Bl_ll = Array{Float64}(undef,0);

        #Lower left corner lower left block, and also the upper right corner of
        #the upper right block
        #    ---------------
        #   \       \   \ x \
        #   \       \---\---\
        #   \       \   \   \
        #   \---------------\
        #   \   \   \       \
        #   \-------\       \
        #   \ x \   \       \
        #    ---------------
        #           (2)

        col_idx = col_idx +treeR.treeL.m;
        Bu_lr, Bl_lr, peakMem = HSS_Expansion_Coeffs_OffDiag(treeL.treeL,treeR.treeR,row_idx,col_idx,N,peakMem);
        bu_lr = treeL.Rl'*Bu_lr*treeR.Wr;
        bl_lr = treeR.Rr'*Bl_lr*treeL.Wl;
        #keep track of peak memory
        count = length(bu_lr)+length(bl_lr);
        peakMem[1] = peakMem[1]+count;
        peakMem[2] = max(peakMem[1],peakMem[2]);

        count = length(Bu_lr)+length(Bl_lr);
        peakMem[1] = peakMem[1]-count;
        peakMem[2] = max(peakMem[1],peakMem[2]);

        Bu_lr = Array{Float64}(undef,0);
        Bl_lr = Array{Float64}(undef,0);

        #upper right corner of lower block, lower left corner of upper block
        #    ---------------
        #   \       \   \   \
        #   \       \---\---\
        #   \       \ x \   \
        #   \---------------\
        #   \   \ x \       \
        #   \-------\       \
        #   \   \   \       \
        #    ---------------
        #           (3)

        col_idx = col_idx -treeR.treeL.m;
        row_idx = row_idx +treeL.treeL.m;
        Bu_rl, Bl_rl, peakMem = HSS_Expansion_Coeffs_OffDiag(treeL.treeR,treeR.treeL,row_idx,col_idx,N,peakMem);
        bu_rl = treeL.Rr'*Bu_rl*treeR.Wl;
        bl_rl = treeR.Rl'*Bl_rl*treeL.Wr;
        #keep track of peak memory
        count = length(bu_rl)+length(bl_rl);
        peakMem[1] = peakMem[1]+count;
        peakMem[2] = max(peakMem[1],peakMem[2]);

        count = length(Bu_rl)+length(Bl_rl);
        peakMem[1] = peakMem[1]-count;
        peakMem[2] = max(peakMem[1],peakMem[2]);

        Bu_rl = Array{Float64}(undef,0);
        Bl_rl = Array{Float64}(undef,0);

        #lower right corner of both lower left and upper left blocks
        #    ---------------
        #   \       \   \   \
        #   \       \---\---\
        #   \       \   \ x \
        #   \---------------\
        #   \   \   \       \
        #   \-------\       \
        #   \   \ x \       \
        #    ---------------
        #           (4)

        col_idx = col_idx +treeR.treeL.m;
        Bu_rr, Bl_rr, peakMem = HSS_Expansion_Coeffs_OffDiag(treeL.treeR,treeR.treeR,row_idx,col_idx,N,peakMem);
        bu_rr = treeL.Rr'*Bu_rr*treeR.Wr;
        bl_rr = treeR.Rr'*Bl_rr*treeL.Wr;
        #keep track of peak memory
        count = length(bu_rr)+length(bl_rr);
        peakMem[1] = peakMem[1]+count;
        peakMem[2] = max(peakMem[1],peakMem[2]);

        count = length(Bu_rr)+length(Bl_rr);
        peakMem[1] = peakMem[1]-count;
        peakMem[2] = max(peakMem[1],peakMem[2]);

        Bu_rr = Array{Float64}(undef,0);
        Bl_rr = Array{Float64}(undef,0);

        #Add b's to form the B's
        Bu = bu_ll+bu_rl+bu_lr+bu_rr;
        Bl = bl_ll+bl_rl+bl_lr+bl_rr;
        #keep track of peak memory
        count = length(Bu)+length(Bl);
        peakMem[1] = peakMem[1]+count;
        peakMem[2] = max(peakMem[1],peakMem[2]);

        count = length(bu_ll)+length(bu_rl)+length(bu_lr)+length(bu_rr)+length(bl_ll)+length(bl_rl)+length(bl_lr)+length(bl_rr);
        peakMem[1] = peakMem[1]-count;
        peakMem[2] = max(peakMem[1],peakMem[2]);

        bu_ll = Array{Float64}(undef,0);
        bu_rl = Array{Float64}(undef,0);
        bu_lr = Array{Float64}(undef,0);
        bu_rr = Array{Float64}(undef,0);
        bl_ll = Array{Float64}(undef,0);
        bl_rl = Array{Float64}(undef,0);
        bl_lr = Array{Float64}(undef,0);
        bl_rr = Array{Float64}(undef,0);

        Bu, Bl, peakMem

    elseif isa(treeL, Leaf) && !isa(treeR, Leaf) #If left node is a leaf, and right node is not
        #COMPUTE EXPANSION COEFFICIENTS

        #left block of upper right corner and upper block of lower left corner
        #(these are computed in the same pass)
        #    ---------------
        #   \       \   \   \
        #   \       \ x \<--\-- MxN
        #   \       \   \   \
        #   \---------------\
        #   \   x<--\-------\-- NxM
        #   \-------\       \
        #   \       \       \
        #    ---------------
        #           (5)

        Bu_l, Bl_l, peakMem = HSS_Expansion_Coeffs_OffDiag(treeL,treeR.treeL,row_idx,col_idx,N,peakMem);
        bu_l = Bu_l*treeR.Wl;
        bl_l = treeR.Rl'*Bl_l;
        #keep track of peak memory
        count = length(bu_l)+length(bl_l);
        peakMem[1] = peakMem[1]+count;
        peakMem[2] = max(peakMem[1],peakMem[2]);

        count = length(Bu_l)+length(Bl_l);
        peakMem[1] = peakMem[1]-count;
        peakMem[2] = max(peakMem[1],peakMem[2]);

        Bu_l = Array{Float64}(undef,0);
        Bl_l = Array{Float64}(undef,0);

        #right block of upper right corner, upper left block of lower left corner.
        #    ---------------
        #   \       \   \   \
        #   \       \   \ x \
        #   \       \   \   \
        #   \---------------\
        #   \       \       \
        #   \-------\       \
        #   \   x   \       \
        #    ---------------
        #           (6)

        col_idx = col_idx + treeR.treeL.m;
        Bu_r, Bl_r, peakMem = HSS_Expansion_Coeffs_OffDiag(treeL,treeR.treeR,row_idx,col_idx,N,peakMem);
        col_idx = col_idx - treeR.treeL.m;
        bu_r = Bu_r*treeR.Wr;
        bl_r = treeR.Rr'*Bl_r;
        #keep track of peak memory
        count = length(bu_r)+length(bl_r);
        peakMem[1] = peakMem[1]+count;
        peakMem[2] = max(peakMem[1],peakMem[2]);

        count = length(Bu_r)+length(Bl_r);
        peakMem[1] = peakMem[1]-count;
        peakMem[2] = max(peakMem[1],peakMem[2]);

        Bu_r = Array{Float64}(undef,0);
        Bl_r = Array{Float64}(undef,0);

        Bu = bu_l + bu_r;
        Bl = bl_l + bl_r;
        #keep track of peak memory
        count = length(Bu)+length(Bl);
        peakMem[1] = peakMem[1]+count;
        peakMem[2] = max(peakMem[1],peakMem[2]);

        count = length(bu_l)+length(bu_r)+length(bl_l)+length(bl_r);
        peakMem[1] = peakMem[1]-count;
        peakMem[2] = max(peakMem[1],peakMem[2]);

        bu_l = Array{Float64}(undef,0);
        bu_r = Array{Float64}(undef,0);
        bl_l = Array{Float64}(undef,0);
        bl_r = Array{Float64}(undef,0);

        Bu, Bl, peakMem

    elseif !isa(treeL, Leaf) && isa(treeR, Leaf) #If right node is a leaf, and left node is not
        #COMPUTE EXPANSION COEFFICIENTS

        #upper block of upper right corner and left block of lower left corner
        #(these are computed in the same pass)
        #    ---------------
        #   \       \   x <-\-- MxN
        #   \       \-------\
        #   \       \       \
        #   \---------------\
        #   \   \   \       \
        #   \ x \<--\-------\-- NxM
        #   \   \   \       \
        #    ---------------
        #           (7)

        Bu_l, Bl_l, peakMem = HSS_Expansion_Coeffs_OffDiag(treeL.treeL,treeR,row_idx,col_idx,N,peakMem);
        bu_l = treeL.Rl'*Bu_l;
        bl_l = Bl_l*treeL.Wl;
        #keep track of peak memory
        count = length(bu_l)+length(bl_l);
        peakMem[1] = peakMem[1]+count;
        peakMem[2] = max(peakMem[1],peakMem[2]);

        count = length(Bu_l)+length(Bl_l);
        peakMem[1] = peakMem[1]-count;
        peakMem[2] = max(peakMem[1],peakMem[2]);

        Bu_l = Array{Float64}(undef,0);
        Bl_l = Array{Float64}(undef,0);

        #upper block of upper right corner and left block of lower left corner
        #(these are computed in the same pass)
        #    ---------------
        #   \       \       \
        #   \       \-------\
        #   \       \   x   \
        #   \---------------\
        #   \   \   \       \
        #   \   \ x \       \
        #   \   \   \       \
        #    ---------------
        #           (8)

        row_idx = row_idx + treeL.treeL.m;
        Bu_r, Bl_r, peakMem = HSS_Expansion_Coeffs_OffDiag(treeL.treeR,treeR,row_idx,col_idx,N,peakMem);
        row_idx = row_idx - treeL.treeL.m;
        bu_r = treeL.Rr'*Bu_r;
        bl_r = Bl_r*treeL.Wr;
        #keep track of peak memory
        count = length(bu_r)+length(bl_r);
        peakMem[1] = peakMem[1]+count;
        peakMem[2] = max(peakMem[1],peakMem[2]);

        count = length(Bu_r)+length(Bl_r);
        peakMem[1] = peakMem[1]-count;
        peakMem[2] = max(peakMem[1],peakMem[2]);

        Bu_r = Array{Float64}(undef,0);
        Bl_r = Array{Float64}(undef,0);

        Bu = bu_l + bu_r;
        Bl = bl_l + bl_r;
        #keep track of peak memory
        count = length(Bu)+length(Bl);
        peakMem[1] = peakMem[1]+count;
        peakMem[2] = max(peakMem[1],peakMem[2]);

        count = length(bu_l)+length(bu_r)+length(bl_l)+length(bl_r);
        peakMem[1] = peakMem[1]-count;
        peakMem[2] = max(peakMem[1],peakMem[2]);

        bu_l = Array{Float64}(undef,0);
        bu_r = Array{Float64}(undef,0);
        bl_l = Array{Float64}(undef,0);
        bl_r = Array{Float64}(undef,0);

        Bu, Bl, peakMem
    end
end
