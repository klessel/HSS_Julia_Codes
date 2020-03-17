cd("C:\\Users\\klessel\\Dropbox\\Julia_home")
push!(LOAD_PATH, pwd())
# push!(LOAD_PATH,"/home/klessel/Dropbox/Julia_home/Polynomials/src")
using HSS_Types
using Polynomials
using LinearAlgebra
# include("/home/klessel/Dropbox/Julia_home/Gen_HSS_Memory_Efficient.jl")
# include("/home/klessel/Dropbox/Julia_home/Build_Matrix_From_HSS.jl")
include("C:\\Users\\klessel\\Dropbox\\Julia_home\\Gen_HSS_Memory_Efficient.jl")
include("C:\\Users\\klessel\\Dropbox\\Julia_home\\Build_Matrix_From_HSS.jl")
include("C:\\Users\\klessel\\Dropbox\\Julia_home\\FMM_Functions.jl")
include("C:\\Users\\klessel\\Dropbox\\Julia_home\\svd_ranktol.jl")
N = Int64;
r = Int64;
d = Int64;
minPart = Int64;
N = 256;
r=3; #hankel block rank
minPart = 3*r;

# d = 1;
# tree, dummy = Gen_Tree_Complete(minPart,N,d);

#tree = Gen_Tree_Mtree(minPart,N);
tree, dummy = label_tree(tree,1);

##
# #complete test tree
# N = 160;
# leaf1 = Leaf_Spine(4,4,-1,-1);
# leaf2 = Leaf_Spine(8,8,-1,-1);
# node1 = Node_Spine(leaf1,leaf2,12,12,-1,-1);
#
# leaf3 = Leaf_Spine(7,7,-1,-1);
# leaf4 = Leaf_Spine(21,21,-1,-1);
# node2 = Node_Spine(leaf3,leaf4,28,28,-1,-1);
#
# leaf5 = Leaf_Spine(14,14,-1,-1);
# leaf6 = Leaf_Spine(16,16,-1,-1);
# node3 = Node_Spine(leaf5,leaf6,30,30,-1,-1);
#
# leaf7 = Leaf_Spine(25,25,-1,-1);
# leaf8 = Leaf_Spine(65,65,-1,-1);
# node4 = Node_Spine(leaf7,leaf8,90,90,-1,-1);
#
# node5 = Node_Spine(node1,node2,40,40,-1,-1);
# node6 = Node_Spine(node3,node4,120,120,-1,-1);
#
# tree = Node_Spine(node5,node6,160,160,-1,-1);
# tree, dummy = label_tree(tree,1);
##
println("N = ", N, ", rank = ", r, ", Minimum Partition = ", minPart);
hss, peakMem = @time Gen_HSS_Memory_Efficient(tree,r);
println("Peak memory count is ", peakMem[2], " Float64 assignments");

#test
A = Array{Float64,2}(undef,N,N);
A = f(1:N,1:N,N);

B, dummy1, dummy2 = Build_Matrix_From_HSS(hss);
rel_err =  norm(A-B)/norm(A);
println("||A-B||/||A|| = ", rel_err);
