# Loading packages
LoadPackage( "repsn" );
LoadPackage("RepnDecomp");

# regular partition of set of size ab a partition of b parts of cardinality a.
testa := [1 .. 4];
testb := PartitionsSet(testa, 2);
# takes a partition and the intended size of each partition for the permutation to be regular.
# returns a function which takes a partition and returns whether or not the partition is an a-regular partition.
isRegularPartition := function(a)
    return function(n)
        local size, x;
        size:= Length(n[1]);
        for x in n do
            if not (Length(x) = a) then
                return false;
            fi;
        od;
        return true;
    end;
end;
# takes a set of partitions and the size each partition should be for the partitions to be regular.
# returns a set of regular partitions with each partition of size a.
regularPartitions := function(n,a)
    return Filtered(n, isRegularPartition(a));
end;
# takes a range and a number of partitions b as input.
# b should divide the cardinality of the range.
# returns the regular partitions of the range of b parts.
createRegularPartitions := function(range, b)
    local a;
    a := Length(range) / b;
    return regularPartitions(PartitionsSet(range, b), a);
end;

testc:= regularPartitions(testb, 2);
testd:= createRegularPartitions([1 .. 6], 2);

# natural action of symmetric group on k element sets of [1 .. n] section
#Combinations(range, k);

# natural action of symmetric group on elements of [1 .. n] section

symmetricGroupTable := function(n)
    return CharacterTable(Action( SymmetricGroup(n), [1 .. n], OnPoints));
end;

teste := symmetricGroupTable(6);
symmetricRegularAction := Action(SymmetricGroup(3), [1 .. 3], OnPoints);  # Permutation action of G on X
# createPermutationMatrix := function(g, X, action)
#     local n, P, i, j;
#     n := Length(X);                        # Dimension of the space
#     P := ZeroMatrix(Integers, n, n);      # Initialize an n x n zero matrix
#     for i in [1..n] do
#         j := Position(X, action( X[i], g)); # Find where g maps X[i]
#         P[i][j] := 1;                      # Set the (i,j)-th entry to 1
#     od;
#     return P;
# end;

# representation := function(g, X, action)
#     return createPermutationMatrix(g, X, action); # Map g to its permutation matrix
# end;

# Example: Get the matrices for all elements of G
# testMatrices := List(GeneratorsOfGroup(SymmetricGroup(3)), g -> representation(g, [1 .. 3], OnPoints));
G := SymmetricGroup(4);
testImages := List(GeneratorsOfGroup(G), g -> PermutationMat(g, 4));
testRho := GroupHomomorphismByImages(G, Group(testImages));
testTensor := TensorProductOfRepresentations(testRho, testRho);
#testIrr := REPN_ComputeUsingSerre(testTensor );
