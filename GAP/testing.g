# Loading packages
# as it turns out we don't actually need these to do character theory, which is all we need, so I've commented these out
# LoadPackage( "repsn" );
# LoadPackage("RepnDecomp");

# Our goal is to calculate the quantum query complexity (of finding the base size) of a given group action on a domain,
# which is equivalent to the number of times you need to tensor the permutation representation of the action
# with itself for its decomposition into irreducibles to contain every irreducible character of the group.

# Note that GAP has an interpreter, so I've been using these functions in the command line to analyze the output.
# to do so, write "Read("testing.g");" when GAP is opened in the directory containing this file, then call functions in this file.

# for a given group and a character of a permutation representation of an action of that group,
# we take the tensor product of the permutation character with itself until the decomposition contains every irreducible of the group.
quantumQueryComplexity := function(group, character)
    local count, tensorChar, constituents, groupTable;
    groupTable := CharacterTable(group);
    tensorChar := character;
    count := 0;
    repeat
        count := count + 1;
        constituents := ConstituentsOfCharacter(groupTable, tensorChar);                                                            
        Print(constituents);
        Print("\n");
        tensorChar := Tensored([character], [tensorChar])[1];
    until IsEqualSet(constituents, Irr(groupTable));
    return count;
end;

# Now we want to define group actions and find their permutation character.

# natural action of symmetric group on elements of [1 .. n] 

regSymActChar := function(n)
    return PermutationCharacter(SymmetricGroup(n), [1..n], OnPoints);
end;

# We create functions to create regular partitions, which will be the domain of our action of the symmetric group on regular partitions.
# a regular partition of set of size ab is a partition of b parts of cardinality a.

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

# natural action of symmetric group on partitions of [1 .. ab] where there are b parts of size a
regSymActPartsChar := function(a,b)
    return PermutationCharacter(SymmetricGroup(a * b), createRegularPartitions([1 .. a * b], b), OnSetsSets);
end;



# a bunch of wasted code where I was working with actual representations instead of their characters, which is inefficent and uneccesary

#symmetricRegularAction := Action(SymmetricGroup(4), [1 .. 4], OnPoints);  # Permutation action of G on X

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
# G := SymmetricGroup(4);

# testImages := List(GeneratorsOfGroup(G), g -> PermutationMat(g, 4));
# testRho := GroupHomomorphismByImages(G, Group(testImages));
# testTensor := TensorProductOfRepresentations(testRho, testRho);
# testIrr := REPN_ComputeUsingSerre(testTensor );
