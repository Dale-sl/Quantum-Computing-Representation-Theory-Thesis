# Loading packages
# as it turns out we don't actually need these to do character theory, which is all we need, so I've commented these out
# LoadPackage( "repsn" );
# LoadPackage("RepnDecomp");

# Our goal is to calculate the quantum query complexity (of finding the base size) of a given group action on a domain,
# which is equivalent to the number of times you need to tensor the permutation representation of the action
# with itself for its decomposition into irreducibles to contain every irreducible character of the group.

# Note that GAP has an interpreter, so I've been using these functions in the command line to analyze the output.
# to do so, write "Read("testing.g");" when GAP is opened in the directory containing this file,
# then call functions in this file.

# for a given group and a character of a permutation representation of an action of that group,
# we take the tensor product of the permutation character with itself until the
# decomposition contains every irreducible of the group.
# quantumQueryComplexity := function(group, character, groupDegree)
#     local queries, tensorChar, constituents, groupTable, probabilityOfSuccess, boundedQueries;
#     groupTable := CharacterTable(group);
#     tensorChar := character;
#     queries := 0;
#     boundedQueries := 0;  
#     repeat
#         queries := queries + 1;
#         constituents := ConstituentsOfCharacter(groupTable, tensorChar);
#         probabilityOfSuccess := Sum(constituents, x -> DegreeOfCharacter(x)^2) / Size(group);
#         if boundedQueries = 0 and probabilityOfSuccess >= 2/3 then
#             boundedQueries := queries;
#         fi;
        
#         # Print(constituents);
#         # Print("\n");
#         # Print("\n");
#         # tensorChar := Tensored([character], [tensorChar])[1];
#         tensorChar := character * tensorChar;
#     until IsEqualSet(constituents, Irr(groupTable));
#     return rec(groupDegree := groupDegree, queries := queries, boundedQueries := boundedQueries);
# end;

# new version that associates a partition with each irreducible character
quantumQueryComplexity := function(tbl, permutationCharacter, groupDegree)
    local queries, tensorChar, constituents, probabilityOfSuccess, boundedQueries,
    irr, cp, mapping, partitionsHistory, groupOrder, i, pos, currentIterationPartitions;
    
    queries := 0;
    boundedQueries := 0;  
    tensorChar := permutationCharacter;
    
    # Get the list of irreducible characters from the table.
    irr := Irr(tbl);
    if Size(ConstituentsOfCharacter(tbl, permutationCharacter)) = 1 then
        Print("Initial character only contains the trivial character, would loop forever.");
        return;
    fi;
    
    # Compute group order as the sum of squares of the degrees.
    groupOrder := Sum(irr, x -> x[1]^2);
    
    # Get the character parameters; each entry isq of the form [ 1, partition ].
    cp := CharacterParameters(tbl);
    # Build a mapping: for each irreducible (by position), store the associated partition.
    mapping := [];
    for i in [1..Length(cp)] do
        mapping[i] := cp[i][2];
    od;
    
    # Initialize the history of partitions encountered.
    partitionsHistory := [];
    
    repeat
        queries := queries + 1;
        constituents := ConstituentsOfCharacter(tbl, tensorChar);
        probabilityOfSuccess := Sum(constituents, x -> x[1]^2) / groupOrder;
        if boundedQueries = 0 and probabilityOfSuccess >= 2/3 then
            boundedQueries := queries;
        fi;
        
        # For each constituent, find its associated partition (by position in irr).
        currentIterationPartitions := [];
        for constituent in constituents do
            pos := Position(irr, constituent);
            if pos = fail then
                Add(currentIterationPartitions, []);
                Print("Something is very wrong, couldn't find partition associated to irreducible!");
            else
                Add(currentIterationPartitions, mapping[pos]);
            fi;
        od;
        Add(partitionsHistory, currentIterationPartitions);
        
        # Tensor the permutationCharacter with itself.
        tensorChar := permutationCharacter * tensorChar;
    until IsEqualSet(constituents, irr);
    
    return rec(
      groupDegree      := groupDegree,
      queries          := queries,
      boundedQueries   := boundedQueries,
      partitionHistory := partitionsHistory
    );
end;

quantumQueryProbabilityHistory := function(tbl, permutationCharacter, groupDegree)
    local tensorChar, constituents, probabilityOfSuccess, groupOrder, irr, probHistory;
    
    tensorChar := permutationCharacter;
    irr := Irr(tbl);
    if Size(ConstituentsOfCharacter(tbl, permutationCharacter)) = 1 then
        Print("Initial character only contains the trivial character.");
        return;
    fi;
    
    # Compute the group order as the sum of the squares of the degrees.
    groupOrder := Sum(irr, x -> x[1]^2);
    
    # Initialize the list that will store the probability of success at each iteration.
    probHistory := [];
    
    repeat
        constituents := ConstituentsOfCharacter(tbl, tensorChar);
        probabilityOfSuccess := Sum(constituents, x -> x[1]^2) / groupOrder;
        Add(probHistory, probabilityOfSuccess);
        
        # Tensor the permutationCharacter with the current tensorChar.
        tensorChar := permutationCharacter * tensorChar;
    until IsEqualSet(constituents, irr);
    
    return rec(
      groupDegree   := groupDegree,
      probabilities := probHistory
    );
end;

quantumQueryCharacterPartitions := function(tbl, permutationCharacter, groupDegree)
    local constituents, irr, cp, mapping, partitionList, i, pos;
    
    # Get the list of irreducible characters from the table.
    irr := Irr(tbl);
    
    # Get the character parameters; each entry is of the form [ 1, partition ].
    cp := CharacterParameters(tbl);
    
    # Build a mapping: for each irreducible (by position), store the associated partition.
    mapping := [];
    for i in [1..Length(cp)] do
        mapping[i] := cp[i][2];
    od;
    
    # Compute the constituents of the given permutation character.
    constituents := ConstituentsOfCharacter(tbl, permutationCharacter);
    
    # For each constituent, find its associated partition.
    partitionList := [];
    for constituent in constituents do
        pos := Position(irr, constituent);
        if pos = fail then
            Add(partitionList, []);
            Print("Something is very wrong, couldn't find partition associated to irreducible!");
        else
            Add(partitionList, mapping[pos]);
        fi;
    od;
    
    return rec(
         groupDegree := groupDegree,
         partitions  := partitionList
    );
end;

TensorPartitions := function(partition1, partition2, n)
    local tbl, cp, irr, mapping, pos1, pos2, chi1, chi2, tensorChar, constituents, result, constituent, pos;
    
    # Create the symmetric group character table for degree n.
    tbl := CharacterTable("symmetric", n);
    
    # Obtain the character parameters; each entry is of the form [ 1, partition ].
    cp := CharacterParameters(tbl);
    
    # Get the list of irreducible characters.
    irr := Irr(tbl);
    
    # Build a mapping so that mapping[i] is the partition corresponding to irr[i].
    mapping := [];
    for i in [1..Length(cp)] do
        mapping[i] := cp[i][2];
    od;
    
    # Find the positions corresponding to partition1 and partition2.
    pos1 := Position(mapping, partition1);
    pos2 := Position(mapping, partition2);
    if pos1 = fail or pos2 = fail then
        Error("One of the given partitions is not valid for S_", n, ".");
    fi;
    
    chi1 := irr[pos1];
    chi2 := irr[pos2];
    
    # Compute the tensor product of the two characters.
    tensorChar := chi1 * chi2;
    
    # Decompose the tensor product into irreducible constituents.
    constituents := ConstituentsOfCharacter(tbl, tensorChar);
    
    # For each constituent, look up its corresponding partition.
    result := [];
    for constituent in constituents do
        pos := Position(irr, constituent);
        if pos = fail then
            Add(result, "Error");
        else
            Add(result, mapping[pos]);
        fi;
    od;
    
    return result;
end;

# manually getting a csv of the partitions of a specific action
# PrintCSV("SnParts16.csv", quantumQueryCharacterPartitions(CharacterTable("symmetric", 16), regActPartsChar(4,4, SymmetricGroup), 16), ["groupDegree", "partitions"]);


# Now we want to define group actions and find their permutation character.

# natural action of given permutation group (constructing function) on elements of [1 .. n] 
regActChar := function(n, groupConstructor)
    return PermutationCharacter(groupConstructor(n), [1..n], OnPoints);
end;

# natural action of given permutation group on k element subsets of [1..n]
regActSubsetsChar := function(n,k, groupConstructor)
    return PermutationCharacter(groupConstructor(n), Combinations([1..n], k), OnSets);
end;

# natural action of given permutation group (constructing function) on partitions of [1 .. ab] where there are b parts of size a
regActPartsChar := function(a,b, groupConstructor)
    return PermutationCharacter(groupConstructor(a * b), RegularPartitionsSet(a * b, b), OnSetsSets);
end;

# The following section contains auxillary functions so we can create the desired actions.

# Main function.
# RegularPartitionsSet( n, b )
#   n: size of the set [1..n]
#   b: number of parts; we require that b divides n.
#
# Returns the list of all (a,b)-regular partitions of [1..n], where a = n/b.
RegularPartitionsSet := function(n, b)
    local a, L;
    if n mod b <> 0 then
        Error("RegularPartitionsOfSet: n must be divisible by b.");
    fi;
    a := n / b;
    L := [1..n];
    return PartitionsOfList(L, b, a, 1);
end;

# Auxiliary recursive function.
# PartitionsOfList( L, b, a, minVal )
#   L: a sorted list (the set to be partitioned)
#   b: number of parts to partition L into
#   a: target size of each part (so that Length(L) should equal a*b)
#   minVal: lower bound for the first element of the part to be chosen.
#
# It returns a list of partitions; each partition is a list of b subsets (each a list of a numbers),
# with the property that the first element of each subset is at least minVal,
# and (for uniqueness) if the partition is read in order the first element of part i is strictly less than that of part i+1.
PartitionsOfList := function(L, b, a, minVal)
    local partitions, comb, L2, subparts, p;
    partitions := [];
    if b = 0 then
        # If no parts are to be chosen, L must be empty.
        if Length(L) = 0 then
            return [ [] ];
        else
            return [];
        fi;
    fi;
    if b = 1 then
        # Only one part remains: it must exactly equal L and have length a,
        # and we require its first element to be at least minVal.
        if Length(L) = a and L[1] >= minVal then
            return [ [ L ] ];
        else
            return [];
        fi;
    fi;
    # For b > 1, choose a combination 'comb' of a elements from L with comb[1] >= minVal.
    for comb in Combinations(L, a) do
        if comb[1] >= minVal then
            # Remove the chosen elements from L to get L2.
            L2 := Difference(L, comb);
            # For the remaining parts, require that the first element is greater than comb[1].
            subparts := PartitionsOfList(L2, b - 1, a, comb[1] + 1);
            for p in subparts do
                Add(partitions, Concatenation([ comb ], p));
            od;
        fi;
    od;
    return partitions;
end;

# The following section contains data collection functions.

# data collection for regular partitions
ProcessQuantumQueryComplexities := function(groupConstructor, groupName, b, filename)
    local results, a, n, ret;
    
    results := [];
    a := 2; # a starts at 2 to insure a nontrivial action
    
    while true do
        n := a * b;
        Print("Starting computation for n = ", n, "\n");
        
        ret := IO_CallWithTimeout( rec(seconds := 30), 
            function()
                return quantumQueryComplexity(
                           CharacterTable(groupName, n),
                           regActPartsChar(a, b, groupConstructor),
                           n
                       );
            end );
        
        if ret[1] = false then
            Print("Timeout encountered for group with n = ", n, ". Stopping loop.\n");
            break;
        elif ret[1] = fail then
            Print("Computation failed for group with n = ", n, ". Stopping loop.\n");
            break;
        elif ret[1] = true then
            Add(results, ret[2]);
        fi;
        
        a := a + 1;
    od;
    
    PrintCSV(filename, results, ["groupDegree", "queries", "boundedQueries", "partitionHistory"]);
    return results;
end;

# data collection for k element subsets
ProcessQuantumQueryComplexitiesSubsets := function(groupConstructor, groupName, k, filename)
    local results, n, ret;
    
    results := [];
    n := k + 1;  # We need at least k points to form k-element subsets, k+1 insures nontrivial action
    
    while true do
        Print("Starting computation for n = ", n, "\n");
        
        ret := IO_CallWithTimeout( rec(seconds := 30), 
            function()
                return quantumQueryComplexity(
                           CharacterTable(groupName, n),
                           regActSubsetsChar(n, k, groupConstructor),
                           n
                       );
            end );
        
        if ret[1] = false then
            Print("Timeout encountered for group with n = ", n, ". Stopping loop.\n");
            break;
        elif ret[1] = fail then
            Print("Computation failed for group with n = ", n, ". Stopping loop.\n");
            break;
        elif ret[1] = true then
            Add(results, ret[2]);
        fi;
        
        n := n + 1;
    od;
    
    PrintCSV(filename, results, ["groupDegree", "queries", "boundedQueries", "partitionHistory"]);
    return results;
end;

ProcessQuantumQueryProbabilityHistoryRegular := function(groupConstructor, groupName, b, filename)
    local results, a, n, ret;
    
    results := [];
    a := 2;  # a starts at 2 to ensure a nontrivial action
    
    while true do
        n := a * b;
        Print("Starting computation for n = ", n, "\n");
        
        ret := IO_CallWithTimeout( rec(seconds := 30), 
            function()
                return quantumQueryProbabilityHistory(
                           CharacterTable(groupName, n),
                           regActPartsChar(a, b, groupConstructor),
                           n
                       );
            end );
        
        if ret[1] = false then
            Print("Timeout encountered for group with n = ", n, ". Stopping loop.\n");
            break;
        elif ret[1] = fail then
            Print("Computation failed for group with n = ", n, ". Stopping loop.\n");
            break;
        elif ret[1] = true then
            Add(results, ret[2]);
        fi;
        
        a := a + 1;
    od;
    
    PrintCSV(filename, results, ["groupDegree", "probabilities"]);
    return results;
end;

ProcessQuantumQueryProbabilityHistorySubsets := function(groupConstructor, groupName, k, filename)
    local results, n, ret;
    
    results := [];
    n := k + 1;  # We need at least k points to form k-element subsets; k+1 ensures nontrivial action
    
    while true do
        Print("Starting computation for n = ", n, "\n");
        
        ret := IO_CallWithTimeout( rec(seconds := 30), 
            function()
                return quantumQueryProbabilityHistory(
                           CharacterTable(groupName, n),
                           regActSubsetsChar(n, k, groupConstructor),
                           n
                       );
            end );
        
        if ret[1] = false then
            Print("Timeout encountered for group with n = ", n, ". Stopping loop.\n");
            break;
        elif ret[1] = fail then
            Print("Computation failed for group with n = ", n, ". Stopping loop.\n");
            break;
        elif ret[1] = true then
            Add(results, ret[2]);
        fi;
        
        n := n + 1;
    od;
    
    PrintCSV(filename, results, ["groupDegree", "probabilities"]);
    return results;
end;


ProcessCharacterPartitionsFactorizations := function(groupConstructor, groupName, filename)
    local results, n, factors, a, b, ret, factorPairs, factorFound, recEntry, usedFactorPairs;
    
    results := [];
    n := 2;
    
    while true do
        Print("Starting computations for n = ", n, "\n");
        
        # Get all factors of n
        factors := FactorsInt(n);
        
        # Initialize list to store used ordered pairs to avoid duplicate calls.
        usedFactorPairs := [];
        factorPairs := [];
        for a in factors do
            b := n / a;
            if a > 1 and b > 1 then
                if not ([a, b] in usedFactorPairs) then
                    Add(factorPairs, [a, b]);
                    Add(usedFactorPairs, [a, b]);
                    # If the factors are distinct, mark both orders as used.
                    if a <> b then
                        Add(usedFactorPairs, [b, a]);
                    fi;
                fi;
            fi;
        od;
        
        # If no nontrivial factorization exists for n, skip to the next n.
        if Length(factorPairs) = 0 then
            Print("No nontrivial factorizations for n = ", n, ". Skipping to next n.\n");
            n := n + 1;
            continue;
        fi;
        
        factorFound := false;
        # Iterate over each factor pair.
        for factors in factorPairs do
            a := factors[1];
            b := factors[2];
            
            if a = b then
                # For a square factorization, call only once.
                Print("Trying factorization: ", a, " * ", b, " = ", n, "\n");
                ret := IO_CallWithTimeout(rec(seconds := 600),
                    function()
                        return quantumQueryCharacterPartitions(
                                   CharacterTable(groupName, n),
                                   regActPartsChar(a, b, groupConstructor),
                                   n
                               );
                    end);
                if ret[1] = true then
                    factorFound := true;
                    recEntry := rec(
                        a           := a,
                        b           := b,
                        groupDegree := n,
                        partitions  := ret[2].partitions
                    );
                    Add(results, recEntry);
                fi;
            else
                # For distinct factors, make two calls: (a, b) and (b, a).
                Print("Trying factorization: ", a, " * ", b, " = ", n, "\n");
                ret := IO_CallWithTimeout(rec(seconds := 600),
                    function()
                        return quantumQueryCharacterPartitions(
                                   CharacterTable(groupName, n),
                                   regActPartsChar(a, b, groupConstructor),
                                   n
                               );
                    end);
                if ret[1] = true then
                    factorFound := true;
                    recEntry := rec(
                        a           := a,
                        b           := b,
                        groupDegree := n,
                        partitions  := ret[2].partitions
                    );
                    Add(results, recEntry);
                fi;
                
                Print("Trying factorization: ", b, " * ", a, " = ", n, "\n");
                ret := IO_CallWithTimeout(rec(seconds := 600),
                    function()
                        return quantumQueryCharacterPartitions(
                                   CharacterTable(groupName, n),
                                   regActPartsChar(b, a, groupConstructor),
                                   n
                               );
                    end);
                if ret[1] = true then
                    factorFound := true;
                    recEntry := rec(
                        a           := b,
                        b           := a,
                        groupDegree := n,
                        partitions  := ret[2].partitions
                    );
                    Add(results, recEntry);
                fi;
            fi;
        od;
        
        # If no factorization for this n produced a successful call, stop looping.
        if not factorFound then
            Print("All calls failed for n = ", n, ". Stopping loop.\n");
            break;
        fi;
        
        n := n + 1;
    od;
    
    PrintCSV(filename, results, ["a", "b", "groupDegree", "partitions"]);
    return results;
end;



# ProcessQuantumQueryComplexities(SymmetricGroup, "symmetric", 3, "SnPartitions3.csv");
# ProcessQuantumQueryComplexities(AlternatingGroup, "alternating", 3, "AnPartitions3.csv");
# ProcessQuantumQueryComplexitiesSubsets(SymmetricGroup, "symmetric", 3, "SnSubsets3.csv");
# ProcessQuantumQueryComplexitiesSubsets(AlternatingGroup, "alternating", 3, "AnSubsets3.csv");


























# code graveyard

# # We create functions to create regular partitions,
# # which will be the domain of our action of the symmetric group on regular partitions.
# # a regular partition of set of size ab is a partition of b parts of cardinality a.

# # takes a partition and the intended size of each partition for the permutation to be regular.
# # returns a function which takes a partition and returns whether or not the partition is an a-regular partition.
# isRegularPartition := function(a)
#     return function(n)
#         local size, x;
#         size:= Length(n[1]);
#         for x in n do
#             if not (Length(x) = a) then
#                 return false;
#             fi;
#         od;
#         return true;
#     end;
# end;

# # takes a set of partitions and the size each partition should be for the partitions to be regular.
# # returns a set of regular partitions with each partition of size a.
# regularPartitions := function(n,a)
#     return Filtered(n, isRegularPartition(a));
# end;

# takes a range and a number of partitions b as input.
# b should divide the cardinality of the range.
# returns the regular partitions of the range of b parts.
# createRegularPartitions := function(range, b)
#     local a;
#     a := Length(range) / b;
#     return regularPartitions(PartitionsSet(range, b), a);
# end;
