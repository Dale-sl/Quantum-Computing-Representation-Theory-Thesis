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
quantumQueryComplexity := function(group, character)
    local queries, tensorChar, constituents, groupTable, probabilityOfSuccess, boundedQueries;
    groupTable := CharacterTable(group);
    tensorChar := character;
    queries := 0;
    boundedQueries := 0;  
    repeat
        queries := queries + 1;
        constituents := ConstituentsOfCharacter(groupTable, tensorChar);
        probabilityOfSuccess := Sum(constituents, x -> DegreeOfCharacter(x)^2) / Size(group);
        if boundedQueries = 0 and probabilityOfSuccess >= 2/3 then
            boundedQueries := queries;
        fi;
        
        # Print(constituents);
        # Print("\n");
        # Print("\n");
        tensorChar := Tensored([character], [tensorChar])[1];
    until IsEqualSet(constituents, Irr(groupTable));
    return rec(group := group, queries := queries, boundedQueries := boundedQueries);
end;

# Now we want to define group actions and find their permutation character.

# natural action of given permutation group (constructing function) on elements of [1 .. n] 

regActChar := function(n, groupConstructor)
    return PermutationCharacter(groupConstructor(n), [1..n], OnPoints);
end;

# natural action of given permutation group on k element subsets of [1..n]

regActSubsetsChar := function(n,k, groupConstructor)
    return PermutationCharacter(groupConstructor(n), Combinations([1..n], k), OnSets);
end;


# We create functions to create regular partitions,
# which will be the domain of our action of the symmetric group on regular partitions.
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

# natural action of given permutation group (constructing function) on partitions of [1 .. ab] where there are b parts of size a
regActPartsChar := function(a,b, groupConstructor)
    return PermutationCharacter(groupConstructor(a * b), createRegularPartitions([1 .. a * b], b), OnSetsSets);
end;

# data collection for regular partitions
ProcessQuantumQueryComplexitiesHPC := function(groupConstructor, b, filename)
    local results, a, n, recResult, task, startTime, elapsed;
    
    results := [];
    a := 1;
    
    # Loop over increasing values of a (with total degree n = a * b)
    while true do
        n := a * b;
        
        # Launch a task to compute quantumQueryComplexity for the current group and action.
        task := RunTask( function()
            return quantumQueryComplexity(
                groupConstructor(n),                      # Construct the group of degree n.
                regActPartsChar(a, b, groupConstructor)     # Compute its permutation character.
            );
        end );
        
        # Poll until the task finishes, but if it takes longer than 60 seconds, cancel it.
        startTime := CPUTime();
        # Initially, recResult is undefined.
        recResult := false;
        while not TaskFinished(task) do
            elapsed := CPUTime() - startTime;
            if elapsed > 60 then
                CancelTask(task);
                recResult := fail;
                break;
            fi;
        od;
        
        if recResult = fail then
            Print("Timeout encountered for group with n = ", n, ". Stopping loop.\n");
            break;
        fi;
        
        recResult := TaskResult(task);
        if recResult = fail then
            Print("Task for group with n = ", n, " failed.\n");
            break;
        fi;
        
        Add(results, recResult);
        a := a + 1;
    od;
    
    # Write the results to a CSV file using PrintCSV.
    PrintCSV(filename, results, ["group", "queries", "boundedQueries"]);
    
    return results;
end;

# data collection for k-element subsets
ProcessQuantumQueryComplexitiesSubsetsHPC := function(groupConstructor, k, filename)
    local results, n, recResult, task, startTime, elapsed;
    
    results := [];
    n := k;  # Start at n = k since we need at least k points to form k-element subsets.
    
    # Loop over increasing values of n until a task exceeds 60 seconds.
    while true do
        # Launch a task to compute quantumQueryComplexity for groupConstructor(n)
        task := RunTask( function()
            return quantumQueryComplexity(
                groupConstructor(n),
                regActSubsetsChar(n, k, groupConstructor)
            );
        end );
        
        # Poll until the task is finished, but break out if it takes more than 60 seconds.
        startTime := CPUTime();
        while not TaskFinished(task) do
            elapsed := CPUTime() - startTime;
            if elapsed > 60 then
                CancelTask(task);
                recResult := fail;
                break;
            fi;
            # Optionally, one might insert a short delay here.
        od;
        
        if recResult = fail then
            Print("Timeout encountered for group with n = ", n, ". Stopping loop.\n");
            break;
        fi;
        
        recResult := TaskResult(task);
        if recResult = fail then
            Print("Task for group with n = ", n, " failed.\n");
            break;
        fi;
        
        Add(results, recResult);
        n := n + 1;
    od;
    
    # Write the results to a CSV file using PrintCSV.
    PrintCSV(filename, results, ["group", "queries", "boundedQueries"]);
    
    return results;
end;

ProcessQuantumQueryComplexitiesHPC(SymmetricGroup, 3, "SnPartitions3");
ProcessQuantumQueryComplexitiesHPC(AlternatingGroup, 3, "AnPartitions3");
ProcessQuantumQueryComplexitiesSubsetsHPC(SymmetricGroup, 3, "SnSubsets3");
ProcessQuantumQueryComplexitiesSubsetsHPC(AlternatingGroup, 3, "AnSubsets3");
