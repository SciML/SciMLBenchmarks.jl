using LinearAlgebra

## compute the percentage banded for a matrix given a bandwidth
function compute_bandedness(A, bandwidth)
    n = size(A, 1)
    total_band_positions = 0
    non_zero_in_band = 0
    bandwidth = bandwidth
    for r in 1:n
        for c in 1:n
            if abs(r - c) <= bandwidth
                total_band_positions += 1  # This position belongs to the band
                if A[r, c] != 0
                    non_zero_in_band += 1  # This element is non-zero in the band
                end
            end
        end
    end

    percentage_filled = non_zero_in_band / total_band_positions * 100
    return percentage_filled
end

## use the compute banded function to compute the bandedness given different bandsizes
function isbanded(A)
    n = size(A, 1)
    count = 5
    potential_band_sizes = []

    while count < n - 5: 
        push!(potential_band_sizes, count += 5)
    end
    
    for bandsize in potential_band_sizes
        percentage_filled = compute_bandedness(A, bandsize)
        if percentage_filled >= 0.9
            return true 
        end
    end

    return false

## compute the sparsity for a given matrix
function compute_sparsity(A)
    n = size(A, 1)
    percentage_sparsity = length(nonzeros(A)) / n^2
    return percentage_sparsity
end

## get the percentage banded for a bandwidth of 5 and percentage sparsity
function getstructure(A::Matrix)::Any
    percentage_banded = compute_bandedness(A, 5)
    percentage_sparsity = compute_sparsity(A)

    return (percentage_banded, percentage_sparsity)
end

## return the best type of matrix for a given sparse matrix
function sparsestructure(A::SparseMatrixCSC)::Any
    sym = issymmetric(A)
    herm = ishermitian(A)
    banded = isbanded(A)
    posdef = isposdef(A)
    lower_triangular = istril(A)
    upper_triangular = istriu(A)


    n = size(A, 1)

    if banded
        return BandedMatrix(A)
    end
    
    if sym 
       return Symmetric(A)
    end 

    if herm
        return Hermitian(A)
    end

    if lower_triangular
        return LowerTriangular(A)
    end

    if upper_triangular
        return UpperTriangular(A)
    end

    return SparseMatrixCSC(A)
end





    