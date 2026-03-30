using DelimitedFiles

function matrix_from_file(file_loc::String, d::Int, n::Int)
    m = readdlm(file_loc);
    real_mat = reshape(m[1:d*n], (d,n));
    comp_mat = reshape(m[d*n+1:2*d*n], (d,n));
    if iszero(comp_mat)
        real_mat
    else
        real_mat + im*comp_mat
    end
end

function real_etf_to_case_O(gram::Matrix{Float64}, d::Int, char::Int)
    n = size(gram)[1];
    (n == size(gram)[2]) || throw(DomainError(size(gram),"gram must be square"));

    ff = GF(char);
    if n == 2*d
        if !is_square(ff(n-1))
            ff = GF(char, 2, "a");
        end
        a = sqrt(ff(n-1));
        c = 2*a;
    else
        a = ff(Int(sqrt(d*(n-1)/(n-d))));
        c = ff(Int((n/d)*sqrt(d*(n-1)/(n-d))));
    end
    
    gram = (1/abs(gram[1,2]))*gram;
    gram[diagind(gram)] .= 0;
    matrix(ff, round.(Int, gram)) + diagonal_matrix(a,n,n)
end

