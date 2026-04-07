export maximal_case_O_etf, real_etf_to_case_O, real_dx2d_etf

function maximal_case_O_etf(d::Int, p::Int)
    (d > 1 && gcd(p, d-7) == p) || throw(DomainError((d, p),"p must divide d-7"));

    n = Int(d*(d+1)/2);
    two_sets = collect(combinations(d+1,2));
    gram = zeros(Int, (n, n))

    for verts in combinations(n,2) 
        i = verts[1];
        j = verts[2];
        i_set = two_sets[i];
        j_set = two_sets[j];

        if i_set[1] == j_set[1] || i_set[1] == j_set[2] || i_set[2] == j_set[1] || i_set[2] == j_set[2]
            gram[i,j] = 1;  ## This might be flipped
        else
            gram[i,j] = -1;
        end
    end
    gram += gram'
    matrix(GF(p), gram + 3*I[1:n,1:n])
end


function real_etf_to_case_O(gram::Matrix{BigFloat}, d::Int, char::Int)::FqMatrix
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


function real_dx2d_etf(prime_power::Int)
    (e, b) = is_power(prime_power);
    (is_prime(b) && (gcd(e,2)==2 || gcd(b-3, 4)==1 )) || throw(DomainError(prime_power, "input must be a prime power: p^k, such that k is even or p is not 3 mod 4"))
    ff = GF(b^e);

    a = sqrt(BigFloat(prime_power))

    conf_mat = ones(BigFloat, (b^e+1, b^e+1));
    for (i,x) in Iterators.enumerate(ff)
        for (j,y) in Iterators.enumerate(ff)
            if i == j 
                conf_mat[i,j]=a;
            elseif !is_square(y-x)
                conf_mat[i,j]=-1;
            end
        end
    end
    conf_mat[b^e+1, b^e+1] = a;

    conf_mat
end