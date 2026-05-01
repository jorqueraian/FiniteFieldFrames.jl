export real_dx2d_etf_from_prime_power


function real_dx2d_etf_from_prime_power(prime_power::Int)::Matrix{BigFloat}
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