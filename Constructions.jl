include("Frames.jl")


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
