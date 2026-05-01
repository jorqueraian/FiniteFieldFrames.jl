export etf_from_triangular_graph, real_etf_to_case_O, etf_from_pmod_diff_set, etf_from_modular_hadamard

function etf_from_triangular_graph(d::Int, p::Int)::FqMatrix
    # This construction comes from Theorem 5.4 in [2].
    # Constructs a d by ((d^2-d)/2) ETF in Case O, as long as p divides d-7.
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


function etf_from_pmod_diff_set(D, n, q)
    # This construction comes from Theorem 5.7 in [1].
    # Constructs a d=|D| by n ETF in Case U, as long as n divides q+1 
    # and D is a p-modular difference set. 
    # This construction does not verify that D is p-modular difference set.
    (e, b) = is_power(q);
    (is_prime(b)) || throw(DomainError(q, "input must be a prime power: q=p^k"));
    (gcd(n, q+1)==n) || throw(DomainError((n,q), "n must divide q+1"));
    
    base_ff = GF(q);
    
    ## In general, gen will not return a multiplicative generator.
    ## But becasue base_ff is defined over the prime subfield, it will!
    base_ff_gen = gen(base_ff)
    Kx, x = base_ff["x"];
    ff = GF(x^2-base_ff_gen, "a");
    
    ## becasue ff is defined as a field extention, gen would give a 
    ## generator with respect to the base_field.
    ff_gen = nothing
    for r in ff
        ff_gen = r
        findall(is_one, [r^i for i in 1:(q^2-1)]).size[1] == 1 && break
    end

    findall(is_one, [ff_gen^i for i in 1:(q^2-1)]).size[1] == 1 || throw(error("Couldnt find multiplicative generator of finite field."))

    w = ff_gen^Int((q^2-1)/n)
    F = matrix(ff, [[w^(i*j) for i in 0:(n-1)] for j in 0:(n-1)])
    F[D.+1,:]
end


function real_etf_to_case_O(gram::Matrix{BigFloat}, d::Int, char::Int)::FqMatrix
    # This construction comes from Proposition 3.2 in [2].
    # Constructs an d'xn ETF in case O from a real dxn ETF. 
    # d'<= d with equality depending on char.
    
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


function etf_from_modular_hadamard(modHadamard::FqMatrix, return_gram::Bool=true)::FqMatrix
    # This construction comes from Theorem 2 in [4].
    # Constructs a dxd^2 ETF in case U, given a modular hadamard defined over 
    # a quadratic extension of a finite field such that x^2+1 does not split in 
    # the base field and d = 8 mod characteristic.

    d = size(modHadamard)[1];
    (d == size(modHadamard)[2]) || throw(DomainError(size(modHadamard),"Modular Hadamard must be square"));
    
    ff = modHadamard.base_ring;
    p = characteristic(ff);
    
    base_ff = base_field(ff);

    ## need root of x^2+1
    ## this means that in the base field -1 is not a square.
    (degree(ff) == 2) || throw(DomainError(ff,"in Case U, the provided field must be finite and must be a degree 2 extension"));
    !is_square(base_ff(-1)) || throw(DomainError(ff,"The base field must have -1 be a non-square"));

    (gcd(d-8, p) == p) || throw(DomainError(characteristic(modHadamard.base_ring),"must be that d = 8 mod p"));
    iszero(transpose(modHadamard)*modHadamard - diagonal_matrix([ff(d) for i in 1:d])) || throw(DomainError(modHadamard,"input must be modular Hadamard, does not satsify H^t*H=dI"));
    iszero(modHadamard.^2 .- ff(1)) || throw(DomainError(modHadamard,"input must be modular Hadamard, entries must square to 1."));


    ffx, x = ff["x"];
    ff_im = roots(x^2+1)[1]
        

    Phi = matrix(ff, zeros(Int, (d,d^2)));
    for i in 1:d 
        thing = matrix(ff, ones(Int, (d,1)));
        thing[i] += (ff(-2)*(ff(1)+ff_im));
        for j in 1:d 
            ind = j+(i-1)*d;
            Phi[:,ind] = matrix(modHadamard[:,j]).*thing;
        end
    end
    if return_gram
        conjugate_transpose(Phi)*Phi
    else
        Phi
    end
end