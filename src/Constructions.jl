export etf_from_triangular_graph, real_etf_to_case_O, etf_from_pmod_diff_set, etf_from_modular_hadamard, etf_from_singer_diff_set

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


function is_conway(ff)
    p = characteristic(ff);
    deg = degree(ff)
    (deg == absolute_degree(ff)) || throw(error("bad"));

    coeffs = [Int(0) for i  in 0:deg];
    # The type that _nmod_poly_conway wants for the first param is a nn_ptr.
    # I dont know what this is though. So i am not sure if this function is writing outside of coeffs.
    res = @ccall Oscar.Nemo.libflint._nmod_poly_conway(coeffs::Ref{Int}, p::UInt, deg::Int)::Int

    if res == 0
        false
    else
        poly_to_check = defining_polynomial(ff);
        conway_poly = poly_to_check.parent(coeffs);
        
        conway_poly == poly_to_check
    end
end


function etf_from_pmod_diff_set(D, n::T, q::T; return_gram::Bool=true, verify_mult_gen::Bool=false) where {T<:Union{Int,BigInt}}
    # This construction comes from Theorem 5.7 in [1].
    # Constructs a d=|D| by n ETF in Case U, as long as n divides q+1 
    # and D is a p-modular difference set. 
    # This construction does not verify that D is p-modular difference set.
    (e, b) = is_power(q);
    (is_prime(b)) || throw(DomainError(q, "input must be a prime power: q=p^k"));
    (gcd(n, q+1)==n) || throw(DomainError((n,q), "n must divide q+1"));
    
    ff = GF(BigInt(q)^2);
    
    # According to Oscar.jl this is not guaranteed to be a multiplicative generator.
    # However since Oscar.jl is using Nemo.jl which is using FLINT
    # calling GF will attempt to instantiate a finite field using a Conways Polynomial, if possible.
    # It does not appear that Nemo.jl has provide wrapper functionality for checking if this is the case.
    # So all that can be done is to hope that the underlying field is using a Conway Polynomial, and gen gives a multiplicative generator
    ff_gen = gen(ff);

    if verify_mult_gen
        is_conway(ff) || throw(error("The field constructed has no easily identifiable multiplicative generator"));
    end

    w = ff_gen^BigInt((BigInt(q)^2-1)/n);
    F = matrix(ff, [[w^(i*j) for i in 0:(n-1)] for j in 0:(n-1)]);

    Phi = F[D.+1,:];
    if return_gram
        conjugate_transpose(Phi)*Phi
    else
        Phi
    end
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


function etf_from_modular_hadamard(modHadamard::FqMatrix; return_gram::Bool=true)::FqMatrix
    # This construction comes from Theorem 2 in [4].
    # Constructs a dxd^2 ETF in case U, given a modular hadamard defined over 
    # a quadratic extension of a finite field such that x^2+1 does not split in 
    # the base field and d = 8 mod characteristic.

    d = size(modHadamard)[1];
    (d == size(modHadamard)[2]) || throw(DomainError(size(modHadamard),"Modular Hadamard must be square"));
    
    ff = modHadamard.base_ring;

    ## need root of x^2+1
    ## this means that in the base field -1 is not a square.
    (gcd(2, degree(ff)) == 2) || throw(DomainError(ff,"in Case U, the provided field must be finite and must be an extension of even degree"));
    
    p = characteristic(ff);
    abs_d = absolute_degree(ff);
    ((gcd(4,p-3)==4) && (gcd(4, abs_d)==2)) || throw(DomainError(ff,"The base field must have -1 be a non-square"));

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


function etf_from_singer_diff_set(D, p::T, k::T, r::T; verify_mult_gen::Bool=false, return_gram::Bool=true) where {T<:Union{Int,BigInt}}
    # This construction comes from Theorem 21 in [3].
    # Constructs a dxd^2 ETF in Case U, as long as p divides r-1 and r^2+r+1 divides p^k+1
    # and D is a Singer difference set. 
    # This construction does not verify that D is Singer difference set.

    (e, b) = is_power(r);
    (is_prime(b)) || throw(DomainError(r, "input r must be a prime power."));
    (is_prime(p)) || throw(DomainError(p, "input p must be a prime."));

    q = BigInt(p)^k;
    d = r^2+r+1;
    (BigInt(d) == BigInt(r)^2+r+1) || throw(DomainError((r), "Int overflow when computing for d=r^2+r+1, pick a (much much) smaller r"));

    (gcd(p, r-1)==p) || throw(DomainError((p,r), "p must divide r-1"));
    (gcd(d, q+1)==d) || throw(DomainError((p,k,r), "r^2+r+1 must divide p^k+1"));

    ff = GF(BigInt(q)^2);
    
    # According to Oscar.jl this is not guaranteed to be a multiplicative generator.
    # However since Oscar.jl is using Nemo.jl which is using FLINT
    # calling GF will attempt to instantiate a finite field using a Conways Polynomial, if possible.
    # It does not appear that Nemo.jl has provide wrapper functionality for checking if this is the case.
    # So all that can be done is to hope that the underlying field is using a Conway Polynomial, and gen gives a multiplicative generator
    ff_gen = gen(ff);

    ## Every the smallest examples of this construction have q^2-1 not being an integer, making this verification impossible.
    if verify_mult_gen
        is_conway(ff) || throw(error("The field constructed has no easily identifiable multiplicative generator"));
    end

    w = ff_gen^BigInt((BigInt(q)^2-1)/d);

    perm = [[i for i in 2:d]; 1];
    modulation = diagonal_matrix([ w^x for x in 0:(d-1)]);

    one_D = zeros(Int,(d,1));
    one_D[(D.+1),:] .= 1;
    one_D = matrix(ff, one_D);

    Phi = matrix(ff, zeros(Int,(d,d^2)));
    
    ind = 1;
    for s in 0:Int(d-1)
        for t in 0:Int(d-1)
            if s == 0 
                if t == 0
                    Phi[:,ind] = one_D
                else
                    Phi[:,ind] = modulation*Phi[:,ind-1]
                end
            else
                Phi[:,ind] = Phi[perm,ind-d]
            end
            ind += 1;
        end
    end

    if return_gram
        conjugate_transpose(Phi)*Phi
    else
        Phi
    end
end
