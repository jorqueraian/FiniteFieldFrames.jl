using Oscar
using LinearAlgebra
using Printf


## Note: I will do everything on the level of gram matrices.

function is_frame(gram, case::String)
    n = size(gram)[1];
    (n == size(gram)[2]) || throw(DomainError(size(gram),"gram must be square"));
    d = rank(gram);

    if case == "O"
        if iszero(gram - transpose(gram))
            return (true, d)
        else
            return (false, nothing)
        end
    elseif case == "U"
        (degree(gram.base_ring) == 2) || throw(DomainError(gram.base_ring,"in Case U, the provided field must be finite and must be a degree 2 extension"));
        if iszero(gram - conjugate_transpose(gram))
            return (true, d)
        else
            return (false, nothing)
        end
    else
        throw(DomainError(case,"I can only do Case O and U at the moment"));
    end
end

function case_O_frame_discr_is_square(gram)
    n = size(gram)[1];
    (n == size(gram)[2]) || throw(DomainError(size(gram),"gram must be square"));
    frame_bool, d = is_frame(gram, "O");
    (frame_bool) || throw(DomainError(gram,"gram must be symmetric"));
    
    for inds in combinations(1:n, d)
        sub_gram = gram[inds,inds];
        if rank(sub_gram) == d
            return is_square(det(sub_gram))
        end
    end
end

function is_equiangular(gram, case::String)
    if case == "O"
        gram_mat_modulus_sqrd = gram.^2;
    elseif case == "U"
        (degree(gram.base_ring) == 2) || throw(DomainError(gram.base_ring,"in Case U, the provided field must be finite and must be a degree 2 extension"));
        gram_mat_modulus_sqrd = norm.(gram);
    else
        throw(DomainError(case,"I can only do Case O and U at the moment"));
    end

    n = size(gram)[1];
    (n == size(gram)[2]) || throw(DomainError(size(gram),"gram must be square"));

    a = gram[1,1];
    a_sqrd = gram_mat_modulus_sqrd[1,1];
    b = gram_mat_modulus_sqrd[1,2];

    test_mat = (gram_mat_modulus_sqrd.-b) - (a_sqrd-b)*matrix(gram_mat_modulus_sqrd.base_ring, I[1:n,1:n])

    if iszero(test_mat)
        return (true, a, b)
    else
        return (false, nothing, nothing)
    end
end


function is_frame_tight(gram)
    # Note that this algorithm assumes that gram is the gram of a frame.
    # if tight, will return tight constant and the rank.
    n = size(gram)[1];
    (n == size(gram)[2]) || throw(DomainError(size(gram),"gram must be square"));

    gram_sqrd = gram * gram;

    d = rank(gram)

    if iszero(gram_sqrd)
        # c=0 iff G^2 = 0*G
        return (true, gram.base_ring(0))
    elseif !iszero(gram[1,1])
        c = gram_sqrd[1,1] * inv(gram[1,1])
    elseif !iszero(gram[2,1])
        c = gram_sqrd[2,1] * inv(gram[2,1])
    else
        throw(DomainError(gram,"gram has to many zeros, and I haven't written cases to deal with that yet"));
    end

    if iszero(gram_sqrd - (c*gram))
        return (true, c) 
    else
        return (false, nothing)
    end
end

function is_ETF(gram, case)
    frame_bool, d = is_frame(gram, case);
    if !frame_bool
        return (false, nothing, nothing, nothing, nothing)
    else
        equiangular_bool, a, b = is_equiangular(gram, case);
        tight_bool, c = is_frame_tight(gram);
        return ((equiangular_bool && tight_bool), a, b, c, d);
    end
end


function reconstruct_frame_from_gram(gram, case)
    # see Thm 3.13 and 3.15 of citation needed.
    n = size(gram)[1];
    (n == size(gram)[2]) || throw(DomainError(size(gram),"gram must be square"));

    frame_bool, d = is_frame(gram, case);
    frame_bool || throw(DomainError(size(gram),"gram must be the Gram matrix of a frame"));

    ff = gram.base_ring;

    if case == "O"
        q = 1;
    elseif case == "U"
        (degree(ff) == 2) || throw(DomainError(ff,"in Case U, the provided field must be finite and must be a degree 2 extension"));
        q = Int(sqrt(size(ff)));
    end

    A = matrix(ff, I[1:n,1:n]);
    for j in 2:n 
        for i in 1:(j-1) 
            A[i,j] = gram[i,j] - sum([A[k,i]^q*A[k,j] for k in 1:(i-1)], init=ff(0));
        end
    end
    
    if case == "O"
        for i in 1:n
            num = gram[i,i]- transpose(A[:,i])*A[:,i]
            if !is_square(num)
                Kx, x = ff["x"];
                ff, ffa = finite_field(x^2-num, "a");
                A = ff.(A)
                gram = ff.(gram)
            end
        end
    end
    B = diagonal_matrix(ff(0), n, n);
    
    for i in 1:n 
        if case == "O"
            B[i,i] = sqrt(gram[i,i]- transpose(A[:,i])*A[:,i])
        else
            Kx, x = ff["x"];
            B[i,i] = roots(x^(q+1)-gram[i,i]+conjugate_transpose(A[:,i])*A[:,i])[1]
        end
    end

    Psi = [A; B]
    rank(Psi) == n || throw(error("Something did not work. Constructed matrix Psi should have rank n, but it does not."));
    
    if case == "O"
        form = symmetric_form(gram);
        m = symmetric_form(diagonal_matrix([[ff(1) for i in 1:d]; [ff(0) for i in (d+1):n]]));
        #(cong_bool, C) =  is_congruent(m, form);
        # C * gram * C' = daig(1 1 1 1... 0 0 0....)
        # D = C^{-1}
        (cong_bool, D) =  is_congruent(form, m);
        cong_bool || throw(error("Frankly I thought this was impossible."));

        Phi = transpose(D)[1:d, :]
        (rank(Phi) == d) || throw(error("Something broke and I dont know why"));
    elseif case == "U"
        form = hermitian_form(gram);
        m = hermitian_form(diagonal_matrix([[ff(1) for i in 1:d]; [ff(0) for i in (d+1):n]]));
        (cong_bool, D) =  is_congruent(form, m);
        Phi = conjugate_transpose(D)[1:d, :]
        (rank(Phi) == d) || throw(error("Something broke and I dont know why"));
    end
    Phi
end