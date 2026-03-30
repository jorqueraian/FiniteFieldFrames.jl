using Oscar
using LinearAlgebra


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

function binder_finder(s::Int, gram, case::String)::Matrix{Int}
    ####################
    # TODO Update documentation
    ###### Binder Finder
    ###### Code to find the embedded simplices (ETFs for their spans with one
    ###### more vector than their rank) in a given ETF.
    ######
    ###### Input: Phi, the synthesis matrix of an equiangular tight frame.
    ######
    ###### Output: Binder, the incidence matrix (rows = simplices, cols = frame
    ###### vectors) of the simplices embedded in the given ETF
    ######
    ###### Citation: "Equiangular tight frames that contain regular simplices,"
    ###### Matthew Fickus, John Jasper, Emily J. King, Dustin G. Mixon, 2017
    ###### Modified for finite fields by Ian Jorquera
    ###### This is not a true finite field binder finder
    ######
    ###### Code created: July 2016
    ###### Last updated: March 27th, 2026
    #####################

    n = size(gram)[1];
    k = s+1;
    (n == size(gram)[2]) || throw(DomainError(size(gram),"gram must be square"));
    ff = gram.base_ring;
    etf_bool, a, b, c, d = is_ETF(gram, case);
    (etf_bool) || throw(DomainError(gram,"gram must be the gram of and ETF"));

    (a^2 == ff(s^2)*b) || throw(DomainError((a^2, ff(s^2)*b),"Must satisfy that a^2=s^2*b which the inputs do not."));

    if !iszero(a)
        (!iszero(ff(s)) && is_square(a*inv(ff(s)))) || throw(DomainError((a,s),"need a/s to be a square"));
        neg_a3bys3 = ff(-1)*(a^3)*inv(ff(s^3))

        TripleCounter = 0;
        mat_list = [];
        for ii in 1:(n-2)
            for jj in (ii+1):(n-1)
                for kk in (jj+1):n
                    if gram[ii,jj]*gram[jj,kk]*gram[kk,ii] == neg_a3bys3
                        TripleCounter += 1;
                        row_lst = zeros(Int, (1,n));
                        row_lst[ii] = 1;
                        row_lst[jj] = 1;
                        row_lst[kk] = 1;
                        append!(mat_list, row_lst)
                    end
                end
            end
        end

        jTuple = Int.(reshape(mat_list, (n, Int(mat_list.size[1]/n)))');
        NewjTuple = nothing;

        for jTupleSize in 3:s  # should this be s? or is my python code wrong??
            for ii in 1:(jTuple.size[1])
                Indices = findall(x->x==1, jTuple[ii,:]);
                Indicator = (sum(Int.(jTuple[:,Indices] * (ones(Int, (jTupleSize,jTupleSize))-I[1:jTupleSize,1:jTupleSize]) .== (jTupleSize-1)),dims=2))' * jTuple
                Indicator[Indices] .= 0
                NewIndices = findall(x -> x == jTupleSize, Indicator[1,:])
                A = kron(ones(Int, (NewIndices.size[1],1)),jTuple[ii,:]')
                A[:,NewIndices] = I[1:NewIndices.size[1],1:NewIndices.size[1]]

                # A will contains set of jTupleSize + 1 vector with equal triple produ cts. if jtupleSize is not 

                if NewjTuple == nothing || ii == 1
                    NewjTuple = A;
                else
                    NewjTuple = [NewjTuple; A];
                end
                jTuple[ii,:] = zeros(Int, (1,jTuple.size[2]))
            end
            jTuple = NewjTuple
        end 
        
        actual_simps = []
        for ii in 1:(jTuple.size[1])
            simp = findall(x->x==1, jTuple[ii,:]);
            fixed = simp[1];
            non_simp = findall(x->x==0, jTuple[ii,:]);
            
            tp_sums = sum(reshape([gram[fixed,kk]*gram[kk,jj]*gram[jj,fixed] for jj in simp for kk in non_simp], (non_simp.size[1], simp.size[1])), dims=2);
            if iszero(tp_sums .- (ff(s+1)*inv(ff(s))*a*b))
                append!(actual_simps, ii);
            end
        end
        jTuple[actual_simps,:]
    else
        # Probably a better way to do this
        print("Warning: this implementation is slow and does not guarantee correct outputs when a != 0. Sorry\n\n")
        flat_binder = [];
        for inds in combinations(1:n, k)
            sub_gram = gram[inds,inds];

            # This is a bit of an issue,
            # it could be the case that the k vectors selected in the frame are LI
            # but are degenerate, so just checking the rank of the gram
            # is not enough. Not sure how to deal with this atm. But for now
            # just assume this isnt an issue.
            if rank(sub_gram) == s
                tight_bool, c = is_frame_tight(sub_gram);
                if tight_bool
                    simplex_row = zeros(Int, (n, 1));
                    simplex_row[inds] .= 1;
                    append!(flat_binder, simplex_row)
                end
            end
        end

        reshape(flat_binder, (n,Int(flat_binder.size[1]/n)))'
    end
end


function binder_finder(gram, case::String)::Matrix{Int}
    ####################
    ###### Binder Finder
    ###### Code to find the embedded simplices (ETFs for their spans with one
    ###### more vector than their rank) in a given ETF.
    ######
    ###### Input: Phi, the synthesis matrix of an equiangular tight frame.
    ######
    ###### Output: Binder, the incidence matrix (rows = simplices, cols = frame
    ###### vectors) of the simplices embedded in the given ETF
    ######
    ###### Citation: "Equiangular tight frames that contain regular simplices,"
    ###### Matthew Fickus, John Jasper, Emily J. King, Dustin G. Mixon, 2017
    ###### Modified for finite fields by Ian Jorquera
    ###### This is not a true finite field binder finder
    ######
    ###### Code created: July 2016
    ###### Last updated: March 27th, 2026
    #####################

    n = size(gram)[1];
    (n == size(gram)[2]) || throw(DomainError(size(gram),"gram must be square"));

    ff = gram.base_ring;
    etf_bool, a, b, c, d = is_ETF(gram, case);
    (etf_bool) || throw(DomainError(gram,"gram must be the gram of and ETF"));

    (!iszero(a)) || throw(DomainError(a,"need a != 0 for this implementation, try specifying a particular value for s as the first parameter"));

    ##possible s's
    # Note that this also depends on 
    allss = [s for s in 2:d if a^2 == ff(s^2)*b];
    ss1 = [s for s in allss if ff(s)==ff(allss[1])];
    ss2 = [s for s in allss if !(s in ss1)];

    simplices = [];
    for ss in [ss1, ss2]
        jTuple = [];
        neg_a3bys3 = ff(-1)*(a^3)*inv(ff(ss[1]^3))

        TripleCounter = 0;
        mat_list = [];
        for ii in 1:(n-2)
            for jj in (ii+1):(n-1)
                for kk in (jj+1):n
                    if gram[ii,jj]*gram[jj,kk]*gram[kk,ii] == neg_a3bys3
                        TripleCounter += 1;
                        row_lst = zeros(Int, (1,n));
                        row_lst[ii] = 1;
                        row_lst[jj] = 1;
                        row_lst[kk] = 1;
                        append!(mat_list, row_lst)
                    end
                end
            end
        end

        jTuple = Int.(reshape(mat_list, (n, Int(mat_list.size[1]/n)))');

        # if s could be 2, we need to check that none of the triples are simplices.
        if 2 in ss
            s = 2
            actual_simps = []
            for ii in 1:(jTuple.size[1])
                simp = findall(x->x==1, jTuple[ii,:]);
                fixed = simp[1];
                non_simp = findall(x->x==0, jTuple[ii,:]);
                
                tp_sums = sum(reshape([gram[fixed,kk]*gram[kk,jj]*gram[jj,fixed] for jj in simp for kk in non_simp], (non_simp.size[1], simp.size[1])), dims=2);
                if iszero(tp_sums .- (ff(s+1)*inv(ff(s))*a*b))
                    append!(actual_simps, ii);
                end
            end
            append!(simplices, jTuple[actual_simps,:]')
            jTuple = jTuple[[i for i in 1:(jTuple.size[1]) if !(i in actual_simps)],:]
        end

        NewjTuple = [];
        for jTupleSize in 3:maximum(ss)  # should this be s? or is my python code wrong??
            for ii in 1:(jTuple.size[1])
                Indices = findall(x->x==1, jTuple[ii,:]);
                Indicator = (sum(Int.(jTuple[:,Indices] * (ones(Int, (jTupleSize,jTupleSize))-I[1:jTupleSize,1:jTupleSize]) .== (jTupleSize-1)),dims=2))' * jTuple
                Indicator[Indices] .= 0
                NewIndices = findall(x -> x == jTupleSize, Indicator[1,:])
                A = kron(ones(Int, (NewIndices.size[1],1)),jTuple[ii,:]')
                A[:,NewIndices] = I[1:NewIndices.size[1],1:NewIndices.size[1]]

                if NewjTuple.size[1] == 0 || ii == 1
                    NewjTuple = A;
                else
                    NewjTuple = [NewjTuple; A];
                end
                jTuple[ii,:] = zeros(Int, (1,jTuple.size[2]))
            end
            jTuple = NewjTuple
            
            if jTupleSize in ss
                s = jTupleSize
                actual_simps = []
                for ii in 1:(jTuple.size[1])
                    simp = findall(x->x==1, jTuple[ii,:]);
                    fixed = simp[1];
                    non_simp = findall(x->x==0, jTuple[ii,:]);
                    
                    tp_sums = sum(reshape([gram[fixed,kk]*gram[kk,jj]*gram[jj,fixed] for jj in simp for kk in non_simp], (non_simp.size[1], simp.size[1])), dims=2);
                    if iszero(tp_sums .- (ff(s+1)*inv(ff(s))*a*b))
                        append!(actual_simps, ii);
                    end
                end
                append!(simplices, jTuple[actual_simps,:]')
                jTuple = jTuple[[i for i in 1:(jTuple.size[1]) if !(i in actual_simps)],:]
            end
        end
    end
    reshape(simplices, (n, Int(simplices.size[1]/n)))'
end