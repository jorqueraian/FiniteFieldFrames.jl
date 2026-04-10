export binder_finder

import Base.Threads: @threads


function binder_finder(gram::FqMatrix, case::String, verbose::Bool=false)::Matrix{Int}
    ####################
    ###### Binder Finder
    ###### Citation: "Equiangular tight frames that contain regular simplices,"
    ###### Matthew Fickus, John Jasper, Emily J. King, Dustin G. Mixon, 2017
    ###### Modified for finite fields by Ian Jorquera
    ###### This is not a true finite field binder finder
    ######
    ###### Code created: July 2016
    ###### Last updated: April 1st, 2026
    #####################
    # Notes: this may not return a complete binder, better writeup coming soon

    n = size(gram)[1];
    (n == size(gram)[2]) || throw(DomainError(size(gram),"gram must be square"));

    ff = gram.base_ring;
    char = Int(characteristic(ff));
    etf_bool, a, b, c, d = is_ETF(gram, case);
    (etf_bool) || throw(DomainError(gram,"gram must be the gram of and ETF"));

    (!iszero(a)) || throw(DomainError(a,"need a != 0 for this implementation, try specifying a particular value for s as the first parameter"));

    ##possible s's
    allss = [s for s in 2:d if a^2 == ff(s^2)*b];
    bad_s_inds1 = findall(s->(gcd(char, s)==char), allss);
    bad_s_inds2 = findall(s->(gcd(char, s+1)==char), allss);
    good_s_inds = findall(s->(gcd(char, s)==1 && gcd(char, s+1)==1), allss);
    if verbose
        if bad_s_inds2.size[1] > 0
            println("Skipping the following s's:")
            println(allss[bad_s_inds2])
            println("This is becasue the characteristic of ff divides s+1, and so no s-simplices can exist. ")
        end
        if bad_s_inds1.size[1] > 0
            # this wont happen, because im already excluding a = 0.
            println("Skipping the following s's:")
            println(allss[bad_s_inds1])
            println("This is becasue the characteristic of ff divides s, and so binder finder does not apply. This means the resulting binder may be incomplete.")
        end
    end
    allss = allss[good_s_inds];

    ss1 = [s for s in allss if ff(s)==ff(allss[1])];
    ss2 = [s for s in allss if !(s in ss1)];

    simplices = [];
    for ss in [ss1, ss2]
        if ss.size[1] == 0
            continue
        end
        if verbose
            print("Checking simplex dims: ");
            println(ss);
        end

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
            if verbose
                @printf("Found %i simplices with s=2\n", actual_simps.size[1])
            end
            append!(simplices, jTuple[actual_simps,:]')
            jTuple = jTuple[[i for i in 1:(jTuple.size[1]) if !(i in actual_simps)],:]
        end

        for jTupleSize in 3:maximum(ss)  # should this be s? or is my python code wrong??
            if verbose
                @printf("Looking for subsets of %i with equal TPs\n", jTupleSize+1)
            end
            
            # creating an array for each iteration of the next loop.
            threadedNewjTuples = [Array{Int}(undef, 0, 0)  for ii in 1:(jTuple.size[1])];
            @threads for ii in 1:(jTuple.size[1])
                Indices = findall(x->x==1, jTuple[ii,:]);
                Indicator = (sum(Int.(jTuple[:,Indices] * (ones(Int, (jTupleSize,jTupleSize))-I[1:jTupleSize,1:jTupleSize]) .== (jTupleSize-1)),dims=2))' * jTuple
                Indicator[Indices] .= 0
                NewIndices = findall(x -> x == jTupleSize, Indicator[1,:])
                A = kron(ones(Int, (NewIndices.size[1],1)),jTuple[ii,:]')
                A[:,NewIndices] = I[1:NewIndices.size[1],1:NewIndices.size[1]]

                @inbounds threadedNewjTuples[ii] = A
                @inbounds jTuple[ii,:] = zeros(Int, (1,jTuple.size[2]))
            end
            jTuple = vcat(threadedNewjTuples...)
            
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
                if verbose
                    @printf("Found %i simplices with s=%i\n", actual_simps.size[1], s)
                end
                append!(simplices, jTuple[actual_simps,:]')
                jTuple = jTuple[[i for i in 1:(jTuple.size[1]) if !(i in actual_simps)],:]
            end
        end
    end
    reshape(simplices, (n, Int(simplices.size[1]/n)))'
end


function binder_finder(s::Int, gram::FqMatrix, case::String)::Matrix{Int}
    # While the other version picks values of s automatically,
    # you can specify a specific s here.
    # This also allows a != 0, BUT it will give you the warning 
    # that there is no way to verify what this code found is correct.

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
        Phi = reconstruct_frame_from_gram(gram, case);

        flat_binder = [];
        for inds in combinations(1:n, k)
            maybe_simplex = Phi[:, inds];
            sub_gram = gram[inds,inds];

            # verify gram of a frame for a s-dim space
            if rank(sub_gram) == s && rank(sub_gram) == rank(maybe_simplex)
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