using FiniteFieldFrames
using Oscar
using Test

function verify_binder(gram, case, binder, non_binder=nothing)
    Phi = reconstruct_frame_from_gram(gram, case=case);

    if case == :O
        if (transpose(Phi)*Phi != Phi.base_ring.(gram))
            return False
        end
    elseif case == :U
        if (conjugate_transpose(Phi)*Phi != Phi.base_ring.(gram))
            return False
        end
    else
    end

    for ii in 1:(binder.size[1])
        simp_ind = findall(x->x==1, binder[ii,:]);
        simplex_in_Phi = Phi[:, simp_ind];
        d = rank(simplex_in_Phi);
        etf_bool, a, b, c, dd = is_ETF(gram[simp_ind, simp_ind], case=case);

        if !etf_bool || d != dd
            return false
        end
    end

    if non_binder != nothing
        for ii in 1:(non_binder.size[1])
            nsimp_ind = findall(x->x==1, non_binder[ii,:]);
            nsimplex_in_Phi = Phi[:, nsimp_ind];
            d = rank(nsimplex_in_Phi);
            etf_bool, a, b, c, dd = is_ETF(gram[nsimp_ind, nsimp_ind], case=case);

            if etf_bool && d != dd && nsimp_ind.size[1] == d+1
                return false
            end
        end 
    end

    return true
end

@testset "FiniteFieldFrames.jl" begin
    gram = etf_from_triangular_graph(12, 5);
    Phi = reconstruct_frame_from_gram(gram, case=:O);
    @test transpose(Phi)*Phi == Phi.base_ring.(gram);

    gram2 = real_etf_to_case_O(real_dx2d_etf_from_prime_power(3^4), Int((3^4+1)/2), 7);
    Phi2 = reconstruct_frame_from_gram(gram2, case=:O);
    @test transpose(Phi2)*Phi2 == Phi2.base_ring.(gram2);

    Psi = etf_from_pmod_diff_set([0,4,6,7,8,11,13], 14, 27);
    gram3 = conjugate_transpose(Psi)*Psi;
    Phi3 = reconstruct_frame_from_gram(gram3, case=:U);
    @test conjugate_transpose(Phi3)*Phi3 == Phi3.base_ring.(gram3);
end

@testset "Constructions.jl" begin
    gram = etf_from_triangular_graph(12, 5);
    etf_bool, a, b, c, d = is_ETF(gram, case=:O);
    @test (etf_bool && d == 12);

    gram2 = real_etf_to_case_O(real_dx2d_etf_from_prime_power(3^4), Int((3^4+1)/2), 7);
    etf_bool2, a2, b2, c2, d2 = is_ETF(gram2, case=:O);
    @test (etf_bool2 && d2 == Int((3^4+1)/2));

    Phi = etf_from_pmod_diff_set([0,4,6,7,8,11,13], 14, 27);
    gram2 = conjugate_transpose(Phi)*Phi;
    etf_bool3, a3, b3, c3, d3 = is_ETF(gram2, case=:U);
    @test (etf_bool3 && d3 == 7);
end

@testset "BinderFinder.jl" begin
    # These serve only to test that the output of binder_finder is a set of simplicies.
    # This does not check that it found the full binder.
    gram = etf_from_triangular_graph(7, 5);
    @test verify_binder(gram, :O, binder_finder(gram, case=:O))

    gram2 = etf_from_triangular_graph(7, 7);
    @test verify_binder(gram2, :O, binder_finder(gram2, case=:O))

    gram3 = real_etf_to_case_O(real_dx2d_etf_from_prime_power(5^2), Int((5^2+1)/2), 7);
    @test verify_binder(gram3, :O, binder_finder(gram3, case=:O))
end
