export load_matrix_from_file

function load_matrix_from_file(file_loc::String, d::Int, n::Int)::Union{Matrix{BigFloat}, Matrix{Complex{BigFloat}}}
    m = readdlm(file_loc, BigFloat);
    real_mat = reshape(m[1:d*n], (d,n));
    comp_mat = reshape(m[d*n+1:2*d*n], (d,n));
    if iszero(comp_mat)
        real_mat
    else
        real_mat + im*comp_mat
    end
end