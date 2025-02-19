dofs(n, indices) = n == 1 ? indices : vec([n * i - j for j = n-1:-1:0, i = indices])

function assemble_kr(m, nf=1)
    N = nf * nnodes(m)
    K = zeros(N, N)
    r = zeros(N)

    for e = entities(m)
        I = dofs(nf, nodeindices(e))
        hasdata(e, :ke_func) && (K[I, I] += data(e, :ke_func)(e))
        hasdata(e, :re_func) && (r[I] += data(e, :re_func)(e))
    end

    return sparse(K), r
end

function assemble_kmr(m, nf=1)
    N = nf * nnodes(m)
    K = zeros(N, N)
    M = zeros(N, N)
    r = zeros(N)

    for e = entities(m)
        I = dofs(nf, nodeindices(e))
        hasdata(e, :ke_func) && (K[I, I] += data(e, :ke_func)(e))
        hasdata(e, :me_func) && (M[I, I] += data(e, :me_func)(e))
        hasdata(e, :re_func) && (r[I] += data(e, :re_func)(e))
    end

    return sparse(K), sparse(M), r
end

function apply_dirichlet_bcs!(K, r, dofs, v=0)
    K[dofs, :] .= 0
    r[dofs] .= v
    K[diagind(K)[dofs]] .= 1
    nothing
end
