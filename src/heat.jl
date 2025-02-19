using MMJMesh
using MMJMesh.Geometries
using MMJMesh.Mathematics

# Gauss points for numerical integration
const gauss_w = ones(4)
const gauss_p = 1 / sqrt(3) * [[-1, -1], [1, -1], [1, 1], [-1, 1]]

# Shape functions and gradients
const N = MappingFromComponents(nodalbasis(makeelement(:lagrange, QHat, k=1))...)
const ∇ξN = MMJMesh.Mathematics.TransposeMapping(jacobian(N))

# Element matrix
function heat_ke(λ)
    function ke_func(e::Face{2, 2, 4})
        Ke = zeros(4, 4)
        jF = jacobian(parametrization(geometry(e)))
        for (ξ, w) ∈ zip(gauss_p, gauss_w)
            J = jF(ξ)
            ∇ₓN = inv(J') * ∇ξN(ξ)
            B = ∇ₓN
            Ke += w * λ * B' * B * det(J)
        end
        return Ke
    end
    return ke_func
end

# Element vector
function heat_re(w)
    function re_func(e)
        re = zeros(4)
        jF = jacobian(parametrization(geometry(e)))
        for ξ ∈ gaussPoints
            J = jF(ξ)
            re += w * N(ξ) * det(J)
        end
        return re
    end
    return re_func
end
