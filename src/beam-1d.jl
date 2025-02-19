function beam_ke(EI)
    function ke_func(e)
        l = length(e)

        return [
            12EI/l^3 6EI/l^2 -12EI/l^3 6EI/l^2
            6EI/l^2 4EI/l -6EI/l^2 2EI/l
            -12EI/l^3 -6EI/l^2 12EI/l^3 -6EI/l^2
            6EI/l^2 2EI/l -6EI/l^2 4EI/l
        ]
    end
    return ke_func
end

function beam_re(q)
    function re_func(e)
        l = length(e)
        
        return [
            0.5 * l * q
            0.0833333333333333 * l^2 * q
            0.5 * l * q
            -0.0833333333333333 * l^2 * q
        ]
    end
    return re_func
end


