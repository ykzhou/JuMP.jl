addToExpression(x, y::Vector{Variable}, z) = addToExpression(x, z, y)

addToExpression(x::AffExpr, y::Number, z::Vector{Variable}) = addToExpression(x, y, convert(Vector{AffExpr},z))
addToExpression(x::QuadExpr, y::Number, z::Vector{Variable}) = addToExpression(x, y, convert(Vector{AffExpr},z))
addToExpression(x, y::Vector{Variable}, z::Vector{Variable}) = error()
addToExpression(x, y::Vector{Variable}, z::OneIndexedArray) = error()
addToExpression(x, y::OneIndexedArray, z::Vector{Variable}) = error()
addToExpression(x, y::OneIndexedArray, z::OneIndexedArray) = error()

addToExpression(x, y, z::Vector{Variable}) = addToExpression(x, y, convert(AffExpr,z))

addToExpression{T}(x, y, z::JuMPArray{T,1,true}) = addToExpression(x, y, z.innerArray)

addToExpression{T}(x, y::JuMPArray{T,1,true}, z) = addToExpression(x, y.innerArray, z)

function addToExpression(x::Number, y::Number, z::Vector{AffExpr})
    v = copy(z)
    for i in 1:length(v)
        t = v[i]
        for j in 1:length(t.vars)
            t.coeffs[j] *= y
        end
        t.constant += x
    end
    return v
end

function addToExpression(x::Number, y::Number, z::Vector{QuadExpr})
    v = copy(z)
    for i in 1:length(v)
        t = v[i]
        for j in 1:length(t.qvars1)
            t.qcoeffs[j] *= y
        end
        for k in 1:length(t.aff.vars)
            t.aff.coeffs *= y
        end
        t.aff.constant += x
    end
    return v
end

# ugly fallback, should find a better solution to make things work
function addToExpression{T<:JuMPTypes}(x::AffExpr, y::Number, z::Vector{T})
    (isempty(x.vars) && isempty(x.coeffs)) || error("Nonempty affine expression")
    return addToExpression(x.constant, y, z)
end

function addToExpression{T<:JuMPTypes}(x::QuadExpr, y::Number, z::Vector{T})
    (isempty(x.qvars1) && isempty(x.aff.vars)) || error("Nonempty affine expression")
    return addToExpression(x.aff.constant, y, z)
end

function addToExpression{R<:Number}(x::Vector{R}, y::Number, z::Vector{AffExpr})
    length(x) == length(z) || error("Incompatible dimensions")
    v = copy(z)
    for i in 1:length(v)
        t = v[i]
        for j in 1:length(t.vars)
            t.coeffs[j] *= y
        end
        t.constant += x[i]
    end
    return v
end

function addToExpression{R<:Number}(x::Vector{R}, y::Number, z::Vector{QuadExpr})
    length(x) == length(z) || error("Incompatible dimensions")
    v = copy(z)
    for i in 1:length(v)
        t = v[i]
        for j in 1:length(t.qvars1)
            t.qcoeffs[j] *= y
        end
        for k in 1:length(t.aff.vars)
            t.aff.coeffs[k] *= y
        end
        t.aff.constant += x[i]
    end
    return v
end

function addToExpression{T<:Union(AffExpr,QuadExpr),R<:Real}(x::Vector{T}, y::Number, z::Vector{R})
    length(x) == length(z) || error("Incompatible dimensions")
    for i in 1:length(x)
        x[i].constant += y*z[i]
    end
    return x
end

function addToExpression(x::Vector{AffExpr}, y::Number, z::Vector{AffExpr})
    length(x) == length(z) || error("Incompatible dimensions")
    for i in 1:length(x)
        append!(x[i].vars, z[i].vars)
        append!(x[i].coeffs, fill(y,length(z[i].vars)))
        x[i].constant += y*z[i].constant
    end
    return x
end

function addToExpression(x::Vector{AffExpr}, y::Number, z::Vector{QuadExpr})
    length(x) == length(z) || error("Incompatible dimensions")
    v = copy(z)
    for i in 1:length(v)
        t = v[i]
        for j in 1:length(v.qvars1)
            t.qcoeffs[j] *= y
        end
        for k in 1:length(v.aff.vars)
            t.aff.coeffs[k] *= y
        end
        append!(t.aff.vars, x[i].vars)
        append!(t.aff.coeffs, x[i].coeffs)
        t.constant += y*z[i].constant
    end
    return x
end

function addToExpression(x::Vector{QuadExpr}, y::Number, z::Vector{AffExpr})
    length(x) == length(z) || error("Incompatible dimensions")
    for i in 1:length(x)
        append!(x[i].aff.vars, z[i].vars)
        append!(x[i].aff.coeffs, fill(y,length(z[i].vars)))
        x[i].aff.constant += y*z[i].constant
    end
    return x
end

function addToExpression(x::Vector{QuadExpr}, y::Number, z::Vector{QuadExpr})
    length(x) == length(z) || error("Incompatible dimensions")
    for i in 1:length(x)
        append!(x[i].qvars1, z[i].qvars1)
        append!(x[i].qvars2, z[i].qvars2)
        st = length(x[i].qcoeffs)
        append!(x[i].qcoeffs, z[i].qcoeffs)
        for j in (st+1):length(x[i].qcoeffs)
            x[i].qcoeffs[j] *= y
        end
        append!(x[i].aff.vars, z[i].aff.vars)
        st = length(x[i].aff.coeffs)
        append!(x[i].aff.coeffs, z[i].aff.coeffs)
        for j in (st+1):length(x[i].aff.coeffs)
            x[i].aff.coeffs[j] *= y
        end
        x[i].aff.constant += y*z[i].constant
    end
    return x
end
