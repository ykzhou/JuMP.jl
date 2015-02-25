#  Copyright 2015, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

typealias VectTypes Union(JuMPTypes,Real)

function addToExpression{T<:VectTypes}(x::AffExpr, y::Real, z::Array{T})
    (isempty(x.vars) && isempty(x.coeffs)) || error("Cannot add an affine expression to Array{$T}")
    return x.constant + y.*z
end
function addToExpression{T<:VectTypes}(x::AffExpr, y::Array{T}, z::Real)
    (isempty(x.vars) && isempty(x.coeffs)) || error("Cannot add an affine expression to Array{$T}")
    return x.constant + y.*z
end
function addToExpression{T<:VectTypes,R<:VectTypes}(x::AffExpr, y::Array{T}, z::Array{R})
    (isempty(x.vars) && isempty(x.coeffs)) || error("Cannot add an affine expression to Array{$T}")
    return x.constant + y*z
end

addToExpression{T<:VectTypes,R<:VectTypes,S<:VectTypes}(x::Array{T},y::R,z::S) = (x .+ y*z)
addToExpression{T<:VectTypes,R<:VectTypes,S<:VectTypes}(x::Array{T},y::Array{R},z::S) = (x .+ y.*z)
addToExpression{T<:VectTypes,R<:VectTypes,S<:VectTypes}(x::Array{T},y::R,z::Array{S}) = (x .+ y.*z)
addToExpression{T<:VectTypes,R<:VectTypes,S<:VectTypes}(x::Array{T},y::R,z::Array{S}) = (x .+ y.*z)
addToExpression{T<:VectTypes,R<:VectTypes,S<:VectTypes}(x::Array{T},y::Array{R},z::Array{S}) = (x .+ y*z)

addToExpression{T<:Real,R<:VectTypes,S<:VectTypes}(x::T,y::Array{R},z::S) = (x .+ y.*z)
addToExpression{T<:Real,R<:VectTypes,S<:VectTypes}(x::T,y::R,z::Array{S}) = (x .+ y.*z)

addToExpression(x, y, z::OneIndexedArray) = addToExpression(x, y, z.innerArray)
addToExpression(x, y::OneIndexedArray, z::OneIndexedArray) = addToExpression(x, y.innerArray, z.innerArray)
addToExpression(x, y::OneIndexedArray, z) = addToExpression(x, y.innerArray, z)

addToExpression(x, y::OneIndexedArray, z::SparseMatrixCSC) = addToExpression(x, y.innerArray, full(z))
addToExpression(x, y, z::SparseMatrixCSC) = addToExpression(x, y, full(z)) # lol
addToExpression(x, y::SparseMatrixCSC, z::SparseMatrixCSC) = addToExpression(x, full(y), full(z))
addToExpression(x, y::SparseMatrixCSC, z::OneIndexedArray) = addToExpression(x, full(y), z.innerArray)
addToExpression(x, y::SparseMatrixCSC, z) = addToExpression(x, full(y), z)
