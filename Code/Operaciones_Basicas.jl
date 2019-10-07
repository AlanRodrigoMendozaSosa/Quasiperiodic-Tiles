function prod_Punto(A,B)
    return A[1]*B[1] + A[2]*B[2]
end

function norma_Vector(A)
    return sqrt(A[1]^2 + A[2]^2)
end

function vector_Ortogonal(A)
    return [A[2], -A[1]]
end