function prod_Punto(A,B)
    return A[1]*B[1] + A[2]*B[2]
end

function norma_Vector(A)
    return sqrt(A[1]^2 + A[2]^2)
end

function vector_Ortogonal(A)
    return [A[2], -A[1]]
end

#Función que nos genera un punto aleatorio en un cuadrado de semilado SL centrado en el origen.
function punto_Arbitrario(SL)
    #Definimos la variable Punto que contendrá las coordenadas del punto de interés
    Punto = zeros(Float64,2);

    #Llaves para determinar el cuadrante donde estará el punto
    x = rand();
    y = rand();

    if (x > 0.5) && (y > 0.5)
        Punto = [rand()*SL, rand()*SL];
    elseif (x > 0.5) && (y < 0.5)
        Punto = [rand()*SL, -rand()*SL];
    elseif (x < 0.5) && (y > 0.5)
        Punto = [-rand()*SL, rand()*SL];
    elseif (x < 0.5) && (y < 0.5)
        Punto = [-rand()*SL, -rand()*SL];
    end

    return Punto
end