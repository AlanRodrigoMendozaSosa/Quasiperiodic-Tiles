#La función determina, fijando los vectores estrella Ej y Ek, fijando los enteros Nj y Nk, y para algún conjunto
#de alfas (uno por cada vector estrella), los cuatros puntos de la red dual
function cuatro_Regiones(J, K, Nj, Nk, Vectores_Estrella, Arreglo_Alfas)
    VecJ = Vectores_Estrella[J];
    VecK = Vectores_Estrella[K];
    
    #No podemos considerar a los mismos vectores
    if J == K
        error("Los índices J y K deben ser diferentes")
    else        
        #Definimos los dos vectores con los que se consigue la intersección en la malla generada por los vectores
        #estrella que estamos considerando.
        Ej = Vectores_Estrella[J]; #Creo que son innecesarios, ya están definidos arriba
        Ek = Vectores_Estrella[K]; #Creo que son innecesarios, ya están definidos arriba
        
        #Obtenemos los vectores ortogonales a estos dos vectores
        EjOrt = vector_Ortogonal(Ej);
        EkOrt = vector_Ortogonal(Ek);
        
        #Definimos los valores reales con los que se crearon las rectas ortogonales a cada vector Ej y Ek para tomar
        #la intersección
        Factor_Ej = Nj + Arreglo_Alfas[J];
        Factor_Ek = Nk + Arreglo_Alfas[K];
        
        #Obtenemos el área que forman los dos vectores Ej y Ek
        AreaJK = Ej[1]*Ek[2] - Ej[2]*Ek[1];
        
        #Definimos lo que será el vértice del Dual denominado t^{0} en la teoría
        Punto_Cero = Nj*Ej + Nk*Ek;
        
        for i in 1:length(Vectores_Estrella)
            if i == J || i == K
                nothing
            else
                Factor_Ei = (Factor_Ej/AreaJK)*(prod_Punto(EkOrt, Vectores_Estrella[i])) - (Factor_Ek/AreaJK)*(prod_Punto(EjOrt, Vectores_Estrella[i]));
                #println(i,floor(Factor_Ei - Arreglo_Alfas[i]))
                Punto_Cero += (floor(Factor_Ei - Arreglo_Alfas[i]))*Vectores_Estrella[i];
            end
        end
        
        #Obtenemos los otros tres puntos asociados al punto t^{0}
        Punto_Uno = Punto_Cero - Ej;
        Punto_Dos = Punto_Cero - Ej - Ek;
        Punto_Tres = Punto_Cero - Ek;
        
        #Agregamos la información correspondiente a cómo se obtuvieron estos vértices
        Info = [J,K,Nj,Nk];
        
        return Punto_Cero, Punto_Uno, Punto_Dos, Punto_Tres, Info
    end
end

#La función determina, fijando los vectores estrella Ej y Ek, y para algún conjunto de alfas (uno por cada vector 
#estrella), los puntos de la red dual corriendo los valores enteros Nj y Nk desde -N hasta N
function puntos_Dual_JK(J, K, N, Vectores_Estrella, Arreglo_Alfas)
    #Definimos el arreglo que contendra los vectores asociados a cada punto
    Puntos_Red_Dual = [];
    Informacion_Puntos_Red_Dual = [];
    
    for Nj in -N:N
        for Nk in -N:N
            try #Ponemos el Try-Catch para que no debamos separar los casos en que J == K, el error lo maneja autom.
                t0, t1, t2, t3, info = cuatro_Regiones(J, K, Nj, Nk, Vectores_Estrella, Arreglo_Alfas)
                push!(Puntos_Red_Dual, t0);
                push!(Puntos_Red_Dual, t1);
                push!(Puntos_Red_Dual, t2);
                push!(Puntos_Red_Dual, t3);
                push!(Informacion_Puntos_Red_Dual, info);
            catch y
                nothing
            end
        end
    end
    
    return Puntos_Red_Dual, Informacion_Puntos_Red_Dual
end

#Ahora hay que barrer todas las posibles combinaciones de los vectores Ej y Ek
function puntos_Dual(N, Vectores_Estrella, Arreglo_Alfas)
    #Definimos el arreglo que contendrá todos los puntos de la red Dual
    Puntos_Duales = [];
    Informacion_Duales = [];
    
    #Corramos los índices J y K de los vectores Ej y Ek
    for J in 1:length(Vectores_Estrella)
        for K in (J+1):length(Vectores_Estrella)
            Puntos_JK = puntos_Dual_JK(J, K, N, Vectores_Estrella, Arreglo_Alfas)
            push!(Puntos_Duales, Puntos_JK);
            push!(Informacion_Duales, Informacion_JK);
        end
    end
    
    #Se utiliza la función flatten para que convierta un arreglo de arreglos en un sólo arreglo grande.
    return collect(Iterators.flatten(Puntos_Duales)), collect(Iterators.flatten(Informacion_Duales))
end
