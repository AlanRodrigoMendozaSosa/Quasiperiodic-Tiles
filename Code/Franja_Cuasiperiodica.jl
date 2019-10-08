#La función determina, fijando los vectores estrella Ej y Ek, y para algún conjunto de alfas (uno por cada vector 
#estrella), los puntos de la red dual dejando fijo el valor Nj y corriendo el valor enteros Nk desde -N hasta N
#"J" y "K" son los índices de los vectores estrella a considerar
#"Nj" es el número entero asociado a la recta ortogonal al vector Ej
#"N_Intervalo_K" es el número entero que determina el rango [-N, N] en el que se barren los enteros de las rectas ortogonales al vector Ek.
#"Vectores_Estrella" es el conjunto de vectores estrella
#"Arreglo_Alfas" es el conjunto con las constantes alfas asociadas a cada vector estrella
function puntos_Dual_JK_Franja(J, K, Nj, N_Intervalo_K, Vectores_Estrella, Arreglo_Alfas)
    #Definimos el arreglo que contendra los vectores asociados a cada punto
    Puntos_Red_Dual = [];
    Informacion_Puntos_Red_Dual = [];
    
    for Nk in -N_Intervalo_K:N_Intervalo_K
        try #Ponemos el Try-Catch para que no debamos separar los casos en que Ej y Ek son paralelos, el error lo maneja autom.
            t0, t1, t2, t3, info = cuatro_Regiones(J, K, Nj, Nk, Vectores_Estrella, Arreglo_Alfas)
            push!(Puntos_Red_Dual, t0);
            push!(Puntos_Red_Dual, t1);
            push!(Puntos_Red_Dual, t2);
            push!(Puntos_Red_Dual, t3);
            push!(Informacion_Puntos_Red_Dual, info);
        catch
            nothing
        end
    end
    
    return Puntos_Red_Dual, Informacion_Puntos_Red_Dual
end

#Funcion que genera los polígonos obtenidos al usar el vector J (fijo) y el vector K (variable).
#"J" es el índice del vector estrella que se mantendrá fijo
#"Nj" es el número entero asociado a la recta ortogonal al vector Ej
#"N_Intervalo_K" es el número entero que determina el rango [-N, N] en el que se barren los enteros de las rectas ortogonales al vector Ek.
#"Vectores_Estrella" es el conjunto de vectores estrella
#"Arreglo_Alfas" es el conjunto con las constantes alfas asociadas a cada vector estrella
function franjas_Cuasiperiodicas(J, Nj, N_Intervalo_K, Vectores_Estrella, Arreglo_Alfas)
    #Definimos el arreglo que contendrá todos los puntos de la red Dual
    Puntos_Duales = [];
    Informacion_Duales = [];
    
    #Corramos los índices K de los vectores Ek
    for K in 1:length(Vectores_Estrella)
        if K == J
            nothing;
        else
            Puntos_JK, Informacion_JK = puntos_Dual_JK_Franja(J, K, Nj, N_Intervalo_K, Vectores_Estrella, Arreglo_Alfas)
            push!(Puntos_Duales, Puntos_JK);
            push!(Informacion_Duales, Informacion_JK);
        end
    end
    
    #Se utiliza la función flatten para que convierta un arreglo de arreglos en un sólo arreglo grande.
    return collect(Iterators.flatten(Puntos_Duales)), collect(Iterators.flatten(Informacion_Duales))
end