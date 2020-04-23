#Función que determina, fijando los vectores estrella Ej y Ek, y para algún conjunto de constante alfa, los puntos de la retícula cuasiperiódica dejando fijo el valor Nj y corriendo 
#el valor entero Nk desde -N hasta N.
#"J" y "K" son los índices de los vectores estrella a considerar.
#"Nj" es el número entero asociado a la recta ortogonal al vector Ej.
#"N_Intervalo_K" es el número entero que determina el rango [-N, N] en el que se barren los enteros de las rectas ortogonales al vector Ek.
#"Vectores_Estrella" es el conjunto de vectores estrella.
#"Arreglo_Alfas" es el arreglo con los valores numéricos de la separación respecto al origen del conjunto de rectas ortogonales a los vectores estrella en el método generalizado dual.
function puntos_Dual_JK_Franja(J, K, Nj, N_Intervalo_K, Vectores_Estrella, Arreglo_Alfas)
    #Paso 1: Definimos el arreglo que contendrá las coordenadas de los vértices de la retícula cuasiperiódica.
    Puntos_Red_Dual = [];
    
    #Paso 2: Iteramos sobre los distintos valores que puede tomar el número entero Nk.
    for Nk in -N_Intervalo_K:N_Intervalo_K
        #Permitimos que el Try---Catch maneje el error que surje cuando los vectores son colineales.
        try
            #Salida: [X,Y]
            t0, t1, t2, t3 = cuatro_Regiones(J, K, Nj, Nk, Vectores_Estrella, Arreglo_Alfas);
            push!(Puntos_Red_Dual, t0);
            push!(Puntos_Red_Dual, t1);
            push!(Puntos_Red_Dual, t2);
            push!(Puntos_Red_Dual, t3);
        catch
            nothing
        end
    end
    
    return Puntos_Red_Dual
end

#Funcion que genera los polígonos obtenidos al usar el vector J (fijo) y el vector K (variable).
#"J" es el índice del vector estrella que se mantendrá fijo.
#"Nj" es el número entero asociado a la recta ortogonal al vector Ej.
#"N_Intervalo_K" es el número entero que determina el rango [-N, N] en el que se barren los enteros de las rectas ortogonales al vector Ek.
#"Vectores_Estrella" es el conjunto de vectores estrella.
#"Arreglo_Alfas" es el arreglo con los valores numéricos de la separación respecto al origen del conjunto de rectas ortogonales a los vectores estrella en el método generalizado dual.
function franjas_Cuasiperiodicas(J, Nj, N_Intervalo_K, Vectores_Estrella, Arreglo_Alfas)
    #Paso 1: Definimos el arreglo que contendrá las coordenadas de los vértices de la retícula cuasiperiódica.
    Puntos_Duales = [];
    
    #Paso 2: Iteramos sobre los diferentes vectores Ek
    for K in 1:length(Vectores_Estrella)
        if K == J
            nothing;
        else
            #Salida: [[X,Y]]
            Puntos_JK = puntos_Dual_JK_Franja(J, K, Nj, N_Intervalo_K, Vectores_Estrella, Arreglo_Alfas);
            push!(Puntos_Duales, Puntos_JK);
        end
    end
    
    #Se utiliza la función flatten para que convierta un arreglo de arreglos en un sólo arreglo grande.
    return collect(Iterators.flatten(Puntos_Duales))
end