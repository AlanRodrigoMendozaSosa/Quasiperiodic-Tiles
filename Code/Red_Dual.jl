#Función que determina, fijando los vectores estrella Ej y Ek, fijando los enteros Nj y Nk, y para algún conjunto de constantes alfa, los cuatros puntos de la red dual.
#"J" y "K" son los índices de los vectores estrella a considerar.
#"Nj" y "Nk" son los enteros con los que se generan las rectas ortogonales a Ej y Ek.
#"Vectores_Estrella" es el conjunto de vectores estrella.
#"Arreglo_Alfas" es el arreglo con los valores numéricos de la separación respecto al origen del conjunto de rectas ortogonales a los vectores estrella en el método generalizado dual.
function cuatro_Regiones(J, K, Nj, Nk, Vectores_Estrella, Arreglo_Alfas)
    #Paso 0: Verifiquemos si los vectores a considerar son colineales, en cuyo caso manda un error.
    if (length(Vectores_Estrella) % 2 == 0) && (K == J + length(Vectores_Estrella)/2)
        error("Los vectores Ej y Ek no pueden ser paralelos")
    else
        #Paso 1: Definimos los dos vectores con los que se consigue la intersección en la malla generada por los vectores estrella que estamos considerando.
        Ej = Vectores_Estrella[J];
        Ek = Vectores_Estrella[K];
        
        #Paso 2: Obtenemos los vectores ortogonales a estos dos vectores.
        EjOrt = vector_Ortogonal(Ej);
        EkOrt = vector_Ortogonal(Ek);
        
        #Paso 3: Definimos los valores reales con los que se crearon las rectas ortogonales a cada vector Ej y Ek para tomar la intersección.
        Factor_Ej = Nj + Arreglo_Alfas[J];
        Factor_Ek = Nk + Arreglo_Alfas[K];
        
        #Paso 4: Obtenemos el área que forman los dos vectores Ej y Ek.
        AreaJK = Ej[1]*Ek[2] - Ej[2]*Ek[1];
        
        #Paso 5: Definimos lo que será el vértice en el espacio real de la retícula cuasiperiódica. Este vértice se denomina t^{0} en la teoría.
        Punto_Cero = Nj*Ej + Nk*Ek;
        
        #Paso 6: Generamos los términos asociados a la proyección del vector Ej y Ek con los demás vectores estrella
        for i in 1:length(Vectores_Estrella)
            if i == J || i == K
                nothing
            else
                Factor_Ei = (Factor_Ej/AreaJK)*(prod_Punto(EkOrt, Vectores_Estrella[i])) - (Factor_Ek/AreaJK)*(prod_Punto(EjOrt, Vectores_Estrella[i]));
                Punto_Cero += (floor(Factor_Ei - Arreglo_Alfas[i]))*Vectores_Estrella[i];
            end
        end
        
        #Paso 7: Obtenemos los otros tres vértices asociados al punto t^{0}.
        Punto_Uno = Punto_Cero - Ej;
        Punto_Dos = Punto_Cero - Ej - Ek;
        Punto_Tres = Punto_Cero - Ek;
        
        #Agregamos la información correspondiente a cómo se obtuvieron estos vértices
        #POSIBLEMENTE UTIL EN UN FUTURO, NO BORRAR EN FUTURAS REVISIONES DE SOFTWARE
        #Info = [J,K,Nj,Nk];
        
        return Punto_Cero, Punto_Uno, Punto_Dos, Punto_Tres
    end
end

#Función que determina, fijando los vectores estrella Ej y Ek, y para algún conjunto de constantes alfa, los puntos de la retícula en el espacio real corriendo los valores enteros 
#Nj y Nk desde -N hasta N.
#"J" y "K" son los índices de los vectores estrella a considerar.
#"N" es el número entero que determina el rango [-N, N] en el que se barren los enteros de las rectas ortogonales a los vectores Ej y Ek.
#"Vectores_Estrella" es el conjunto de vectores estrella.
#"Arreglo_Alfas" es el arreglo con los valores numéricos de la separación respecto al origen del conjunto de rectas ortogonales a los vectores estrella en el método generalizado dual.
function puntos_Dual_JK(J, K, N, Vectores_Estrella, Arreglo_Alfas)
    #Paso 1: Definimos el arreglo que contendra las coordenadas [X,Y] de cada vértice de la retícula cuasiperiódica.
    Puntos_Red_Dual = [];
    
    #Paso 2: Barrermos el conjunto de números enteros Nj y Nk para formar las intersecciones de las rectas ortogonales en el mallado.
    for Nj in -N:N
        for Nk in -N:N
            #Dejamos que el try---catch detecte el error que surje cuando los vectores considerados son colineales.
            try
                #Salida: [X,Y]
                t0, t1, t2, t3 = cuatro_Regiones(J, K, Nj, Nk, Vectores_Estrella, Arreglo_Alfas)
                push!(Puntos_Red_Dual, t0);
                push!(Puntos_Red_Dual, t1);
                push!(Puntos_Red_Dual, t2);
                push!(Puntos_Red_Dual, t3);
            catch
                nothing
            end
        end
    end
    
    return Puntos_Red_Dual
end

#Función que genera una vecindad alrededor del origen de una retícula cuasiperiódica por el método generalizado dual.
#"N" es el número entero que determina el rango [-N, N] en el que se barren los enteros de las rectas ortogonales a los vectores Ej y Ek.
#"Vectores_Estrella" es el conjunto de vectores estrella.
#"Arreglo_Alfas" es el arreglo con los valores numéricos de la separación respecto al origen del conjunto de rectas ortogonales a los vectores estrella en el método generalizado dual.
function puntos_Dual(N, Vectores_Estrella, Arreglo_Alfas)
    #Paso 1: Definimos el arreglo que contendrá todos los vértices de la retícula cuasiperiódica generada.
    Puntos_Duales = [];
    
    #Paso 2: Consideramos todos las posibles parejas de vectores Ej y Ek.
    for J in 1:length(Vectores_Estrella)
        for K in (J+1):length(Vectores_Estrella)
            #Salida: [[X,Y]]
            Puntos_JK = puntos_Dual_JK(J, K, N, Vectores_Estrella, Arreglo_Alfas)
            push!(Puntos_Duales, Puntos_JK);
        end
    end
    
    #Se utiliza la función flatten para que convierta un arreglo de arreglos en un sólo arreglo grande.
    return collect(Iterators.flatten(Puntos_Duales))
end