#Función que determina si el Vértice dado está dentro de un círculo de radio Radio, la función regresa un true o false.
#"Vertice" es un arreglo dos dimensional con las coordenadas de un vértice de un polígono.
#"Centro" es un arreglo dos dimensional con las coordenadas del centro del círculo.
#"Radio" es un número flotante asociado al radio del círculo.
function dentro_Radio(Vertice, Centro, Radio)
    if ((Vertice[1]-Centro[1])^2 + (Vertice[2] - Centro[2])^2 <= Radio^2)
        return true
    else
        return false
    end
end

#MODIFICACIÓN QUE NO GUARDA LA INFORMACIÓN DE LOS POLÍGONOS QUE CAEN FUERA DE UN CÍRCULO.
#Función que genera los vértices de un arreglo cuasiperiódico asociados a la vecindad de un Punto arbitrario.
#"Proyecciones" es un arreglo de números asociados a las proyecciones del punto de interés con los vectores estrella.
#"Vectores_Estrella" es el conjunto de vectores estrella.
#"Arreglo_Alfas" es el conjunto con las constantes alfas asociadas a cada vector estrella.
#"N" es el margen de error numérico al determinar el número entero de la proyección del punto de interés con los vectores estrella.
#"Radio" es un número flotante asociado al radio del círculo.
#"Centro" es un arreglo dos dimensional con las coordenadas del centro del círculo.
function generador_Vecindades_Acotado(Proyecciones, Vectores_Estrella, Arreglo_Alfas, N, Radio, Centro)
    #Definimos el arreglo que contendra los vértices asociados a cada combinación
    Puntos_Red_Dual = [];
    
    for i in 1:length(Vectores_Estrella)
        for j in i+1:length(Vectores_Estrella)
            
            for n in -N:N
                for m in -N:N
                    #Vamos a dejar que el try ... catch se encargue de los vectores estrella paralelos
                    try
                        #Obtengamos los vértices del arreglo considerando los vectores Ei y Ej con sus respectivos números enteros
                        t0, t1, t2, t3 = cuatro_Regiones(i, j, round(Proyecciones[i])+n, round(Proyecciones[j])+m, Vectores_Estrella, Arreglo_Alfas);
                        
                        #|| dentro_Radio(t1, Centro, Radio) || dentro_Radio(t2, Centro, Radio) || dentro_Radio(t3, Centro, Radio)
                        if dentro_Radio(t0, Centro, Radio) || dentro_Radio(t1, Centro, Radio) || dentro_Radio(t2, Centro, Radio) || dentro_Radio(t3, Centro, Radio)
                            push!(Puntos_Red_Dual, t0);
                            push!(Puntos_Red_Dual, t1);
                            push!(Puntos_Red_Dual, t2);
                            push!(Puntos_Red_Dual, t3);
                        end
                    catch
                        nothing
                    end                    
                end
            end
            
        end
    end
    
    return Puntos_Red_Dual
end

#Función que genera los vértices del arreglo cuasiperiódico alrededor de un punto arbitrario (dado o no por el usuario). El arreglo cuasiperiódico tendrá únicamente
#polígonos con al menos un vértice dentro de un radio alrededor del punto de interés.
#"N" es el margen de error numérico al determinar el número entero de la proyección del punto de interés con los vectores estrella.
#"Radio" es un número flotante asociado al radio del círculo.
#"Promedios_Distancia" es el arreglo con la separación entre las franjas cuasiperiódicas.
#"Vectores_Estrella" es el conjunto de vectores estrella.
#"Arreglo_Alfas" es el conjunto con las constantes alfas asociadas a cada vector estrella.
#"Punto" es el punto alrededor del cual se genera la vecindad cuasiperiódica.
function region_Local_Radio(N, Radio, Promedios_Distancia, Vectores_Estrella, Arreglo_Alfas, Punto)
    #Obtenemos las proyecciones del punto deseado con cada uno de los vectores estrella para generar los posibles números enteros
    Proyecciones = proyecciones_Pto_Direccion_Franjas(Punto, Promedios_Distancia, Vectores_Estrella);
    #Con las proyecciones y el punto deseado como centro del círculo discriminante, generamos la vecindad cuasiperiódica
    Puntos_Duales = generador_Vecindades_Acotado(Proyecciones, Vectores_Estrella, Arreglo_Alfas, N, Radio, Punto);

    return Puntos_Duales
end

#Función que determina el polígono que contiene al punto de interés.
#"Coordenadas_X" es un arreglo con las coordenadas en X de los vértices de los polígonos. Cada 4 entradas corresponden a un mismo polígono.
#"Coordenadas_Y" es un arreglo con las coordenadas en Y de los vértices de los polígonos. Cada 4 entradas corresponden a un mismo polígono.
#"Punto" es un arreglo dos dimensional con las coordenadas de un punto en el espacio 2D.
function poligono_Contenedor_Radio(Coordenadas_X, Coordenadas_Y, Punto)
    #Generamos los segmentos de los polígonos de la retícula cuasiperiódica.
    Poligonos = obtener_Segmentos_Vertices(Coordenadas_X, Coordenadas_Y);
    #Obtenemos los vértices del polígono contenedor
    Vertices_Contenedor = encontrar_Poligono(Punto, Poligonos);
    
    return Vertices_Contenedor
end