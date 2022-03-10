#Función que dado un Punto, obtiene aproximadamente los enteros que forman el polígono que lo contiene.
#"Punto" es un arreglo [X,Y] con las coordenadas de un punto en el espacio 2D.
#"Promedios_Distancia" es el arreglo con la separación entre las franjas cuasiperiódicas.
#"Vectores_Estrella" es el arreglo con los vectores estrella que generan la retícula cuasiperiódica deseada.
function proyecciones_Pto_Direccion_Franjas(Punto, Promedios_Distancia, Vectores_Estrella)
    #Paso 1: Generemos un arreglo en donde irán los números reales resultado de proyectar el Punto con los vectores estrella.
    Arreglo_Proyeccion_Pto_Direccion_Franjas = [];
    
    #Paso 2: Para cada vector estrella, proyectamos el Punto sobre dicho vector y reescalamos con la separación entre las franjas cuasiperiódicas.
    for i in 1:length(Vectores_Estrella)
        Proyeccion = prod_Punto(Punto, Vectores_Estrella[i]/norm(Vectores_Estrella[i]))/Promedios_Distancia[i];
        #Al parecer no hay que restarle "- Arreglo_Distancia_Origen_Primera_Recta[i]" que es la DOPR
        push!(Arreglo_Proyeccion_Pto_Direccion_Franjas, Proyeccion);
    end
    
    return Arreglo_Proyeccion_Pto_Direccion_Franjas
end

#Función que genera los vértices de un arreglo cuasiperiódico asociados a la vecindad de un Punto arbitrario.
#"Proyecciones" es el conjunto de números enteros candidatos a ser los que generan el polígono que contiene al punto.
#"Vectores_Estrella" es el conjunto de vectores estrella.
#"Arreglo_Alfas" es el arreglo con los valores numéricos de la separación respecto al origen del conjunto de rectas ortogonales a los vectores estrella en el método generalizado dual.
#"N" es el margen de error asociado a los números enteros generados por la proyección del punto sobre los vectores estrella.
function generador_Vecindades_Vertices(Proyecciones, Vectores_Estrella, Arreglo_Alfas, N)
    #Paso 1: Definimos el arreglo que contendrá a los vértices asociados a cada combinación de vectores estrella (incluídos los generados al considerar el margen de error)
    Puntos_Red_Dual = [];
    
    #Paso 2: Consideramos todas las posibles combinaciones de vectores estrella con los posibles números enteros correspondientes (y su margen de error)
    for i in 1:length(Vectores_Estrella)
        for j in i+1:length(Vectores_Estrella)
            #Consideramos el margen de error a cada número entero
            for n in -N:N
                for m in -N:N
                    #Vamos a dejar que el try ... catch se encargue de los casos en que los vectores estrella sean paralelos
                    try
                        #Obtengamos los vértices del arreglo considerando los vectores Ei y Ej con sus respectivos números enteros (y su margen de error)
                        #Salida: [X,Y]
                        t0, t1, t2, t3 = cuatro_Regiones(i, j, round(Proyecciones[i])+n, round(Proyecciones[j])+m, Vectores_Estrella, Arreglo_Alfas);
                        push!(Puntos_Red_Dual, t0);
                        push!(Puntos_Red_Dual, t1);
                        push!(Puntos_Red_Dual, t2);
                        push!(Puntos_Red_Dual, t3);
                    catch
                        nothing;
                    end
                end
            end
            
        end
    end
    
    return Puntos_Red_Dual
end