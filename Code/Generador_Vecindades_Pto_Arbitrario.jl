#Función que dado un Punto, obtiene aproximadamente los enteros que forman el polígono que lo contiene
#"Punto" es un arreglo dos dimensional con las coordenadas de un punto en el espacio 2D
#"Arreglo_Distancias_Promedio_Rectas_V" es un arreglo con la separación promedio entre rectas ortogonales a cada vector estrella
#"Vectores_Estrella" es el conjunto de vectores estrella
function proyecciones_Pto_Direccion_Franjas(Punto, Arreglo_Distancias_Promedio_Rectas_V, Vectores_Estrella)

    #Obtengamos ahora la proyección del punto con la línea de cada una de las franjas generadas por cada vector
    Arreglo_Proyeccion_Pto_Direccion_Franjas = [];
    
    for i in 1:length(Vectores_Estrella)
        Proyeccion = prod_Punto(Punto, Vectores_Estrella[i])/Arreglo_Distancias_Promedio_Rectas_V[i];
        #Al parecer no hay que restarle "- Arreglo_Distancia_Origen_Primera_Recta[i]" que es la DOPR
        push!(Arreglo_Proyeccion_Pto_Direccion_Franjas, Proyeccion);
    end
    
    return Arreglo_Proyeccion_Pto_Direccion_Franjas
end

#Función que genera los vértices de un arreglo cuasiperiódico asociados a la vecindad de un Punto arbitrario
#"Proyecciones" es el conjunto de números enteros candidatos a ser los que generan las rectas ortogonales a los vectores Ej y Ek que contienen al punto
#"Vectores_Estrella" es el conjunto de vectores estrella
#"Arreglo_Alfas" es el conjunto con las constantes alfas asociadas a cada vector estrella
#El parámetro "N" nos sirve para "corregir" el posible error numérico al determinar el número entero asociado a cada vector
function generador_Vecindades_Vertices(Proyecciones, Vectores_Estrella, Arreglo_Alfas, N)
    #Definimos el arreglo que contendra los vértices asociados a cada combinación
    Puntos_Red_Dual = [];
    Informacion_Puntos_Red_Dual = [];
    
    for i in 1:length(Vectores_Estrella)
        for j in i+1:length(Vectores_Estrella)
            
            for n in -N:N
                for m in -N:N
                    try #Vamos a dejar que el try ... catch se encargue de los casos en que los vectores estrella sean paralelos
                        #println(i,j)
                        #Obtengamos los vértices del arreglo considerando los vectores Ei y Ej con sus respectivos números enteros
                        t0, t1, t2, t3, info = cuatro_Regiones(i, j, round(Proyecciones[i])+n, round(Proyecciones[j])+m, Vectores_Estrella, Arreglo_Alfas);
                        push!(Puntos_Red_Dual, t0);
                        push!(Puntos_Red_Dual, t1);
                        push!(Puntos_Red_Dual, t2);
                        push!(Puntos_Red_Dual, t3);
                        push!(Informacion_Puntos_Red_Dual, info);
                    catch
                        nothing;
                    end
                end
            end
            
        end
    end
    
    return Puntos_Red_Dual, Informacion_Puntos_Red_Dual
end