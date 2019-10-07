#Función que dado un Punto, obtiene aproximadamente los enteros que forman el polígono que lo contiene
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
#El parámetro N nos sirve para "corregir" el posible error numérico al determinar el número entero asociado a 
#cada vector
function generador_Vecindades_Vertices(Proyecciones, Vectores_Estrella, Arreglo_Alfas, N)
    #Definimos el arreglo que contendra los vértices asociados a cada combinación
    Puntos_Red_Dual = [];
    Informacion_Puntos_Red_Dual = [];
    
    for i in 1:length(Vectores_Estrella)
        for j in i+1:length(Vectores_Estrella)
            
            for n in -N:N
                for m in -N:N
                    #println(i,j)
                    #Obtengamos los vértices del arreglo considerando los vectores Ei y Ej con sus respectivos números enteros
                    t0, t1, t2, t3, info = cuatro_Regiones(i, j, round(Proyecciones[i])+n, round(Proyecciones[j])+m, Vectores_Estrella, Arreglo_Alfas);
                    push!(Puntos_Red_Dual, t0);
                    push!(Puntos_Red_Dual, t1);
                    push!(Puntos_Red_Dual, t2);
                    push!(Puntos_Red_Dual, t3);
                    push!(Informacion_Puntos_Red_Dual, info);
                end
            end
            
        end
    end
    
    return Puntos_Red_Dual, Informacion_Puntos_Red_Dual
end