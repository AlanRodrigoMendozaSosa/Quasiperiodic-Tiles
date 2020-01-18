#Función que determina si el Vértice dado está dentro de un círculo de radio Radio, la función regresa un true o false
#"Vertice" es un arreglo dos dimensional con las coordenadas de un vértice de un polígono.
#"Centro" es un arreglo dos dimensional con las coordenadas del centro del círculo
#"Radio" es un número flotante asociado al radio del círculo.
function dentro_Radio(Vertice, Centro, Radio)
    if ((Vertice[1]-Centro[1])^2 + (Vertice[2] - Centro[2])^2 <= Radio^2)
        return true
    else
        return false
    end
end

#MODIFICACIÓN QUE NO GUARDA LA INFORMACIÓN DE LOS POLÍGONOS QUE CAEN FUERA DE UN CÍRCULO
#Función que genera los vértices de un arreglo cuasiperiódico asociados a la vecindad de un Punto arbitrario
#"Proyecciones" es un arreglo de números asociados a las proyecciones del punto de interés con los vectores estrella
#"Vectores_Estrella" es el conjunto de vectores estrella
#"Arreglo_Alfas" es el conjunto con las constantes alfas asociadas a cada vector estrella
#"N" es el margen de error numérico al determinar el número entero de la proyección del punto de interés con los vectores estrella
#"Radio" es un número flotante asociado al radio del círculo
#"Centro" es un arreglo dos dimensional con las coordenadas del centro del círculo
function generador_Vecindades_Acotado(Proyecciones, Vectores_Estrella, Arreglo_Alfas, N, Radio, Centro)
    #Definimos el arreglo que contendra los vértices asociados a cada combinación
    Puntos_Red_Dual = [];
    Informacion_Puntos_Red_Dual = [];
    
    for i in 1:length(Vectores_Estrella)
        for j in i+1:length(Vectores_Estrella)
            
            for n in -N:N
                for m in -N:N
                    try #Vamos a dejar que el try ... catch se encargue de los vectores estrella paralelos
                        #println(i,j)
                        #Obtengamos los vértices del arreglo considerando los vectores Ei y Ej con sus respectivos números enteros
                        t0, t1, t2, t3, info = cuatro_Regiones(i, j, round(Proyecciones[i])+n, round(Proyecciones[j])+m, Vectores_Estrella, Arreglo_Alfas);
                        
                        #|| dentro_Radio(t1, Centro, Radio) || dentro_Radio(t2, Centro, Radio) || dentro_Radio(t3, Centro, Radio)
                        if dentro_Radio(t0, Centro, Radio) || dentro_Radio(t1, Centro, Radio) || dentro_Radio(t2, Centro, Radio) || dentro_Radio(t3, Centro, Radio)
                            push!(Puntos_Red_Dual, t0);
                            push!(Puntos_Red_Dual, t1);
                            push!(Puntos_Red_Dual, t2);
                            push!(Puntos_Red_Dual, t3);
                            #push!(Puntos_Red_Dual, t0);
                            push!(Informacion_Puntos_Red_Dual, info);
                        end
                    catch
                        nothing
                    end                    
                end
            end
            
        end
    end
    
    return Puntos_Red_Dual, Informacion_Puntos_Red_Dual
end

#Función que genera los vértices del arreglo cuasiperiódico alrededor de un punto arbitrario (dado o no por el usuario). El arreglo cuasiperiódico tendrá únicamente
#polígonos con al menos un vértice dentro de un radio alrededor del punto de interés.
#"N" es el margen de error numérico al determinar el número entero de la proyección del punto de interés con los vectores estrella
#"Cota" es el Semilado de la caja centrada en el origen dentro de la cual se va a calcular el punto arbitrario si el usuario no da uno
#"Radio" es un número flotante asociado al radio del círculo
function region_Local_Radio(N, Cota, Radio, Promedios_Distancia, Vectores_Estrella, Arreglo_Alfas, Punto = true)
    #Si el usuario no proporciona un punto, generamos uno
    if Punto == true
        x = rand();
        y = rand();

        if (x > 0.5) && (y > 0.5)
            Punto = [rand()*Cota, rand()*Cota];
        elseif (x > 0.5) && (y < 0.5)
            Punto = [rand()*Cota, -rand()*Cota];
        elseif (x < 0.5) && (y > 0.5)
            Punto = [-rand()*Cota, rand()*Cota];
        elseif (x < 0.5) && (y < 0.5)
            Punto = [-rand()*Cota, -rand()*Cota];
        end
    end
    
    Proyecciones = proyecciones_Pto_Direccion_Franjas(Punto, Promedios_Distancia, Vectores_Estrella);
    
    Puntos_Duales, Informacion_Duales = generador_Vecindades_Acotado(Proyecciones, Vectores_Estrella, Arreglo_Alfas, N, Radio, Punto);
    
    Coordenadas_X, Coordenadas_Y = separacion_Arreglo_de_Arreglos_2D(Puntos_Duales);
    
    return Coordenadas_X, Coordenadas_Y, Punto, Informacion_Duales, Proyecciones
end

#Función que determina el polígono que contiene al punto de interés.
#"Coordenadas_X" es un arreglo con las coordenadas en X de los vértices de los polígonos. Cada 4 entradas corresponden a un mismo polígono.
#"Coordenadas_Y" es un arreglo con las coordenadas en Y de los vértices de los polígonos. Cada 4 entradas corresponden a un mismo polígono.
#"Punto" es un arreglo dos dimensional con las coordenadas de un punto en el espacio 2D.
#"Informacion_Duales" es un arreglo con la información con la que se generó cada polígono del arreglo cuasiperiódico
function poligono_Contenedor_Radio(Coordenadas_X, Coordenadas_Y, Punto, Informacion_Duales, Vectores_Estrella, Arreglo_Alfas)
    
    Poligonos = obtener_Segmentos_Vertices(Coordenadas_X, Coordenadas_Y);
    
    Informacion_Poligono_Contenedor = encontrar_Poligono(Punto, Poligonos, Informacion_Duales);
    
    Punto1, Punto2, Punto3, Punto4, Info = cuatro_Regiones(Int(Informacion_Poligono_Contenedor[1]), Int(Informacion_Poligono_Contenedor[2]), Int(Informacion_Poligono_Contenedor[3]), Int(Informacion_Poligono_Contenedor[4]), Vectores_Estrella, Arreglo_Alfas);
    
    return Punto1, Punto2, Punto3, Punto4, Informacion_Poligono_Contenedor
    
end