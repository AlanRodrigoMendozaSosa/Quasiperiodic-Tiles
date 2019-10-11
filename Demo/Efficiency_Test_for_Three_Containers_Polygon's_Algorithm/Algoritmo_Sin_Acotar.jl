#Función que genera los vértices de un arreglo cuasiperiódico asociados a la vecindad de un Punto arbitrario (dado o no por el usuario)
#"N" es el margen de error numérico al determinar el número entero de la proyección del punto de interés con los vectores estrella
#"Cota" es el Semilado de la caja centrada en el origen dentro de la cual se va a calcular el punto arbitrario si el usuario no da uno
function region_Local(N, Cota, Promedios_Distancia, Vectores_Estrella, Arreglo_Alfas, Punto = true)
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
    
    Puntos_Duales, Informacion_Duales = generador_Vecindades_Vertices(Proyecciones, Vectores_Estrella, Arreglo_Alfas, N);
    
    Coordenadas_X, Coordenadas_Y = separacion_Arreglo_de_Arreglos_2D(Puntos_Duales);
    
    return Coordenadas_X, Coordenadas_Y, Punto, Informacion_Duales, Proyecciones
end

#Función que determina el polígono que contiene al punto de interés.
#"Coordenadas_X" es un arreglo con las coordenadas en X de los vértices de los polígonos. Cada 4 entradas corresponden a un mismo polígono.
#"Coordenadas_Y" es un arreglo con las coordenadas en Y de los vértices de los polígonos. Cada 4 entradas corresponden a un mismo polígono.
#"Punto" es un arreglo dos dimensional con las coordenadas de un punto en el espacio 2D.
#"Informacion_Duales" es un arreglo con la información con la que se generó cada polígono del arreglo cuasiperiódico
function poligono_Contenedor(Coordenadas_X, Coordenadas_Y, Punto, Informacion_Duales, Vectores_Estrella, Arreglo_Alfas)
    
    Poligonos = obtener_Segmentos_Vertices(Coordenadas_X, Coordenadas_Y);
    
    #return Poligonos
    
    Informacion_Poligono_Contenedor = encontrar_Poligono(Punto, Poligonos, Informacion_Duales);
    
    Punto1, Punto2, Punto3, Punto4, Info = cuatro_Regiones(Int(Informacion_Poligono_Contenedor[1]), Int(Informacion_Poligono_Contenedor[2]), Int(Informacion_Poligono_Contenedor[3]), Int(Informacion_Poligono_Contenedor[4]), Vectores_Estrella, Arreglo_Alfas);
    
    return Punto1, Punto2, Punto3, Punto4, Informacion_Poligono_Contenedor
    
end