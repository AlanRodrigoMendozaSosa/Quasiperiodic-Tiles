function region_Local(N, Cota, Punto = true)
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
    
    println(length(Coordenadas_X)/4)
    
    return Coordenadas_X, Coordenadas_Y, Punto, Informacion_Duales, Proyecciones
end

function poligono_Contenedor(Coordenadas_X, Coordenadas_Y, Punto, Informacion_Duales)
    
    Poligonos = obtener_Segmentos_Vertices(Coordenadas_X, Coordenadas_Y);
    
    #return Poligonos
    
    Informacion_Poligono_Contenedor = encontrar_Poligono(Punto, Poligonos, Informacion_Duales);
    
    Punto1, Punto2, Punto3, Punto4, Info = cuatro_Regiones(Int(Informacion_Poligono_Contenedor[1]), Int(Informacion_Poligono_Contenedor[2]), Int(Informacion_Poligono_Contenedor[3]), Int(Informacion_Poligono_Contenedor[4]), Vectores_Estrella, Arreglo_Alfas);
    
    plot()    
    plot([Punto1[1], Punto2[1], Punto3[1], Punto4[1], Punto1[1]], [Punto1[2], Punto2[2], Punto3[2], Punto4[2], Punto1[2]], markersize = 0.2, key = false, aspect_ratio=:equal, grid = false, color =:black)
    scatter!([Punto[1]], [Punto[2]], markersize = 5, markeralpha = 0.5, markerstrokewidth = 0, markercolor = :red)
    
    return Informacion_Poligono_Contenedor, Punto
    
end
