#Función que determina si el Vértice dado está dentro de un círculo de radio Radio, la función regresa un true o false
function dentro_Radio(Vertice, Centro, Radio)
    if ((Vertice[1]-Centro[1])^2 + (Vertice[2] - Centro[2])^2 <= Radio^2)
        return true
    else
        return false
    end
end

#MODIFICACIÓN QUE NO GUARDA LA INFORMACIÓN DE LOS POLÍGONOS QUE CAEN FUERA DE UN CÍRCULO
#Función que genera los vértices de un arreglo cuasiperiódico asociados a la vecindad de un Punto arbitrario
#El parámetro N nos sirve para "corregir" el posible error numérico al determinar el número entero asociado a 
#cada vector
function generador_Vecindades_Acotado(Proyecciones, Vectores_Estrella, Arreglo_Alfas, N, Radio, Centro)
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
                    
                    #|| dentro_Radio(t1, Centro, Radio) || dentro_Radio(t2, Centro, Radio) || dentro_Radio(t3, Centro, Radio)
                    if dentro_Radio(t0, Centro, Radio) || dentro_Radio(t1, Centro, Radio) || dentro_Radio(t2, Centro, Radio) || dentro_Radio(t3, Centro, Radio)
                        push!(Puntos_Red_Dual, t0);
                        push!(Puntos_Red_Dual, t1);
                        push!(Puntos_Red_Dual, t2);
                        push!(Puntos_Red_Dual, t3);
                        #push!(Puntos_Red_Dual, t0);
                        push!(Informacion_Puntos_Red_Dual, info);
                    end
                    
                end
            end
            
        end
    end
    
    return Puntos_Red_Dual, Informacion_Puntos_Red_Dual
end

function region_Local_Radio(N, Cota, Radio, Punto = true)
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
    
    println(length(Coordenadas_X)/4)
    
    return Coordenadas_X, Coordenadas_Y, Punto, Informacion_Duales, Proyecciones
end

function poligono_Contenedor_Radio(Coordenadas_X, Coordenadas_Y, Punto, Informacion_Duales)
    
    Poligonos = obtener_Segmentos_Vertices(Coordenadas_X, Coordenadas_Y);
    
    Informacion_Poligono_Contenedor = encontrar_Poligono(Punto, Poligonos, Informacion_Duales);
    
    Punto1, Punto2, Punto3, Punto4, Info = cuatro_Regiones(Int(Informacion_Poligono_Contenedor[1]), Int(Informacion_Poligono_Contenedor[2]), Int(Informacion_Poligono_Contenedor[3]), Int(Informacion_Poligono_Contenedor[4]), Vectores_Estrella, Arreglo_Alfas);
    
    plot()    
    plot([Punto1[1], Punto2[1], Punto3[1], Punto4[1], Punto1[1]], [Punto1[2], Punto2[2], Punto3[2], Punto4[2], Punto1[2]], markersize = 0.2, key = false, aspect_ratio=:equal, grid = false, color =:black)
    scatter!([Punto[1]], [Punto[2]], markersize = 5, markeralpha = 0.5, markerstrokewidth = 0, markercolor = :red)
    
    return Informacion_Poligono_Contenedor, Punto
    
end
