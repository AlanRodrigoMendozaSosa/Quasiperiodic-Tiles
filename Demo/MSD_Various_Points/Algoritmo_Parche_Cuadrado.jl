#Misma función que "Centroides_a_Vertices" sólo que esta está pensada para que, en caso de que no encuentre en un diccionario al centroide en cuestión, el código no
#se rompa, si no que sólo ignore ese paso y continúe; eventualmente el centroide de interés estará en alguno de los cuatro diccionarios a revisar.
#"Centroides" es un arreglo con los centroides de interés.
#"Diccionario_Centroides" es un diccionario que relaciona los centroides con los vértices del polígono que lo generó.
function centroides_A_Vertices_2(Centroides, Diccionario_Centroides)
    Coordenadas_X = [];
    Coordenadas_Y = [];

    for i in Centroides
        try
            push!(Coordenadas_X, Diccionario_Centroides[i][1])
            push!(Coordenadas_Y, Diccionario_Centroides[i][2])
        catch
            nothing
        end
    end
    
    return collect(Iterators.flatten(Coordenadas_X)), collect(Iterators.flatten(Coordenadas_Y))
end

#Función para generar un parche cuadrado con vértices de los polígonos
#"Informacion_Parche" es el arreglo con todos los datos generados para conocer el radio de los parches circulares
#"Semilado" es el semilado del cuadrado centrado en el origen donde, en caso de no tener un punto dado, el programa generará uno aleatoriamente.
#"Factor_Reduccion" es el factor por el que se multiplica el radio promedio de los parches circulares, se toma generalmente como 0.8 y ha mostrado buenos resultados.
#"Punto" es el parámetro que en caso de ser cierto, implica que queremos que el programa nos genere un punto aleatorio
function parche_Cuadrado(Informacion_Parche, Semilado, Factor_Reduccion, Promedios_Distancia, Vectores_Estrella, Arreglo_Alfas, Punto = true)
    #Si el usuario no proporciona un punto, generamos uno
    if Punto == true
        x = rand();
        y = rand();

        if (x > 0.5) && (y > 0.5)
            Punto = [rand()*Semilado, rand()*Semilado];
        elseif (x > 0.5) && (y < 0.5)
            Punto = [rand()*Semilado, -rand()*Semilado];
        elseif (x < 0.5) && (y > 0.5)
            Punto = [-rand()*Semilado, rand()*Semilado];
        elseif (x < 0.5) && (y < 0.5)
            Punto = [-rand()*Semilado, -rand()*Semilado];
        end
    end
    
    #Paso 1: Generemos los cuatro centros para los parches circulares.
    P_1 = Punto + [(sqrt(2)/2)*(Factor_Reduccion)*(Informacion_Parche[6]), (sqrt(2)/2)*(Factor_Reduccion)*(Informacion_Parche[6])];
    P_2 = Punto + [-(sqrt(2)/2)*(Factor_Reduccion)*(Informacion_Parche[6]), (sqrt(2)/2)*(Factor_Reduccion)*(Informacion_Parche[6])];
    P_3 = Punto - [(sqrt(2)/2)*(Factor_Reduccion)*(Informacion_Parche[6]), (sqrt(2)/2)*(Factor_Reduccion)*(Informacion_Parche[6])];
    P_4 = Punto + [(sqrt(2)/2)*(Factor_Reduccion)*(Informacion_Parche[6]), -(sqrt(2)/2)*(Factor_Reduccion)*(Informacion_Parche[6])];
    
    #Paso 2: Generamos las vecindades de los arreglos cuasiperiódicos centrados en esos puntos dados.
    #Generación del arreglo cuasiperiódico con margen de error (1era entrada), la segunda entrada es irrelevante
    #dado que estamos dado el punto en el cual centraremos el arreglo cuasiperiódico; el punto es la tercera entrada
    Coordenadas_X1, Coordenadas_Y1, P_1 = region_Local_Voronoi(Informacion_Parche[1], Informacion_Parche[2], Promedios_Distancia, Vectores_Estrella, Arreglo_Alfas, P_1);
    Coordenadas_X2, Coordenadas_Y2, P_2 = region_Local_Voronoi(Informacion_Parche[1], Informacion_Parche[2], Promedios_Distancia, Vectores_Estrella, Arreglo_Alfas, P_2);
    Coordenadas_X3, Coordenadas_Y3, P_3 = region_Local_Voronoi(Informacion_Parche[1], Informacion_Parche[2], Promedios_Distancia, Vectores_Estrella, Arreglo_Alfas, P_3);
    Coordenadas_X4, Coordenadas_Y4, P_4 = region_Local_Voronoi(Informacion_Parche[1], Informacion_Parche[2], Promedios_Distancia, Vectores_Estrella, Arreglo_Alfas, P_4);
    
    #Paso 3: Obtenemos los centroides asociados a estos polígonos del arreglo cuasiperiódico.
    #Obtención de los centroides de los poligonos del arreglo cuasiperiódico en: 
    #coordenada X, coordenada Y, Float64(X,Y) y su diccionario Centroides <-> Vertices.
    Centroides_X1, Centroides_Y1, Centroides1, Diccionario_Centroides1 = centroides(Coordenadas_X1, Coordenadas_Y1);
    Centroides_X2, Centroides_Y2, Centroides2, Diccionario_Centroides2 = centroides(Coordenadas_X2, Coordenadas_Y2);
    Centroides_X3, Centroides_Y3, Centroides3, Diccionario_Centroides3 = centroides(Coordenadas_X3, Coordenadas_Y3);
    Centroides_X4, Centroides_Y4, Centroides4, Diccionario_Centroides4 = centroides(Coordenadas_X4, Coordenadas_Y4);
    
    #Paso 4: Depuramos los arreglos cuasiperiódicos obtenidos para quedarnos con el clúster principal.
    #Obtención de los centroides de los buenos polígonos tras iterar N veces (3era entrada) algoritmo de áreas.
    Centroides_Iterados1 = centroides_Area_Acotada_Iterada(Centroides1, Informacion_Parche[4], Informacion_Parche[3]);
    Centroides_Iterados2 = centroides_Area_Acotada_Iterada(Centroides2, Informacion_Parche[4], Informacion_Parche[3]);
    Centroides_Iterados3 = centroides_Area_Acotada_Iterada(Centroides3, Informacion_Parche[4], Informacion_Parche[3]);
    Centroides_Iterados4 = centroides_Area_Acotada_Iterada(Centroides4, Informacion_Parche[4], Informacion_Parche[3]);
    
    #Paso 5: Agrupamos todos los buenos centroides en un único arreglo, quedándonos únicamente con los centroides sin
    #repetir.
    Arreglo_Centroides = vcat(Centroides_Iterados1, Centroides_Iterados2, Centroides_Iterados3, Centroides_Iterados4);
    #println("There is a total of $(length(Arreglo_Centroides)) polygons with possible duplicates.")
    unique!(Arreglo_Centroides);
    #println("There is a total of $(length(Arreglo_Centroides)) uniques polygons.")
    
    #Paso5.5: Generamos un único diccionario al unir los 4 diccionarios generados previamente
    Diccionario_Centroides_Unico = merge(Diccionario_Centroides1, Diccionario_Centroides2, Diccionario_Centroides3, Diccionario_Centroides4)

    #Paso 6: Recuperamos los polígonos del arrego cuasiperiódico que corresponden a los clúster principales.
    #Visualizamos ahora el arreglo cuasiperiódico tras quitarle las n capas.
    Coordenadas_X_CP, Coordenadas_Y_CP = centroides_A_Vertices(Arreglo_Centroides, Diccionario_Centroides_Unico);
    
    return Coordenadas_X_CP, Coordenadas_Y_CP, Punto, Arreglo_Centroides, Diccionario_Centroides_Unico
end
