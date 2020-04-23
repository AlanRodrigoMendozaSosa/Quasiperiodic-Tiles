#Función para generar un parche cuadrado con vértices de las celdas unitarias del arreglo cuasiperiódico.
#"Informacion_Parche" es el arreglo con todos los datos asociados a los parches circulares.
#"Factor_Reduccion" es el factor por el que se multiplica el radio promedio de los parches circulares para generar el radio seguro de los parches.
#"Promedios_Distancia" es el arreglo con la separación entre las franjas cuasiperiódicas.
#"Vectores_Estrella" es el arreglo con los vectores estrella que generan la retícula cuasiperiódica deseada.
#"Arreglo_Alfas" es el arreglo con los valores numéricos de la separación respecto al origen del conjunto de rectas ortogonales a los vectores estrella en el método generalizado dual.
#"Punto" es el punto alrededor de donde se va a generar la vecindad.
function parche_Cuadrado(Informacion_Parche, Factor_Reduccion, Promedios_Distancia, Vectores_Estrella, Arreglo_Alfas, Punto)    
    #Paso 1: Generemos los cuatro centros para los parches circulares.
    #Informacion_Parche[6] es el radio promedio del parche circular.
    #Salida: [X,Y]
    P_1 = Punto + [(sqrt(2)/2)*(Factor_Reduccion)*(Informacion_Parche[6]), (sqrt(2)/2)*(Factor_Reduccion)*(Informacion_Parche[6])];
    P_2 = Punto + [-(sqrt(2)/2)*(Factor_Reduccion)*(Informacion_Parche[6]), (sqrt(2)/2)*(Factor_Reduccion)*(Informacion_Parche[6])];
    P_3 = Punto - [(sqrt(2)/2)*(Factor_Reduccion)*(Informacion_Parche[6]), (sqrt(2)/2)*(Factor_Reduccion)*(Informacion_Parche[6])];
    P_4 = Punto + [(sqrt(2)/2)*(Factor_Reduccion)*(Informacion_Parche[6]), -(sqrt(2)/2)*(Factor_Reduccion)*(Informacion_Parche[6])];
    
    #Paso 2: Generamos las vecindades de los arreglos cuasiperiódicos centrados en esos puntos dados.
    #Informacion_Parche[1] es el margen de error a los números enteros "n_j" empleado al generar los parches circulares.
    #Salida: [[X,Y]]
    Coordenadas_1 = region_Local_Voronoi(Informacion_Parche[1], Promedios_Distancia, Vectores_Estrella, Arreglo_Alfas, P_1);
    Coordenadas_2 = region_Local_Voronoi(Informacion_Parche[1], Promedios_Distancia, Vectores_Estrella, Arreglo_Alfas, P_2);
    Coordenadas_3 = region_Local_Voronoi(Informacion_Parche[1], Promedios_Distancia, Vectores_Estrella, Arreglo_Alfas, P_3);
    Coordenadas_4 = region_Local_Voronoi(Informacion_Parche[1], Promedios_Distancia, Vectores_Estrella, Arreglo_Alfas, P_4);
    
    #Paso 3: Obtenemos los centroides asociados a estos polígonos del arreglo cuasiperiódico.
    #Obtención de los centroides de los poligonos del arreglo cuasiperiódico en: 
    #Salida: [Float64(X,Y)] y diccionario Centroides <-> Vertices.
    Centroides1, Diccionario_Centroides1 = centroides(Coordenadas_1);
    Centroides2, Diccionario_Centroides2 = centroides(Coordenadas_2);
    Centroides3, Diccionario_Centroides3 = centroides(Coordenadas_3);
    Centroides4, Diccionario_Centroides4 = centroides(Coordenadas_4);
    
    #Paso 4: Depuramos los arreglos cuasiperiódicos obtenidos para quedarnos con el clúster principal.
    #Informacion_Parche[4] es la cota del área empleada como discriminante en el algoritmo de Voronoi.
    #Informacion_Parche[3] es el número de iteraciones del algoritmo de Voronoi para obtener el cluster principal.
    #Salida: [Float64(X,Y)]
    Centroides_Iterados1 = centroides_Area_Acotada_Iterada(Centroides1, Informacion_Parche[4], Informacion_Parche[3]);
    Centroides_Iterados2 = centroides_Area_Acotada_Iterada(Centroides2, Informacion_Parche[4], Informacion_Parche[3]);
    Centroides_Iterados3 = centroides_Area_Acotada_Iterada(Centroides3, Informacion_Parche[4], Informacion_Parche[3]);
    Centroides_Iterados4 = centroides_Area_Acotada_Iterada(Centroides4, Informacion_Parche[4], Informacion_Parche[3]);
    
    #Paso 5: Agrupamos todos los buenos centroides en un único arreglo, quedándonos únicamente con los centroides sin repetir.
    Arreglo_Centroides = vcat(Centroides_Iterados1, Centroides_Iterados2, Centroides_Iterados3, Centroides_Iterados4);
    unique!(Arreglo_Centroides);
    #println("There is a total of $(length(Arreglo_Centroides)) polygons with possible duplicates.")
    #println("There is a total of $(length(Arreglo_Centroides)) uniques polygons.")
    
    #Paso 6: Generamos un único diccionario al unir los 4 diccionarios generados previamente
    Diccionario_Centroides_Unico = merge(Diccionario_Centroides1, Diccionario_Centroides2, Diccionario_Centroides3, Diccionario_Centroides4)

    #Paso 7: Recuperamos los polígonos del arrego cuasiperiódico que corresponden a los clúster principales.
    Coordenadas_X_CP, Coordenadas_Y_CP = centroides_A_Vertices(Arreglo_Centroides, Diccionario_Centroides_Unico);
    
    return Coordenadas_X_CP, Coordenadas_Y_CP
end
