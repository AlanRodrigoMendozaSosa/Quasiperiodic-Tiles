#Función para generar un parche cuadrado con vértices de las celdas unitarias del arreglo cuasiperiódico.
#"Informacion_Parche" es el arreglo con todos los datos asociados a los parches circulares.
#"Factor_Reduccion" es el factor por el que se multiplica el radio promedio de los parches circulares para generar el radio seguro de los parches.
#"Promedios_Distancia" es el arreglo con la separación entre las franjas cuasiperiódicas.
#"Vectores_Estrella" es el arreglo con los vectores estrella que generan la retícula cuasiperiódica deseada.
#"Arreglo_Alfas" es el arreglo con los valores numéricos de la separación respecto al origen del conjunto de rectas ortogonales a los vectores estrella en el método generalizado dual.
#"Punto" es el punto alrededor de donde se va a generar la vecindad.
function parche_Cuadrado(Informacion_Parche, Factor_Reduccion, Promedios_Distancia, Vectores_Estrella, Arreglo_Alfas, Punto)    
    #Paso 1: Generemos los cuatro centros para los parches circulares.
    #Informacion_Parche[5] es el radio promedio del parche circular.
    #Salida: [X,Y]
    P_1 = Punto + [(sqrt(2)/2)*(Factor_Reduccion)*(Informacion_Parche[5]), (sqrt(2)/2)*(Factor_Reduccion)*(Informacion_Parche[5])];
    P_2 = Punto + [-(sqrt(2)/2)*(Factor_Reduccion)*(Informacion_Parche[5]), (sqrt(2)/2)*(Factor_Reduccion)*(Informacion_Parche[5])];
    P_3 = Punto - [(sqrt(2)/2)*(Factor_Reduccion)*(Informacion_Parche[5]), (sqrt(2)/2)*(Factor_Reduccion)*(Informacion_Parche[5])];
    P_4 = Punto + [(sqrt(2)/2)*(Factor_Reduccion)*(Informacion_Parche[5]), -(sqrt(2)/2)*(Factor_Reduccion)*(Informacion_Parche[5])];
    
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

    #Paso 4: Generamos los teselados de Voronoi asociados a los centroides.
    sites1 = [(Centroides1[i][1], Centroides1[i][2]) for i in 1:length(Centroides1)]
    sites2 = [(Centroides2[i][1], Centroides2[i][2]) for i in 1:length(Centroides2)]
    sites3 = [(Centroides3[i][1], Centroides3[i][2]) for i in 1:length(Centroides3)]
    sites4 = [(Centroides4[i][1], Centroides4[i][2]) for i in 1:length(Centroides4)]
    Initial_Voronoi1 = getVoronoiDiagram(sites1);
    Initial_Voronoi2 = getVoronoiDiagram(sites2);
    Initial_Voronoi3 = getVoronoiDiagram(sites3);
    Initial_Voronoi4 = getVoronoiDiagram(sites4);
    
    #Paso 5: Depuramos los arreglos cuasiperiódicos obtenidos para quedarnos con el clúster principal.
    #Informacion_Parche[3] es la cota del área empleada como discriminante en el algoritmo de Voronoi.
    #Salida: [Float64(X,Y)]
    Centroides_Principal1 = centroides_Area_Acotada(Initial_Voronoi1, Informacion_Parche[3], P_1);
    Centroides_Principal2 = centroides_Area_Acotada(Initial_Voronoi2, Informacion_Parche[3], P_2);
    Centroides_Principal3 = centroides_Area_Acotada(Initial_Voronoi3, Informacion_Parche[3], P_3);
    Centroides_Principal4 = centroides_Area_Acotada(Initial_Voronoi4, Informacion_Parche[3], P_4);
    
    #Paso 6: Agrupamos todos los buenos centroides en un único arreglo, quedándonos únicamente con los centroides sin repetir.
    Arreglo_Centroides = vcat(Centroides_Principal1, Centroides_Principal2, Centroides_Principal3, Centroides_Principal4);
    unique!(Arreglo_Centroides);
    #println("There is a total of $(length(Arreglo_Centroides)) polygons with possible duplicates.")
    #println("There is a total of $(length(Arreglo_Centroides)) uniques polygons.")
    
    #Paso 7: Generamos un único diccionario al unir los 4 diccionarios generados previamente
    Diccionario_Centroides_Unico = merge(Diccionario_Centroides1, Diccionario_Centroides2, Diccionario_Centroides3, Diccionario_Centroides4)

    #Paso 8: Recuperamos los polígonos del arreglo cuasiperiódico que corresponden a los clúster principales.
    Coordenadas_X_CP, Coordenadas_Y_CP = centroides_A_Vertices(Arreglo_Centroides, Diccionario_Centroides_Unico);
    
    return Coordenadas_X_CP, Coordenadas_Y_CP
end
