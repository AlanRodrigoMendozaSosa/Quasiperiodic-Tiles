#Algoritmo para conocer el radio promedio del cluster principal de diferentes configuraciones.
function radio_Seguro(Margen_Error, Semilado, Iteraciones_Areas, Cota_Area, Iteraciones_Algoritmo)
    Radio_Seguro = 0; #Definimos la variable que nos dará el radio del parche
    Contador_Iteraciones = 0; #Contador de las iteraciones que se han realizado

    for i in 1:Iteraciones_Algoritmo
        #Generación del arreglo cuasiperiódico con margen de error (1era entrada) y semilado del cuadrado donde 
        #colocar un pto arbitrario.
        Coordenadas_X, Coordenadas_Y, Punto = region_Local_Voronoi(Margen_Error, Semilado);

        #Obtención de los centroides de los poligonos del arreglo cuasiperiódico en: 
        #coordenada X, coordenada Y, Float64(X,Y) y su diccionario Centroides <-> Vertices.
        Centroides_X, Centroides_Y, Centroides, Diccionario_Centroides = centroides(Coordenadas_X, Coordenadas_Y);
        
        #Obtención de los centroides de los buenos polígonos tras iterar N veces (3era entrada) algoritmo de áreas.
        Centroides_Iterados = centroides_Area_Acotada_Iterada(Centroides, Cota_Area, Iteraciones_Areas);
        
        #Visualizamos ahora el arreglo cuasiperiódico tras quitarle las n capas.
        Coordenadas_X_2, Coordenadas_Y_2 = Centroides_a_Vertices(Centroides_Iterados, Diccionario_Centroides);
        
        #Algoritmo para encontrar el radio seguro de los polígonos de Voronoi en un arreglo cuadrado.
        MaximoX = maximum(Coordenadas_X_2)
        MaximoY = maximum(Coordenadas_Y_2)
        MinimoX = minimum(Coordenadas_X_2)
        MinimoY = minimum(Coordenadas_Y_2)

        Distancias = [abs(MaximoX-Punto[1]), abs(MinimoX-Punto[1]), abs(MaximoY-Punto[2]), abs(MinimoY-Punto[2])]

        Radio_Seguro += minimum(Distancias)
        
        Contador_Iteraciones += 1;
        println("Se han realizado ", Contador_Iteraciones, " iteraciones.")
    end
    
    return Radio_Seguro/Iteraciones_Algoritmo
end