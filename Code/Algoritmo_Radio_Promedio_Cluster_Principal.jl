#Función que calcula el radio promedio del cluster principal de diferentes configuraciones.
#"Margen_Error" es el número entero que consideramos como error para las proyecciones del punto de interés con los vectores estrella.
#"Semilado" es el número flotante asociado al semilado de un cuadrado centrado en el origen dentro del cual pondremos un punto aleatorio.
#"Iteraciones_Areas" es el número entero de veces que aplicaremos el algoritmo de quitar polígonos por el discriminante del área a cada arreglo de centroides.
#"Cota_Area" es un número flotante que sirve como área discrminante de polígonos del cluster principal de los que no.
#"Iteraciones_Algoritmo" es el número entero de puntos y sus vecindades del arreglo cuasiperiódico que consideraremos.
function radio_Seguro(Margen_Error, Semilado, Iteraciones_Areas, Cota_Area, Iteraciones_Algoritmo, Promedios_Distancia, Vectores_Estrella, Arreglo_Alfas)
    Radio_Seguro = 0; #Definimos la variable que nos dará el radio del parche

    for i in 1:Iteraciones_Algoritmo
        #Generamos un punto arbitrario en un cuadrado centrado en el origen de semilado "Semilado"
        Punto_Arbitrario = punto_Arbitrario(Semilado);

        #Generación del arreglo cuasiperiódico con margen de error (1era entrada)
        Coordenadas_Vertices = region_Local_Voronoi(Margen_Error, Promedios_Distancia, Vectores_Estrella, Arreglo_Alfas, Punto_Arbitrario);

        #Obtención de los centroides de los poligonos del arreglo cuasiperiódico en: 
        #coordenada X, coordenada Y, Float64(X,Y) y su diccionario Centroides <-> Vertices.
        Centroides, Diccionario_Centroides = centroides(Coordenadas_Vertices);
        
        #Obtención de los centroides de los buenos polígonos tras iterar N veces (3era entrada) algoritmo de áreas.
        Centroides_Iterados = centroides_Area_Acotada_Iterada(Centroides, Cota_Area, Iteraciones_Areas);
        
        #Visualizamos ahora el arreglo cuasiperiódico tras quitarle las n capas.
        Coordenadas_X, Coordenadas_Y = centroides_A_Vertices(Centroides_Iterados, Diccionario_Centroides);
        
        #Algoritmo para encontrar el radio seguro de los polígonos de Voronoi en un arreglo cuadrado.
        MaximoX = maximum(Coordenadas_X)
        MaximoY = maximum(Coordenadas_Y)
        MinimoX = minimum(Coordenadas_X)
        MinimoY = minimum(Coordenadas_Y)

        Distancias = [abs(MaximoX-Punto_Arbitrario[1]), abs(MinimoX-Punto_Arbitrario[1]), abs(MaximoY-Punto_Arbitrario[2]), abs(MinimoY-Punto_Arbitrario[2])]

        Radio_Seguro += minimum(Distancias)
        
        println("Se han realizado $(i) iteraciones.")
    end
    
    return Radio_Seguro/Iteraciones_Algoritmo
end