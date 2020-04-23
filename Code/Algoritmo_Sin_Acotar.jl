#Función que genera los vértices de un arreglo cuasiperiódico asociados a la vecindad de un Punto arbitrario (dado o no por el usuario).
#"N" es el margen de error numérico al determinar el número entero de la proyección del punto de interés con los vectores estrella.
#"Promedios_Distancia" es el arreglo con la separación entre las franjas cuasiperiódicas.
#"Vectores_Estrella" es el arreglo con los vectores estrella que generan la retícula cuasiperiódica deseada.
#"Arreglo_Alfas" es el arreglo con los valores numéricos de la separación respecto al origen del conjunto de rectas ortogonales a los vectores estrella en el método generalizado dual.
#"Punto" es el punto alrededor de donde se va a generar la vecindad.
function region_Local(N, Promedios_Distancia, Vectores_Estrella, Arreglo_Alfas, Punto)
    #Paso 1: Dado el Punto proyectamos este con los Vectores Estrella para obtener los enteros aproximados asociados al polígono contenedor.
    #Salida: [n1,n2,n3,...]    
    Proyecciones = proyecciones_Pto_Direccion_Franjas(Punto, Promedios_Distancia, Vectores_Estrella);
    
    #Paso 2: A partir de los valores enteros aproximados, generamos la vecindad del arreglo cuasiperiódico que contenga al punto.
    #Salida: [[X,Y]]
    Puntos_Duales = generador_Vecindades_Vertices(Proyecciones, Vectores_Estrella, Arreglo_Alfas, N);
    
    return Puntos_Duales
end

#Función que determina el polígono que contiene al punto de interés.
#"Coordenadas_X" es un arreglo con las coordenadas en X de los vértices de los polígonos. Cada 4 entradas corresponden a un mismo polígono.
#"Coordenadas_Y" es un arreglo con las coordenadas en Y de los vértices de los polígonos. Cada 4 entradas corresponden a un mismo polígono.
#"Punto" es un arreglo dos dimensional con las coordenadas de un punto en el espacio 2D.
function poligono_Contenedor(Coordenadas_X, Coordenadas_Y, Punto)
    #Generamos los segmentos de los polígonos de la retícula cuasiperiódica.
    Poligonos = obtener_Segmentos_Vertices(Coordenadas_X, Coordenadas_Y);
    #Obtenemos los vértices del polígono contenedor
    Vertices_Contenedor = encontrar_Poligono(Punto, Poligonos);
    
    return Vertices_Contenedor
end