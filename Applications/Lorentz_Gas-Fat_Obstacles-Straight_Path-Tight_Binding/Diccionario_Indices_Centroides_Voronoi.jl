#Definimos una función que nos genere un diccionario de las duplas de los centroides en la lista dada por el usuario al índice del polígono asociado 
#en la red de Voronoi.
#Arreglo_Duplas_Centroides: Arreglo con las duplas (X,Y) asociadas a las coordenadas de los centroides de los polígonos del arreglo cuasiperiódico
#Voronoi: Estructura de datos generada por el algoritmo para los polígonos de Voronoi
function diccionario_Centroides_Indice_Voronoi(Arreglo_Duplas_Centroides, Voronoi)
    Diccionario_Centroides_Indice = Dict(); #Definimos el diccionario vacío que iremos llenando

    for i in 1:length(Arreglo_Duplas_Centroides)
        Diccionario_Centroides_Indice[Voronoi.faces[i].site] = i; #Añadimos la llave de la dupla del centroide y le asociamos su índice en los polígonos de Voronoi
    end

    return Diccionario_Centroides_Indice
end