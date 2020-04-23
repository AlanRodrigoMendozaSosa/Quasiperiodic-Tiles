#Función que elimina la capa más externa de uno o varios clústers de polígonos.
#"Voronoi" es la estructura generada por Enrique con getVoronoiDiagram().
#"Cota_Area" es el área que servirá como discriminante parar separar a los polígonos del clúster principal de los que no están ahí empleando el área de los polígonos de Voronoi.
function centroides_Area_Acotada(Voronoi, Cota_Area)
    #Paso 1: Definimos un arreglo donde irán los centroides cuya área del polígono de Voronoi asociado cae debajo de la cota dada.
    Arreglo_Centroides = [];

    #Paso 2: Iteramos sobre todos los polígonos asociados a centroides
    for Face in Voronoi.faces
        if Face.area < Cota_Area #Compara el área del polígono en turno contra el área de la cota dada
            push!(Arreglo_Centroides, Face.site) #Agregamos al arreglo correspondiente las coordenadas (X,Y) del centroide que da lugar al polígono
        end
    end

    return Arreglo_Centroides
end

#Función que itera varias veces el algoritmo de centroides_Area_Acotada(). Esto sirve para ir quitando varias capas usando como discriminante el área de los polígonos de Voronoi.
#"Centroides" es un arreglo con los centroides de los polígonos del arreglo cuasiperiódico de interés.
#"Cota_Area" es el área que servirá como discriminante parar separar a los polígonos del clúster principal de los que no están ahí empleando el área de los polígonos de Voronoi.
#"Num_Iteraciones" es un número entero asociado al número de veces que queremos implementar el algoritmo.
function centroides_Area_Acotada_Iterada(Centroides, Cota_Area, Num_Iteraciones)
    for i in 1:Num_Iteraciones
        #Definimos las duplas requeridas por el algoritmo "Voronoi" de Enrique con las coordenadas de los centroides.
        sites = [(Centroides[j][1], Centroides[j][2]) for j in 1:length(Centroides)]

        #Generamos la estructura de datos para los polígonos de Voronoi de los centroides de los polígonos del arreglo cuasiperiódico.
        voronoi_inicial = getVoronoiDiagram(sites);

        #Reescribimos el conjunto "Centroides" al eliminar aquellos que, en la iteración i-ésima, obtuvieron un polígono de Voronoi con área mayor a "Cota_Area"
        #Salida: [Float64(X,Y)]
        Centroides = centroides_Area_Acotada(voronoi_inicial, Cota_Area);
    end
    
    return Centroides
end