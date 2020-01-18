#=
#Polygon es una variable que contiene los vértices de forma ordenada. Podrías remplazarla por un arreglo de vértices 
#ordenados. 
mutable struct Polygon
    Vertices::Array{Array{Number,1}, 1}
end

#Función de Ata para calcular áreas de polígonos.
#"P" es un polígono.
function area(P::Polygon)
    A = 0; #Variable que eventualmente nos dará el área del polígono
    Vertices = P.Vertices; #Arreglo con los vértices del polígono (presuponiendo que están ordenados)
    N = length(Vertices); #Número de vértices que conforman al polígono
    for i in 1:N-1
        A += Vertices[i][1]*Vertices[i+1][2]-Vertices[i][2]*Vertices[i+1][1]; #Algoritmo para determinar el área.
    end
    A += Vertices[N][1]*Vertices[1][2]- Vertices[1][1]*Vertices[N][2]; #Iteramos el algoritmo con los puntos 1 y N
    A *= 1/2; #Terminamos dividiendo entre dos.
    if A < 0
        return -A
    else
        return A
    end
end
=#

#Función que itera sobre todos los centroides de los polígonos del arreglo cuasiperiódico. Si el centroide, al pasar a polígonos de Voronoi tiene asociado un 
#polígono con área menor a una "Cota_Area", entonces guardamos dicho centroides caso contrario lo eliminamos.
#"Voronoi" es la estructura generada por Enrique con getVoronoiDiagram()-
#"Cota_Area" es el área que servirá como discriminante de los polígonos de Voronoi del clúster principal de los que no están ahí.
function centroides_Area_Acotada(Voronoi, Cota_Area)
    Arreglo_Centroides = []; #Arreglo donde irán los centroides cuya área cae debajo de la cota dada.

    for Face in Voronoi.faces #Iteramos sobre todos los polígonos asociados a centroides
        if Face.area < Cota_Area #Compara el área del polígono en turno contra el área de la cota dada
            push!(Arreglo_Centroides, Face.site) #Agregamos al arreglo correspondiente las coordenadas (X,Y) del centroide que da lugar al polígono
        end
    end

    return Arreglo_Centroides
end

#Función que itera varias veces el algoritmo de centroides_Area_Acotada(). Esto sirve para ir quitando varias capas usando como discriminante el área.
#"Centroides" es un arreglo con los centroides de los polígonos del arreglo cuasiperiódico de interés.
#"Cota_Area" es el área que servirá como discriminante de los polígonos de Voronoi del clúster principal de los que no están ahí.
#"Num_Iteraciones" es un número entero asociado al número de veces que queremos implementar el algoritmo.
function centroides_Area_Acotada_Iterada(Centroides, Cota_Area, Num_Iteraciones)
    for i in 1:Num_Iteraciones
        #Definimos las duplas requeridas por el algoritmo "Voronoi" de Enrique con las coordenadas de los centroides.
        sites = [(Centroides[i][1], Centroides[i][2]) for i in 1:length(Centroides)]

        #Generamos la estructura de datos para los polígonos de Voronoi de los centroides de los polígonos del arreglo cuasipe.
        voronoi_inicial = getVoronoiDiagram(sites);

        #Arreglo con los centroides del Cluster Principal
        Centroides = centroides_Area_Acotada(voronoi_inicial, Cota_Area);
    end
    
    return Centroides
end