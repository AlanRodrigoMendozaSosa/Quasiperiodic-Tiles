#Función que mantiene los centroides de las teselas del cluster principal
#"Voronoi" es la estructura generada por Enrique con getVoronoiDiagram().
#"Cota_Area" es el área que servirá como discriminante parar separar a los polígonos dentro de clúster de los de frontera
#"APoint" es el punto en el espacio alrededor del cual se genera la vecindad cuasiperiódica
function centroides_Area_Acotada(Voronoi, Cota_Area, APoint)
    Arreglo_Centroides = []; #Arreglo con Centroides de Celdas Voronoi del cluster principal
    Distancia_Frontera_Cercana = Inf;

    #Iteramos sobre todos los polígonos asociados a centroides
    for Face in Voronoi.faces
        if Face.area > Cota_Area #Compara el área del polígono en turno contra el área de la cota dada
            Distancia_Celda = norm([Face.site[1], Face.site[2]] - APoint);
            Distancia_Celda  < Distancia_Frontera_Cercana ? (Distancia_Frontera_Cercana = Distancia_Celda) : (nothing);
        end
    end
    
    for Face in Voronoi.faces
        if Face.area < Cota_Area #Compara el área del polígono en turno contra el área de la cota dada
            Distancia_Celda = norm([Face.site[1], Face.site[2]] - APoint);
            Distancia_Celda  < Distancia_Frontera_Cercana ? (push!(Arreglo_Centroides, Face.site)) : (nothing)
        end
    end

    return Arreglo_Centroides
end