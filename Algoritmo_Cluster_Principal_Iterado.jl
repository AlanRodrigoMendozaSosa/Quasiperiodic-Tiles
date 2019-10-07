#Función que itera sobre todos los centroides de los polígonos del arreglo cuasiperiódico. Si el centroide, al pasar a
#polígonos de Voronoi tiene asociado un polígono con área menor a una "Cota_Area", entonces guardamos dicho centroides
#caso contrario lo eliminamos.
function centroides_Area_Acotada(Voronoi, Cota_Area)
    Arreglo_Centroides = []; #Arreglo donde irán los centroides cuya área cae debajo de la cota dada.
    
    for Face in Voronoi.faces #Iteramos sobre todos los centroides.
        
        Halfedge = Face.outerComponent; #Elemento de los pol. Voronoi, si es igual a "nothing" el pol. no se cierra
        Cerrado = true; #Condicional para cortar el "while" una vez que se ha encontrado que el pol. no cierra
        
        if Halfedge == nothing #Condición indicadora de que el pol. no cierra
            Cerrado = false; #Ponemos el condicional como falso para que el While no ejecute
        end
        
        while Cerrado && Halfedge.next != Face.outerComponent #Condiciones para seguir buscando si el pol. se cierra

            if Halfedge.next == nothing #Condición indicadora de que el pol. no cierra
                Cerrado = false; #Ponemos el condicional como falso para que el While no ejecute
            end

            Halfedge = Halfedge.next; #Probamos con el siguiente halfedge del polígono
            
        end
        
        if Cerrado
            
            Halfedge = Face.outerComponent; #Elemento de los pol. Voronoi, si es igual a "nothing" el pol. no se cierra
            Vertices_Poligono = [[Halfedge.origin.coordinates[1], Halfedge.origin.coordinates[2]]];
            
            while Halfedge.next != Face.outerComponent
                Halfedge = Halfedge.next
                push!(Vertices_Poligono, [Halfedge.origin.coordinates[1], Halfedge.origin.coordinates[2]]); 
            end
            
            Poligono = Polygon(Vertices_Poligono); #Generamos la estructura adecuada para la función área de Ata
            
            if area(Poligono) < Cota_Area
                push!(Arreglo_Centroides, Face.site)
            end
            
        end
        
    end
    
    return Arreglo_Centroides
end

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