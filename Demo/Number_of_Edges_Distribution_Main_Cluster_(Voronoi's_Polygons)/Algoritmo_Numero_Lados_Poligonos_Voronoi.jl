#Función que nos regresa el número de lados de los polígonos de Voronoi de los clúster principales (antes de remover las capas externas para evitar deformar los 
#polígonos cercanos a la última capa del cluster). La función se itera un cierto número de veces para tener una mayor estadística.
#"Iteraciones" es el número de clúster principales que vamos a considerar para el número de lados.
#"Iteraciones_Cluster" es el número de iteraciones del algoritmo de acotar polígonos de Voronoi por su área que se empleará para determinar un clúster principal
#"Area_Cota" es el área que servirá como discriminante de buenos y malos polígonos en el algoritmo de acotar poligonos
#"Margen_Error" es el parámetro que determina qué margen de error se considera al generar los posibles números enteros asociados al polígono de Voronoi que contiene 
#al punto arbitrario.
#"Semilado_Caja" es el semilado de un cuadrado en el cual se genera un punto arbitrario, alrededor del cual se genera la vecindad del arreglo cuasiperiódico.
function arreglo_Numero_Lados_Buenos_Poligonos(Iteraciones, Iteraciones_Cluster, Area_Cota, Margen_Error, Semilado_Caja, Promedios_Distancia, Vectores_Estrella, Arreglo_Alfas)
    Arreglo_Numero_Lados_Total = []; #Arreglo donde irán el número de lados de los polígonos tras todas las iteraciones
    Contador_Iteraciones = 0; #Contador de las iteraciones que se han realizado
    Arreglo_Arreglos_Vertices = []; #Arreglo que contendrá los arreglos con los vértices de cada polígono
    
    for i in 1:Iteraciones
        #Arreglo donde irán los vértices de los polígonos que viven en el cluster principal de cada una de las vecindades generadas
        Arreglo_Poligonos_Cluster_Principal = [];

        #Generación del arreglo cuasiperiódico con margen de error (1era entrada) y semilado del cuadrado donde colocar un pto arbitrario.
        Coordenadas_X, Coordenadas_Y, Punto = region_Local_Voronoi(Margen_Error, Semilado_Caja, Promedios_Distancia, Vectores_Estrella, Arreglo_Alfas);

        #Obtención de los centroides de los poligonos del arreglo cuasiperiódico en: 
        #coordenada X, coordenada Y, Float64(X,Y) y su diccionario Centroides <-> Vertices.
        Centroides_X, Centroides_Y, Centroides, Diccionario_Centroides = centroides(Coordenadas_X, Coordenadas_Y);
        
        Centroides_Iterados = centroides_Area_Acotada_Iterada(Centroides, Area_Cota, Iteraciones_Cluster);
        
        #Definimos las duplas requeridas por el algoritmo "Voronoi" de Enrique con las coordenadas de los centroides iniciales
        sites_Inicial = [(Centroides[i][1], Centroides[i][2]) for i in 1:length(Centroides)];

        #Obtenemos el arreglo de Voronoi con los centroides iniciales
        voronoi_Inicial = getVoronoiDiagram(sites_Inicial);

        for Face in voronoi_Inicial.faces #Iteramos sobre todos los centroides
    
            if Face.site in Centroides_Iterados #Si es cierto, el polígono es uno de nuestros polígonos de interés
                Arreglo_Coordenadas_Vertices = []; #Arreglo donde iremos poniendo las coordenadas de los vértices de los pol. Voronoi
                Halfedge = Face.outerComponent; #Elemento de los pol. Voronoi, si es igual a "nothing" el pol. no se cierra
                Vertice_Inicial = Halfedge.origin.coordinates; #Coordenadas del primer vértice
                push!(Arreglo_Coordenadas_Vertices, Vertice_Inicial)
                Contador_Vertices = 1; #Contador que llevará el número de vértices que tiene el polígono

                while Halfedge.next != Face.outerComponent
                    Halfedge = Halfedge.next; #Nos movemos al siguiente elemento del polígono de Voronoi
                    Vertice = Halfedge.origin.coordinates; #Coordenadas del vértice

                    if Vertice in Arreglo_Coordenadas_Vertices
                        nothing
                    else
                        #Comparamos con el vértice previo para ver si el lado que se añade es demasiado pequeño, si es el caso lo quitamos.
                        if (Vertice[1] - Arreglo_Coordenadas_Vertices[end][1])^2 + (Vertice[2] - Arreglo_Coordenadas_Vertices[end][2])^2 > 1e-6 && (Vertice[1] - Arreglo_Coordenadas_Vertices[1][1])^2 + (Vertice[2] - Arreglo_Coordenadas_Vertices[1][2])^2 > 1e-6
                            push!(Arreglo_Coordenadas_Vertices, Vertice);
                            Contador_Vertices += 1; #Sumamos uno en el contador de vértices
                        end
                    end
                end

                push!(Arreglo_Numero_Lados_Total, Contador_Vertices);
                push!(Arreglo_Poligonos_Cluster_Principal, Arreglo_Coordenadas_Vertices);
                #push!(Arreglo_Arreglos_Vertices, Arreglo_Coordenadas_Vertices)

            end

        end

        #Agregamos al arreglo donde irán todos los vértices de los polígonos de los clusteres principales, el arreglo con los polígonos del cluster principal de
        #esta iteración
        push!(Arreglo_Arreglos_Vertices, Arreglo_Poligonos_Cluster_Principal);
        
        Contador_Iteraciones += 1;
        println("Se han realizado ", Contador_Iteraciones, " iteraciones.")
        
    end
    
    return Arreglo_Numero_Lados_Total, Arreglo_Arreglos_Vertices
end