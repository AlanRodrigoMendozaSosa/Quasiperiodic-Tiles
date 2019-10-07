#Función que regresa las coordenadas de los centroides de aquellos polígonos que no se cierran al pasar a Voronoi.
function centroides_Sin_Poligono(Voronoi)
    Sitios = []; #Arreglo donde irán las coordenadas de los sitios de los polígonos que no se cierran
    for Face in Voronoi.faces #Iteramos sobre todos los centroides.

        Halfedge = Face.outerComponent; #Elemento de los pol. Voronoi, si es igual a "nothing" el pol. no se cierra
        llave = true; #Condicional para cortar el "while" una vez que se ha encontrado que el pol. no cierra

        if Halfedge == nothing #Condición indicadora de que el pol. no cierra
            push!(Sitios, Face.site); #Si no cierra mandamos las coordenadas del centroide al arreglo Sitios
            llave = false; #Ponemos el condicional como falso para que el While no ejecute
        end

        while llave && Halfedge.next != Face.outerComponent #Condiciones para seguir buscando si el pol. se cierra

            if Halfedge.next == nothing #Condición indicadora de que el pol. no cierra
                push!(Sitios, Face.site); #Si no cierra mandamos las coordenadas del centroide al arreglo Sitios
                llave = false; #Ponemos el condicional como falso para que el While no ejecute
            end

            Halfedge = Halfedge.next; #Probamos con el siguiente halfedge del polígono
        end

    end
    
    return Sitios
end

#Función que regresa los Centroides de los polígonos que quedan tras ir quitando las capas más externas del arreglo Cuasiperiódico (La cebolla)
function algoritmo_Sir_Davos!(Centroides, Num_Capas_Cebolla)
    
    for i in 1:Num_Capas_Cebolla
        #Calculamos los centroides que se usarán para generar Voronoi
        sites = [(Centroides[i][1], Centroides[i][2]) for i in 1:length(Centroides)]

        #Aplicamos el algoritmo de Voronoi
        voronoi = getVoronoiDiagram(sites);

        #Obtenemos los centroides que no cierran en la iteración de Voronoi
        Sitios = centroides_Sin_Poligono(voronoi);

        #Removemos los centroides que quedan fuera en la capa de la cebolla
        for i in Sitios
            deleteat!(Centroides,findfirst(x -> x == i, Centroides))
        end
    end
    
    return Centroides
end

#Polygon es una variable que contiene los vértices de forma ordenada. Podrías remplazarla por un arreglo de vértices 
#ordenados. 
mutable struct Polygon
    Vertices::Array{Array{Number,1}, 1}
end

#Función de Ata para calcular áreas de polígonos
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

#Obtengamos un arreglo con las áreas de los polígonos que están bien definidos en Voronoi.
function area_Poligonos_Voronoi_Config_Inicial(Centroides_Sin_Cerrar, Voronoi)
    Areas = []; #arreglo que contendrá las áreas de los polígonos
    
    for Face in Voronoi.faces #Iteramos sobre todos los centroides

        if Face.site in Centroides_Sin_Cerrar #Si es "true" no hacemos nada, el polígono no tiene área
            nothing
        else
            Halfedge = Face.outerComponent; #Elemento de los pol. Voronoi, si es igual a "nothing" el pol. no se cierra
            Vertices_Poligono = [[Halfedge.origin.coordinates[1], Halfedge.origin.coordinates[2]]];
            
            while Halfedge.next != Face.outerComponent
                Halfedge = Halfedge.next
                push!(Vertices_Poligono, [Halfedge.origin.coordinates[1], Halfedge.origin.coordinates[2]]); 
            end
            
            Poligono = Polygon(Vertices_Poligono); #Generamos la estructura adecuada para la función área de Ata
            
            push!(Areas,area(Poligono))
        end
        
    end
    
    return Areas
end

#Obtengamos un arreglo con las áreas de los polígonos que están bien definidos en Voronoi.
function area_Poligonos_Voronoi_Ultima_Capa(Centroides_Refinados, Voronoi_Config_Inicial)
    Areas = []; #arreglo que contendrá las áreas de los polígonos
    
    for Face in Voronoi_Config_Inicial.faces #Iteramos sobre todos los centroides

        if Face.site in Centroides_Refinados #Si es "true" es uno de los buenos poligonos
            Halfedge = Face.outerComponent; #Elemento de los pol. Voronoi, si es igual a "nothing" el pol. no se cierra
            Vertices_Poligono = [[Halfedge.origin.coordinates[1], Halfedge.origin.coordinates[2]]];
            
            while Halfedge.next != Face.outerComponent
                Halfedge = Halfedge.next
                push!(Vertices_Poligono, [Halfedge.origin.coordinates[1], Halfedge.origin.coordinates[2]]); 
            end
            
            Poligono = Polygon(Vertices_Poligono); #Generamos la estructura adecuada para la función área de Ata
            
            push!(Areas,area(Poligono))
        end
        
    end
    
    return Areas
end

#Función que nos regresa las áreas de los polígonos de Voronoi tras remover una "Num_Capas_Remover" cantidad de capas externas a dicho arreglo.
#A priori, conocer el número de capas a remover para quedarnos con el cluster principal es prácticamente imposible (no conozco un método de hacerlo analíticamente)
#Por ello es que al implementarse este algoritmo, se hacen previamente pruebas con distintos valores para la variable "Num_Capas_Remover".
function arreglo_Areas_Buenos_Poligonos(Iteraciones, Num_Capas_Remover, Margen_Error, Semilado_Caja)
    Arreglo_Buenas_Areas = []; #Arreglo donde meteremos las buenas áreas de los polígonos de Voronoi.
    Contador_Iteraciones = 0;
    
    for i in 1:Iteraciones
    
        #Generación del arreglo cuasiperiódico con margen de error (1era entrada) y semilado del cuadrado donde colocar un 
        #pto arbitrario.
        Coordenadas_X, Coordenadas_Y, Punto = region_Local_Voronoi(Margen_Error, Semilado_Caja);

        #Obtención de los centroides de los poligonos del arreglo cuasiperiódico en: 
        #coordenada X, coordenada Y, Float64(X,Y) y su diccionario Centroides <-> Vertices.
        Centroides_X, Centroides_Y, Centroides, Diccionario_Centroides = centroides(Coordenadas_X, Coordenadas_Y);

        #Definimos las duplas requeridas por el algoritmo "Voronoi" de Enrique con las coordenadas de los centroides.
        sites = [(Float64(Centroides_X[i]), Float64(Centroides_Y[i])) for i in 1:length(Centroides_X)];

        #Generamos la estructura de datos para los polígonos de Voronoi de los centroides de los polígonos del arreglo cuasipe.
        voronoi_inicial = getVoronoiDiagram(sites);

        #Arreglo con los centroides de los polígonos a quitar
        Sitios_Centroides_Sin_Cerrar = centroides_Sin_Poligono(voronoi_inicial);
        
        #Obtenemos las área de todos los polígonos de Voronoi (exceptuando únicamente los "polígonos" que no se cierran).
        #Es decir, iteramos sobre todos los centroides, si el centroide corresponde a un polígono que no se cierra en Voronoi
        #lo ignoramos, caso contrario obtenemos el área del polígono.
        Areas_Config_Inicial = area_Poligonos_Voronoi_Config_Inicial(Sitios_Centroides_Sin_Cerrar, voronoi_inicial);
        
        #Hacemos una copia del arreglo de Centroides para poder conservar la configuración inicial y poder hacer distintas 
        #pruebas con esa misma configuración.
        Centroides_Copia = copy(Centroides);

        #Removemos capa externa tras capa externa de los polígonos de Voronoi hasta que mantengamos únicamente la capa central
        #conformada por buenos polígonos de Voronoi. Esta parte es visual, con lo cual el número de capas a quitar puede ser
        #variable.
        Centroides_Copia = algoritmo_Sir_Davos!(Centroides_Copia, Num_Capas_Remover);
        
        #Obtengamos el área de los buenos polígonos de Voronoi. Es decir, iteramos sobre la configuración inicial (esto es 
        #importante ya que, incluso si son buenos polígonos, al aplicarles Voronoi a estos siempre quedarán polígonos abiertos
        #y con áreas grandes, por eso es que el cálculo de las áreas debe ser sobre los buenos polígonos pero en el arreglo de
        #polígonos de Voronoi inicial) con todos los centroides, si el centroide es el asociado a un buen polígono se cálcula
        #su área, caso contrario no.
        Areas_Ultima_Capa = area_Poligonos_Voronoi_Ultima_Capa(Centroides_Copia, voronoi_inicial);
        
        #Área máxima y mínima de los buenos polígonos
        #println("El área máxima de los buenos polígonos es: ", maximum(Areas_Ultima_Capa))
        #println("El área mínima de los buenos polígonos es: ", minimum(Areas_Ultima_Capa))
        
        push!(Arreglo_Buenas_Areas, Areas_Ultima_Capa);
        Contador_Iteraciones += 1;
        println("Se han realizado ", Contador_Iteraciones, " iteraciones.")
    end
    
    return collect(Iterators.flatten(Arreglo_Buenas_Areas))
end
