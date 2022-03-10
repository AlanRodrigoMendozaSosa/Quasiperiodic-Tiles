#Función que regresa las coordenadas de los centroides de aquellos polígonos en la capa externa del teselado de Voronoi.
#"Voronoi" es la estructura generada por Enrique con la función getVoronoiDiagram().
function centroides_Sin_Poligono(Voronoi)
    Sitios = []; #Arreglo donde irán las coordenadas de los sitios de los polígonos que no se cierran
    for Face in Voronoi.faces #Iteramos sobre todos los centroides.
        if Face.area == Inf #Indicador de que el polígono está asociado a un centroide de la capa externa
            push!(Sitios, Face.site)
        end
    end    
    return Sitios
end

#Función que regresa los Centroides de los polígonos que quedan tras ir quitando las capas más externas del arreglo Cuasiperiódico (La cebolla).
#NOTA: Esta función modifica el arreglo Centroides, por lo que en caso de querer conservarlo es necesario trabajar con una copia.
#"Centroides" es un arreglo con los centroides de todos los polígonos generados en la vecindad del arreglo cuasiperiódico.
#"Num_Capas_Cebolla" es un entero que indica cuántas iteraciones del algoritmo centroides_Sin_Poligono() se van a realizar.
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

#Función que calcula las áreas de los polígonos que están en el clúster principal de Voronoi.
#"Centroides_Refinados" es un arreglo con los centroides tras quitarle un cierto número de capas al arreglo de Voronoi.
#"Diccionario_Centroides_Indices" es un diccionario que relaciona las coordenadas del centroide de un polígono con el índice de su polígono de Voronoi
#"Voronoi_Config_Inicial" es la estructura generada por Enrique con la función getVoronoiDiagram() considerando TODOS los centroides.
function area_Poligonos_Voronoi_Ultima_Capa(Centroides_Refinados, Diccionario_Centroides_Indices, Voronoi_Config_Inicial)
    Areas_Cluster_Ultima_Capa = []; #Arreglo donde irán las área de los polígonos que corresponden al cluster principal
    for i in Centroides_Refinados
        push!(Areas_Cluster_Ultima_Capa, Voronoi_Config_Inicial.faces[Diccionario_Centroides_Indices[i]].area);
    end

    return Areas_Cluster_Ultima_Capa
end

#Función que nos regresa las áreas de los polígonos de Voronoi tras remover una "Num_Capas_Remover" cantidad de capas externas a dicho arreglo.
#A priori, conocer el número de capas a remover para quedarnos con el cluster principal es prácticamente imposible (no conozco un método de hacerlo analíticamente)
#Por ello es que al implementarse este algoritmo, se hacen previamente pruebas con distintos valores para la variable "Num_Capas_Remover".
#"Iteraciones" es el número de veces que se va a calcular el área de los polígonos de Voronoi tras retirar capas, cada vez para un punto arbitrario distinto.
#"Num_Capas_Remover" es el número de capas que se van a remover de los polígonos de Voronoi en cada iteración.
#"Margen_Error" es el número entero que consideramos como posible error a los enteros obtenidos al proyectar el punto arbitrario con los vectores estrella
#"Semilado_Caja" es el semilado del cuadrado centrado en el origen dentro del cual se obtiene un punto arbitrario alrededor del que se obtendrá el arrego cuasiper
#"Promedios_Distancia" es el arreglo con la separación entre las franjas cuasiperiódicas.
#"Vectores_Estrella" es el arreglo con los vectores estrella que generan la retícula cuasiperiódica deseada.
#"Arreglo_Alfas" es el arreglo con los valores numéricos de la separación respecto al origen del conjunto de rectas ortogonales a los vectores estrella en el método generalizado dual.
function arreglo_Areas_Buenos_Poligonos(Iteraciones, Num_Capas_Remover, Margen_Error, Semilado_Caja, Promedios_Distancia, Vectores_Estrella, Arreglo_Alfas)
    Arreglo_Buenas_Areas = []; #Arreglo donde meteremos las buenas áreas de los polígonos de Voronoi.
    
    for i in 1:Iteraciones
        #Generamos un punto arbitrario dentro de un cuadrado centrado en el origen de semilado "Semilado_Caja"
        APoint = Float64[]; #Arreglo donde colocaremos las coordenadas del punto arbitrario generado

        #Generamos dos parámetros que determinarán donde se localiza el punto arbitrario dentro de la caja
        x = rand();
        y = rand();

        if (x > 0.5) && (y > 0.5)
            APoint = [rand()*Semilado_Caja, rand()*Semilado_Caja];
        elseif (x > 0.5) && (y < 0.5)
            APoint = [rand()*Semilado_Caja, -rand()*Semilado_Caja];
        elseif (x < 0.5) && (y > 0.5)
            APoint = [-rand()*Semilado_Caja, rand()*Semilado_Caja];
        elseif (x < 0.5) && (y < 0.5)
            APoint = [-rand()*Semilado_Caja, -rand()*Semilado_Caja];
        end

        #Generación del arreglo cuasiperiódico con margen de error (1era entrada) y semilado del cuadrado donde colocar un pto arbitrario.
        Puntos_Duales = region_Local_Voronoi(Margen_Error, Promedios_Distancia, Vectores_Estrella, Arreglo_Alfas, APoint);

        #Obtención de los centroides de los poligonos del arreglo cuasiperiódico en: 
        #Float64(X,Y) y su diccionario Centroides <-> Vertices.
        Centroides, Diccionario_Centroides = centroides(Puntos_Duales);

        #Definimos las duplas requeridas por el algoritmo "Voronoi" de Enrique con las coordenadas de los centroides.
        sites = [(Float64(Centroides[i][1]), Float64(Centroides[i][2])) for i in 1:length(Centroides)];

        #Generamos la estructura de datos para los polígonos de Voronoi de los centroides de los polígonos del arreglo cuasipe.
        voronoi_inicial = getVoronoiDiagram(sites);

        #Generamos el diccionario "Centroides -> Indices Voronoi" de la configuración inicial
        Diccionario_Centroides_Indices = diccionario_Centroides_Indice_Voronoi(sites, voronoi_inicial);

        #Arreglo con los centroides de los polígonos a quitar
        Sitios_Centroides_Sin_Cerrar = centroides_Sin_Poligono(voronoi_inicial);
        
        #Hacemos una copia del arreglo de Centroides para poder conservar la configuración inicial y poder hacer distintas pruebas con esa misma configuración.
        Centroides_Copia = copy(Centroides);

        #Removemos capa externa tras capa externa de los polígonos de Voronoi hasta que mantengamos únicamente la capa central conformada por buenos polígonos de 
        #Voronoi. Esta parte es visual, con lo cual el número de capas a quitar puede ser variable.
        Centroides_Copia = algoritmo_Sir_Davos!(Centroides_Copia, Num_Capas_Remover);
        
        #Obtengamos el área de los buenos polígonos de Voronoi. Es decir, iteramos sobre la configuración inicial (esto es importante ya que, incluso si son buenos 
        #polígonos, al aplicarles Voronoi a estos siempre quedarán polígonos abiertos y con áreas grandes, por eso es que el cálculo de las áreas debe ser sobre 
        #los buenos polígonos pero en el arreglo de polígonos de Voronoi inicial) con todos los centroides, si el centroide es el asociado a un buen polígono se 
        #cálcula su área, caso contrario no.
        Areas_Ultima_Capa = area_Poligonos_Voronoi_Ultima_Capa(Centroides_Copia, Diccionario_Centroides_Indices, voronoi_inicial);
        
        push!(Arreglo_Buenas_Areas, Areas_Ultima_Capa);
        println("Se han realizado $(i) iteraciones.")
    end
    
    return collect(Iterators.flatten(Arreglo_Buenas_Areas))
end
