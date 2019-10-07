#Función que nos regresa el centroide de un conjunto de rombos (paralelepípedo de cuatro vértices) dados sus vértices
#"Coordenadas_X" es un arreglo con las coordenadas en X de los vértices de los polígonos. Cada 4 entradas corresponden a un mismo polígono.
#"Coordenadas_Y" es un arreglo con las coordenadas en Y de los vértices de los polígonos. Cada 4 entradas corresponden a un mismo polígono.
function centroides(Coordenadas_X, Coordenadas_Y)
    #Generamos los arreglos que contendrán las coordenadas de los centroides
    Centroide_X = [];
    Centroide_Y = [];
    Centroides = [];
    
    #Definimos un diccionario que nos servirá después para, dado un centroide, nos regrese sus vértices
    Diccionario_Centroides = Dict()
    
    #Para cada cuatro entradas, calculamos el centroide que es el promedio de los cuatro vértices
    for i in 1:4:length(Coordenadas_X)
        X = (Coordenadas_X[i]+Coordenadas_X[i+1]+Coordenadas_X[i+2]+Coordenadas_X[i+3])/4;
        Y = (Coordenadas_Y[i]+Coordenadas_Y[i+1]+Coordenadas_Y[i+2]+Coordenadas_Y[i+3])/4;
        
        push!(Centroide_X, X)
        push!(Centroide_Y, Y)
        push!(Centroides, (Float64(X),Float64(Y))) #Se pone en Float64 debido a que el algoritmo de Enrique sólo trabaja con ese formato
        
        Diccionario_Centroides[(Float64(X),Float64(Y))] = ([Coordenadas_X[i], Coordenadas_X[i+1], Coordenadas_X[i+2], Coordenadas_X[i+3]], [Coordenadas_Y[i], Coordenadas_Y[i+1], Coordenadas_Y[i+2], Coordenadas_Y[i+3]])
    end
    
    return Centroide_X, Centroide_Y, Centroides, Diccionario_Centroides
end

#Definimos una función que, dado el índice asociado a un polígono de Voronoi dentro de todos los generados, nos regresa el índice de todos los polígonos vecinos.
#"Indice" es el índice del polígono de Voronoi en el que estamos interesados.
#"Voronoi" es la estructura generada por Enrique con la función getVoronoiDiagram().
function vecinos_Voronoi(Indice, Voronoi)
    #Definimos un arreglo donde iremos guardando los centroides de los polígonos vecinos
    Vecinos_Centroides = [];
    
    #Asociado a cada polígono de Voronoi existe una clase llamada Halfedge, la cual tiene un punto de origen (un vér-
    #tice de los polígonos de Voronoi) y a partir de ahí va recorriendo los demás vértices.
    Vertice_Inicial = Voronoi.faces[Indice].outerComponent; #Dupla con las coordenadas del vértice
    Siguiente_Vertice = Voronoi.faces[Indice].outerComponent.next; #Dupla del siguiente vértice
    
    #Estudiemos el caso del gemelo asociado al inicial
    Vertice_Inicial_Gemelo = Vertice_Inicial.twin
    
    #Contador_Vertices_Recorridos = 1;
    #Si nos arroja un "nothing" es porque dicho vértice no es el que posee la información del polígono, por lo que
    #habrá que pasarnos a otro hasta que hallemos el que contenga esa información.
    while Vertice_Inicial_Gemelo.incidentFace == nothing
        #Contador_Vertices_Recorridos += 1;
        #Recorremos al siguiente vértice en el nuevo poligono de Voronoi
        Vertice_Inicial_Gemelo = Vertice_Inicial_Gemelo.next;
    end
    #println(Contador_Vertices_Recorridos)
    
    #En cuanto ubicamos el vértice que posee la información, le pedimos las coordenadas de su centroide
    Centroide = Vertice_Inicial_Gemelo.incidentFace.site;
    push!(Vecinos_Centroides, Centroide)
    
    #Si el siguiente vértice no es el mismo con el que empezamos, realizamos el algoritmo para ubicar su centroide
    #asociado distinto al polígono inicial.
    while Vertice_Inicial.origin.coordinates != Siguiente_Vertice.origin.coordinates
        Siguiente_Vertice_Gemelo = Siguiente_Vertice.twin; #Esto nos cambia el polígono al que "pertenece" el vértice
        
        #Contador_Vertices_Recorridos = 1;
        #Si nos arroja un "nothing" es porque dicho vértice no es el que posee la información del polígono, por lo que
        #habrá que pasarnos a otro hasta que hallemos el que contenga esa información.
        while Siguiente_Vertice_Gemelo.incidentFace == nothing
            #Contador_Vertices_Recorridos += 1;
            #Recorremos al siguiente vértice en el nuevo poligono de Voronoi
            Siguiente_Vertice_Gemelo = Siguiente_Vertice_Gemelo.next;
        end
        #println(Contador_Vertices_Recorridos)
        
        #En cuanto ubicamos el vértice que posee la información, le pedimos las coordenadas de su centroide
        Centroide = Siguiente_Vertice_Gemelo.incidentFace.site;
        push!(Vecinos_Centroides, Centroide)
        
        #Recorremos en uno el vértice del polígono original
        Siguiente_Vertice = Siguiente_Vertice.next
    end    
    
    return Vecinos_Centroides
end

#Función que localiza el índice del polígono de Voronoi asociado a algún centroide dado. Lo realiza por fuerza bruta.
#"Dupla_Centroide" es una dupla con las coordenadas del Centroide de interés
#"Voronoi" es la estructura generada por Enrique con la función getVoronoiDiagram().
function indice_Voronoi_Centroide(Dupla_Centroide, Voronoi)
    #Comencemos tomando el primero de los posibles polígonos de Voronoi
    i = 1;
    
    #Mientras en el índice "i" no encontremos la dupla, pasamos a otro índice
    while Dupla_Centroide != Voronoi.faces[i].site
        i += 1;
    end
    
    #Cuando encontremos el índice, lo regresamos al usuario.
    return i
end

#Función que nos regresa, dado un conjunto de centroides, los vértices de los polígonos que los generan. Emplea un diccionario para ello.
#"Centroides" es un arreglo con las duplas de los centroides de interés
#"Diccionario_Centroides" es un diccionario que relaciona los centroides con los vértices del polígono que los genera.
function Centroides_a_Vertices(Centroides, Diccionario_Centroides)
    Coordenadas_X = [];
    Coordenadas_Y = [];

    for i in Centroides
        push!(Coordenadas_X, Diccionario_Centroides[i][1])
        push!(Coordenadas_Y, Diccionario_Centroides[i][2])
    end
    
    return collect(Iterators.flatten(Coordenadas_X)), collect(Iterators.flatten(Coordenadas_Y))
end

#Función que itera sobre los distintos polígonos hasta encontrar el que contiene al punto "Punto".
#"Punto" es el punto que queremos checar si está dentro o fuera. Es un arreglo con dos entradas.
#"Poligonos" es un arreglo de arreglos, cada uno con los segmentos que conforman al polígono en cuestión.
function encontrar_Poligono_Voronoi(Punto, Poligonos)
    #Iteremos sobre los posibles polígonos para encontrar el que contenga al punto
    for i in 1:length(Poligonos)
        #Si el polígono i-ésimo contiene al punto, regresa la info de ese polígono
        if dentro(Poligonos[i], Punto)
            return i
        end
        
    end
    #Si no encuentra polígono, que nos mande impresión con dicha información
    println("Error: No hay polígono que contenga al punto")
end

#Función que genera una vecindad del arreglo cuasiperiódico alrededor de un punto dado (o arbitrario dentro de un cuadrado)
#"N" es el margen de error que vamos a permitir en los posibles números enteros generados por la proyección del punto sobre los vectores estrella.
#"Cota" es el Semilado de la caja centrada en el origen dentro de la cual se va a calcular el punto arbitrario si el usuario no da uno.
function region_Local_Voronoi(N, Cota, Punto = true)
    #Si el usuario no proporciona un punto, generamos uno
    if Punto == true
        x = rand();
        y = rand();

        if (x > 0.5) && (y > 0.5)
            Punto = [rand()*Cota, rand()*Cota];
        elseif (x > 0.5) && (y < 0.5)
            Punto = [rand()*Cota, -rand()*Cota];
        elseif (x < 0.5) && (y > 0.5)
            Punto = [-rand()*Cota, rand()*Cota];
        elseif (x < 0.5) && (y < 0.5)
            Punto = [-rand()*Cota, -rand()*Cota];
        end
    end
    
    #Paso 1: Dado el punto proyectamos para obtener los enteros aproximados que pueden dar lugar al polígono contenedor.
    Proyecciones = proyecciones_Pto_Direccion_Franjas(Punto, Promedios_Distancia, Vectores_Estrella);

    #Paso 2: A partir de los valores enteros aproximados, generamos la vecindad del arreglo cuasiperiódico que contenga al
    #punto.
    Puntos_Duales, Informacion_Duales = generador_Vecindades_Vertices(Proyecciones, Vectores_Estrella, Arreglo_Alfas, N);

    #Paso 3: Separamos en coordenadas X y coordenadas Y los vértices de los polígonos que conforman la vecindad del punto
    #anterior.
    Coordenadas_X, Coordenadas_Y = separacion_Arreglo_de_Arreglos_2D(Puntos_Duales);
    
    return Coordenadas_X, Coordenadas_Y, Punto
end

#Función que busca el polígono del arreglo cuasiperiódico que contiene al punto de interés empleando polígonos de Voronoi.
#"Coordenadas_X" es un arreglo con las coordenadas en X de los vértices de los polígonos. Cada 4 entradas corresponden a un mismo polígono.
#"Coordenadas_Y" es un arreglo con las coordenadas en Y de los vértices de los polígonos. Cada 4 entradas corresponden a un mismo polígono.
#"Punto" es un arreglo dos dimensional con las coordenadas de un punto en el espacio 2D.
function poligono_Contenedor_Voronoi(Coordenadas_X, Coordenadas_Y, Punto)
    #Paso 4: Obtenemos el conjunto de centroides de nuestros futuros polígonos.
    Centroides_X, Centroides_Y, Centroides, Diccionario_Centroides = centroides(Coordenadas_X, Coordenadas_Y);

    #Paso 5: Agregamos al conjunto de Centroides_X y Centroides_Y las respectivas coordenadas del punto arbitrario
    push!(Centroides_X, Punto[1])
    push!(Centroides_Y, Punto[2])

    #Definimos las duplas con las coordenadas de los centroides
    sites = [(Float64(Centroides_X[i]), Float64(Centroides_Y[i])) for i in 1:length(Centroides_X)]

    voronoi = getVoronoiDiagram(sites);
    
    #Obtenemos el índice del polígono correspondiente a nuestro punto arbitrario
    Indice = indice_Voronoi_Centroide((Float64(Centroides_X[end]), Float64(Centroides_Y[end])), voronoi);

    #Obtenemos los vecinos al polígono de Voronoi asociado a nuestro punto arbitrario
    Vecinos_Centroides = vecinos_Voronoi(Indice, voronoi);

    #Recuperamos la información de los vértices de los polígonos de Voronoi candidatos
    X,Y = Centroides_a_Vertices(Vecinos_Centroides, Diccionario_Centroides);
    
    #Transformamos los vertices con coordenadas (X,Y) a estructura de polígono
    Poligonos = obtener_Segmentos_Vertices(X,Y);
    
    #Iteramos sobre los polígonos candidatos para obtener el polígono contenedor
    Indice_Poligono_Contenedor = encontrar_Poligono_Voronoi(Punto, Poligonos);
    
    plot()
    plot([Poligonos[Indice_Poligono_Contenedor][1].inicio[1], Poligonos[Indice_Poligono_Contenedor][2].inicio[1], Poligonos[Indice_Poligono_Contenedor][3].inicio[1], Poligonos[Indice_Poligono_Contenedor][4].inicio[1], Poligonos[Indice_Poligono_Contenedor][1].inicio[1]], [Poligonos[Indice_Poligono_Contenedor][1].inicio[2], Poligonos[Indice_Poligono_Contenedor][2].inicio[2], Poligonos[Indice_Poligono_Contenedor][3].inicio[2], Poligonos[Indice_Poligono_Contenedor][4].inicio[2], Poligonos[Indice_Poligono_Contenedor][1].inicio[2]], markersize = 0.2, key = false, aspect_ratio=:equal, grid = false, color =:black)
    scatter!([Punto[1]], [Punto[2]], markersize = 5, markeralpha = 0.5, markerstrokewidth = 0, markercolor = :red)
end
