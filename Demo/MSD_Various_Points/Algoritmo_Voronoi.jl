#Función que nos regresa el centroide de un conjunto de rombos (paralelepípedo de cuatro vértices) dados sus vértices.
#"Vertices" es un arreglo con las coordenadas (X,Y) de los vértices de los polígonos. Cada 4 entradas corresponden a un mismo polígono.
function centroides(Vertices)
    #Paso 1: Generamos el arreglo que contendrá las coordenadas de los centroides.
    Centroides = [];
    
    #Paso 2: Definimos un diccionario que nos servirá después para, dado un centroide, nos regrese sus vértices
    Diccionario_Centroides = Dict();
    
    #Paso 3: Para cada cuatro entradas, calculamos el centroide que es el promedio de los cuatro vértices
    for i in 1:4:length(Vertices)
        X = (Vertices[i][1] + Vertices[i+1][1] + Vertices[i+2][1] + Vertices[i+3][1])/4;
        Y = (Vertices[i][2] + Vertices[i+1][2] + Vertices[i+2][2] + Vertices[i+3][2])/4;
        
        push!(Centroides, (Float64(X),Float64(Y))); #Se pone en Float64 debido a que el algoritmo de Enrique (Voronoi) sólo trabaja con ese formato
        
        Diccionario_Centroides[(Float64(X),Float64(Y))] = ([Vertices[i][1], Vertices[i+1][1], Vertices[i+2][1], Vertices[i+3][1]], [Vertices[i][2], Vertices[i+1][2], Vertices[i+2][2], Vertices[i+3][2]]);
    end
    
    return Centroides, Diccionario_Centroides
end

#Definimos una función que, dado el índice asociado a un polígono de Voronoi dentro de todos los generados, nos regresa el índice de todos los polígonos vecinos.
#"Indice" es el índice del polígono de Voronoi en el que estamos interesados.
#"Voronoi" es la estructura generada por Enrique con la función getVoronoiDiagram().
function vecinos_Voronoi(Indice, Voronoi)
    #Definimos un arreglo donde iremos guardando los centroides de los polígonos vecinos
    Vecinos_Centroides = [];
    
    #Partimos de un lado del polígono de Voronoi que es de nuestro interés
    Lado_Poligono_Voronoi = Voronoi.faces[Indice].outerComponent;

    #Iniciamos el proceso while para recorrer todos los lados del polígono de Voronoi y en cada uno hallar el polígono vecino
    while true
        #Encontramos el vecino asociado al lado que estamos considerando
        Coordenadas_Vecino = Lado_Poligono_Voronoi.twin.incidentFace.site;
        push!(Vecinos_Centroides, Coordenadas_Vecino);

        #Recorremos al siguiente lado del polígono
        Lado_Poligono_Voronoi = Lado_Poligono_Voronoi.next;

        #Checamos si hemos ya concluido de revisar todos los lados del polígono de Voronoi
        if Lado_Poligono_Voronoi == Voronoi.faces[Indice].outerComponent
            break
        end
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
#"Centroides" es un arreglo con las duplas (X,Y) de los centroides de interés
#"Diccionario_Centroides" es un diccionario que relaciona los centroides con los vértices del polígono que los genera.
function centroides_A_Vertices(Centroides, Diccionario_Centroides)
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

#Función que genera una vecindad de la retícula cuasiperiódica alrededor de un punto dado
#"N" es el margen de error asociado a los números enteros generados por la proyección del punto sobre los vectores estrella.
#"Promedios_Distancia" es el arreglo con la separación entre las franjas cuasiperiódicas.
#"Vectores_Estrella" es el arreglo con los vectores estrella que generan la retícula cuasiperiódica deseada.
#"Arreglo_Alfas" es el arreglo con los valores numéricos de la separación respecto al origen del conjunto de rectas ortogonales a los vectores estrella en el método generalizado dual.
#"Punto" es el punto alrededor de donde se va a generar la vecindad.
function region_Local_Voronoi(N, Promedios_Distancia, Vectores_Estrella, Arreglo_Alfas, Punto)
    #Paso 1: Dado el Punto proyectamos este con los Vectores Estrella para obtener los enteros aproximados asociados al polígono contenedor.
    #Salida: [n1,n2,n3,...]
    Proyecciones = proyecciones_Pto_Direccion_Franjas(Punto, Promedios_Distancia, Vectores_Estrella);

    #Paso 2: A partir de los valores enteros aproximados, generamos la vecindad del arreglo cuasiperiódico que contenga al punto.
    #Salida: [[X,Y]]
    Puntos_Duales = generador_Vecindades_Vertices(Proyecciones, Vectores_Estrella, Arreglo_Alfas, N);
    
    return Puntos_Duales
end

#Función que busca el polígono del arreglo cuasiperiódico que contiene al punto de interés empleando polígonos de Voronoi.
#"Coordenadas_Vertices" es un arreglo con las coordenadas en (X,Y) de los vértices de los polígonos. Cada 4 entradas corresponden a un mismo polígono.
#"Punto" es un arreglo dos dimensional con las coordenadas de un punto en el espacio 2D.
function poligono_Contenedor_Voronoi(Coordenadas_Vertices, Punto)
    #Paso 4: Obtenemos el conjunto de centroides de nuestros futuros polígonos.
    Centroides, Diccionario_Centroides = centroides(Coordenadas_Vertices);

    #Paso 5: Agregamos al conjunto de Centroides las coordenadas del punto arbitrario
    push!(Centroides, (Float64(Punto[1]), Float64(Punto[2])))

    #Definimos las duplas con las coordenadas de los centroides
    sites = [(Float64(Centroides[i][1]), Float64(Centroides[i][2])) for i in 1:length(Centroides)]

    voronoi = getVoronoiDiagram(sites);
    
    #Obtenemos el índice del polígono correspondiente a nuestro punto arbitrario
    Indice = indice_Voronoi_Centroide((Float64(Punto[1]), Float64(Punto[2])), voronoi);

    #Obtenemos los vecinos al polígono de Voronoi asociado a nuestro punto arbitrario
    Vecinos_Centroides = vecinos_Voronoi(Indice, voronoi);

    #Recuperamos la información de los vértices de los polígonos de Voronoi candidatos
    X,Y = centroides_A_Vertices(Vecinos_Centroides, Diccionario_Centroides);
    
    #Transformamos los vertices con coordenadas (X,Y) a estructura de polígono
    Poligonos = obtener_Segmentos_Vertices(X,Y);
    
    #Iteramos sobre los polígonos candidatos para obtener el polígono contenedor
    Indice_Poligono_Contenedor = encontrar_Poligono_Voronoi(Punto, Poligonos);
    
    return Poligonos[Indice_Poligono_Contenedor]
end
