#Función que busca los vertices del arreglo cuasiperiódico más cercanos al punto de interés empleando polígonos de Voronoi.
#"Punto" es un arreglo dos dimensional con las coordenadas de un punto en el espacio 2D.
#"Vertices_Unicos" es un arreglo con las coordenadas (X,Y) de los vértices sin repetir.
function vertices_Cercanos_Voronoi(Punto, Vertices_Unicos)
    #Definimos las duplas con las coordenadas de los centroides
    #sites = [(Float64(Vertices_Unicos[i][1]), Float64(Vertices_Unicos[i][2])) for i in 1:length(Vertices_Unicos)];

    #Agregamos el punto de interés a la lista de vertices
    push!(Vertices_Unicos, (Float64(Punto[1]), Float64(Punto[2])))

    #Generamos la estructura de Voronoi asociada a las coordenadas de los vertices de los polígonos en el arreglo cuasiperiódico
    voronoi = getVoronoiDiagram(Vertices_Unicos);
    
    #Obtenemos el índice del polígono correspondiente a nuestro punto arbitrario
    Indice = indice_Voronoi_Centroide(Vertices_Unicos[end], voronoi);

    #Obtenemos los vecinos al polígono de Voronoi asociado a nuestro punto arbitrario
    Vecinos_Vertices = vecinos_Voronoi(Indice, voronoi);

    #Eliminamos el punto que agregamos a nuestra lista con los vertices del polígono de voronoi de los centros de los obstáculos
    pop!(Vertices_Unicos);

    return Vecinos_Vertices
end

#Función que dado un conjunto de posibles vértices a ser el vértice más cercano al pto de interes
#nos regresa el vertice más cercano.
function vertice_Cercano(Candidatos, Punto)
    N = length(Candidatos);
    Distancias = zeros(BigFloat,N); #Arreglo que tendrá las distancias euclidianas del pto con los vertices
    Diccionario_Norma_Indice = Dict(); #Diccionario que relaciona "Norma -> Indice"
    for i in 1:N
        Norma = norm([Candidatos[i][1], Candidatos[i][2]] - Punto); #Calculamos la separación del pto con el i-ésimo vértice
        Distancias[i] += Norma
        Diccionario_Norma_Indice[Norma] = i; #Asociamos la norma con el índice
    end
    
    return Candidatos[Diccionario_Norma_Indice[minimum(Distancias)]]
end

#Posicion_Inicial: Coordenadas (X,Y) de la posición de la partícula antes de la colisión
#Velocidad_Inicial: Coordenadas (Vx, Vy) de la velocidad de la partícula antes de la colisión
#Radio_Obstaculo: Radio de los obstáculos
#Posicion_Obstaculo: Coordenadas (Ox, Oy) de la posición del obstáculo.
#funcion_Posicion: Función que determina la posición tras la interacción Partícula-Obstáculo
#funcion_Velocidad: Función que determina la velocidad tras la interacción Partícula-Obstáculo
function colision_Obstaculo(Posicion_Inicial, Velocidad_Inicial, Radio_Obstaculo, Posicion_Obstaculo, funcion_Posicion::Function, funcion_Velocidad::Function)
    Posicion_Colision, Tiempo_Colision = funcion_Posicion(Posicion_Inicial, Posicion_Obstaculo, Radio_Obstaculo, Velocidad_Inicial);
    
    if Tiempo_Colision < Inf
        Velocidad_Final = funcion_Velocidad(Posicion_Colision, Posicion_Obstaculo, Velocidad_Inicial);
        
        return Posicion_Colision, Velocidad_Final, Tiempo_Colision
    else
        return Posicion_Inicial, Velocidad_Inicial, Tiempo_Colision
    end
end

#Función que, dados los puntos A y B del segmento AB y los puntos C y D del segmento CD nos 
#regresa el punto de intersección de ambas rectas
#A,B,C,D: Arreglos de la forma (x,y) 
function interseccion_rectas(A, B, C, D)
    #Obtengamos los coeficientes "a1", "b1" y "c1" asociados a la línea AB:
    #a1x + b1y = c1
    a1 = B[2] - A[2];
    b1 = A[1] - B[1];
    c1 = a1*A[1] + b1*A[2];
    
    #Obtengamos los coeficientes "a2", "b2" y "c2" asociados a la línea CD
    #a2x + b2y = c2
    a2 = D[2] - C[2];
    b2 = C[1] - D[1];
    c2 = a2*C[1] + b2*C[2];
    
    Determinante = a1*b2 - a2*b1;
    
    #Chequemos si las líneas son paralelas o no
    if Determinante == 0 #Líneas paralelas
        return [Inf, Inf] #Punto de colisión
    else
        return [(b2*c1 - b1*c2)/Determinante, (a1*c2 - a2*c1)/Determinante] #Pto Colision
    end
end

#Función que, dada la posición inicial y la posición final de una partícula, con velocidad V, 
#regresa el tiempo que tarda en desplazarse.
#Posicion_Inicial: Arreglo [X,Y] con las coordenadas de la posición inicial de la partícula.
#Posicion_Inicial: Arreglo [X,Y] con las coordenadas de la posición final de la partícula.
#Rapidez: Escalar asociado a la rapidez de la partícula (Norma de su velocidad).
function tiempo_Vuelo(Posicion_Inicial, Posicion_Final, Rapidez)
    #Obtenemos la distancia euclidiana entre los dos puntos
    Distancia = norm(Posicion_Final - Posicion_Inicial);
    
    #El tiempo de vuelo está determinado por t = Distancia/Rapidez
    return Distancia/Rapidez;
end

#Función que dado el índice del polígono de Voronoi que contiene a la partícula, su posición y
#su velocidad, nos regresa el tiempo de vuelo en salir del poligono, así como el índice del
#polígono al cual llega.
#Posicion_Particula: Posición de la partícula.
#Velocidad_Particula: Velocidad de la partícula.
#Voronoi: Estructura de Voronoi de los vértices.
#Indice_Poligono_Contenedor: Indice del polígono que contiene a la partícula.
#Diccionario_Vertice_Indice: Diccionario "Vertice (X,Y) -> Indice_Voronoi"
function cambio_Poligono(Posicion_Particula, Velocidad_Particula, Voronoi, Indice_Poligono_Contenedor, Extremos_Entrada, Diccionario_Vertices_Indice, Diccionario_Vertices_Cuadrado)
    if Diccionario_Vertices_Cuadrado[Voronoi.faces[Indice_Poligono_Contenedor].site] #Si la condición se cumple, estamos en un vértice dentro del parche cuadrado seguro
        Arreglo_Tiempos_Vuelo = Float64[]; #Arreglo donde iremos guardando los tiempos de vuelo
        Diccionario_Tiempos_Vuelo_Indice = Dict(); #Diccionario que relaciona "tiempo -> Indice_Voronoi"
        Diccionario_Tiempos_Vuelo_Extremos = Dict(); #Diccionario que relaciona "tiempo -> Vertices del lado por el que entra la partícula"

        #Determinamos el punto final de la partícula tras moverse un tiempo t = 1.
        Posicion_Particula_2 = Posicion_Particula + Velocidad_Particula;
        
        #Partimos de un lado del polígono de Voronoi que es de nuestro interés
        Lado_Poligono_Voronoi = Voronoi.faces[Indice_Poligono_Contenedor].outerComponent;
        
        #Iniciamos el proceso while para recorrer todos los lados del polígono de Voronoi y en cada uno hallar el polígono vecino
        while true
            #Encontramos el indice del vecino asociado al lado que estamos considerando
            Indice_Vecino = Diccionario_Vertices_Indice[Lado_Poligono_Voronoi.twin.incidentFace.site];
            #println("El vecino es $(Indice_Vecino)")
            #Obtenemos los puntos que conforman el lado del polígono de Voronoi
            Punto1 = Lado_Poligono_Voronoi.origin.coordinates;
            Punto2 = Lado_Poligono_Voronoi.next.origin.coordinates;

            if [Punto1, Punto2] == Extremos_Entrada || [Punto2, Punto1] == Extremos_Entrada #Si es la puerta por donde entró, no la consideramos
                nothing
            elseif (Punto1[1] - Punto2[1])^2 + (Punto1[2] - Punto2[2])^2 < 1e-6 #El lado es muy pequeño, es un falso lado del polígono Voronoi
            	nothing
            else            
                #Calculamos el punto de intersección entre la partícula y la línea descrita por ese lado
                Punto_Colision = interseccion_rectas(Posicion_Particula, Posicion_Particula_2, Punto1, Punto2);
                
                if dot(Punto_Colision - Posicion_Particula, Velocidad_Particula) < 0.0 #La colisión se da a tiempos negativos
                    push!(Arreglo_Tiempos_Vuelo, Inf); #Agregamos el tiempo de vuelo al arreglo
                else
                    Tiempo_Vuelo = Float64(tiempo_Vuelo(Posicion_Particula, Punto_Colision, 1)); #Tiempo que tardaría la partícula en colisionar
                    push!(Arreglo_Tiempos_Vuelo, Tiempo_Vuelo); #Añadimos el tiempo de vuelo al arreglo
                    Diccionario_Tiempos_Vuelo_Indice[Tiempo_Vuelo] = Indice_Vecino; #Relacionamos el tiempo de vuelo con el índice del vecino
                    Diccionario_Tiempos_Vuelo_Extremos[Tiempo_Vuelo] = [Punto1, Punto2]; #Relacionamos el tiempo de vuelo con la puerta de entrada a la partícula
                end
            end
            
            #Recorremos al siguiente lado del polígono
            Lado_Poligono_Voronoi = Lado_Poligono_Voronoi.next;

            #Checamos si hemos ya concluido de revisar todos los lados del polígono de Voronoi
            if Lado_Poligono_Voronoi === Voronoi.faces[Indice_Poligono_Contenedor].outerComponent
                break
            end
        end
        
        #Obtenemos el mínimo de los tiempos de vuelo
        Tiempo_Vuelo_Minimo = minimum(Arreglo_Tiempos_Vuelo);
        
        if Tiempo_Vuelo_Minimo == Inf
            println("Contiene $(Indice_Poligono_Contenedor) y Falla")
            return Indice_Poligono_Contenedor, 0.0, [Inf, Inf], true
        end
        
        #Obtenemos el índice del polígono al cual entrará la partícula
        Indice_Poligono_Entrante = Diccionario_Tiempos_Vuelo_Indice[Tiempo_Vuelo_Minimo]
        Extremos_Entrada = Diccionario_Tiempos_Vuelo_Extremos[Tiempo_Vuelo_Minimo]
        #println("Contiene $(Indice_Poligono_Contenedor)")
        #println("Entra a $(Indice_Poligono_Entrante)")
        return Indice_Poligono_Entrante, Tiempo_Vuelo_Minimo, Extremos_Entrada, false
    else
        return Indice_Poligono_Contenedor, 0.0, Extremos_Entrada, true
    end
end

#Función que dadas las condiciones iniciales de la partícula y un tiempo de vuelo dado, nos regresa la trayectoria de 
#la partícula.
#Posicion: Posicion inicial de la partícula
#Velocidad: Velocidad inicial de la partícula
#Tiempo_Vuelo: Tiempo de vuelo deseado que recorra la partícula
#Indice_Obstaculo_Cercano: Indice en la lista de Voronoi del obstáculo más cercano a la partícula
#Radio_Obstaculo: Radio de los obstáculos
#Voronoi: Estructura de voronoi de los vertices
#Diccionario_Vertices_Indices: Diccionario "Vertices (X,Y) -> Indices_Voronoi"
#Diccionario_Vertices_Cuadrado: Diccionario "Vertices (X,Y) -> Booleano" que nos indica si el polígono es un polígono dentro de algún parche cuadrado o no
#Lado_Poligono_Bloqueado: Información sobre el lado del polígono por el que la partícula entró.
function vuelo_Una_Particula_Grafico(Posicion, Velocidad, Tiempo_Vuelo, Indice_Obstaculo_Cercano, Radio_Obstaculo, Voronoi, Diccionario_Vertices_Indices, Diccionario_Vertices_Cuadrado, Lado_Poligono_Bloqueado)
    #plot(); #VISUALIZATION PURPOSE ONLY
    while Tiempo_Vuelo > 0.0
        #Guardamos las coordenadas X,Y del obstáculo más cercano a nuestra partícula
        Posicion_Obstaculo = [Voronoi.faces[Indice_Obstaculo_Cercano].site[1], Voronoi.faces[Indice_Obstaculo_Cercano].site[2]];
        
        #=
        x(t) = Posicion_Obstaculo[1] + Obstacle_Radius*cos(t);
        y(t) = Posicion_Obstaculo[2] + Obstacle_Radius*sin(t);
        plot!(x, y, 0, 2π, leg=false, color="blue")
        =#
        
        #Obtenemos la posición y la velocidad tras la colisión (incluso si no hay colisión).
        #Si no hay una colisión, entonces el tiempo de colisión es infinito.
        #NOTA: Esta parte del algoritmo presupone que, si hay colisión con obstáculo, esta se da el obstáculo
        #centrado en el polígono de Voronoi que contiene a la partícula: esto es cierto ÚNICAMENTE PARA RADIOS PEQUEÑOS
        Posicion_Colision, Velocidad_Colision, Tiempo_Vuelo_Colision = colision_Obstaculo(Posicion, Velocidad, Radio_Obstaculo, Posicion_Obstaculo, colision_Esfera_Dura, velocidad_Esfera_Dura);
        
        #Obtengamos el índice del polígono que contendrá a la partícula y el tiempo que le tomará cambiar de polígono
        Indice_Proximo_Obstaculo, Tiempo_Vuelo_Cambio, Lado_Poligono_Bloqueado_Cambio, Frontier_Polygon = cambio_Poligono(Posicion, Velocidad, Voronoi, Indice_Obstaculo_Cercano, Lado_Poligono_Bloqueado, Diccionario_Vertices_Indices, Diccionario_Vertices_Cuadrado);
        #@show Indice_Proximo_Obstaculo
        #Chequemos si el polígono contenedor de la partícula está en los bordes exteriores del parche cuadrado seguro.
        #Si lo está, no calculamos más, nos salimos de la función para generar otro parche centrado en esta posición
        if Frontier_Polygon
            return Posicion, Velocidad, Tiempo_Vuelo, Indice_Obstaculo_Cercano, Lado_Poligono_Bloqueado
        end
        
        if Tiempo_Vuelo_Colision < Tiempo_Vuelo_Cambio #La partícula colisiona antes de moverse de polígono
            if (Tiempo_Vuelo - Tiempo_Vuelo_Colision) > 0.0 #Hay tiempo para la colisión completa
                #plot!([Posicion[1], Posicion_Colision[1]], [Posicion[2], Posicion_Colision[2]], color = "black")
                #Actualizamos las variables
                Tiempo_Vuelo -= Tiempo_Vuelo_Colision; #Se reduce al tiempo de vuelo restante el tiempo en colisionar
                Posicion = Posicion_Colision; #La nueva posición de la partícula es la posición de la colisión
                Velocidad = Velocidad_Colision; #La nueva velocidad de la partícula es la velocidad tras la colisión
                #scatter!([Posicion[1]],[Posicion[2]], markersize = 0.5); #VISUALIZATION PURPOSE ONLY
                
                #Si hay colisión permitimos que la partícula salga por donde entró
                Lado_Poligono_Bloqueado = [Inf, Inf];
            else #Se agota el tiempo de vuelo antes de la colisión completa
                #plot!([Posicion[1], (Posicion + Tiempo_Vuelo*Velocidad)[1]], [Posicion[2], (Posicion + Tiempo_Vuelo*Velocidad)[2]], color = "black")
                #La posición final de la partícula será la última posición válida más el tiempo restante de vuelo por la última velocidad
                Posicion = Posicion + Tiempo_Vuelo*Velocidad;
                #scatter!([Posicion[1]],[Posicion[2]], markersize = 0.5); #VISUALIZATION PURPOSE ONLY
                return Posicion, Velocidad, 0.0, Indice_Obstaculo_Cercano, Lado_Poligono_Bloqueado
            end
        else #La partícula se mueve de polígono antes de colisionar con el obstáculo de su polígono contenedor
            if (Tiempo_Vuelo - Tiempo_Vuelo_Cambio) > 0.0 #Hay tiempo para que la partícula salga del polígono  
                #plot!([Posicion[1], (Posicion + (Tiempo_Vuelo_Cambio)*Velocidad)[1]], [Posicion[2], (Posicion + (Tiempo_Vuelo_Cambio)*Velocidad)[2]], color = "black")
                Tiempo_Vuelo -= (Tiempo_Vuelo_Cambio); #Actualizamos el tiempo de vuelo restante
                Posicion = Posicion + (Tiempo_Vuelo_Cambio)*Velocidad; #Avancemos la posición de la partícula al nuevo polígono
                Indice_Obstaculo_Cercano = Indice_Proximo_Obstaculo; #Actualizamos el índice del polígono contenedor
                Lado_Poligono_Bloqueado = Lado_Poligono_Bloqueado_Cambio; #El lado del polígono bloqueado cambia
                #scatter!([Posicion[1]],[Posicion[2]], markersize = 0.5); #VISUALIZATION PURPOSE ONLY
            else #No hay tiempo de vuelo restante para salir del polígono
                #plot!([Posicion[1], (Posicion + Tiempo_Vuelo*Velocidad)[1]], [Posicion[2], (Posicion + Tiempo_Vuelo*Velocidad)[2]], color = "black")
                #La posición final de la partícula será la última posición válida más el tiempo restante de vuelo por la última velocidad
                Posicion = Posicion + Tiempo_Vuelo*Velocidad;
                #scatter!([Posicion[1]],[Posicion[2]], markersize = 0.5); #VISUALIZATION PURPOSE ONLY
                return Posicion, Velocidad, 0.0, Indice_Obstaculo_Cercano, Lado_Poligono_Bloqueado
            end
        end
    end
end

#Function that, given a initial position in the Quasiperiodic Array, calculates the MSD of the system iterating over
#a given number of flies (each one conformed as a set of short-time flies), each one with a random velocity.
#Patch_Information: An array that contains information about the semi-circular patches of the Quasiperiodic Array
#Reduction_Factor: The factor with which we multiple the Average Radius to generate the Safe Radius
#Average_Distance_Stripes: Array with the average distance between stripes
#Star_Vectors: Array wich will contain the Star Vectors
#Alphas_Array: Array of the alphas constant
#APoint_Initial: A fixed position in the Quasiperiodic array that will be the Initial Position of all the particles,
#each one with a different random velocity
#Number_Velocities: It's the number of different random velocities with which the user desires to obtain the MSD
#Number_Flights: The number of short flies that conform a long one
#Short_Fly_Time: The flying time in which the user desires to divide the long-fly time
#Obstacle_Radius: The radius of the obstacles centered in the vertices of the quasiperiodic array
#MSD_Array: An array that have the MSD as a function of the fly time
function MSD_Varying_Velocities(Patch_Information, Reduction_Factor, Average_Distance_Stripes, Star_Vectors, Alphas_Array, APoint_Initial, Number_Velocities, Number_Flights, Short_Fly_Time, Obstacle_Radius, MSD_Array)
    SL = Patch_Information[2]; #Size of a half side of the square in which the algorithm generate a random point inside it
    Average_Radius = Patch_Information[end]; #The value of the previously calculated "Average Radius"
    Safe_Radius = Reduction_Factor*Average_Radius; #The Safe Radius value

    #STEP 1: Generate a square patch of the Quasiperiodic Array around the Initial Position
    X,Y,APoint_Initial,Centroids,Dictionary_Centroids = parche_Cuadrado(Patch_Information,SL,Reduction_Factor,Average_Distance_Stripes,Star_Vectors,Alphas_Array,APoint_Initial);

    #STEP 2: Get the structure of the Voronoi's Polygons of the Vertices that conform the Quasiperiodic Array
    #Obtain the vertices of all the polygons as duples.
    Sites_Vertices = [(Float64(X[i]), Float64(Y[i])) for i in 1:length(X)];
    unique!(Sites_Vertices); #Eliminate all the copies of a vertex

    #Let's generate a Dictionary with the coordinates of the vertices that lay inside the square cluster.
    #The Dictionary relates "Vertex (X,Y) -> "true" or "false" depending if the vertex is inside the square or not
    Dictionary_Vertices_Inside_Square = Dict();
    for i in Sites_Vertices
        if (APoint_Initial[1]-sqrt(2)*Safe_Radius <= i[1] <= APoint_Initial[1]+sqrt(2)*Safe_Radius) && (APoint_Initial[2]-sqrt(2)*Safe_Radius <= i[2] <= APoint_Initial[2]+sqrt(2)*Safe_Radius)
            Dictionary_Vertices_Inside_Square[i] = true;
        else
            Dictionary_Vertices_Inside_Square[i] = false;
        end
    end

    #Get the Voronoi structure with only the non-repeated vertices
    Voronoi_Vertices = getVoronoiDiagram(Sites_Vertices);

    #Let's get a dictionary that relates "Vertex (X,Y) -> Index Voronoi's Polygon"
    Dictionary_Vertex_Index = diccionario_Centroides_Indice_Voronoi(Sites_Vertices, Voronoi_Vertices);
    
    #STEP 2.5: Defining some arrays and variables
    Array_Position_Patches = [APoint_Initial]; #Array that will contain the position of all the centers of the patches
    
    #=
    Number_Flights = Int64(Long_Fly_Time/Short_Fly_Time); #The number of short flies that conform a long one
    
    #ARRAYS USED IN THE MSD ALGORITHM
    Accumulated_Fly_Time = zeros(Number_Flights); #An array with the accumulated fly time
    for i in 1:Number_Flights
        Accumulated_Fly_Time[i] = i*Short_Fly_Time;
    end
    MSD_Array = zeros(Number_Flights); #An array that will contain the MSD of the system
    =#
    
    #STEP 3: Iterates over the desired number of velocities with the fixed initial position
    for α in 1:Number_Velocities
        println("The algorithm is calculating the particle $(α).")

        #STEP 4: Get a Initial Velocity
        θ = 2*π*rand(); #Get a random angle between 0 and 2*pi
        AVelocity_Initial = [cos(θ), sin(θ)];

        #STEP 5: Find the container polygon of the Initial Position
        APoint = copy(APoint_Initial); #A copy of Initial Position in order to preserve the original one

        #Let's obtain the vertices coordinates of the possible nearest vertex to the point
        Nearest_Vertex_Candidates = vertices_Cercanos_Voronoi(APoint, Sites_Vertices);

        #Get the index of the nearest vertex from the list of candidates
        Nearest_Vertex = vertice_Cercano(Nearest_Vertex_Candidates, APoint);

        #Get the index in the Vertex's Voronoi of the nearest vertex
        Nearest_Vertex_Index_Voronoi = Dictionary_Vertex_Index[Nearest_Vertex];

        #At the begin, all the sides of the container polygons are valid
        Forbidden_Side_Polygon = [Inf, Inf];

        #STEP 6: Let's the particle fly
        for β in 1:Number_Flights
            #Make a little fly
            APoint,AVelocity_Initial,Remaining_Fly_Time,Nearest_Vertex_Index_Voronoi,Forbidden_Side_Polygon = vuelo_Una_Particula_Grafico(APoint,AVelocity_Initial,Short_Fly_Time,Nearest_Vertex_Index_Voronoi,Obstacle_Radius,Voronoi_Vertices,Dictionary_Vertex_Index,Dictionary_Vertices_Inside_Square,Forbidden_Side_Polygon);
            Nearest_Vertex_End_Patch = Voronoi_Vertices.faces[Nearest_Vertex_Index_Voronoi].site; #X,Y coordinates of the obstacle that correspond to the container polygon

            if Remaining_Fly_Time > 0.0
                #STEP 7: Generate a new square patch if the previous one wasn't enough
                X,Y,APoint,Centroids,Dictionary_Centroids = parche_Cuadrado(Patch_Information,SL,Reduction_Factor,Average_Distance_Stripes,Star_Vectors,Alphas_Array,APoint);
                push!(Array_Position_Patches, APoint);

                #STEP 8: Add the new patch vertices information to the previous one patch vertices information
                Sites_Vertices_New_Patch = [(Float64(X[i]), Float64(Y[i])) for i in 1:length(X)];
                unique!(Sites_Vertices_New_Patch); #Eliminate all the copies of a vertex
                Sites_Vertices = vcat(Sites_Vertices, Sites_Vertices_New_Patch); #A new array with the vertices of the bigger patch
                unique!(Sites_Vertices); #Eliminate all the copies of a vertex

                #Let's generate a Dictionary with the coordinates of the vertices that lay inside the squares clusters.
                #The Dictionary relates "Vertex (X,Y) -> "true" or "false" depending if the vertex is inside the square or not
                Dictionary_Vertices_Inside_Square = Dict();
                for k in Sites_Vertices
                    for l in Array_Position_Patches
                        if (l[1]-sqrt(2)*Safe_Radius <= k[1] <= l[1]+sqrt(2)*Safe_Radius) && (l[2]-sqrt(2)*Safe_Radius <= k[2] <= l[2]+sqrt(2)*Safe_Radius)
                            Dictionary_Vertices_Inside_Square[k] = true; break
                        else
                            Dictionary_Vertices_Inside_Square[k] = false;
                        end
                    end
                end

                #Get the Voronoi structure with only the non-repeated vertices
                Voronoi_Vertices = getVoronoiDiagram(Sites_Vertices);

                #Let's get a dictionary that relates "Vertex (X,Y) -> Index Voronoi's Polygon"
                Dictionary_Vertex_Index = diccionario_Centroides_Indice_Voronoi(Sites_Vertices, Voronoi_Vertices);

                #Get the index in the Vertex's Voronoi of the nearest vertex
                Nearest_Vertex_Index_Voronoi = Dictionary_Vertex_Index[Nearest_Vertex_End_Patch];

                #STEP 9: Re-fly the particle
                APoint,AVelocity_Initial,Remaining_Fly_Time,Nearest_Vertex_Index_Voronoi,Forbidden_Side_Polygon = vuelo_Una_Particula_Grafico(APoint,AVelocity_Initial,Remaining_Fly_Time,Nearest_Vertex_Index_Voronoi,Obstacle_Radius,Voronoi_Vertices,Dictionary_Vertex_Index,Dictionary_Vertices_Inside_Square,Forbidden_Side_Polygon);
                Nearest_Vertex_End_Patch = Voronoi_Vertices.faces[Nearest_Vertex_Index_Voronoi].site;
            end

            #Calculate the MSD of the system
            MSD_Actual = (APoint[1] - APoint_Initial[1])^2 + (APoint[2]- APoint_Initial[2])^2;
            MSD_Array[β] = (MSD_Array[β]*(α-1) + MSD_Actual)/α;

            #=
            if α == 1
                MSD_Array[β] = (APoint[1] - APoint_Initial[1])^2 + (APoint[2]- APoint_Initial[2])^2
            else
                MSD_Actual = (APoint[1] - APoint_Initial[1])^2 + (APoint[2]- APoint_Initial[2])^2
                MSD_Array[β] = (MSD_Array[β]*(α-1) + MSD_Actual)/α;
            end
            =#
        end
    end
    return MSD_Array
end