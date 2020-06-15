struct DeadTrajectory <: Exception end

#Función que busca los vertices del arreglo cuasiperiódico más cercanos al punto de interés empleando polígonos de Voronoi.
#"Punto" es un arreglo dos dimensional con las coordenadas de un punto en el espacio 2D.
#"Vertices_Unicos" es un arreglo con las coordenadas (X,Y) de los vértices sin repetir.
function vertices_Cercanos_Voronoi(Punto, Vertices_Unicos)
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

#Función que dado un conjunto de posibles vértices a ser el vértice más cercano al punto de interes nos regresa el vertice más cercano.
#"Candidatos" es un arreglo con las coordenadas (X,Y) de los vértices cercanos al punto de interés en el plano.
#"Punto" es un arreglo dos dimensional con las coordenadas de un punto en el espacio 2D.
function vertice_Cercano(Candidatos, Punto)
    N = length(Candidatos);
    Distancias = zeros(BigFloat,N); #Arreglo que tendrá las distancias euclidianas del punto con los vertices
    Diccionario_Norma_Indice = Dict(); #Diccionario que relaciona "Norma -> Indice"
    for i in 1:N
        Norma = norm([Candidatos[i][1], Candidatos[i][2]] - Punto); #Calculamos la separación del punto con el i-ésimo vértice
        Distancias[i] += Norma
        Diccionario_Norma_Indice[Norma] = i; #Asociamos la norma con el índice
    end
    
    return Candidatos[Diccionario_Norma_Indice[minimum(Distancias)]]
end

#Función que dadas las condiciones iniciales de la partícula y un tiempo de vuelo dado, nos regresa la trayectoria de la partícula considerando colisiones con obstáculos.
#Posicion: Posicion inicial de la partícula.
#Velocidad: Velocidad inicial de la partícula.
#Tiempo_Vuelo: Tiempo de vuelo deseado que recorra la partícula.
#Indice_Obstaculo_Cercano: Indice en la lista de Voronoi del obstáculo más cercano a la partícula.
#Radio_Obstaculo: Radio de los obstáculos.
#Voronoi: Estructura de voronoi de los vertices.
#Diccionario_Vertices_Indices: Diccionario "Vertices (X,Y) -> Indices_Voronoi".
#Diccionario_Vertices_Cuadrado: Diccionario "Vertices (X,Y) -> Booleano" que nos indica si el polígono es un polígono dentro de algún parche cuadrado o no.
#Lado_Poligono_Bloqueado: Información sobre el lado del polígono por el que la partícula entró.
function vuelo_Una_Particula_Obstaculos(Posicion, Velocidad, Tiempo_Vuelo, Indice_Obstaculo_Cercano, Voronoi, Diccionario_Vertices_Indices, Diccionario_Vertices_Cuadrado, Lado_Poligono_Bloqueado, funcion_Cambio_Celda::Function, funcion_Estado_Tras_Tiempo::Function)
    while Tiempo_Vuelo > 0.0
        #Chequemos si el polígono contenedor de la partícula está en los bordes exteriores del parche cuadrado seguro.
        #Si lo está, no calculamos más, nos salimos de la función para generar otro parche centrado en esta posición
        if Diccionario_Vertices_Cuadrado[Voronoi.faces[Indice_Obstaculo_Cercano].site] == false
            return Posicion, Velocidad, Tiempo_Vuelo, Indice_Obstaculo_Cercano, Lado_Poligono_Bloqueado
        end

        #Guardamos las coordenadas X,Y del obstáculo más cercano a nuestra partícula
        Posicion_Obstaculo = [Voronoi.faces[Indice_Obstaculo_Cercano].site[1], Voronoi.faces[Indice_Obstaculo_Cercano].site[2]];
        Posicion_Inicial_Cambio = copy(Posicion);
        Velocidad_Inicial_Cambio = copy(Velocidad);
        Lado_Poligono_Bloqueado_Cambio = copy(Lado_Poligono_Bloqueado);

        #Generamos las variables que se requieren para "cambio_Poligono" y potencialmente para "funcion_Estado_Tras_Tiempo"
        Diccionario_Segmento_Indice = Dict(); #Diccionario que relaciona los extremos [[X1,Y1], [X2,Y2]] con el índice del polígono exterior que comparte ese lado
        Diccionario_Segmentos_Centro_Vecino = Dict(); #Diccionario que relaciona los los extremos [[X1,Y1], [X2,Y2]] con la posición del obstáculo de dicha celda vecina
        Arreglo_Vertices_Contenedor = []; #Arreglo donde se guardarán las coordenadas [X,Y] de los vértices de la celda contenedora

        #Partimos de un lado del polígono de Voronoi que es de nuestro interés
        Lado_Poligono_Voronoi = Voronoi.faces[Indice_Obstaculo_Cercano].outerComponent;
            
        #Iniciamos el proceso while para recorrer todos los lados del polígono de Voronoi y en cada uno hallar el polígono vecino
        while true
            #Encontramos el indice del vecino asociado al lado que estamos considerando
            Obstaculo_Vecino = Lado_Poligono_Voronoi.twin.incidentFace.site;
            Indice_Vecino = Diccionario_Vertices_Indices[Lado_Poligono_Voronoi.twin.incidentFace.site];

            #Obtenemos los puntos que conforman el lado del polígono de Voronoi
            Punto1 = Lado_Poligono_Voronoi.origin.coordinates;
            Punto2 = Lado_Poligono_Voronoi.next.origin.coordinates;

            push!(Arreglo_Vertices_Contenedor, [Punto1[1], Punto1[2]]);
            Diccionario_Segmento_Indice[[Punto1[1], Punto1[2]], [Punto2[1], Punto2[2]]] = Indice_Vecino;
            Diccionario_Segmentos_Centro_Vecino[[Punto1[1], Punto1[2]], [Punto2[1], Punto2[2]]] = Obstaculo_Vecino;
                
            #Recorremos al siguiente lado del polígono
            Lado_Poligono_Voronoi = Lado_Poligono_Voronoi.next;

            #Checamos si hemos ya concluido de revisar todos los lados del polígono de Voronoi
            if Lado_Poligono_Voronoi === Voronoi.faces[Indice_Obstaculo_Cercano].outerComponent
                break
            end
        end

        #Empleamos la función del usuario para que nos regrese la posición y velocidad de la partícula al salir del polígono contenedor actual, así como el tiempo de vuelo de en salir del polígono
        #y las coordenadas [[X1, Y1], [X2, Y2]] del segmento por el que sale la partícula
        Posicion_Cambio_Celda, Velocidad_Cambio_Celda, Tiempo_Vuelo_Cambio, Lado_Poligono_Bloqueado_Cambio = funcion_Cambio_Celda(Posicion_Inicial_Cambio, Velocidad_Inicial_Cambio, Arreglo_Vertices_Contenedor, Lado_Poligono_Bloqueado_Cambio, Posicion_Obstaculo, Diccionario_Segmentos_Centro_Vecino);

        #Definimos la variable del índice de centro de Voronoi de la celda receptora
        Indice_Proximo_Obstaculo = 0;

        #Obtengamos el índice de la celda contenedora receptora
        if Tiempo_Vuelo_Cambio < Inf
            try
                Indice_Proximo_Obstaculo = Diccionario_Segmento_Indice[Lado_Poligono_Bloqueado_Cambio[1], Lado_Poligono_Bloqueado_Cambio[2]];
            catch
                Indice_Proximo_Obstaculo = Diccionario_Segmento_Indice[Lado_Poligono_Bloqueado_Cambio[2], Lado_Poligono_Bloqueado_Cambio[1]];
            end
        else
            #Dependiendo de la trayectoria puede ocurrir que la partícula no colisione con obstáculos ni cambie de polígono, en ese caso la partícula se encuentra en
            #una trayectoria muerta, por ejemplo una trayectoria circular cerrada, en ese caso el algoritmo debe terminar ahí.
            Indice_Proximo_Obstaculo = Indice_Obstaculo_Cercano;
            throw(DeadTrajectory())
        end
        
        if (Tiempo_Vuelo - Tiempo_Vuelo_Cambio) > 0.0 #Hay tiempo para que la partícula salga del polígono
            Tiempo_Vuelo -= (Tiempo_Vuelo_Cambio); #Actualizamos el tiempo de vuelo restante
            Posicion = Posicion_Cambio_Celda; #Avancemos la posición de la partícula al nuevo polígono
            Velocidad = Velocidad_Cambio_Celda; #Actualizamos la velocidad de la partícula a la velocidad que tiene cuando cambia de celda
            Indice_Obstaculo_Cercano = Indice_Proximo_Obstaculo; #Actualizamos el índice del polígono contenedor
            Lado_Poligono_Bloqueado = Lado_Poligono_Bloqueado_Cambio; #El lado del polígono bloqueado cambia
        else #No hay tiempo de vuelo restante para salir del polígono
            Posicion, Velocidad, Lado_Poligono_Bloqueado = funcion_Estado_Tras_Tiempo(Posicion, Velocidad, Arreglo_Vertices_Contenedor, Lado_Poligono_Bloqueado, Posicion_Obstaculo, Diccionario_Segmentos_Centro_Vecino, Tiempo_Vuelo);
            return Posicion, Velocidad, 0.0, Indice_Obstaculo_Cercano, Lado_Poligono_Bloqueado
        end
    end
end

#Function that, given a initial position in the Quasiperiodic Array, calculates the MSD of the system iterating over a given number of flies (each one conformed as 
#a set of short-time flies), each one with a random velocity.
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
#funcion_Cambio_Celda: Function that calculates the position and velocity of a particle when it leaves the container polygon, as well as the time it requires to
#leave the polygon and the coordinates of the vertices of the segment where it leaves.
#funcion_Estado_Tras_Tiempo: Function that calculates the position and velocity of a particle, as well as the vertices of the segment from which the particle enter
#in the latest move.
function MSD_Varying_Velocities(Patch_Information, Reduction_Factor, Average_Distance_Stripes, Star_Vectors, Alphas_Array, APoint_Initial, Number_Velocities, Number_Flights, Short_Fly_Time, MSD_Array, funcion_Cambio_Celda::Function, funcion_Estado_Tras_Tiempo::Function)
    Average_Radius = Patch_Information[end]; #The value of the previously calculated "Average Radius"
    Safe_Radius = Reduction_Factor*Average_Radius; #The Safe Radius value

    #STEP 1: Generate a square patch of the Quasiperiodic Array around the Initial Position
    X,Y = parche_Cuadrado(Patch_Information,Reduction_Factor,Average_Distance_Stripes,Star_Vectors,Alphas_Array,APoint_Initial);

    #STEP 2: Get the structure of the Voronoi's Polygons of the Vertices that conform the Quasiperiodic Array
    Sites_Vertices = [(Float64(X[i]), Float64(Y[i])) for i in 1:length(X)]; #Obtain the vertices of all the polygons as duples.
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

    #STEP 3: Find the container polygon of the Initial Position
    APoint = copy(APoint_Initial); #A copy of Initial Position in order to preserve the original one

    #Let's obtain the vertices coordinates of the possible nearest vertex to the point
    Nearest_Vertex_Candidates = vertices_Cercanos_Voronoi(APoint, Sites_Vertices);

    #Get the index of the nearest vertex from the list of candidates
    Nearest_Vertex = vertice_Cercano(Nearest_Vertex_Candidates, APoint);

    Live_Flies = 0; #Contador para el número de buenas trayectorias realizadas (Trayectorias "vivas", es decir, no-cerradas)
    
    #STEP 4: Iterates over the desired number of velocities with the fixed initial position
    for α in 1:Number_Velocities
        println("The algorithm is calculating the particle $(α).")
        Dead_Trajectory = false; #We suppose that the trajectory won't be a closed one

        #STEP 5: Get a Initial Velocity
        θ = 2*π*rand(); #Get a random angle between 0 and 2*pi
        AVelocity_Initial = [cos(θ), sin(θ)]; #Generate a unitary velocity
        AVelocity = copy(AVelocity_Initial); #A copy of Initial Velocity in order to preserve the original one
        APoint = copy(APoint_Initial); #A copy of Initial Position in order to preserve the original one

        #Get the index in the Vertex's Voronoi of the nearest vertex
        Nearest_Vertex_Index_Voronoi = Dictionary_Vertex_Index[Nearest_Vertex];

        #At the begin, all the sides of the container polygons are valid
        Forbidden_Side_Polygon = [[Inf, Inf], [Inf, Inf]];

        #STEP 6: Let's the particle fly
        for β in 1:Number_Flights
            #Make a little fly
            Remaining_Fly_Time = Short_Fly_Time;

            try
                APoint,AVelocity,Remaining_Fly_Time,Nearest_Vertex_Index_Voronoi,Forbidden_Side_Polygon = vuelo_Una_Particula_Obstaculos(APoint,AVelocity,Short_Fly_Time,Nearest_Vertex_Index_Voronoi,Voronoi_Vertices,Dictionary_Vertex_Index,Dictionary_Vertices_Inside_Square,Forbidden_Side_Polygon,funcion_Cambio_Celda,funcion_Estado_Tras_Tiempo);
            catch Err
            	if isa(Err, DeadTrajectory)
            		println("The particle $(α) is a dead one")
            		println("La posicion inicial es: $(APoint_Initial)")
            		println("La velocidad inicial es: $(AVelocity_Initial)")
            		Dead_Trajectory = true;
            		break
            	else
            		println("Check the initial conditions, something went wrong")
            		println("La posicion inicial es: $(APoint_Initial)")
            		println("La velocidad inicial es: $(AVelocity_Initial)")
            		throw(ErrorException())
            	end
            end

            Nearest_Vertex_End_Patch = Voronoi_Vertices.faces[Nearest_Vertex_Index_Voronoi].site; #X,Y coordinates of the obstacle that correspond to the container polygon

            while Remaining_Fly_Time > 0.0
                println("Estamos generando un parche nuevo")
                #STEP 7: Generate a new square patch if the previous one wasn't enough
                X,Y = parche_Cuadrado(Patch_Information,Reduction_Factor,Average_Distance_Stripes,Star_Vectors,Alphas_Array,APoint);
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
                APoint,AVelocity,Remaining_Fly_Time,Nearest_Vertex_Index_Voronoi,Forbidden_Side_Polygon = vuelo_Una_Particula_Obstaculos(APoint,AVelocity,Remaining_Fly_Time,Nearest_Vertex_Index_Voronoi,Voronoi_Vertices,Dictionary_Vertex_Index,Dictionary_Vertices_Inside_Square,Forbidden_Side_Polygon,funcion_Cambio_Celda,funcion_Estado_Tras_Tiempo);
                Nearest_Vertex_End_Patch = Voronoi_Vertices.faces[Nearest_Vertex_Index_Voronoi].site;
            end

            #Calculate the MSD of the system
            MSD_Actual = (APoint[1] - APoint_Initial[1])^2 + (APoint[2]- APoint_Initial[2])^2;
            MSD_Array[β] = (MSD_Array[β]*(Live_Flies) + MSD_Actual)/(Live_Flies+1);
        end

        #If the trajectory is not a dead one, we add this one in the counter of the validate trajectories to calculate the MSD
        if Dead_Trajectory == false
            Live_Flies += 1;
        end
    end
    return MSD_Array, Live_Flies
end