#A function that, given the coordinates P1=(X1,Y1) and P2=(X2,Y2), and given a Point=(X,Y), return the distance between
#the Point and the line that passes by the points P1 and P2
#Punto: The coordinates (X,Y) of a given point.
#P1: The coordinates (X1,Y1) of a given point P1.
#P2: The coordinates (X2,Y2) of a given point P2.
function distancia_Punto_Segmento(Punto, P1, P2)
    return abs((P2[2] - P1[2])*Punto[1] - (P2[1] - P1[1])*Punto[2] + P2[1]*P1[2] - P2[2]*P1[1])/(sqrt((P2[2] - P1[2])^2 + (P2[1] - P1[1])^2));
end

#A function that, given the coordinates of the vertices (X,Y) of the Quasiperiodic Array, determines the minimum
#distance between each vertex and the sides of his container polygon IF the vertex is inside of the Square Patch.
#The function return the minimum of this minimum distances.
#Voronoi: The Voronoi's Polygons' structure of the vertices of the Quasiperiodic Array
#Dictionario_Vertices_Parche_Cuadrado: A dictionary that relates Vertex (X,Y) -> "true" or "false" depending if the vertex is inside the square or not
function radio_Eficiente_Seguro(Voronoi, Dictionario_Vertices_Parche_Cuadrado)
    #Definimos el arreglo que contendrá la mínima distancia de un vértice a sus lados del polígono contenedor
    Arreglo_Minimas_Distancias = [];
    
    for Poligono_Voronoi in Voronoi.faces #Iteramos sobre todos los polígonos de Voronoi
        Vertice = Poligono_Voronoi.site; #Obtenemos el vértice del arreglo cuasiperiódico asociado al polígono
        if Dictionario_Vertices_Parche_Cuadrado[Vertice] #Verificamos si el vértice está dentro del parche cuadrado
            #Partimos de un lado del polígono de Voronoi que es de nuestro interés
            Lado_Poligono_Voronoi = Poligono_Voronoi.outerComponent;
            #Definimos el arreglo que contendrá las distancias del vértice a los lados del polígono
            Arreglo_Distancias = [];
            
            #Iniciamos el proceso while para recorrer todos los lados del polígono de Voronoi
            while true
                #Obtenemos los puntos que conforman el lado del polígono de Voronoi
                Punto1 = Lado_Poligono_Voronoi.origin.coordinates;
                Punto2 = Lado_Poligono_Voronoi.next.origin.coordinates;
                
                if (Punto1[1] - Punto2[1])^2 + (Punto1[2] - Punto2[2])^2 < 1e-6 #El lado es muy pequeño, es un falso lado del polígono Voronoi
                    nothing
                else
                    #Obtenemos la distancia entre el vértice y el segmento dado por [Punto1, Punto2]
                    Distancia = distancia_Punto_Segmento(Vertice, Punto1, Punto2);
                    push!(Arreglo_Distancias, Distancia)
                end
                
                #Recorremos al siguiente lado del polígono
                Lado_Poligono_Voronoi = Lado_Poligono_Voronoi.next;
                
                #Checamos si hemos ya concluido de revisar todos los lados del polígono de Voronoi
                if Lado_Poligono_Voronoi === Poligono_Voronoi.outerComponent
                    break
                end
            end
            
            #Obtenemos la mínima distancia del vértice a los lados de su polígono
            Minima_Distancia = minimum(Arreglo_Distancias);
            push!(Arreglo_Minimas_Distancias, Minima_Distancia);
        end
    end
    
    return minimum(Arreglo_Minimas_Distancias)
end

#A function that obtain the safe (and efficient) radius of the obstacles in a Quasiperiodic Array, looking in a given
#number of Square Patches.
#Number_Square_Patch: The number of Square Patches that the users wants to review to obtain the Safe Obstacle Radius
#Patch_Information: An array that contains information about the semi-circular patches of the Quasiperiodic Array
#Reduction_Factor: The factor with which we multiple the Average Radius to generate the Safe Radius
#Average_Distance_Stripes: Array with the average distance between stripes
#Star_Vectors: Array wich will contain the Star Vectors
#Alphas_Array: Array of the alphas constant
function radio_Eficiente_Seguro_Iterado(Number_Square_Patch, Patch_Information, Reduction_Factor, Average_Distance_Stripes, Star_Vectors, Alphas_Array)
    Array_Minimum_Distance_Patch = zeros(Float64, Number_Square_Patch);

    for j in 1:Number_Square_Patch
        println("The algorithm is calculating the Minimum Distance in the Square Patch $(j)")
        
        #Generate a square patch of the Quasiperiodic Array around the Initial Position
        X,Y,APoint_Initial,Centroids,Dictionary_Centroids = parche_Cuadrado(Patch_Information,Patch_Information[2],Reduction_Factor,Average_Distance_Stripes,Star_Vectors,Alphas_Array);

        #Obtain the vertices of all the polygons as duples.
        Sites_Vertices = [(Float64(X[i]), Float64(Y[i])) for i in 1:length(X)];
        unique!(Sites_Vertices); #Eliminate all the copies of a vertex

        #Get the Voronoi structure with only the non-repeated vertices
        Voronoi_Vertices = getVoronoiDiagram(Sites_Vertices);

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

        Safe_Obstacle_Radius = radio_Eficiente_Seguro(Voronoi_Vertices, Dictionary_Vertices_Inside_Square);
        Array_Minimum_Distance_Patch[j] = Safe_Obstacle_Radius;
    end

    return minimum(Array_Minimum_Distance_Patch)
end