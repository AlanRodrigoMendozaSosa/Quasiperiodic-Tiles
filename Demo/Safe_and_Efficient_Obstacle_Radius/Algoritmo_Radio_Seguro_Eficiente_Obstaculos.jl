#A function that, given the Voronoi structure of the vertices of the Quasiperiodic lattice, return an array with the minimum
#distance between each vertex and his neighbors IF the vertex is inside of the Square Patch.
#Voronoi: The Voronoi's Polygons' structure of the vertices of the Quasiperiodic Array
#Dictionario_Vertices_Parche_Cuadrado: A dictionary that relates Vertex (X,Y) -> "true" or "false" depending if the vertex is inside the square or not
function arreglo_Minimas_Distancias(Voronoi, Diccionario_Vertices_Parche_Cuadrado)
    #Definimos el arreglo que contendrá la mínima distancia de un vértice a sus lados del polígono contenedor
    Arreglo_Minimas_Distancias = [];
    
    for i in 1:length(Voronoi.faces)
        Arreglo_Distancias_Vertice_Vecinos = Float64[];
        Vertice = Voronoi.faces[i].site;
        if Diccionario_Vertices_Parche_Cuadrado[Vertice] #Verificamos si el vértice está dentro del parche cuadrado
            Vertices_Vecinos = vecinos_Voronoi(i, Voronoi)
            for j in Vertices_Vecinos
                Distancia = sqrt((j[1] - Vertice[1])^2 + (j[2] - Vertice[2])^2)
                push!(Arreglo_Distancias_Vertice_Vecinos, Distancia)
            end
            push!(Arreglo_Minimas_Distancias, minimum(Arreglo_Distancias_Vertice_Vecinos))
        end
    end
    
    return Arreglo_Minimas_Distancias
end

#A function that, given the coordinates of the vertices (X,Y) of the Quasiperiodic Array, determines the minimum
#radius of the obstacles that ensures that the obstacle will be contained in the corresponding Voronoi's polygon
#Voronoi: The Voronoi's Polygons' structure of the vertices of the Quasiperiodic Array
#Dictionario_Vertices_Parche_Cuadrado: A dictionary that relates Vertex (X,Y) -> "true" or "false" depending if the vertex is inside the square or not
function radio_Eficiente_Seguro(Voronoi, Diccionario_Vertices_Parche_Cuadrado)
    Arreglo_Candidatos_Radio = arreglo_Minimas_Distancias(Voronoi, Diccionario_Vertices_Parche_Cuadrado);

    return minimum(Arreglo_Candidatos_Radio)/2
end

#A function that obtain the safe (and efficient) radius of the obstacles in a Quasiperiodic Array, looking in a given number of Square Patches.
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
        
        #Generate an arbitrary point in the 2D plane
        APoint_Initial = punto_Arbitrario(Patch_Information[2]);

        #Generate a square patch of the Quasiperiodic Array around the Initial Position
        X,Y = parche_Cuadrado(Patch_Information,Reduction_Factor,Average_Distance_Stripes,Star_Vectors,Alphas_Array,APoint_Initial);

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
        Array_Minimum_Distance_Patch[j] += Safe_Obstacle_Radius;
    end

    return minimum(Array_Minimum_Distance_Patch)
end