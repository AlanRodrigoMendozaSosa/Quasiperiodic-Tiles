function areas_Poligonos_Voronoi(Error_Margin, SL, Bounded_Area, Reduction_Factor, Average_Distance_Stripes, Star_Vectors, Alphas_Array)
    APoint = punto_Arbitrario(SL); #We generate an arbitrary point inside the square centered on the (0,0)
    
    #We generate the local neighborhood around the arbitrary point
    Dual_Points = region_Local_Voronoi(Error_Margin, Average_Distance_Stripes, Star_Vectors, Alphas_Array, APoint);

    #Let's get the coordinates as tuples of the centroids and the dictionary that relates the centroid's coordinates 
    #with the polygons vertices' coordinates of the polygon that generate the centroid.
    Centroids, Centroids_Dictionary = centroides(Dual_Points);

    #Let's get the initial Voronoi's Lattice
    Sites = [(Centroids[i][1], Centroids[i][2]) for i in 1:length(Centroids)]
    Initial_Voronoi = getVoronoiDiagram(Sites);

    #Let's get the centroids that remains after the iterations of the areas algorithm
    Inside_Clusters_Centroids, Radius = centroides_Area_Acotada(Initial_Voronoi, Bounded_Area, APoint);

    #Let's get the X and Y coordinates of the vertices of the retained polygons in the quasiperiodic lattice
    X,Y = centroides_A_Vertices(Inside_Clusters_Centroids, Centroids_Dictionary);
    
    #Get the structure of the Voronoi's Polygons of the Vertices that conform the Quasiperiodic Array
    Sites_Vertices = [(Float64(X[i]), Float64(Y[i])) for i in 1:length(X)]; #Obtain the vertices of all the polygons as duples.
    unique!(Sites_Vertices); #Eliminate all the copies of a vertex
    
    #Let's generate a Dictionary with the coordinates of the vertices that lay inside the circle.
    #The Dictionary relates "Vertex (X,Y) -> "true" or "false" depending if the vertex is inside the circle or not
    Dictionary_Vertices_Inside_Circle = Dict();
    for i in Sites_Vertices
        if norm([i[1], i[2]] - APoint) > Reduction_Factor*Radius
            Dictionary_Vertices_Inside_Circle[i] = false;
        else
            Dictionary_Vertices_Inside_Circle[i] = true;
        end
    end
    
    #Get the Voronoi structure with the non-repeated vertices
    Voronoi_Vertices = getVoronoiDiagram(Sites_Vertices);

    #Let's get a dictionary that relates "Vertex (X,Y) -> Index Voronoi's Polygon"
    Dictionary_Vertex_Index = diccionario_Centroides_Indice_Voronoi(Sites_Vertices, Voronoi_Vertices);
    
    #Arreglo donde irán las áreas de los polígonos de Voronoi que corresponden a los sitios del círculo de seguridad
    Areas_Voronoi = [];
    
    for i in Sites_Vertices
        #Si el vértice está dentro del círculo de seguridad, lo consideramos. Caso contrario lo ignoramos
        if Dictionary_Vertices_Inside_Circle[i] == true
            push!(Areas_Voronoi, Voronoi_Vertices.faces[Dictionary_Vertex_Index[i]].area);
        else
            nothing
        end
    end

    return Areas_Voronoi
end