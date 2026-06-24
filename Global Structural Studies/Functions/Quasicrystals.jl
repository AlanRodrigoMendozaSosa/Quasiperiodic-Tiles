using LinearAlgebra #Paquetería para emplear funciones de álgebra lineal
using LsqFit #Paquetería para realizar ajuste de mínimos cuadrados
include("voronoi.jl") #Conjuntos de funciones y algoritmos programados por Enrique para generar teselados de Voronoi de manera eficiente.
########################################################################################################################################################################
########################################################################## Operaciones Basicas #########################################################################
########################################################################################################################################################################
#Función que calcula el producto punto entre dos vectores 2D.
function prod_Punto(A, B)
    return A[1]*B[1] + A[2]*B[2]
end

#Función que calcula la norma Euclideana de un vector en 2D.
function norma_Vector(A)
    return sqrt(A[1]^2 + A[2]^2)
end

#Función que genera el vector ortogonal de un vector en 2D.
function vector_Ortogonal(A)
    return [A[2], -A[1]]
end

#Función que nos genera un punto aleatorio en un cuadrado de semilado SL centrado en el origen.
function punto_Arbitrario(SL)
    #Definimos la variable Punto que contendrá las coordenadas del punto de interés
    Punto = zeros(Float64,2);

    #Llaves para determinar el cuadrante donde estará el punto
    x = rand();
    y = rand();

    if (x > 0.5) && (y > 0.5)
        Punto = [rand()*SL, rand()*SL]; #Primer cuadrante
    elseif (x > 0.5) && (y < 0.5)
        Punto = [rand()*SL, -rand()*SL]; #Cuarto cuadrante
    elseif (x < 0.5) && (y > 0.5)
        Punto = [-rand()*SL, rand()*SL]; #Segundo cuadrante
    elseif (x < 0.5) && (y < 0.5)
        Punto = [-rand()*SL, -rand()*SL]; #Tercer cuadrante
    end

    return Punto
end
########################################################################################################################################################################
##################################################################### Generación del teselado Cuasiperiódico ###########################################################
########################################################################################################################################################################
#Función que determina, fijando los vectores estrella Ej y Ek, fijando los enteros Nj y Nk, y para algún conjunto de constantes alfa, los cuatros puntos de la red dual.
#"J" y "K" son los índices de los vectores estrella a considerar.
#"Nj" y "Nk" son los enteros con los que se generan las rectas ortogonales a Ej y Ek.
#"Vectores_Estrella" es el conjunto de vectores estrella.
#"Arreglo_Alfas" es el arreglo con los valores numéricos de la separación respecto al origen del conjunto de rectas ortogonales a los vectores estrella en el método generalizado dual.
function cuatro_Regiones(J, K, Nj, Nk, Vectores_Estrella, Arreglo_Alfas)
    #Paso 0: Verifiquemos si los vectores a considerar son colineales, en cuyo caso manda un error.
    if (length(Vectores_Estrella) % 2 == 0) && (K == J + length(Vectores_Estrella)/2)
        error("Los vectores Ej y Ek no pueden ser paralelos")
    else
        #Paso 1: Definimos los dos vectores con los que se consigue la intersección en la malla generada por los vectores estrella que estamos considerando.
        Ej = Vectores_Estrella[J];
        Ek = Vectores_Estrella[K];
        
        #Paso 2: Obtenemos los vectores ortogonales a estos dos vectores.
        EjOrt = vector_Ortogonal(Ej);
        EkOrt = vector_Ortogonal(Ek);
        
        #Paso 3: Definimos los valores reales con los que se crearon las rectas ortogonales a cada vector Ej y Ek para tomar la intersección.
        Factor_Ej = Nj + Arreglo_Alfas[J];
        Factor_Ek = Nk + Arreglo_Alfas[K];
        
        #Paso 4: Obtenemos el área que forman los dos vectores Ej y Ek.
        AreaJK = Ej[1]*Ek[2] - Ej[2]*Ek[1];
        
        #Paso 5: Definimos lo que será el vértice en el espacio real de la retícula cuasiperiódica. Este vértice se denomina t^{0} en la teoría.
        Punto_Cero = Nj*Ej + Nk*Ek;
        
        #Paso 6: Generamos los términos asociados a la proyección del vector Ej y Ek con los demás vectores estrella
        for i in 1:length(Vectores_Estrella)
            if i == J || i == K
                nothing
            else
                Factor_Ei = (Factor_Ej/AreaJK)*(prod_Punto(EkOrt, Vectores_Estrella[i])) - (Factor_Ek/AreaJK)*(prod_Punto(EjOrt, Vectores_Estrella[i]));
                Punto_Cero += (floor(Factor_Ei - Arreglo_Alfas[i]))*Vectores_Estrella[i];
            end
        end
        
        #Paso 7: Obtenemos los otros tres vértices asociados al punto t^{0}.
        Punto_Uno = Punto_Cero - Ej;
        Punto_Dos = Punto_Cero - Ej - Ek;
        Punto_Tres = Punto_Cero - Ek;
        
        return Punto_Cero, Punto_Uno, Punto_Dos, Punto_Tres
    end
end

#Función que dado un Punto, obtiene aproximadamente los enteros que forman el polígono que lo contiene.
#"Punto" es un arreglo [X,Y] con las coordenadas de un punto en el espacio 2D.
#"Promedios_Distancia" es el arreglo con la separación entre las franjas cuasiperiódicas.
#"Vectores_Estrella" es el arreglo con los vectores estrella que generan la retícula cuasiperiódica deseada.
function proyecciones_Pto_Direccion_Franjas(Punto, Promedios_Distancia, Vectores_Estrella)
    #Paso 1: Generemos un arreglo en donde irán los números reales resultado de proyectar el Punto con los vectores estrella.
    Arreglo_Proyeccion_Pto_Direccion_Franjas = [];
    
    #Paso 2: Para cada vector estrella, proyectamos el Punto sobre dicho vector y reescalamos con la separación entre las franjas cuasiperiódicas.
    for i in 1:length(Vectores_Estrella)
        Proyeccion = prod_Punto(Punto, Vectores_Estrella[i]/norm(Vectores_Estrella[i]))/Promedios_Distancia[i];
        #Al parecer no hay que restarle "- Arreglo_Distancia_Origen_Primera_Recta[i]" que es la DOPR
        push!(Arreglo_Proyeccion_Pto_Direccion_Franjas, Proyeccion);
    end
    
    return Arreglo_Proyeccion_Pto_Direccion_Franjas
end

#Función que genera los vértices de un arreglo cuasiperiódico asociados a la vecindad de un Punto arbitrario.
#"Proyecciones" es el conjunto de números enteros candidatos a ser los que generan el polígono que contiene al punto.
#"Vectores_Estrella" es el conjunto de vectores estrella.
#"Arreglo_Alfas" es el arreglo con los valores numéricos de la separación respecto al origen del conjunto de rectas ortogonales a los vectores estrella en el método generalizado dual.
#"N" es el margen de error asociado a los números enteros generados por la proyección del punto sobre los vectores estrella.
function generador_Vecindades_Vertices(Proyecciones, Vectores_Estrella, Arreglo_Alfas, N)
    #Paso 1: Definimos el arreglo que contendrá a los vértices asociados a cada combinación de vectores estrella (incluídos los generados al considerar el margen de error)
    Puntos_Red_Dual = [];
    
    #Paso 2: Consideramos todas las posibles combinaciones de vectores estrella con los posibles números enteros correspondientes (y su margen de error)
    for i in 1:length(Vectores_Estrella)
        for j in i+1:length(Vectores_Estrella)
            #Consideramos el margen de error a cada número entero
            for n in -N:N
                for m in -N:N
                    #Vamos a dejar que el try ... catch se encargue de los casos en que los vectores estrella sean paralelos
                    try
                        #Obtengamos los vértices del arreglo considerando los vectores Ei y Ej con sus respectivos números enteros (y su margen de error)
                        #Salida: [X,Y]
                        t0, t1, t2, t3 = cuatro_Regiones(i, j, round(Proyecciones[i])+n, round(Proyecciones[j])+m, Vectores_Estrella, Arreglo_Alfas);
                        push!(Puntos_Red_Dual, t0);
                        push!(Puntos_Red_Dual, t1);
                        push!(Puntos_Red_Dual, t2);
                        push!(Puntos_Red_Dual, t3);
                    catch
                        nothing;
                    end
                end
            end
            
        end
    end
    
    return Puntos_Red_Dual
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

#Función que nos regresa, dado un conjunto de centroides, los vértices de los polígonos que los generan. Emplea un diccionario para ello.
#"Centroides" es un arreglo con las duplas (X,Y) de los centroides de interés.
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

#Función que mantiene los centroides de las teselas del cluster principal
#"Voronoi" es la estructura generada por Enrique con getVoronoiDiagram().
#"Cota_Area" es el área que servirá como discriminante parar separar a los polígonos dentro de clúster de los de frontera.
#"APoint" es el punto en el espacio alrededor del cual se genera la vecindad cuasiperiódica.
function centroides_Area_Acotada(Voronoi, Cota_Area, APoint)
    Arreglo_Centroides = []; #Arreglo con Centroides de Celdas Voronoi del cluster principal
    Distancia_Frontera_Cercana = Inf;

    #Iteramos sobre todos los polígonos asociados a centroides
    for Face in Voronoi.faces
        if Face.area > Cota_Area #Compara el área del polígono en turno contra el área de la cota dada
            Distancia_Celda = norm([Face.site[1], Face.site[2]] - APoint);
            Distancia_Celda  < Distancia_Frontera_Cercana ? (Distancia_Frontera_Cercana = Distancia_Celda) : (nothing);
        end
    end
    
    for Face in Voronoi.faces
        if Face.area < Cota_Area #Compara el área del polígono en turno contra el área de la cota dada
            Distancia_Celda = norm([Face.site[1], Face.site[2]] - APoint);
            Distancia_Celda  < Distancia_Frontera_Cercana ? (push!(Arreglo_Centroides, Face.site)) : (nothing)
        end
    end

    return Arreglo_Centroides, Distancia_Frontera_Cercana
end

#Función que estima el margen de error necesario para construir una vecindad de radio R de un sistema cuasiperiódico de simetría
#rotacional N.
#"NSides" es la simetría rotacional del sistema cuasiperiódico.
#"Desired_Radius" es el radio deseado para la vecindad circular a construir.
#"Alpha_Value" es el valor del parámetro "alpha" a usar en la construcción del sistema cuasiperiódico por el GDM.
#"Beta_Start" es el margen de error inicial para estimar la recta que mejor describa la relación Beta vs R.
#"Beta_End" es el margen de error final para estimar la recta que mejor describa la relación Beta vs R.
function beta_Calculation(NSides, Desired_Radius; Alpha_Value = 0.0, Beta_Start = 5, Beta_End = 10)
    Star_Vectors = [[BigFloat(1),0]]; #Arrangement that will contain the star vectors
    for i in 1:(NSides-1)
        push!(Star_Vectors, [cos((2*i)*pi/NSides), sin((2*i)*pi/NSides)]); #Vertices of the polygon with "NSides" sides
    end
    Alphas_Array = fill(Alpha_Value, NSides); #Array with the alpha constants of the GDM
    Average_Distance_Stripes = fill(NSides/2, NSides); #Array with the average spacing between stripes

    βArray = []; #Array that will contain the beta value used to generate a Circular Cluster of radius R
    RArray = []; #Array that will contain the value of the Circular Cluster generated with the beta value

    for β in Beta_Start:Beta_End
        Error_Margin = β; #Margin of error in the integers "n" associated with the star vectors of the GDM
        SL = 1e6; #Size of the half side of the square centered on the (0,0) where a starting position is sought.
        APoint = punto_Arbitrario(SL); #We generate an arbitrary point inside the square centered on the (0,0)
        
        #We generate the local neighborhood around the Test_Point point
        Dual_Points = region_Local_Voronoi(Error_Margin, Average_Distance_Stripes, Star_Vectors, Alphas_Array, APoint);

        #Let's get the coordinates as tuples of the centroids and the dictionary that relates the centroid's coordinates 
        #with the polygons vertices' coordinates of the polygon that generate the centroid.
        Centroids, Centroids_Dictionary = centroides(Dual_Points);

        #Let's get the initial Voronoi's Lattice
        Sites = [(Centroids[i][1], Centroids[i][2]) for i in 1:length(Centroids)]
        Initial_Voronoi = getVoronoiDiagram(Sites);

        #The value for the area of the polygons that will be a discriminator value in the areas algorithm
        Bounded_Area = 1.2;

        #Let's get the centroids that remains after the iterations of the areas algorithm
        Inside_Clusters_Centroids, Radius = centroides_Area_Acotada(Initial_Voronoi, Bounded_Area, APoint);
        
        push!(βArray, β); push!(RArray, Radius)
    end

    @. model(x, p) = p[1]*x + p[2]; #Linear Model for the Radius of the quasiperiodic circular patch as a function of the beta parameter
    xdata = βArray; #Independent variable
    ydata = RArray; #Dependent variable
    p0 = [0.5, 0.5]; #Starting guess values for the linear model parameters

    fit = curve_fit(model, xdata, ydata, p0);
    a, b = fit.param;

    Requiredβ = Int(ceil((Desired_Radius - b)/a)); #Estimated beta value to produce a quasiperiodic circular patch of radius 'Desired_Radius'
    println("The estimated beta value to produce a Quasiperiodic Circular neighbourhood with R = $(Desired_Radius) and N = $(NSides) is β = $(Requiredβ)")
end

#Función que genera una vecindad circular de un sistema cuasiperiódico de simetría rotacional N alrededor de un punto arbitrario.
#Esta función requiere que el punto arbitrario esté a lo más dentro de un cuadrado centrado en el origen de lado 2e6, esto debido a que hace uso del algoritmo
#de voronoi.jl para eliminar las teselas basura de la vecindad central. Su uso es principalmente para validar la estimación del parámetro Beta realizado por
#la función "beta_Calculation" y se desaconseja su uso en simulaciones largas.
#"NSides" es la simetría rotacional del sistema cuasiperiódico.
#"Error_Margin" es el margen de error asociado a los enteros Nj correspondientes a los vectores estrella. También determina el tamaño de la vecindad.
#"APoint" es el punto alrededor del cual se va a generar el sistema cuasiperiódico.
#"Alpha_Value" es el valor del parámetro "alpha" a usar en la construcción del sistema cuasiperiódico por el GDM.
function quasiperiodic_Neighbourhood_Voronoi_Version(NSides, Error_Margin, APoint; Alpha_Value = 0.0)
    Star_Vectors = [[BigFloat(1),0]]; #Arrangement that will contain the star vectors
    for i in 1:(NSides-1)
        push!(Star_Vectors, [cos((2*i)*pi/NSides), sin((2*i)*pi/NSides)]); #Vertices of the polygon with "NSides" sides
    end
    Alphas_Array = fill(Alpha_Value, NSides); #Array with the alpha constants of the GDM
    Average_Distance_Stripes = fill(NSides/2, NSides); #Array with the average spacing between stripes
    ###################################################################################################################
    #SL = 1e6; #Size of the half side of the square centered on the (0,0) where a starting position is sought.
    #APoint = punto_Arbitrario(SL); #We generate an arbitrary point inside the square centered on the (0,0)

    #We generate the local neighborhood around the Test_Point point
    Dual_Points = region_Local_Voronoi(Error_Margin, Average_Distance_Stripes, Star_Vectors, Alphas_Array, APoint);

    #Let's get the coordinates as tuples of the centroids and the dictionary that relates the centroid's coordinates with the polygons vertices' coordinates of the polygon that generate the centroid.
    Centroids, Centroids_Dictionary = centroides(Dual_Points);

    #Let's get the initial Voronoi's Lattice
    Sites = [(Centroids[i][1], Centroids[i][2]) for i in 1:length(Centroids)]
    Initial_Voronoi = getVoronoiDiagram(Sites);

    #The value for the area of the polygons that will be a discriminator value in the areas algorithm
    Bounded_Area = 1.2;

    #Let's get the centroids that remains after the iterations of the areas algorithm
    Inside_Clusters_Centroids, Radius = centroides_Area_Acotada(Initial_Voronoi, Bounded_Area, APoint);

    #Let's get the X and Y coordinates of the vertices of the retained polygons in the quasiperiodic lattice
    X, Y = centroides_A_Vertices(Inside_Clusters_Centroids, Centroids_Dictionary);
    ###################################################################################################################
    #Get the structure of the Voronoi's Polygons of the Vertices that conform the Quasiperiodic Array
    Sites_Vertices = [(Float64(X[i]), Float64(Y[i])) for i in 1:length(X)]; #Obtain the vertices of all the polygons as duples.
    unique!(Sites_Vertices); #Eliminate all the copies of a vertex
    #unique!(x->(round(x[1],digits = 5), round(x[2],digits = 5)), Sites_Vertices); #Versión alterna del unique!() para F64
    #Sites_Vertices = [round.(x, digits = 5) for x in Sites_Vertices]; #Redondeamos realmente y no sólo quitamos copias
    println("The quasiperiodic neighbourhood has a radius of R = $(Radius)")
    return X, Y, Sites_Vertices
end

#Función que genera una vecindad circular de radio "Radius" alrededor del punto "APoint" de un sistema cuasiperiódico de simetría rotacional "NSides". Para este algoritmo es necesario conocer
#a priori el valor del margen de error "Error_Margin" asociado al valor del radio "Radius", para lo cual es posible usar la función beta_Calculation.
#"NSides" es la simetría rotacional del sistema cuasiperiódico.
#"Error_Margin" es el margen de error asociado a los enteros Nj correspondientes a los vectores estrella. También determina el tamaño de la vecindad.
#"Radius" es el radio que tendrá la vecindad circular del sistema cuasiperiódico.
#"APoint" es el punto alrededor del cual se va a generar el sistema cuasiperiódico.
#"Alpha_Value" es el valor del parámetro "alpha" a usar en la construcción del sistema cuasiperiódico por el GDM.
function quasiperiodic_Neighbourhood(NSides, Error_Margin, Radius, APoint; Alpha_Value = 0.0)
    Star_Vectors = [[BigFloat(1),0]]; #Arrangement that will contain the star vectors
    for i in 1:(NSides-1)
        push!(Star_Vectors, [cos((2*i)*pi/NSides), sin((2*i)*pi/NSides)]); #Vertices of the polygon with "NSides" sides
    end
    Alphas_Array = fill(Alpha_Value, NSides); #Array with the alpha constants of the GDM
    Average_Distance_Stripes = fill(NSides/2, NSides); #Array with the average spacing between stripes

    #We generate the local neighborhood around the Test_Point point
    Dual_Points = region_Local_Voronoi(Error_Margin, Average_Distance_Stripes, Star_Vectors, Alphas_Array, APoint);

    #Let's get the coordinates as tuples of the centroids and the dictionary that relates the centroid's coordinates with the polygons vertices' coordinates of the polygon that generate the centroid.
    Centroids, Centroids_Dictionary = centroides(Dual_Points);

    Cluster_Centroids = []; #Array that will contain the coordinates of the centroid of the tiles inside the circular region
    for e in Centroids
        if norm([e[1], e[2]] - APoint) < Radius
            push!(Cluster_Centroids, e)
        end
    end

    #Let's get the X and Y coordinates of the vertices of the retained polygons in the quasiperiodic lattice
    X, Y = centroides_A_Vertices(Cluster_Centroids, Centroids_Dictionary);

    #Remove all the copies of a same vertex and keep only the lattice sites that are inside the circular neighborhood
    unique!(Dual_Points);

    Cluster_Sites = []; #Array that will contain only the lattice sites that are inside the circular neighbourhood
    for e in Dual_Points
        if norm(e - APoint) < Radius
            push!(Cluster_Sites, e)
        end
    end

    return X, Y, Cluster_Sites
end
########################################################################################################################################################################
################################################################################ Randomizar Teselados  #################################################################
########################################################################################################################################################################
function ina(x, y)
    for ys in y
        if norm(x .- ys) ≤ 1e-8
            return true
        end
    end
    return false
end

function inarhombous(x, rhombi)
    n = length(rhombi);
    rombos = [];
    for i in 1:n
        test = ina(x, rhombi[i]);
        if test
            push!(rombos, i);
        end
    end
    return rombos
end

function angulo(x, rhombus)
    i = findfirst(v -> norm(v .- x) ≤ 1e-8, rhombus);
    i₋ = mod1(i-1, 4);
    i₊ = i+1;
    v₋ = rhombus[i] .- rhombus[i₋];
    v₊ = rhombus[i₊] .- rhombus[i];
    
    θ1 = mod(atan(v₋[2], v₋[1]), 2π);
    θ2 = mod(atan(v₊[2], v₊[1]), 2π);
    θ = min(mod(π-θ1+θ2, 2π), mod(π-θ2+θ1, 2π));
    return θ, -1 .* v₋, v₊
end  

function inarhombous2(x, rhombi)
    n = length(rhombi);
    rombos = [];
    angulos = [];
    vs = [];
    for i in 1:n
        test = ina(x, rhombi[i]);
        if test
            θ, v₋, v₊ = angulo(x, rhombi[i]);
            push!(rombos, i);
            push!(angulos, θ);
            push!(vs, [v₋, v₊]);
        end
    end
    return rombos, round(180*sum(angulos)/π, digits = 5), vs
end
    
function approx(v1, v2)
    if norm(v1 .- v2) ≤ 1e-8 || norm(v1 .+ v2) ≤ 1e-8
        return true
    else
        return false
    end
end
    
function edges_vertices(v, vs, rhombi)   
    edges = [];
    centros = [v];
    vertices1 = [];
    for i in 1:3
        k = findfirst(x -> norm(v .- x) ≤ 1e-8, rhombi[i]); #Vértice del rombo i que es vértice de los otros 3 rombos
        j = findfirst(x -> !(ina(x, vs[i])), vs[mod1(i+1, 3)]);
        push!(vertices1, rhombi[i][mod1(k+2, 4)]);
        push!(edges, vs[mod1(i+1, 3)][j]);
    end
    v2 = vertices1[1] .+ edges[1];
    push!(centros, v2);
    return edges, vertices1, centros
end 

#Está horrible, lo sé, en algún momento reduciré esta función para que sea menos horrible. 
function orienta_vertices(rhombi, v)
    k1 = findfirst(x -> norm(v .- x) ≤ 1e-8, rhombi[1]);
    vertices = typeof(rhombi[1][1])[];
    push!(vertices, rhombi[1][mod1(k1+1, 4)]);
    push!(vertices, rhombi[1][mod1(k1+2, 4)]);
    push!(vertices, rhombi[1][mod1(k1+3, 4)]);
    k2 = findfirst(x -> norm(v .- x) ≤ 1e-8, rhombi[2]);
    k3 = findfirst(x -> norm(v .- x) ≤ 1e-8, rhombi[3]);
    if approx(rhombi[2][mod1(k2+1, 4)], vertices[3])
        push!(vertices, rhombi[2][mod1(k2+2, 4)]);
        push!(vertices, rhombi[2][mod1(k2+3, 4)]);
        push!(vertices, rhombi[3][mod1(k3+2, 4)]);
    elseif approx(rhombi[2][mod1(k2+3, 4)], vertices[3])
        push!(vertices, rhombi[2][mod1(k2+2, 4)]);
        push!(vertices, rhombi[2][mod1(k2+1, 4)]);
        push!(vertices, rhombi[3][mod1(k3+2, 4)]);
    elseif approx(rhombi[3][mod1(k3+1, 4)], vertices[3])
        push!(vertices, rhombi[3][mod1(k3+2, 4)]);
        push!(vertices, rhombi[3][mod1(k3+3, 4)]);
        push!(vertices, rhombi[2][mod1(k2+2, 4)]);
    elseif approx(rhombi[3][mod1(k3+3, 4)], vertices[3])
        push!(vertices, rhombi[3][mod1(k3+2, 4)]);
        push!(vertices, rhombi[3][mod1(k3+1, 4)]);
        push!(vertices, rhombi[2][mod1(k2+2, 4)]);
    end
    push!(vertices, vertices[1]);
    return vertices
end

function flip!(rhombi, sites; apoint = [0, 0.], rmax = 15)
    i = rand(1:length(sites));
    v = sites[i];
    while norm(v .- apoint) > rmax
        i = rand(1:length(sites));
        v = sites[i];
    end
    I, θs, vs = inarhombous2(v, rhombi);
    if length(I) ≠ 3 || θs ≠ 360
        return 1, v, sites, rhombi[I], rhombi
    end
    edges, vertices1, centros = edges_vertices(v, vs, rhombi[I]);
    vertices = orienta_vertices(rhombi[I], v);
    r1 = [centros[2], vertices[4], vertices[5], vertices[6], centros[2]];
    r2 = [centros[2], vertices[6], vertices[1], vertices[2], centros[2]];
    r3 = [centros[2], vertices[2], vertices[3], vertices[4], centros[2]];
    rhombi[I] = [r1, r2, r3];
    sites[i] = centros[2];
    return i, sites[i], sites, [r1, r2, r3], rhombi
end

function plot_rhombi!(rhombi, apoint, radius)
    for r in rhombi  
        centro = [sum(rs[1] for rs in r[1:4])/4, sum(rs[2] for rs in r[1:4])/4];
        rd = [0.95 .* (rs[1] - centro[1], rs[2] - centro[2]) for rs in r];
        r2 = [(rs[1] + centro[1], rs[2] + centro[2]) for rs in rd];
        a = min(1, area(r2));
        θ = asin(min(1, a));
        d = 2*sin(θ/2);
        plot!(r2, color = RGB(0.0, d/sqrt(2), d/sqrt(2)), lw = 0, fill = true, alpha = 1);
        plot!(r2, color = :black, lw = 0.5);
    end
    plot!(xlim = (apoint[1]-radius-1, apoint[1]+radius+1), ylim = (apoint[2]-radius-1, apoint[2]+radius+1));
end 

function randomiza!(sites, rhombi, paso, n, apoint, rmax; dibuja = false)
    pasos = 0;
    for k in 1:(n+1)
        if k == 1
            println("Ha comenzado el ciclo");
        end
        if k > 1
            kk, centro, sites, rs, rhombi = flip!(rhombi, sites, rmax =  rmax, apoint = apoint);
            for k2 in 1:(paso-1)
                kk, centro, sites, rs, rhombi = flip!(rhombi, sites, rmax =  rmax, apoint = apoint);
            end
        end
        if dibuja 
            plot(aspect_ratio = 1, key = false, axis = false, ticks = false, grid = false);
            plot_rhombi!(rhombi, apoint, rmax);
            plot!(title = "steps = $(pasos)", show = :ijulia);
        end
        pasos += paso;
        println("Lleva $(round(100*(k-1)/n, digits = 3)) %");
    end
    return sites, rhombi
end
########################################################################################################################################################################
################################################################################ Análisis estructural  #################################################################
########################################################################################################################################################################
#Función que calcula la escala de longitud estimada para los sistemas cuasiperiódicos perfectos de simetría rotacional "NSides".
#"NSides" es la simetría rotacional del sistema cuasiperiódico.
AproxLambda(NSides) = (2π / (1 - cos(2π / NSides)));

#Función que calcula la distribución de la distancia al vecino más cercano de los sitios que conforman a una vecindad circular del sistema cuasiperiódico.
#"Error_Margin" es el parámetro asociado al tamaño de la vecindad circular del sistema cuasiperiódico (margen de error de los enteros n_j asociados a los vectores estrella e_j).
#"SL" semilado del cuadrado centrado en el origen dentro del cual se va a seleccionar al azar el centro de la vecindad circular del sistema cuasiperiódico.
#"Bounded_Area" cota superior para el valor del área de las celdas de Voronoi generadas con los centroides de las teselas, empleado para discriminar teselas del interior de la
#vecindad circular de las que están fuera de ella.
#"Reduction_Factor" porcentaje del tamaño de la vecindad circular que se considera completamente interior, con las teselas que caen dentro de esta vecindad poseyendo todas sus
#teselas vecinas.
#"Average_Distance_Stripes" separación promedio entre las franjas cuasiperiódicas de teselas que conforman al teselado.
#"Star_Vectors" conjunto de vectores estrella con los que se construyen los teselados cuasiperiódicos.
#"Alphas_Array" conjunto de valores del parámetro alpha que determinan el desplazamiento con respecto al origen de las rectas ortogonales del mallado en el GDM.
function first_Neighbor_Distance(Error_Margin, SL, Bounded_Area, Reduction_Factor, Average_Distance_Stripes, Star_Vectors, Alphas_Array)
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
    
    #Array that will contain the index of the Voronoi's tiles associated to the vertices inside the circle.
    First_Neighbor_Distance_Array = [];

    for i in Sites_Vertices
        if Dictionary_Vertices_Inside_Circle[i] == true
            Voronoi_Index = Dictionary_Vertex_Index[i];
            Neighbors_Array = vecinos_Voronoi(Voronoi_Index, Voronoi_Vertices);
            Minimum_Distance = Inf;
            for j in Neighbors_Array
                x = norm([j[1], j[2]] - [i[1], i[2]])
                if x < Minimum_Distance
                    Minimum_Distance = x
                end
            end
            push!(First_Neighbor_Distance_Array, Minimum_Distance)
        end
    end
    
    return First_Neighbor_Distance_Array
end

#Función que calcula la distribución de la distancia a los tres vecinos más cercano de los sitios que conforman a una vecindad circular del sistema cuasiperiódico.
#"Error_Margin" es el parámetro asociado al tamaño de la vecindad circular del sistema cuasiperiódico (margen de error de los enteros n_j asociados a los vectores estrella e_j).
#"SL" semilado del cuadrado centrado en el origen dentro del cual se va a seleccionar al azar el centro de la vecindad circular del sistema cuasiperiódico.
#"Bounded_Area" cota superior para el valor del área de las celdas de Voronoi generadas con los centroides de las teselas, empleado para discriminar teselas del interior de la
#vecindad circular de las que están fuera de ella.
#"Reduction_Factor" porcentaje del tamaño de la vecindad circular que se considera completamente interior, con las teselas que caen dentro de esta vecindad poseyendo todas sus
#teselas vecinas.
#"Average_Distance_Stripes" separación promedio entre las franjas cuasiperiódicas de teselas que conforman al teselado.
#"Star_Vectors" conjunto de vectores estrella con los que se construyen los teselados cuasiperiódicos.
#"Alphas_Array" conjunto de valores del parámetro alpha que determinan el desplazamiento con respecto al origen de las rectas ortogonales del mallado en el GDM.
function third_Neighbor_Distance(Error_Margin, SL, Bounded_Area, Reduction_Factor, Average_Distance_Stripes, Star_Vectors, Alphas_Array)
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
    
    #Array that will contain the 3-tuples with the nearest neighbor distance
    Neighbor_Distance_Array = [];

    for i in Sites_Vertices
        if Dictionary_Vertices_Inside_Circle[i] == true
            Voronoi_Index = Dictionary_Vertex_Index[i];
            Neighbors_Array = vecinos_Voronoi(Voronoi_Index, Voronoi_Vertices);
            AD = []; #AD = Array_Distances
            for j in Neighbors_Array
                x = norm([j[1], j[2]] - [i[1], i[2]])
                push!(AD, x)
            end
            sort!(AD)
            push!(Neighbor_Distance_Array, (AD[1], AD[2], AD[3]))
        end
    end
    
    return Neighbor_Distance_Array
end

#Función que calcula el área de las celdas de Voronoi generadas al usar como centros de Voronoi a los sitios de la retícula cuasiperiódica.
#"Error_Margin" es el parámetro asociado al tamaño de la vecindad circular. Se relaciona con el margen de error de los enteros "Nj" de los vectores estrella "ej".
#"SL" es el tamaño del semilado del cuadrado centrado en el origen dentro del cual se genera aleatoriamente un centro para la vecindad circular.
#"Bounded_Area" es la cota superior para el área de las celdas de Voronoi generadas por los centroides de las teselas cuasiperiódicas que están dentro de la
#vecindad circular.
#"Reduction_Factor" es el porcentaje del radio de la vecindad circular que se considera completamente dentro de la misma, asegurando que los efectos de borde
#no afecten a los cálculos.
#"Average_Distance_Stripes" es el arreglo con la separación promedio entre las franjas cuasiperiódicas en el GDM.
#"Star_Vectors" es el arreglo con los vectores estrella que determinan al sistema cuasiperiódico.
#"Alphas_Array" es el arreglo con los valores de los parámetros alpha asociados a cada vector estrella, este parámetro determina el desplazamiento respecto
#al origen de las rectas ortogonales en el GDM.
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

#Función que calcula el número de sitios que caen dentro de una ventana circular como función de su radio.
#"NSides" es la simetría rotacional del sistema cuasiperiódico a generar.
#"Error_Margin" es el parámetro relacionado al tamaño de la vecindad circular, asociado con el margen de error al calcular los enteros "Nj" de los vectores "Ej".
#"ΔR" es el salto de valor entre Radio y Radio al incrementar el tamaño de la ventana circular.
#"Radius" es el radio de las vecindades circulares a generar.
#"SL" es el tamaño del semilado del cuadrado centrado en el origen dentro del cual se genera aleatoriamente un centro para la vecindad circular.
#"Iteraciones" es el número de vecindades circulares que vamos a considerar para calcular la estadística de N(R).
#"Areas_Values" es el arreglo con las diferentes áreas, redondeadas a 6 dígitos en Float64, de las teselas del teselado cuasiperiódico, ORDENADAS de MENOR a MAYOR.
#"Progress_File_Name" es una cadena de caracteres con el nombre a usar para verificar el progreso de los cálculos.
#"Data_Name_CPAG" es una cadena de caracteres con el nombre a usar para guardar los datos asociados al caso de Decoración en Centroides separados por área.
#"Data_Name_CAT" es una cadena de caracteres con el nombre a usar para guardar los datos asociados al caso de Decoración en Centroides todas las teselas.
#"Data_Name_CATwo1" es una cadena de caracteres con el nombre a usar para guardar los datos asociados al caso de Decoración en Centroides todas las teselas, excepto la menor.
#"Data_Name_CATwo2" es una cadena de caracteres con el nombre a usar para guardar los datos asociados al caso de Decoración en Centroides todas las teselas, excepto la segunda menor.
#"Data_Name_Vertices" es una cadena de caracteres con el nombre a usar para guardar los datos asociados al caso de Decoración en Vértices de las teselas.
#"Alpha_Value" es el valor del parámetro alpha asociado al desplazamiento de las rectas ortogonales en el GDM.

#*Nb es el control interno de Notebooks o más en general, hilo, que se emplea.
#*Progress_File_Name sugerido = "Estado_NR_N$(NSides)_Alfa0P0_$(Nb).txt"
#*Data_Name_CPAG sugerido = "NR_CentroidsDeco_N$(NSides)_Alfa0P0_R848_Step1_Notebook$(Nb)_Area$(AreaTile)_$(Iteracion).csv"
#*Data_Name_CAT sugerido = "NR_AllTilesCentroid_N$(NSides)_Alfa0P0_R848_Step1_Notebook$(Nb)_$(Iteracion).csv"
#*Data_Name_CATwo1 sugerido = "NR_AllTilesExceptSmallestCentroid_N$(NSides)_Alfa0P0_R848_Step1_Notebook$(Nb)_$(Iteracion).csv"
#*Data_Name_CATwo2 sugerido = "NR_AllTilesExcept2ndSmallestCentroid_N$(NSides)_Alfa0P0_R848_Step1_Notebook$(Nb)_$(Iteracion).csv"
#*Data_Name_Vertices sugerido = "NR_VerticesDeco_N$(NSides)_Alfa0P0_R848_Step1_Notebook$(Nb)_$(Iteracion).csv"
function number_Sites_Circular_Window(NSides, Error_Margin, ΔR, Radius, SL, Iteraciones, Areas_Values,
                                      Progress_File_Name, DataFolderPath, Data_Name_CPAG, Data_Name_CAT, Data_Name_CATwo1, Data_Name_CATwo2, Data_Name_Vertices;
                                      Alpha_Value = 0.0)
    Star_Vectors = [[BigFloat(1),0]]; #Arreglo que contendrá a los vectores estrella
    for i in 1:(NSides-1)
        push!(Star_Vectors, [cos((2*i)*pi/NSides), sin((2*i)*pi/NSides)]); #Vertices de un polígono regular con NSides lados
    end
    Alphas_Array = fill(Alpha_Value, NSides); #Arreglo con los valores del parámetro alpha definido en el GDM
    Average_Distance_Stripes = fill(NSides/2, NSides); #Arreglo con la separación promedio entre las franjas cuasiperiódicas

    R_Array = ΔR:ΔR:Radius #Arreglo con los diferentes radios de la ventana circular a considerar
    for Iteracion in 1:Iteraciones
        #Escribimos un archivo .txt que lleve un conteo del avance de la función
        open(Progress_File_Name, "w") do file
            write(file, "Estamos generando el cluster $(Iteracion).")
        end

        APoint = punto_Arbitrario(SL); #Elegimos un sitio arbitrario para el centro de la vecindad cuasiperiódica

        #Generamos los sitios de la retícula cuasiperiódica alrededor del punto APoint
        Dual_Points = region_Local_Voronoi(Error_Margin, Average_Distance_Stripes, Star_Vectors, Alphas_Array, APoint);
        
        ##################################ARREGLOS DECORACIÓN EN CENTROIDES#############################################
        #Centroids Matrix Value X
        CMVPAG = [[] for i in 1:length(Areas_Values)]; #Arreglo centroides separadas por área de teselas
        CMVAT = []; #Arreglo centroides todas las teselas
        CMVATwo1 = []; #Arreglo centroides todas las teselas excepto la más pequeña
        CMVATwo2 = []; #Arreglo centroides todas las teselas excepto la segunda más pequeña
        
        for i in 1:4:length(Dual_Points)
            #Arreglos con los vértices de una tesela cerrada, esto es V1->V2->V3->V4->V1, separada en X y Y.
            XX = [Dual_Points[i][1], Dual_Points[i+1][1], Dual_Points[i+2][1], Dual_Points[i+3][1], Dual_Points[i][1]];
            YY = [Dual_Points[i][2], Dual_Points[i+1][2], Dual_Points[i+2][2], Dual_Points[i+3][2], Dual_Points[i][2]];

            A = round(Float64(Area(XX, YY)), digits = 6); #Area de la tesela redondeada a seis dígitos
            Index = findfirst(x -> x == A, Areas_Values); #Índice del correspondiente valor del área, de la menor área a la mayor área
            
            #[Cx, Cy] arreglo con las coordenadas X y Y del centroide del rombo
            Centroid = [(Dual_Points[i][1] + Dual_Points[i+1][1] + Dual_Points[i+2][1] + Dual_Points[i+3][1])/4,
                        (Dual_Points[i][2] + Dual_Points[i+1][2] + Dual_Points[i+2][2] + Dual_Points[i+3][2])/4]
            
            #Enviamos a la lista de centroides correspondientes si están dentro de la vecindad circular
            if norm(Centroid - APoint) < Radius
                push!(CMVPAG[Index], Centroid) #Centroide separada por área
                push!(CMVAT, Centroid) #Centroide todas las teselas
                
                if A != Areas_Values[1]
                    push!(CMVATwo1, Centroid) #Centroide todas las teselas excepto la más pequeña
                end
                
                if A != Areas_Values[2]
                    push!(CMVATwo2, Centroid) #Centroide todas las teselas excepto la segunda más pequeña
                end
            end
        end

        ##################################ARREGLOS DECORACIÓN EN VÉRTICES#############################################
        unique!(Dual_Points); #Eliminamos las múltiples copias de un mismo vértice en la retícula cuasiperiódica
        
        #Vertices Matrix Value
        VMV = []; #Arreglo vertices únicos de todas las teselas
        
        #Enviamos a la lista de vertices los puntos que se encuentran dentro de clúster principal
        for Vertex in Dual_Points
            if norm(Vertex - APoint) < Radius
                push!(VMV, Vertex)
            end
        end
        
        Dual_Points = nothing; #Eliminamos la lista de Dual_Points para liberar memoria de la computadora
    
        ##################################CÁLCULOS DE N(R) CASO CENTROIDES###########################################
        #Cálculo de las NR por área de teselas
        for AreaTile in 1:length(Areas_Values)
            #Ordenamos la lista de centroides por su distancia al centro de la vecindad circular
            Cluster_Sites = sort(CMVPAG[AreaTile], by = x -> norm(x - APoint))

            #Generamos el arreglo que va a contener la información de N(R)
            N_R = zeros(length(R_Array));

            #Obtenemos el primer elemento de los centroides que satisface que su distancia al centro de la vecindad circular es mayor que R_Array[i]
            FirstIndex = 1;
            for i in 1:length(R_Array)
                for j in FirstIndex:length(Cluster_Sites)
                    if j == length(Cluster_Sites)
                        N_R[i] = j; First_Index = j;
                    else
                        (norm(Cluster_Sites[j] - APoint) > R_Array[i]) ? (N_R[i] = j-1; FirstIndex = j; break) : (nothing)
                    end
                end
            end

            writedlm(DataFolderPath * Data_Name_CPAG, N_R, ',')
        end
        
        #Cálculo de la NR para todas las teselas
        #Ordenamos la lista de centroides por su distancia al centro de la vecindad circular
        Cluster_Sites = sort(CMVAT, by = x -> norm(x - APoint))

        #Generamos el arreglo que va a contener la información de N(R)
        N_R = zeros(length(R_Array));

        #Obtenemos el primer elemento de los centroides que satisface que su distancia al centro de la vecindad circular es mayor que R_Array[i]
        FirstIndex = 1;
        for i in 1:length(R_Array)
            for j in FirstIndex:length(Cluster_Sites)
                if j == length(Cluster_Sites)
                    N_R[i] = j; First_Index = j;
                else
                    (norm(Cluster_Sites[j] - APoint) > R_Array[i]) ? (N_R[i] = j-1; FirstIndex = j; break) : (nothing)
                end
            end
        end

        writedlm(DataFolderPath * Data_Name_CAT, N_R, ',')
        
        #Cálculo de la NR para todas las teselas excepto la más pequeña
        #Ordenamos la lista de centroides por su distancia al centro de la vecindad circular
        Cluster_Sites = sort(CMVATwo1, by = x -> norm(x - APoint))

        #Generamos el arreglo que va a contener la información de N(R)
        N_R = zeros(length(R_Array));

        #Obtenemos el primer elemento de los centroides que satisface que su distancia al centro de la vecindad circular es mayor que R_Array[i]
        FirstIndex = 1;
        for i in 1:length(R_Array)
            for j in FirstIndex:length(Cluster_Sites)
                if j == length(Cluster_Sites)
                    N_R[i] = j; First_Index = j;
                else
                    (norm(Cluster_Sites[j] - APoint) > R_Array[i]) ? (N_R[i] = j-1; FirstIndex = j; break) : (nothing)
                end
            end
        end

        writedlm(DataFolderPath * Data_Name_CATwo1, N_R, ',')
        
        #Cálculo de la NR para todas las teselas excepto la segunda más pequeña
        #Ordenamos la lista de centroides por su distancia al centro de la vecindad circular
        Cluster_Sites = sort(CMVATwo2, by = x -> norm(x - APoint))

        #Generamos el arreglo que va a contener la información de N(R)
        N_R = zeros(length(R_Array));

        #Obtenemos el primer elemento de los centroides que satisface que su distancia al centro de la vecindad circular es mayor que R_Array[i]
        FirstIndex = 1;
        for i in 1:length(R_Array)
            for j in FirstIndex:length(Cluster_Sites)
                if j == length(Cluster_Sites)
                    N_R[i] = j; First_Index = j;
                else
                    (norm(Cluster_Sites[j] - APoint) > R_Array[i]) ? (N_R[i] = j-1; FirstIndex = j; break) : (nothing)
                end
            end
        end

        writedlm(DataFolderPath * Data_Name_CATwo2, N_R, ',')

        #########################################CÁLCULOS DE N(R) CASO VÉRTICES###########################################
        #Cálculo de la NR para los vértices de todas las teselas
        #Ordenamos la lista de centroides por su distancia al centro de la vecindad circular
        Cluster_Sites = sort(VMV, by = x -> norm(x - APoint))
        
        #Generamos el arreglo que va a contener la información de N(R)
        N_R = zeros(length(R_Array));
        
        #Obtenemos el primer elemento de los centroides que satisface que su distancia al centro de la vecindad circular es mayor que R_Array[i]
        FirstIndex = 1;
        for i in 1:length(R_Array)
            for j in FirstIndex:length(Cluster_Sites)
                if j == length(Cluster_Sites)
                    N_R[i] = j; First_Index = j;
                else
                    (norm(Cluster_Sites[j] - APoint) > R_Array[i]) ? (N_R[i] = j-1; FirstIndex = j; break) : (nothing)
                end
            end
        end
        
        writedlm(DataFolderPath * Data_Name_Vertices, N_R, ',')
    end
end
########################################################################################################################################################################
########################################################################### Conexiones Aristas Teselados ###############################################################
########################################################################################################################################################################

########################################################################################################################################################################
############################################################################## Definición de polígonos #################################################################
########################################################################################################################################################################
#Estructura para trabajar con las aristas de los polígonos.
#"inicio" es el arreglo de coordenadas [x1, y1] del vértice inicial de la arista.
#"fin" es el arreglo de coordenadas [x2, y2] del vértice final de la arista.
struct segmento
    inicio; #Arreglo [X1,Y1]
    fin; #Arreglo [X2,Y2]
end

#Función que calcula la intersección o no de dos segmentos a partir del cálculo de la intersección de dos rectas.
#"S1" es el primero de los segmentos a intersectar.
#"S2" es el segundo de los segmentos a intersectar.
function interseccion(S1::segmento, S2::segmento)
    m1 = (S1.fin[2] - S1.inicio[2])/(S1.fin[1] - S1.inicio[1]); #Pendiente del primer segmento
    m2 = (S2.fin[2] - S2.inicio[2])/(S2.fin[1] - S2.inicio[1]); #Pendiente del segundo segmento
    b1 = S1.inicio[2] - m1*S1.inicio[1]; #Ordenada al origen del primer segmento
    b2 = S2.inicio[2] - m2*S2.inicio[1]; #Ordenada al origen del segundo segmento
    A = [-m1 1; -m2 1]; b = [b1; b2];
    X = inv(A)*b; #Punto de intersección de las dos rectas definidas por los segmentos S1 y S2
    
    #Verificación de que el punto de intersección X cae, en su primer coordenada, "sobre" el segmento S1
    if (X[1] ≤ S1.fin[1]) && (X[1] ≥ S1.inicio[1]) 
        coordx1 = true
    elseif (X[1] ≥ S1.fin[1]) && (X[1] ≤ S1.inicio[1])
        coordx1 = true
    else
        coordx1 = false
    end
    
    #Verificación de que el punto de intersección X cae, en su primer coordenada, "sobre" el segmento S2
    if (X[1] ≤ S2.fin[1]) && (X[1] ≥ S2.inicio[1])
        coordx2 = true
    elseif (X[1] ≥ S2.fin[1]) && (X[1] ≤ S2.inicio[1])
        coordx2 = true
    else
        coordx2 = false
    end
    
    #Verificación de que el punto de intersección X cae, en su segunda coordenada, "sobre" el segmento S1
    if (X[2] ≤ S1.fin[2]) && (X[2] ≥ S1.inicio[2]) 
        coordy1 = true
    elseif (X[2] ≥ S1.fin[2]) && (X[2] ≤ S1.inicio[2])
        coordy1 = true
    else
        coordy1 = false
    end
    
    #Verificación de que el punto de intersección X cae, en su segunda coordenada, "sobre" el segmento S2
    if (X[2] ≤ S2.fin[2]) && (X[2] ≥ S2.inicio[2])
        coordy2 = true
    elseif (X[2] ≥ S2.fin[2]) && (X[2] ≤ S2.inicio[2])
        coordy2 = true
    else
        coordy2 = false
    end
    
    #Verificación de que el punto de intersección X cae, en su primer coordenada, "sobre" los segmentos S1 y S2
    if (coordx1 == true) && (coordx2 == true)
        coordx = true
    else
        coordx = false
    end
    
    #Verificación de que el punto de intersección X cae, en su segunda coordenada, "sobre" los segmentos S1 y S2
    if (coordy1 == true) && (coordy2 == true)
        coordy = true
    else
        coordy = false
    end
    
    #Verificación de que el punto de intersección X cae "sobre" los segmentos S1 y S2
    if (coordx == true) && (coordy == true)
        return true
    else
         return false
    end
end

#La función regresa "true" si el punto está dentro del polígono y "false" si no.
#"Poligono" es un arreglo con los segmentos que conforman al polígono en cuestión.
#"Punto" es el punto que queremos checar si está dentro o fuera. Es un arreglo con dos entradas.
function dentro(Poligono, Punto)
    #Generamos una recta con dirección "arbitraria".
    Recta = segmento(Punto, [1e6*cos(rand(0:1e-6:2*pi)),1e6*sin(rand(0:1e-6:2*pi))]);
    
    Contenedor = [];
    
    for i in 1:length(Poligono)
        if interseccion(Recta, Poligono[i])
            push!(Contenedor, true)
        end
    end
    
    #Calculemos el número de intersecciones que hay entre los segmentos del arreglo prueba.
    #Si el número de intersecciones es cero, manda "false"
    if length(Contenedor) == 0
        return false;
    elseif iseven(length(Contenedor))
        return false;
    else
        return true;
    end
end

#Función que regresa los segmentos que conforman cada uno de los polígonos del arreglo cuasiperiódico.
#"Coordenadas_X" es un arreglo con las coordenadas en X de los vértices de los polígonos. Cada 4 entradas corresponden a un mismo polígono.
#"Coordenadas_Y" es un arreglo con las coordenadas en Y de los vértices de los polígonos. Cada 4 entradas corresponden a un mismo polígono.
function obtener_Segmentos_Vertices(Coordenadas_X, Coordenadas_Y)
    #Arreglo donde irán los segmentos asociados a cada polígono
    Poligonos = [];
        
    for i in 1:4:length(Coordenadas_X)
        Segmento1 = segmento([Coordenadas_X[i], Coordenadas_Y[i]], [Coordenadas_X[i+1], Coordenadas_Y[i+1]]);
        Segmento2 = segmento([Coordenadas_X[i+1], Coordenadas_Y[i+1]], [Coordenadas_X[i+2], Coordenadas_Y[i+2]]);
        Segmento3 = segmento([Coordenadas_X[i+2], Coordenadas_Y[i+2]], [Coordenadas_X[i+3], Coordenadas_Y[i+3]]);
        Segmento4 = segmento([Coordenadas_X[i+3], Coordenadas_Y[i+3]], [Coordenadas_X[i], Coordenadas_Y[i]]);
        
        push!(Poligonos, [Segmento1, Segmento2, Segmento3, Segmento4]);
    end
    
    return Poligonos
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
########################################################################################################################################################################
####################################################################### Celda contenedora y primeros vecinos ###########################################################
########################################################################################################################################################################
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

#Definimos una función que nos genere un diccionario de las duplas de los centroides en la lista dada por el usuario al índice del polígono asociado 
#en la red de Voronoi.
#Arreglo_Duplas_Centroides: Arreglo con las duplas (X,Y) asociadas a las coordenadas de los centroides de los polígonos del arreglo cuasiperiódico
#Voronoi: Estructura de datos generada por el algoritmo para los polígonos de Voronoi
function diccionario_Centroides_Indice_Voronoi(Arreglo_Duplas_Centroides, Voronoi)
    Diccionario_Centroides_Indice = Dict(); #Definimos el diccionario vacío que iremos llenando

    for i in 1:length(Arreglo_Duplas_Centroides)
        Diccionario_Centroides_Indice[Voronoi.faces[i].site] = i; #Añadimos la llave de la dupla del centroide y le asociamos su índice en los polígonos de Voronoi
    end

    return Diccionario_Centroides_Indice
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
########################################################################################################################################################################
################################################################################### Area Teselas #######################################################################
########################################################################################################################################################################
#Función que calcula el área de un cuadrilátero a partir de las coordenadas X y Y de sus vértices.
#"XX" es el arreglo con las coordenadas X de sus cuatro vértices, repitiendo la información del primer vértice al final del arreglo.
#"YY" es el arreglo con las coordenadas Y de sus cuatro vértices, repitiendo la información del primer vértice al final del arreglo.
function Area(XX,YY)
    A = 0 #Initial area of a tile set as zero
    for i in 1:4
        A += XX[i]*YY[i+1]-YY[i]*XX[i+1]  
    end
    A *= 1/2
    return abs(A)
end

#Función que obtiene el área de todas las teselas que conforman a una vecindad del teselado cuasiperiódico, posteriormente mantiene únicamente los valores únicos de
#dichas áreas y las almacena en un arreglo ordenado de menor a mayor área.
#"X_Tiles_Coord" es un arreglo con las coordenadas X de las teselas, cada cuatro valores corresponden a una misma tesela.
#"Y_Tiles_Coord" es un arreglo con las coordenadas Y de las teselas, cada cuatro valores corresponden a una misma tesela.
function areas_Tiles(X_Tiles_Coord, Y_Tiles_Coord)
    Area_Tiles_Array = []; #Array that will contain all the different areas values for the tiles of the quasiperiodic tiling

    for i in 1:4:length(X_Tiles_Coord)
        XX = [X_Tiles_Coord[i],X_Tiles_Coord[i+1],X_Tiles_Coord[i+2],X_Tiles_Coord[i+3],X_Tiles_Coord[i]];
        YY = [Y_Tiles_Coord[i],Y_Tiles_Coord[i+1],Y_Tiles_Coord[i+2],Y_Tiles_Coord[i+3],Y_Tiles_Coord[i]];
        
        A = Area(XX, YY);
        push!(Area_Tiles_Array, A)
    end

    unique!(Area_Tiles_Array); #Kepp only the different values for the areas
    sort!(Area_Tiles_Array); #Sort the area values from smallest to biggest

    return Area_Tiles_Array
end
########################################################################################################################################################################
############################################################################### Lados Celdas Voronoi ###################################################################
########################################################################################################################################################################
#Función que nos regresa el número de lados de los polígonos de Voronoi de los clúster principales (antes de remover las capas externas para evitar deformar los 
#polígonos cercanos a la última capa del cluster). La función se itera un cierto número de veces para tener una mayor estadística.
#"Iteraciones" es el número de clúster principales que vamos a considerar para el número de lados.
#"Area_Cota" es el área que servirá como discriminante para los polígonos de la frontera.
#"Margen_Error" es el parámetro que determina el margen de error considerado al generar los posibles números enteros asociados a los vectores estrella ei, ej.
#"Semilado_Caja" es el semilado de un cuadrado en el cual se genera un punto arbitrario, alrededor del cual se genera la vecindad del arreglo cuasiperiódico.
#"Promedios_Distancia" es el arreglo con la separación entre las franjas cuasiperiódicas.
#"Vectores_Estrella" es el arreglo con los vectores estrella que generan la retícula cuasiperiódica deseada.
#"Arreglo_Alfas" es el arreglo con los parámetros de la separación respecto al origen del conjunto de rectas ortogonales a los vectores estrella.
function numero_Lados_Poligonos(Iteraciones, Area_Cota, Margen_Error, Semilado_Caja, Promedios_Distancia, Vectores_Estrella, Arreglo_Alfas)
    Arreglo_Numero_Lados_Total = []; #Arreglo donde irá el número de lados de los polígonos tras todas las iteraciones
    Arreglo_Arreglos_Vertices = []; #Arreglo que contendrá los arreglos con los vértices de cada polígono
    
    for i in 1:Iteraciones
        #Arreglo donde irán los vértices de los polígonos que viven en el cluster principal de cada una de las vecindades generadas
        Arreglo_Poligonos_Cluster_Principal = [];

        #Generamos un punto arbitrario dentro de un cuadrado centrado en el origen de semilado "Semilado_Caja"
        APoint = punto_Arbitrario(Semilado_Caja)

        #Generación del arreglo cuasiperiódico.
        Puntos_Duales = region_Local_Voronoi(Margen_Error, Promedios_Distancia, Vectores_Estrella, Arreglo_Alfas, APoint);

        #Obtención de los centroides de los poligonos del arreglo cuasiperiódico en: 
        Centroides, Diccionario_Centroides = centroides(Puntos_Duales);

        #Definimos las duplas requeridas por el algoritmo "Voronoi" de Enrique con las coordenadas de los centroides iniciales
        sites_Inicial = [(Centroides[i][1], Centroides[i][2]) for i in 1:length(Centroides)];

        #Obtenemos el arreglo de Voronoi con los centroides iniciales
        voronoi_Inicial = getVoronoiDiagram(sites_Inicial);
        
        #Nos quedamos únicamente con los centroides correspondientes a los polígonos del cluster principal
        Centroides_Cluster_Principal = centroides_Area_Acotada(voronoi_Inicial, Area_Cota, APoint);

        #Obtenemos el diccionario "Centroides -> Indices Poligonos Voronoi" de la configuración inicial
        Diccionario_Centroides_Indice = diccionario_Centroides_Indice_Voronoi(sites_Inicial, voronoi_Inicial)

        for j in Centroides_Cluster_Principal #Checamos únicamente los polígonos del cluster principal
            Face = voronoi_Inicial.faces[Diccionario_Centroides_Indice[j]]; #Polígono del clúster principal que estamos evaluando
            Arreglo_Coordenadas_Vertices = []; #Arreglo donde iremos poniendo las coordenadas de los vértices de los pol. Voronoi
            Halfedge = Face.outerComponent; #Elemento de los pol. Voronoi, si es igual a "nothing" el pol. no se cierra
            Vertice_Inicial = Halfedge.origin.coordinates; #Coordenadas del primer vértice
            push!(Arreglo_Coordenadas_Vertices, Vertice_Inicial)

            while Halfedge.next != Face.outerComponent
                Halfedge = Halfedge.next; #Nos movemos al siguiente elemento del polígono de Voronoi
                Vertice = Halfedge.origin.coordinates; #Coordenadas del vértice

                if (Vertice[1] - Arreglo_Coordenadas_Vertices[end][1])^2 + (Vertice[2] - Arreglo_Coordenadas_Vertices[end][2])^2 > 1e-6 && (Vertice[1] - Arreglo_Coordenadas_Vertices[1][1])^2 + (Vertice[2] - Arreglo_Coordenadas_Vertices[1][2])^2 > 1e-6
                    push!(Arreglo_Coordenadas_Vertices, Vertice);
                end 
            end

            push!(Arreglo_Numero_Lados_Total, length(Arreglo_Coordenadas_Vertices));
            push!(Arreglo_Poligonos_Cluster_Principal, Arreglo_Coordenadas_Vertices);
        end
        
        #Agregamos al arreglo donde irán todos los vértices de los polígonos de los clusteres principales, el arreglo con los polígonos del cluster principal de
        #esta iteración
        push!(Arreglo_Arreglos_Vertices, Arreglo_Poligonos_Cluster_Principal);
        
        println("Se han realizado $(i) iteraciones.")
    end

    return Arreglo_Numero_Lados_Total, Arreglo_Arreglos_Vertices
end
########################################################################################################################################################################
############################################################### Algoritmo_Barrido.jl ##############################################################
########################################################################################################################################################################
#Función que genera, dado un arreglo de centros, un arreglo con los polos asociados a cada obstáculo y su centro.
#Centros: Arreglo con las coordenadas (X,Y) de los centros de los obstáculos a considerar.
#r: Radio de los obstáculos.
#n: número de obstáculos que se están considerando.
function genera_Cola(Centros::Array, r::Float64, n::Int64)
    Q = []; #Arreglo donde irá la información de los polos y el centro de los obstáculos.
    for i in 1:n
        push!(Q, [Centros[i][1]+r, Centros[i][2], -1, i]); #Polo de la derecha
        push!(Q, [Centros[i][1]-r, Centros[i][2], 1, i]); #Polo de la izquierda
        push!(Q, [Centros[i][1], Centros[i][2], 0, i]); #Centro
    end
    sort!(Q); #Ordenamos el arreglo con la info
end

#Función que calcula el número de centros que caen, dado un radio, dentro de cada obstáculo. (VMU = Vecinos_Menores_R)
#Centros: Arreglo con las coordenadas (X,Y) de los centros de los obstáculos a considerar.
#r: Radio de los obstáculos.
function numero_VMR(Centros::Array, r::Float64)
    n = length(Centros); #Número de obstáculos a considerar
    Q = genera_Cola(Centros, r, n); #Puntos relevantes al análisis
    V = zeros(n); #Arreglo con contadores de centros que caen dentro para cada uno de los obstáculos
    A = Int[]; #Arreglo que contendrá los índices de los obstáculos activos
    for j in 1:length(Q)
        Q1 = Q[j]; #Tomamos cada uno de los puntos importantes tras ser ordenados respecto al eje X
        if Q1[3] == 1 #El punto es un polo de la izquierda, con lo cual activamos dicho obstáculo
            push!(A, Int(Q1[4])); #Agregamos el índice asociado al obstáculo que vamos a activar al arreglo A
        elseif Q1[3] == -1 #El punto es un polo de la derecha, con lo cual desactivamos dicho obstáculo
            filter!(e -> e ≠ Q1[4], A); #Quitamos el índice asociado al obstáculo que vamos a desactivar del arreglo A
        elseif Q1[3] == 0 #El punto es un centro, con lo cual revisamos si los obstáculos activos tienen su centro dentro
            for i in A #Revisamos con cada uno de los obstáculos activos
                #Si el centro del obstáculo activo dista menos que "r" del centro que estamos considerando, entonces sumamos un "1" al contador del obstáculo
                if 0 < norm([Centros[i][1], Centros[i][2]] - [Centros[Int(Q1[4])][1], Centros[Int(Q1[4])][2]]) < r
                    V[i] += 1; #Sumamos un uno al contador asociado al obstáculo
                end
            end
        end
    end 
    return V #Regresa el número de centros que caen dentro de un obstáculo dado un radio "r"
end
########################################################################################################################################################################
####################################################################### M I S C E L A N E A / O U T D A T E D ##########################################################
########################################################################################################################################################################
#Función que determina, fijando los vectores estrella Ej y Ek, y para algún conjunto de constantes alfa, los puntos de la retícula en el espacio real corriendo los valores enteros 
#Nj y Nk desde -N hasta N.
#"J" y "K" son los índices de los vectores estrella a considerar.
#"N" es el número entero que determina el rango [-N, N] en el que se barren los enteros de las rectas ortogonales a los vectores Ej y Ek.
#"Vectores_Estrella" es el conjunto de vectores estrella.
#"Arreglo_Alfas" es el arreglo con los valores numéricos de la separación respecto al origen del conjunto de rectas ortogonales a los vectores estrella en el método generalizado dual.
function puntos_Dual_JK(J, K, N, Vectores_Estrella, Arreglo_Alfas)
    #Paso 1: Definimos el arreglo que contendra las coordenadas [X,Y] de cada vértice de la retícula cuasiperiódica.
    Puntos_Red_Dual = [];
    
    #Paso 2: Barrermos el conjunto de números enteros Nj y Nk para formar las intersecciones de las rectas ortogonales en el mallado.
    for Nj in -N:N
        for Nk in -N:N
            #Dejamos que el try---catch detecte el error que surje cuando los vectores considerados son colineales.
            try
                #Salida: [X,Y]
                t0, t1, t2, t3 = cuatro_Regiones(J, K, Nj, Nk, Vectores_Estrella, Arreglo_Alfas)
                push!(Puntos_Red_Dual, t0);
                push!(Puntos_Red_Dual, t1);
                push!(Puntos_Red_Dual, t2);
                push!(Puntos_Red_Dual, t3);
            catch
                nothing
            end
        end
    end
    
    return Puntos_Red_Dual
end

#Función que genera una vecindad alrededor del origen de una retícula cuasiperiódica por el método generalizado dual.
#"N" es el número entero que determina el rango [-N, N] en el que se barren los enteros de las rectas ortogonales a los vectores Ej y Ek.
#"Vectores_Estrella" es el conjunto de vectores estrella.
#"Arreglo_Alfas" es el arreglo con los valores numéricos de la separación respecto al origen del conjunto de rectas ortogonales a los vectores estrella en el método generalizado dual.
function puntos_Dual(N, Vectores_Estrella, Arreglo_Alfas)
    #Paso 1: Definimos el arreglo que contendrá todos los vértices de la retícula cuasiperiódica generada.
    Puntos_Duales = [];
    
    #Paso 2: Consideramos todos las posibles parejas de vectores Ej y Ek.
    for J in 1:length(Vectores_Estrella)
        for K in (J+1):length(Vectores_Estrella)
            #Salida: [[X,Y]]
            Puntos_JK = puntos_Dual_JK(J, K, N, Vectores_Estrella, Arreglo_Alfas)
            push!(Puntos_Duales, Puntos_JK);
        end
    end
    
    #Se utiliza la función flatten para que convierta un arreglo de arreglos en un sólo arreglo grande.
    return collect(Iterators.flatten(Puntos_Duales))
end

#La función determina, fijando los vectores estrella Ej y Ek, y para algún conjunto de alfas (uno por cada vector 
#estrella), los puntos de la red dual dejando fijo el valor Nj y corriendo el valor enteros Nk desde -N hasta N
#"J" y "K" son los índices de los vectores estrella a considerar
#"Nj" es el número entero asociado a la recta ortogonal al vector Ej
#"N_Intervalo_K" es el número entero que determina el rango [-N, N] en el que se barren los enteros de las rectas ortogonales al vector Ek.
#"Vectores_Estrella" es el conjunto de vectores estrella
#"Arreglo_Alfas" es el conjunto con las constantes alfas asociadas a cada vector estrella
function puntos_Dual_JK_Franja(J, K, Nj, N_Intervalo_K, Vectores_Estrella, Arreglo_Alfas)
    #Definimos el arreglo que contendra los vectores asociados a cada punto
    Puntos_Red_Dual = [];
    Informacion_Puntos_Red_Dual = [];
    
    for Nk in -N_Intervalo_K:N_Intervalo_K
        try #Ponemos el Try-Catch para que no debamos separar los casos en que Ej y Ek son paralelos, el error lo maneja autom.
            t0, t1, t2, t3 = cuatro_Regiones(J, K, Nj, Nk, Vectores_Estrella, Arreglo_Alfas)
            push!(Puntos_Red_Dual, t0);
            push!(Puntos_Red_Dual, t1);
            push!(Puntos_Red_Dual, t2);
            push!(Puntos_Red_Dual, t3);
            push!(Informacion_Puntos_Red_Dual, info);
        catch
            nothing
        end
    end
    
    return Puntos_Red_Dual, Informacion_Puntos_Red_Dual
end

#Funcion que genera los polígonos obtenidos al usar el vector J (fijo) y el vector K (variable).
#"J" es el índice del vector estrella que se mantendrá fijo
#"Nj" es el número entero asociado a la recta ortogonal al vector Ej
#"N_Intervalo_K" es el número entero que determina el rango [-N, N] en el que se barren los enteros de las rectas ortogonales al vector Ek.
#"Vectores_Estrella" es el conjunto de vectores estrella
#"Arreglo_Alfas" es el conjunto con las constantes alfas asociadas a cada vector estrella
function franjas_Cuasiperiodicas(J, Nj, N_Intervalo_K, Vectores_Estrella, Arreglo_Alfas)
    #Definimos el arreglo que contendrá todos los puntos de la red Dual
    Puntos_Duales = [];
    Informacion_Duales = [];
    
    #Corramos los índices K de los vectores Ek
    for K in 1:length(Vectores_Estrella)
        if K == J
            nothing;
        else
            Puntos_JK, Informacion_JK = puntos_Dual_JK_Franja(J, K, Nj, N_Intervalo_K, Vectores_Estrella, Arreglo_Alfas)
            push!(Puntos_Duales, Puntos_JK);
            push!(Informacion_Duales, Informacion_JK);
        end
    end
    
    #Se utiliza la función flatten para que convierta un arreglo de arreglos en un sólo arreglo grande.
    return collect(Iterators.flatten(Puntos_Duales)), collect(Iterators.flatten(Informacion_Duales))
end