struct segmento
    inicio; #Arreglo [X1,Y1]
    fin; #Arreglo [X2,Y2]
end

function interseccion(S1::segmento,S2::segmento)
    
    m1 = (S1.fin[2] - S1.inicio[2])/(S1.fin[1] - S1.inicio[1]); #Pendiente del primer segmento
    m2 = (S2.fin[2] - S2.inicio[2])/(S2.fin[1] - S2.inicio[1]); #Pendiente del segundo segmento
    b1 = S1.inicio[2] - m1*S1.inicio[1]; #Ordenada al origen del primer segmento
    b2 = S2.inicio[2] - m2*S2.inicio[1]; #Ordenada al origen del segundo segmento
    A = [-m1 1; -m2 1]; b = [b1; b2];
    X = inv(A)*b;
    
    if (X[1] ≤ S1.fin[1]) && (X[1] ≥ S1.inicio[1]) 
        coordx1 = true
    elseif (X[1] ≥ S1.fin[1]) && (X[1] ≤ S1.inicio[1])
        coordx1 = true
    else
        coordx1 = false
    end
    
    if (X[1] ≤ S2.fin[1]) && (X[1] ≥ S2.inicio[1])
        coordx2 = true
    elseif (X[1] ≥ S2.fin[1]) && (X[1] ≤ S2.inicio[1])
        coordx2 = true
    else
        coordx2 = false
    end
    
    if (X[2] ≤ S1.fin[2]) && (X[2] ≥ S1.inicio[2]) 
        coordy1 = true
    elseif (X[2] ≥ S1.fin[2]) && (X[2] ≤ S1.inicio[2])
        coordy1 = true
    else
        coordy1 = false
    end
    
    if (X[2] ≤ S2.fin[2]) && (X[2] ≥ S2.inicio[2])
        coordy2 = true
    elseif (X[2] ≥ S2.fin[2]) && (X[2] ≤ S2.inicio[2])
        coordy2 = true
    else
        coordy2 = false
    end
    
    if (coordx1 == true) && (coordx2 == true)
        coordx = true
    else
        coordx = false
    end
    
    if (coordy1 == true) && (coordy2 == true)
        coordy = true
    else
        coordy = false
    end
    
    if (coordx == true) && (coordy == true)
        return true
    else
         return false
    end
    
end

#La función regresa "true" si el punto está dentro del polígono y "false" si no.
#Poligono es un arreglo con los segmentos que conforman al polígono en cuestión.
#Punto es el punto que queremos checar si está dentro o fuera.
function dentro(Poligono, Punto)
    #Generamos una recta con dirección "arbitraria". Recomiendo pasar a una distribución circular
    Recta = segmento(Punto, [100*cos(rand(0:1e-6:2*pi)),100*sin(rand(0:1e-6:2*pi))]);
    
    Contenedor = [];
    
    for i in 1:length(Poligono)
        if interseccion(Recta, Poligono[i])
            push!(Contenedor, true)
        end
    end
    
    #Calculemos el número de intersecciones que hay entre los segmentos del arreglo prueba.
    #Hay que regresar true o false, en lugar del texto, para que esto sea funcional.
    if iseven(length(Contenedor))
        return false;
    else
        return true;
    end
end

#Función que regresa los segmentos que conforman cada uno de los polígonos del arreglo cuasiperiódico.
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

#Funciones que itera sobre los distintos polígonos hasta encontrar el que contiene al punto "Punto".
function encontrar_Poligono(Punto, Poligonos, Informacion_Duales)
    #Iteremos sobre los posibles polígonos para encontrar el que contenga al punto
    for i in 1:length(Poligonos)
        #Si el polígono i-ésimo contiene al punto, regresa la info de ese polígono
        if dentro(Poligonos[i], Punto)
            return Informacion_Duales[i]
        end
        
    end
    #Si no encuentra polígono, que nos mande impresión con dicha información
    println("Error: No hay polígono que contenga al punto")
end
