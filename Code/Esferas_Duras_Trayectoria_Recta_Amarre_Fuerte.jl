####################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################
#################################################################### PARAMETROS TRAYECTORIA Y OBSTACULOS ###########################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################

#Definimos el radio del obstáculo a considerar en este caso particular
#const Radio_Obstaculo = 0.36;

####################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################
######################################################################### VALIDACION POSICION INICIAL ##############################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################

function posicion_Valida_Obstaculos_Circulares(Posicion_Prueba, Vertices_Lattice_Periferia)
    #Barremos todos los vertices que están cercanos a la posición propuesta para posición inicial de la partícula
    for i in Vertices_Lattice_Periferia
        #Checamos con cada uno de estos vertices (que son los "centros" de nuestros obstáculos) si al insertar nuestro obstáculo, la posición propuesta cae dentro
        #del obstáculo o no
        if norm(Posicion_Prueba - i) < Radio_Obstaculo #Si el condicional es cierto, la posición está dentro de un obstáculo
            return false #Regresamos como output de la función que la posición propuesta como posición válida es falsa
        end
    end
    return true #Si la posición propuesta no cae dentro de ningún obstáculo, entonces regresar que la posición propuesta como posición válida es verdadera 
end

####################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################
################################################################################# COLISION OBSTACULOS ##############################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################

#This is a function made by Atahualpa that gives the position of the collition and the time of collition of a particle and a sphere
"collision(x1,x2,r,v), dado la posicion de la particula en Rn, x1, la posición del centro de la esfera x2, el radio
de la esfera r, y la velocidad de la particula v, regresa el punto de colisión x y el tiempo que tarda en alcanzarlo
t"
function colision_Esfera_Dura(x1,x2,r,v)
    v2=dot(v,v)
    b=dot((x1-x2),v)/v2
    c=(dot((x1-x2),(x1-x2))-r^2)/v2
    if(BigFloat(b^2-c)<0.)
        return false, Inf
    end
    t=-b-sqrt(b^2-c)
    #=
    if(t<0) #La colisión ocurre a tiempos negativos (la esfera está "detrás" de la trayectoria)
        t=-b+sqrt(b^2-c) #Corrije el tiempo de vuelo para colisionar, sigue negativo
    end
    =#
    if(t<0) #Si tras corregir el tiempo de vuelo es negativo, la colisión no se da en tiempo >0
        t = Inf
    end
    
    x=v*t+x1
    return x, t
end

"velo_col(x1,x2,v,p), dado el punto de colisión x1, el centro de la esfera x2, v la velocidad de la particula y p la posición de la particula
regresa la velocidad v después de la colisión si consideramos que se refleja (sigue la ley de Snell de la reflexión)"
function velocidad_Esfera_Dura(x1,x2,v,p)
    rapidez=norm(v);
    n=(x1-x2)
    n=n/norm(n)
    vn=dot(n,v)*n
    v=v-2*vn
    v=v/norm(v)
    return v*rapidez
end

####################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################
################################################################################## COLISION SEGMENTOS ##############################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################

#Función que, dados los puntos A y B del segmento AB y los puntos C y D del segmento CD nos regresa el punto de intersección de ambas rectas.
#A,B,C,D: Arreglos de la forma (x,y).
function interseccion_Rectas(A, B, C, D)
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

#Función que calcula el punto de intersección entre la recta que describe la trayectoria de la partícula y la recta que generan dos punto.
#Posicion_Particula: Posición inicial de la partícula [X,Y].
#Velocidad_Particula: Velocidad inicial de la partícula [Vx,Vy].
#Punto1: Coordenadas [X,Y] de uno de los puntos a considerar.
#Punto2: Coordenadas [X,Y] de uno de los puntos a considerar.
########################################## FUNCION = interseccion_Particula_Celda ########################################################################
function interseccion_Particula_Celda_Trayectorias_Rectas(Posicion_Particula, Velocidad_Particula, Punto1, Punto2)
    #Determinamos el punto final de la partícula tras moverse un tiempo t = 1.
    Posicion_Particula_2 = Posicion_Particula + Velocidad_Particula;

    #Punto en el que se intersecta la trayectoria de la partícula con la recta generada por los extremos de la pared de la celda
    Punto_Interseccion = interseccion_Rectas(Posicion_Particula, Posicion_Particula_2, Punto1, Punto2);

    #Revisemos si la intersección de la trayectoria de la partícula con la recta del lado de la celda ocurre en el sentido de la trayectoria
    if dot(Punto_Interseccion - Posicion_Particula, Velocidad_Particula) < 0.0 #La colisión se da a tiempos negativos
        Tiempo_Vuelo = Inf;
        return Punto_Interseccion, Tiempo_Vuelo
    else
        Tiempo_Vuelo = norm(Punto_Interseccion - Posicion_Particula)/norm(Velocidad_Particula);
        return Punto_Interseccion, Tiempo_Vuelo
    end
end

function velocidad_Cambio_Celda_Trayectorias_Rectas(Posicion_Particula, Punto_Cambio_Celda, Velocidad_Particula)
    return Velocidad_Particula
end

####################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################
##################################################################################### CAMBIO DE CELDA ##############################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################

#Función que regresa la posición y velocidad de la partícula tras la colisión con el obstáculo, así como el tiempo requerido para la colisión. 
#Si no hay colisión, regresa la posición y velocidad inicial y un tiempo de colición infinito.
#Posicion_Inicial: Coordenadas (X,Y) de la posición de la partícula antes de la colisión.
#Velocidad_Inicial: Coordenadas (Vx, Vy) de la velocidad de la partícula antes de la colisión.
#Radio_Obstaculo: Radio de los obstáculos.
#Posicion_Obstaculo: Coordenadas (Ox, Oy) de la posición del obstáculo.
#funcion_Posicion: Función que determina la posición tras la interacción Partícula-Obstáculo.
#funcion_Velocidad: Función que determina la velocidad tras la interacción Partícula-Obstáculo.
function colision_Obstaculo(Posicion_Inicial, Velocidad_Inicial, Radio_Obstaculo, Posicion_Obstaculo, funcion_Posicion::Function, funcion_Velocidad::Function)
    Posicion_Colision, Tiempo_Colision = funcion_Posicion(Posicion_Inicial, Posicion_Obstaculo, Radio_Obstaculo, Velocidad_Inicial);
    
    if Tiempo_Colision < Inf
        Velocidad_Final = funcion_Velocidad(Posicion_Colision, Posicion_Obstaculo, Velocidad_Inicial, Posicion_Inicial);
        return Posicion_Colision, Velocidad_Final, Tiempo_Colision
    else
        return Posicion_Inicial, Velocidad_Inicial, Tiempo_Colision
    end
end

function colision_Segmento(Posicion_Inicial, Velocidad_Inicial, Arreglo_Vertices_Celda, Extremos_Entrada, Posicion_Obstaculo, Diccionario_Segmentos_Centro_Vecino, funcion_Interseccion_Particula_Celda::Function, funcion_Velocidad_Cambio_Celda::Function)
	Arreglo_Tiempos_Vuelo = Float64[]; #Arreglo donde iremos guardando los tiempos de vuelo
	Diccionario_Tiempos_Vuelo_Extremos = Dict(); #Diccionario que relaciona "tiempo -> Vertices del lado por el que entra la partícula"
	Diccionario_Tiempos_Vuelo_Punto_Colision = Dict(); #Diccionario que relaciona "tiempo -> Punto donde la partícula cambia de celda"

	#Probamos con cada uno de los segmentos de la celda contenedora
	for i in 1:(length(Arreglo_Vertices_Celda)-1)
		P1 = Arreglo_Vertices_Celda[i]; #Primer punto del segmento a considerar
		P2 = Arreglo_Vertices_Celda[i+1]; #Segundo punto del segmento a considerar

		if [P1, P2] == Extremos_Entrada || [P2, P1] == Extremos_Entrada #Si es la puerta por donde entró, no la consideramos
			nothing
        elseif norm(P1 - P2) < 1e-6 #Si el segmento es muy pequeño, es un segmento falso (errores introducidos al generar el Voronoi)
            nothing
		else
			#Obtenemos el punto de colisión de la partícula con ese segmento y el tiempo de vuelo que le toma
			Punto_Colision, Tiempo_Vuelo = funcion_Interseccion_Particula_Celda(Posicion_Inicial, Velocidad_Inicial, P1, P2);
			Tiempo_Vuelo = Float64(Tiempo_Vuelo);

        	push!(Arreglo_Tiempos_Vuelo, Tiempo_Vuelo); #Añadimos el tiempo de vuelo al arreglo

        	if Tiempo_Vuelo < Inf
            	Diccionario_Tiempos_Vuelo_Extremos[Tiempo_Vuelo] = [P1, P2]; #Relacionamos el tiempo de vuelo con la puerta de entrada a la partícula
            	Diccionario_Tiempos_Vuelo_Punto_Colision[Tiempo_Vuelo] = Punto_Colision; #Relacionamos el tiempo de vuelo con el punto donde la partícula cambia de celda
        	end
        end
    end

    #Probamos con el último de los segmentos, el formado por el vértice N y el vértice 1
    P1 = Arreglo_Vertices_Celda[end]; #Primer punto del segmento a considerar
	P2 = Arreglo_Vertices_Celda[1]; #Segundo punto del segmento a considerar

	if [P1, P2] == Extremos_Entrada || [P2, P1] == Extremos_Entrada #Si es la puerta por donde entró, no la consideramos
		nothing
    elseif norm(P1 - P2) < 1e-6 #Si el segmento es muy pequeño, es un segmento falso (errores introducidos al generar el Voronoi)
        nothing
	else
		#Obtenemos el punto de colisión de la partícula con ese segmento y el tiempo de vuelo que le toma
		Punto_Colision, Tiempo_Vuelo = funcion_Interseccion_Particula_Celda(Posicion_Inicial, Velocidad_Inicial, P1, P2);
		Tiempo_Vuelo = Float64(Tiempo_Vuelo);

        push!(Arreglo_Tiempos_Vuelo, Tiempo_Vuelo); #Añadimos el tiempo de vuelo al arreglo

        if Tiempo_Vuelo < Inf
            Diccionario_Tiempos_Vuelo_Extremos[Tiempo_Vuelo] = [P1, P2]; #Relacionamos el tiempo de vuelo con la puerta de entrada a la partícula
            Diccionario_Tiempos_Vuelo_Punto_Colision[Tiempo_Vuelo] = Punto_Colision; #Relacionamos el tiempo de vuelo con el punto donde la partícula cambia de celda
        end
    end

    #Obtenemos el mínimo de los tiempos de vuelo
    Tiempo_Vuelo_Minimo = minimum(Arreglo_Tiempos_Vuelo);

    #EL PRESENTE WHILE FUE LA SOLUCIÓN A UN ERROR QUE DURÓ MESES EN EL PROYECTO... PRESS "F" TO PAY RESPECT.
    while true
        Segmento_Salida = Diccionario_Tiempos_Vuelo_Extremos[Tiempo_Vuelo_Minimo]
        Centro_Poligono_Contenedor = Posicion_Obstaculo; #Coordenadas del centro del obstáculo del polígono contenedor
        Centro_Poligono_Receptor = [0.0, 0.0];
        try
        	Centro_Poligono_Receptor = Diccionario_Segmentos_Centro_Vecino[Segmento_Salida[1], Segmento_Salida[2]]; #Coordenadas del centro del obstáculo del polígono receptor
        catch
        	Centro_Poligono_Receptor = Diccionario_Segmentos_Centro_Vecino[Segmento_Salida[2], Segmento_Salida[1]]; #Coordenadas del centro del obstáculo del polígono receptor
        end
        Vector_Orientado_Contenedor_Receptor = [Centro_Poligono_Receptor[1] - Centro_Poligono_Contenedor[1], Centro_Poligono_Receptor[2] - Centro_Poligono_Contenedor[2]]; #Vector que va del polígono contenedor al polígono receptor

        #Calculamos la velocidad de la partícula en el punto en el cuál supuestamente sale de la celda
        Punto_Cambio_Celda = Diccionario_Tiempos_Vuelo_Punto_Colision[Tiempo_Vuelo_Minimo];
        Velocidad_Cambio_Celda = funcion_Velocidad_Cambio_Celda(Posicion_Inicial, Punto_Cambio_Celda, Velocidad_Inicial);

        if dot(Velocidad_Cambio_Celda, Vector_Orientado_Contenedor_Receptor) > 0 #Revisamos si la velocidad de la partícula coincide con el supuesto cambio de polígono (Solución al error en las esquinas)
            return Punto_Cambio_Celda, Velocidad_Cambio_Celda, Tiempo_Vuelo_Minimo, Segmento_Salida
        else
            filter!(x->x != Tiempo_Vuelo_Minimo, Arreglo_Tiempos_Vuelo); #Eliminamos el tiempo de vuelo mínimo que genera errores de la lista de tiempos de vuelo
            Tiempo_Vuelo_Minimo = minimum(Arreglo_Tiempos_Vuelo); #Volvemos a obtener el mínimo del tiempo de vuelo
        end
    end
end

function cambio_Celda_Esferas_Duras_Trayectoria_Recta(Posicion_Inicial, Velocidad_Inicial, Arreglo_Vertices_Celda, Coordenadas_Segmento_Entrada, Posicion_Obstaculo, Diccionario_Segmentos_Centro_Vecino)
	Tiempo_Acumulado_Vuelo = 0; #Variable que contendrá el tiempo de vuelo que llevamos acumulado

	Posicion_Colision, Velocidad_Colision, Tiempo_Vuelo_Colision = colision_Obstaculo(Posicion_Inicial, Velocidad_Inicial, Radio_Obstaculo, Posicion_Obstaculo, colision_Esfera_Dura, velocidad_Esfera_Dura);

	Posicion_Cambio_Celda, Velocidad_Cambio_Celda, Tiempo_Vuelo_Cambio_Celda, Extremos_Cambio_Celda = colision_Segmento(Posicion_Inicial, Velocidad_Inicial, Arreglo_Vertices_Celda, Coordenadas_Segmento_Entrada, Posicion_Obstaculo, Diccionario_Segmentos_Centro_Vecino, interseccion_Particula_Celda_Trayectorias_Rectas, velocidad_Cambio_Celda_Trayectorias_Rectas);

	if Tiempo_Vuelo_Colision < Tiempo_Vuelo_Cambio_Celda #La partícula colisiona antes de moverse de polígono
        #Actualizamos las variables
        Tiempo_Acumulado_Vuelo += Tiempo_Vuelo_Colision; #Se reduce al tiempo de vuelo restante el tiempo en colisionar
        Posicion_Inicial = Posicion_Colision; #La nueva posición de la partícula es la posición de la colisión
        Velocidad_Inicial = Velocidad_Colision; #La nueva velocidad de la partícula es la velocidad tras la colisión
                
        #Si hay colisión permitimos que la partícula salga por donde entró
        Coordenadas_Segmento_Entrada = [[Inf, Inf], [Inf, Inf]];

        #Calcula ahora la única posibilidad que queda que es salir de la celda
        Posicion_Cambio_Celda, Velocidad_Cambio_Celda, Tiempo_Vuelo_Cambio_Celda, Extremos_Cambio_Celda = colision_Segmento(Posicion_Inicial, Velocidad_Inicial, Arreglo_Vertices_Celda, Coordenadas_Segmento_Entrada, Posicion_Obstaculo, Diccionario_Segmentos_Centro_Vecino, interseccion_Particula_Celda_Trayectorias_Rectas, velocidad_Cambio_Celda_Trayectorias_Rectas);

        Tiempo_Acumulado_Vuelo += Tiempo_Vuelo_Cambio_Celda; #Actualizamos el tiempo de vuelo restante
        return Posicion_Cambio_Celda, Velocidad_Cambio_Celda, Tiempo_Acumulado_Vuelo, Extremos_Cambio_Celda, Posicion_Colision
    else #La partícula se mueve de polígono antes de colisionar con el obstáculo de su polígono contenedor
        Tiempo_Acumulado_Vuelo += Tiempo_Vuelo_Cambio_Celda; #Actualizamos el tiempo de vuelo restante
        return Posicion_Cambio_Celda, Velocidad_Cambio_Celda, Tiempo_Acumulado_Vuelo, Extremos_Cambio_Celda, [Inf, Inf]
    end
end

function estado_Tras_Tiempo_Rectas_Parcial(Posicion_Inicial, Velocidad_Inicial, Tiempo_Vuelo)
    Posicion_Final = Posicion_Inicial + Velocidad_Inicial*Tiempo_Vuelo;
    Velocidad_Final = Velocidad_Inicial;
    return Posicion_Final, Velocidad_Final
end

function estado_Tras_Tiempo_Rectas(Posicion_Inicial, Velocidad_Inicial, Arreglo_Vertices_Celda, Coordenadas_Segmento_Entrada, Posicion_Obstaculo, Diccionario_Segmentos_Centro_Vecino, Tiempo_Vuelo_Parcial)
	Posicion_Colision, Velocidad_Colision, Tiempo_Vuelo_Colision = colision_Obstaculo(Posicion_Inicial, Velocidad_Inicial, Radio_Obstaculo, Posicion_Obstaculo, colision_Esfera_Dura, velocidad_Esfera_Dura);

	Posicion_Cambio_Celda, Velocidad_Cambio_Celda, Tiempo_Vuelo_Cambio_Celda, Extremos_Cambio_Celda = colision_Segmento(Posicion_Inicial, Velocidad_Inicial, Arreglo_Vertices_Celda, Coordenadas_Segmento_Entrada, Posicion_Obstaculo, Diccionario_Segmentos_Centro_Vecino, interseccion_Particula_Celda_Trayectorias_Rectas, velocidad_Cambio_Celda_Trayectorias_Rectas);

	if Tiempo_Vuelo_Colision < Tiempo_Vuelo_Cambio_Celda #La partícula colisiona antes de moverse de polígono
        if (Tiempo_Vuelo_Parcial - Tiempo_Vuelo_Colision) > 0.0 #Hay tiempo para la colisión completa
            #Actualizamos las variables
            Tiempo_Vuelo_Parcial -= Tiempo_Vuelo_Colision; #Se reduce al tiempo de vuelo restante el tiempo en colisionar
            Posicion_Inicial = Posicion_Colision; #La nueva posición de la partícula es la posición de la colisión
            Velocidad_Inicial = Velocidad_Colision; #La nueva velocidad de la partícula es la velocidad tras la colisión
            Coordenadas_Segmento_Entrada = [[Inf, Inf], [Inf, Inf]]; #Si hay colisión permitimos que la partícula salga por donde entró

           	Posicion_Inicial, Velocidad_Inicial = estado_Tras_Tiempo_Rectas_Parcial(Posicion_Inicial, Velocidad_Inicial, Tiempo_Vuelo_Parcial);
			return Posicion_Inicial, Velocidad_Inicial, Coordenadas_Segmento_Entrada, Posicion_Colision
        else #Se agota el tiempo de vuelo antes de la colisión completa
            Posicion_Inicial, Velocidad_Inicial = estado_Tras_Tiempo_Rectas_Parcial(Posicion_Inicial, Velocidad_Inicial, Tiempo_Vuelo_Parcial);
			return Posicion_Inicial, Velocidad_Inicial, Coordenadas_Segmento_Entrada, [Inf, Inf]
        end
    else #La partícula se mueve de polígono antes de colisionar con el obstáculo de su polígono contenedor
        #No hay tiempo de vuelo restante para salir del polígono
        Posicion_Inicial, Velocidad_Inicial = estado_Tras_Tiempo_Rectas_Parcial(Posicion_Inicial, Velocidad_Inicial, Tiempo_Vuelo_Parcial);
        return Posicion_Inicial, Velocidad_Inicial, Coordenadas_Segmento_Entrada, [Inf, Inf]
    end
end
