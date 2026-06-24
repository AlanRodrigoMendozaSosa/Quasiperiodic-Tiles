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
const Radio_Obstaculo = 0.29;
#Definimos el radio de la trayectoria que sigue la partícula
const Radio_Trayectoria = 0.4;

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

function posicion_Valida_Obstaculos_Circulares(Posicion_Prueba, Vertices_Celda_Lattice)
	#Barremos todos los vertices que están cercanos a la posición propuesta para posición inicial de la partícula
	for i in Vertices_Celda_Lattice
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
############################################################################# TRAYECTORIAS CIRCULARES ##############################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################

#Función que nos regresa el centro de la circunferencia de radio "Radio", que pasa por "Posicion" y con velocidad lineal en "Posicion" igual a "Velocidad".
#Velocidad: Arreglo con coordenadas (Vx,Vy) del vector tangente a la circunferencia paralelo a la velocidad lineal de la partícula en el punto "Posicion".
#Posicion: Arreglo con coordenadas (X,Y) de la posición por la cual pasa la circunferencia deseada.
#Radio: Radio de la circunferencia que describe la partícula.
function centro_Circunferencia(Velocidad::Array, Posicion::Array, Radio)
	Velocidad = Velocidad/norm(Velocidad); #Vector unitario tangente a la circunferencia con dirección paralela a la velocidad lineal de la partícula
	Velocidad_Ortogonal = [Velocidad[2], -Velocidad[1]]; #Vector unitario que apunta, desde la posición de la partícula "Posicion", al centro de la circunferencia
	Centro = Velocidad_Ortogonal*Radio + Posicion; #Arreglo con las coordenadas del centro de la circunferencia
	return Centro
end

#Función que nos regresa el ángulo de un punto sobre una circunferencia a partir del eje X en el sentido de las manecillas del reloj
#Punto_Circulo: Arreglo con las coordenadas [X,Y] del punto sobre la circunferencia.
#Centro_Circulo: Arreglo con las coordenadas [X,Y] del centro de la circunferencia.
function angulo_Posicion_Circulo(Punto_Circulo, Centro_Circulo)
	#Posición relativa al centro del círculo
	Posicion_Relativa = Punto_Circulo - Centro_Circulo;

	#Ángulo de [-Pi, Pi] de la posición relativa respecto al eje X en el sentido de a las manecillas del reloj
	Angulo = -atan(Posicion_Relativa[2], Posicion_Relativa[1]);

	if Angulo < 0.
		Angulo += 2*pi;
	end

	return Angulo
end

#Función que calcula la diferencia angular entre dos posiciones. El tiempo de vuelo asociado a ir de la posición inicial a la final es la longitud de arco entre la rapidez de la partícula.
#"Posicion_Inicial" es un arreglo [X,Y] con las coordenadas de la posición inicial de la partícula.
#"Posicion_Final" es un arreglo [X,Y] con las coordenadas de la posición final de la partícula.
#"Velocidad" es un arreglo [Vx,Vy] con las componentes de la velocidad inicial de la partícula.
#"Centro": Arreglo con las coordenadas [X,Y] del centro de la circunferencia.
#Radio: Radio de la circunferencia que describe la partícula.
function tiempo_Trayectoria(Posicion_Inicial::Array, Posicion_Final::Array, Velocidad::Array, Centro::Array, Radio)
    Angulo_Posicion_Inicial = angulo_Posicion_Circulo(Posicion_Inicial, Centro);
    Angulo_Posicion_Final = angulo_Posicion_Circulo(Posicion_Final, Centro);

    Diferencia_Angulo = Angulo_Posicion_Final - Angulo_Posicion_Inicial;
    if Diferencia_Angulo >= 0.0
        return Diferencia_Angulo*Radio/norm(Velocidad)
    else
        return (2*pi - abs(Diferencia_Angulo))*Radio/norm(Velocidad)
    end
end

#Funcion que nos da un arreglo con las coordenadas (Vx,Vy) de la velocidad lineal de la partícula en el punto "Posicion" considerando que sigue una trayectoria circular en el sentido de las manecillas del reloj
#centrada en "Centro".
#Posicion: Arreglo con las coordenadas (X,Y) asociadas a la posición de la partícula en una trayectoria circular.
#Centro: Arreglo con las coordenadas (Cx, Cy) asociadas al centro de la circunferencia que describe la partícula al moverse.
function velocidad_Particula(Posicion::Array, Velocidad::Array, Centro::Array)
	Direccion = Posicion - Centro; #Vector que apunta del centro a la posición de la partícula
	Direccion = Direccion/norm(Direccion); #Volvemos unitario el vector previo
	Rapidez = norm(Velocidad); #Rapidez con la que se mueve la partícula
	return [Direccion[2]*Rapidez, -Direccion[1]*Rapidez] #Regresamos al usuario el vector unitario que apunta en la dirección de la velocidad lineal de la partícula
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

#Función que calcula la matriz de rotación V1 -> V2, siendo V1 un vector y V2 otro vector en un espacio 2D. De igual forma regresa la "Distancia" entre los puntos.
#V1: Arreglo con las coordenadas (X,Y) del vector V1.
#V2: Arreglo con las coordenadas (X,Y) del vector V2.
function rotacion(V1::Array, V2::Array)
	Distancia = norm(V1-V2); #Distancia entre los puntos V1 y V2

	CC = dot(V2 - V1, [1.0; 0.0])/Distancia; #Cos(Theta), siendo Theta el ángulo entre el eje X y el vector V2-V1
	SS = dot(V2 - V1, [0.0; 1.0])/Distancia; #Sin(Theta), siendo Theta el ángulo entre el eje X y el vector V2-V1

	return [[CC -SS]; [SS CC]], Distancia
end

#Función que calcula los dos puntos de intersección entre dos círculos, uno centrado en (0,0) y otro centrado en (X,0).
#Radio1: Radio del círculo centrado en (0,0).
#Radio2: Radio del círculo centrado en (X,0).
#X: Coordenada X del segundo círculo.
function interseccion_1(Radio1, Radio2, X)	
	X1 = (Radio1^2 - Radio2^2 + X^2)/(2*X); #Coordenada en X de la intersección de los dos círculos (si es que existe)
	Discriminante = -Radio1^4 + 2*(Radio1^2)*(Radio2^2) - Radio2^4 + 2*(Radio1^2)*(X^2) + 2*(Radio2^2)*(X^2) - X^4; #Valor que nos indica si los círculos se intersectan o no

	if Discriminante >= 0.0 #En caso de ser cierta la desigualdad, hay intersección de los círculos
		Y1 = -sqrt(Discriminante)/(2*X); #Coordenada en Y de la intersección de los dos círculos
		return [X1; Y1], [X1; -Y1] #Regresamos las coordenadas de los dos puntos de intersección de los círculos
	else
		return [Inf; Inf], [Inf; Inf] #Caso contrario, no hay intersección y regresamos como puntos de intersección el infinito
	end
end

#Función que calcula los dos puntos de intersección entre dos círculos, uno con centro en Centro1 y Radio1, otro con centro en Centro2 y Radio2.
#Centro1: Arreglo con las coordenadas (X1,Y1) del centro del círculo 1.
#Radio1: Radio del círculo centrado en (X1,Y1).
#Centro2: Arreglo con las coordenadas (X2,Y2) del centro del círculo 2.
#Radio2: Radio del círculo centrado en (X2,Y2).
function interseccion_2(Radio1, Radio2, Centro1::Array, Centro2::Array)
	#Solucionamos la intersección como si los círculos estuvieran en (0,0) y en (d,0).
	Matriz_Rotacion, Distancia = rotacion(Centro1, Centro2);
	Prev_Colision1, Prev_Colision2 = interseccion_1(Radio1, Radio2, Distancia);

	if Prev_Colision1[1] == Inf #Si la colisión se da en infinito, eso no se modificará tras rotar y trasladar
		return [Inf; Inf], [Inf; Inf]
	end

	#Con la matriz de rotación regresamos las soluciones al caso en el que no están centrados los círculos en el origen y en el (d,0).
	Colision1 = Matriz_Rotacion*Prev_Colision1 + Centro1;
	Colision2 = Matriz_Rotacion*Prev_Colision2 + Centro1;

	return Colision1, Colision2
end

#Función que calcula el punto de colisión entre una partícula que sigue una trayectoria circular recorrida en el sentido de las manecillas del reloj y un obstáculo circular.
#"Posicion_Inicial" es un arreglo [X,Y] con las coordenadas de la posición inicial de la partícula.
#"Posicion_Obstaculo" es un arreglo [X,Y] con las coordenadas del centro del obstáculo circular.
#"Radio_Obstaculo" es el radio de los obstáculos circulares a considerar.
#"Velocidad_Inicial" es un arreglo [Vx, Vy] con las componentes de la velocidad inicial de la partícula.
########################################## FUNCION = funcion_Posicion_Colision_Obstaculo ################################################################################
function colision_Esferas_Duras_Campo_Magnetico(Posicion_Inicial, Posicion_Obstaculo, Radio_Obstaculo, Velocidad_Inicial)
	#Calculamos el centro de la circunferencia que describe la partícula al moverse
	Centro_Trayectoria = centro_Circunferencia(Velocidad_Inicial, Posicion_Inicial, Radio_Trayectoria);

	#Obtenemos las coordenadas de las posibles colisiones entre la partícula y el obstáculo
	Colision1, Colision2 = interseccion_2(Radio_Trayectoria, Radio_Obstaculo, Centro_Trayectoria, Posicion_Obstaculo);

	if Colision2[1] == Inf || Colision2[1] == -Inf || Colision2[1] == NaN
		return [Inf; Inf], Inf #Si la colisión no ocurre, regresa al usuario que la colisión se da en infinito en un tiempo infinito
	end

	Tiempo = tiempo_Trayectoria(Posicion_Inicial, Colision2, Velocidad_Inicial, Centro_Trayectoria, Radio_Trayectoria); #Obtenemos el tiempo que tarda en producirse la colisión.
	return Colision2, Tiempo
end

#Función que calcula la velocidad lineal de la partícula tras colisionar con algún obstáculo.
#Centro_Obstaculo: Arreglo con las coordenadas (Cx, Cy) asociadas al centro del obstáculo con el que la partícula colisionó.
#Posicion_Colision: Arreglo con las coordenadas (X,Y) asociadas a la posición de la partícula al colisionar con el obstáculo.
#Velocidad_Colision: Arreglo con las coordenadas (Vx, Vy) asociadas a la velocidad de la partícula al colisionar con el objetivo.
function velocidad_Tras_Colision(Centro_Obstaculo::Array, Posicion_Colision::Array, Velocidad_Colision::Array)
	Rapidez = norm(Velocidad_Colision); #Rapidez de la partícula al desplazarse
	Normal_Colision = (Centro_Obstaculo - Posicion_Colision); #Vector que apunta de la posición "Posicion_Colision" al "Centro_Obstaculo"
	Normal_Colision = Normal_Colision/norm(Normal_Colision); #Normalizamos el vector previo
	Velocidad_Normal_Colision = dot(Normal_Colision, Velocidad_Colision)*Normal_Colision; #Componente normal a la recta tangente al obstáculo en el punto de colisión de la velocidad de la partícula al colisionar
	Velocidad_Colision = Velocidad_Colision - 2*Velocidad_Normal_Colision; #Velocidad de la partícula tras colisionar con el obstáculo
	Velocidad_Colision = Velocidad_Colision/norm(Velocidad_Colision); #Normalizamos la velocidad de la partícula tras la colisión
	Velocidad_Colision = Velocidad_Colision*Rapidez; #Nos aseguramos que la rapidez de la partícula se mantenga constante
	return Velocidad_Colision
end

#Función que calcula la velocidad tras la colisión entre una partícula que sigue una trayectoria circular recorrida en el sentido de las manecillas del reloj y un obstáculo circular.
#Posicion_Colision: Arreglo con las coordenadas (X,Y) asociadas a la posición de la partícula al colisionar con el obstáculo.
#Posicion_Obstaculo: Arreglo con las coordenadas (Cx, Cy) asociadas al centro del obstáculo con el que la partícula colisionó.
#Velocidad_Inicial: Arreglo con las componentes de la velocidad inicial de la partícula.
#Posicion_Inicial: Arreglo con las coordenadas [X,Y] de la posición inicial de la partícula.
########################################## FUNCION = funcion_Velocidad_Colision_Obstaculo ################################################################################
function velocidad_Esferas_Duras_Campo_Magnetico(Posicion_Colision, Posicion_Obstaculo, Velocidad_Inicial, Posicion_Inicial)
	#Obtenemos el centro de la circunferencia que describe la partícula al moverse
	Centro_Trayectoria = centro_Circunferencia(Velocidad_Inicial, Posicion_Inicial, Radio_Trayectoria);

	#Obtenemos la velocidad de la partícula al colisionar
	Velocidad_Colision = velocidad_Particula(Posicion_Colision, Velocidad_Inicial, Centro_Trayectoria);

	#Obtenemos la velocidad de la partícula tras colisionar
	Velocidad_Final = velocidad_Tras_Colision(Posicion_Obstaculo, Posicion_Colision, Velocidad_Colision);

	return Velocidad_Final
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

#Función que obtiene las coordenadas de los puntos donde se intersecta una recta y una esfera (si no se intersectan arroja coordenadas [Inf, Inf]).
#Punto_Inicial: Arreglo con las coordenadas [X,Y] de un punto sobre la recta.
#Centro_Obstaculo: Arreglo con las coordenadas [X,Y] del centro de la esfera.
#Radio_Obstaculo: Radio de la esfera.
#Direccion_Recta: Arreglo [X,Y] correspondiente a la dirección en la que se desplaza la recta.
function colision_Esfera_Recta(Punto_Inicial,Centro_Obstaculo,Radio_Obstaculo,Direccion_Recta)
	#Parémtros que aún no entiendo
    v2=dot(Direccion_Recta,Direccion_Recta);
    b=dot((Punto_Inicial - Centro_Obstaculo),Direccion_Recta)/v2;
    c=(dot((Punto_Inicial - Centro_Obstaculo),(Punto_Inicial - Centro_Obstaculo)) - Radio_Obstaculo^2)/v2;

    #Veamos si la recta intersecta a la esfera, si no mandamos punto de intersección a "falso" y tiempo de intersección a "infinito"
    if(BigFloat(b^2-c)<0.)
        return [Inf, Inf], [Inf, Inf]
    end

    #Obtenemos los dos puntos donde se intersectan la recta y la esfera
    t1=-b-sqrt(b^2-c);
    x1=Direccion_Recta*t1 + Punto_Inicial;

    t2=-b+sqrt(b^2-c);
    x2=Direccion_Recta*t2 + Punto_Inicial;
    
    return x1, x2
end

#Función que calcula el primer punto de intersección de una partícula que se desplaza en una trayectoria circular en el sentido de las manecillas del reloj y una recta.
#Posicion_Particula: Arreglo con las coordenadas [X,Y] de la posición inicial de la partícula.
#Velocidad_Particula: Arreglo con las componentes [Vx,Vy] de la velocidad inicial de la partícula.
#Punto1: Arreglo con las coordenadas [X,Y] de uno de los puntos de la recta.
#Punto2: Arreglo con las coordenadas [X,Y] de uno de los puntos de la recta.
################################################# FUNCION = interseccion_Particula_Celda #########################################
function interseccion_Particula_Celda_Trayectoria_Circular(Posicion_Particula, Velocidad_Particula, Punto1, Punto2)
	#Obtenemos el centro de la trayectoria que describe la partícula
	Centro_Trayectoria_Particula = centro_Circunferencia(Velocidad_Particula, Posicion_Particula, Radio_Trayectoria);

	#Obtenemos la dirección de la recta que genera el lado de la celda de Voronoi
	Direccion_Lado_Celda = Punto2 - Punto1;

	#Obtenemos los sitios donde se intersecta la recta y la trayectoria de la partícula
	Interseccion1, Interseccion2 = colision_Esfera_Recta(Punto1, Centro_Trayectoria_Particula, Radio_Trayectoria, Direccion_Lado_Celda);

	if Interseccion1[1] == Inf
		return [Inf, Inf], Inf
	end

	#Obtengamos los ángulos de las posibles intersecciones y de la posición de la partícula
	Angulo_Posicion_Particula = angulo_Posicion_Circulo(Posicion_Particula, Centro_Trayectoria_Particula);
	Angulo_Interseccion1 = angulo_Posicion_Circulo(Interseccion1, Centro_Trayectoria_Particula);
	Angulo_Interseccion2 = angulo_Posicion_Circulo(Interseccion2, Centro_Trayectoria_Particula);

	#Fijamos el cero de nuestra escala en el ángulo de la posición de la partícula y recalculamos los ángulos
	Angulo_Interseccion1 -= Angulo_Posicion_Particula;
	Angulo_Interseccion2 -= Angulo_Posicion_Particula;

	#Generamos el diccionario que asocia el ángulo con el punto en el círculo
	Diccionario_Angulo_Interseccion = Dict();
	Diccionario_Angulo_Interseccion[Angulo_Interseccion1] = Interseccion1;
	Diccionario_Angulo_Interseccion[Angulo_Interseccion2] = Interseccion2;

	#Definimos las variables que corresponderán al punto de intersección y al tiempo de vuelo
	Punto_Interseccion = [Inf, Inf];
	Tiempo_Vuelo = Inf;

	#Obtengamos ahora cuál es el punto de intersección con el que colisionará primero la partícula considerando su trayectoria circular a partir del punto inicial
	if (Angulo_Interseccion1 > 0 && Angulo_Interseccion2 > 0) || (Angulo_Interseccion1 < 0 && Angulo_Interseccion2 < 0)
		Punto_Interseccion = Diccionario_Angulo_Interseccion[minimum([Angulo_Interseccion1, Angulo_Interseccion2])];
		Tiempo_Vuelo = tiempo_Trayectoria(Posicion_Particula, Punto_Interseccion, Velocidad_Particula, Centro_Trayectoria_Particula, Radio_Trayectoria);
		#Tiempo_Vuelo = abs(minimum([Angulo_Interseccion1, Angulo_Interseccion2])*Radio_Trayectoria/norm(Velocidad_Particula)); 
		if abs(Tiempo_Vuelo) < 1e-6
			Punto_Interseccion = Diccionario_Angulo_Interseccion[maximum([Angulo_Interseccion1, Angulo_Interseccion2])];
			Tiempo_Vuelo = tiempo_Trayectoria(Posicion_Particula, Punto_Interseccion, Velocidad_Particula, Centro_Trayectoria_Particula, Radio_Trayectoria);
			#Tiempo_Vuelo = abs(maximum([Angulo_Interseccion1, Angulo_Interseccion2])*Radio_Trayectoria/norm(Velocidad_Particula));
			#println("El Tiempo Vuelo caso 1 es $(Tiempo_Vuelo)"); ##V-ONLY##
		end
	else
		Punto_Interseccion = Diccionario_Angulo_Interseccion[maximum([Angulo_Interseccion1, Angulo_Interseccion2])];
		Tiempo_Vuelo = tiempo_Trayectoria(Posicion_Particula, Punto_Interseccion, Velocidad_Particula, Centro_Trayectoria_Particula, Radio_Trayectoria);
		#Tiempo_Vuelo = abs(maximum([Angulo_Interseccion1, Angulo_Interseccion2])*Radio_Trayectoria/norm(Velocidad_Particula));
		if abs(Tiempo_Vuelo) < 1e-6
			Punto_Interseccion = Diccionario_Angulo_Interseccion[minimum([Angulo_Interseccion1, Angulo_Interseccion2])];
			Tiempo_Vuelo = tiempo_Trayectoria(Posicion_Particula, Punto_Interseccion, Velocidad_Particula, Centro_Trayectoria_Particula, Radio_Trayectoria);
			#Tiempo_Vuelo = abs(minimum([Angulo_Interseccion1, Angulo_Interseccion2])*Radio_Trayectoria/norm(Velocidad_Particula)); 
			#println("El Tiempo Vuelo caso 2 es $(Tiempo_Vuelo)"); ##V-ONLY##
		end
	end

	return Punto_Interseccion, Tiempo_Vuelo
end

#Función que calcula la velocidad de la partícula en un punto determinado de una circunferencia.
#Posicion_Particula: Arreglo con las coordenadas [X,Y] de la posición inicial de la partícula.
#Punto_Interseccion: Arreglo con las coordenadas [X,Y] de la posición final de la partícula.
#Velocidad_Particula: Arreglo con las componentes [Vx,Vy] de la velocidad inicial de la partícula.
function velocidad_Cambio_Celda_Trayectoria_Circular(Posicion_Particula, Punto_Interseccion, Velocidad_Particula)
	#Obtenemos el centro de la trayectoria que describe la partícula
	Centro_Trayectoria_Particula = centro_Circunferencia(Velocidad_Particula, Posicion_Particula, Radio_Trayectoria);

	#Obtengamos la velocidad de la partícula justo cuando está en la posición "Punto_Interseccion"
	Velocidad_Punto_Interseccion = velocidad_Particula(Punto_Interseccion, Velocidad_Particula, Centro_Trayectoria_Particula);

	return Velocidad_Punto_Interseccion
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

			#Corroboremos si el punto cae dentro del segmento adecuado
        	if abs(norm(P1 - Punto_Colision) + norm(P2 - Punto_Colision) - norm(P1 - P2)) > 1e-6
        		Tiempo_Vuelo = Inf;
        	end

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

		#Corroboremos si el punto cae dentro del segmento adecuado
        if abs(norm(P1 - Punto_Colision) + norm(P2 - Punto_Colision) - norm(P1 - P2)) > 1e-6
        	Tiempo_Vuelo = Inf;
        end

        push!(Arreglo_Tiempos_Vuelo, Tiempo_Vuelo); #Añadimos el tiempo de vuelo al arreglo

        if Tiempo_Vuelo < Inf
            Diccionario_Tiempos_Vuelo_Extremos[Tiempo_Vuelo] = [P1, P2]; #Relacionamos el tiempo de vuelo con la puerta de entrada a la partícula
            Diccionario_Tiempos_Vuelo_Punto_Colision[Tiempo_Vuelo] = Punto_Colision; #Relacionamos el tiempo de vuelo con el punto donde la partícula cambia de celda
        end
    end

    #Obtenemos el mínimo de los tiempos de vuelo
    Tiempo_Vuelo_Minimo = minimum(Arreglo_Tiempos_Vuelo);

    #Dependiendo de la trayectoria de la partícula puede ocurrir que no haya intersección con ninguna pared de la celda (cambio de celda Voronoi)
    if Tiempo_Vuelo_Minimo == Inf && Extremos_Entrada[1] == [Inf, Inf] #Tras una colisión el círculo ya no va a colisionar con ninguna pared.
        Punto_Cambio_Celda = [Inf, Inf];
        Velocidad_Cambio_Celda = [Inf, Inf];
        Tiempo_Vuelo_Minimo = Inf;
        Extremos_Entrada_Cambio = [[Inf, Inf], [Inf, Inf]];
        return Punto_Cambio_Celda, Velocidad_Cambio_Celda, Tiempo_Vuelo_Minimo, Extremos_Entrada_Cambio
    elseif Tiempo_Vuelo_Minimo == Inf
        #Intentamos probar con la puerta de entrada
        P1 = Extremos_Entrada[1];
        P2 = Extremos_Entrada[2];
        Punto_Colision, Tiempo_Vuelo = funcion_Interseccion_Particula_Celda(Posicion_Inicial, Velocidad_Inicial, P1, P2);
        Tiempo_Vuelo = Float64(Tiempo_Vuelo);
        push!(Arreglo_Tiempos_Vuelo, Tiempo_Vuelo); #Añadimos el tiempo de vuelo al arreglo

        Diccionario_Tiempos_Vuelo_Extremos[Tiempo_Vuelo] = Extremos_Entrada; #Relacionamos el tiempo de vuelo con la puerta de entrada a la partícula
        Diccionario_Tiempos_Vuelo_Punto_Colision[Tiempo_Vuelo] = Punto_Colision; #Relacionamos el tiempo de vuelo con el punto donde la partícula cambia de celda

        Tiempo_Vuelo_Minimo = Tiempo_Vuelo;
    end

    #EL PRESENTE WHILE FUE LA SOLUCIÓN A UN ERROR QUE DURÓ MESES EN EL PROYECTO... PRESS "F" TO PAY RESPECT.
    while true
        #println(Tiempo_Vuelo_Minimo) ##V-ONLY##
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
            Extremos_Entrada_Cambio = Diccionario_Tiempos_Vuelo_Extremos[Tiempo_Vuelo_Minimo];
            return Punto_Cambio_Celda, Velocidad_Cambio_Celda, Tiempo_Vuelo_Minimo, Extremos_Entrada_Cambio
        else
            filter!(x->x != Tiempo_Vuelo_Minimo, Arreglo_Tiempos_Vuelo); #Eliminamos el tiempo de vuelo mínimo que genera errores de la lista de tiempos de vuelo
            Tiempo_Vuelo_Minimo = minimum(Arreglo_Tiempos_Vuelo); #Volvemos a obtener el mínimo del tiempo de vuelo
        end
    end
end

#Función que calcula la posición y la velocidad final de una partícula que se mueve en una trayectoria circular en el sentido de las manecillas del reloj tras un
#tiempo de vuelo dado.
#Posicion_Inicial: Arreglo con las coordenadas [X,Y] de la posición inicial de la partícula.
#Velocidad_Inicial: Arreglo con las componentes [Vx,Vy] de la velocidad inicial de la partícula.
#Tiempo_Vuelo: Número real positivo con el tiempo de vuelo de la partícula.
function estado_Tras_Tiempo_Circulares_Parcial(Posicion_Inicial, Velocidad_Inicial, Tiempo_Vuelo)
	#El tiempo de vuelo está dado en radianes y suponemos que su valor cae en [0, 2*pi]

	#Obtenemos el centro de la trayectoria que describe la partícula
	Centro_Trayectoria_Particula = centro_Circunferencia(Velocidad_Inicial, Posicion_Inicial, Radio_Trayectoria);

	#Obtenemos la longitud de arco que recorrerá la partícula en el tiempo de vuelo asignado
	Longitud_Arco = Tiempo_Vuelo*norm(Velocidad_Inicial);

	#Obtenemos los ángulos en radianes que va a recorrer la partícula en el tiempo de vuelo asignado
	Radianes = Longitud_Arco/Radio_Trayectoria;

	#Obtenemos los radianes que hay desde el eje X a la posición inicial de la partícula cuando se recorre la trayectoria en el sentido de las manecillas del reloj
	Angulo_Posicion_Inicial = angulo_Posicion_Circulo(Posicion_Inicial, Centro_Trayectoria_Particula);

	#Obtenemos los radianes de la partícula tras el tiempo de vuelo, luego hacemos modulo 2*pi para que esté determinado de manera única cada punto en el círculo
	Angulo_Posicion_Final = mod(Angulo_Posicion_Inicial + Radianes, 2*pi);

	#Generamos la posición final de la partícula empleando el ángulo final
	Posicion_Final = Centro_Trayectoria_Particula + Radio_Trayectoria*[cos(Angulo_Posicion_Final), -sin(Angulo_Posicion_Final)];

	#Obtenemos la velocidad final de la partícula
	Velocidad_Final = velocidad_Particula(Posicion_Final, Velocidad_Inicial, Centro_Trayectoria_Particula);

	return Posicion_Final, Velocidad_Final
end

function cambio_Celda_Esferas_Duras_Trayectoria_Circular(Posicion_Inicial, Velocidad_Inicial, Arreglo_Vertices_Celda, Coordenadas_Segmento_Entrada, Posicion_Obstaculo, Diccionario_Segmentos_Centro_Vecino)
	Tiempo_Acumulado_Vuelo = 0; #Variable que contendrá el tiempo de vuelo que llevamos acumulado

	while true
		Posicion_Colision, Velocidad_Colision, Tiempo_Vuelo_Colision = colision_Obstaculo(Posicion_Inicial, Velocidad_Inicial, Radio_Obstaculo, Posicion_Obstaculo, colision_Esferas_Duras_Campo_Magnetico, velocidad_Esferas_Duras_Campo_Magnetico);

		Posicion_Cambio_Celda, Velocidad_Cambio_Celda, Tiempo_Vuelo_Cambio_Celda, Extremos_Cambio_Celda = colision_Segmento(Posicion_Inicial, Velocidad_Inicial, Arreglo_Vertices_Celda, Coordenadas_Segmento_Entrada, Posicion_Obstaculo, Diccionario_Segmentos_Centro_Vecino, interseccion_Particula_Celda_Trayectoria_Circular, velocidad_Cambio_Celda_Trayectoria_Circular);

		if Tiempo_Vuelo_Colision == Tiempo_Vuelo_Cambio_Celda == Inf
			return Posicion_Cambio_Celda, Velocidad_Cambio_Celda, Tiempo_Vuelo_Cambio_Celda, Extremos_Cambio_Celda
		end

		if Tiempo_Vuelo_Colision < Tiempo_Vuelo_Cambio_Celda #La partícula colisiona antes de moverse de polígono
            #Actualizamos las variables
            Tiempo_Acumulado_Vuelo += Tiempo_Vuelo_Colision; #Se reduce al tiempo de vuelo restante el tiempo en colisionar
            Posicion_Inicial = Posicion_Colision; #La nueva posición de la partícula es la posición de la colisión
            Velocidad_Inicial = Velocidad_Colision; #La nueva velocidad de la partícula es la velocidad tras la colisión
            push!(Central_Obstacles, Posicion_Obstaculo); ##V-ONLY##
        	push!(Positions_Particle, Posicion_Inicial); ##V-ONLY##
        	push!(Velocities_Particle, Velocidad_Inicial); ##V-ONLY##
        	push!(Vertices_Door, [[Inf, Inf], [Inf, Inf]]); ##V-ONLY##
        	push!(Vertices_Voronoi, Arreglo_Vertices_Celda); ##V-ONLY##
                
            #Si hay colisión permitimos que la partícula salga por donde entró
            Coordenadas_Segmento_Entrada = [[Inf, Inf], [Inf, Inf]];
        else #La partícula se mueve de polígono antes de colisionar con el obstáculo de su polígono contenedor
        	Tiempo_Acumulado_Vuelo += Tiempo_Vuelo_Cambio_Celda; #Actualizamos el tiempo de vuelo restante
            return Posicion_Cambio_Celda, Velocidad_Cambio_Celda, Tiempo_Acumulado_Vuelo, Extremos_Cambio_Celda 
        end
    end
end

function estado_Tras_Tiempo_Circular(Posicion_Inicial, Velocidad_Inicial, Arreglo_Vertices_Celda, Coordenadas_Segmento_Entrada, Posicion_Obstaculo, Diccionario_Segmentos_Centro_Vecino, Tiempo_Vuelo_Parcial)
	while Tiempo_Vuelo_Parcial > 0.0
		Posicion_Colision, Velocidad_Colision, Tiempo_Vuelo_Colision = colision_Obstaculo(Posicion_Inicial, Velocidad_Inicial, Radio_Obstaculo, Posicion_Obstaculo, colision_Esferas_Duras_Campo_Magnetico, velocidad_Esferas_Duras_Campo_Magnetico);

		Posicion_Cambio_Celda, Velocidad_Cambio_Celda, Tiempo_Vuelo_Cambio_Celda, Extremos_Cambio_Celda = colision_Segmento(Posicion_Inicial, Velocidad_Inicial, Arreglo_Vertices_Celda, Coordenadas_Segmento_Entrada, Posicion_Obstaculo, Diccionario_Segmentos_Centro_Vecino, interseccion_Particula_Celda_Trayectoria_Circular, velocidad_Cambio_Celda_Trayectoria_Circular);

		if Tiempo_Vuelo_Colision < Tiempo_Vuelo_Cambio_Celda #La partícula colisiona antes de moverse de polígono
            if (Tiempo_Vuelo_Parcial - Tiempo_Vuelo_Colision) > 0.0 #Hay tiempo para la colisión completa
                #Actualizamos las variables
                Tiempo_Vuelo_Parcial -= Tiempo_Vuelo_Colision; #Se reduce al tiempo de vuelo restante el tiempo en colisionar
                Posicion_Inicial = Posicion_Colision; #La nueva posición de la partícula es la posición de la colisión
                Velocidad_Inicial = Velocidad_Colision; #La nueva velocidad de la partícula es la velocidad tras la colisión
                Coordenadas_Segmento_Entrada = [[Inf, Inf], [Inf, Inf]]; #Si hay colisión permitimos que la partícula salga por donde entró
            else #Se agota el tiempo de vuelo antes de la colisión completa
                Posicion_Inicial, Velocidad_Inicial = estado_Tras_Tiempo_Circulares_Parcial(Posicion_Inicial, Velocidad_Inicial, Tiempo_Vuelo_Parcial);
                return Posicion_Inicial, Velocidad_Inicial, Coordenadas_Segmento_Entrada
            end
        else #La partícula se mueve de polígono antes de colisionar con el obstáculo de su polígono contenedor
            if (Tiempo_Vuelo_Parcial - Tiempo_Vuelo_Cambio_Celda) > 0.0 #Hay tiempo para que la partícula salga del polígono
            	throw(ErrorException())
            else #No hay tiempo de vuelo restante para salir del polígono
                Posicion_Inicial, Velocidad_Inicial = estado_Tras_Tiempo_Circulares_Parcial(Posicion_Inicial, Velocidad_Inicial, Tiempo_Vuelo_Parcial);
                return Posicion_Inicial, Velocidad_Inicial, Coordenadas_Segmento_Entrada
            end
        end
	end
end