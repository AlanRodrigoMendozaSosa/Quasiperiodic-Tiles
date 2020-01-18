#This is a function made by Atahualpa that gives the position of the collition and the time of
#collition of a particle and a sphere
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

" velo_col(x1,x2,v), dado el punto de colisión x1, el centro de la esfera x2 y la velocidad de la particula
regresa la velocidad v después de la colisión si consideramos que se refleja (sigue la ley de Snell de la reflexión)"
function velocidad_Esfera_Dura(x1,x2,v)
    n=(x1-x2)
    n=n/norm(n)
    vn=dot(n,v)*n
    v=v-2*vn
    v=v/norm(v)
    return v
end