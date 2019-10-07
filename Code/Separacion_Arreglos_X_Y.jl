function separacion_Arreglo_de_Arreglos_2D(Arreglo)
    #Definimos los contenedores para cada una de las coordenadas
    Coordenadas_X = [];
    Coordenadas_Y = [];
    
    #Para cada arreglo dentro del arreglo de arreglos, extraemos su coordenada X y su coordenada Y
    for i in Arreglo
        push!(Coordenadas_X, i[1]);
        push!(Coordenadas_Y, i[2]);
    end
    
    return Coordenadas_X, Coordenadas_Y
end