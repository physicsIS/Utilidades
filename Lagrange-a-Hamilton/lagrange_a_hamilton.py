# Elmer Barquero
def cambio_L_H(Lagrangiano,coordenadas_generalizadas:list):
    import sympy as smp
    ###########################################################
    #         Generar simbolos velocidades generalizadas
    ###########################################################
    velocidades_generalizadas = []
    for i in range(len(coordenadas_generalizadas)):
        aux = smp.diff(coordenadas_generalizadas[i])
        velocidades_generalizadas.append(aux)

    ###############################################################
    #          Generar simbolos para los momentos canonicos
    ###############################################################
    simbolos_momentos_canonicos = []
    for i in range(len(coordenadas_generalizadas)):
        aux = smp.Symbol(f'p_{{{coordenadas_generalizadas[i]}}}')
        simbolos_momentos_canonicos.append(aux)

    ################################################################
    #                 Generar momentos canonicos
    ################################################################
    momentos_canonicos_en_coordenadas_generalizadas = []
    for i in range(len(coordenadas_generalizadas)):
        aux = smp.diff(Lagrangiano, velocidades_generalizadas[i])
        momentos_canonicos_en_coordenadas_generalizadas.append(aux)
    ################################################################

    ############################################################################
    #  Despejar las velocidades generalizadas en terminos de momentos canonicos
    ############################################################################

    # Sistema de ecuaciones a resolver para poner las velocidades generalizas
    # en terminos de los momentos canonicos

    ecuaciones = []
    for i in range(len(coordenadas_generalizadas)):
        aux = smp.Eq(simbolos_momentos_canonicos[i] , momentos_canonicos_en_coordenadas_generalizadas[i])
        ecuaciones.append(aux)

    # Resolver el sistema de ecuaciones para encontrar cada velocidad generalizada
    # en terminos de los momentos canonicos

    incognitas = tuple(velocidades_generalizadas)

    despeje_velocidades_generalizadas_en_momentos_canonico = smp.solve(ecuaciones, incognitas, dict=True)[0]

    claves_del_despeje = velocidades_generalizadas

    ##################################################################################
    # Realizar suma de la multiplicacion pq_punto de la transformada del Hamiltoneano
    ##################################################################################
    HH = 0
    for i in range(len(velocidades_generalizadas)):
        HH += simbolos_momentos_canonicos[i]*despeje_velocidades_generalizadas_en_momentos_canonico[claves_del_despeje[i]]

    ##########################################################################################################
    #   Realizar el cambio de variable en el Lagragiano, cambiando q_punto por el respectivo momento canonico
    ##########################################################################################################
    L_cambio_variables_q_a_p = Lagrangiano.subs([(velocidades_generalizadas[i],despeje_velocidades_generalizadas_en_momentos_canonico[claves_del_despeje[i]]) for i in range(len(coordenadas_generalizadas))])

    ##########################################################################################################
    #                                           EL HAMILTONEANO
    ##########################################################################################################

    H = (HH - L_cambio_variables_q_a_p).expand()

    return H ,simbolos_momentos_canonicos, momentos_canonicos_en_coordenadas_generalizadas
