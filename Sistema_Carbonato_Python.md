
# Resolviendo el sistema carbonato usando Python
****
Dr. Fernando J. Gomez (CICTERRA-CONICET-UNC)

## Introducción

Las ecuaciones que expresan el equilibrio del sistema carbonato, tal como se lo detalla a continuación, puede manipularse de forma de obtener un polinomio que, al resolverlo, nos permite calcular el pH y por consiguiente las concentraciones relativas de otros compuestos como ser $CO_2(aq)$, $HCO_3^-$ y $CO_3^{-2}$. Esto, junto a la concentración de $Ca^{+2}$ nos permite calcular el estado de saturación ($\Omega$) para diferentes minerales de carbonato de calcio como ser calcita ($\Omega_c$) y aragonita ($\Omega_a$).
En este ejemplo utilizaremos la función de **numpy** llamada *numpy.roots* (que en el código abreviamos como *np.roots*) para encontrar las raíces del polinomio de quinto orden que resulta de resolver el equilibrio del sistema carbonato en función del carbono inorngánico disuelto ($DIC$) y la alcalinidad ($ALK$). Crearemos una función que tiene como entrada o "input" la temperatura ($TK$) (en grados Kelvin), la salinidad ($S$), la profundidad ($z$), la alcalinidad ($ALK$) y el carbono inorgánico disuelto ($DIC$) y que nos permite calcular ("output") $fCO_{2}$, $pH$, $CO_2$, $HCO3^{-}$, $CO3^{-2}$, $H^{+}$, $OH^{-}$, $B(OH)_3$, $B(OH)_4^-$, y el *factor de Revelle* ($R$).

El ejemplo esta tomado de **Numerical Analysis with R: Solutions to ODEs and PDEs** de Graham W. Griffiths. En este libro el autor crea esta función en **R**, aqui simplemente trasladaremos dicha función (con ligeras modificaciones) en el lenguaje de **Python**.
La idea es mostrar como crear y utilizar funciones en **Python** y librerías de **numpy** para obtener las raíces (reales) del polinomio que nos permite obtener la concentración de $H^-$ y por consiguiente el $pH$ y los demás componentes del sistema carbonato ($CO_2$, $HCO3^{-}$, $CO3^{-2}$).

Las ecuaciones que dan bases a nuestro ejemplo son las siguientes (ver mas detalles en Griffiths, 2016):

Segun la Ley de Henry podemos calcular las concentraciones en equilibrio entre la presión parcial de $CO_2$ y la concentración de dióxido de carbono disuelto en el agua $[CO_2(aq)]$, donde $K_H$ es la constante de Henry (mol/kg/atm):

$$pCO_2=\frac{[CO_2(aq)]}{K_H}$$

Cuando el $CO_2(aq)$ se disuelve en el agua se generan las siguientes reacciones de equilibrio:

$$CO_2(aq) + H_2O \rightleftharpoons H_2CO_3$$

$$H_2CO_3 \rightleftharpoons HCO_3^- + H^+$$

$$HCO_3^- \rightleftharpoons CO_3^{-2} + H^+$$

En condiciones de equilibrio las concentraciones relativas estan controladas por las correspondientes constantes de equilibrio:

$$K_{CO_2(aq)}=\frac{[CO_2(aq)]}{H_2CO_3}$$

$$K_{H_2CO_3}=\frac{[H^+][HCO_3^-]}{H_2CO_3}$$

Usando $H_2CO_3$ para igualar las ecuaciones anteriores obtenemos la primera y segunda constante de disociación para el acido carbónico:

$$K_{a_1}=\frac{[H^+][HCO_3^-]}{CO_2(aq)}$$

$$K_{a_2}=\frac{[H^+][CO_3^{-2}]}{HCO_3^-}$$

Tambien es necesario considerar la disociación del agua y su consiguiente constante de disociación:

$$H^+ + OH^- \rightleftharpoons H_2O$$

$$K_W=[OH^-][H^+]$$


Otra variable necesaria es el carbono inorgánico disuelto ($DIC$):

$$DIC=[HCO_3^-]+[CO_3^{-2}]+[CO_2(aq)]$$

Al combinar las ecuaciones anteriores podemos obtener $[CO_2(aq)]$, $[HCO_3^-]$ y $[CO_3^{-2}]$:

$$[CO_2(aq)]= \frac{DIC}{1+\frac{K_{a1}}{[H^+]}+\frac{K_{a1}K_{a2}}{[H^+]^2}}$$


$$[HCO_3]= \frac{DIC}{\frac{[H^+]}{K_{a1}}+1+\frac{K_{a2}}{[H^+]}}$$

$$[CO_3^{-2}]= \frac{DIC}{1+\frac{[H^+]^2}{K_{a1}K_{a2}}+\frac{[H^+]}{K_{a2}}}$$

El ácido bórico tambien debe ser considerado, el mismo se disocia mediante la siguiente reacción:

$$B(OH)_3+H_2O \rightleftharpoons B(OH)_4^-+H^+$$

 y su equilibrio se encuentra controlado por la siguiente expresión donde K_B es la constante de equilibrio del ácido bórico:
 
 $$K_B=\frac{[B(OH)_4^-][H^+]}{[B(OH)_3]}$$
 
 La cantidad todal de boro $B_T$ es conservativa y depende de la salinidad:
 
 $$B_T=0.000416\frac{s}{35.0}$$
 
 Esto es igual a:
 
 $$B_T=[B(OH)_3]+[B(OH)_4^-]$$
 
 Ahora si estamos en condiciones de definir la alcalinidad total ($TA$) (esto despreciando otros compuestos menores que tambien contribuyen a la $TA$):
 
 $$TA=[HCO_3^-]+2[CO_3^{-2}]+[B(OH)_4^-]+[OH^-]-[H^+]$$


Ahora estamos en condiciones de obtener el polinomio de quinto orden (ver Bacastow y Keeling, 1973) que nos permite calcular la concentración de $H^+$ (y por consiguiente el $pH$):

$$0=\sum_{i=0}^5 a_i[H^+]^i$$

cuyos coeficientes son los siguientes:

Para $a_6$:

$$a_6=K_WK_{a1}K_{a2}K_{B}$$

Para $a_5$:

$$a_5=-2K_{a1} K_{a2} DIC  K_{B}-K_B B_T K_{a1} K_{a2} - K_W K_{a1} K_B+TA K_{a1} K_{a2}K_B-TAK_{a1}K_{a2}K_B-K_WK_{a1}K_{a2}$$

Para $a_4$:

$$a_4= -K_{a1}DICK_{B}-2K_{a1}K_{a2}DIC-K_BB_TK_{a1}-K_WK_{a1}+TAK_{a1}K_{a2}+K_{a1}K_{a2}K_{B}+TAK_{a1}K_{B}-K_WK_B$$

Para $a_3$:

$$a_3= -K_{a1}DICK_{B}-2K_{a1}K_{a2}DIC-K_BB_TK_{a1}-K_WK_{a1}+TAK_{a1}K_{a2}+K_{a1}K_{a2}K_{B}+TAK_{a1}K_{B}-K_WK_B$$

Para $a_2$:

$$a_2=K_{a1}+TA+K_B$$

Para $a_1$:

$$a_1=1$$


Para calcular el grado de saturación de la calcita o aragonita la reacción y la constante de equilibrio son (respectivamente):

$$CaCO_3 \rightleftharpoons Ca^{+2}+CO_3^{-2}$$

$$K'_{sp}=[Ca^{+2}][CO_3^{-2}]$$

El grado de saturacion esta representado por $\Omega$:

$$\Omega=\frac{[Ca^{+2}][CO_3^{-2}]}{K'_{sp}}$$

Si $\Omega$ es mayor a $1$ el sistema esta sobresaturado y si es menor a $1$ el sistema esta subsaturado en esos minerales. $K'_{sp}(cal)$ y $K'_{sp}(ar)$ son las constantes de equilibrio para la calcita y aragonita respectivamente.


## Implementando las ecuaciones con Python 
El primer paso sera crear la función, luego comprobar que nos da los resultados esperados y luego calcular los parámetros del sistema carbonato para un *set* de valores de $DIC$.


```python
#Este codigo esta basado en el codigo titulado CarbEq.R el libro Numerical Analysis with R.
#En el trabajo original de Emerson and Hedges (2008) p129/130 el codigo estaba en Matlab,
#luego Griffiths lo escribio en R y aqui lo trasladamos a Python.
# La funcion calcula fCO2, HCO3, and CO3 a partir de TK, s, z, ALK, DIC
#Aqui agregamos Ca+2 como input para luego poder calcular el grado de saturacion
#en carbonatos (Calcita y Aragonita).
# Ejemplo:
# Input:    carbEq (293.15, 35, 0, 2.300, 1.970)
# Output:   fco2 <- 0.00028007
#           pH <- 8.1715
#           co2 <- 9.0764e-006
#           hco3 <- 0.0017
#           co3 <- 0.00020373
#------------------------------------------

#Llamamos las librerías a utilizar
import numpy as np
import matplotlib.pyplot as plt
import math
#%matplotlib inline #en caso queremos graficar

#Aqui creamos la función que hará los cálculos
def carbEq(TK, s, z, alk, dic, Ca):
    alk=alk/1000 # alk now in (mol/m3)
    dic=dic/1000 # dic now in (mol/m3)
    #------------------------------------------
    R=83.131
    
    # Calculamos boratos totales (Tbor) a partir de clorinidad
    Tbor=.000416 * s / 35.0
    
    # Calculamos KH, el coeficiente de Henry (cf. a Weiss, 1974)
    U1= -60.2409+93.4517 * (100/TK)+23.3585 * math.log(TK/100)
    U2= s * (.023517 -0.023656 * (TK/100)+.0047036 * (TK/100)**2)
    KH=np.exp(U1+U2)
    
    # Calculamos KB de la TK y S (cf. a Dickson, 1990)
    KB=np.exp((-8966.9-2890.53 * s**0.5-77.942 * s+1.728 * s**1.5 -0.0996*s**2)/TK+\
              148.0248+137.1942 * s**0.5+1.62142 * s -(24.4344+25.085 * s**0.5+0.2474 * s)\
              * np.log(TK)+0.053105 * s**0.5 * TK)
    
    # Calculamos K1 y K2 (cf. a Luecker et al., 2000)
    K1=10**(-(3633.86/TK-61.2172+9.67770 * math.log(TK)-0.011555*s+0.0001152 * s**2))
    K2=10**(-(471.78/TK+25.9290-3.16967 * math.log(TK)-0.01781*s+0.0001122 * s**2))
    
    #-----------------------------------------------
    # Dependencia de Kw en función de la temperatura, TK (cf. a DoE, 1994)
    KW1=148.96502-13847.26/TK-23.65218*math.log(TK)
    KW2=(118.67/TK-5.977+1.0495*math.log(TK))*s**0.5-0.01615*s
    KW=np.exp(KW1+KW2)
    
    #-----------------------------------------------
    # Resolvemos para obtener H+ (cf. a Zeebe and Wolf-Gladrow, 2000)
    a1=1
    a2=(alk+KB+K1)
    a3=(alk*KB-KB*Tbor-KW+alk*K1+K1*KB+K1*K2-dic*K1)
    a4=(-KW*KB+alk*KB*K1-KB*Tbor*K1-KW*K1+alk*K1*K2+KB*K1*K2-dic*KB*K1-2*dic*K1*K2)
    a5=(-KW*KB*K1+alk*KB*K1*K2-KW*K1*K2-KB*Tbor*K1*K2-2*dic*KB*K1*K2)
    a6= -KB*KW*K1*K2
    p= [a1, a2, a3, a4, a5, a6]
    r = np.roots(p)
    h= max(np.real(r)) # Esto es para seleccionar la raiz real mas grande
    #
    # Calculamos HCO3, CO3 and CO2aq, usando DIC, AlK y H+
    hco3=dic/(1+h/K1+K2/h)
    co3=dic/(1+h/K2+h*h/(K1*K2))
    co2=dic/(1+K1/h+K1*K2/(h*h))
    fco2=co2 / KH
    pH=-math.log10(h)
    
    # calculamos B(OH)4- y OH-
    BOH4=KB*Tbor/(h+KB)
    oh=KW/h
    BOH3=Tbor-BOH4

    # Calculamos el factor de Revelle, ver detalles en
    # Egleston, et al, (2010). Revelle revisited: Buffer factors that
    # quantify the response of ocean chemistry to changes in DIC and
    # alkalinity, 'Global Biogeochemical Cycles, vol 24, 1-9.

    ATc= (2 * co3 + hco3) # alkalinidad debido a las especies de C.
    Gdic=(hco3 + 4*co3 + h*BOH4/(KB+h) + h + oh)
    R=dic/(dic - ATc**2/Gdic)
    
    #Constantes de equilibrio Calcita y Aragonita
    Ksp_ar=6.76e-7
    Ksp_cal=4.37e-7
    
    #Saturacion de Calcita y Aragonita
    omega_ar=Ca*co3/Ksp_ar
    omega_cal=Ca*co3/Ksp_cal
    
    return fco2, pH, co2, hco3, co3, h, oh, BOH3, BOH4, R, omega_ar, omega_cal

#-----------------------------------------------------------------------

```

## Probando con un par de ejemplos

Una vez creada la función primero veremos si nos da valores
similares a los de Griffiths, para *input* carbEq ($TK$=293.15, $s$=35, $z$=0, $ALK$=2.300, $DIC$=1.970)
*Output* esperado: $fCO_2$ = 0.00028007, $pH$ = 8.1715, $CO_2$ = 9.0764e-006, $HCO_3^-$ = 0.0017, $CO_3^{-2}$ = 0.00020373

Pueden compararse los 5 primeros valores detallados a continuación con el output esperado.       


```python
print(carbEq (293.15, 35, 0, 2.300, 1.970, 0.0103))
```

    (0.00028007007970335009, 8.171545623563299, 9.0763555320580811e-06, 0.0017301985592125143, 0.00023072508525542773, 6.7368112011194785e-09, 5.6957740282050387e-06, 0.00031333776694037336, 0.00010266223305962662, 9.1849481830157149, 3.5154857664658365, 5.4381427417183188)


Luego la usaremos la función creada para valores de DIC entre 1000 y 2500 y generaremos una tabla con los valores. En este caso tambien incluiremos el calculo de saturacion de la calcita y aragonita expresado como $\Omega$:


```python
%matplotlib inline
dic_values=[1.000, 1.100, 1.200, 1.300, 1.400, 1.500, 1.600,\
            1.700, 1.800, 1.900, 2.000, 2.100, 2.200, 2.300,\
            2.400, 2.500]

for dic in dic_values:
    fco2, pH, co2, hco3, co3, h, oh, BOH3, BOH4, R, omega_ar, omega_cal= carbEq(TK=293.15,\
                                        s=35, z=0, alk=2.300, dic=dic, Ca=0.0103)
    print ('|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    print ('|pH= %-3g CO2= %-3g HCO3-= %-3f CO3-2= %-3f' % (pH, co2, hco3, co3))
    print ('|omega_ar= %-3g omega_cal= %-3g' % (omega_ar, omega_cal))
    print ('|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
```

    |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    |pH= 9.59079 CO2= 4.4382e-08 HCO3-= 0.000222 CO3-2= 0.000778
    |omega_ar= 11.8513 omega_cal= 18.3329
    |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    |pH= 9.38775 CO2= 1.09814e-07 HCO3-= 0.000344 CO3-2= 0.000756
    |omega_ar= 11.5114 omega_cal= 17.8072
    |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    |pH= 9.21142 CO2= 2.33237e-07 HCO3-= 0.000487 CO3-2= 0.000712
    |omega_ar= 10.8547 omega_cal= 16.7912
    |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    |pH= 9.05831 CO2= 4.36406e-07 HCO3-= 0.000641 CO3-2= 0.000659
    |omega_ar= 10.0346 omega_cal= 15.5227
    |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    |pH= 8.92103 CO2= 7.47205e-07 HCO3-= 0.000800 CO3-2= 0.000599
    |omega_ar= 9.13016 omega_cal= 14.1235
    |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    |pH= 8.79321 CO2= 1.20591e-06 HCO3-= 0.000962 CO3-2= 0.000537
    |omega_ar= 8.17935 omega_cal= 12.6527
    |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    |pH= 8.66978 CO2= 1.87459e-06 HCO3-= 0.001125 CO3-2= 0.000473
    |omega_ar= 7.20198 omega_cal= 11.1408
    |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    |pH= 8.54625 CO2= 2.85473e-06 HCO3-= 0.001290 CO3-2= 0.000408
    |omega_ar= 6.2095 omega_cal= 9.60554
    |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    |pH= 8.41799 CO2= 4.32377e-06 HCO3-= 0.001454 CO3-2= 0.000342
    |omega_ar= 5.20981 omega_cal= 8.05912
    |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    |pH= 8.27923 CO2= 6.62001e-06 HCO3-= 0.001617 CO3-2= 0.000276
    |omega_ar= 4.21017 omega_cal= 6.51276
    |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    |pH= 8.1216 CO2= 1.04648e-05 HCO3-= 0.001778 CO3-2= 0.000211
    |omega_ar= 3.22049 omega_cal= 4.98181
    |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    |pH= 7.93149 CO2= 1.76331e-05 HCO3-= 0.001934 CO3-2= 0.000148
    |omega_ar= 2.2609 omega_cal= 3.49741
    |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    |pH= 7.68792 CO2= 3.31633e-05 HCO3-= 0.002076 CO3-2= 0.000091
    |omega_ar= 1.38508 omega_cal= 2.14259
    |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    |pH= 7.38631 CO2= 6.9823e-05 HCO3-= 0.002182 CO3-2= 0.000048
    |omega_ar= 0.727106 omega_cal= 1.12477
    |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    |pH= 7.10354 CO2= 0.000137255 HCO3-= 0.002237 CO3-2= 0.000026
    |omega_ar= 0.388687 omega_cal= 0.601264
    |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    |pH= 6.89671 CO2= 0.000223301 HCO3-= 0.002261 CO3-2= 0.000016
    |omega_ar= 0.243947 omega_cal= 0.377364
    |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

