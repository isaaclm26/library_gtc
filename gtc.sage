#
#
#	Author: Isaac Lacort Magan
#
#	
#
#
#
#

########## GTC Practica 0 ##########

# 1. Distancia entre dos puntos

def dist(A,B):
    return sqrt((A[0]-B[0])**2+(A[1]-B[1])**2)

# 2. Cuadrado de la distancia entre dos puntos

def dist2(A,B):
    return (A[0]-B[0])**2+(A[1]-B[1])**2

# 3. Area signada

def sarea(A,B,C):
    if A == B or A == C or B == C: return 0
    area = (1/2)*((B[0] - A[0])*(C[1] - B[1]) - (C[0] - B[0])*(B[1] - A[1]))
    return area

# 4. Orientacion de un triangulo

def orientation(A,B,C):  # la salida sera 1, -1 o 0
    return sign(sarea(A,B,C))

# 5. Punto medio de dos puntos

def midPoint(A,B):
    return [(A[0]+B[0])/2,(A[1]+B[1])/2]

# 6. Punto en segmento

def inSegment(P,s):
    return orientation(P,s[0],s[1]) == 0 and min(s[0][0],s[1][0]) <= P[0] <= max(s[0][0],s[1][0]) and min(s[0][1],s[1][1]) <= P[1] <= max(s[0][1],s[1][1])

# 7. Punto en triangulo

def inTriangle(P,t):
    ori = orientation(t[0],t[1],t[2]); 
    return ori*orientation(P,t[1],t[2]) > -1 and ori*orientation(t[0],P,t[2]) > -1 and ori*orientation(t[0],t[1],P) > -1

# 8. Interseccion de dos segmentos

def segmentIntersectionTest(a,b):
    if orientation(a[0],a[1],b[0]) == orientation(a[0],a[1],b[1]) == 0: return min(a) <= max(b) and max(a) >= min(b)
    return orientation(a[0],a[1],b[0])*orientation(a[0],a[1],b[1]) <= 0 and orientation(b[0],b[1],a[0])*orientation(b[0],b[1],a[1]) <= 0

# 9. Interseccion de dos rectas

def lineIntersection(r,s):
    # Comprobamos que no son paralelas
    aR = r[1][1] - r[0][1];  bR = - (r[1][0] - r[0][0]);  
    aS = s[1][1] - s[0][1];  bS = - (s[1][0] - s[0][0]); 
    
    # Primero compruebo si son paralelas o coincidentes
    c = (aR * bS) - (bR * aS);
    if c == 0:
        if sarea(r[0],r[1],s[0]) == 0:
            # Las rectas son coincidentes, devuelvo r, ya que su interseccion es toda la recta
            return r
        else: 
            return []
        
    cR = (r[0][1]*r[1][0]) - (r[1][1]*r[0][0])
    cS = (s[0][1]*s[1][0]) - (s[1][1]*s[0][0])
    return [((bR * cS) - (cR * bS))/c, ((cR * aS) - (aR * cS))/c]

# 10. Mediatriz de dos puntos

def perpendicularBisector(A,B):
    p = midPoint(A,B)
    return [p,[p[0] + A[1] - B[1], p[1] + B[0] - A[0]]]

# 11. Circuncentro de tres puntos

def circumcenter(a,b,c):
    
    if sarea(a,b,c) == 0: 
        if a == b == c: 
            return a
        elif a == b != c or a == c != b: 
            return midPoint(c,b)
        elif a != b == c: 
            return midPoint(a,b)
        else: 
            return []
            
    return lineIntersection(perpendicularBisector(a,b),perpendicularBisector(b,c))

# 12. Inclusion de un punto en circulo

# Fuente: Manuel abellanas oar
def svolume(a,b,c,d):
    M = [a,b,c,d]
    for i in M:
        i.insert(0,1)
    m = matrix(M).transpose()
    #m = traspose(M)
    return m.det()/6
#Fuente: Manuel Abellanas Oar
def incircleVol(a,b,c,d):
    ar = sarea(a,b,c)
    if ar == 0:
        return -1
    
    A,B,C,D = copy(a), copy(b), copy(c), copy(d)
    A.append(a[0]**2+a[1]**2)
    B.append(b[0]**2+b[1]**2)
    C.append(c[0]**2+c[1]**2)
    D.append(d[0]**2+d[1]**2)
    v = svolume(A,B,C,D)
    return -sign(ar*v)

# Este metodo calcula el circumcentro del circulo
# circunscrito llamando al metodo circumcenter
# y evalua las distancias al punto d, para saber si
# se encuentra en el circulo

def incircleInt(a,b,c,d):
    if sarea(a,b,c) == 0: return -1

    p = circumcenter(a,b,c)
    # No necesito hacer la raiz cuadrada
    l = dist2(p,d); r = dist2(p,a)
    
    if    l >  r: return -1
    elif  l == r: return  0
    else        : return  1

# Este metodo hace uso del arco capaz para calcular el
# circumcentro. 

def incircle(a,b,c,d):
    sa = sarea(a,b,c)
    if sa == 0: return -1
    # Trabajo con los vertices A, B y C que deben tener
    # orientacion positiva
    A = a
    if sa > 0:
        B = b
        C = c
    else:
        B = c
        C = b
    # Creo los vectores del triangulo
    AB = [B[0]-A[0],B[1]-A[1]]
    AC = [C[0]-A[0],C[1]-A[1]]
    BC = [C[0]-B[0],C[1]-B[1]]
    # Calculo el coseno y seno del angulo que forman los vectores AB y AC
    denom = dist2(A,B)*dist2(A,C)
    num = (AB[0]*AC[0] + AB[1]*AC[1])
    cos_a = num / sqrt(denom)
    sin_a = sqrt(1 - (num*num / denom))
    # Este vector es el resultado de aplicar la aplicacion
    # de giro en sentido negativo con el seno y el coseno
    # que acabo de calcular.
    v = [cos_a*BC[0] + sin_a*BC[1], - sin_a*BC[0] + cos_a*BC[1]]
    # Giro el vector anterios 90 en sentido positivo
    u = [-v[1],v[0]]
    # calculo la interseccion entre la mediatriz de BC
    # y la recta que pasa por B con direccion u
    p = lineIntersection(perpendicularBisector(B,C),[B,[B[0] + u[0], B[1] + u[1]]])
    # Utilizo las distancias al cuadrado para ahorrar
    # coste computacional
    l = dist2(p,d); r = dist2(p,a)
    
    if    l >  r: return -1
    elif  l == r: return  0
    else        : return  1

    
########## GTC Practica 1 ##########

# 1. Hallar el maximo/minimo en abscisa/ordenada

def max_min(p):
    if len(p) == 0: return []
    xmin = xmax = ymin = ymax = p[0]
    for i in p:
        if i[0] < xmin[0]: xmin = i
        if i[1] < ymin[1]: ymin = i
        if i[0] > xmax[0]: xmax = i
        if i[1] > ymax[1]: ymax = i
    return xmin,xmax,ymin,ymax

# 2. Hallar el maximo/minimo segun la direccion dada por un vector

def max_min_v(p,v):
    if len(p) == 0: return []
    pmin = pmax = p[0]
    pvmin = pvmax = p[0][0]*v[0] + p[0][1]*v[1]
    for i in p:
        pv = i[0]*v[0] + i[1]*v[1]
        if pv < pvmin: pvmin = pv; pmin = i;
        elif pv > pvmax: pvmax = pv; pmax = i
    return pmin,pmax

# 3. Si los puntos estan en el primer cuadrante, hallar aquel cuyo vector de posicion forma un angulo maximo / minimo con el eje de abscisas.

def max_min_a(p):
    if len(p) == 0: return []
    pmin = pmax = [-1,-1]
    amin = -2; amax = 2 #A mayor angulo menor coseno
    for i in p:
        if min(i) >= 0:
            an = i[0]/sqrt(i[0]**2 + i[1]**2)
            if an < amax: amax = an; pmax = i
            if an > amin: amin = an; pmin = i
    return pmin,pmax

# 4. Dividir la lista en dos listas del mismo tamanyo separadas por una recta vertical / horizontal / de pendiente dada.

def list2(p,m):
    l = len(p)
    if l == 0 or l == 1: return p 
    max = (l/2).ceil()
    list = sorted(p, key = lambda x: x[0]*m - x[1]);
    return [list[0:max],list[max::]]

# 5. Hallar el menor rectangulo isotetico (de lados paralelos a los ejes) que los contenga. ("bounding box").

def minRectangle(p):
    min_x = max_x = p[0][0]; min_y = max_y = p[0][1]
    for i in p:
        if i[0] < min_x: min_x = i[0]
        elif i[0] > max_x: max_x = i[0]
        if i[1] < min_y: min_y = i[1]
        elif i[1] > max_y: max_y = i[1]
    return [[min_x,min_y],[min_x,max_y],[max_x,max_y],[max_x,min_y]]

# 6. Hallar el menor rectangulo de lados paralelos o perpendiculares a un vector v que los contenga. ("bounding box orientado").

def minRectangleVector(p,v):
    max_min = max_min_v(p,v); min_v = max_min[0]; max_v = max_min[1]
    max_min = max_min_v(p,[v[1],-v[0]]); min_pv = max_min[0]; max_pv = max_min[1]
    point(p,color="blue",size=20)+line([min_v,[min_v[0]+v[1],min_v[1]-v[0]]],color="red")
    return [lineIntersection([min_v,[min_v[0]+v[1],min_v[1]-v[0]]],[min_pv,[min_pv[0] + v[0],min_pv[1] + v[1]]]),
            lineIntersection([min_v,[min_v[0]+v[1],min_v[1]-v[0]]],[max_pv,[max_pv[0] + v[0],max_pv[1] + v[1]]]),
            lineIntersection([max_v,[max_v[0]+v[1],max_v[1]-v[0]]],[max_pv,[max_pv[0] + v[0],max_pv[1] + v[1]]]),
            lineIntersection([max_v,[max_v[0]+v[1],max_v[1]-v[0]]],[min_pv,[min_pv[0] + v[0],min_pv[1] + v[1]]])]

# 7. Ordenar la lista por abscisas / ordenadas

def sort_abs_ord(p):
    return sorted(p, key = lambda x: [x[1],x[0]])

# 8. Ordenar la lista segun una direccion dada por un vector.

def sortv(p,v):
    return sorted(p, key = lambda x: x[0]*v[0] + x[1]*v[1])

# 9. Ordenar la lista angularmente respecto de un punto dado.

def angularSort(p,A):
    p_t = []; p_d = [];
    def clave(x):
        return (x[0]-A[0])/dist(x,A)
    for i in p:
        if i[1] > A[1]: p_t.append(i)
        elif i[1] < A[1]: p_d.append(i)
    p_t = sorted(p_t, key = clave); p_d = sorted(p_d, key = clave)
    p_t.reverse()
    if A in p:
        if len(p_t) != 0 and len(p_d) != 0 and (clave(p_t[-1]) < 0 or clave(p_d[0]) < 0):
            return [A]+p_t+p_d
        return [A]+p_d+p_t
    return p_d+p_t

########## GTC Practica 2 ##########

# poligonizacion X-monotona

def polygonization(p):
    L = []; R = []
    pmax = max(p); pmin = min(p)
    for i in p:
        if orientation(pmin,pmax,i) >= 0: L.append(i)
        else: R.append(i)
    L = sorted(L)
    R = sorted(R)
    L.reverse()
    return L+R
    
# poligonizacion estrellada

def starPolygonization(p):    
    return angularSort(p,midPoint(p[0],p[1]))

# Halla el poligono interseccion de P con el semiplano izquierdo de r

def clipping(P,r):
    if len(P) == 0: return []
    p = []; newori = 0; lastori = orientation(r[0],r[1],P[len(P)-1])
    if lastori >= 0: p.append(P[len(P)-1])
    for i in range(len(P)):
        newori = orientation(r[0],r[1],P[i])
        if newori != 0 and newori != lastori:
            p.append(lineIntersection(r,[P[i-1],P[i]]))
        if newori >= 0: 
            p.append(P[i])
        lastori = newori
    return p

# Poligono cuyos puntos se unen a todos los vertices mediante segmentos
# completamente contenidos en el poligono

def kernel(p):
    C=copy(p)
    for i in range(len(p)):
        C=clipping(C,[p[i-1],p[i]])

    return C

### funcion auxiliar que detecta si un vertice es oreja en tiempo O(n)    

def earTest(P,i):
    l = len(P)
    if orientation(P[i-1],P[i],P[(i+1)%l]) <= 0: return False
    j = 0
    while j < l:
        if j in [(i-1)%l,i,(i+1)%l]: j += 1
        elif inTriangle(P[j],[P[i-1],P[i],P[(i+1)%l]]): return False
        else: j += 1
    return True
    

# devolver el indice del vertice oreja

def ear(P):
    i = 0;
    while (i < len(P))and(earTest(P,i) == False):
        i += 1
    return i

# Metodo que devuelve la triangulacion calculada a partir de sus orejas.

def otectomyTriangulation(P):
    p = copy(P); t = [];
    while len(p) > 3:
        i = ear(p); n = len(p)
        t.append([p[i-1],p[i],p[(i+1)%n]])
        del p[i]
    t.append(p)
    return t

# area signada de un poligono

def sareaPolygon(p):
    a=0
    for i in range(len(p)):
        a+= sarea(p[0],p[i-1],p[i])
    return a

# Metodo que te dice si un punto es interior a un poligono

def pointInPolygon(Q,c):
    cont = 0;
    # Genero el punto p que me asegura que el segmento [p,c] 
    # tocal en al menos un punto, si c es punto interior.
    xmin, xmax, ymin, ymax = max_min(Q)
    p = [xmin[0],ymin[1]]
    
    for i in range(len(Q)):
        if segmentIntersectionTest([p,c],[Q[i-1],Q[i]]):
            cont += 1
            
    # Si el numero de aristas con las que corta es par o 0 entonces se encuentra fuera del poligono
    if cont % 2 == 0:
        return false
    else:
        return true
    
# Diagonal interna de un poligono

def diagonal(p):
    i = 0; n = len(p)
    while sarea(p[i-1],p[i],p[(i+1)%n]) <= 0:
        i+=1
    inT = []; j = 0
    while j < n:
        if j in [(i-1)%n,i,(i+1)%n]:
            j+=1
        else:
            if inTriangle(p[j],[p[i-1],p[i],p[(i+1)%n]]):
                inT.append(j)
            j+=1
    if len(inT) == 0:
        return [(i-1)%n,(i+1)%n]
    
    sa_max = abs(sarea(p[i-1],p[(i+1)%n],p[inT[0]]));
    w = inT[0];
    aux_sa = 0;
    for k in inT:
        aux_sa = abs(sarea(p[i-1],p[(i+1)%n],p[k]))
        if aux_sa > sa_max:
            w = k;
            sa_max = aux_sa;
    
    return [i,w]

### Triangulacion recursiva mediante el calculo de diagonales

def splitPolygon(p,i,j):
    k = min(i,j)
    l = max(i,j)
    return [p[k:l+1],p[l::]+p[0:k+1]]

def polygonTriangulation(p):
    if len(P) < 3: return p
    if len(P) == 3: return [p]
    T = []
    pendientes = [p]
    while pendientes:
        Q = pendientes.pop()
        d = diagonal(Q)
        S = splitPolygon(Q,d[0],d[1])
        if len(S[0]) <= 3: T.append(S[0])
        else: pendientes.append(S[0])
        if len(S[1]) <= 3: T.append(S[1])
        else: pendientes.append(S[1])
    return T


############# GTC 3 #################

# Calcula el cierre convexo comprobando vertice por vertice si es exterior
def fuerzaBrutaVertices(P):
    Q = copy(P)
    for i in P:
        for j in P:
            for k in P:
                w = 0
                if i != j != k and i != k and orientation(i,j,k) != 0:
                    while w < len(Q):
                        if Q[w] != i and Q[w] != j and Q[w] != k and inTriangle(Q[w],[i,j,k]): del Q[w]
                        else: w+=1
    return angularSort(Q,max(Q))

# Calcula el cierre convexo comprobando cada par de aristas

def fuerzaBrutaAristas(P):
    a = []
    esbuena = True
    for i in range(len(P)):
        for j in range(i+1,len(P)):
            esbuena = True
            lastor = orientation(P[i],P[j],P[i-1])
            for k in range(len(P)):
                if i != k and k != j:
                    newor = orientation(P[i],P[j],P[k])
                    if inSegment(P[k],[P[i],P[j]]) or (lastor != 0 and newor != 0 and newor != lastor):
                        esbuena = False
                        break
                    if newor != 0:
                        lastor = newor
            if esbuena and [P[i],P[j]] not in a:
                a.append([P[i],P[j]])
    A = a[0]
    del a[0]
    for i in range (len(a)):
        for j in range(len(a)):
            if A[-1] == a[j][0]: 
                A.append(a[j][1])
                del a[j]
                break;
            elif A[-1] == a[j][1]:
                A.append(a[j][0])
                del a[j]
                break;
             
    return A

# Calculo del Cierre Convexo:

#     - Primero calcula los puntos con los maximos y minimos en los ejes 'x' e 'y'
#     y los mete en Una Lista De Forma ordenada en el sentido de las agujas del reloj.
#     Para ello utilizo los metodo auxiliares 'AbscisasOrdenadas' e 'InsertarOrdenadamente'.

#     - En segundo lugar meto en 4 listas los puntos que optan a estar en el cierre
#     convexo, es decir, si tengo las listas de puntos con el valor maximo para 'x' e 'y',
#     llamadas 'xmax' y 'ymax', creo la lista 'c1' para los elementos que se encuentren
#     entre los puntos maximos de las listas y cuya orientacion 'orientation(xmax,i,ymax)'
#     sea postiva.

#     - En el ultimo paso anyado a las listas 'arc' (arcos) los puntos con coseno maximo o
#     minimo dependiendo del caso, respecto el eje de abscisas.
 
# Con este engorroso trabajo consigo un metodo para calcular cierres convexos con complejidad
# O(n^2) en el peor de los casos y O(n) en el mejor de estos. Aun siendo asi, salvo en el caso
# de los puntos en la circunferencia, este algoritmo ha superado en tiempos al resto. En el
# caso de la circunferencia se ve bastante desfavorecido quedando a la altura de Jarvis por
# ser O(n^2).

# i es el indice de cada punto que hay que comprobar
# m == 0, de menor a mayor; m == 1, de mayor a menor;
def InsertarOrdenadamente(L,p,i,m):
    j = len(L);
    if m == 0:
        for k in range(j):
            if L[k][i] > p[i]: 
                j = k
                break;
    else: 
        for k in range(j):
            if L[k][i] < p[i]: 
                j = k
                break;
    return L[:j]+[p]+L[j:]

# Esta funcion devuelve los puntos maximos y minimos en abscisas ordenadas. 
# Si son varios devuelve una lista ordenada en sentido positivo de ellos
def AbscisasOrdenadas(P):
    if P == []: return []
    xmax = [P[0]]; xmin = [P[0]]; ymax = [P[0]]; ymin = [P[0]]
    for i in range(1,len(P)):
        if P[i][0] > xmax[0][0]:
            xmax = [P[i]]
        elif P[i][0] == xmax[0][0]:
            xmax = InsertarOrdenadamente(xmax,P[i],1,0)
        
        if P[i][0] < xmin[0][0]:
            xmin = [P[i]]
        elif P[i][0] == xmin[0][0]:
            xmin = InsertarOrdenadamente(xmin,P[i],1,1)
        
        if P[i][1] > ymax[0][1]:
            ymax = [P[i]]
        elif P[i][1] == ymax[0][1]:
            ymax = InsertarOrdenadamente(ymax,P[i],0,1)
        
        if P[i][1] < ymin[0][1]:
            ymin = [P[i]]
        elif P[i][1] == ymin[0][1]:
            ymin = InsertarOrdenadamente(ymin,P[i],0,0)
            
    return xmax,ymax,xmin,ymin

def CierreConvexo(P):
    xmax,ymax,xmin,ymin = AbscisasOrdenadas(P)
    c1 = [ymax[0]]; c2 = [xmin[0]]; c3 = [ymin[0]]; c4 = [xmax[0]]
    for i in P:
        if xmax[-1] != ymax[0] and i[0] > ymax[0][0] and i[1] > xmax[-1][1] and sarea(xmax[-1],i,ymax[0]) > 0:
            c1.append(i)
        elif ymax[-1] != xmin[0] and i[0] < ymax[-1][0] and i[1] > xmin[0][1] and sarea(ymax[-1],i,xmin[0]) > 0:
            c2.append(i)
        elif xmin[-1] != ymin[0] and i[0] < ymin[0][0] and i[1] < xmin[-1][1] and sarea(xmin[-1],i,ymin[0]) > 0:
            c3.append(i)
        elif ymin[-1] != xmax[0] and i[0] > ymin[-1][0] and i[1] < xmax[0][1] and sarea(ymin[-1],i,xmax[0]) > 0:
            c4.append(i)
        
    arc1 = [xmax[-1]]; arc2 = [ymax[-1]]; arc3 = [xmin[-1]]; arc4 = [ymin[-1]]
    arc = arc1
    def keyCuadrante (X):
        return (X[0]-arc[-1][0])/dist(X,arc[-1])
    
        
    while arc1[-1] != ymax[0]:
        arc1.append(max(c1, key = keyCuadrante))
        i = 0
        while i < len(c1):
            if c1[i] == arc1[-1] or sarea(arc1[-1],c1[i],ymax[0]) < 0:
                del c1[i]
            else: i+=1
    
    arc = arc2
    while arc2[-1] != xmin[0]:
        arc2.append(min(c2, key = keyCuadrante))
        i = 0
        while i < len(c2):
            if c2[i] == arc2[-1] or sarea(arc2[-1],c2[i],ymin[0]) < 0:
                del c2[i]
            else: i+=1
    
    arc = arc3
    while arc3[-1] != ymin[0]:
        arc3.append(min(c3, key = keyCuadrante))
        i = 0
        while i < len(c3):
            if c3[i] == arc3[-1] or sarea(arc3[-1],c3[i],ymin[0]) < 0:
                del c3[i]
            else: i+=1
    
    arc = arc4
    while arc4[-1] != xmax[0]:
        arc4.append(max(c4, key = keyCuadrante))
        i = 0
        while i < len(c4):
            if c4[i] == arc4[-1] or sarea(arc4[-1],c4[i],xmax[0]) < 0:
                del c4[i]
            else: i+=1
    
    return xmax[1::]+arc1[1::]+ymax+arc2[1::]+xmin[1::]+arc3[1::]+ymin[1::]+arc4[1::]

# Gift wrapping (Algoritmo de Jarvis)

def Jarvis(P):
    CC = [max(P)]; 
    while True:
        B = P[0]
        if CC[-1] == B:
            B = P[1]
        for i in range(1,len(P)):
            if sarea(CC[-1],B,P[i]) < 0:
                B = P[i]
        if B != CC[0]: 
            CC.append(B)
        else: 
            return CC

# Mediante la poligonizacion estrellada (Algorimo de Graham)

def Graham(P):
    Q = angularSort(P,max(P))
    i = 1;
    while Q[i] != Q[-1]:
        if sarea(Q[i-1],Q[i],Q[(i+1)%len(Q)]) <= 0:
            del Q[i]; i -= 1
        else: i += 1
        
    return Q

# Metodo para preprocesar los datos para el algoritmo de Graham

# Fuente: Manuel Abellanas Oar

def preproccessCC(p):
    def maxv(p,v):
        return max(p, key = lambda x : x[0]*v[0] + x[1]*v[1])
    
    P = [maxv(p,i) for i in [[0,-1],[1,-1],[1,0],[1,1],[0,1],[-1,1],[-1,0],[-1,-1]]]
    
    for i in p:
        for j in range(8):
            if sarea(P[j-1],P[j],i) < 0:
                P.append(i)
    return P


# QuickHull. Metodo de triangulacion recursiva

# Fuente: Manuel Abellanas Oar

#A es el maximo y B el minimo
def QuickHull(P):
    A, B = max(P), min(P)
    def envolvente(P,A,B):
        if len(P) == 0: return [A]
        R = []; L = []
        m = 0 
        for i in range(len(P)):
            if sarea(A,P[i],B) > m:
                m = i
        M = P[m]
        for i in range (len(P)):
            if sarea(M,A,P[i]) > 0:
                R.append(P[i])
            if sarea(B,M,P[i]) > 0:
                L.append(P[i])
        return envolvente(R,A,M)+envolvente(L,M,B)
    return envolvente(P,A,B)+envolvente(P,B,A)

# Calcula el mayor diametro de un conjunto de puntos probando todas
# las combinaciones de puntos posibles

def DiametroFuerzaBruta(P):
    if P == []: return 0
    maxd = 0; dia = []
    for i in range(len(P)):
        for j in range(i,len(P)):
            auxd = dist(P[i],P[j])
            if maxd < auxd:
                maxd = auxd
                dia = [P[i],P[j]]
    return maxd

# Calcula el mayor diametro de un conjunto de puntos
# calculando su cierre convexo y comprobando cada par de
# vertices de este.

def DiametroCierreConvexo(P):
    if P == []: return 0
    Q = Graham(P)
    maxd = 0; dia = []
    for i in range(len(Q)):
        for j in range(i,len(Q)):
            auxd = dist(Q[i],Q[j])
            if maxd < auxd:
                maxd = auxd
                dia = [Q[i],Q[j]]
    return maxd

# Calcula el diametro por comparacion de los angulos
# que forman las aristas del cierre convexo

# fuente: David Munuera Mazarro

def cosUV(u,v):
    return (u[0]*v[0] + u[1]*v[1])/(sqrt(u[0]**2 + u[1]**2)*sqrt(v[0]**2 + v[1]**2))

def diametroRotacion(P):
    if len(P) == 0 or len(P) == 1:
        return 0
    elif len(P) == 2:
        return dist(P[0], P[1])
    else:
        dia = 0
        Q = Graham(P)
        imax, imin = 0, 0
        vmin = [0, 1]
        vmax = [0, -1]
        n = len(Q)
        for k in range(len(Q)):
            if Q[k][0] > Q[imax][0]: imax = k 
            if Q[k][0] < Q[imin][0]: imin = k 
        xmin = Q[imin]
        it = True
        while it or Q[imax%n] != xmin:
            dia = max(dia, dist(Q[imax%n], Q[imin%n]))
            v_max = [Q[(imax+1)%n][0] - Q[imax%n][0], Q[(imax+1)%n][1] - Q[imax%n][1]]
            v_min = [Q[(imin+1)%n][0] - Q[imin%n][0], Q[(imin+1)%n][1] - Q[imin%n][1]]
            c_max = cosUV(vmax, v_max)
            c_min = cosUV(vmin, v_min)
            if c_max >= c_min:
                vmax = v_max
                vmin = [-vmax[0], -vmax[1]]
                imax += 1
            else:
                vmin = v_min
                vmax = [-vmin[0], -vmin[1]]
                imin += 1
                if it: it = False
        return dia

# Metodo de triangulacion basado en el calculo del cierre
# convexo de graham

def GrahamTriangulation(P):
    Q = angularSort(P,max(P))
    T = []
    
    for i in range(1, n-1):
        T.append([Q[0],Q[i],Q[i+1]])
    
    i = 1;
    while Q[i] != Q[-1]:
        if sarea(Q[i-1],Q[i],Q[(i+1)%len(Q)]) < 0:
            T.append([Q[i-1],Q[i],Q[(i+1)%len(Q)]])
            del Q[i]; 
            i -= 1
        else: i += 1
        
    return Q


############# GTC 4 #################

# funcion que crea un DCEL, para un poligono

# Fuente: Manuel Abellanas Oar

def dcel(P):
    n=len(P)
    V=[[P[i],i] for i in range(len(P))]
    e=[[i,n+i,(i-1)%n,(i+1)%n,1]for i in range(n)]+[[(i+1)%n,i,n+(i+1)%n,n+(i-1)%n,0]for i in range(n)]
    f=[n,0]
    return [V,e,f]

# funciones para referirse a los elementos asociados a un elemento del DCEL

# indice del origen de una arista e

# Fuente: Manuel Abellanas Oar
def origin(e,D):
    return D[1][e][0]

# coordenadas del origen de la arista e

# Fuente: Manuel Abellanas Oar
def originCoords(e,D):
    return D[0][origin(e,D)][0]

# arista gemela de la arista e

# Fuente: Manuel Abellanas Oar    
def twin(e,D):
    return D[1][e][1]

# arista previa de la arista e

# Fuente: Manuel Abellanas Oar    
def prev(e,D):
    return D[1][e][2]

# arista siguiente de la arista e

# Fuente: Manuel Abellanas Oar
def next(e,D):
    return D[1][e][3] 

# indice de la cara de cuyo borde forma parte la arista e

# Fuente: Manuel Abellanas Oar
def face(e,D):
    return D[1][e][4]               

# indice de una de las aristas del borde de la cara c

# Fuente: Manuel Abellanas Oar
def edge(c,D):
    return D[2][c]
        
# funcion para dibujar las aristas de un DCEL

# Fuente: Manuel Abellanas Oar
def tercio(A,B):
    return [.7*A[0] + .3*B[0], .7*A[1] + .3*B[1]]

# Funcion para dibujar un DCEL en el que se indican los
# indices de las aristas y vertices

# Fuente: Manuel Abellanas Oar

def plotDCELN(D):
    aristas = sum(line([originCoords(i,D),originCoords(twin(i,D),D)],aspect_ratio=1,thickness=.3) for i in range(len(D[1])))
    aristas += sum(text(i, tercio(originCoords(i,D),originCoords(twin(i,D),D)),color="red") for i in range(len(D[1])))
    vertices = point([i[0] for i in D[0]]) + sum(text(i, vector(D[0][i][0])+vector([.05,.05])) for i in range(len(D[0])))
    return aristas + vertices

# Funcion para dibujar un DCEL

# Fuente: Manuel Abellanas Oar

def plotDCEL(D):
    return sum(line([originCoords(i,D),originCoords(twin(i,D),D)],aspect_ratio=1) for i in range(len(D[1])))  
    
# funcion para colorear una cara de un DCEL

# Fuente: Manuel Abellanas Oar

def plotFace(c,D,col):

    f=D[2][c]
    C=[f]
    f=next(f,D)
    while f <> C[0]:
        C.append(f)
        f=next(f,D)
    
    P=[originCoords(j,D) for j in C]
    return polygon(P,color=col, alpha=.5)     

# funcion para colorear las caras de un DCEL

# Fuente: Manuel Abellanas Oar
    
def colorDCEL(D):
    return sum(plotFace(i,D,(random(),random(),random())) for i in range(1,len(D[2])))

# funcion para dividir una cara del DCEL D por una diagonal
# e1 y e2 son las aristas cuyos origenes son los extremos de la diagonal que divide la cara

# Fuente: Manuel Abellanas Oar

def splitFace(e1,e2,D):
    # si no son aristas de la misma cara o si son adyacentes sus origenes no definen una diagonal
    if face(e1,D) <> face(e2,D) or origin(e2,D) == origin(twin(e1,D),D) or origin(e1,D) == origin(twin(e2,D),D):
        print "no diagonal"
        return
    
    nv, ne, nf = len(D[0]), len(D[1]), len(D[2])
    preve1 = prev(e1,D)
    preve2 = prev(e2,D)
    k=face(e1,D)
    
    # anyadimos las aristas nuevas
    D[1].append([origin(e1,D),ne+1,preve1,e2,k])
    D[1].append([origin(e2,D),ne,preve2,e1,nf])
    
    # modificamos aristas afectadas
    D[1][preve1][3]=ne
    D[1][e1][2]=ne+1
    D[1][preve2][3]=ne+1
    D[1][e2][2]=ne
    i=e1
    while i<>ne+1:
        D[1][i][4]=nf
        i=next(i,D)
    
    #modificamos la cara afectada
    D[2][k]=ne
    
    # anyadimos la nueva cara
    D[2].append(ne+1)

# Obtiene la lista de los indices de las aristas de la cara c del DCEL
# en el orden en el que se recorren
    
def faceEdges (c,D):
    C = [edge (c,D)]
    aux = next(C[-1],D)
    while aux <> C[0]:
        C.append(aux)
        aux = next (aux,D)
    return C

# Obtiene la lista de los indices de los vertices de la cara c del DCEL
# en el orden en el que se recorren

def faceVertices(c,D):
    V = faceEdges(c,D)
    return [origin(i,D) for i in V]

# Obtiene el poligono borde de la cara c del DCEL D conservando la orientacion

def faceVerticesCoords(c,D):
    V = faceVertices(c,D)
    return[D[0][i][0] for i in V]

# Obtiene los indices de las caras adyacentes a la cara c del DCEL

def faceNeighbors(c,D):
    return [face(twin(e,D),D) for e in faceEdges(c,D)]

# Obtiene la lista de los indices de las aristas del DCEL cuyo
# origen es v

def vertexEdges(v,D):
    E = [D[0][v][1]]
    aux = twin(prev(E[-1],D),D)
    while aux <> E[0]:
        E.append(aux)
        aux = twin(prev(E[-1],D),D)
    return E

# Obtiene la lista de los indices de las caras del DCEL
# que contienen al vertice v

def vertexFan (v,D):
    C = vertexEdges (v,D)
    return [face(i,D) for i in C]

# Triangula el exterior del borde de la cara externa
# iterando la funcion splitFace

def convexHullDCEL(p):
    P = angularSort(p,min(p))
    D = dcel(P)
    e0 = len(D[1]) - 1
    e = next(e0,D)
    while e != e0:
        if sarea(originCoords(e,D), originCoords(next(e,D),D), originCoords(next(next(e,D),D),D)) > 0:
            aux = prev(e,D)
            splitFace(e,next(next(e,D),D),D)
            if aux != e0:
                e = aux
            else:
                e = next(e0,D)
        else: e = next(e,D)
    return D

# Obtiene el DCEL de la triangulacion de un conjunto de puntos p

def triangulation(p):
    D = convexHullDCEL(p)
    for i in range(len(p)-2,1,-1):
        splitFace(0,i,D)
    return D

# Varia minimamente los parametros

def shake(p):
    return [[i[0]+0.0000001*random(),i[1]+.0000001] for i in p]

# Dado el DCEL de una triangulacion, elimina la arita especificada
# y anyade la arista que une los vertices adyacentes a los vertices
# de la arista eliminada

def flip(a,D):
    oa=D[1][a][0]
    ga=D[1][a][1]
    ca=D[1][a][4]
    aa=D[1][a][2]
    pa=D[1][a][3]
    cb=D[1][ga][4]
    ab=D[1][ga][2]
    pb=D[1][ga][3]
    oga=D[1][ga][0]
    D[1][a]=[D[1][aa][0],ga,pa,ab,ca]
    D[1][ga]=[D[1][ab][0],a,pb,aa,cb]
    D[1][pa][2]=ab
    D[1][pa][3]=a
    D[1][aa][2]=ga
    D[1][aa][3]=pb
    D[1][aa][4]=cb
    D[1][pb][2]=aa
    D[1][pb][3]=ga
    D[1][ab][2]=a
    D[1][ab][3]=pa
    D[1][ab][4]=ca
    D[2][ca]=a
    D[2][cb]=ga
    D[0][oa][1]=pb
    D[0][oga][1]=pa

    return

# Analiza si la arista a es flipable

def flipable (a,D):
    if face(a,D) == 0 or face(twin(a,D),D) == 0:
        return False
    if segmentIntersectionTest([originCoords(a,D),originCoords(twin(a,D),D)],[originCoords(prev(a,D),D),originCoords(prev(twin(a,D),D),D)]):
        return True
    return False

# Analiza si la arista a es legal

def legal(a,D):
    p1 = originCoords(a,D)
    p2 = originCoords(twin(a,D),D)
    p4 = originCoords(prev(a,D),D)
    p3 = originCoords(prev(twin(a,D),D),D)
    if sarea (p1,p2,p3) != 0 and flipable(a,D) and incircle(p1,p2,p3,p4) == 1: #or incircle(p1,p2,p4,p3) == 1 or inclircle(p1,p4,p3,p2) == 1 or incircle(p4,p2,p3,p1) == 1:
        return False
    return True

# Comprueba si una arista es legal y en el caso
# de ser flipable hace un flip con ella, de lo
# contrario se marca como legal

def legalize(T):
    pendientes = [0..len(T[1])-1]
    
    while pendientes:
        e = pendientes.pop()
        if not legal(e,T) :
            flip(e,T)
            pendientes += [prev(e,T),next(e,T),prev(twin(e,T),T),next(twin(e,T),T)]    
    return 

# Calcula la triangulacion de Delaunay de
# un conjunto de puntos p

def Delaunay(p):
    T = triangulation(p)
    legalize(T)
    return T

# Calcula la triangulacion de Delaunay por
# fuerza bruta.

def delone(p):
    T = []
    for i in Combinations(p,3):
        de = True
        for j in p:
            if incircle(i[0],i[1],i[2],j) == 1 and j != i[0] and j != i[1] and j != i[2]:
                de = False
                break
        if de:
            T.append(i)
    return T


########## GTC Practica 6 ##########

# Calcula la region de Voronoi a partir de la
# triangulacion de Delaunay

def VoronoiClipping (p):
    
    D = Delaunay(p)
    
    # Esta funcion devuelve la region de Voronoi
    # del vertice v
    def VoronoiRegion (p, i):
        E = vertexEdges (i,D)
        
        [A , B] = perpendicularBisector(originCoords(E[0],D), originCoords(twin(E[0],D),D))
        v = vector(B) - vector (A)
        A0 = list(vector(A) -1000000*v)
        A1 = list(vector(A) +1000000*v)
        A2 = list(vector(A1) +1000000*vector([-v[1],v[0]]))
        A3 = list(vector(A0) +1000000*vector([-v[1],v[0]]))
        R = [A0, A1, A2, A3]
    
        for k in range(1,len(E)):
            R = clipping(R,perpendicularBisector(originCoords(E[k],D),originCoords(twin(E[k],D),D)))
        
        return R
    
    
    V = [VoronoiRegion(p,i) for i in range(len(p))]
    return V

# Calcula la region de Voronoi de un conjunto de puntos
# p recortado por el rectangulo R

def Voronoi(p):

    def infinityPoint(e,D):
        A = originCoords(e,D)
        B = originCoords(twin(e,D),D)
        u = vector(B)-vector(A)
        v = vector([-u[1], u[0]])
        F = originCoords(prev(twin(e,D),D),D)
        m = vector(circumcenter(A,B,F))
        return list(m+10000*v)
    
    D = Delaunay(p)
    V = []
    for i in range(len(D[0])):
        E = vertexEdges(i,D)
        R = []
        for j in E:
            if face(j,D) == 0:
                R.append(infinityPoint(j,D))
                R.append(infinityPoint(prev(j,D),D))
            else:
                R.append(circumcenter(originCoords(j,D), originCoords(next(j,D),D), originCoords(prev(j,D),D)))
        V.append(R)
        
    return V

# Este metodo reduce el poligono dado acercando las aristas al centro la distancia desp

def reducePoligono(P, desp):
    T = copy(P)
    for i in range(len(P)):
        v = vector(P[i]) - vector(P[i-1])
        v = vector([-v[1],v[0]])
        v /= v.norm()
        p1 = [P[i-1][0] + desp*(v[0]), P[i-1][1] + desp*(v[1])]
        p2 = [P[i][0] + desp*(v[0]), P[i][1] + desp*(v[1])]
        T = clipping(T,[p1,p2])
        p1 = [p2[0],p2[1]]
    return T


########## Metodo de Amenta de reconstruccion de curvas ##########

# Funcion que te devuelve los indices de las aristas mas
# cortas que la distancia r

def shortEdges(D, r):
    E = []
    for i in range(len(D[1])):
        if dist(originCoords(i,D),originCoords(twin(i,D),D)) <= r:
            E.append(i)
    return E

# Funcion que devuelve la lista de vertices del diagrama
# de Voronoi de una triangulacion de Delaunay

def VoronoiVertices(D):
    vertices = []
    for i in range(1, len(D[1])):
        vertices.append(circumcenter(originCoords(i,D), originCoords(prev(i,D),D), originCoords(next(i,D),D)))
    return vertices

# Me devuelve el punto en el infinito que junto con el
# circumcentro determina la region de Voronoi de una
# cara de la triangulacion de Delaunay

def infinityPoint(e,D):
        A = originCoords(e,D)
        B = originCoords(twin(e,D),D)
        u = vector(B)-vector(A)
        v = vector([-u[1], u[0]])
        F = originCoords(prev(twin(e,D),D),D)
        m = vector(circumcenter(A,B,F))
        return list(m+10000*v)

# Devuelve la arista dual de la arista i,
# es decir la arista que determina la region
# de Voronoi de uno de los posibles triangulos
# del que forma parte esa arista

def dual(e,D):
    A = originCoords(e,D)
    B = originCoords(twin(e,D),D)
    if face(e,D) == 0:
        C = infinityPoint(e,D)
        F = circumcenter(A,B,originCoords(prev(twin(e,D),D),D))
    elif face(twin(e,D),D) == 0:
        C = circumcenter(A,B,originCoords(prev(e,D),D))
        F = infinityPoint(twin(e,D),D)
    else:
        C = circumcenter(A,B,originCoords(prev(e,D),D))
        F = circumcenter(A,B,originCoords(prev(twin(e,D),D),D))
    return [C,F]    

# Determina si la arista e forma parte del crust 

def crustTest(e,D):
    A = originCoords(e,D)
    B = originCoords(twin(e,D),D)
    C = dual(e,D)
    segmentIntersectionTest([A,B],[A,B])
    b1 = segmentIntersectionTest([A,B],C)
    return  b1 and incircle(A,B,C[0],C[1]) == -1

# Calcula el crust de la triangulacion de Delaunay D.

def crust(D):
    return [i for i in range(len(D[1])) if crustTest(i,D)]

# Dibuja las aristas en el color especificado

def plotEdges(C, D, col):
    return sum(line([originCoords(c,D),originCoords(twin(c,D),D)], color=col) for c in C)



########## GTC Practica 8 ##########

# Crea las listas de vertices, aristas y caras estas
# dos ultimas estaran ordenadas por su longitud o la
# de su radio de la circumferencia circunscrita a
# partir del DCEL D
def chainComplex (D):
    vertices = [v for [v, i] in D[0]]
    
    # e1 y e2 almacenaran los indices de los vertices de la arista actual
    # almaceno los indices de menor a mayor indice
    e1 = origin(i,D)
    e2 = origin(twin(i,D),D)
    # contiene los pares de indices de los vertices que contienen las aristas
    edges = [[min(e1,e2),max(e1,e2)]]
    # lista que contiene el indice de la arista en el DCEL
    index = [0]
    # lon almacena las longitudes de las aristas, en una lista que sige el orden de edges
    lon = [dist(vertices[e1],vertices[e2])]
    # En este bucle se insertaran de forma ordenada las aristas
    for i in range(1,len(D[1])):
        e1 = origin(i,D)
        e2 = origin(twin(i,D),D)
        # Si la arista gemela esta en edges, pasa a la siguiente arista
        if [min(e1,e2),max(e1,e2)] in edges:
            continue;
        # lon_e almacena la longitud de la arista e
        lon_e = dist(vertices[e1],vertices[e2])
        
        edges.append([e1,e2])
        index.append(i)
        lon.append(lon_e)
        j = len(edges)-1
        while j > 0 and lon[j-1] > lon_e:
            edges[j] = edges[j-1]
            index[j] = index[j-1]
            lon[j] = lon[j-1]
            j -= 1
        edges[j] = [min(e1,e2),max(e1,e2)]
        index[j] = i
        lon[j] = lon_e
    
    # Estos son los puntos de la primera cara
    A = origin(edge(1,D),D)
    B = origin(twin(edge(1,D),D),D)
    C = origin(prev(edge(1,D),D),D)
    # f guarda las caras
    f = [edges.index([min(A,B),max(A,B)]),edges.index([min(B,C),max(B,C)]),edges.index([min(C,A),max(C,A)])]
    # rad es el array de radios
    rad = [dist(D[0][A][0],circumcenter(D[0][A][0],D[0][B][0],D[0][C][0]))]
    faces = [f]
    for i in range(2,len(D[2])):
        A = origin(edge(i,D),D)
        B = origin(twin(edge(i,D),D),D)
        C = origin(prev(edge(i,D),D),D)
        f = [edges.index([min(A,B),max(A,B)]),edges.index([min(B,C),max(B,C)]),edges.index([min(C,A),max(C,A)])]
        rad_f = dist(D[0][A][0],circumcenter(D[0][A][0],D[0][B][0],D[0][C][0]))
        faces.append(f)
        rad.append(rad_f)
        j = len(faces)-1
        while j > 0 and rad[j-1] > rad_f:
            faces[j] = faces[j-1]
            rad[j] = rad[j-1]
            j -= 1
        faces[j] = f
        rad[j] = rad_f
        
    return [vertices, edges, faces]

# Crea las matrices de las aplicaciones borde

# operadores borde de dimensiones 3,2,1,0 del complejo simplicial C
# devuelve la matriz de la aplicacion lineal
def borderOperator(i,C):
    if i == 0:
        return matrix(GF(2),1,len(C[0]),[0 for i in range(len(C[0]))])
    elif i == 1:
        return matrix(GF(2),len(C[0]),len(C[1]),[[1 if i in C[1][j] else 0 for j in range(len(C[1]))] for i in range(len(C[0]))])
    elif i == 2:
        return matrix(GF(2),len(C[1]),len(C[2]),[[1 if i in C[2][j] else 0 for j in range(len(C[2]))] for i in range(len(C[1]))])
    elif i == 3:
        return matrix(GF(2),len(C[2]),1,[[0] for i in range(len(C[2]))])
    return


# C[0] lista de vertices
# C[1] lista de tuplas de indices de los vertices
# C[2] lista de triplas de indices de las aristas
# C[0][] vertice en cuestion
# C[1][] tupla en cuestion
# C[2][] tripla en cuestion
# C[1][][] indice de un vertice de la arista -> C[0][C[1][][]] uno de los vertices de una arista
# C[2][][] indice de una de las tuplas -> C[1][C[2][][]] la tupla en cuestion
# C[1][C[2][][]][] indice de un vertice de la arista -> C[0][C[1][C[2][][]][]] un vertice de un triangulo

# Esta funcion devuelve los indices de los vertices de un triangulo i de un alfa-complejo C
def getFace (i, C):
    X = C[1][C[2][i][0]][0]
    Y = C[1][C[2][i][0]][1]
    Z = C[1][C[2][i][1]][0]
    if Z == X or Z == Y:
        Z = C[1][C[2][i][1]][1]
    return [X,Y,Z]

# vamos a definir los alpha-complejos dando las listas de aristas y de triangulos
# por sus indices en las listas correspondientes del complejo de cadenas

# En este metodo he anyadido el parametro permiteCuerdas que voy a utilizar para 
# indicar si quiero que el alfa complejo contenga los indices de las aristas
# que no se encuentran en ninguna cara. Esta opcion me ayuda a visualizar y 
# las figuras de ejemplo con mayor precision y a obtener los numeros de Betty 
# que estas deben tener
def alphaComplex(a, C, permiteCuerdas):
    E = []
    i = 0;
    while i < len(C[1]) and dist(C[0][C[1][i][0]],C[0][C[1][i][1]]) <= a*2:
        E.append(i)
        i += 1
    
    j = 0;
    F = []
    while j < len(C[2]):
        if C[2][j][0] >= i or C[2][j][1] >= i or C[2][j][2] >= i:
            j += 1
            continue;
            
        t = getFace (j, C)
        if dist(C[0][t[0]], circumcenter(C[0][t[0]], C[0][t[1]], C[0][t[2]])) <= a:
            F.append(j)
            j += 1
        else:
            break;
    
    # Si el parametro permiteCuerdas == false
    # Construyo una nueva lista de aristas que no contiene los indices de las aristas
    # que no forman parte de ningun triangulo
    if not permiteCuerdas:
        E1 = []
        for i in F:
            for j in range(3):
                if not C[2][i][j] in E1:
                    E1.append(C[2][i][j])
        E = E1
        # Esta parte desordena las edges, pero no afectara
    
    return E,F

# recibo el parametro permiteCuerdas para pasarselo al metodo alphaComplex
def plotAlphaComplex(a,p,permiteCuerdas):
    D = Delaunay(p)
    C = chainComplex(D)
    E,F = alphaComplex(a,C,permiteCuerdas)
    
    g = point(C[0])
    g += sum(line([C[0][C[1][E[i]][0]], C[0][C[1][E[i]][1]]]) for i in range(len(E)))

    for j in F:
        t = getFace(j,C)
        g += polygon([C[0][t[0]],C[0][t[1]],C[0][t[2]]], alpha=.3)

    return g

# Pinta las caras en colores
def plotAlphaComplexColor(a,p,permiteCuerdas):
    D = Delaunay(p)
    C = chainComplex(D)
    E,F = alphaComplex(a,C,permiteCuerdas)
    
    g = point(C[0])
    g += sum(line([C[0][C[1][E[i]][0]], C[0][C[1][E[i]][1]]]) for i in range(len(E)))

    for j in F:
        t = getFace(j,C)
        g += polygon([C[0][t[0]],C[0][t[1]],C[0][t[2]]], color=(random(),random(),random()), alpha=.6)

    return g

# operador borde del alpha-complejo A del complejo de cadenas C
def borderOperatorAlpha(i,C,A): # i = 0 ,1 ,2 ,3
    if i == 0:
        return matrix(GF(2),1,len(C[0]),[0 for i in range(len(C[0]))])
    elif i == 1:
        return matrix(GF(2),len(C[0]),len(A[0]), [[1 if i in C[1][j] else 0 for j in A[0]] for i in range(len(C[0]))])
    elif i == 2:
        return matrix(GF(2),len(A[0]),len(A[1]), [[1 if i in C[2][j] else 0 for j in A[1]] for i in A[0]])
    elif i == 3:
        return matrix(GF(2),len(A[1]),1, [[0] for i in range(len(A[1]))])
    return 

# Calcula el alfa complejo a partir de un conjunto de puntos p
# y devuelve sus numeros de Betti

def BettiNumbersAlphaComplex(p,a,permiteCuerdas):
    D = Delaunay(p)
    C = chainComplex(D)
    A = alphaComplex(a,C,permiteCuerdas)
    
    Im0 = borderOperatorAlpha(0,C,A).rank(); ker0 = len(C[0])
    Im1 = borderOperatorAlpha(1,C,A).rank(); ker1 = len(A[0]) - Im1;
    Im2 = borderOperatorAlpha(2,C,A).rank(); ker2 = len(A[1]) - Im2;
    Im3 = borderOperatorAlpha(3,C,A).rank(); ker3 = Im3;
    
    b0 = ker0 - Im1
    b1 = ker1 - Im2
    b2 = ker2 - Im3
    
    print "B0 = ", b0, ". Hay ", b0, "componentes conexas."
    print "B1 = ", b1, ". Hay ", b1, "agujeros."
    print "B2 = ", b2, ". Hay ", b2, "cavidades."
    
    return [b0,b1,b2]
