import math
import array

def absoluteError (a : float, b: float) :
    c : float = abs(a - b)
    return c

def relativeError(a : float, b : float) :
    c : float = (absoluteError(a,b) / abs(compareAfterComa(a, b)))
    return c

def vectorNorma(a : array) :
    nv1 : any = 0

    for i in a :
        nv1 += pow(i,2)

    return math.sqrt(nv1)

def VectScallMult(v1 : array, v2 : array) :
    sdv : any = 0
    if len(v1) != len(v2) :
        return "Неможна шукати скалярний добуток векторів різної розмірності "

    for i in range(len(v1)) :
        sdv += (v1[i] * v2[i])

    return sdv

def higestLen(a : float, hLen: float) :
    count : int
    tStr = str(hLen)

    if "." in tStr :
        count = len(tStr.split(".")[1])
    else :
        count = 0
    return round(a, count)

def compareAfterComa(a : float, b : float) :
    lenA : int = len(str(a).split(".")[1])
    lenB : int = len(str(b).split(".")[1])

    if lenA > lenB :
        return a
    elif lenB > lenA :
        return b
    else :
        return a


#Похибки
p1 : float = 3.14159265
p2 : float = 3.14

asnAblolute : float = higestLen(absoluteError(p1,p2), compareAfterComa(p1, p2))
asnRelative : float = higestLen(relativeError(p1,p2), compareAfterComa(p1, p2))

print(f" Абсолютна похибка : {asnAblolute}\n Відносна похибка : {asnRelative}")


#Норма вектора
vect1 = array.array('f', [2,2,2])
print(f" Норма вектора ({', '.join(str(x) for x in vect1)}) = {vectorNorma(vect1)}")


#Скалярний добуток векторів
vectX = array.array('f', [1,2,3])
vectY = array.array('f', [2,3,4])
print(f" Скалярний добуток векторів x({', '.join(str(x) for x in vectX)}) та y({', '.join(str(x) for x in vectY)}) = {VectScallMult(vectX, vectY)}")