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
        return "Неправильна довжина векторів"

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

def main():

    print(
    "Меню:\n" \
    "1.Обчислення Абсолютної та Відносної похибки\n" \
    "2.Обчислення порми n-вимірного вектора\n" \
    "3.Обчислення Скалярного добутку n-вимірного ветора\n" \
    "0.Вихід\n" \
    ">> ", end="")

    val = input()

    match str(val):
        case "1":
            print("Введіть 1 число: ", end="")
            p1 = float(input())
            print("Введіть 2 число: ", end="")
            p2 = float(input())

            asnAblolute = higestLen(absoluteError(p1,p2), compareAfterComa(p1,p2))
            asnRelative = higestLen(relativeError(p1,p2), compareAfterComa(p1,p2))

            print(f" Абсолютна похибка : {asnAblolute}\n Відносна похибка : {asnRelative}")

            return
        case "2":
            print("Введіть вектор (Наприклад : 1, 1, 1): ", end= "")
            vectInput = input()

            vect = array.array("f",[float(x) for x in vectInput.split(',')])
            print(f" Норма вектора ({', '.join(str(x) for x in vect)}) = {vectorNorma(vect)}")

            return
        case "3":
            print("Введіть вектор X (Наприклад : 1, 1, 1): ", end= "")
            vectX = array.array('f',[float(x) for x in input().split(',')])

            print("Введіть вектор Y (Наприклад : 1, 1, 1): ", end= "")
            vectY = array.array('f',[float(x) for x in input().split(',')])

            if len(vectX) != len(vectY) :
                print("Помилка!!! Неможна шукати скалярний добуток векторів різної розмірності")
            else :
                print(f"Скалярний добуток векторів x({', '.join(str(x) for x in vectX)}) та y({', '.join(str(x) for x in vectY)}) = {VectScallMult(vectX, vectY)}")
            return
        case "0":
            return 0
        case _:
           return print("Wrong choise!!!")

if __name__ == "__main__" :
    main()  