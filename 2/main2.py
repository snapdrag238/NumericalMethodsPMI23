import math
import array
import re

EPSILON = 1e-12
MAX_ITERATIONS = 50
INITIAL_GUESS = [0.0, 0.0, 0.0]

#варіант 15
# -7x1 +3x2 +2x3 = 13
# -3x1 -7x2 -2x3 = 25       
# -x1 -2x2 -6x3 = 13

class Slar :
    def __init__(self, nInp: int) :
        self.n = nInp
        self.matrixA = [[0.0] * self.n for _ in range(self.n)]
        self.matrixB = [0.0] * self.n
        self.matrixAns = [[0.0] * self.n for _ in range(self.n)]

    def fillMatrix(self,inpType : str) -> bool:
        match inpType:
            case "file" :
                print("Введіть ім'я файлу:", end="")
                filename = str(input())

                while(not filename.endswith(".txt")) :
                    print("Такий формат файлу не підтримується , спробуйте відкрити файл формату .txt")
                    filename = str(input())

                self.matrixFileRead(filename)

                return True
            case "cmd":
                for i in range(self.n) :
                    for j in range(self.n) :
                        print(f"Введіть префіксну частину невідомого x{i+1}{j+1} >> ",end="")
                        self.matrixA[i][j] = float(input())

                for i in range(self.n) :
                    print(f"Введіть вільний член b{i+1} >> ", end="")
                    self.matrixB[i] = float(input())
                
                return True
            case _ :
                print(f"Варіант вибору {inpType} недійсний")
                return False

    def calcVectResidual(self, x : list[float]) -> tuple[list[float], float] :  #обч вектор невязки r = b - Ax, його неск форму ||r||inf
        r = [.0] * self.n
        for i in range(self.n) : 
            Axi = sum(self.matrixA[i][j] * x[j] for j in range(self.n))
            r[i] = self.matrixB[i] - Axi
        normR = max(abs(val) for val in r)
        return r, normR

    def checkSuffConditionConv(self) -> bool :  #перевірка достатньої умови збіжності (|a_ii| > sum_{j!=i} |a_ij|)
        striclDomm = True
        for i in range(self.n) :
            diag = abs(self.matrixA[i][i])
            offDiag = sum(abs(self.matrixA[i][j]) for j in range(self.n) if j != i)
            if diag <= offDiag:
                striclDomm = False
        return striclDomm


    def methodGaus(
        self, eps: float = EPSILON
    ) -> tuple[bool, list[float] | None]:
        n = self.n
        A = [row[:] for row in self.matrixA]
        b = self.matrixB[:]

        swapCount = 0
        det = 1.0

        for k in range(n):
            #пошук макс ел за модулем у k-му стовпці
            maxRow = k
            maxVal = abs(A[k][k])
            for i in range(k + 1, n):
                if abs(A[i][k]) > maxVal:
                    maxVal = abs(A[i][k])
                    maxRow = i

            #перевірка на виродженість
            if maxVal < eps:
                print(
                    f"Помилка!!! Матриця вироджена або близька до виродженої (|a_kk| < {eps:.1e})."
                )
                return False, None

            #перестановка рядків
            if maxRow != k:
                A[k], A[maxRow] = A[maxRow], A[k]
                b[k], b[maxRow] = b[maxRow], b[k]
                swapCount += 1

            #накопичення визначника
            det *= A[k][k]

            #виключення невідомих
            for i in range(k + 1, n):
                fct = A[i][k] / A[k][k]
                A[i][k] = 0.0
                for j in range(k + 1, n):
                    A[i][j] -= fct * A[k][j]
                b[i] -= fct * b[k]

        #обч знаку визначника
        det *= (-1) ** swapCount
        print(
            "\nРезультат за методом Гауса:\n"
            f"Кількість перестановок рядків: {swapCount}\n"
            f"Визначник матриці (det A): {det:.6f}"
        )

        #зворотний хід
        x = [0.0] * n
        for i in range(n - 1, -1, -1):
            sumAx = sum(A[i][j] * x[j] for j in range(i + 1, n))
            x[i] = (b[i] - sumAx) / A[i][i]

        print("\nРозвязок:")
        for i in range(n):
            print(f"x{i + 1} = {x[i]:.6f}")

        r, normr = self.calcVectResidual(x)
        print(f"Норма невязки ||r||inf = {normr:.6e}")

        return True, x
            

    def methodYacobi(self) -> bool:
        pass

    def methodZeidel(self) -> bool:
        pass

    def simpleMatrixOutput(self):
        print("Тип matrixA:", type(self.matrixA))
        print("Вміст matrixA:", self.matrixA)
        print("Тип matrixB:", type(self.matrixB))
        print("Вміст matrixB:", self.matrixB)
        for i in range(self.n):
            print(" | ", end="")
            for j in range(self.n):
                print(f"{self.matrixA[i][j]:.2f}", end=" ")
            print(f" | {self.matrixB[i]:.2f}")




    def finalMatrixOutput(self):
        pass

    def matrixFileRead(self,fname : str) -> bool:
        try:
            with open(fname, "r", encoding="utf-8") as f:
                        lines = f.readlines()
        except FileNotFoundError:
            print(f"Помилка!!! Файл {fname} не знайдено.")
            return False

        for i in range(self.n):
            row, self.matrixB[i] = parseEquastion(lines[i])
            self.matrixA[i] = [row[var] for var in sorted(row.keys())]

        return True
            
def parseEquastion(eqStr : str):
        leftPart, rightpart = eqStr.split("=")

        pattern = r"([+-]?\s*\d*)\s*(x\d+)"
        rhs = float(rightpart.strip())

        coefs = {}
        for rCoefs, varName in re.findall(pattern, leftPart):
            clear = rCoefs.replace(" ","")

            if clear in ("", "+"):
                value = 1.0
            elif clear == "-":
                value = -1.0
            else:
                value = float(clear)

            coefs[varName] = value

        return coefs, rhs

def main() :
    slar = Slar(3)
    print("Виберіть як вихочети заповнити СЛАР" \
            "1. .txt файд введіть 'file'" \
            "2. Вручну через консоль 'cmd'" \
            ">>", end="")
    slar.fillMatrix(str(input()))
    #slar.simpleMatrixOutput()
    slar.methodGaus()



if __name__ == "__main__" :
    main()  