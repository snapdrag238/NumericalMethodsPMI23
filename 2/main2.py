import re

DEF_DIM = 3
EPSILON_GAUS = 1e-12
EPSILON_ITER = 1e-3
MAX_ITERATIONS = 100

#варіант 15
# -7x1 +3x2 +2x3 = 13
# -3x1 -7x2 -2x3 = 25       
# -x1 -2x2 -6x3 = 13

class Slar :
    def __init__(self, nInp: int = None) :
        if nInp is None :
            self.n = DEF_DIM
        else:
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

                return self.matrixFileRead(filename)

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
            

    def calcVectResidual(self, x : list[float]) -> tuple[list[float], float] :
        r = [.0] * self.n
        for i in range(self.n) : 
            Axi = sum(self.matrixA[i][j] * x[j] for j in range(self.n))
            r[i] = self.matrixB[i] - Axi
        normR = max(abs(val) for val in r)
        return r, normR
    

    def checkSuffConditionConv(self) -> bool :
        striclDomm = True
        for i in range(self.n) :
            diag = abs(self.matrixA[i][i])
            offDiag = sum(abs(self.matrixA[i][j]) for j in range(self.n) if j != i)
            if diag <= offDiag:
                striclDomm = False
        return striclDomm


    def methodGaus(
        self, eps: float = EPSILON_GAUS
    ) -> tuple[bool, list[float] | None]:
        n = self.n
        A = [row[:] for row in self.matrixA]
        b = self.matrixB[:]

        swapCount = 0
        det = 1.0

        for k in range(n):
            maxRow = k
            maxVal = abs(A[k][k])
            for i in range(k + 1, n):
                if abs(A[i][k]) > maxVal:
                    maxVal = abs(A[i][k])
                    maxRow = i

            if maxVal < eps:
                print(
                    f"Помилка!!! Матриця вироджена або близька до виродженої (|a_kk| < {eps:.1e})."
                )
                return False, None

            if maxRow != k:
                A[k], A[maxRow] = A[maxRow], A[k]
                b[k], b[maxRow] = b[maxRow], b[k]
                swapCount += 1

            det *= A[k][k]

            for i in range(k + 1, n):
                fct = A[i][k] / A[k][k]
                A[i][k] = 0.0
                for j in range(k + 1, n):
                    A[i][j] -= fct * A[k][j]
                b[i] -= fct * b[k]

        det *= (-1) ** swapCount
        print(
            "\nРезультат за методом Гауса:\n"
            f"Кількість перестановок рядків: {swapCount}\n"
            f"Визначник матриці (det A): {det:.6f}"
        )

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
            

    def methodYacobi(self, eps : float = EPSILON_ITER, maxIter : int = MAX_ITERATIONS, x0 : list[float] | None = None) -> tuple[int, list[float]] :
        n = self.n
        xPrev = [0.0] * n if x0 is None else x0[:]
        xCurr = [0.0] * n

        if self.checkSuffConditionConv():
            print(
                "\n[Якобі] Перевірка збіжності: матриця МАЄ суворе діагональне переважання. Збіжність гарантована."
            )
        else:
            print(
                "\n[Якобі] Матриця НЕ МАЄ суворого діагонального переважання! Збіжність не гарантується."
            )

        print("\nМетод Якобі:")
        colNames = " | ".join(f"{f'x{i+1}':>12}" for i in range(n))
        header = f" Ітер | {colNames} | {'||x(k)-x(k-1)||':>16}"
        print("-" * len(header))
        print(header)
        print("-" * len(header))

        iters = 0
        conv = False
        for it in range(1, maxIter + 1) :
            iters = it
            for i in range(n) :
                s = sum(self.matrixA[i][j] * xPrev[j] for j in range(n) if j != i)
                xCurr[i] = (self.matrixB[i] - s) / self.matrixA[i][i]

            diffNorm = max(abs(curr - prev) for curr, prev in zip(xCurr, xPrev))

            rowVals = " | ".join(f"{val:12.6f}" for val in xCurr)
            print(f"{it:5d} | {rowVals} | {diffNorm:16.6e}")

            if diffNorm < eps:
                conv = True
                break
            xPrev = xCurr[:]

        print("-" * len(header))

        if not conv:
            print(
                f"[Якобі] Досягнуто ліміту {maxIter} ітерацій без досягнення заданої точності eps={eps}!"
            )

        print(f"Виконано ітреацій: {iters}")
        _, normr = self.calcVectResidual(xCurr)
        print(f"Норма невязки ||r||inf = {normr:.6e}")

        return iters, xCurr
    

    def methodZeidel(self, eps : float = EPSILON_ITER, maxIter : int = MAX_ITERATIONS, x0 : list[float] | None = None) -> tuple[int, list[float]]:
        n = self.n
        x = [0.0] * n if x0 is None else x0[:]

        if self.checkSuffConditionConv():
            print(
                "\n[Зейдель] Перевірка збіжності: матриця МАЄ суворе діагональне переважання. Збіжність гарантована."
            )
        else:
            print(
                "\n[Зейдель] Матриця НЕ МАЄ суворого діагонального переважання! Збіжність не гарантується.")

        print("\nМетод Зейдаля:")
        colNames = " | ".join(f"{f'x{i+1}':>12}" for i in range(n))
        header = f" Ітер | {colNames} | {'||x(k)-x(k-1)||':>16}"
        print("-" * len(header))
        print(header)
        print("-" * len(header))

        iters = 0
        conv = False
        for it in range(1, maxIter + 1) :
            iters = it
            maxDiff = 0.0

            for i in range(n) :
                s = sum(self.matrixA[i][j] * x[j] for j in range(n) if j != i)
                newXi = (self.matrixB[i] - s) / self.matrixA[i][i]

                diff = abs(newXi - x[i])
                if diff > maxDiff:
                    maxDiff = diff
                x[i] = newXi

            rowVals = " | ".join(f"{val:12.6f}" for val in x)
            print(f"{it:5d} | {rowVals} | {maxDiff:16.6e}")

            if maxDiff < eps :
                conv = True
                break

        print("-" * len(header))

        if not conv:
            print(
                f"[Зейдель] Досягнуто ліміту {maxIter} ітерацій без досягнення заданої точності eps={eps}!"
            )
        
        print(f"Виконано ітерацій : {iters}")
        _, normr = self.calcVectResidual(x)
        print(f"Норма невязки ||r||inf = {normr:.6e}")

        return iters, x
    

    def simpleMatrixOutput(self):
       # print("Тип matrixA:", type(self.matrixA))
        #print("Вміст matrixA:", self.matrixA)
        #print("Тип matrixB:", type(self.matrixB))
        #print("Вміст matrixB:", self.matrixB)
        for i in range(self.n):
            print(" | ", end="")
            for j in range(self.n):
                print(f"{self.matrixA[i][j]:8.2f}", end=" ")
            print(f" | {self.matrixB[i]:8.2f}")


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


def printMenu() :
    print("\n\nГоловне меню:\n" \
    "1. Зчитати матрицю\n" \
    "2. Вивести поточну матрицю\n" \
    "3. Розвязання методом Гауса\n" \
    "4. Розвязання методом Якобі\n" \
    "5. Розвязання методом Зейдаля\n" \
    "6. Послідовне розвязання трьома методами\n" \
    "0. Вийти\n" \
    ">>", end=" ")


def getInitialGuess(n : int) -> list[float] | None :
    ans = input("Використати нульове початкове наближення x^(0)? (y/n) [за замовч. y]:").strip().lower()

    if ans == "n" :
        print(f"Введіть {n} через пробіл для x^(0) :")
        while True :
            try :
                raw = input(">> ").split()
                if len(raw) != n:
                    print(f"Потрібно ввести рівно {n} чисел!")
                    continue
                return [float(val) for val in raw]
            except ValueError :
                print("Помилка: введіть коректні дійсні числа!")
    return None

def main() :

    slar = Slar()
    m_isMatrixLoaded = False
    m_initChoise = input("Перед початком роботи потрібно завантажити матрицю"
    "Оберіть джерело даних ('file' - завантаження з .txt файлу, 'cmd' ручний ввід з консолі)\n >> ")

    if slar.fillMatrix(m_initChoise) :
        m_isMatrixLoaded = True
        slar.simpleMatrixOutput()
    
    while True :
        printMenu()
        m_choise = input()

        if m_choise == "0" :
            print("Завершення роботи...")
            break
        
        if m_choise == "1" :
            m_inpMode = input("Оберіть джерело даних ('file' - завантаження з .txt файлу, 'cmd' ручний ввід з консолі)\n >> ")
            if slar.fillMatrix(m_inpMode) :
                m_isMatrixLoaded = True
                slar.simpleMatrixOutput()

        if not m_isMatrixLoaded:
            print("\nУвага!!! Матриця ще не завантажена! Спочатку оберіть пункт 1.")
            continue
                
        match str(m_choise) :
            case "2" :
                slar.simpleMatrixOutput()

            case "3" :
                slar.methodGaus()

            case "4" :
                x0 = getInitialGuess(slar.n)
                slar.methodYacobi(x0 = x0)

            case "5" :
                x0 = getInitialGuess(slar.n)
                slar.methodZeidel(x0 = x0)

            case "6" :
                slar.methodGaus()
                slar.methodYacobi(x0 = None)
                slar.methodZeidel(x0 = None)

            case _:
                print("\nНевірний пункт меню! Введіть цифру від 0 до 6 включно.")


if __name__ == "__main__" :
    main()