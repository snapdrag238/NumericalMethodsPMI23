import math
import array
import re

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

    def methodGaus(self) -> bool:
        pass

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
                print(f"{self.matrixA[i][j]:6.2f}", end="")
            print(f" | {self.matrixB[i]:6.2f}")




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
            self.matrixA = [row[var] for var in sorted(row.keys())]

        return True
            
def parseEquastion(eqStr : str) :
        leftPart, rightpart = eqStr.split("=")

        pattern = r"([+-]?\s*\d*)\s*(x\d+)"
        rhs = float(rightpart.strip())

        coefs = {}
        for rCoefs, varName in re.findall(pattern, leftPart):
            clear = rCoefs.replace(" ","")

            if clear in ("", "+"):
                value = 1.0
            elif clear in ("", "-"):
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
    slar.simpleMatrixOutput()



if __name__ == "__main__" :
    main()  