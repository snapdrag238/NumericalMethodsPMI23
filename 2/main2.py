import math
import array

#варіант 15
# -7x1 + 3x2 + 2x3 = 13
#          
#

class Slar :
    def __init__(self, nInp: int) :
        self.n = nInp
        self.matrixA = [[0] * self.n in range(self.n)]
        self.matrixB = [self.n]

    def fillMatrix(inpType : str) -> bool:
        #print("Виберіть як вихочети заповнити СЛАР" \
        #"1. .txt файд введіть 'file'" \
        #"2. Вручну через консоль 'cmd'" \
        #">>", end="")

        match inpType:
            case "file" :
                print("Введіть ім'я файлу:", end="")
                filename = str(input())
                while(not filename.endswith(".txt")) :
                    print("Такий формат файлу не підтримується , спробуйте відкрити файл формату .txt")
                    filename = str(input())
                # читання слар з файлу

                return True
            case "cmd":
                for i in range(Slar.n) :
                    for j in range(Slar.n) :
                        print(f"Введіть префіікс невідомого x{i}{j} >> ",end="")
                        Slar.matrixA[i][j] = float[input()]
                for i in range(Slar.n) :
                    print(f"Введіть вільний член b{i} >> ", end="")
                    Slar.matrixB[i] = float(input())
                
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



    def main() :
        pass

    if __name__ == "__main2__" :
        main()  