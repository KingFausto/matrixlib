import math


# TODO: decomposicion LU

class Matrix:

    def __init__(self, matrix=None, rows=None, cols=None, mtype=None):
        if matrix:
            self.matrix = matrix
        elif rows and cols:
            if mtype == "i" and rows == cols:
                for i in range(rows):
                    for j in range(cols):
                        self.matrix[i][j] = 1 if i == j else 0
            else:
                self.matrix = self.__fill(rows, cols, mtype)

        self.rows = len(self.matrix)
        self.cols = len(self.matrix[0])
        self.size = self.rows, self.cols

    def __fill(self, rows, cols, val):
        return [[val for j in range(cols)] for i in range(rows)]

    def __distance(self, vector_a, vector_b):
        total = 0
        for i in range(vector_a.rows):
            total += math.pow(vector_b.matrix[i][0] - vector_a.matrix[i][0], 2)
        return math.sqrt(total)

    def add(self, val):
        inserted = False
        for i in range(self.rows):
            for j in range(self.rows):
                if self.matrix[i][j] is None:
                    self.matrix[i][j] = val
                    inserted = True
                    break
            if inserted:
                break

    def splice(self, other, col):
        """Insert a column in a specific column of the matrix."""
        if len(other.matrix) == self.cols:
            out = Matrix(rows=self.rows, cols=self.cols)
            for i in range(self.rows):
                for j in range(self.cols):
                    if j == col:
                        out.matrix[i][j] = other.matrix[i][0]
                    else:
                        out.matrix[i][j] = self.matrix[i][j]
            return out

    def resize(self, row, col):
        """Removes a row and column from matrix."""
        out = Matrix(rows=self.rows-1, cols=self.cols-1)
        for i in range(self.rows):
            for j in range(self.cols):
                if i != row and j != col:
                    out.add(self.matrix[i][j])
        return out

    def solve(self, ans, method="determinant", error=3):
        """Solves a system of equations using the determinant method."""
        if method == "determinant":
            out = Matrix(rows=ans.rows, cols=1)
            D = self.determinant()
            for i in range(self.cols):
                di = self.splice(ans, i).determinant()
                out.matrix[i][0] = di / D
            return out
        elif method == "jacobi":
            """jacobi method."""
            error = math.pow(10, -error)
            prev_x = Matrix(matrix=[[0] for i in range(self.rows)])

            r = self.strict_lower() + self.strict_upper()
            dinv = self.diagonal().inverse()
            x = dinv * (ans - r * prev_x)
            d = self.__distance(x, prev_x)
            i = 0

            while d > error:
                prev_x = x
                x = dinv * (ans - (r * prev_x))
                d = self.__distance(x, prev_x)
                i += 1
            return x, i

        elif method == "gauss":
            """gauss-seidol method."""
            error = math.pow(10, -error)
            prev_x = Matrix(matrix=[[0] for i in range(self.rows)])

            dlinv = (self.diagonal() + self.strict_lower()).inverse()
            x = dlinv * (ans - self.strict_upper() * prev_x)
            d = self.__distance(x, prev_x)
            i = 0

            while d > error:
                prev_x = x
                x = dlinv * (ans - self.strict_upper() * prev_x)
                d = self.__distance(x, prev_x)
                i += 1
            return x, i

    def transpose(self):
        """Swaps the rows and columns of a matrix."""
        out = Matrix(rows=self.cols, cols=self.rows)
        for i in range(out.rows):
            for j in range(out.cols):
                out.matrix[i][j] = self.matrix[j][i]
        return out

    def determinant(self):
        """Dynamic method for calculating the determinant of a matrix."""
        if self.size == (1, 1):
            determinant = self.matrix[0][0]
        else:
            determinant = 0
            for i in range(self.cols):
                partial_determinant = self.resize(0, i).determinant()
                determinant += self.matrix[0][i] * \
                    partial_determinant * math.pow(-1, i)
        return determinant

    def inverse(self):
        out = Matrix(rows=self.rows, cols=self.cols)

        for i in range(self.rows):
            for j in range(self.cols):
                partial_determinant = self.resize(i, j).determinant()
                out.matrix[i][j] = partial_determinant * math.pow(-1, i+j)
        return out.transpose() * (1/self.determinant())

    # f: lower, upper, diagonal
    def __triangle(self, f):
        if self.rows == self.cols:
            out = Matrix(rows=self.rows, cols=self.cols, mtype=0)
            for i in range(out.rows):
                for j in range(out.cols):
                    if f(i, j):
                        out.matrix[i][j] = self.matrix[i][j]
            return out

    def lower(self):
        return self.__triangle(lambda x, y: x >= y)

    def strict_lower(self):
        return self.__triangle(lambda x, y: x > y)

    def upper(self):
        return self.__triangle(lambda x, y: x <= y)

    def strict_upper(self):
        return self.__triangle(lambda x, y: x < y)

    def diagonal(self):
        return self.__triangle(lambda x, y: x == y)

    # f: add, sub
    def __addersubtractor(self, other, f):
        if self.size == other.size:
            out = Matrix(rows=self.rows, cols=self.cols)
            for i in range(out.rows):
                for j in range(out.cols):
                    out.matrix[i][j] = f(self.matrix[i][j], other.matrix[i][j])
            return out

    def __add__(self, other):
        return self.__addersubtractor(other, lambda x, y: x + y)

    def __sub__(self, other):
        return self.__addersubtractor(other, lambda x, y: x - y)

    def __mul__(self, other):
        if isinstance(other, int) or isinstance(other, float):
            out = Matrix(rows=self.rows, cols=self.cols)
            for i in range(out.rows):
                for j in range(out.cols):
                    out.matrix[i][j] = self.matrix[i][j] * other
            return out

        elif isinstance(other, Matrix) and self.cols == other.rows:
            out = Matrix(rows=self.rows, cols=other.cols, mtype=0)
            for i in range(self.rows):
                for j in range(other.cols):
                    for k in range(self.cols):
                        out.matrix[i][j] += self.matrix[i][k] * \
                            other.matrix[k][j]
            return out

    def __repr__(self):
        return '\n'.join([' '.join(map(str, row)) for row in self.matrix])


def main():
    mat_a = Matrix(
        [
            [10, 12, 15],
            [6, 8, 12],
            [12, 12, 18],
        ]
    )

    sol_a = Matrix(
        [
            [16*60],
            [11*60],
            [18*60],
        ]
    )

    print(mat_a.solve(sol_a, method="determinant"))


if __name__ == "__main__":
    main()
