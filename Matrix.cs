using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace OperatResAndOptimMet
{
    class Matrix
    {
        protected int rows, columns;     //Переменные для строк и столбцов матрицы
        protected double[,] data;        //Двумерный массив, т.е. массив с 2-мя измерениями (строки и столбцы); [,,] - трехмерный массив...

        public Matrix(int r, int c)
        {
            this.rows = r; this.columns = c;
            data = new double[rows, columns];
            for (int i = 0; i < rows; i++)
                for (int j = 0; j < columns; j++) data[i, j] = 0;
        }

        public Matrix(double[,] mm)
        {
            this.rows = mm.GetLength(0); this.columns = mm.GetLength(1);
            data = new double[rows, columns];
            for (int i = 0; i < rows; i++)
                for (int j = 0; j < columns; j++)
                    data[i, j] = mm[i, j];
        }

        public int GetCountRows()
        {
            return rows;
        }

        public int GetCountColumns()
        {
            return columns;
        }

        public double this[int i, int j]
        {
            get
            {
                if (i < 0 && j < 0 && i >= rows && j >= columns)
                {
                    // Console.WriteLine(" Индексы вышли за пределы матрицы ");
                    return Double.NaN;
                }
                else
                    return data[i, j];
            }
            set
            {
                if (i < 0 && j < 0 && i >= rows && j >= columns)
                {
                    //Console.WriteLine(" Индексы вышли за пределы матрицы ");
                }
                else
                    data[i, j] = value;
            }
        }

        public Vector GetRow(int r)  //Получить строки
        {
            if (r >= 0 && r < rows)
            {
                Vector row = new Vector(columns);
                for (int j = 0; j < columns; j++) row[j] = data[r, j];
                return row;
            }
            return null;
        }

        public void SetRow(int r, Vector rr) //Набор строк
        {
            if (r >= 0 && r < rows)
            {
                if (columns == rr.Size)
                    for (int j = 0; j < columns; j++) data[r, j] = rr[j];
            }
        }

        public Vector GetColumn(int c)   //Получить столбцы
        {
            if (c >= 0 && c < columns)
            {
                Vector column = new Vector(rows);
                for (int i = 0; i < rows; i++) column[i] = data[i, c];
                return column;
            }
            return null;
        }

        public void SetColumn(int c, Vector cc)  //Набор столбцов
        {
            if (c >= 0 && c < columns)
            {
                if (rows == cc.Size)
                    for (int i = 0; i < rows; i++) data[i, c] = cc[i];
            }
        }

        //Корирование матрицы
        public Matrix Copy()
        {
            Matrix matrix_copy = new Matrix(columns, rows);
            for (int i = 0; i < rows; i++)
                for (int j = 0; j < columns; j++)
                    matrix_copy.data[i, j] = data[i, j];
            return matrix_copy;

        }

        public void Print() //Mетод, который печатает матрицу
        {
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < columns; j++)
                {
                    Console.Write($"{data[i, j]} \t");             //t - горизонтальная табуляция
                }
                Console.WriteLine();
            }
            Console.WriteLine("\n");
        }

        public override string ToString()
        {
            string s = "{";
            for (int i = 0; i < rows - 1; i++)
            {
                Vector v = this.GetRow(i);
                s += v.ToString() + ";";
            }
            Vector vc = this.GetRow(rows - 1);
            s += vc.ToString() + "}";
            return s;
        }

        public static Vector operator *(Matrix a, Vector b)    //умножаем матрицу на вектор
        {
            if (a.columns != b.Size) return null;
            Vector r = new Vector(a.rows);
            for (int i = 0; i < a.rows; i++)
            {
                r[i] = a.GetRow(i) * b;
            }
            return r;
        }

        public Matrix Transp() //Метод, который транспонирует матрицу
        {
            Matrix tr = new Matrix(columns, rows);
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < columns; j++)
                {
                    tr.data[j, i] = data[i, j];
                    //Console.Write(data[j, i] + "\t"); 
                }
            }
            return tr;
        }

        public static Matrix Umnch(Matrix a, double ch)   //Умножаем матрицу на число
        {
            Matrix resMass = new Matrix(a.rows, a.columns);
            for (int i = 0; i < a.rows; i++)
            {
                for (int j = 0; j < a.columns; j++)
                {
                    resMass[i, j] = a[i, j] * ch;
                }
            }
            return resMass;
        }

        public static Matrix UmnMatrix(Matrix a, Matrix b)     //Умножение матриц
        {
            /*if (a.GetCountColumns() != b.GetCountRows())
            {
                throw new Exception("Умножение не возможно! Количество столбцов первой матрицы не равно количеству строк второй матрицы.");
            }*/

            Matrix rezMatrix = new Matrix(a.rows, b.columns);

            try
            {
                if (a.GetCountColumns() != b.GetCountRows()) { throw new OverflowException(); }
                else
                {
                    for (int i = 0; i < a.rows; i++)
                    {
                        for (int j = 0; j < b.columns; j++)
                        {
                            rezMatrix[i, j] = 0;
                            for (int k = 0; k < a.columns; k++)
                            {
                                rezMatrix[i, j] += a[i, k] * b[k, j];
                            }
                        }
                    }
                }
            }
            catch (OverflowException)
            {
                Console.WriteLine("Умножение невозможно! Количество столбцов первой матрицы не равно количеству строк второй матрицы!");
            }
            return rezMatrix;
        }

        public static Matrix SumMatrix(Matrix a, Matrix b)     //Сложение матриц
        {
            Matrix rezMatrix = new Matrix(a.rows, a.columns);

            try
            {                                                                                               //Генерируем ошибку
                if (a.GetCountColumns() != b.GetCountColumns() || a.GetCountRows() != b.GetCountRows()) { throw new OverflowException(); }
                else
                {
                    for (int i = 0; i < a.rows; i++)
                    {
                        for (int j = 0; j < b.columns; j++)
                        {
                            rezMatrix[i, j] += a[i, j] + b[i, j];
                        }
                    }
                    //return rezMatrix;
                }
            }
            catch (OverflowException)
            {
                Console.WriteLine("Сложение невозможно, размеры матриц должны быть одинаковыми!");
            }

            return rezMatrix;
        }

        public static Matrix SubtractionMatrix(Matrix a, Matrix b)     //Вычитание матриц
        {
            Matrix rezMatrix = new Matrix(a.rows, a.columns);

            try
            {                                                                                               //Генерируем ошибку
                if (a.GetCountColumns() != b.GetCountColumns() || a.GetCountRows() != b.GetCountRows()) { throw new OverflowException(); }
                else
                {
                    for (int i = 0; i < a.rows; i++)
                    {
                        for (int j = 0; j < b.columns; j++)
                        {
                            rezMatrix[i, j] += a[i, j] - b[i, j];
                        }
                    }
                }
            }
            catch (OverflowException)
            {
                Console.WriteLine("Вычитание невозможно, размеры матриц должны быть одинаковыми!");
            }

            return rezMatrix;
        }

        //СЛУ с нижней треугольной матрицей
        public static Vector SLU_DOWN(Matrix a, Vector b)
        {
            //Контроль
            if (a.rows != b.Size) return null;
            for (int i = 0; i < a.rows; i++)
                if (a[i, i] == 0) return null;

            //for (int j = i + 1; j < a.columns; j++) if (a[i, j] != 0.0) return null;
            //    for (int i = 0; i < a.columns; i++) if (a[i, i] == 0.0) return null;

            Vector x = new Vector(b.Size);
            //Forward
            double c;
            x[0] = b[0] / a[0, 0];
            for (int i = 1; i < a.rows; i++)
            {
                c = 0;
                for (int j = 0; j < i; j++) { c += a[i, j] * x[j]; }
                x[i] = (b[i] - c) / (a[i, i]);
            }
            return x;

        }

        //СЛУ с верхней треугольной матрицей
        public static Vector SLU_UP(Matrix a, Vector b)
        {
            //Контроль
            if (a.rows != b.Size) return null;
            for (int i = 0; i < a.rows; i++)
                if (a[i, i] == 0) return null;

            //for (int j = i - 1; j < a.columns; j++) if (a[i, j] != 0.0) return null;// ошибка
            //    for (int i = 0; i < a.columns; i++) if (a[i, i] == 0.0) return null;

            Vector x = new Vector(b.Size);
            x[b.Size - 1] = b[b.Size - 1] / a[a.rows - 1, a.columns - 1];

            for (int i = a.rows - 2; i >= 0; i--)
            {
                double c = 0;
                for (int j = b.Size - 1; j > i; j--) { c += a[i, j] * x[j]; }
                x[i] = (b[i] - c) / (a[i, i]);
            }
            return x;

        }

        //Метод Гаусса
        public static Vector Gauss(Matrix aa, Vector bb)
        {
            //Цель: привести матрицу к верхней треугольной матрице, обнуляя элементы под главной диагональю, а далее решить СЛУ
            Matrix a = aa.Copy();
            Vector b = bb.Copy();
            Vector x;

            if (a.columns != a.rows) return null;
            if (a.columns != b.Size) return null;

            for (int j = 0; j < a.rows - 1; j++)
            {
                //Ищем макс элемент в j-ом столбце
                double maxc = 0;
                int jmax = j;
                for (int i = j; i < a.rows; i++)
                {
                    if (Math.Abs(a[i, j]) > maxc) { maxc = Math.Abs(a[i, j]); jmax = i; }
                }

                if (maxc == 0) return null;

                //Меняем строки j и jmax
                Vector rowj = a.GetRow(j);
                Vector rowjmax = a.GetRow(jmax);
                //Меняем строки j и jmax
                for (int k = 0; k < a.columns; k++)
                {
                    a[j, k] = rowjmax[k];
                    a[jmax, k] = rowj[k];
                }
                double c = b[j]; b[j] = b[jmax]; b[jmax] = c;

                for (int i = j + 1; i < a.rows; i++)
                {
                    double aij = a[i, j] / a[j, j];
                    //умножаем строку и добавляем строку //b

                    //обнуляем ниже диагонали
                    for (int k = 0; k < a.columns; k++)
                    {
                        a[i, k] = a[i, k] - aij * a[j, k];

                    }
                    b[i] = b[i] - aij * b[j];
                }
            }
            x = SLU_UP(a, b);

            return x;
        }

        //Метод Гивенса
        public static Vector MethodGivens(Matrix aa, Vector bb)
        {
            // Метод вращений ( Гивенса ) основан на преобразовании исходной матрицы А путем
            //последовательного умножения на матрицы вращений к треугольной матрице

            Matrix a = aa.Copy();
            //Matrix_NumMet givens = new Matrix_NumMet(a.columns, a.rows); //Суммарная или общая матрца вращений
            Vector b = bb.Copy();
            Vector x;

            if (a.columns != a.rows) return null;   //Проверка! матрица д/б квадратной
            if (a.columns != b.Size) return null;   //Провекра вектора на размерность

            double[,] Ematrix = new double[a.rows, a.columns];
            for (int k = 0; k < a.rows; k++)
                for (int l = 0; l < a.columns; l++)
                    if (k == l)
                        Ematrix[k, l] = 1;
                    else
                        Ematrix[k, l] = 0;

            for (int i = 0; i < a.rows - 1; i++)
                for (int j = i + 1; j < a.rows; j++)
                {
                    Matrix Ed = new Matrix(Ematrix);

                    double cos = a[i, i] / Math.Sqrt(a[i, i] * a[i, i] + a[j, i] * a[j, i]);    // cos угла поворота
                    double sin = -a[j, i] / Math.Sqrt(a[i, i] * a[i, i] + a[j, i] * a[j, i]);   // sin угла поворота

                    Ed[i, i] = cos;
                    Ed[i, j] = sin;
                    Ed[j, i] = -sin;
                    Ed[j, j] = cos;
                    a = Ed.Transp() * a;
                    b = Ed.Transp() * b;
                }

            x = SLU_UP(a, b);
            return x;
        }

        //Метод прогонки (полностью не работает)
        //вариант 1 через векторы
        public static Vector Progonka_PolnostNeRabotaet(Vector d, Vector m, Vector u, Vector xx)
        {
            Vector a = d.Copy();
            Vector b = m.Copy();
            Vector c = u.Copy();
            Vector x = xx.Copy();
            Vector otv = new Vector(b.size);
            Vector alpha = new Vector(b.size);
            Vector beta = new Vector(b.size);
            Vector y = new Vector(b.size);

            if (a.size + 1 != b.size || b.size != c.size + 1 || a.size != c.size) return null;

            y[0] = b[0];
            alpha[0] = -c[0] / y[0];
            beta[0] = x[0] / y[0];

            for (int i = 1; i < b.size - 1; i++)
            {

                y[i] = b[i] + ((a[i - 1]) * alpha[i - 1]);
                alpha[i] = -c[i] / y[i];
                beta[i] = (x[i] - a[i - 1] * beta[i - 1]) / y[i];
            }

            y[b.size - 1] = b[b.size - 1] + a[b.size - 2] * alpha[b.size - 2];
            beta[b.size - 1] = (x[b.size - 1] - a[b.size - 2] * beta[b.size - 2]) / y[b.size - 1];

            otv[b.size - 1] = beta[b.size - 1];

            for (int i = b.size - 2; i > -1; i--)
            {
                otv[i] = alpha[i] * otv[i + 1] + beta[i];
            }

            return otv;
        }

        //Метод прогонки (работает)
        public static Vector ProgonkaRabotaet(Vector c, Vector d, Vector e, Vector b)
        {
            //Дано: трехдиагональная матрица
            int n = d.Size;
            if (b.Size != n) return null;

            Vector x = new Vector(n);
            Vector alpha = new Vector(n);
            Vector betta = new Vector(n);

            for (int i = 0; i < n; i++)
            {
                if (d[i] == 0) return null;
            }

            //Прямой ход
            alpha[1] = -e[0] / d[0];
            betta[1] = b[0] / d[0];
            for (int i = 1; i < n - 1; i++)
            {
                double zn = d[i] + c[i] * alpha[i];
                alpha[i + 1] = -e[i] / zn;
                betta[i + 1] = (-c[i] * betta[i] + b[i]) / zn;

            }

            //Обратный ход
            x[n - 1] = (-c[n - 1] * betta[n - 1] + b[n - 1]) / (d[n - 1] + c[n - 1] * alpha[n - 1]);
            for (int i = n - 2; i >= 0; i--)
            {
                x[i] = alpha[i + 1] * x[i + 1] + betta[i + 1];
            }
            return x;

        }

        //Метод последовательных приближений
        public static Vector PosledPr(Matrix aa, Vector bb, double eps, int kmax)
        {
            Matrix a = aa.Copy();

            Vector b = bb.Copy();
            Vector x = bb.Copy();
            Vector xs = bb.Copy();
            Vector r = bb.Copy();

            if (a.columns != b.Size || a.rows != a.columns) return null;

            for (int i = 0; i < a.rows; i++)
            {
                if (a[i, i] != 0)
                {
                    for (int j = 0; j < a.rows; j++)
                    {
                        if (i != j)
                            //if (a[i, i] != 0.0)
                            a[i, j] = -a[i, j] / a[i, i];
                        //else Console.WriteLine("Ошибка! Деление на ноль!");

                    }
                    b[i] = b[i] / a[i, i];
                }
                a[i, i] = 0;
            }

            x = b;
            int k = 0;
            do
            {
                k++;
                xs = a * x + b;
                r = xs - x;
                x = xs.Copy();
                if (k > kmax) return null;
            }
            while (r.NormaE() > eps);

            return x;
        }

        public void View()
        {
            for (int i = 0; i < this.rows; i++)
            {
                for (int j = 0; j < this.columns; j++)
                    Console.Write("{0} ", this[i, j]);
                Console.WriteLine(" ");
            }
        }

        public static Matrix operator *(Matrix a, Matrix b)
        {
            return Matrix.UmnMatrix(a, b);
        }

        public static Matrix operator *(Matrix a, double b)
        {
            return Matrix.Umnch(a, b);
        }

        public static Matrix operator *(double b, Matrix a)
        {
            return Matrix.Umnch(a, b);
        }

        public static Matrix operator +(Matrix a, Matrix b)
        {
            return Matrix.SumMatrix(a, b);
        }

        public static Matrix operator -(Matrix a, Matrix b)
        {
            return Matrix.SubtractionMatrix(a, b);
        }

    }
}
