using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace OperatResAndOptimMet
{
    class Vector
    {
        public double[] vector;
        public int size = 0;    //Размерность

        public Vector(int size)
        {
            this.size = size;
            vector = new double[size];

        }

        public Vector() : base()
        {
            this.size = 0;
            vector = new double[0];
        }

        public Vector(double[] m)
        {
            this.size = m.Length;
            vector = new double[size];
            for (int i = 0; i < size; i++)
                vector[i] = m[i];
        }

        public Vector(Vector m)
        {
            this.size = m.size;
            vector = new double[size];
            for (int i = 0; i < size; i++)
                vector[i] = m[i];
        }

        public double this[int i]
        {
            get
            {
                if (i < 0 && i >= size)
                {
                    Console.WriteLine("Индексы вышли за пределы матрицы");
                    return 0;
                }
                else
                    return vector[i];
            }
            set
            {
                if (i < 0 && i >= size)
                {
                    Console.WriteLine("Индексы вышли за пределы матрицы");
                }
                else
                    vector[i] = value;
            }
        }

        public int Size { get { return size; } }    //Получение размера

        //public int GetSize() {  return size;  }   

        public bool SetElement(double el, int index)    //Установить значение по индексу
        {
            if (index < 0 || index >= size) return false;
            vector[index] = el;
            return true;
        }

        public double GetElevent(int index)     //Получить значение по индексу
        {
            if (index < 0 || index >= size) return default(double); //С плавающей точкой
            return vector[index];
        }

        public Vector Copy()
        {
            Vector rez = new Vector(vector);
            return rez;
        }

        public void Clear()
        {
            for (int i = 0; i < this.size; i++) vector[i] = 0.0;
        }

        /*
        public override string ToString()
        {
            string s = "(";
            for (int i = 0; i < this.size; i++)
                s += "" + this[i] + " ";
            s += ")";
            return s;
        }
        */

        public override string ToString()
        => $"{{{string.Join(" ; ", this.vector)}}}";

        public void View()
        {
            Console.Write("( ");
            for (int i = 0; i < this.size; i++)
                Console.Write("{0} ", this[i]);
            Console.WriteLine(")");
        }

        public Vector Addition(Vector a)      //Сложение векторов
        {
            if (size == a.size)
            {
                Vector rez = new Vector(size);
                for (int i = 0; i < size; i++)
                    rez[i] = vector[i] + a[i];
                return rez;
            }
            return null;
        }

        public Vector Subtraction(Vector a)       //Вычитание векторов
        {
            if (size == a.size)
            {
                Vector rez = new Vector(size);
                for (int i = 0; i < size; i++)
                    rez[i] = vector[i] - a[i];
                return rez;
            }
            return null;
        }

        public double Multiplication(Vector a)   //Скалярное умножение векторов
        {
            if (size == a.size)                         //Формула: a*b = ∑(a(i)*b(i)) = a(x)*b(x) + a(y)*b(y)
            {
                double s = 0;
                for (int i = 0; i < size; i++)
                    s += vector[i] * a[i];
                return s;
            }
            return 0;
        }

        public Vector Multiplication_x(double x)       //Умножение вектора на число 
        {
            Vector rez = new Vector(size);
            for (int i = 0; i < size; i++)
                rez[i] = vector[i] * x;
            return rez;
        }

        public double Len()     //Определение длины вектора
        {                                                    //|a| = sqrt(a^2(x) + a^2(y))
            double x = 0;
            for (int i = 0; i < size; i++)
                x += Math.Pow(vector[i], 2);                 //Pow - Возвращает указанное число, возведенное в указанную степень
            x = Math.Sqrt(x);
            return x;
        }

        public Vector Normalization()     //Нормализация вектора
        {
            Vector rez = new Vector(vector);
            double x = Len();
            if (x > 0)                                          //Формула: a(норм) = a * 1/|a|, т.е. делим вектор на его длину;
            {
                for (int i = 0; i < size; i++)
                    rez[i] = rez[i] / x;
                return rez;
            }
            return null;
        }

        /* Второй вариант нормализации */
        public double NormaE()
        {
            return Math.Sqrt(this * this);
        }


        /* ************************************************************************************************************** */

        /* Сложение_векторов
         * Вычитание_векторов
         * Скалярное_проиведение_вектора_на_число
         * Скалярное_произведение_векоров
         * ПРИ ПОМОЩИ "operator"        
         */

        public static Vector operator +(Vector a, Vector b)
        {
            if (a.Size == b.Size)
            {
                Vector c = new Vector(a.Size);
                for (int i = 0; i < a.Size; i++)
                    c[i] += a[i] + b[i];
                return c;
            }
            return null;
        }

        public static Vector operator -(Vector a, Vector b)
        {
            if (a.Size == b.Size)
            {
                Vector c = new Vector(a.Size);
                for (int i = 0; i < a.Size; i++)
                    c[i] += a[i] - b[i];
                return c;
            }
            return null;
        }

        public static Vector operator *(Vector a, double c)
        {
            Vector r = new Vector(a.Size);
            for (int i = 0; i < a.Size; i++)
                r[i] = a[i] * c;
            return r;
        }

        public static Vector operator *(double c, Vector a)
        {
            Vector r = new Vector(a.Size);
            for (int i = 0; i < a.Size; i++)
                r[i] = a[i] * c;
            return r;
        }

        public static double operator *(Vector a, Vector b)   //Скалярное произведение 2-ух векторов
        {
            if (a.Size == b.Size)                                           //Формула: a*b = ∑(a(i)*b(i)) = a(x)*b(x) + a(y)*b(y)
            {
                double s = 0.0;
                for (int i = 0; i < a.Size; i++)
                    s += a[i] * b[i];
                return s;
            }
            return Double.NaN;
        }

        /* ************************************************************************************************************** */


        private Vector ForEach(Action2 action, Vector vector)
        {
            for (int index = 0; index < this.Size; ++index)
                action(ref this.vector[index], ref vector.vector[index]);

            return this;
        }

        public double ScalarProduct(Vector vector)     //Скалярное произведение 2-ух векторов при помощи делегатов (и лямбда-выражения); 
        {
            if (vector.Size != this.Size)
                throw new InvalidOperationException("Size of both vectors should be equals");

            var result = 0.0;
            this.ForEach((ref double f, ref double l) => result += f * l, vector);
            return result;
        }
        private delegate void Action2(ref double v, ref double l);
        private delegate void Action(ref double v);

        //Случайный вектор единичной длины (нужен для метода случайного поиска)
        public Vector RandMSP(int n, Random rnd)
        {
            double[] x = new double[n];
            for (int i = 0; i < n; i++)
            {
                x[i] = rnd.NextDouble() - 0.5;
            }
            Vector psi = new Vector(x);
            double norm_psi = psi.NormaE();
            for (int i = 0; i < n; i++)
            {
                psi[i] = psi[i] / norm_psi;
            }
            return psi;

        }

    }
}
