using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace OperatResAndOptimMet
{
    delegate double func(double x);  //Для одномерной функции
    delegate double funcVec(Vector x);
    delegate double[] funcVecMas(double[] x);
    delegate Vector funcVecKriteriy(Vector x);  //Вектор частных критериев

    class Optimization
    {
        //Численный метод поиска локального экстремума (Пошаговый метод)
        public static double StepByStep(double xn, double h, double eps, func f)
        {
            int k = 2;
            double fn = f(xn);
            double xs = xn + h;
            double fs = fn + f(xs);

            while (Math.Abs(h) > eps)
            {
                if (fs > fn)
                {
                    h = -h / 2;
                }
                else
                {
                    h = h * 1.2;
                }
                xn = xs;
                fn = fs;
                xs = xn + h;
                fs = f(xs);
                k++; 
            }
            //Console.WriteLine("k = {0} - кол-во итераций в пошаговом методе", k);
            return xs;

        }

        //Метод золотого сечения
        public static double GoldenRatio(double a, double b, double eps, func f)
        {
            int k = 2;
            double l = b - a;
            double fi = (Math.Sqrt(5.0) - 1) / 2.0;
            double v = a + l * (1.0 - fi);
            double w = a + l * fi;
            double fv = f(v);
            double fw = f(w);

            while (l > eps)
            {
                if (fv < fw)
                {
                    b = w;
                    w = v;
                    fw = fv;
                    l = b - a;
                    v = a + l * (1.0 - fi);
                    fv = f(v);

                }
                else
                {
                    a = v;
                    v = w;
                    fv = fw;
                    l = b - a;
                    w = a + l * fi; ;
                    fw = f(w);
                }
                k++;
            }
            Console.WriteLine("\nk = {0} - кол-во итераций в методе 'Золотое сечение'", k);
            double x = (a + b) / 2.0;
            return x;
            
        }

        //Метод квадратичной аппроскимации
        public static double MethodKvApproksim(double x1, double x2, double x3, double eps, func fun)
        {
            double xpred = DopMethod(x1, x2, x3, fun);
            double xtek = 0;
            
            while (Math.Abs(xpred - xtek) < eps)
            {
                if ((x1 > x2) && (x1 > x3))
                {
                    x1 = xpred;
                    xtek = DopMethod(x1, x2, x3, fun);
                }

                if ((x2 > x1) && (x2 > x3))
                {
                    x2 = xpred;
                    xtek = DopMethod(x1, x2, x3, fun);
                }

                if ((x3 > x1) && (x3 > x2))
                {
                    x3 = xpred;
                    xtek = DopMethod(x1, x2, x3, fun);
                }
            }

            return xpred;
        }

        //Для квадратичной аппроксимации
        public static double DopMethod(double x1, double x2, double x3, func fun)
        {
            double f1 = fun(x1);
            double f2 = fun(x2);
            double f3 = fun(x3);
            double c = ((f3 - f1) * (x2 - x1) - (f2 - f1) * (x3 - x1)) / ((x3 * x3 - x1 * x1) * (x2 - x1) - (x2 * x2 - x1 * x1) * (x3 - x1));
            double b = ((f2 - f1) - c * (x2 * x2 - x1 * x1)) / (x2 - x1);
            double a = (f1 - b * x1 - c * x1 * x1);
            double xpred = -b / 2 * c;
            double fm = xpred * xpred * c + b * xpred + a;
            return xpred;

        }

        //Градиентный метод поиска минимума
        public static Vector Grad(Vector xn, double h, double eps, funcVec f)
        {
            Vector xnach = new Vector(xn);
            int k = 0;
            int n = xn.size;
            double fn = f(xn);
            Vector xs = new Vector(xn);
            Vector gr;  //Градиент
            Vector dx;  //Изменение
            double delta = 0.5 * eps;
            //double fs = f(xs);
            do
            {
                gr = new Vector(n);
                for (int i = 0; i < n; i++)
                {
                    Vector xg = xn.Copy();
                    xg[i] = xg[i] + delta;
                    gr[i] = (f(xg) - fn) / delta;
                }
                dx = gr * -h;
                xs = xn + dx;  //След точка
                double fs = f(xs);

                if (fs > fn)
                {
                    h = h / 2;

                }
                else
                {
                    h = 1.2 * h;
                }
                xn = xs;
                fn = fs;
                k++;

            } while (Math.Abs(dx.NormaE()) > eps);

            Console.WriteLine("\nk = {0} - кол-во итераций в методе градиента", k);
            return xs;

        }

        // Метод наискорейшего градиентного спуска (МНГС)
        public static Vector MNGS(Vector xn, double eps, funcVec f)
        {
            Vector xnach = new Vector(xn);
            int k = 0;  //Переменная для подсчета итераций
            int n = xn.size;  //Размерность
            Vector xs = new Vector(xn);  //Следующее значение вектора
            Vector grad = new Vector(n);  //Вектор градиента
            Vector dx;  //Вектор производной (изменение)
            double delta = 0.5 * eps;  //Достаточно малое число
            double fn = f(xn);
            //double fs = fn + f(xs);
            double fs = f(xs);
            double h_opt = 0.1;

            do
            {
                //Расчитываем градиент - начало
                for (int i = 0; i < n; i++)
                {
                    Vector xg = xn.Copy();
                    xg[i] = xg[i] + delta;
                    grad[i] = (f(xg) - fn) / delta;
                }
                //Расчитываем градиент - конец

                dx = -h_opt * grad;  //Изменение точки
                h_opt = StepByStep(h_opt, 0.5, eps, h_tmp => f(xn - h_tmp * grad));
                xs = xn + dx;  //X^(k+1) = X^(k) - h_opt * S^(k)  --- именно отсюда берется ответ
                fs = f(xs);  //То, что мы получили в xs (т.е. мы нашли x1 и x2) подставили в исходную функцию
                
                xn = xs;  //Для дальнейших итераций/вычислений
                fn = fs;
                k++;

            } while (Math.Abs(dx.NormaE()) > eps);

            Console.WriteLine("\nk = {0} - кол-во итераций в методе наискорейшего градиентного спуска", k);
            return xs;
        }

        //Метод сопряженных направлений (МСН)  (Проверить)
        public static Vector MSN(Vector xn, double eps, funcVec f)
        {
            Vector xnach = new Vector(xn);
            int k = 0;  //Переменная для подсчета итераций
            int n = xn.size;  //Размерность
            Vector xs = new Vector(xn);  //Следующее значение вектора
            Vector grad = new Vector(n);  //Вектор градиента
            Vector dx;  //Вектор производной
            double delta = 0.5 * eps;  //Достаточно малое число
            double fn = f(xn);
            double fs = fn + f(xs);
            double h_opt = 0.1;

            Vector grad_old = new Vector(n);
            double gamma;

            //Высчитываем градиент - начало
            for (int i = 0; i < n; i++)
            {
                Vector xg = xn.Copy();
                xg[i] = xg[i] + delta;
                grad[i] = (f(xg) - fn) / delta;
            }
            //Высчитываем градиент - конец

            dx = -h_opt * grad;
            h_opt = StepByStep(h_opt, 0.5, eps, h_tmp => f(xn - h_tmp * grad));
            xs = xn + dx;  //X^(k+1) = X^(k) - h_opt * S^(k)  --- именно отсюда берется ответ
            fs = f(xs);  //То, что мы получили в xs (т.е. мы нашли x1 и x2) подставили в исходную функцию
            
            grad_old = grad.Copy();
            xn = xs;  //Для дальнейших итераций/вычислений
            fn = fs;
            k++;

            do
            {
                //Расчитываем градиент - начало
                for (int i = 0; i < n; i++)
                {
                    Vector xg = xn.Copy();
                    xg[i] = xg[i] + delta;
                    grad[i] = (f(xg) - fn) / delta;
                }
                //Расчитываем градиент - конец

                h_opt = StepByStep(h_opt, 0.5, eps, h_tmp => f(xn - h_tmp * grad));

                gamma = (grad.NormaE() * grad.NormaE()) / (grad_old.NormaE() * grad_old.NormaE());
                grad = grad + gamma * grad_old;

                dx = -h_opt * grad;  //?
                xs = xn + dx;  //X^(k+1) = X^(k) - h_opt * S^(k)  --- именно отсюда берется ответ
                fs = f(xs);  //То, что мы получили в xs (т.е. мы нашли x1 и x2) подставили в исходную функцию

                grad_old = grad.Copy();
                xn = xs;  //Для дальнейших итераций/вычислений
                fn = fs;
                k++;

            } while (Math.Abs(dx.NormaE()) > eps);

            Console.WriteLine("\nk = {0} - кол-во итераций в методе сопряженных направлений", k);
            return xs;
        }
                        
        //Метод наилучшей пробы (МСП)
        public static Vector MSP_best_sampling(Vector xn, int n, double h, double eps, funcVec f)
        {
            Random rnd = new Random();             
            int ni = xn.Size;
            Vector xs = new Vector(ni);
            int M = 3 * n;  //Число испытаний на текущей итерации;  n - размерность функции
            Vector[] x_prob = new Vector[M];
            int k = 0;
                       
            do
            {
                //xs = xn - xn.RandMSP(2, rnd);
                double fn = f(xn);
                double fs = double.MaxValue;

                for (int i = 0; i < M; i++)
                {
                    //Формируем вокруг xn набор случайных точек:
                    Vector psi = xn.RandMSP(2, rnd); //Вектор единичной длины [пси - Ψ]
                    Vector Xj_prob = xn - h * psi;  //Следующая точка
                    x_prob[i] = Xj_prob;
                    double fj = f(x_prob[i]);
                    
                    if (fs > fj)
                    {
                        fs = fj;
                        xs = x_prob[i];
                    }
                }

                if (fs < fn)  //Шаг удачный
                {
                    xn = xs;
                    h = h * 1.2;
                }

                if (fs >= fn)  //Шаг неудачный
                {
                    h = h * 0.5;
                }
                k++;

            } while (h > eps);

            Console.WriteLine("\nk = {0} - кол-во итераций в методе случайного поиска", k);
            return xn;
        }

        
        //Метод деформируемого многогранника (м. Нелдера-Мида)
        public static Vector Nelder_Mead(Vector[] xn, double eps, funcVec func)
        {
            int n = xn.Length;
            double m = n + 1;


            Vector xotr = new Vector(xn[0].Size);
            double cps = 1.0 / (n - 1.0);

            double k = 0;
            double eeeps = 10;
            double[] fmas = new double[n];

            do
            {
                Vector xc = new Vector(xn[0].Size);
                double kof = 1;
                for (int i = 0; i < n; i++)
                {
                    fmas[i] = func(xn[i]);
                }
                fmas = ShellSort(fmas);
                for (int i = 0; i < n; i++)
                {
                    if (func(xn[i]) == fmas[n - 1]) { }
                    else
                    {
                        xc += xn[i] * cps;
                    }
                }
                xotr = xc + kof * (xc - FindXN(xn, fmas[n - 1], func));
                double fotr = func(xotr);

                if (fotr < fmas[0])
                {
                    double newfotr = fotr;

                    eeeps = Math.Abs(func(FindXN(xn, fmas[n - 1], func)) - func(xotr));
                    if (fotr < fmas[1])
                    {
                        kof = kof * 2;
                        xotr = xc + kof * (xotr - xc);
                    }
                    for (int i = 0; i < n; i++)
                    {
                        if (xn[i] == FindXN(xn, fmas[n - 1], func))
                        {
                            xn[i] = xotr.Copy();
                        }
                    }
                }
                else if (fotr > fmas[0])
                {
                    kof = kof / 2;
                    xotr = xc + kof * (FindXN(xn, fmas[n - 1], func) - xc);
                    eeeps = Math.Abs(-func(FindXN(xn, fmas[n - 1], func)) + func(xotr));
                    for (int i = 0; i < n; i++)
                    {
                        if (xn[i] == FindXN(xn, fmas[n - 1], func))
                        {
                            xn[i] = xotr.Copy();
                        }
                    }

                }
                else
                {
                    for (int i = 0; i < n; i++)
                    {
                        if (xn[i] != FindXN(xn, fmas[0], func))
                        {
                            for (int j = 0; j < xn[i].Size; j++)
                            {
                                xn[i][j] = xn[i][j] / 2;
                            }
                        }
                    }
                }
                k++;

            } while (eps < eeeps);  // https://habr.com/ru/post/332092/
            Console.WriteLine("\nk = {0} - кол-во итераций в методе деформируемого многогранника:", k);
            return FindXN(xn, fmas[0], func);

        }

        private static Vector FindXN(Vector[] xn, double res, funcVec func)
        {
            Vector n = new Vector(xn[0].Size);
            for (int i = 0; i < xn.Length; i++)
            {
                if (func(xn[i]) == res)
                {
                    n = xn[i];
                }


            }
            if (n == new Vector(xn[0].Size)) throw new Exception();
            return n;
        }

        private static double[] ShellSort(double[] list)  //Сортировка Шелла (для метода деформир многогранника)
        {
            //расстояние между элементами, которые сравниваются
            var d = list.Length / 2;
            while (d >= 1)
            {
                for (var i = d; i < list.Length; i++)
                {
                    var j = i;
                    while ((j >= d) && (list[j - d].CompareTo(list[j]) > 0))
                    {
                        Swap(ref list[j], ref list[j - d]);
                        j = j - d;
                    }
                }

                d = d / 2;
            }
            return list;
        }

        private static void Swap(ref double a, ref double b)
        {
            double c = a; a = b; b = c;
        }


        static Random rnd = new Random();

        //Метод ОЗУ (с векторами)
        public static Vector OZU(Vector xn, int n, double h, double eps, Vector Fv, Vector Fn, int[] tip_ogr, funcVecKriteriy ff)
        {
            int k = 0;
            int m = xn.Size;
            int pmin = int.MinValue;
            Matrix xpTemp = new Matrix(m, n);

            Vector xp = new Vector(n);
            Vector fp = new Vector(m);
            double ft = Limitation(xn, Fn, Fv, tip_ogr, ff);
            double fpmin = double.MaxValue;
            double temp, len;

            do
            {
                for (int i = 0; i < m; i++)
                {
                    for (int j = 0; j < n; j++)
                    {
                        temp = rnd.NextDouble() - 0.5;
                        xpTemp[i, j] = temp;
                    }

                    //Нормализация массива
                    len = 0;
                    for (int j = 0; j < n; j++)
                    {
                        len += (xpTemp[i, j] - xn[j]) * (xpTemp[i, j] - xn[j]);
                    }

                    len = Math.Sqrt(len);
                    for (int j = 0; j < n; j++)
                    {
                        xpTemp[i, j] = xn[j] + h * (xpTemp[i, j] - xn[i]) / len;
                    }

                }

                for (int i = 0; i < m; i++)
                {
                    for (int j = 0; j < n; j++)
                    {
                        xp[j] = xpTemp[i, j];
                    }
                    fp[i] = Limitation(xn, Fn, Fv, tip_ogr, ff);
                    if (fp[i] < fpmin)
                    {
                        fpmin = fp[i];
                        pmin = i;
                    }
                }

                if (fpmin < ft)
                {
                    for (int j = 0; j < n; j++)
                    {
                        xn[j] = xpTemp[pmin, j];
                    }
                    ft = fpmin;
                    h = h * 1.2;
                }
                else
                    h = h / 2.0;

                k++;

            } while (h > eps);

            Console.WriteLine("\nk = {0} - кол-во итераций в методе ОЗУ", k);
            return xn;

        }

        //Для ОЗУ (Расчет ограничений)
        public static double Limitation(Vector x, Vector jv, Vector jn, int[] tip_ogr, funcVecKriteriy ff)
        {
            Vector fc = ff(x);
            Vector phi = new Vector(fc.Size);
            int max_phi = 0;
            double phi_1;  //Для двухстороннего ограничения
            double phi_2;  //Для двухстороннего ограничения

            for (int i = 0; i < fc.Size; i++)
            {
                if (tip_ogr[i] == 0) //Ограничение снизу
                {
                    if (jn[i] > 0)
                    {
                        phi[i] = 2 - fc[i] / jn[i];
                    }

                    if (jn[i] == 0)
                    {
                        phi[i] = 1 - fc[i] / jn[i];
                    }

                    if (jn[i] < 0)
                    {
                        phi[i] = fc[i] / jn[i];
                    }

                }

                if (tip_ogr[i] == 1) //Ограничение сверху
                {
                    if (jv[i] > 0)
                    {
                        phi[i] = fc[i] / jn[i];
                    }

                    if (jv[i] == 0)
                    {
                        phi[i] = 1 - fc[i] / jv[i];
                    }

                    if (jv[i] < 0)
                    {
                        phi[i] = 2 - fc[i] / jv[i];
                    }

                }

                if (tip_ogr[i] == 2) //Двухстороннее ограничение 
                {
                    phi_1 = (fc[i] - jn[i]) / (jv[i] - jn[i]);
                    phi_2 = (jv[i] - fc[i]) / (jv[i] - jn[i]);

                    if (phi_1 > phi_2)
                    {
                        phi[i] = phi_1;
                    }
                    else
                        phi[i] = phi_2;

                }

                if (phi[i] > phi[max_phi])
                {
                    max_phi = i;
                }

            }
            return phi[max_phi];

        }

        


        //Метод ОЗУ (С помощью массивов)
        public double[] OZU_mas(double[] x, double eps, double h, double[] Fn, double[] Fv, int[] tipf, funcVecMas f)
        {
            int k = 1;
            int n = x.Length;
            int m = 3 * x.Length;
            int pmin = int.MinValue;
            double[,] xpTemp = new double[m, n];
            double[] xp = new double[m];
            double[] fp = new double[m];
            double ft = Limitation_mas(x, Fn, Fv, tipf, f);
            double fpmin = double.MaxValue;
            double temp, len;

            do
            {
                for (int i = 0; i < m; i++)
                {
                    for (int j = 0; j < n; j++)
                    {
                        temp = rnd.NextDouble() - 0.5;
                        xpTemp[i, j] = temp;
                    }

                    //нормализация массива
                    len = 0;
                    for (int j = 0; j < n; j++)
                    {
                        len += (xpTemp[i, j] - x[j]) * (xpTemp[i, j] - x[j]);
                    }

                    len = Math.Sqrt(len);
                    for (int j = 0; j < n; j++)
                    {
                        xpTemp[i, j] = x[j] + h * (xpTemp[i, j] - x[i]) / len;
                    }

                }

                for (int i = n; i < m; i++)
                {
                    for (int j = 0; j < n; j++)
                    {
                        xp[j] = xpTemp[i, j];
                    }
                    fp[i] = Limitation_mas(xp, Fn, Fv, tipf, f);
                    if (fp[i] < fpmin)
                    {
                        fpmin = fp[i];
                        pmin = i;
                    }
                }

                if (fpmin < ft)
                {
                    for (int j = 0; j < n; j++)
                    {
                        x[j] = xpTemp[pmin, j];
                    }
                    ft = fpmin;
                    h = h * 1.2;
                }
                else
                    h = h / 2.0;

                k++;

            } while (h > eps);
            return x;

        }

        //Для ОЗУ
        public static double Limitation_mas(double[] x, double[] jn, double[] jv, int[] tipf, funcVecMas f)
        {
            double[] jx = f(x);
            int sizeCriteria = jx.Length;
            double[] phi = new double[sizeCriteria];
            int max_phi = 0;
            double g1;
            double g2;

            for (int i = 0; i < sizeCriteria; i++)
            {
                if (tipf[i] == 0) //Ограничение снизу
                {
                    if (jn[i] > 0)
                    {
                        phi[i] = 2 - jx[i] / jn[i];
                    }

                    if (jn[i] == 0)
                    {
                        phi[i] = 1 - jx[i] / jn[i];
                    }

                    if (jn[i] < 0)
                    {
                        phi[i] = jx[i] / jn[i];
                    }

                }

                if (tipf[i] == 1) //Ограничение сверху
                {
                    if (jv[i] > 0)
                    {
                        phi[i] = jx[i] / jn[i];
                    }

                    if (jv[i] == 0)
                    {
                        phi[i] = 1 - jx[i] / jv[i];
                    }

                    if (jv[i] < 0)
                    {
                        phi[i] = 2 - jx[i] / jv[i];
                    }

                }

                if (tipf[i] == 2) //Двухстороннее ограничение 
                {
                    g1 = (jx[i] - jn[i]) / (jv[i] - jn[i]);
                    g2 = (jv[i] - jx[i]) / (jv[i] - jn[i]);

                    if (g1 > g2)
                    {
                        phi[i] = g1;
                    }
                    else
                        phi[i] = g2;

                }

                if (phi[i] > phi[max_phi])
                {
                    max_phi = i;
                }

            }
            return phi[max_phi];

        }



    }
}
