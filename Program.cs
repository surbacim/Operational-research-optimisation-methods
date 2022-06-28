using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace OperatResAndOptimMet
{
    class Program
    {
        public static Vector ffff(Vector x)
        {
            Vector fc = new Vector(3);
            fc[0] = x[0] * x[0] + x[0] * x[1] + x[1] * x[1];
            fc[1] = x[0] * x[0] + x[1];
            fc[2] = x[1] + 2 * x[1] * x[1];
            return fc;
        }

        static void Main(string[] args)
        {

            Console.WriteLine("Численный метод поиска локального экстремума: {0}", Optimization.StepByStep(1, 0.5, 0.0001, x => x + 8 / Math.Pow(x, 4)));
            Console.WriteLine("Золотое сечение: {0}", Optimization.GoldenRatio(1, 5, 0.0001, x => x + 8 / Math.Pow(x, 4)));
            
            Vector xn = new Vector(2);
            xn[0] = 1;
            xn[1] = 1;

            int[] tipogr = { 0, 1, 2 };
            Vector fv = new Vector(3);
            fv[0] = 4; fv[1] = 1; fv[2] = 5;

            Vector fn = new Vector(3);
            fn[0] = 8; fn[1] = 3; fn[2] = 4;

            Console.WriteLine("\nМетод квадратичной аппроксимации: {0}", Optimization.MethodKvApproksim(-5, 0, 5, 0.0001, x => x * x + 5 * x));

            Console.WriteLine("Градиентный метод поиска минимума: {0}", Optimization.Grad(xn, 0.1, 0.0001, x => x[0] * x[0] + x[1] * x[0] + 2 * x[1] * x[1]));

            Console.WriteLine("Метод наискорейшего градиентного спуска (МНГС): {0}", Optimization.MNGS(xn, 0.0001, x => x[0] * x[0] + x[1] * x[0] + 2 * x[1] * x[1]));

            Console.WriteLine("Метод сопряженных направлений (МСН): {0}", Optimization.MSN(xn, 0.0001, x => x[0] * x[0] + x[1] * x[0] + 2 * x[1] * x[1]));

            Console.WriteLine("Метод случайного поиска (МСП): {0}", Optimization.MSP_best_sampling(xn, 2, 0.1, 0.0001, x => x[0] * x[0] + x[0] * x[1] + 2 * x[1] * x[1]));

            Console.WriteLine("ОЗУ {0}", Optimization.OZU(xn, 2, 0.1, 0.001, fv, fn, tipogr, ffff));

            Vector[] xnn = new Vector[3];
            xnn[0] = xn;
            double[] xnnn = new double[2];
            xnnn[0] = 0; xnnn[1] = 0;
            xnn[1] = new Vector(xnnn);
            xnnn[0] = 1; xnnn[1] = 0;
            xnn[2] = new Vector(xnnn);

            Console.WriteLine("Метод деформируемого многогранника: {0}", Optimization.Nelder_Mead(xnn, 0.001, x => x[0] * x[0] + x[1] * x[0] + x[1] * x[1] - 9 * x[1] - 6 * x[0]));



            Console.ReadKey();
            

        }
    }
}
