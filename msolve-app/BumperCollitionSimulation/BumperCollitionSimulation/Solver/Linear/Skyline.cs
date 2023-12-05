using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BumperCollitionSimulation
{
    class Skyline : LinearSolution
    {
        public static (double[], int[]) Create_Skyline(double[,] A)
        {
            int[] band_per_column = new int[A.GetLength(1)];
            int nzv = 0;
            int[] diag_offsets = new int[A.GetLength(1) + 1];
            for (int j = 0; j < A.GetLength(1); j++)
            {
                int i = 0;
                bool flag = true;
                while (i <= j & flag)
                {
                    if (A[i, j] != 0)
                    {
                        flag = false;
                        band_per_column[j] = j - i;
                        diag_offsets[j] = nzv;
                        nzv = nzv + (j - i + 1);
                    }
                    i++;
                }
            }
            diag_offsets[A.GetLength(1)] = nzv;
            // non zero values
            double[] values = new double[nzv];
            int k = 0;
            for (int j = 0; j < A.GetLength(1); j++)
            {
                for (int i = j; i >= (j - (diag_offsets[j + 1] - diag_offsets[j] - 1)); i--)
                {
                    values[k] = A[i, j];
                    k++;
                }
            }
            return (values, diag_offsets);
        }
        private static void Cholesky_Skyline(double[] values, int[] diag_offsets, int n)
        {
            for (int i = 0; i < n; i++)
            {
                int hi = diag_offsets[i + 1] - diag_offsets[i] - 1;
                int mi = i - hi;
                for (int j = mi; j < i; j++)
                {
                    int hj = diag_offsets[j + 1] - diag_offsets[j] - 1;
                    int mj = j - hj;
                    for (int k = Math.Max(mi, mj); k < j; k++)
                    {
                        values[diag_offsets[i] + i - j] = values[diag_offsets[i] + i - j] - (values[diag_offsets[i] + i - k] * values[diag_offsets[j] + j - k]);
                    }
                    values[diag_offsets[i] + i - j] = values[diag_offsets[i] + i - j] / values[diag_offsets[j]];
                }
                for (int k = mi; k < i; k++)
                {
                    values[diag_offsets[i]] = values[diag_offsets[i]] - (values[diag_offsets[i] + i - k] * values[diag_offsets[i] + i - k]);
                }
                values[diag_offsets[i]] = Math.Sqrt(values[diag_offsets[i]]);
            }
        }
        private static double[] Solve_system(double[] L_skyline, int[] diag_offsets, double[] b, int n)
        {
            // Forward Substitution
            double[] temp = new double[b.Length];
            for (int i = 0; i < n; i++)
            {
                temp[i] = b[i];
                int hi = diag_offsets[i + 1] - diag_offsets[i] - 1;
                for (int j = 0; j < i; j++)
                {
                    if (Math.Abs(i - j) <= hi)
                    {
                        temp[i] = temp[i] - (L_skyline[diag_offsets[i] + i - j] * temp[j]);
                    }
                }
                temp[i] = temp[i] / L_skyline[diag_offsets[i]];
            }
            //Backward Subtitution with lower triangular L
            double[] solution = new double[temp.Length];
            for (int i = (n - 1); i >= 0; i--)
            {
                solution[i] = temp[i];
                for (int j = i + 1; j < n; j++)
                {
                    int hj = diag_offsets[j + 1] - diag_offsets[j] - 1;
                    if (Math.Abs(i - j) <= hj)
                    {
                        solution[i] = solution[i] - L_skyline[diag_offsets[j] + j - i] * solution[j]; // pay attetion to the order of i,j
                    }
                }
                solution[i] = solution[i] / L_skyline[diag_offsets[i]];
            }
            return solution;
        }
        public static void Test()
        {
            double[,] A = { { 21, 1, 0, 4, 0 }, { 1, 22, 2, 0, 0 }, { 0, 2, 23, 1, 3 }, { 4, 0, 1, 24, 2 }, { 0, 0, 3, 2, 25 } };
            (double[] values, int[] diag_offsets) = Create_Skyline(A);
            Cholesky_Skyline(values, diag_offsets, A.GetLength(0));
            double[] b = { 1, 2, 3, 4, 5 };
            double[] solution = Solve_system(values, diag_offsets, b, A.GetLength(0));
            double[] b_calc = VectorOperations.MatrixVectorProduct(A, solution);
            CheckEquality.ManualCheckArrays(b, b_calc);
        }
        public override double[] Solve(double[,] stiffnessMatrix, double[] rhsVector)
        {
            (double[] values, int[] diag_offsets) = Create_Skyline(stiffnessMatrix);
            Cholesky_Skyline(values, diag_offsets, stiffnessMatrix.GetLength(0));
            double[] solution = Solve_system(values, diag_offsets, rhsVector, stiffnessMatrix.GetLength(0));
            return solution;
        }
    }
}