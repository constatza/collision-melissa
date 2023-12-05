using System;
using System.Collections.Generic;
using System.Text;
using System.Diagnostics;


namespace BumperCollitionSimulation
{
    static class CheckEquality
    {
        public static void CheckArrays(double[] y1, double[] y2)
        {
            for (int i = 0; i < y1.Length; i++)
            {
                Debug.Assert(y1[i] == y2[i]);
            }
        }
        public static void CheckMatrices(double[,] M1, double[,] M2)
        {
            for (int i = 0; i < M1.GetLength(0); i++)
            {
                for (int j = 0; j < M1.GetLength(1); j++)
                {
                    Debug.Assert(M1[i, j] == M2[i, j]);
                }
            }
        }
        public static void ManualCheckArrays(double[] y1, double[] y2)
        {
            for (int i = 0; i < y1.Length; i++)
            {
                if (Math.Abs(y1[i] - y2[i]) > 1e-7)
                {
                    throw new Exception("Not equal elements");
                }
            }
        }
        public static void ManualCheckMatrices(double[,] M1, double[,] M2)
        {
            for (int i = 0; i < M1.GetLength(0); i++)
            {
                for (int j = 0; j < M1.GetLength(1); j++)
                {
                    if (Math.Abs(M1[i, j] - M2[i, j]) > 1e-8)
                    {
                        throw new Exception("Not equal elements");
                    }
                }
            }
        }
    }
}
