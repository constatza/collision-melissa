using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace BumperCollitionSimulation
{
    public interface ILinearSolution
    {
        double[] Solve(double[,] stiffnessMatrix, double[] forceVector);
        double[] Solve(double[,] stiffnessMatrixK, double[,] stiffnessMatrixB, double[,] stiffnessMatrixC,
            double[] forceVectorM, double[] forceVectorC, int [,] fillInLevels, int fillLevel);

        Tuple<double[], int> Solve(double[,] stiffnessMatrixK, double[,] stiffnessMatrixB, double[,] stiffnessMatrixC,
            double[] forceVectorM, double[] forceVectorC, int[,] fillInLevels, int fillLevel, bool countIterations);

        Tuple<double[], int> Solve(double[,] stiffnessMatrix, double[] forceVector, bool countIterations);

    }
}