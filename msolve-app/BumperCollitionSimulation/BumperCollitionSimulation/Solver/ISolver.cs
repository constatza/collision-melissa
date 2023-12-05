﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace BumperCollitionSimulation
{
    public interface ISolver
    {
        ILinearSolution LinearScheme { get; set; }
        IAssembly AssemblyData { get; set; }
        void Solve(double[] rhsVector);
        void Solve(Dictionary<int, double[]> rhsVectors, Dictionary<int, int> numberOfLoadStepsForEachLoad);
        bool ActivateNonLinearSolver { get; set; }
        INonLinearSolution NonLinearScheme { get; set; }
        void PrintSolution();
        double[,] CustomStiffnessMatrix { get; set; }
        double[] GetSolution();
        Dictionary<int, double[]> GetInternalForces();
        Dictionary<int, double[]> GetAllStepsSolutions();
        Dictionary<int, double[]> GetAllStepsFullSolutionVectors();
    }
}