﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace BumperCollitionSimulation
{
    class StaticSolver : ISolver
    {
        private double[] staticSolutionVector;
        public IAssembly AssemblyData { get; set; }
        public ILinearSolution LinearScheme { get; set; }
        public INonLinearSolution NonLinearScheme { get; set; }
        public bool ActivateNonLinearSolver { get; set; }
        public double[,] CustomStiffnessMatrix { get; set; }

        public void Solve(double[] rhsVector)
        {
            if (ActivateNonLinearSolver == true)
            {
                staticSolutionVector = NonLinearScheme.Solve(AssemblyData, LinearScheme, rhsVector);
            }
            else
            {
                if (CustomStiffnessMatrix == null)
                {
                    double[,] coefMatrix = AssemblyData.CreateTotalStiffnessMatrix();
                    staticSolutionVector = LinearScheme.Solve(coefMatrix, rhsVector);
                }
                //double[,] coefMatrix = AssemblyData.CreateTotalStiffnessMatrix();
                //staticSolutionVector = LinearScheme.Solve(coefMatrix, rhsVector);
                else
                {
                    double[,] coefMatrix = CustomStiffnessMatrix;
                    staticSolutionVector = LinearScheme.Solve(coefMatrix, rhsVector);
                }
            }
        }
        public void Solve(Dictionary<int, double[]> rhsVectors, Dictionary<int, int> numberOfLoadStepsForEachLoad)
        {
            if (ActivateNonLinearSolver == true)
            {
                staticSolutionVector = NonLinearScheme.Solve(AssemblyData, LinearScheme, rhsVectors, numberOfLoadStepsForEachLoad);
            }
            else
            {
                if (CustomStiffnessMatrix == null)
                {
                    double[,] coefMatrix = AssemblyData.CreateTotalStiffnessMatrix();
                    staticSolutionVector = LinearScheme.Solve(coefMatrix, rhsVectors[1]);
                }
                //double[,] coefMatrix = AssemblyData.CreateTotalStiffnessMatrix();
                //staticSolutionVector = LinearScheme.Solve(coefMatrix, rhsVector);
                else
                {
                    double[,] coefMatrix = CustomStiffnessMatrix;
                    staticSolutionVector = LinearScheme.Solve(coefMatrix, rhsVectors[1]);
                }
            }
        }
        public void PrintSolution()
        {
            VectorOperations.PrintVector(staticSolutionVector);
        }

        public double[] GetSolution()
        {
            return staticSolutionVector;
        }

        public Dictionary<int, double[]> GetInternalForces()
        {
            return NonLinearScheme.InternalForces;
        }

        public Dictionary<int, double[]> GetAllStepsSolutions()
        {
            return NonLinearScheme.Solutions;
        }

        public Dictionary<int, double[]> GetAllStepsFullSolutionVectors()
        {

            Dictionary<int, double[]> allStepsSolutions = GetAllStepsSolutions();
            Dictionary<int, double[]> result = new Dictionary<int, double[]>();
            foreach (KeyValuePair<int, double[]> solutionVector in allStepsSolutions)
            {
                double[] fullSolutionVector = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(solutionVector.Value, AssemblyData.BoundedDOFsVector);
                result.Add(solutionVector.Key, fullSolutionVector);
            }
            return result;
        }
    }
}