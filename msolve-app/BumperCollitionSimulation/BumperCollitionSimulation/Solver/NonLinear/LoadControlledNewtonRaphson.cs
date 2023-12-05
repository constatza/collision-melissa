using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading;

namespace BumperCollitionSimulation
{
    public class LoadControlledNewtonRaphson : NonLinearSolution
    {
        private double[] localSolutionVector;

        public LoadControlledNewtonRaphson()
        {

        }
        public LoadControlledNewtonRaphson(double[] exSolution)
        {
            localSolutionVector = exSolution;
        }
        private double[] LoadControlledNR(double[] forceVector)
        {
            double[] incrementDf = VectorOperations.VectorScalarProductNew(forceVector, lambda);
            double[] solutionVector = localSolutionVector;
            double[] incrementalExternalForcesVector = new double[forceVector.Length];
            double[] tempSolutionVector = new double[solutionVector.Length];
            double[] deltaU = new double[solutionVector.Length];
            double[] internalForcesTotalVector;
            double[] dU;
            double[] residual;
            double residualNorm;
            discretization.InitializeContactTangentialProperties();
            discretization.InitializeContactSurfaceVectors();
            double[,] stiffnessMatrixLinearPart = discretization.CreateStiffnessMatrixLinearPart();
            for (int i = 0; i < numberOfLoadSteps; i++)
            {
                incrementalExternalForcesVector = VectorOperations.VectorVectorAddition(incrementalExternalForcesVector, incrementDf);
                discretization.UpdateDisplacements(solutionVector);
                discretization.UpdateEASParameters(solutionVector);//added EAS
                discretization.UpdateContactSurfaceVectors();
                internalForcesTotalVector = discretization.CreateTotalInternalForcesVector();
                double[,] UpdatedStiffnessMatrix = discretization.UpdateStifnessMatrix(stiffnessMatrixLinearPart);
                dU = linearSolver.Solve(UpdatedStiffnessMatrix, incrementDf);
                discretization.UpdateElementsIncrementalDisplacements(dU);
                solutionVector = VectorOperations.VectorVectorAddition(solutionVector, dU);
                residual = VectorOperations.VectorVectorSubtraction(internalForcesTotalVector, incrementalExternalForcesVector);
                residualNorm = VectorOperations.VectorNorm2(residual);
                int iteration = 0;
                Array.Clear(deltaU, 0, deltaU.Length);
                //Array.Clear(deltaU, 0, deltaU.Length);
                while (residualNorm > Tolerance && iteration < MaxIterations)
                {
                    UpdatedStiffnessMatrix = discretization.UpdateStifnessMatrix(stiffnessMatrixLinearPart);
                    deltaU = VectorOperations.VectorVectorSubtraction(deltaU, linearSolver.Solve(UpdatedStiffnessMatrix, residual));
                    tempSolutionVector = VectorOperations.VectorVectorAddition(solutionVector, deltaU);
                    discretization.UpdateElementsIncrementalDisplacements(VectorOperations.VectorVectorAddition(dU, deltaU));
                    discretization.UpdateDisplacements(tempSolutionVector);
                    discretization.UpdateEASParameters(tempSolutionVector);//added EAS
                    internalForcesTotalVector = discretization.CreateTotalInternalForcesVector();
                    residual = VectorOperations.VectorVectorSubtraction(internalForcesTotalVector, incrementalExternalForcesVector);
                    residualNorm = VectorOperations.VectorNorm2(residual);
                    if (residualNorm <= Tolerance)
                    {
                        OnConvergenceResult(new ConvergenceValues() { ResidualNorm = residualNorm, LoadStep = i, Iteration = iteration, Tolerance = Tolerance, ConvergenceResult = true });
                    }
                    else
                    {
                        OnConvergenceResult(new ConvergenceValues() { ResidualNorm = residualNorm, LoadStep = i, Iteration = iteration, Tolerance = Tolerance, ConvergenceResult = false });
                    }
                    iteration += 1;
                }
                InternalForces.Add(i + 1, internalForcesTotalVector);
                solutionVector = VectorOperations.VectorVectorAddition(solutionVector, deltaU);
                Solutions.Add(i + 1, solutionVector);
                discretization.UpdateContactTangentialProperties();
                //VectorOperations.PrintVectorToFile(solutionVector, @"C:\Users\Public\Documents\" + "Solution37.dat");

                //discretization.NextStepFrictionUpdate();
                if (iteration >= MaxIterations)
                {
                    OnConvergenceResult(new ConvergenceValues() { ResidualNorm = residualNorm, LoadStep = i, Iteration = iteration, Tolerance = Tolerance, ConvergenceResult = true });
                    LoadStepConvergence.Add("Solution not converged.");
                    break;
                }
                LoadStepConvergence.Add("Solution converged.");
            }
            return solutionVector;
        }

        public override double[] Solve(IAssembly assembly, ILinearSolution linearScheme, double[] forceVector)
        {
            InternalForces = new Dictionary<int, double[]>();
            Solutions = new Dictionary<int, double[]>();
            LoadStepConvergence = new List<string>();
            if (localSolutionVector == null)
            {
                localSolutionVector = new double[forceVector.Length];
            }
            discretization = assembly;
            linearSolver = linearScheme;
            lambda = 1.0 / numberOfLoadSteps;
            //double[] solution = null;

            //Thread tcore1 = new Thread(() => LoadControlledNR(forceVector));
            //tcore1.Start();
            //tcore1.Join();
            double[] solution = LoadControlledNR(forceVector);
            return solution;
        }

    }
}