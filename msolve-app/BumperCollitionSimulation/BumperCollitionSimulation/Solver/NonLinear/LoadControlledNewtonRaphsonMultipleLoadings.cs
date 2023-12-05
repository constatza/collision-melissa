using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading;

namespace BumperCollitionSimulation
{
    public class LoadControlledNewtonRaphsonMultipleLoadings : NonLinearSolution
    {
        private double[] localSolutionVector;

        public LoadControlledNewtonRaphsonMultipleLoadings()
        {

        }
        public LoadControlledNewtonRaphsonMultipleLoadings(double[] exSolution)
        {
            localSolutionVector = exSolution;
        }
        private double[] LoadControlledNRMultipleLoadings(Dictionary<int, double[]> forceVectors, Dictionary<int, int> numberOfLoadStepsForEachLoad)
        {
            var solutionsCount = 1;
            double[] solutionVector = localSolutionVector;
            discretization.InitializeContactTangentialProperties();
            discretization.InitializeContactSurfaceVectors();
            for (var loadingCounter = 1; loadingCounter <= forceVectors.Count; loadingCounter++)
            {
                lambda = 1.0 / numberOfLoadStepsForEachLoad[loadingCounter];
                double[] incrementDf = VectorOperations.VectorScalarProductNew(forceVectors[loadingCounter], lambda);
                double[] incrementalExternalForcesVector = new double[forceVectors[1].Length];
                if (loadingCounter > 1)
                {
                    for (var localCounter = 1; localCounter < loadingCounter; localCounter++)
                    {
                        incrementalExternalForcesVector = VectorOperations.VectorVectorAddition(incrementalExternalForcesVector, forceVectors[localCounter]);
                    }
                }
                double[] tempSolutionVector = new double[solutionVector.Length];
                double[] deltaU = new double[solutionVector.Length];
                double[] internalForcesTotalVector;
                double[] dU;
                double[] residual;
                double residualNorm;
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
                    InternalForces.Add(solutionsCount, internalForcesTotalVector);
                    solutionVector = VectorOperations.VectorVectorAddition(solutionVector, deltaU);
                    Solutions.Add(solutionsCount, solutionVector);
                    solutionsCount += 1;
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
            }
            return solutionVector;
        }

        public override double[] Solve(IAssembly assembly, ILinearSolution linearScheme, Dictionary<int, double[]> forceVectors,
            Dictionary<int, int> numberOfLoadStepsForEachLoad)
        {
            InternalForces = new Dictionary<int, double[]>();
            Solutions = new Dictionary<int, double[]>();
            LoadStepConvergence = new List<string>();
            if (localSolutionVector == null)
            {
                localSolutionVector = new double[forceVectors[1].Length];
            }
            discretization = assembly;
            linearSolver = linearScheme;
            lambda = 1.0 / numberOfLoadSteps;
            double[] solution = LoadControlledNRMultipleLoadings(forceVectors, numberOfLoadStepsForEachLoad);
            return solution;
        }

    }
}