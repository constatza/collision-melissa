using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading;

namespace BumperCollitionSimulation
{
    public class ModifiedNewtonRaphson : NonLinearSolution
    {
        private double[] localSolutionVector;

        public ModifiedNewtonRaphson()
        {

        }
        public ModifiedNewtonRaphson(double[] exSolution)
        {
            localSolutionVector = exSolution;
        }
        private double[] LoadControlledModifiedNR(double[] forceVector)
        {
            double[] incrementDf = VectorOperations.VectorScalarProductNew(forceVector, lambda);
            double[] solutionVector = localSolutionVector;
            double[] incrementalExternalForcesVector = new double[forceVector.Length];
            double[] tempSolutionVector = new double[solutionVector.Length];
            double[] deltaU = new double[solutionVector.Length];
            double[] internalForcesTotalVector;
            double[] dU;
            double[] residual;
            //
            //double[] oldInternalForces;
            //double relativeInternalForcesNorm;
            //
            //List<int> iterationsToUpdateStiffnessMatrix = new List<int>();
            //for (int i = 0; i <= 7; i ++)
            //{
            //    iterationsToUpdateStiffnessMatrix.Add(i);
            //}
            //for (int i = 15; i <= MaxIterations - 5; i += 20)
            //{
            //    iterationsToUpdateStiffnessMatrix.Add(i);
            //    iterationsToUpdateStiffnessMatrix.Add(i + 1);
            //    iterationsToUpdateStiffnessMatrix.Add(i + 2);

            //}
            double residualNorm;
            discretization.InitializeContactTangentialProperties();
            discretization.InitializeContactSurfaceVectors();
            for (int i = 0; i < numberOfLoadSteps; i++)
            {
                int count = 0;
                int countIncrements = 1;
                //int countIncrements2 = 0;
                //int countIncrements3 = 0;
                //int countIncrements0 = 1;
                //int countIncrements1 = 1;
                //int countIncrements2 = 1;
                //int countIncrements3 = 1;
                //int countIncrements4 = 1;

                incrementalExternalForcesVector = VectorOperations.VectorVectorAddition(incrementalExternalForcesVector, incrementDf);
                discretization.UpdateDisplacements(solutionVector);
                discretization.UpdateContactSurfaceVectors();
                internalForcesTotalVector = discretization.CreateTotalInternalForcesVector();
                double[,] stiffnessMatrix = discretization.CreateTotalStiffnessMatrix();
                //OnConvergenceResult("Newton-Raphson: Solution not converged at load step" + i); 
                dU = linearSolver.Solve(stiffnessMatrix, incrementDf);
                discretization.UpdateElementsIncrementalDisplacements(dU);
                solutionVector = VectorOperations.VectorVectorAddition(solutionVector, dU);
                residual = VectorOperations.VectorVectorSubtraction(internalForcesTotalVector, incrementalExternalForcesVector);
                residualNorm = VectorOperations.VectorNorm2(residual);
                int iteration = 0;
                //
                //oldInternalForces = internalForcesTotalVector;
                //
                Array.Clear(deltaU, 0, deltaU.Length);
                while (residualNorm > Tolerance && iteration < MaxIterations)
                {
                    if (residualNorm >= 10.00)
                    {
                        stiffnessMatrix = discretization.CreateTotalStiffnessMatrix();
                        count += 1;
                        countIncrements += 1;
                    }
                    else if (residualNorm >= 1.00)
                    {
                        if (countIncrements % 2 == 0)
                        {
                            stiffnessMatrix = discretization.CreateTotalStiffnessMatrix(); count += 1;
                        }
                        countIncrements += 1;
                    }
                    else if (1000.0 * Tolerance <= residualNorm && residualNorm < 1.00)
                    {
                        if (countIncrements % 5 == 0)
                        {
                            stiffnessMatrix = discretization.CreateTotalStiffnessMatrix(); count += 1;
                        }
                        countIncrements += 1;
                    }
                    else if (10.0 * Tolerance <= residualNorm && residualNorm < 1000.0 * Tolerance)
                    {
                        if (countIncrements % 10 == 0)
                        {
                            stiffnessMatrix = discretization.CreateTotalStiffnessMatrix(); count += 1;
                        }
                        countIncrements += 1;
                    }

                    deltaU = VectorOperations.VectorVectorSubtraction(deltaU, linearSolver.Solve(stiffnessMatrix, residual));
                    tempSolutionVector = VectorOperations.VectorVectorAddition(solutionVector, deltaU);
                    discretization.UpdateElementsIncrementalDisplacements(VectorOperations.VectorVectorAddition(dU, deltaU));
                    discretization.UpdateDisplacements(tempSolutionVector);
                    //oldInternalForces = internalForcesTotalVector;
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
                    iteration = iteration + 1;
                    //(Application.Current.Windows[0] as MainWindow).LogTool.Text = "ok"; 
                    //OnConvergenceResult("Newton-Raphson: Solution not converged at load step" + iteration);
                }
                InternalForces.Add(i + 1, internalForcesTotalVector);
                solutionVector = VectorOperations.VectorVectorAddition(solutionVector, deltaU);
                Solutions.Add(i + 1, solutionVector);
                discretization.UpdateContactTangentialProperties();
                //discretization.NextStepFrictionUpdate();
                if (iteration >= MaxIterations)
                {
                    OnConvergenceResult(new ConvergenceValues() { ResidualNorm = residualNorm, LoadStep = i, Iteration = iteration, Tolerance = Tolerance, ConvergenceResult = false });
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
            double[] solution = LoadControlledModifiedNR(forceVector);
            return solution;
        }

    }
}

