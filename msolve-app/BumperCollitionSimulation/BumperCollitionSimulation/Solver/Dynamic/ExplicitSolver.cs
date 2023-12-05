﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace BumperCollitionSimulation
{
    public class ExplicitSolver
    {
        public IAssembly Assembler { get; set; }
        public ILinearSolution LinearSolver { get; set; }
        public double[] ExternalForcesVector { get; set; }
        public InitialConditions InitialValues { get; set; }
        private int numberOfLoadSteps = 20;
        private double tolerance = 1.0e-4;
        private int maxIterations = 1000;
        private double lambda;
        private double totalTime;
        private int timeStepsNumber;
        private double timeStep;
        private double a0, a1, a2, a3;
        public Dictionary<int, double[]> explicitSolution = new Dictionary<int, double[]>();
        public Dictionary<int, double[]> displacement = new Dictionary<int, double[]>();
        private Dictionary<int, double[]> explicitVelocity = new Dictionary<int, double[]>();
        private Dictionary<int, double[]> explicitAcceleration = new Dictionary<int, double[]>();
        public bool ActivateNonLinearSolution { get; set; }
        public double[,] CustomStiffnessMatrix { get; set; }
        public double[,] CustomMassMatrix { get; set; }
        public double[,] CustomDampingMatrix { get; set; }
        public Dictionary<int, double> TimeAtEachStep { get; set; }

        public ExplicitSolver(double totalTime, int timeStepsNumber)
        {
            this.totalTime = totalTime;
            this.timeStepsNumber = timeStepsNumber;
            timeStep = totalTime / timeStepsNumber;
            a0 = 1.0 / (timeStep * timeStep);
            a1 = 1.0 / (2.0 * timeStep);
            a2 = 2.0 * a0;
            a3 = 1.0 / a2;
            TimeAtEachStep = new Dictionary<int, double>();
        }

        private double[] LoadControlledNR(double[] forceVector)
        {

            lambda = 1.0 / numberOfLoadSteps;
            double[] incrementDf = VectorOperations.VectorScalarProductNew(forceVector, lambda);
            double[] solutionVector = new double[forceVector.Length];
            double[] incrementalExternalForcesVector = new double[forceVector.Length];
            double[] tempSolutionVector = new double[solutionVector.Length];
            double[] deltaU = new double[solutionVector.Length];
            double[] internalForcesTotalVector;
            double[] dU;
            double[] residual;
            double residualNorm;
            //Assembler.UpdateAccelerations(CalculateAccelerations(InitialValues.InitialAccelerationVector));

            for (int i = 0; i < numberOfLoadSteps; i++)
            {
                incrementalExternalForcesVector = VectorOperations.VectorVectorAddition(incrementalExternalForcesVector, incrementDf);
                Assembler.UpdateDisplacements(solutionVector);

                Assembler.UpdateAccelerations(explicitAcceleration.Values.Last());

                internalForcesTotalVector = Assembler.CreateTotalInternalForcesVector();

                double[,] tangentMatrix = CalculateHatMMatrix();
                dU = LinearSolver.Solve(tangentMatrix, incrementDf);
                solutionVector = VectorOperations.VectorVectorAddition(solutionVector, dU);

                Assembler.UpdateDisplacements(solutionVector);
                tangentMatrix = CalculateHatMMatrix();
                internalForcesTotalVector = Assembler.CreateTotalInternalForcesVector();

                residual = VectorOperations.VectorVectorSubtraction(internalForcesTotalVector, incrementalExternalForcesVector);
                residualNorm = VectorOperations.VectorNorm2(residual);
                int iteration = 0;
                Array.Clear(deltaU, 0, deltaU.Length);
                while (residualNorm > tolerance && iteration < maxIterations)
                {
                    tangentMatrix = CalculateHatMMatrix();
                    deltaU = VectorOperations.VectorVectorSubtraction(deltaU, LinearSolver.Solve(tangentMatrix, residual));
                    tempSolutionVector = VectorOperations.VectorVectorAddition(solutionVector, deltaU);
                    Assembler.UpdateDisplacements(tempSolutionVector);

                    //Assembler.UpdateAccelerations(CalculateAccelerations(solutionVector));

                    internalForcesTotalVector = Assembler.CreateTotalInternalForcesVector();
                    residual = VectorOperations.VectorVectorSubtraction(internalForcesTotalVector, incrementalExternalForcesVector);
                    residualNorm = VectorOperations.VectorNorm2(residual);
                    iteration = iteration + 1;
                }
                solutionVector = VectorOperations.VectorVectorAddition(solutionVector, deltaU);
                if (iteration >= maxIterations) Console.WriteLine("Newton-Raphson: Solution not converged at current iterations");
            }

            return solutionVector;
        }

        private double[] NewtonIterations(double[] forceVector)
        {

            lambda = 1.0 / numberOfLoadSteps;
            double[] incrementDf = VectorOperations.VectorScalarProductNew(forceVector, lambda);
            double[] solutionVector = new double[forceVector.Length];
            double[] deltaU = new double[solutionVector.Length];
            double[] internalForcesTotalVector;
            double[] dU;
            double[] residual;
            double residualNorm;

            solutionVector = explicitSolution.Values.Last();

            Assembler.UpdateDisplacements(solutionVector);

            Assembler.UpdateAccelerations(explicitAcceleration.Values.Last());
            residual = VectorOperations.VectorVectorSubtraction(forceVector, Assembler.CreateTotalInternalForcesVector());
            int iteration = 0;
            Array.Clear(deltaU, 0, deltaU.Length);

            for (int i = 0; i < maxIterations; i++)
            {
                double[,] tangentMatrix = CalculateHatMMatrix();
                deltaU = LinearSolver.Solve(tangentMatrix, residual);
                solutionVector = VectorOperations.VectorVectorAddition(solutionVector, deltaU);
                Assembler.UpdateDisplacements(solutionVector);
                //Assembler.UpdateAccelerations(CalculateAccelerations());

                internalForcesTotalVector = Assembler.CreateTotalInternalForcesVector();
                residual = VectorOperations.VectorVectorSubtraction(forceVector, internalForcesTotalVector);
                residualNorm = VectorOperations.VectorNorm2(residual);
                if (residualNorm < tolerance)
                {
                    break;
                }
                iteration = iteration + 1;
            }
            Console.WriteLine(iteration);
            if (iteration >= maxIterations) Console.WriteLine("Newton-Raphson: Solution not converged at current iterations");

            return solutionVector;
        }

        private double[] NewtonIterationsExplicit(int timeStep, double[] forceVector, double[,] tangentMatrix)
        {
            //lambda = 1.0 / numberOfLoadSteps;
            //double[] incrementDf = VectorOperations.VectorScalarProductNew(forceVector, lambda);
            double[] solutionVector = new double[forceVector.Length];
            double[] deltaU = new double[solutionVector.Length];
            //double[] internalForcesTotalVector;
            double[] residual;
            double residualNorm;
            double[] hatRPrevious;
            double[] hatRNext;

            solutionVector = explicitSolution.Values.Last();

            Assembler.UpdateDisplacements(solutionVector);

            Assembler.UpdateAccelerations(explicitAcceleration.Values.Last());
            hatRPrevious = CalculateHatRVectorNL(timeStep);
            hatRNext = hatRPrevious;
            //internalForcesTotalVector = Assembler.CreateTotalInternalForcesVector();
            //residual = VectorOperations.VectorVectorSubtraction(hatR, internalForcesTotalVector);
            residual = hatRPrevious;
            int iteration = 0;
            Array.Clear(deltaU, 0, deltaU.Length);

            for (int i = 0; i < maxIterations; i++)
            {
                if (i == 0)
                {
                    //deltaU = LinearSolver.Solve(tangentMatrix, residual);
                    //solutionVector = VectorOperations.VectorVectorAddition(solutionVector, deltaU);
                    solutionVector = LinearSolver.Solve(tangentMatrix, residual);
                    Assembler.UpdateDisplacements(solutionVector);
                    //Assembler.UpdateAccelerations(CalculateAccelerations());
                    hatRNext = CalculateHatRVectorNL(timeStep);
                    //internalForcesTotalVector = Assembler.CreateTotalInternalForcesVector();
                    residual = VectorOperations.VectorVectorSubtraction(hatRNext, hatRPrevious);
                    residualNorm = VectorOperations.VectorNorm2(residual);
                    if (residualNorm < tolerance)
                    {
                        break;
                    }
                    iteration = iteration + 1;
                    hatRPrevious = hatRNext;
                }
                else
                {
                    deltaU = LinearSolver.Solve(tangentMatrix, residual);
                    solutionVector = VectorOperations.VectorVectorAddition(solutionVector, deltaU);
                    //solutionVector = LinearSolver.Solve(tangentMatrix, residual);
                    Assembler.UpdateDisplacements(solutionVector);
                    //Assembler.UpdateAccelerations(CalculateAccelerations());
                    hatRNext = CalculateHatRVectorNL(timeStep);
                    //internalForcesTotalVector = Assembler.CreateTotalInternalForcesVector();
                    residual = VectorOperations.VectorVectorSubtraction(hatRNext, hatRPrevious);
                    residualNorm = VectorOperations.VectorNorm2(residual);
                    if (residualNorm < tolerance)
                    {
                        break;
                    }
                    iteration = iteration + 1;
                    hatRPrevious = hatRNext;
                }
            }
            //Console.WriteLine(iteration);
            if (iteration >= maxIterations) Console.WriteLine("Newton-Raphson: Solution not converged at current iterations");

            return solutionVector;
        }
        //private double[] UpdatedF(double[] forceVector)
        //{


        //    lambda = 1.0 / numberOfLoadSteps;
        //    double[] incrementDf = VectorOperations.VectorScalarProductNew(forceVector, lambda);
        //    double[] solutionVector = new double[forceVector.Length];
        //    double[] incrementalExternalForcesVector = new double[forceVector.Length];
        //    double[] tempSolutionVector = new double[solutionVector.Length];
        //    double[] deltaU = new double[solutionVector.Length];
        //    double[] internalForcesTotalVector;
        //    double[] dU;





        //    double[] currentU = explicitSolution.Values.Last();
        //    Assembler.UpdateDisplacements(currentU);
        //    internalForcesTotalVector = Assembler.CreateTotalInternalForcesVector();
        //    double[,] TotalMassMatrix = Assembler.CreateTotalMassMatrix();
        //    double[] a2_M_ut = VectorOperations.MatrixVectorProduct(
        //                        MatrixOperations.ScalarMatrixProductNew(a2, TotalMassMatrix), currentU);

        //    double[] a0_M_ut = VectorOperations.MatrixVectorProduct(
        //                        MatrixOperations.ScalarMatrixProductNew(a0, TotalMassMatrix), currentU);

        //    double[,] hatMMatrix = CalculateHatMMatrix();
        //    double[] hatM_utprevious = VectorOperations.MatrixVectorProduct(hatMMatrix, )




        //double[,] hatK = MatrixOperations.MatrixSubtraction(TotalStiffnessMatrix,
        //                    MatrixOperations.ScalarMatrixProductNew(a2, TotalMassMatrix));
        //double[] hatCurrentU = VectorOperations.MatrixVectorProduct(hatKMatrix, explicitSolution[i - 1]);
        //    double[] hatPreviousU = VectorOperations.MatrixVectorProduct(hatMMatrix, explicitSolution[i - 2]);

        //    double[] hatR = VectorOperations.VectorVectorSubtraction(ExternalForcesVector,
        //                    VectorOperations.VectorVectorAddition(hatCurrentU, hatPreviousU));



        //    double[,] tangentMatrix = CalculateHatMMatrix();
        //        dU = LinearSolver.Solve(tangentMatrix, incrementDf);
        //        solutionVector = VectorOperations.VectorVectorAddition(solutionVector, dU);

        //        Assembler.UpdateDisplacements(solutionVector);
        //        tangentMatrix = CalculateHatMMatrix();
        //        internalForcesTotalVector = Assembler.CreateTotalInternalForcesVector();

        //        residual = VectorOperations.VectorVectorSubtraction(internalForcesTotalVector, incrementalExternalForcesVector);
        //        residualNorm = VectorOperations.VectorNorm2(residual);
        //        int iteration = 0;
        //        Array.Clear(deltaU, 0, deltaU.Length);



        //    return solutionVector;
        //}

        /// <summary>
        /// Calculates accelerations for time t
        /// </summary>
        /// <returns></returns>
        private double[] CalculateAccelerations() //Bathe page 771
        {
            int steps = explicitSolution.Count;
            double[] aCurrent =
                VectorOperations.VectorScalarProductNew(
                    VectorOperations.VectorVectorAddition(explicitSolution[steps - 4],
                        VectorOperations.VectorVectorAddition(
                            VectorOperations.VectorScalarProductNew(explicitSolution[steps - 3], -2.0), explicitSolution[steps - 2])), a0);

            return aCurrent;
        }

        private double[] CalculateInitialAccelerations() //Bathe page 771
        {
            if (CustomStiffnessMatrix != null) return InitialValues.InitialAccelerationVector;
            int step = explicitSolution.Count - 2;
            Assembler.UpdateDisplacements(explicitSolution[step]);
            double[,] stiffness = Assembler.CreateTotalStiffnessMatrix();
            double[,] mass = Assembler.CreateTotalMassMatrix();

            double[] Ku = VectorOperations.MatrixVectorProduct(stiffness, explicitSolution[step]);
            double[] RHS = VectorOperations.VectorVectorSubtraction(ExternalForcesVector, Ku);

            double[] acceleration = LinearSolver.Solve(mass, RHS);

            return acceleration;
        }

        private double[] CalculatePreviousDisplacementVector()
        {
            double[] previousDisp = VectorOperations.VectorVectorAddition(
                                    VectorOperations.VectorVectorSubtraction(InitialValues.InitialDisplacementVector,
                                    VectorOperations.VectorScalarProductNew(InitialValues.InitialVelocityVector, timeStep)),
                                    VectorOperations.VectorScalarProductNew(InitialValues.InitialAccelerationVector, a3));
            return previousDisp;
        }

        private double[,] CalculateHatMMatrix()
        {
            double[,] TotalMassMatrix;
            double[,] TotalDampingMatrix;
            if (CustomMassMatrix != null)
            {
                TotalMassMatrix = CustomMassMatrix;
                TotalDampingMatrix = CustomDampingMatrix;
            }
            else
            {
                TotalMassMatrix = Assembler.CreateTotalMassMatrix();
                TotalDampingMatrix = Assembler.CreateTotalDampingMatrix();
            }
            double[,] a0M = MatrixOperations.ScalarMatrixProductNew(a0, TotalMassMatrix);
            double[,] a1C = MatrixOperations.ScalarMatrixProductNew(a1, TotalDampingMatrix);
            double[,] hutM = MatrixOperations.MatrixAddition(a0M, a1C);
            return hutM;
        }

        private double[,] CalculateHatKMatrix()
        {
            double[,] TotalMassMatrix;
            double[,] TotalStiffnessMatrix;
            if (CustomStiffnessMatrix != null)
            {
                TotalMassMatrix = CustomMassMatrix;
                TotalStiffnessMatrix = CustomStiffnessMatrix;
            }
            else
            {
                TotalMassMatrix = Assembler.CreateTotalMassMatrix();
                TotalStiffnessMatrix = Assembler.CreateTotalStiffnessMatrix();
            }
            double[,] hatK = MatrixOperations.MatrixSubtraction(TotalStiffnessMatrix,
                                MatrixOperations.ScalarMatrixProductNew(a2, TotalMassMatrix));
            return hatK;
        }

        private double[] CalculateHatRVector(int i)
        {
            double[,] hatKMatrix = CalculateHatKMatrix();
            double[,] hatMMatrix = CalculateHatMMatrix();
            double[] hatCurrentU = VectorOperations.MatrixVectorProduct(hatKMatrix, explicitSolution[i - 1]);
            double[] hatPreviousU = VectorOperations.MatrixVectorProduct(hatMMatrix, explicitSolution[i - 2]);

            double[] hatR = VectorOperations.VectorVectorSubtraction(ExternalForcesVector,
                            VectorOperations.VectorVectorAddition(hatCurrentU, hatPreviousU));
            return hatR;
        }

        private double[] CalculateHatRVectorNL(int i)
        {
            double[] solutionVector = explicitSolution.Values.Last();
            Assembler.UpdateDisplacements(solutionVector);
            double[,] totalMassMatrix = Assembler.CreateTotalMassMatrix();
            double[,] totalDampingMatrix = Assembler.CreateTotalDampingMatrix();
            double[,] a2M = MatrixOperations.ScalarMatrixProductNew(a2, totalMassMatrix);
            double[,] a0M = MatrixOperations.ScalarMatrixProductNew(a0, totalMassMatrix);
            double[,] a1C = MatrixOperations.ScalarMatrixProductNew(a1, totalDampingMatrix);
            double[,] hutM = MatrixOperations.MatrixAddition(MatrixOperations.ScalarMatrixProductNew(-1.0, a0M), a1C);

            double[] F = Assembler.CreateTotalInternalForcesVector();
            double[] hatPreviousU = VectorOperations.MatrixVectorProduct(hutM, explicitSolution[i - 2]);
            double[] a2Mut = VectorOperations.MatrixVectorProduct(a2M, explicitSolution[i - 1]);


            double[] hatR1 = VectorOperations.VectorVectorSubtraction(ExternalForcesVector, F);
            double[] hatR2 = VectorOperations.VectorVectorAddition(a2Mut, hatPreviousU);
            double[] hatRtotal = VectorOperations.VectorVectorAddition(hatR1, hatR2);
            return hatRtotal;
        }


        public void SolveExplicit()
        {
            double[,] hatMassMatrix = CalculateHatMMatrix();
            explicitSolution.Add(-1, CalculatePreviousDisplacementVector());
            explicitSolution.Add(0, InitialValues.InitialDisplacementVector);
            explicitAcceleration.Add(0, CalculateInitialAccelerations());
            TimeAtEachStep.Add(-1, -1 * timeStep + InitialValues.InitialTime);
            TimeAtEachStep.Add(0, 0.0);
            double[] nextSolution;
            //Assembler.CalculateEASMatrices();//added EAS
            for (int i = 1; i < timeStepsNumber; i++)
            {
                double time = i * timeStep + InitialValues.InitialTime;

                if (ActivateNonLinearSolution == false)
                {
                    double[] hatRVector = CalculateHatRVector(i);
                    nextSolution = LinearSolver.Solve(hatMassMatrix, hatRVector);
                    Console.WriteLine("Solution for time Step {0} is:", i);
                    VectorOperations.PrintVector(nextSolution);
                }
                else
                {
                    double[] hatRVectorNL = CalculateHatRVectorNL(i);
                    nextSolution = LinearSolver.Solve(hatMassMatrix, hatRVectorNL);
                    //nextSolution = NewtonIterations(hatRVector);
                    Console.WriteLine("Solution for time Step {0} is:", i);
                    VectorOperations.PrintVector(nextSolution);
                }
                explicitSolution.Add(i, nextSolution);
                explicitAcceleration.Add(i, CalculateAccelerations());
                TimeAtEachStep.Add(i, time);
            }
            ExportToFile.ExportExplicitResults(explicitSolution, TimeAtEachStep, 1, 1);
            var tempTimeAtEachStep = TimeAtEachStep;
            TimeAtEachStep = new Dictionary<int, double>();
            int k = 0;
            foreach (var step in tempTimeAtEachStep)
            {
                TimeAtEachStep.Add(k, step.Value);
                k = k + 1;
            }
            
            
            k = 0;
            foreach (var step in explicitSolution)
            {
                displacement.Add(k, step.Value);
                k = k + 1;
            }
        }

        public void PrintExplicitSolution()
        {
            foreach (KeyValuePair<int, double[]> element in explicitSolution)
            {
                int step = element.Key;
                double[] solutionInStep = element.Value;
                Console.WriteLine("Solution for Load Step {0} is:", step);
                VectorOperations.PrintVector(solutionInStep);
            }
        }

        #region Newmark_Method
        private List<double> CalculateIntegrationConstantsNewmark()
        {
            double a = 0.25 * 1.21;
            double delta = 0.5 + 0.1;
            List<double> aConst = new List<double>();
            //{
            //    [0] = 1.0 / (a * Math.Pow(timeStep, 2.0)),
            //    [1] = delta / (a * timeStep),
            //    [2] = 1.0 / (a * timeStep),
            //    [3] = (1.0 / (2.0 * a)) - 1.0,
            //    [4] = (delta / a) - 1.0,
            //    [5] = (timeStep / 2.0) * ((delta / a) - 2.0),
            //    [6] = timeStep * (1.0 - delta),
            //    [7] = delta * timeStep
            //};
            aConst.Add(1.0 / (a * Math.Pow(timeStep, 2.0)));
            aConst.Add(delta / (a * timeStep));
            aConst.Add(1.0 / (a * timeStep));
            aConst.Add((1.0 / (2.0 * a)) - 1.0);
            aConst.Add((delta / a) - 1.0);
            aConst.Add((timeStep / 2.0) * ((delta / a) - 2.0));
            aConst.Add(timeStep * (1.0 - delta));
            aConst.Add(delta * timeStep);
            return aConst;
        }

        private double[,] CalculateHatKMatrixNewmark(List<double> aConstants)
        {
            if (Assembler.ActivateParallelCalculations)
            {
                if (CustomMassMatrix != null)
                {
                    var hutK = MatrixOperations.MatrixAdditionParallelNew(CustomStiffnessMatrix,
                                    MatrixOperations.MatrixAdditionParallelNew(MatrixOperations.ScalarMatrixProductInParallel(aConstants[0], CustomMassMatrix),
                                    MatrixOperations.ScalarMatrixProductInParallel(aConstants[1], CustomDampingMatrix)));
                    return hutK;
                }
                else
                {
                    var hutK = MatrixOperations.MatrixAdditionParallelNew(Assembler.CreateTotalStiffnessMatrix(),
                                    MatrixOperations.MatrixAdditionParallelNew(MatrixOperations.ScalarMatrixProductInParallel(aConstants[0], Assembler.CreateTotalMassMatrix()),
                                    MatrixOperations.ScalarMatrixProductInParallel(aConstants[1], Assembler.CreateTotalDampingMatrix())));
                    return hutK;
                }
            }
            else
            {
                if (CustomMassMatrix != null)
                {
                    var hutK = MatrixOperations.MatrixAddition(CustomStiffnessMatrix,
                                    MatrixOperations.MatrixAddition(MatrixOperations.ScalarMatrixProductNew(aConstants[0], CustomMassMatrix),
                                    MatrixOperations.ScalarMatrixProductNew(aConstants[1], CustomDampingMatrix)));
                    return hutK;
                }
                else
                {
                    var hutK = MatrixOperations.MatrixAddition(Assembler.CreateTotalStiffnessMatrix(),
                                    MatrixOperations.MatrixAddition(MatrixOperations.ScalarMatrixProductNew(aConstants[0], Assembler.CreateTotalMassMatrix()),
                                    MatrixOperations.ScalarMatrixProductNew(aConstants[1], Assembler.CreateTotalDampingMatrix())));
                    return hutK;
                }
            }
        }

        private double[] CalculateHatRVectorNewmark(int i, List<double> aConstants)
        {
            if (Assembler.ActivateParallelCalculations)
            {
                if (CustomMassMatrix != null)
                {
                    double[] currentU = explicitSolution[i - 1];
                    double[] currentdU = explicitVelocity[i - 1];
                    double[] currentddU = explicitAcceleration[i - 1];

                    double[] a0U = VectorOperations.VectorScalarProductParallel(currentU, aConstants[0]);
                    double[] a2dU = VectorOperations.VectorScalarProductParallel(currentdU, aConstants[2]);
                    double[] a3ddU = VectorOperations.VectorScalarProductParallel(currentddU, aConstants[3]);
                    double[] a1U = VectorOperations.VectorScalarProductParallel(currentU, aConstants[1]);
                    double[] a4dU = VectorOperations.VectorScalarProductParallel(currentdU, aConstants[4]);
                    double[] a5ddU = VectorOperations.VectorScalarProductParallel(currentddU, aConstants[5]);

                    var hatR = VectorOperations.VectorVectorAdditionParallel(ExternalForcesVector,
                                    VectorOperations.VectorVectorAdditionParallel(VectorOperations.MatrixVectorProduct(CustomMassMatrix,
                                    VectorOperations.VectorVectorAdditionParallel(a0U, VectorOperations.VectorVectorAdditionParallel(a2dU, a3ddU))),
                                    VectorOperations.MatrixVectorProduct(CustomDampingMatrix,
                                    VectorOperations.VectorVectorAdditionParallel(a1U, VectorOperations.VectorVectorAdditionParallel(a4dU, a5ddU)))));
                    return hatR;
                }
                else
                {
                    double[] currentU = explicitSolution[i - 1];
                    double[] currentdU = explicitVelocity[i - 1];
                    double[] currentddU = explicitAcceleration[i - 1];

                    double[] a0U = VectorOperations.VectorScalarProductParallel(currentU, aConstants[0]);
                    double[] a2dU = VectorOperations.VectorScalarProductParallel(currentdU, aConstants[2]);
                    double[] a3ddU = VectorOperations.VectorScalarProductParallel(currentddU, aConstants[3]);
                    double[] a1U = VectorOperations.VectorScalarProductParallel(currentU, aConstants[1]);
                    double[] a4dU = VectorOperations.VectorScalarProductParallel(currentdU, aConstants[4]);
                    double[] a5ddU = VectorOperations.VectorScalarProductParallel(currentddU, aConstants[5]);

                    var hatR = VectorOperations.VectorVectorAdditionParallel(ExternalForcesVector,
                                    VectorOperations.VectorVectorAdditionParallel(VectorOperations.MatrixVectorProduct(Assembler.CreateTotalMassMatrix(),
                                    VectorOperations.VectorVectorAdditionParallel(a0U, VectorOperations.VectorVectorAdditionParallel(a2dU, a3ddU))),
                                    VectorOperations.MatrixVectorProduct(Assembler.CreateTotalDampingMatrix(),
                                    VectorOperations.VectorVectorAdditionParallel(a1U, VectorOperations.VectorVectorAdditionParallel(a4dU, a5ddU)))));
                    return hatR;
                }
            }
            else
            {
                if (CustomMassMatrix != null)
                {
                    double[] currentU = explicitSolution[i - 1];
                    double[] currentdU = explicitVelocity[i - 1];
                    double[] currentddU = explicitAcceleration[i - 1];

                    double[] a0U = VectorOperations.VectorScalarProductNew(currentU, aConstants[0]);
                    double[] a2dU = VectorOperations.VectorScalarProductNew(currentdU, aConstants[2]);
                    double[] a3ddU = VectorOperations.VectorScalarProductNew(currentddU, aConstants[3]);
                    double[] a1U = VectorOperations.VectorScalarProductNew(currentU, aConstants[1]);
                    double[] a4dU = VectorOperations.VectorScalarProductNew(currentdU, aConstants[4]);
                    double[] a5ddU = VectorOperations.VectorScalarProductNew(currentddU, aConstants[5]);

                    var hatR = VectorOperations.VectorVectorAddition(ExternalForcesVector,
                                    VectorOperations.VectorVectorAddition(VectorOperations.MatrixVectorProduct(CustomMassMatrix,
                                    VectorOperations.VectorVectorAddition(a0U, VectorOperations.VectorVectorAddition(a2dU, a3ddU))),
                                    VectorOperations.MatrixVectorProduct(CustomDampingMatrix,
                                    VectorOperations.VectorVectorAddition(a1U, VectorOperations.VectorVectorAddition(a4dU, a5ddU)))));
                    return hatR;
                }
                else
                {
                    double[] currentU = explicitSolution[i - 1];
                    double[] currentdU = explicitVelocity[i - 1];
                    double[] currentddU = explicitAcceleration[i - 1];

                    double[] a0U = VectorOperations.VectorScalarProductNew(currentU, aConstants[0]);
                    double[] a2dU = VectorOperations.VectorScalarProductNew(currentdU, aConstants[2]);
                    double[] a3ddU = VectorOperations.VectorScalarProductNew(currentddU, aConstants[3]);
                    double[] a1U = VectorOperations.VectorScalarProductNew(currentU, aConstants[1]);
                    double[] a4dU = VectorOperations.VectorScalarProductNew(currentdU, aConstants[4]);
                    double[] a5ddU = VectorOperations.VectorScalarProductNew(currentddU, aConstants[5]);

                    var hatR = VectorOperations.VectorVectorAddition(ExternalForcesVector,
                                    VectorOperations.VectorVectorAddition(VectorOperations.MatrixVectorProduct(Assembler.CreateTotalMassMatrix(),
                                    VectorOperations.VectorVectorAddition(a0U, VectorOperations.VectorVectorAddition(a2dU, a3ddU))),
                                    VectorOperations.MatrixVectorProduct(Assembler.CreateTotalDampingMatrix(),
                                    VectorOperations.VectorVectorAddition(a1U, VectorOperations.VectorVectorAddition(a4dU, a5ddU)))));
                    return hatR;
                }
            }
        }

        private double[] CalculateAccelerationNewmark(int step, List<double> aConstants)
        {
            double[] dUt = explicitVelocity[step - 1];
            double[] ut = explicitSolution[step - 1];
            double[] utplusDt = explicitSolution[step];
            double[] ddUt = explicitAcceleration[step - 1];

            double[] part1 = VectorOperations.VectorScalarProductNew(
                                VectorOperations.VectorVectorSubtraction(utplusDt, ut), aConstants[0]);
            double[] part2 = VectorOperations.VectorScalarProductNew(dUt, aConstants[2]);
            double[] part3 = VectorOperations.VectorScalarProductNew(ddUt, aConstants[3]);

            double[] ddUtplusDt = VectorOperations.VectorVectorSubtraction(part1,
                                    VectorOperations.VectorVectorAddition(part2, part3));
            return ddUtplusDt;
        }

        private double[] CalculateVelocityNewmark(int step, List<double> aConstants)
        {
            double[] dUt = explicitVelocity[step - 1];
            double[] ddUt = explicitAcceleration[step - 1];
            double[] ddUtplusDt = explicitAcceleration[step];

            double[] part1 = VectorOperations.VectorScalarProductNew(ddUt, aConstants[6]);
            double[] part2 = VectorOperations.VectorScalarProductNew(ddUtplusDt, aConstants[7]);

            double[] dUtplusDt = VectorOperations.VectorVectorAddition(dUt,
                                    VectorOperations.VectorVectorAddition(part1, part2));
            return dUtplusDt;
        }

        public void SolveNewmark()
        {
            List<double> aConstants = CalculateIntegrationConstantsNewmark();
            double[,] hatStiffnessMatrixNewmark = CalculateHatKMatrixNewmark(aConstants);
            explicitSolution.Add(0, InitialValues.InitialDisplacementVector);
            explicitVelocity.Add(0, InitialValues.InitialVelocityVector);
            explicitAcceleration.Add(0, CalculateInitialAccelerationsNewmark());
            TimeAtEachStep.Add(0, 0.0);
            double[] nextSolution;
            //Assembler.InitializeEASParameters();//Added EAS
            //Assembler.CalculateEASMatrices();//added EAS
            for (int i = 1; i <= timeStepsNumber; i++)
            {
                double time = i * timeStep + InitialValues.InitialTime;

                if (ActivateNonLinearSolution == false)
                {
                    double[] hatRVectorNewmark = CalculateHatRVectorNewmark(i, aConstants);
                    nextSolution = LinearSolver.Solve(hatStiffnessMatrixNewmark, hatRVectorNewmark);
                    //Console.WriteLine("Solution for Load Step {0} is:", i);
                    VectorOperations.PrintVector(nextSolution);
                }
                else
                {
                    //Assembler.StoreFinalStepDisplacementVector(explicitSolution.Values.Last());//Added EAS
                    nextSolution = NewtonIterationsNewmark(ExternalForcesVector, i, aConstants);
                    //Console.WriteLine("Solution for Load Step {0} is:", i);
                    //VectorOperations.PrintVector(nextSolution);
                }
                explicitSolution.Add(i, nextSolution);
                double[] fullDynamicSol = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(nextSolution, Assembler.BoundedDOFsVector);
                VectorOperations.PrintVectorToFile(fullDynamicSol, @"C:\Users\Public\Documents\DynamicSol" + i.ToString() + ".dat");                
                explicitAcceleration.Add(i, CalculateAccelerationNewmark(i, aConstants));
                explicitVelocity.Add(i, CalculateVelocityNewmark(i, aConstants));

                double[] fullVelocitySol = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(explicitVelocity[i], Assembler.BoundedDOFsVector);
                double[] fullAccelerationSol = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(explicitAcceleration[i], Assembler.BoundedDOFsVector);
                VectorOperations.PrintVectorToFile(fullVelocitySol, @"C:\Users\Public\Documents\Velocity" + i.ToString() + ".dat");
                VectorOperations.PrintVectorToFile(fullAccelerationSol, @"C:\Users\Public\Documents\Acceleration" + i.ToString() + ".dat");
                TimeAtEachStep.Add(i, time);
            }
            ExportToFile.ExportExplicitResults(explicitSolution, TimeAtEachStep, 1, 1);
            //ShowToGUI.ShowResults(explicitSolution, TimeAtEachStep, 1, 1);
        }

        private double[] CalculateHatRVectorNewmarkNL(int i, List<double> aConstants, double[] previousIterationSolution)
        {
            if (Assembler.ActivateParallelCalculations)
            {
                if (CustomMassMatrix != null)
                {
                    double[] a0U = VectorOperations.VectorScalarProductParallel(
                    VectorOperations.VectorVectorSubtractionParallel(explicitSolution[i - 1], previousIterationSolution), aConstants[0]);
                    double[] a2dU = VectorOperations.VectorScalarProductParallel(explicitVelocity[i - 1], aConstants[2]);
                    double[] a3ddU = VectorOperations.VectorScalarProductParallel(explicitAcceleration[i - 1], aConstants[3]);
                    double[] a1U = VectorOperations.VectorScalarProductParallel(explicitSolution[i - 1], aConstants[1]);
                    double[] a4dU = VectorOperations.VectorScalarProductParallel(explicitVelocity[i - 1], aConstants[4]);
                    double[] a5ddU = VectorOperations.VectorScalarProductParallel(explicitAcceleration[i - 1], aConstants[5]);

                    var hatR = VectorOperations.VectorVectorAdditionParallel(ExternalForcesVector,
                                VectorOperations.VectorVectorAdditionParallel(VectorOperations.MatrixVectorProduct(CustomMassMatrix,
                                VectorOperations.VectorVectorAdditionParallel(a0U, VectorOperations.VectorVectorAdditionParallel(a2dU, a3ddU))),
                                VectorOperations.MatrixVectorProduct(CustomDampingMatrix,
                                VectorOperations.VectorVectorAdditionParallel(a1U, VectorOperations.VectorVectorAdditionParallel(a4dU, a5ddU)))));
                    return hatR;
                }
                else
                {
                    double[] a0U = VectorOperations.VectorScalarProductParallel(
                       VectorOperations.VectorVectorSubtractionParallel(explicitSolution[i - 1], previousIterationSolution), aConstants[0]);
                    double[] a2dU = VectorOperations.VectorScalarProductParallel(explicitVelocity[i - 1], aConstants[2]);
                    double[] a3ddU = VectorOperations.VectorScalarProductParallel(explicitAcceleration[i - 1], aConstants[3]);
                    double[] a1U = VectorOperations.VectorScalarProductParallel(explicitSolution[i - 1], aConstants[1]);
                    double[] a4dU = VectorOperations.VectorScalarProductParallel(explicitVelocity[i - 1], aConstants[4]);
                    double[] a5ddU = VectorOperations.VectorScalarProductParallel(explicitAcceleration[i - 1], aConstants[5]);

                    var hatR = VectorOperations.VectorVectorAdditionParallel(ExternalForcesVector,
                                VectorOperations.VectorVectorAdditionParallel(VectorOperations.MatrixVectorProduct(Assembler.CreateTotalMassMatrix(),
                                VectorOperations.VectorVectorAdditionParallel(a0U, VectorOperations.VectorVectorAdditionParallel(a2dU, a3ddU))),
                                VectorOperations.MatrixVectorProduct(Assembler.CreateTotalDampingMatrix(),
                                VectorOperations.VectorVectorAdditionParallel(a1U, VectorOperations.VectorVectorAdditionParallel(a4dU, a5ddU)))));
                    return hatR;
                }
            }
            else
            {
                if (CustomMassMatrix != null)
                {
                    double[] a0U = VectorOperations.VectorScalarProductNew(
                    VectorOperations.VectorVectorSubtraction(explicitSolution[i - 1], previousIterationSolution), aConstants[0]);
                    double[] a2dU = VectorOperations.VectorScalarProductNew(explicitVelocity[i - 1], aConstants[2]);
                    double[] a3ddU = VectorOperations.VectorScalarProductNew(explicitAcceleration[i - 1], aConstants[3]);
                    double[] a1U = VectorOperations.VectorScalarProductNew(explicitSolution[i - 1], aConstants[1]);
                    double[] a4dU = VectorOperations.VectorScalarProductNew(explicitVelocity[i - 1], aConstants[4]);
                    double[] a5ddU = VectorOperations.VectorScalarProductNew(explicitAcceleration[i - 1], aConstants[5]);

                    var hatR = VectorOperations.VectorVectorAddition(ExternalForcesVector,
                                VectorOperations.VectorVectorAddition(VectorOperations.MatrixVectorProduct(CustomMassMatrix,
                                VectorOperations.VectorVectorAddition(a0U, VectorOperations.VectorVectorAddition(a2dU, a3ddU))),
                                VectorOperations.MatrixVectorProduct(CustomDampingMatrix,
                                VectorOperations.VectorVectorAddition(a1U, VectorOperations.VectorVectorAddition(a4dU, a5ddU)))));
                    return hatR;
                }
                else
                {
                    double[] a0U = VectorOperations.VectorScalarProductNew(
                       VectorOperations.VectorVectorSubtraction(explicitSolution[i - 1], previousIterationSolution), aConstants[0]);
                    double[] a2dU = VectorOperations.VectorScalarProductNew(explicitVelocity[i - 1], aConstants[2]);
                    double[] a3ddU = VectorOperations.VectorScalarProductNew(explicitAcceleration[i - 1], aConstants[3]);
                    double[] a1U = VectorOperations.VectorScalarProductNew(explicitSolution[i - 1], aConstants[1]);
                    double[] a4dU = VectorOperations.VectorScalarProductNew(explicitVelocity[i - 1], aConstants[4]);
                    double[] a5ddU = VectorOperations.VectorScalarProductNew(explicitAcceleration[i - 1], aConstants[5]);

                    var hatR = VectorOperations.VectorVectorAddition(ExternalForcesVector,
                                VectorOperations.VectorVectorAddition(VectorOperations.MatrixVectorProduct(Assembler.CreateTotalMassMatrix(),
                                VectorOperations.VectorVectorAddition(a0U, VectorOperations.VectorVectorAddition(a2dU, a3ddU))),
                                VectorOperations.MatrixVectorProduct(Assembler.CreateTotalDampingMatrix(),
                                VectorOperations.VectorVectorAddition(a1U, VectorOperations.VectorVectorAddition(a4dU, a5ddU)))));
                    return hatR;
                }
            }
        }

        private double[] CalculateInitialAccelerationsNewmark() //Bathe page 771
        {
            if (CustomStiffnessMatrix != null) return InitialValues.InitialAccelerationVector;
            int step = explicitSolution.Count - 1;
            Assembler.UpdateDisplacements(explicitSolution[step]);
            double[,] stiffness = Assembler.CreateTotalStiffnessMatrix();
            double[,] mass = Assembler.CreateTotalMassMatrix();

            double[] Ku = VectorOperations.MatrixVectorProduct(stiffness, explicitSolution[step]);
            double[] RHS = VectorOperations.VectorVectorSubtraction(ExternalForcesVector, Ku);

            double[] acceleration = LinearSolver.Solve(mass, RHS);

            return acceleration;
        }

        private double[] NewtonIterationsNewmark(double[] forceVector, int stepNumber, List<double> aConstants)
        {

            lambda = 1.0 / numberOfLoadSteps;
            //double[] incrementDf = VectorOperations.VectorScalarProductNew(forceVector, lambda);
            double[] solutionVector = new double[forceVector.Length];
            double[] deltaU = new double[solutionVector.Length];
            double[] internalForcesTotalVector;
            double[] residual;
            double residualNorm;
            double[] hatR;

            solutionVector = explicitSolution.Values.Last();

            Assembler.UpdateDisplacements(solutionVector);
            Assembler.UpdateEASParameters(solutionVector);//added EAS

            //Assembler.UpdateAccelerations(explicitAcceleration.Values.Last());
            hatR = CalculateHatRVectorNewmarkNL(stepNumber, aConstants, solutionVector);
            internalForcesTotalVector = Assembler.CreateTotalInternalForcesVector();
            residual = VectorOperations.VectorVectorSubtraction(hatR, internalForcesTotalVector);
            //residual = VectorOperations.VectorVectorSubtraction(forceVector, Assembler.CreateTotalInternalForcesVector());
            residualNorm = VectorOperations.VectorNorm2(residual);
            int iteration = 0;
            Array.Clear(deltaU, 0, deltaU.Length);

            for (int i = 0; i < maxIterations; i++)
            {
                double[,] tangentMatrix = CalculateHatKMatrixNewmark(aConstants);
                deltaU = LinearSolver.Solve(tangentMatrix, residual);
                solutionVector = VectorOperations.VectorVectorAddition(solutionVector, deltaU);
                Assembler.UpdateDisplacements(solutionVector);
                Assembler.UpdateEASParameters(solutionVector);//added EAS
                //Assembler.UpdateAccelerations(CalculateAccelerations());
                hatR = CalculateHatRVectorNewmarkNL(stepNumber, aConstants, solutionVector);
                internalForcesTotalVector = Assembler.CreateTotalInternalForcesVector();
                residual = VectorOperations.VectorVectorSubtraction(hatR, internalForcesTotalVector);
                residualNorm = VectorOperations.VectorNorm2(residual);
                if (residualNorm < 1.0e-4)
                {
                    break;
                }
                iteration = iteration + 1;
                //Assembler.CalculateEASMatrices();//added EAS
            }
            //Console.WriteLine(iteration);
            if (iteration >= maxIterations && residualNorm > tolerance) throw new Exception("Newton-Raphson: Solution not converged at current iterations");

            return solutionVector;
        }
        #endregion
        public Tuple<Dictionary<int, double[]>, Dictionary<int, double>> GetResults()
        {
            return new Tuple<Dictionary<int, double[]>, Dictionary<int, double>>(explicitSolution, TimeAtEachStep);
        }
    }
}
