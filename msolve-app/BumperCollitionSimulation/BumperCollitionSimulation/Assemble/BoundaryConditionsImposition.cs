using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BumperCollitionSimulation
{
    class BoundaryConditionsImposition
    {
        public static double[,] ReducedTotalStiff(double[,] totalstiff, int[] boundaryDof)
        {
            int rows = totalstiff.GetLength(0);
            int cols = totalstiff.GetLength(1);
            int dofValues = boundaryDof.GetLength(0);
            int newDim = rows - dofValues;
            double[,] reducedMatrix = new double[newDim, newDim];
            int m;
            int n;
            m = 0;
            for (int i = 0; i < rows; i++)
            {

                if (boundaryDof.Contains(i + 1))   //i+1 because C# is zero based//
                    continue;
                else
                    n = 0;
                for (int j = 0; j < cols; j++)
                {

                    if (boundaryDof.Contains(j + 1))
                        continue;
                    else
                    {
                        reducedMatrix[m, n] = totalstiff[i, j];
                        n = n + 1;
                    }
                }
                m = m + 1;
            }
            return reducedMatrix;
        }
        public static Tuple<double[,], double[,], double[,], double[,]> MMCPCGReducedTotalStiffMatrices(double[,] totalstiff, int[] boundaryDof,
            int[] contactDof, int[] NoContactDof)
        {
            int rows = totalstiff.GetLength(0);
            int cols = totalstiff.GetLength(1);
            int dofValuesK = NoContactDof.GetLength(0);
            int dofValuesC = contactDof.GetLength(0);
            int rowsK = new int();
            int rowsC = new int();
            for (int i = 0; i < dofValuesC; i++)
            {
                if (!boundaryDof.Contains(contactDof[i]))
                {
                    rowsC += 1;
                }
            }
            for (int i = 0; i < dofValuesK; i++)
            {
                if (!boundaryDof.Contains(NoContactDof[i]))
                {
                    rowsK += 1;
                }
            }
            double[,] K = new double[rowsK, rowsK];
            double[,] B = new double[rowsK, rowsC];
            double[,] C = new double[rowsC, rowsC];
            double[,] TotalStifnessMatrix = new double[rowsK + rowsC, rowsK + rowsC];

            int m1 = 0;
            int m2 = 0;
            int m3 = 0;
            bool crow = false;
            bool krow = false;
            bool brow = false;
            for (int i = 0; i < rows; i++)
            {

                if (boundaryDof.Contains(i + 1))
                    continue;
                else
                {
                    int n1 = 0;
                    int n2 = 0;
                    int n3 = 0;
                    for (int j = 0; j < cols; j++)
                    {

                        if (boundaryDof.Contains(j + 1))
                            continue;
                        else
                        {
                            if (NoContactDof.Contains(i + 1) && NoContactDof.Contains(j + 1)) 
                            {
                                K[m3, n3] = totalstiff[i, j];
                                n3 += 1;
                                krow = true;
                            }
                            else if (NoContactDof.Contains(i + 1) && contactDof.Contains(j + 1))
                            {
                                B[m2, n2] = totalstiff[i, j];
                                n2 += 1;
                                brow = true;
                            }
                            else if (contactDof.Contains(i + 1) && contactDof.Contains(j + 1))
                            {
                                C[m1, n1] = totalstiff[i, j];
                                n1 += 1;
                                crow = true;
                            }
                        }
                    }
                }
                if (crow)
                {
                    m1 += 1;
                    crow = false;
                }
                if (brow)
                {
                    m2 += 1;
                    brow = false;
                }
                if(krow)
                {
                    m3 += 1;
                    krow = false;
                }
            }
            for(int i =0; i< rowsK + rowsC; i++)
            {
                for(int j = 0; j < rowsK + rowsC; j++)
                {
                    if(i< rowsK && j < rowsK)
                    {
                        TotalStifnessMatrix[i, j] = K[i, j];
                    }
                    else if (i < rowsK && j >= rowsK)
                    {
                        TotalStifnessMatrix[i, j] = B[i, j - rowsK];
                    }
                    else if(i >= rowsK && j < rowsK)
                    {
                        TotalStifnessMatrix[i, j] = B[j, i - rowsK];
                    }
                    else if (i >= rowsK && j >= rowsK)
                    {
                        TotalStifnessMatrix[i, j] = C[i - rowsK, j - rowsK];
                    }
                }
            }
            return new Tuple<double[,], double[,], double[,], double[,]>(TotalStifnessMatrix, K, B, C);
        }

        public static double[] ReducedVector(double[] vectorToReduce, int[] boundaryDOF)
        {
            int reducedVectorLength = vectorToReduce.Length - boundaryDOF.Length;
            double[] reducedVector = new double[reducedVectorLength];
            int newRow = 0;
            for (int oldRow = 0; oldRow < vectorToReduce.Length; oldRow++)
            {
                if (boundaryDOF.Contains(oldRow + 1))
                {
                    continue;
                }
                else
                {
                    reducedVector[newRow] = vectorToReduce[oldRow];
                    newRow = newRow + 1;
                }
            }
            return reducedVector;
        }

        public static Tuple<double[], double[], double[]>  MMCPCGReducedVector(double[] vectorToReduce, int[] boundaryDOF, int[] contactDof, int[] NoContactDof)
        {
            int contactDofValues = contactDof.GetLength(0);
            int noContactDofValues = NoContactDof.GetLength(0);

            int reducedVectorLengthM = new int();
            int reducedVectorLengthC = new int();

            for (int i = 0; i< contactDofValues; i++)
            {
                if (!boundaryDOF.Contains(contactDof[i]))
                {
                    reducedVectorLengthC += 1;
                }
            }
            for (int i = 0; i < noContactDofValues; i++)
            {
                if (!boundaryDOF.Contains(NoContactDof[i]))
                {
                    reducedVectorLengthM += 1;
                }
            }
            double[] reducedVectorM = new double[reducedVectorLengthM];
            double[] reducedVectorC = new double[reducedVectorLengthC];
            double[] reducedTotalVector = new double[reducedVectorLengthM + reducedVectorLengthC];
            int newRowC = 0;
            int newRowM = 0;
            for (int oldRow = 0; oldRow < vectorToReduce.Length; oldRow++)
            {
                if (boundaryDOF.Contains(oldRow + 1))
                {
                    continue;
                }
                else
                {
                    if (NoContactDof.Contains(oldRow + 1))
                    {
                        reducedVectorM[newRowM] = vectorToReduce[oldRow];
                        newRowM += 1;
                    }
                    else
                    {
                        reducedVectorC[newRowC] = vectorToReduce[oldRow];
                        newRowC += 1;
                    }
                }
            }
            for(int i =0; i< reducedVectorLengthM + reducedVectorLengthC; i++)
            {
                if(i< reducedVectorLengthM)
                {
                    reducedTotalVector[i] = reducedVectorM[i];
                }
                else
                {
                    reducedTotalVector[i] = reducedVectorC[i - reducedVectorLengthM];
                }
            }

            return new Tuple<double[], double[], double[]>(reducedVectorM, reducedVectorC, reducedTotalVector);
        }
        public static double[] CreateFullVectorFromReducedVector(double[] reducedVector, int[] boundaryDOF)
        {
            int fullVectorLength = reducedVector.Length + boundaryDOF.Length;
            double[] fullVector = new double[fullVectorLength];
            int reducedVectorRow = 0;
            for (int fullVectorRow = 0; fullVectorRow < fullVectorLength; fullVectorRow++)
            {
                if (boundaryDOF.Contains(fullVectorRow + 1))
                {
                    fullVector[fullVectorRow] = 0;
                }
                else
                {
                    fullVector[fullVectorRow] = reducedVector[reducedVectorRow];
                    reducedVectorRow = reducedVectorRow + 1;
                }
            }
            return fullVector;
        }

        public static Dictionary<int, double[]> CreateFullVectorListFromReduced(Dictionary<int, double[]> reducedVectorList, int[] boundaryDOF)
        {
            Dictionary<int, double[]> fullVectorList = new Dictionary<int, double[]>();
            foreach (KeyValuePair<int, double[]> reducedVector in reducedVectorList)
            {
                double[] fullVector = CreateFullVectorFromReducedVector(reducedVector.Value, boundaryDOF);
                fullVectorList.Add(reducedVector.Key, fullVector);
            }
            return fullVectorList;
        }
        public static double[] MMCPCGCreateFullVectorFromReducedVector(double[] reducedVector, int[] boundaryDOF, int[] contactDof, int[] NoContactDof)
        {
            int fullVectorLength = reducedVector.Length + boundaryDOF.Length;
            double[] fullVector = new double[fullVectorLength];
            int countM = 0;
            foreach(var dof in NoContactDof)
            {
                if (!boundaryDOF.Contains(dof))
                {
                    countM += 1;
                }
            }
            int reducedVectorRow1 = 0;
            int reducedVectorRow2 = countM;
            for (int fullVectorRow = 0; fullVectorRow < fullVectorLength; fullVectorRow++)
            {
                if (boundaryDOF.Contains(fullVectorRow + 1))
                {
                    fullVector[fullVectorRow] = 0;
                }
                else if(NoContactDof.Contains(fullVectorRow + 1))
                {
                    fullVector[fullVectorRow] = reducedVector[reducedVectorRow1];
                    reducedVectorRow1 += 1;
                }
                else if (contactDof.Contains(fullVectorRow + 1))
                {
                    fullVector[fullVectorRow] = reducedVector[reducedVectorRow2];
                    reducedVectorRow2 += 1;
                }
            }
            return fullVector;
        }
        public static double[] MMCPCGRearrangeReducedVector(double[] reducedVector, int[] boundaryDOF, int[] contactDof, int[] NoContactDof)
        {
            int rearrangedReducedVectorLength = reducedVector.Length;
            int fullVectorLength = reducedVector.Length + boundaryDOF.Length;
            double [] rearrangedFullVector = new double[fullVectorLength];
            double[] rearrangedReducedVector = new double[rearrangedReducedVectorLength];
            int countM = 0;
            foreach (var dof in NoContactDof)
            {
                if (!boundaryDOF.Contains(dof))
                {
                    countM += 1;
                }
            }
            int reducedVectorRow1 = 0;
            int reducedVectorRow2 = countM;
            for (int fullVectorRow = 0; fullVectorRow < fullVectorLength; fullVectorRow++)
            {
                if(boundaryDOF.Contains(fullVectorRow + 1))
                {
                    rearrangedFullVector[fullVectorRow] = 0.0;
                    continue;
                }
                else
                {
                    if (NoContactDof.Contains(fullVectorRow + 1))
                    {
                        rearrangedFullVector[fullVectorRow] = reducedVector[reducedVectorRow1];
                        reducedVectorRow1 += 1;
                    }
                    else if (contactDof.Contains(fullVectorRow + 1))
                    {
                        rearrangedFullVector[fullVectorRow] = reducedVector[reducedVectorRow2];
                        reducedVectorRow2 += 1;
                    }
                }

            }
            rearrangedReducedVector = ReducedVector(rearrangedFullVector, boundaryDOF);
            return rearrangedReducedVector;
        }
        public static Tuple<double[], double[]> MMCPCGSeperateReducedVectorMatrices(double[] reducedTotalVector, int[] boundaryDOF, int[] contactDof, int[] NoContactDof)
        {
            int contactDofValues = contactDof.GetLength(0);
            int noContactDofValues = NoContactDof.GetLength(0);

            int reducedVectorLengthM = new int();
            int reducedVectorLengthC = new int();

            for (int i = 0; i < contactDofValues; i++)
            {
                if (!boundaryDOF.Contains(contactDof[i]))
                {
                    reducedVectorLengthC += 1;
                }
            }
            for (int i = 0; i < noContactDofValues; i++)
            {
                if (!boundaryDOF.Contains(NoContactDof[i]))
                {
                    reducedVectorLengthM += 1;
                }
            }
            double[] reducedVectorM = new double[reducedVectorLengthM];
            double[] reducedVectorC = new double[reducedVectorLengthC];
            int newRowC = 0;
            int newRowM = 0;
            for (int oldRow = 0; oldRow < reducedTotalVector.GetLength(0); oldRow++)
            {
                if (oldRow < reducedVectorLengthM)
                {
                    reducedVectorM[newRowM] = reducedTotalVector[oldRow];
                    newRowM += 1;
                }
                else
                {
                    reducedVectorC[newRowC] = reducedTotalVector[oldRow];
                    newRowC += 1;
                }
            }

            return new Tuple<double[], double[]>(reducedVectorM, reducedVectorC);
        }
    }
}

