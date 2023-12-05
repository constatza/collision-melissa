using System;
using System.IO;
using System.Threading.Tasks;

namespace BumperCollitionSimulation
{
    static class MatrixOperations
    {
        public static double[,] TempVariable;
        public static double[,] TempVariable2;
        public static double[,] TempVariable3;
        public static double[,] TempVariable4;
        public static bool ParallelCalculations { get; set; } = false;
        public static void PrintMatrix(double[,] matrix)
        {
            int matrixRows = matrix.GetLength(0);
            int matrixCols = matrix.GetLength(1);
            for (int row = 0; row < matrixRows; row++)
            {
                for (int col = 0; col < matrixCols; col++)
                {
                    Console.Write(String.Format("{0}\t", matrix[row, col]));
                }

                Console.WriteLine();
            }

        }

        public static double[,] CreateRandomMatrix(int rows, int columns)
        {
            double[,] randomMatrix = new double[rows, columns];
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < columns; j++)
                {
                    randomMatrix[i, j] = new Random().NextDouble();
                }
            }
            return randomMatrix;
        }

        public static double[,] CreateDiagonalMatrix(int dimension, double diagonalNumber)
        {
            double[,] diagMatrix = new double[dimension, dimension];
            for (int i = 0; i < dimension; i++)
            {
                diagMatrix[i, i] = diagonalNumber;
            }
            return diagMatrix;
        }

        public static double[,] Transpose(double[,] matrix)
        {
            int matrixRows = matrix.GetLength(0);
            int matrixCols = matrix.GetLength(1);
            double[,] m = new double[matrixCols, matrixRows];
            for (int row = 0; row < matrixRows; row++)
            {
                for (int col = 0; col < matrixCols; col++)
                {
                    m[col, row] = matrix[row, col];
                }
            }
            return m;
        }
        public static double[,] TransposeSquareMatrix(double[,] matrix)
        {
            for (int row = 0; row <= matrix.GetLength(0) - 2; row++)
            {
                for (int col = row + 1; col <= matrix.GetLength(0) - 1; col++)
                {
                    double local = matrix[row, col];
                    matrix[row, col] = matrix[col, row];
                    matrix[col, row] = local;
                }
            }
            return matrix;
        }
        public static double[,] MatrixProduct(double[,] matrix1, double[,] matrix2)
        {
            
            if (ParallelCalculations == true)
            {
                TempVariable2 = matrix1;
                TempVariable3 = matrix2;
                TempVariable4 = new double[matrix1.GetLength(0), matrix2.GetLength(1)];
                int rowsForEachThread = 400;
                Task first = Task.Run(() => ParallelMatrixProductCalculations(0, rowsForEachThread));
                Task second = Task.Run(() => ParallelMatrixProductCalculations(rowsForEachThread, rowsForEachThread * 2));
                Task third = Task.Run(() => ParallelMatrixProductCalculations(rowsForEachThread * 2, rowsForEachThread * 3));
                Task fourth = Task.Run(() => ParallelMatrixProductCalculations(rowsForEachThread * 3, rowsForEachThread * 4));
                Task fifth = Task.Run(() => ParallelMatrixProductCalculations(rowsForEachThread * 4, rowsForEachThread * 5));
                Task.WaitAll(first, second, third, fourth, fifth);
                //Task.WaitAll(first);
                return TempVariable4;
            }
            else
            {
                int matrix1rows = matrix1.GetLength(0);
                int matrix2cols = matrix2.GetLength(1);
                double[,] productMatrix = new double[matrix1rows, matrix2cols];
                for (int i = 0; i < matrix1rows; i++)
                {
                    for (int j = 0; j < matrix2cols; j++)
                    {
                        double sum = 0;
                        for (int k = 0; k < matrix2.GetLength(0); k++)
                        {
                            sum = sum + matrix1[i, k] * matrix2[k, j];
                        }
                        productMatrix[i, j] = sum;
                    }
                }
                return productMatrix;
            }
        }

        private static void ParallelMatrixProductCalculations(int min, int max)
        {
            for (int i = min; i < max; i++)
            {
                for (int j = 0; j < TempVariable3.GetLength(1); j++)
                {
                    double sum = 0;
                    for (int k = 0; k < TempVariable3.GetLength(0); k++)
                    {
                        sum = sum + TempVariable2[i, k] * TempVariable3[k, j];
                    }
                    TempVariable4[i, j] = sum;
                }
            }
        }

        public static double[,] MatrixProductParallel(double[,] matrix1, double[,] matrix2)
        {
            int matrix1rows = matrix1.GetLength(0);
            int matrix2cols = matrix2.GetLength(1);
            double[,] productMatrix = new double[matrix1rows, matrix2cols];
            Parallel.For(0, matrix1rows, i =>
            {
                for (int j = 0; j < matrix2cols; j++)
                {
                    double sum = 0;
                    for (int k = 0; k < matrix2.GetLength(0); k++)
                    {
                        sum = sum + matrix1[i, k] * matrix2[k, j];
                    }
                    productMatrix[i, j] = sum;
                }
            });
            return productMatrix;
        }

        public static double[,] MatrixAddition(double[,] matrix1, double[,] matrix2)
        {
            int threads, threadsAsync;
            threads = 4;
            //ThreadPool.GetAvailableThreads(out threads, out threadsAsync);
            if (matrix1.GetLength(0) != matrix2.GetLength(0) || matrix1.GetLength(1) != matrix2.GetLength(1))
            {
                throw new Exception("Not equally sized matrices");
            }
            if (ParallelCalculations == true)
            {
                TempVariable = (double[,])matrix1.Clone();

                int totalRows = matrix1.GetLength(0);
                int rowsForEachThread = 400;//totalRows / threads;
                Task[] tasks = new Task[threads];
                int k = 0;
                //for (int i = 0; i < threads; i++)
                //{
                //    tasks[i] = Task.Run(() => MatrixAdditionParallel2Calculations(matrix2, i * rowsForEachThread, i * rowsForEachThread + rowsForEachThread));
                //    //k = k + rowsForEachThread;
                //}

                Task first = Task.Run(() => MatrixAdditionParallel2Calculations(matrix2, 0, rowsForEachThread));
                Task second = Task.Run(() => MatrixAdditionParallel2Calculations(matrix2, rowsForEachThread, rowsForEachThread*2));
                Task third = Task.Run(() => MatrixAdditionParallel2Calculations(matrix2, rowsForEachThread*2, rowsForEachThread*3));
                Task fourth = Task.Run(() => MatrixAdditionParallel2Calculations(matrix2, rowsForEachThread*3, rowsForEachThread*4));
                Task fifth = Task.Run(() => MatrixAdditionParallel2Calculations(matrix2, rowsForEachThread * 4, rowsForEachThread * 5));
                Task.WaitAll(first, second, third, fourth, fifth);
                //Task.WaitAll(tasks);
                return TempVariable;
            }
            else
            {
                int matrixrows = matrix1.GetLength(0);
                int matrixcols = matrix1.GetLength(1);
                double[,] sumMatrix = new double[matrixrows, matrixcols];

                for (int row = 0; row < matrixrows; row++)
                {
                    for (int col = 0; col < matrixcols; col++)
                    {
                        sumMatrix[row, col] = matrix1[row, col] + matrix2[row, col];
                    }
                }
                return sumMatrix;
            }
            
        }

        public static double[,] MatrixAdditionParallelNew(double[,] matrix1, double[,] matrix2)
        {
            if (matrix1.GetLength(0) != matrix2.GetLength(0) || matrix1.GetLength(1) != matrix2.GetLength(1))
            {
                throw new Exception("Not equally sized matrices");
            }
            var sumMatrix = new double[matrix1.GetLength(0), matrix1.GetLength(1)];
            Parallel.For(0, matrix1.GetLength(0), i =>
            {
                Parallel.For(0, matrix1.GetLength(1), j =>
                {
                    sumMatrix[i, j] = matrix1[i, j] + matrix2[i, j];
                });
            });
            return sumMatrix;
        }

        private static void DoAdditionCalculations(int row, double[,] matrix1, double[,] matrix2)
        {
            if (matrix1.GetLength(0) != matrix2.GetLength(0) || matrix1.GetLength(1) != matrix2.GetLength(1))
            {
                throw new Exception("Not equally sized matrices");
            }
            int matrixrows = matrix1.GetLength(0);
            int matrixcols = matrix1.GetLength(1);
            double[,] sumMatrix = new double[matrixrows, matrixcols];

            for (int col = 0; col < matrixcols; col++)
            {
                sumMatrix[row, col] = matrix1[row, col] + matrix2[row, col];
            }

            TempVariable = sumMatrix;
        }

        public static void MatrixAdditionParallel(double[,] matrix1, double[,] matrix2)
        {
            Parallel.For(0, 2000, row => DoAdditionCalculations(row, matrix1, matrix2));
        }

        public static double[,] MatrixSubtraction(double[,] matrix1, double[,] matrix2)
        {
            int matrixrows = matrix1.GetLength(0);
            int matrixcols = matrix1.GetLength(1);
            double[,] sumMatrix = new double[matrixrows, matrixcols];

            for (int row = 0; row < matrixrows; row++)
            {
                for (int col = 0; col < matrixcols; col++)
                {
                    sumMatrix[row, col] = matrix1[row, col] - matrix2[row, col];
                }
            }
            return sumMatrix;
        }

        public static double[,] MatrixAdditionParallel2(double[,] matrix1, double[,] matrix2)
        {
            int totalRows = matrix1.GetLength(0);
            TempVariable = matrix1;
            Task first = Task.Run(() => MatrixAdditionParallel2Calculations(matrix2, 0, totalRows/8));
            Task second = Task.Run(() => MatrixAdditionParallel2Calculations(matrix2, totalRows / 8, 2*totalRows / 8));
            Task third= Task.Run(() => MatrixAdditionParallel2Calculations(matrix2, 2*totalRows / 8, 3* totalRows / 8));
            Task fourth = Task.Run(() => MatrixAdditionParallel2Calculations(matrix2, 3* totalRows / 8, totalRows/2));
            Task fifth = Task.Run(() => MatrixAdditionParallel2Calculations(matrix2, totalRows/2, 5*totalRows / 8));
            Task sixth = Task.Run(() => MatrixAdditionParallel2Calculations(matrix2, 5*totalRows / 8, 6*totalRows / 8));
            Task seventh = Task.Run(() => MatrixAdditionParallel2Calculations(matrix2, 6*totalRows / 8, 7 * totalRows / 8));
            Task eightieth = Task.Run(() => MatrixAdditionParallel2Calculations(matrix2, 7 * totalRows / 8, totalRows));
            Task.WaitAll(first, second, third, fourth, fifth, sixth, seventh, eightieth);
            return TempVariable;
        }
        private static void MatrixAdditionParallel2Calculations(double[,] matrix2, int min, int max)
        {
            //if (matrix1.GetLength(0) != matrix2.GetLength(0) || matrix1.GetLength(1) != matrix2.GetLength(1))
            //{
            //    throw new Exception("Not equally sized matrices");
            //}
            //int matrixrows = matrix1.GetLength(0);
            int matrixcols = matrix2.GetLength(1);
            //double[,] sumMatrix = new double[matrixrows, matrixcols];

            for (int row = min; row < max; row++)
            {
                for (int col = 0; col < matrixcols; col++)
                {
                    TempVariable[row, col] = TempVariable[row, col] + matrix2[row, col];
                }
            }

        }

        public static double[,] ScalarMatrixProduct(double scalar, double[,] matrix)
        {
            int matrixrows = matrix.GetLength(0);
            int matrixcols = matrix.GetLength(1);
            for (int row = 0; row < matrixrows; row++)
            {
                for (int col = 0; col < matrixcols; col++)
                {
                    matrix[row, col] = scalar * matrix[row, col];
                }
            }
            return matrix;
        }

        public static double[,] ScalarMatrixProductNew(double scalar, double[,] matrix)
        {
            int matrixrows = matrix.GetLength(0);
            int matrixcols = matrix.GetLength(1);
            double[,] resultMatrix = new double[matrixrows, matrixcols];
            for (int row = 0; row < matrixrows; row++)
            {
                for (int col = 0; col < matrixcols; col++)
                {
                    resultMatrix[row, col] = scalar * matrix[row, col];
                }
            }
            return resultMatrix;
        }

        public static double[,] ScalarMatrixProductInParallel(double scalar, double[,] matrix)
        {
            double[,] resultMatrix = new double[matrix.GetLength(0), matrix.GetLength(1)];
            Parallel.For(0, matrix.GetLength(0), row =>
            {
                Parallel.For(0, matrix.GetLength(1), col =>
                {
                    resultMatrix[row, col] = scalar * matrix[row, col];
                });
            });
            return resultMatrix;
        }

        public static double[] ScalarByVectorProduct(double scalarFactor, double[] initialVector)
        {
            int dimension = initialVector.GetLength(0);
            double[] finalVector = new double[dimension];

            for (int row = 0; row < dimension; row++)
            {
                finalVector[row] = scalarFactor * initialVector[row];
            }
            return finalVector;
        }

        public static double[,] PutZerosInDiag(double[,] matrix)
        {
            for (int i = 0; i < matrix.GetLength(0); i++)
            {
                matrix[i, i] = 0;
            }
            return matrix;
        }

        public static double[,] GetDiagMatrix(double[,] matrix)
        {
            double[,] diagMatrix = new double[matrix.GetLength(0), matrix.GetLength(1)];
            for (int i = 0; i < matrix.GetLength(0); i++)
            {
                diagMatrix[i, i] = matrix[i, i];
            }
            return diagMatrix;
        }

        public static double[,] InvertDiagMatrix(double[,] matrix)
        {
            for (int i = 0; i < matrix.GetLength(0); i++)
            {
                matrix[i, i] = 1 / matrix[i, i];
            }
            return matrix;
        }

        public static void PrintMatrixToFile(double[,] matrix, string path)
        {
            int rows = matrix.GetLength(0);
            int columns = matrix.GetLength(1);
            string[] dataToPrint = new string[rows];
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < columns; j++)
                {
                    dataToPrint[i] = dataToPrint[i] + "\t" + matrix[i, j];
                }
            }
            File.WriteAllLines(path, dataToPrint);
        }
        public static void PrintMatrixToFile2(double[,] matrix, string path)
        {
            int rows = matrix.GetLength(0);
            int columns = matrix.GetLength(1);
            int elements = rows * columns;
            string[] dataToPrint = new string[elements];
            int count = 0;
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < columns; j++)
                {
                    dataToPrint[count] = matrix[i, j].ToString();
                    count += 1;
                }
            }
            File.WriteAllLines(path, dataToPrint);
        }
        public static double Trace(double[,] matrix)
        {
            int matrixrows = matrix.GetLength(0);
            double trace = new double();

            for (int row = 0; row < matrixrows; row++)
            {
                trace += matrix[row, row];
            }
            return trace;
        }
        public static double[,] IChol(double[,] matrix, int[,] fillInLevels, int fillLevel)
        {
            int rows = matrix.GetLength(0);
            int cols = matrix.GetLength(1);

            double[,] lowerPart = new double[rows, cols];
            double sumr;
            double sumc;

            for (int j = 0; j < cols; j++)
            {
                sumr = 0;
                for (int k = 0; k <= j - 1; k++)
                {
                    sumr = sumr + Math.Pow(lowerPart[j, k], 2);
                }
                if (matrix[j, j] - sumr < 0)
                {
                    throw new Exception("Cholesky: Negative number in square root");
                }
                lowerPart[j, j] = Math.Sqrt(matrix[j, j] - sumr);
                for (int i = j + 1; i < rows; i++)
                {
                    //if (matrix[i, j] == 0.0)
                    //{
                    //    lowerPart[i, j] = 0.0;
                    //}
                    //else
                    //{
                    //    sumc = 0;
                    //    for (int k = 0; k <= j - 1; k++)
                    //    {
                    //        sumc = sumc + lowerPart[i, k] * lowerPart[j, k];
                    //    }
                    //    lowerPart[i, j] = (matrix[i, j] - sumc) / lowerPart[j, j];
                    //}
                    if(fillInLevels[i,j]<= fillLevel)
                    {
                        sumc = 0;
                        for (int k = 0; k <= j - 1; k++)
                        {
                            sumc = sumc + lowerPart[i, k] * lowerPart[j, k];
                        }
                        lowerPart[i, j] = (matrix[i, j] - sumc) / lowerPart[j, j];
                    }
                }
            }
            return lowerPart;
        }
        public static int[,] ICholLevels(double[,] matrix, int fillInDegree)
        {
            int rows = matrix.GetLength(0);
            int cols = matrix.GetLength(1);

            int[,] level = new int[rows, cols];
            level = FillMatrixWithConstantValue(level, 1000);

            for (int j = 0; j < cols; j++)
            {
                for (int i = 0; i < rows; i++)
                {
                    if (matrix[i, j] != 0 || i == j)
                    {
                        level[i, j] = 0;
                    }
                }
            }

            for (int i = 1; i < rows; i++)
            {
                for (int k = 0; k <= i - 1; k++)
                {
                    for (int j = k + 1; j < cols; j++)
                    {
                        level[i, j] = Math.Min(level[i, j], level[i, k] + level[k, j] + 1);
                    }
                }
            }
            return level;
        }
        public static double[,] CChol(double[,] matrix)
        {
            int rows = matrix.GetLength(0);
            int cols = matrix.GetLength(1);

            double[,] lowerPart = new double[rows, cols];
            double sumr;
            double sumc;

            for (int j = 0; j < cols; j++)
            {
                sumr = 0;
                for (int k = 0; k <= j - 1; k++)
                {
                    sumr = sumr + Math.Pow(lowerPart[j, k], 2);
                }
                if (matrix[j, j] - sumr < 0)
                {
                    throw new Exception("Cholesky: Negative number in square root");
                }
                lowerPart[j, j] = Math.Sqrt(matrix[j, j] - sumr);
                for (int i = j + 1; i < rows; i++)
                {
                    sumc = 0;
                    for (int k = 0; k <= j - 1; k++)
                    {
                        sumc = sumc + lowerPart[i, k] * lowerPart[j, k];
                    }
                    lowerPart[i, j] = (matrix[i, j] - sumc) / lowerPart[j, j];
                }
            }
            return lowerPart;
        }
        public static double[,] CreateBlockMatrix(double[,] A11, double[,] A12, double[,] A21, double[,] A22)
        {
            double[,] A = new double[A11.GetLength(0) + A22.GetLength(0), A11.GetLength(0) + A22.GetLength(0)];
            for(int i = 0; i< A.GetLength(0); i++)
            {
                for(int j = 0; j < A.GetLength(0); j++)
                {
                    if(i< A11.GetLength(0) && j < A11.GetLength(0))
                    {
                        A[i, j] = A11[i, j];
                    }
                    else if (i < A11.GetLength(0) && j >= A11.GetLength(0))
                    {
                        A[i, j] = A12[i, j - A11.GetLength(0)];
                    }
                    else if (i >= A11.GetLength(0) && j < A11.GetLength(0))
                    {
                        A[i, j] = A21[i - A11.GetLength(0), j];
                    }
                    else if (i >= A11.GetLength(0) && j >= A11.GetLength(0))
                    {
                        A[i, j] = A22[i - A11.GetLength(0), j - A11.GetLength(0)];
                    }
                }
            }
            return A;
        }
        public static int[,] FillMatrixWithConstantValue(int[,] matrix, int value)
        {
            int[,] result = new int[matrix.GetLength(0), matrix.GetLength(1)];
            for (int i = 0; i < matrix.GetLength(0); i++)
            {
                for (int j = 0; j < matrix.GetLength(0); j++)
                {
                    result[i, j] = value;
                }
            }
            return result;
        }

        public static void FillMatrixWithDoubleValue(double[,] matrix, double value)
        {
            for (int i = 0; i < matrix.GetLength(0); i++)
            {
                for (int j = 0; j < matrix.GetLength(0); j++)
                {
                    matrix[i, j] = value;
                }
            }
           
        }

        public static Tuple<double[,], double[,]> LUfactorizedMatrices(double[,] Matrix)
        {
            int rows = Matrix.GetLength(0);
            int cols = Matrix.GetLength(1);

            double[,] lowerPart = new double[rows, cols];
            double[,] upperPart = new double[rows, cols];
            double sumu;
            double suml;

            for (int j = 0; j < cols; j++)
            {
                upperPart[0, j] = Matrix[0, j];
            }

            for (int k = 0; k < rows; k++)
            {
                lowerPart[k, k] = 1.0;
            }

            for (int j = 0; j < cols; j++)
            {
                for (int i = 0; i < rows; i++)
                {
                    if (i <= j)
                    {
                        sumu = 0.0;
                        for (int k = 0; k < i; k++)
                        {
                            sumu = sumu + lowerPart[i, k] * upperPart[k, j];
                        }
                        upperPart[i, j] = (Matrix[i, j] - sumu) / lowerPart[i, i];
                    }
                    else if (i > j)
                    {
                        suml = 0.0;
                        for (int k = 0; k < j; k++)
                        {
                            suml = suml + lowerPart[i, k] * upperPart[k, j];
                        }
                        lowerPart[i, j] = (Matrix[i, j] - suml) / upperPart[j, j];
                    }
                }
            }
            return new Tuple<double[,], double[,]>(lowerPart, upperPart);
        }
        public static bool CheckIfSymmetric(double[,] matrix)
        {
            int n = matrix.GetLength(0);
            int m = matrix.GetLength(1);
            double tol = 0.001;
            bool value = true;
            if (!(n == m))
            {
                throw new Exception("The matrix is not square");
            }
            else
            {
                for(int i = 0; i < n; i++)
                {
                    for(int j = i + 1; j < n; j++)
                    {
                        if(Math.Abs(matrix[i, j] - matrix[j, i])>tol)
                        {
                            value = false;
                            break;
                        }
                    }
                    if (!value) break; else continue;
                }
            }
            return value;
        }
        public static double[,] CreateDiagonalMatrix(int dimension, double[,] diagonalMatrix)
        {
            if(diagonalMatrix.GetLength(0)!=diagonalMatrix.GetLength(1)||
                dimension % diagonalMatrix.GetLength(0) != 0)
            {
                throw new Exception("Dimensions mismach");
            }
            
            double[,] diagMatrix = new double[dimension, dimension];
            int numberOfNodes = dimension / diagonalMatrix.GetLength(0);
            for (int k =0; k<numberOfNodes; k++)
            {
                for (int i = 0; i < diagonalMatrix.GetLength(0); i++)
                {
                    for (int j = 0; j < diagonalMatrix.GetLength(0); j++)
                    {
                        diagMatrix[k* diagonalMatrix.GetLength(0) + i, k * diagonalMatrix.GetLength(0) + j] = diagonalMatrix[i, j];
                    }
                }
            }
            return diagMatrix;
        }
        //public static double[,] CorrectSymmetry(double[,] matrix)
        //{
        //    int n = matrix.GetLength(0);
        //    int m = matrix.GetLength(1);
        //    double tol = 0.000001;
        //    bool value = true;
        //    if (!(n == m))
        //    {
        //        throw new Exception("The matrix is not square");
        //    }
        //    else
        //    {
        //        for (int i = 0; i < n; i++)
        //        {
        //            for (int j = i + 1; j < n; j++)
        //            {
        //                if (Math.Abs(matrix[i, j] - matrix[j, i]) > tol)
        //                {
        //                    value = false;
        //                    break;
        //                }
        //                else if(!(matrix[i, j] == matrix[j, i]))
        //                {
        //                    matrix[j, i] = matrix[i, j];
        //                }
        //            }
        //            if (!value) break; else continue;
        //        }
        //    }
        //    return matrix;
        //}
        public static double[,] CalculateInverse3X3Matrix(double[,] matrix)
        {
            double[,] inverseMatrix = new double[3, 3];

            inverseMatrix[0, 0] = matrix[1, 1] * matrix[2, 2] - matrix[1, 2] * matrix[2, 1];
            inverseMatrix[1, 1] = matrix[2, 2] * matrix[0, 0] - matrix[2, 0] * matrix[2, 2];
            inverseMatrix[2, 2] = matrix[0, 0] * matrix[1, 1] - matrix[0, 1] * matrix[1, 0];

            inverseMatrix[0, 1] = matrix[1, 2] * matrix[2, 0] - matrix[1, 0] * matrix[2, 2];
            inverseMatrix[1, 2] = matrix[2, 0] * matrix[0, 1] - matrix[2, 1] * matrix[0, 0];
            inverseMatrix[2, 0] = matrix[0, 1] * matrix[1, 2] - matrix[0, 2] * matrix[1, 1];

            inverseMatrix[1, 0] = matrix[2, 1] * matrix[0, 2] - matrix[0, 1] * matrix[2, 2];
            inverseMatrix[2, 1] = matrix[0, 2] * matrix[1, 0] - matrix[1, 2] * matrix[0, 0];
            inverseMatrix[0, 2] = matrix[1, 0] * matrix[1, 1] - matrix[2, 0] * matrix[1, 1];

            double det = matrix[0, 0] * inverseMatrix[0, 0] + matrix[0, 1] * inverseMatrix[1, 0] + matrix[0, 2] * inverseMatrix[2, 0];

            inverseMatrix[0, 0] = inverseMatrix[0, 0] / det;
            inverseMatrix[1, 1] = inverseMatrix[1, 1] / det;
            inverseMatrix[2, 2] = inverseMatrix[2, 2] / det;

            inverseMatrix[0, 1] = inverseMatrix[0, 1] / det;
            inverseMatrix[1, 2] = inverseMatrix[1, 2] / det;
            inverseMatrix[2, 0] = inverseMatrix[2, 0] / det;

            inverseMatrix[1, 0] = inverseMatrix[1, 0] / det;
            inverseMatrix[2, 1] = inverseMatrix[2, 1] / det;
            inverseMatrix[0, 2] = inverseMatrix[0, 2] / det;

            return inverseMatrix;
        }
        public static double[,] BlockMatrixInversion6X6(double[,] matrix)
        {
            double[,] inverseMatrix = new double[6, 6];
            double[,] A = new double[3, 3];
            double[,] B = new double[3, 3];
            double[,] C = new double[3, 3];
            double[,] D = new double[3, 3];

            for (int i = 0; i < A.GetLength(0); i++)
            {
                for (int j = 0; j < A.GetLength(0); j++)
                {
                    A[i, j] = matrix[i, j];
                    B[i, j] = matrix[i, j + 3];
                    C[i, j] = matrix[i + 3, j];
                    D[i, j] = matrix[i + 3, j + 3];

                }
            }
            double[,] AInv = CalculateInverse3X3Matrix(A);
            double[,] ASchur = MatrixSubtraction(D,
                MatrixProduct(MatrixProduct(C, AInv), B));
            double[,] ASchurInv = CalculateInverse3X3Matrix(ASchur);
            double[,] inv11 = MatrixAddition(AInv, MatrixProduct(
                MatrixProduct(MatrixProduct(MatrixProduct(AInv,B),
                ASchurInv),C), AInv));
            double[,] inv12 = ScalarMatrixProductNew(-1.0,
                MatrixProduct(MatrixProduct(AInv, B),
                ASchurInv));
            double[,] inv21 = ScalarMatrixProductNew(-1.0,
                MatrixProduct(MatrixProduct(ASchurInv, C),
                AInv));
            double[,] inv22 = ASchurInv;
            for (int i = 0; i < inverseMatrix.GetLength(0); i++)
            {
                for (int j = 0; j < inverseMatrix.GetLength(0); j++)
                {
                    if(i < 3 && j < 3)
                    {
                        inverseMatrix[i, j] = inv11[i, j];
                    }
                    else if(i < 3 && j >= 3)
                    {
                        inverseMatrix[i, j] = inv12[i, j - 3];
                    }
                    else if (i >= 3 && j < 3)
                    {
                        inverseMatrix[i, j] = inv21[i - 3, j];
                    }
                    else if (i >= 3 && j >= 3)
                    {
                        inverseMatrix[i, j] = inv22[i - 3, j - 3];
                    }
                }
            }
            return inverseMatrix;
        }
        public static double[,] CalculateInverse2X2Matrix(double[,] matrix)
        {
            double[,] inverseMatrix = new double[2, 2];

            double detj = matrix[0, 0] * matrix[1, 1] - matrix[0, 1] * matrix[1, 0];

            inverseMatrix[0, 0] = matrix[1, 1] / detj;
            inverseMatrix[0, 1] = -matrix[0, 1] / detj;
            inverseMatrix[1, 0] = -matrix[1, 0] / detj;
            inverseMatrix[1, 1] = matrix[0, 0] / detj;

            return inverseMatrix;
        }
        public static double[,] BlockMatrixInversion4X4(double[,] matrix)
        {
            double[,] inverseMatrix = new double[4, 4];
            double[,] A = new double[2, 2];
            double[,] B = new double[2, 2];
            double[,] C = new double[2, 2];
            double[,] D = new double[2, 2];

            for (int i = 0; i < A.GetLength(0); i++)
            {
                for (int j = 0; j < A.GetLength(0); j++)
                {
                    A[i, j] = matrix[i, j];
                    B[i, j] = matrix[i, j + 2];
                    C[i, j] = matrix[i + 2, j];
                    D[i, j] = matrix[i + 2, j + 2];
                }
            }
            double[,] AInv = CalculateInverse2X2Matrix(A);
            double[,] ASchur = MatrixSubtraction(D,
                MatrixProduct(MatrixProduct(C, AInv), B));
            double[,] ASchurInv = CalculateInverse2X2Matrix(ASchur);
            double[,] inv11 = MatrixAddition(AInv, MatrixProduct(
                MatrixProduct(MatrixProduct(MatrixProduct(AInv, B),
                ASchurInv), C), AInv));
            double[,] inv12 = ScalarMatrixProductNew(-1.0,
                MatrixProduct(MatrixProduct(AInv, B),
                ASchurInv));
            double[,] inv21 = ScalarMatrixProductNew(-1.0,
                MatrixProduct(MatrixProduct(ASchurInv, C),
                AInv));
            double[,] inv22 = ASchurInv;
            for (int i = 0; i < inverseMatrix.GetLength(0); i++)
            {
                for (int j = 0; j < inverseMatrix.GetLength(0); j++)
                {
                    if (i < 2 && j < 2)
                    {
                        inverseMatrix[i, j] = inv11[i, j];
                    }
                    else if (i < 2 && j >= 2)
                    {
                        inverseMatrix[i, j] = inv12[i, j - 2];
                    }
                    else if (i >= 2 && j < 2)
                    {
                        inverseMatrix[i, j] = inv21[i - 2, j];
                    }
                    else if (i >= 2 && j >= 2)
                    {
                        inverseMatrix[i, j] = inv22[i - 2, j - 2];
                    }
                }
            }
            return inverseMatrix;
        }
        public static double[,] ΕliminateMatrixLine(double[,] matrix, int line)
        {
            double[,] resultMatrix = matrix;
            for (int i = 0; i < matrix.GetLength(1); i++)
            {
                resultMatrix[line, i] = 0.0;
            }
            return resultMatrix;
        }
        public static double[,] CalculateInverseMatrixGaussJordanMethod(double[,] matrix)
        {
            if (!(matrix.GetLength(0) == matrix.GetLength(1)))
            {
                throw new Exception("The matrix must be square");
            }
            double[,] augmentedMatrix = new double[matrix.GetLength(0), 2 * matrix.GetLength(0)];
            double[,] inverseMatrix = new double[matrix.GetLength(0), matrix.GetLength(0)];

            for (int i = 0; i < matrix.GetLength(0); i++)
            {
                for (int j = 0; j < matrix.GetLength(0); j++)
                {
                    augmentedMatrix[i, j] = matrix[i, j];
                    if (i == j)
                    {
                        augmentedMatrix[i, j + matrix.GetLength(0)] = 1.0;
                    }
                    else
                    {
                        augmentedMatrix[i, j + matrix.GetLength(0)] = 0.0;
                    }
                }
            }
            for (int i = 0; i < matrix.GetLength(0); i++)
            {
                if (augmentedMatrix[i, i] == 0.0)
                {
                    throw new Exception("Zero diagonal element!");
                }
                for (int j = 0; j < matrix.GetLength(0); j++)
                {
                    if (i != j)
                    {
                        double ratio = augmentedMatrix[j, i] / augmentedMatrix[i, i];
                        for (int k = 0; k < 2 * matrix.GetLength(0); k++)
                        {
                            augmentedMatrix[j, k] = augmentedMatrix[j, k] - ratio * augmentedMatrix[i, k];
                        }
                    }
                }
            }
            for (int i = 0; i < matrix.GetLength(0); i++)
            {
                for (int j = matrix.GetLength(0); j < 2 * matrix.GetLength(0); j++)
                {
                    augmentedMatrix[i, j] = augmentedMatrix[i, j] / augmentedMatrix[i, i];
                }
            }
            for (int i = 0; i < matrix.GetLength(0); i++)
            {
                for (int j = 0; j < matrix.GetLength(0); j++)
                {
                    inverseMatrix[i, j] = augmentedMatrix[i, j + matrix.GetLength(0)];
                }
            }
            return inverseMatrix;
        }
        public static bool CheckDiagonalElements(double[,] matrix)
        {
            bool allPositive = true;
            for (int i = 0; i< matrix.GetLength(0); i++)
            {
                if (matrix[i, i] <= 0)
                {
                    allPositive = false;
                    break;
                }
            }
            return allPositive;
        }
        public static double GetMatrixMaxValue(double[,] matrix)
        {
            var max = new double();
            for(var i = 0; i< matrix.GetLength(0); i++)
            {
                var rowMax = new double();
                for(var j =0; j<matrix.GetLength(1); j++)
                {
                    if (matrix[i, j] > rowMax)
                    {
                        rowMax = matrix[i, j];
                    }
                }
                if (rowMax > max)
                {
                    max = rowMax;
                }
            }
            return max;
        }

        public static double GetMatrixMinValue(double[,] matrix)
        {
            var min = new double();
            for (var i = 0; i < matrix.GetLength(0); i++)
            {
                var rowMin = new double();
                for (var j = 0; j < matrix.GetLength(1); j++)
                {
                    if (matrix[i, j] < rowMin)
                    {
                        rowMin = matrix[i, j];
                    }
                }
                if (rowMin < min)
                {
                    min = rowMin;
                }
            }
            return min;
        }
    }
}
