﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace BumperCollitionSimulation
{
    class Hex8 : IElement
    {
        public Dictionary<int, INode> Nodes { get; }
        public IElementProperties Properties { get; set; }
        public Dictionary<int, bool[]> ElementFreedomSignature { get; } = new Dictionary<int, bool[]>();
        public List<int> ElementFreedomList { get; set; }
        public double[] DisplacementVector { get; set; }
        public double[] AccelerationVector { get; set; }
        public double poisson { get; set; }
        public void InitializeTangentialProperties()
        {
            throw new Exception("Needs to be removed. Has beeb used only for testing purposes");
        }
        public void UpdateTangentialProperties()
        {
            throw new Exception("Needs to be removed. Has beeb used only for testing purposes");
        }
        public void InitializeContactSurfaceGeometry()
        {
            throw new Exception("Needs to be removed. Has beeb used only for testing purposes");
        }
        public void UpdateContactSurfaceGeometry()
        {
            throw new Exception("Needs to be removed. Has beeb used only for testing purposes");
        }
        public void UpdateIncrementalDisplacements(double[] deltaU)
        {
            throw new Exception("Needs to be removed. Has beeb used only for testing purposes");
        }
        public Hex8(IElementProperties properties, Dictionary<int, INode> nodes)
        {
            Properties = properties;
            this.Nodes = nodes;
            ElementFreedomSignature[1] = new bool[] { true, true, true, false, false, false };
            ElementFreedomSignature[2] = new bool[] { true, true, true, false, false, false };
            ElementFreedomSignature[3] = new bool[] { true, true, true, false, false, false };
            ElementFreedomSignature[4] = new bool[] { true, true, true, false, false, false };
            ElementFreedomSignature[5] = new bool[] { true, true, true, false, false, false };
            ElementFreedomSignature[6] = new bool[] { true, true, true, false, false, false };
            ElementFreedomSignature[7] = new bool[] { true, true, true, false, false, false };
            ElementFreedomSignature[8] = new bool[] { true, true, true, false, false, false };
            DisplacementVector = new double[24];
        }
        public void CalculateElementEASMatrices()
        {
            throw new Exception("This method is to be used only for EAS method elements");
        }
        public void InitializeElementEASParameters()
        {
            throw new Exception("Needs to be removed. Has beeb used only for testing purposes");
        }
        public void UpdateElementEASParameters(double[] solutionVector)
        {
            throw new Exception("Needs to be removed. Has beeb used only for testing purposes");
        }
        public void StoreElementFinalStepDisplacementVector(double[] solutionVector)
        {
            throw new Exception("Needs to be removed. Has beeb used only for testing purposes");
        }
        public double ClosestPointProjection()
        {
            throw new Exception("Needs to be removed. Has beeb used only for testing purposes");
        }

        public Dictionary<int, INode> NodesAtFinalState()
        {
            throw new Exception("Method not implemenented");
        }
        public List<double[]> GetStressVector()
        {
            List<double[]> GpointsStress = new List<double[]>();
            double[,] E = CalculateStressStrainMatrix(Properties.YoungMod, Properties.PoissonRatio);
            for (int i = 0; i < 2; i++)
            {
                for (int j = 0; j < 2; j++)
                {
                    for (int k = 0; k < 2; k++)
                    {
                        double[] gP = GaussPoints(i, j, k).Item1;
                        double[] gW = GaussPoints(i, j, k).Item2;
                        Dictionary<string, double[]> localdN = CalculateShapeFunctionsLocalDerivatives(gP);
                        double[,] J = CalculateJacobian(localdN);
                        double[,] invJ = CalculateInverseJacobian(J).Item1;
                        double detJ = CalculateInverseJacobian(J).Item2;
                        Dictionary<int, double[]> globaldN = CalculateShapeFunctionsGlobalDerivatives(localdN, invJ);
                        double[,] B = CalculateBMatrix(globaldN);
                        double[] strainVector = CalculateStrainsVector(B);
                        double[] stressVector = CalculateStressVector(E, strainVector);
                        GpointsStress.Add(stressVector);
                    }
                }
            }
            return GpointsStress;
        }
        public List<double[]> GetStrainVector()
        {
            List<double[]> GpointsDeformation = new List<double[]>();
            for (int i = 0; i < 2; i++)
            {
                for (int j = 0; j < 2; j++)
                {
                    for (int k = 0; k < 2; k++)
                    {
                        double[] gP = GaussPoints(i, j, k).Item1;
                        double[] gW = GaussPoints(i, j, k).Item2;
                        Dictionary<string, double[]> localdN = CalculateShapeFunctionsLocalDerivatives(gP);
                        double[,] J = CalculateJacobian(localdN);
                        double[,] invJ = CalculateInverseJacobian(J).Item1;
                        double detJ = CalculateInverseJacobian(J).Item2;
                        Dictionary<int, double[]> globaldN = CalculateShapeFunctionsGlobalDerivatives(localdN, invJ);
                        double[,] B = CalculateBMatrix(globaldN);
                        double[] strainVector = CalculateStrainsVector(B);
                        GpointsDeformation.Add(strainVector);
                    }
                }
            }
            return GpointsDeformation;
        }
        public List<double[]> GetGaussPointsInPhysicalSpace()
        {
            List<double[]> GpointsPhysicalCoordinates = new List<double[]>();
            double[] xUpdated = UpdateNodalCoordinates(DisplacementVector);
            for (int i = 0; i < 2; i++)
            {
                for (int j = 0; j < 2; j++)
                {
                    for (int k = 0; k < 2; k++)
                    {
                        double[] gP = GaussPoints(i, j, k).Item1;
                        double[] gaussPoint = VectorOperations.MatrixVectorProduct(CalculateShapeFunctionMatrix(gP[0], gP[1], gP[2]), xUpdated);
                        GpointsPhysicalCoordinates.Add(gaussPoint);
                    }
                }
            }
            return GpointsPhysicalCoordinates;
        }
        public List<double[]> GetStressFromElementsNodes()
        {
            List<double[]> l = new List<double[]>();
            l.Add(new double[] { 0.0, 0.0, 0.0 });
            return l;
        }
        public List<double[]> GetStrainFromElementsNodes()
        {
            List<double[]> l = new List<double[]>();
            l.Add(new double[] { 0.0, 0.0, 0.0 });
            return l;
        }
        public List<double[]> GetStressFromElements(List<double[]> parametricCoordinatesVector)
        {
            List<double[]> StessVectorsList = new List<double[]>();
            double[,] E = CalculateStressStrainMatrix(Properties.YoungMod, Properties.PoissonRatio);
            //double[] xUpdated = UpdateNodalCoordinates(DisplacementVector);
            for (int i = 0; i < parametricCoordinatesVector.Count; i++)
            {
                Dictionary<string, double[]> localdN = CalculateShapeFunctionsLocalDerivatives(parametricCoordinatesVector[i]);
                double[,] J = CalculateJacobian(localdN);
                double[,] invJ = CalculateInverseJacobian(J).Item1;
                Dictionary<int, double[]> globaldN = CalculateShapeFunctionsGlobalDerivatives(localdN, invJ);
                double[,] B = CalculateBMatrix(globaldN);
                double[] strainVector = CalculateStrainsVector(B);
                double[] stressVector = CalculateStressVector(E, strainVector);
                StessVectorsList.Add(stressVector);
            }
            return StessVectorsList;
        }
        public List<double[]> GetphysicalCoordinatesFromElements(List<double[]> parametricCoordinatesVector)
        {
            List<double[]> PositionVectorsList = new List<double[]>();
            double[] xUpdated = UpdateNodalCoordinates(DisplacementVector);
            //PositionVectorsList.Add(new double[] { 0.0, 0.0 });
            int count = parametricCoordinatesVector.Count;
            for (int i = 0; i < count; i++)
            {
                double[] parametricCoordinatesVec = parametricCoordinatesVector[i];
                double[] positionVector = VectorOperations.MatrixVectorProduct(CalculateShapeFunctionMatrix(parametricCoordinatesVec[0], parametricCoordinatesVec[1], parametricCoordinatesVec[2]), xUpdated);
                PositionVectorsList.Add(positionVector);
            }
            return PositionVectorsList;
        }
        private double[] UpdateNodalCoordinates(double[] displacementVector)
        {
            double[] updatedCoor = new double[24];
            for (int i = 1; i <= 8; i++)
            {
                updatedCoor[3 * i - 3] = Nodes[i].XCoordinate + displacementVector[3 * i - 3];
                updatedCoor[3 * i - 2] = Nodes[i].YCoordinate + displacementVector[3 * i - 2];
                updatedCoor[3 * i - 1] = Nodes[i].ZCoordinate + displacementVector[3 * i - 1];
            }
            return updatedCoor;
        }
        private double[] InitialNodalCoordinates()
        {
            double[] initialCoor = new double[24];
            for (int i = 1; i <= 8; i++)
            {
                initialCoor[3 * i - 3] = Nodes[i].XCoordinate;
                initialCoor[3 * i - 2] = Nodes[i].YCoordinate;
                initialCoor[3 * i - 1] = Nodes[i].ZCoordinate;
            }
            return initialCoor;
        }
        private Dictionary<int, double> CalculateShapeFunctions(double ksi, double ihta, double mhi)
        {
            Dictionary<int, double> shapeFunctions = new Dictionary<int, double>();
            double N1 = 1.0 / 8.0 * (1 - ksi) * (1 - ihta) * (1 - mhi); shapeFunctions.Add(1, N1);
            double N2 = 1.0 / 8.0 * (1 + ksi) * (1 - ihta) * (1 - mhi); shapeFunctions.Add(2, N2);
            double N3 = 1.0 / 8.0 * (1 + ksi) * (1 + ihta) * (1 - mhi); shapeFunctions.Add(3, N3);
            double N4 = 1.0 / 8.0 * (1 - ksi) * (1 + ihta) * (1 - mhi); shapeFunctions.Add(4, N4);
            double N5 = 1.0 / 8.0 * (1 - ksi) * (1 - ihta) * (1 + mhi); shapeFunctions.Add(5, N5);
            double N6 = 1.0 / 8.0 * (1 + ksi) * (1 - ihta) * (1 + mhi); shapeFunctions.Add(6, N6);
            double N7 = 1.0 / 8.0 * (1 + ksi) * (1 + ihta) * (1 + mhi); shapeFunctions.Add(7, N7);
            double N8 = 1.0 / 8.0 * (1 - ksi) * (1 + ihta) * (1 + mhi); shapeFunctions.Add(8, N8);

            return shapeFunctions;
        }
        private double[,] CalculateShapeFunctionMatrix(double ksi, double ihta, double zita)
        {
            double N1 = 1.0 / 8.0 * (1 - ksi) * (1 - ihta) * (1 - zita);
            double N2 = 1.0 / 8.0 * (1 + ksi) * (1 - ihta) * (1 - zita);
            double N3 = 1.0 / 8.0 * (1 + ksi) * (1 + ihta) * (1 - zita);
            double N4 = 1.0 / 8.0 * (1 - ksi) * (1 + ihta) * (1 - zita);
            double N5 = 1.0 / 8.0 * (1 - ksi) * (1 - ihta) * (1 + zita);
            double N6 = 1.0 / 8.0 * (1 + ksi) * (1 - ihta) * (1 + zita);
            double N7 = 1.0 / 8.0 * (1 + ksi) * (1 + ihta) * (1 + zita);
            double N8 = 1.0 / 8.0 * (1 - ksi) * (1 + ihta) * (1 + zita);
            double[,] shapeFunctionsMat = new double[,] {
                {N1, 0.0, 0.0, N2, 0.0, 0.0, N3, 0.0, 0.0, N4, 0.0, 0.0, N5, 0.0, 0.0, N6, 0.0, 0.0, N7, 0.0, 0.0, N8, 0.0, 0.0 },
                {0.0, N1, 0.0, 0.0, N2, 0.0, 0.0, N3, 0.0, 0.0, N4, 0.0, 0.0, N5, 0.0, 0.0, N6, 0.0, 0.0, N7, 0.0, 0.0, N8, 0.0 },
                {0.0, 0.0, N1, 0.0, 0.0, N2, 0.0, 0.0, N3, 0.0, 0.0, N4, 0.0, 0.0, N5, 0.0, 0.0, N6, 0.0, 0.0, N7, 0.0, 0.0, N8 },
            };
            return shapeFunctionsMat;
        }

        private Dictionary<string, double[]> CalculateShapeFunctionsLocalDerivatives(double[] naturalCoordinates)
        {
            double ksi = naturalCoordinates[0];
            double ihta = naturalCoordinates[1];
            double mhi = naturalCoordinates[2];

            double[] dN_ksi = new double[]
            {
                (-1.0/8.0*(1-ihta)*(1-mhi)),
                (1.0/8.0*(1-ihta)*(1-mhi)),
                (1.0/8.0*(1+ihta)*(1-mhi)),
                (-1.0/8.0*(1+ihta)*(1-mhi)),
                (-1.0/8.0*(1-ihta)*(1+mhi)),
                (1.0/8.0*(1-ihta)*(1+mhi)),
                (1.0/8.0*(1+ihta)*(1+mhi)),
                (-1.0/8.0*(1+ihta)*(1+mhi))
            };

            double[] dN_ihta = new double[]
            {
                (-1.0/8.0*(1-ksi)*(1-mhi)),
                (-1.0/8.0*(1+ksi)*(1-mhi)),
                (1.0/8.0*(1+ksi)*(1-mhi)),
                (1.0/8.0*(1-ksi)*(1-mhi)),
                (-1.0/8.0*(1-ksi)*(1+mhi)),
                (-1.0/8.0*(1+ksi)*(1+mhi)),
                (1.0/8.0*(1+ksi)*(1+mhi)),
                (1.0/8.0*(1-ksi)*(1+mhi))
            };

            double[] dN_mhi = new double[]
            {
                (-1.0/8.0*(1-ksi)*(1-ihta)),
                (-1.0/8.0*(1+ksi)*(1-ihta)),
                (-1.0/8.0*(1+ksi)*(1+ihta)),
                (-1.0/8.0*(1-ksi)*(1+ihta)),
                (1.0/8.0*(1-ksi)*(1-ihta)),
                (1.0/8.0*(1+ksi)*(1-ihta)),
                (1.0/8.0*(1+ksi)*(1+ihta)),
                (1.0/8.0*(1-ksi)*(1+ihta))
            };

            Dictionary<string, double[]> dN = new Dictionary<string, double[]>();
            dN.Add("ksi", dN_ksi);
            dN.Add("ihta", dN_ihta);
            dN.Add("mhi", dN_mhi);
            return dN;
        }

        private double[,] CalculateJacobian(Dictionary<string, double[]> dN)
        {
            double[,] jacobianMatrix = new double[3, 3];
            //DisplacementVector = new double[24];
            //double[] xUpdated = UpdateNodalCoordinates(DisplacementVector);
            double[] xNodalInitial = InitialNodalCoordinates();

            int k = 0;
            for (int i = 0; i < 8; i++)
            {
                jacobianMatrix[0, 0] = jacobianMatrix[0, 0] + xNodalInitial[k] * dN["ksi"][i];
                k = k + 3;
            }
            k = 1;
            for (int i = 0; i < 8; i++)
            {
                jacobianMatrix[0, 1] = jacobianMatrix[0, 1] + xNodalInitial[k] * dN["ksi"][i];
                k = k + 3;
            }
            k = 2;
            for (int i = 0; i < 8; i++)
            {
                jacobianMatrix[0, 2] = jacobianMatrix[0, 2] + xNodalInitial[k] * dN["ksi"][i];
                k = k + 3;
            }

            k = 0;
            for (int i = 0; i < 8; i++)
            {
                jacobianMatrix[1, 0] = jacobianMatrix[1, 0] + xNodalInitial[k] * dN["ihta"][i];
                k = k + 3;
            }
            k = 1;
            for (int i = 0; i < 8; i++)
            {
                jacobianMatrix[1, 1] = jacobianMatrix[1, 1] + xNodalInitial[k] * dN["ihta"][i];
                k = k + 3;
            }
            k = 2;
            for (int i = 0; i < 8; i++)
            {
                jacobianMatrix[1, 2] = jacobianMatrix[1, 2] + xNodalInitial[k] * dN["ihta"][i];
                k = k + 3;
            }

            k = 0;
            for (int i = 0; i < 8; i++)
            {
                jacobianMatrix[2, 0] = jacobianMatrix[2, 0] + xNodalInitial[k] * dN["mhi"][i];
                k = k + 3;
            }
            k = 1;
            for (int i = 0; i < 8; i++)
            {
                jacobianMatrix[2, 1] = jacobianMatrix[2, 1] + xNodalInitial[k] * dN["mhi"][i];
                k = k + 3;
            }
            k = 2;
            for (int i = 0; i < 8; i++)
            {
                jacobianMatrix[2, 2] = jacobianMatrix[2, 2] + xNodalInitial[k] * dN["mhi"][i];
                k = k + 3;
            }

            return jacobianMatrix;
        }

        private Tuple<double[,], double> CalculateInverseJacobian(double[,] jacobianMatrix)
        {
            double[,] jacobianInverseMatrix = new double[3, 3];

            jacobianInverseMatrix[0, 0] = jacobianMatrix[1, 1] * jacobianMatrix[2, 2] - jacobianMatrix[1, 2] * jacobianMatrix[2, 1];
            jacobianInverseMatrix[1, 1] = jacobianMatrix[2, 2] * jacobianMatrix[0, 0] - jacobianMatrix[2, 0] * jacobianMatrix[2, 2];
            jacobianInverseMatrix[2, 2] = jacobianMatrix[0, 0] * jacobianMatrix[1, 1] - jacobianMatrix[0, 1] * jacobianMatrix[1, 0];

            jacobianInverseMatrix[0, 1] = jacobianMatrix[1, 2] * jacobianMatrix[2, 0] - jacobianMatrix[1, 0] * jacobianMatrix[2, 2];
            jacobianInverseMatrix[1, 2] = jacobianMatrix[2, 0] * jacobianMatrix[0, 1] - jacobianMatrix[2, 1] * jacobianMatrix[0, 0];
            jacobianInverseMatrix[2, 0] = jacobianMatrix[0, 1] * jacobianMatrix[1, 2] - jacobianMatrix[0, 2] * jacobianMatrix[1, 1];

            jacobianInverseMatrix[1, 0] = jacobianMatrix[2, 1] * jacobianMatrix[0, 2] - jacobianMatrix[0, 1] * jacobianMatrix[2, 2];
            jacobianInverseMatrix[2, 1] = jacobianMatrix[0, 2] * jacobianMatrix[1, 0] - jacobianMatrix[1, 2] * jacobianMatrix[0, 0];
            jacobianInverseMatrix[0, 2] = jacobianMatrix[1, 0] * jacobianMatrix[1, 1] - jacobianMatrix[2, 0] * jacobianMatrix[1, 1];

            double detj = jacobianMatrix[0, 0] * jacobianInverseMatrix[0, 0] + jacobianMatrix[0, 1] * jacobianInverseMatrix[1, 0] + jacobianMatrix[0, 2] * jacobianInverseMatrix[2, 0];

            jacobianInverseMatrix[0, 0] = jacobianInverseMatrix[0, 0] / detj;
            jacobianInverseMatrix[1, 1] = jacobianInverseMatrix[1, 1] / detj;
            jacobianInverseMatrix[2, 2] = jacobianInverseMatrix[2, 2] / detj;

            jacobianInverseMatrix[0, 1] = jacobianInverseMatrix[0, 1] / detj;
            jacobianInverseMatrix[1, 2] = jacobianInverseMatrix[1, 2] / detj;
            jacobianInverseMatrix[2, 0] = jacobianInverseMatrix[2, 0] / detj;

            jacobianInverseMatrix[1, 0] = jacobianInverseMatrix[1, 0] / detj;
            jacobianInverseMatrix[2, 1] = jacobianInverseMatrix[2, 1] / detj;
            jacobianInverseMatrix[0, 2] = jacobianInverseMatrix[0, 2] / detj;

            return new Tuple<double[,], double>(jacobianInverseMatrix, detj);
        }

        private Dictionary<int, double[]> CalculateShapeFunctionsGlobalDerivatives(Dictionary<string, double[]> dN, double[,] Jinv)
        {
            Dictionary<int, double[]> dNg = new Dictionary<int, double[]>();

            for (int i = 0; i < 8; i++)
            {
                double[] dNlocal = new double[] { dN["ksi"][i], dN["ihta"][i], dN["mhi"][i] };
                double[] dNglobal = VectorOperations.MatrixVectorProduct(Jinv, dNlocal);
                dNg.Add(i, dNglobal);
            }
            return dNg;
        }

        private double[] CalculateStrainsVector(double[,] Bmatrix)
        {
            double[] strains = VectorOperations.MatrixVectorProduct(Bmatrix, DisplacementVector);
            return strains;
        }

        private double[,] CalculateBMatrix(Dictionary<int, double[]> dNglobal)
        {
            double[,] Bmatrix = new double[6, 24];

            for (int i = 0; i < 8; i++)
            {
                Bmatrix[0, i * 3] = dNglobal[i][0];
                Bmatrix[1, i * 3 + 1] = dNglobal[i][1];
                Bmatrix[2, i * 3 + 2] = dNglobal[i][2];
                Bmatrix[3, i * 3] = dNglobal[i][1];
                Bmatrix[3, i * 3 + 1] = dNglobal[i][0];
                Bmatrix[4, i * 3 + 1] = dNglobal[i][2];
                Bmatrix[4, i * 3 + 2] = dNglobal[i][1];
                Bmatrix[5, i * 3] = dNglobal[i][2];
                Bmatrix[5, i * 3 + 2] = dNglobal[i][0];
            }
            return Bmatrix;
        }

        private double[,] CalculateStressStrainMatrix(double E, double v)
        {
            double[,] Ematrix = new double[6, 6];
            double Ehat = E / ((1.0 - 2.0 * v) * (1.0 + v));
            double G = (1.0 / 2.0) * (E / (1.0 + v));

            Ematrix[0, 0] = Ehat * (1.0 - v);
            Ematrix[0, 1] = Ehat * v;
            Ematrix[0, 2] = Ehat * v;
            Ematrix[1, 0] = Ehat * v;
            Ematrix[1, 1] = Ehat * (1.0 - v);
            Ematrix[1, 2] = Ehat * v;
            Ematrix[2, 0] = Ehat * v;
            Ematrix[2, 1] = Ehat * v;
            Ematrix[2, 2] = Ehat * (1.0 - v);
            Ematrix[3, 3] = G;
            Ematrix[4, 4] = G;
            Ematrix[5, 5] = G;
            return Ematrix;
        }

        private Tuple<double[], double[]> GaussPoints(int i, int j, int k)
        {
            double[] gaussPoints = new double[] { -1.0 / Math.Sqrt(3), 1.0 / Math.Sqrt(3) };
            double[] gaussWeights = new double[] { 1.0, 1.0 };

            double[] vectorWithPoints = new double[] { gaussPoints[i], gaussPoints[j], gaussPoints[k] };
            double[] vectorWithWeights = new double[] { gaussWeights[i], gaussWeights[j], gaussWeights[k] };
            return new Tuple<double[], double[]>(vectorWithPoints, vectorWithWeights);
        }
        //private Tuple<double[], double[]> GaussPoints(int i, int j, int k)
        //{
        //    double[] gaussPoints = new double[] { 0.0 };
        //    double[] gaussWeights = new double[] { 2.0 };

        //    double[] vectorWithPoints = new double[] { gaussPoints[i], gaussPoints[j], gaussPoints[k] };
        //    double[] vectorWithWeights = new double[] { gaussWeights[i], gaussWeights[j], gaussWeights[k] };
        //    return new Tuple<double[], double[]>(vectorWithPoints, vectorWithWeights);
        //}
        private Tuple<double[], double[]> LobattoPoints(int i, int j, int k)
        {
            double[] gaussPoints = new double[] { -1.0, 1.0 };
            double[] gaussWeights = new double[] { 1.0, 1.0 };

            double[] vectorWithPoints = new double[] { gaussPoints[i], gaussPoints[j], gaussPoints[k] };
            double[] vectorWithWeights = new double[] { gaussWeights[i], gaussWeights[j], gaussWeights[k] };
            return new Tuple<double[], double[]>(vectorWithPoints, vectorWithWeights);
        }
        public double[,] CreateGlobalStiffnessMatrix()
        {
            double[,] K = new double[24, 24];
            double[,] E = CalculateStressStrainMatrix(Properties.YoungMod, Properties.PoissonRatio);

            for (int i = 0; i < 2; i++)
            {
                for (int j = 0; j < 2; j++)
                {
                    for (int k = 0; k < 2; k++)
                    {
                        double[] gP = GaussPoints(i, j, k).Item1;
                        double[] gW = GaussPoints(i, j, k).Item2;
                        Dictionary<string, double[]> localdN = CalculateShapeFunctionsLocalDerivatives(gP);
                        double[,] J = CalculateJacobian(localdN);
                        double[,] invJ = CalculateInverseJacobian(J).Item1;
                        double detJ = CalculateInverseJacobian(J).Item2;
                        Dictionary<int, double[]> globaldN = CalculateShapeFunctionsGlobalDerivatives(localdN, invJ);
                        double[,] B = CalculateBMatrix(globaldN);
                        K = MatrixOperations.MatrixAddition(K, MatrixOperations.ScalarMatrixProductNew(detJ * gW[0] * gW[1] * gW[2],
                            MatrixOperations.MatrixProduct(MatrixOperations.Transpose(B), MatrixOperations.MatrixProduct(E, B))));

                    }
                }
            }
            return K;
        }

        public double[,] CreateMassMatrix()
        {
            double[,] M = new double[24, 24];
            for (int i = 0; i < 2; i++)
            {
                for (int j = 0; j < 2; j++)
                {
                    for (int k = 0; k < 2; k++)
                    {
                        double[] lP = LobattoPoints(i, j, k).Item1;
                        double[] lW = LobattoPoints(i, j, k).Item2;
                        Dictionary<string, double[]> localdN = CalculateShapeFunctionsLocalDerivatives(lP);
                        double[,] J = CalculateJacobian(localdN);
                        double[,] invJ = CalculateInverseJacobian(J).Item1;
                        double detJ = CalculateInverseJacobian(J).Item2;
                        double[,] N = CalculateShapeFunctionMatrix(lP[0], lP[1], lP[2]);
                        M = MatrixOperations.MatrixAddition(M, MatrixOperations.ScalarMatrixProductNew(detJ * lW[0] * lW[1] * lW[2] * Properties.Density,
                            MatrixOperations.MatrixProduct(MatrixOperations.Transpose(N), N)));
                    }
                }
            }
            return M;
        }

        public double[,] CreateDampingMatrix()
        {
            return new double[24, 24];
        }

        private double[] CalculateStressVector(double[,] E, double[] strain)
        {
            double[] stressVector = VectorOperations.MatrixVectorProduct(E, strain);
            return stressVector;
        }

        public double[] CreateInternalGlobalForcesVector()
        {
            double[] fInt = new double[24];
            double[,] E = CalculateStressStrainMatrix(Properties.YoungMod, Properties.PoissonRatio);

            for (int i = 0; i < 2; i++)
            {
                for (int j = 0; j < 2; j++)
                {
                    for (int k = 0; k < 2; k++)
                    {
                        double[] gP = GaussPoints(i, j, k).Item1;
                        double[] gW = GaussPoints(i, j, k).Item2;
                        Dictionary<string, double[]> localdN = CalculateShapeFunctionsLocalDerivatives(gP);
                        double[,] J = CalculateJacobian(localdN);
                        double[,] invJ = CalculateInverseJacobian(J).Item1;
                        double detJ = CalculateInverseJacobian(J).Item2;
                        Dictionary<int, double[]> globaldN = CalculateShapeFunctionsGlobalDerivatives(localdN, invJ);
                        double[,] B = CalculateBMatrix(globaldN);
                        double[] strainVector = CalculateStrainsVector(B);
                        double[] stressVector = CalculateStressVector(E, strainVector);
                        fInt = VectorOperations.VectorVectorAddition(fInt, VectorOperations.VectorScalarProductNew(
                            VectorOperations.MatrixVectorProduct(MatrixOperations.Transpose(B), stressVector), detJ * gW[0] * gW[1] * gW[2]));
                    }
                }
            }
            //double[,] Kstiff = CreateGlobalStiffnessMatrix();
            //double[] uDisp = DisplacementVector;
            //double[] fInt = VectorOperations.MatrixVectorProduct(Kstiff, uDisp);
            //fInt = VectorOperations.VectorScalarProductNew(fInt, 1.0);
            return fInt;
        }
    }
}