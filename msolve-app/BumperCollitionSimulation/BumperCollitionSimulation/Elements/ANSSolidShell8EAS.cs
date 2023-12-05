using System;
using System.Collections.Generic;

namespace BumperCollitionSimulation
{
    class ANSSolidShell8EAS : IElement
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
        public ANSSolidShell8EAS(IElementProperties properties, Dictionary<int, INode> nodes)
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
            Properties.EASVector = new double[7];
            Properties.EASFEnhancedVector = new double[7];
            Properties.EASLMatrix = new double[7, 24];
            Properties.EASDMatrix = new double[7, 7];
            Properties.DisplacementVectorPreviousIncrement = new double[24];
        }
        public void CalculateElementEASMatrices()
        {
            throw new Exception("This method is to be used only for EAS method elements");
        }
        public void InitializeElementEASParameters()
        {
            Properties.EASVector[0] = 0.0;
            Properties.EASVector[1] = 0.0;
            Properties.EASVector[2] = 0.0;
            Properties.EASVector[3] = 0.0;
            Properties.EASVector[4] = 0.0;
            Properties.EASVector[5] = 0.0;
            Properties.EASVector[6] = 0.0;
        }
        public void UpdateElementEASParameters(double[] totalU)
        {
            double[] deltaU = VectorOperations.VectorVectorSubtraction(totalU, Properties.DisplacementVectorPreviousIncrement);
            double[] aPrevious = Properties.EASVector;
            double[] Pe = Properties.EASFEnhancedVector;
            double[,] Le = Properties.EASLMatrix;
            double[,] DeInv = Properties.EASDMatrix;
            Properties.EASVector = VectorOperations.VectorVectorAddition(aPrevious,
                VectorOperations.VectorScalarProductNew(
                VectorOperations.MatrixVectorProduct(DeInv, VectorOperations.VectorVectorAddition(
                VectorOperations.MatrixVectorProduct(Le, deltaU), Pe)), -1.0));
            Properties.DisplacementVectorPreviousIncrement = totalU;
        }
        public void StoreElementFinalStepDisplacementVector(double[] solutionVector)
        {
            Properties.DisplacementVectorPreviousIncrement = solutionVector;
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
            throw new Exception("Method to be added later in element ANSSolidShell8LEAS1RI");
        }
        public List<double[]> GetStrainVector()
        {
            throw new Exception("Method to be added later in element ANSSolidShell8LEAS1RI");
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
            throw new Exception("Method to be added later in element ANSSolidShell8LEAS1RI");

        }
        public List<double[]> GetStrainFromElementsNodes()
        {
            throw new Exception("Method to be added later in element ANSSolidShell8LEAS1RI");

        }
        public List<double[]> GetStressFromElements(List<double[]> parametricCoordinatesVector)
        {
            throw new Exception("Method to be added later in element ANSSolidShell8LEAS1RI");
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
            double[] initialCoordinates = new double[24];
            for (int i = 1; i <= 8; i++)
            {
                initialCoordinates[3 * i - 3] = Nodes[i].XCoordinate;
                initialCoordinates[3 * i - 2] = Nodes[i].YCoordinate;
                initialCoordinates[3 * i - 1] = Nodes[i].ZCoordinate;
            }
            return initialCoordinates;
        }
        Dictionary<int, double[]> NodalCartesianCoordinatesCurrent(double[] displacementVector)
        {
            Dictionary<int, double[]> updatedCoor = new Dictionary<int, double[]>();
            for (int i = 0; i < 8; i++)
            {
                double[] nodalPositionVector = new double[]
                {
                    Nodes[i + 1].XCoordinate + displacementVector[3 * i],
                    Nodes[i + 1].YCoordinate + displacementVector[3 * i + 1],
                    Nodes[i + 1].ZCoordinate + displacementVector[3 * i + 2]
                };
                updatedCoor.Add(i, nodalPositionVector);
            }
            return updatedCoor;
        }
        Dictionary<int, double[]> NodalCartesianCoordinatesInitial()
        {
            Dictionary<int, double[]> updatedCoor = new Dictionary<int, double[]>();
            for (int i = 0; i < 8; i++)
            {
                double[] nodalPositionVector = new double[]
                {
                    Nodes[i + 1].XCoordinate,
                    Nodes[i + 1].YCoordinate,
                    Nodes[i + 1].ZCoordinate
                };
                updatedCoor.Add(i, nodalPositionVector);
            }
            return updatedCoor;
        }
        Dictionary<int, double[]> NodalNaturalCoordinates()
        {
            Dictionary<int, double[]> naturalCoor = new Dictionary<int, double[]>();
            naturalCoor.Add(0, new double[] { 1.0, 1.0, 1.0 });
            naturalCoor.Add(1, new double[] { -1.0, 1.0, 1.0 });
            naturalCoor.Add(2, new double[] { -1.0, -1.0, 1.0 });
            naturalCoor.Add(3, new double[] { 1.0, -1.0, 1.0 });
            naturalCoor.Add(4, new double[] { 1.0, 1.0, -1.0 });
            naturalCoor.Add(5, new double[] { -1.0, 1.0, -1.0 });
            naturalCoor.Add(6, new double[] { -1.0, -1.0, -1.0 });
            naturalCoor.Add(7, new double[] { 1.0, -1.0, -1.0 });
            return naturalCoor;
        }
        Dictionary<int, double[]> PositionVectors()
        {
            Dictionary<int, double[]> positionVectors = new Dictionary<int, double[]>();
            positionVectors.Add(0, new double[] { 1.0, 1.0, -1.0, -1.0, -1.0, -1.0, 1.0, 1.0 });
            positionVectors.Add(1, new double[] { 1.0, -1.0, -1.0, 1.0, -1.0, 1.0, 1.0, -1.0 });
            positionVectors.Add(2, new double[] { 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0 });
            positionVectors.Add(3, new double[] { 1.0, -1.0, 1.0, -1.0, -1.0, 1.0, -1.0, 1.0 });
            return positionVectors;
        }
        private Dictionary<int, double> CalculateShapeFunctions(double ksi, double ihta, double mhi)
        {
            Dictionary<int, double> shapeFunctions = new Dictionary<int, double>();
            double N1 = 1.0 / 8.0 * (1 + ksi) * (1 + ihta) * (1 + mhi); shapeFunctions.Add(1, N1);
            double N2 = 1.0 / 8.0 * (1 - ksi) * (1 + ihta) * (1 + mhi); shapeFunctions.Add(2, N2);
            double N3 = 1.0 / 8.0 * (1 - ksi) * (1 - ihta) * (1 + mhi); shapeFunctions.Add(3, N3);
            double N4 = 1.0 / 8.0 * (1 + ksi) * (1 - ihta) * (1 + mhi); shapeFunctions.Add(4, N4);
            double N5 = 1.0 / 8.0 * (1 + ksi) * (1 + ihta) * (1 - mhi); shapeFunctions.Add(5, N5);
            double N6 = 1.0 / 8.0 * (1 - ksi) * (1 + ihta) * (1 - mhi); shapeFunctions.Add(6, N6);
            double N7 = 1.0 / 8.0 * (1 - ksi) * (1 - ihta) * (1 - mhi); shapeFunctions.Add(7, N7);
            double N8 = 1.0 / 8.0 * (1 + ksi) * (1 - ihta) * (1 - mhi); shapeFunctions.Add(8, N8);
            return shapeFunctions;
        }
        private double[,] CalculateShapeFunctionMatrix(double ksi, double ihta, double zita)
        {
            double N1 = 1.0 / 8.0 * (1 + ksi) * (1 + ihta) * (1 + zita);
            double N2 = 1.0 / 8.0 * (1 - ksi) * (1 + ihta) * (1 + zita);
            double N3 = 1.0 / 8.0 * (1 - ksi) * (1 - ihta) * (1 + zita);
            double N4 = 1.0 / 8.0 * (1 + ksi) * (1 - ihta) * (1 + zita);
            double N5 = 1.0 / 8.0 * (1 + ksi) * (1 + ihta) * (1 - zita);
            double N6 = 1.0 / 8.0 * (1 - ksi) * (1 + ihta) * (1 - zita);
            double N7 = 1.0 / 8.0 * (1 - ksi) * (1 - ihta) * (1 - zita);
            double N8 = 1.0 / 8.0 * (1 + ksi) * (1 - ihta) * (1 - zita);
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
                (1.0/8.0*(1+ihta)*(1+mhi)),
                (-1.0/8.0*(1+ihta)*(1+mhi)),
                (-1.0/8.0*(1-ihta)*(1+mhi)),
                (1.0/8.0*(1-ihta)*(1+mhi)),
                (1.0/8.0*(1+ihta)*(1-mhi)),
                (-1.0/8.0*(1+ihta)*(1-mhi)),
                (-1.0/8.0*(1-ihta)*(1-mhi)),
                (1.0/8.0*(1-ihta)*(1-mhi))
            };

            double[] dN_ihta = new double[]
            {
                (1.0/8.0*(1+ksi)*(1+mhi)),
                (1.0/8.0*(1-ksi)*(1+mhi)),
                (-1.0/8.0*(1-ksi)*(1+mhi)),
                (-1.0/8.0*(1+ksi)*(1+mhi)),
                (1.0/8.0*(1+ksi)*(1-mhi)),
                (1.0/8.0*(1-ksi)*(1-mhi)),
                (-1.0/8.0*(1-ksi)*(1-mhi)),
                (-1.0/8.0*(1+ksi)*(1-mhi))
            };

            double[] dN_mhi = new double[]
            {
                (1.0/8.0*(1+ksi)*(1+ihta)),
                (1.0/8.0*(1-ksi)*(1+ihta)),
                (1.0/8.0*(1-ksi)*(1-ihta)),
                (1.0/8.0*(1+ksi)*(1-ihta)),
                (-1.0/8.0*(1+ksi)*(1+ihta)),
                (-1.0/8.0*(1-ksi)*(1+ihta)),
                (-1.0/8.0*(1-ksi)*(1-ihta)),
                (-1.0/8.0*(1+ksi)*(1-ihta))
            };

            Dictionary<string, double[]> dN = new Dictionary<string, double[]>();
            dN.Add("ksi", dN_ksi);
            dN.Add("ihta", dN_ihta);
            dN.Add("mhi", dN_mhi);
            return dN;
        }
        private double CalculateShapeFunctionsLocalDerivatives(double[] naturalCoordinates, int parametricAxis,
            int node)
        {
            double ksi = naturalCoordinates[0];
            double ihta = naturalCoordinates[1];
            double mhi = naturalCoordinates[2];
            if(parametricAxis == 0)
            {
                double dN_ksi = new double();
                switch (node)
                {
                    case 0:
                        dN_ksi = 1.0 / 8.0 * (1 + ihta) * (1 + mhi);
                        break;
                    case 1:
                        dN_ksi = -1.0 / 8.0 * (1 + ihta) * (1 + mhi);
                        break;
                    case 2:
                        dN_ksi = -1.0 / 8.0 * (1 - ihta) * (1 + mhi);
                        break;
                    case 3:
                        dN_ksi = 1.0 / 8.0 * (1 - ihta) * (1 + mhi);
                        break;
                    case 4:
                        dN_ksi = 1.0 / 8.0 * (1 + ihta) * (1 - mhi);
                        break;
                    case 5:
                        dN_ksi = -1.0 / 8.0 * (1 + ihta) * (1 - mhi);
                        break;
                    case 6:
                        dN_ksi = -1.0 / 8.0 * (1 - ihta) * (1 - mhi);
                        break;
                    case 7:
                        dN_ksi = 1.0 / 8.0 * (1 - ihta) * (1 - mhi);
                        break;
                }
                return dN_ksi;
            }
            else if (parametricAxis == 1)
            {
                double dN_ihta = new double();
                switch (node)
                {
                    case 0:
                        dN_ihta = 1.0 / 8.0 * (1 + ksi) * (1 + mhi);
                        break;
                    case 1:
                        dN_ihta = 1.0 / 8.0 * (1 - ksi) * (1 + mhi);
                        break;
                    case 2:
                        dN_ihta = -1.0 / 8.0 * (1 - ksi) * (1 + mhi);
                        break;
                    case 3:
                        dN_ihta = -1.0 / 8.0 * (1 + ksi) * (1 + mhi);
                        break;
                    case 4:
                        dN_ihta = 1.0 / 8.0 * (1 + ksi) * (1 - mhi);
                        break;
                    case 5:
                        dN_ihta = 1.0 / 8.0 * (1 - ksi) * (1 - mhi);
                        break;
                    case 6:
                        dN_ihta = -1.0 / 8.0 * (1 - ksi) * (1 - mhi);
                        break;
                    case 7:
                        dN_ihta = -1.0 / 8.0 * (1 + ksi) * (1 - mhi);
                        break;
                }
                return dN_ihta;
            }
            else
            {

                double dN_mhi = new double();
                switch (node)
                {
                    case 0:
                        dN_mhi = 1.0 / 8.0 * (1 + ksi) * (1 + ihta);
                        break;
                    case 1:
                        dN_mhi = 1.0 / 8.0 * (1 - ksi) * (1 + ihta);
                        break;
                    case 2:
                        dN_mhi = 1.0 / 8.0 * (1 - ksi) * (1 - ihta);
                        break;
                    case 3:
                        dN_mhi = 1.0 / 8.0 * (1 + ksi) * (1 - ihta);
                        break;
                    case 4:
                        dN_mhi = -1.0 / 8.0 * (1 + ksi) * (1 + ihta);
                        break;
                    case 5:
                        dN_mhi = -1.0 / 8.0 * (1 - ksi) * (1 + ihta);
                        break;
                    case 6:
                        dN_mhi = -1.0 / 8.0 * (1 - ksi) * (1 - ihta);
                        break;
                    case 7:
                        dN_mhi = -1.0 / 8.0 * (1 + ksi) * (1 - ihta);
                        break;
                }
                return dN_mhi;
            }
        }
        private double[,] CalculateUpdatedJacobianMatrix(double[] xUpdated, Dictionary<string, double[]> dN)
        {
            double[,] jacobianMatrix = new double[3, 3];

            int k = 0;
            for (int i = 0; i < 8; i++)
            {
                jacobianMatrix[0, 0] = jacobianMatrix[0, 0] + xUpdated[k] * dN["ksi"][i];
                k = k + 3;
            }
            k = 1;
            for (int i = 0; i < 8; i++)
            {
                jacobianMatrix[0, 1] = jacobianMatrix[0, 1] + xUpdated[k] * dN["ksi"][i];
                k = k + 3;
            }
            k = 2;
            for (int i = 0; i < 8; i++)
            {
                jacobianMatrix[0, 2] = jacobianMatrix[0, 2] + xUpdated[k] * dN["ksi"][i];
                k = k + 3;
            }

            k = 0;
            for (int i = 0; i < 8; i++)
            {
                jacobianMatrix[1, 0] = jacobianMatrix[1, 0] + xUpdated[k] * dN["ihta"][i];
                k = k + 3;
            }
            k = 1;
            for (int i = 0; i < 8; i++)
            {
                jacobianMatrix[1, 1] = jacobianMatrix[1, 1] + xUpdated[k] * dN["ihta"][i];
                k = k + 3;
            }
            k = 2;
            for (int i = 0; i < 8; i++)
            {
                jacobianMatrix[1, 2] = jacobianMatrix[1, 2] + xUpdated[k] * dN["ihta"][i];
                k = k + 3;
            }

            k = 0;
            for (int i = 0; i < 8; i++)
            {
                jacobianMatrix[2, 0] = jacobianMatrix[2, 0] + xUpdated[k] * dN["mhi"][i];
                k = k + 3;
            }
            k = 1;
            for (int i = 0; i < 8; i++)
            {
                jacobianMatrix[2, 1] = jacobianMatrix[2, 1] + xUpdated[k] * dN["mhi"][i];
                k = k + 3;
            }
            k = 2;
            for (int i = 0; i < 8; i++)
            {
                jacobianMatrix[2, 2] = jacobianMatrix[2, 2] + xUpdated[k] * dN["mhi"][i];
                k = k + 3;
            }

            return jacobianMatrix;
        }
        private double[,] CalculateJacobian(Dictionary<string, double[]> dN)
        {
            double[,] jacobianMatrix = new double[3, 3];
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
        private double CalculateJacobianDet(double[,] jacobianMatrix)
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
            return detj;
        }
        private double[,] CalculateTransformationMatrix(double[,] jacobianInverseMatrix)
        {
            double[,] T = new double[6, 6];
            T[0, 0] = jacobianInverseMatrix[0, 0] * jacobianInverseMatrix[0, 0];
            T[0, 1] = jacobianInverseMatrix[0, 1] * jacobianInverseMatrix[0, 1];
            T[0, 2] = jacobianInverseMatrix[0, 2] * jacobianInverseMatrix[0, 2];
            T[0, 3] = jacobianInverseMatrix[0, 0] * jacobianInverseMatrix[0, 1];
            T[0, 4] = jacobianInverseMatrix[0, 1] * jacobianInverseMatrix[0, 2];
            T[0, 5] = jacobianInverseMatrix[0, 2] * jacobianInverseMatrix[0, 0];

            T[1, 0] = jacobianInverseMatrix[1, 0] * jacobianInverseMatrix[1, 0];
            T[1, 1] = jacobianInverseMatrix[1, 1] * jacobianInverseMatrix[1, 1];
            T[1, 2] = jacobianInverseMatrix[1, 2] * jacobianInverseMatrix[1, 2];
            T[1, 3] = jacobianInverseMatrix[1, 0] * jacobianInverseMatrix[1, 1];
            T[1, 4] = jacobianInverseMatrix[1, 1] * jacobianInverseMatrix[1, 2];
            T[1, 5] = jacobianInverseMatrix[1, 2] * jacobianInverseMatrix[1, 0];

            T[2, 0] = jacobianInverseMatrix[2, 0] * jacobianInverseMatrix[2, 0];
            T[2, 1] = jacobianInverseMatrix[2, 1] * jacobianInverseMatrix[2, 1];
            T[2, 2] = jacobianInverseMatrix[2, 2] * jacobianInverseMatrix[2, 2];
            T[2, 3] = jacobianInverseMatrix[2, 0] * jacobianInverseMatrix[2, 1];
            T[2, 4] = jacobianInverseMatrix[2, 1] * jacobianInverseMatrix[2, 2];
            T[2, 5] = jacobianInverseMatrix[2, 2] * jacobianInverseMatrix[2, 0];

            T[3, 0] = 2 * jacobianInverseMatrix[0, 0] * jacobianInverseMatrix[1, 0];
            T[3, 1] = 2 * jacobianInverseMatrix[0, 1] * jacobianInverseMatrix[1, 1];
            T[3, 2] = 2 * jacobianInverseMatrix[0, 2] * jacobianInverseMatrix[1, 2];

            T[3, 3] = jacobianInverseMatrix[0, 0] * jacobianInverseMatrix[1, 1] +
                                jacobianInverseMatrix[0, 1] * jacobianInverseMatrix[1, 0];

            T[3, 4] = jacobianInverseMatrix[0, 1] * jacobianInverseMatrix[1, 2] +
                                jacobianInverseMatrix[0, 2] * jacobianInverseMatrix[1, 1];

            T[3, 5] = jacobianInverseMatrix[0, 2] * jacobianInverseMatrix[1, 0] +
                                jacobianInverseMatrix[0, 0] * jacobianInverseMatrix[2, 1];

            T[4, 0] = 2 * jacobianInverseMatrix[1, 0] * jacobianInverseMatrix[2, 0];
            T[4, 1] = 2 * jacobianInverseMatrix[1, 1] * jacobianInverseMatrix[2, 1];
            T[4, 2] = 2 * jacobianInverseMatrix[1, 2] * jacobianInverseMatrix[2, 2];

            T[4, 3] = jacobianInverseMatrix[1, 0] * jacobianInverseMatrix[2, 1] +
                                jacobianInverseMatrix[1, 1] * jacobianInverseMatrix[2, 0];

            T[4, 4] = jacobianInverseMatrix[1, 1] * jacobianInverseMatrix[2, 2] +
                                jacobianInverseMatrix[1, 2] * jacobianInverseMatrix[2, 1];

            T[4, 5] = jacobianInverseMatrix[1, 2] * jacobianInverseMatrix[2, 0] +
                                jacobianInverseMatrix[1, 0] * jacobianInverseMatrix[2, 2];

            T[5, 0] = 2 * jacobianInverseMatrix[2, 0] * jacobianInverseMatrix[0, 0];
            T[5, 1] = 2 * jacobianInverseMatrix[2, 1] * jacobianInverseMatrix[0, 1];
            T[5, 2] = 2 * jacobianInverseMatrix[2, 2] * jacobianInverseMatrix[0, 2];

            T[5, 3] = jacobianInverseMatrix[2, 0] * jacobianInverseMatrix[0, 1] +
                                jacobianInverseMatrix[2, 1] * jacobianInverseMatrix[0, 0];

            T[5, 4] = jacobianInverseMatrix[2, 1] * jacobianInverseMatrix[0, 2] +
                                jacobianInverseMatrix[2, 2] * jacobianInverseMatrix[0, 1];

            T[5, 5] = jacobianInverseMatrix[2, 2] * jacobianInverseMatrix[0, 0] +
                                jacobianInverseMatrix[2, 0] * jacobianInverseMatrix[0, 2];
            return T;
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
        private double[] CalculateJacobianPart1(double[] parametricCoordinates)
        {
            double[] jacobianPart1 = new double[3];
            double[] xNodalInitial = InitialNodalCoordinates();

            int k = 0;
            for (int i = 0; i < 8; i++)
            {
                jacobianPart1[0] = jacobianPart1[0] + xNodalInitial[k] * CalculateShapeFunctionsLocalDerivatives(parametricCoordinates, 0, i); 
                k = k + 3;
            }
            k = 1;
            for (int i = 0; i < 8; i++)
            {
                jacobianPart1[1] = jacobianPart1[1] + xNodalInitial[k] * CalculateShapeFunctionsLocalDerivatives(parametricCoordinates, 0, i);
                k = k + 3;
            }
            k = 2;
            for (int i = 0; i < 8; i++)
            {
                jacobianPart1[2] = jacobianPart1[2] + xNodalInitial[k] * CalculateShapeFunctionsLocalDerivatives(parametricCoordinates, 0, i);
                k = k + 3;
            }
            return jacobianPart1;
        }
        private double[] CalculateJacobianPart2(double[] parametricCoordinates)
        {
            double[] jacobianPart2 = new double[3];
            double[] xNodalInitial = InitialNodalCoordinates();

            int k = 0;
            for (int i = 0; i < 8; i++)
            {
                jacobianPart2[0] = jacobianPart2[0] + xNodalInitial[k] * CalculateShapeFunctionsLocalDerivatives(parametricCoordinates, 1, i);
                k = k + 3;
            }
            k = 1;
            for (int i = 0; i < 8; i++)
            {
                jacobianPart2[1] = jacobianPart2[1] + xNodalInitial[k] * CalculateShapeFunctionsLocalDerivatives(parametricCoordinates, 1, i);
                k = k + 3;
            }
            k = 2;
            for (int i = 0; i < 8; i++)
            {
                jacobianPart2[2] = jacobianPart2[2] + xNodalInitial[k] * CalculateShapeFunctionsLocalDerivatives(parametricCoordinates, 1, i);
                k = k + 3;
            }
            return jacobianPart2;
        }
        private double[] CalculateJacobianPart3(double[] parametricCoordinates)
        {
            double[] jacobianPart3 = new double[3];
            double[] xNodalInitial = InitialNodalCoordinates();

            int k = 0;
            for (int i = 0; i < 8; i++)
            {
                jacobianPart3[0] = jacobianPart3[0] + xNodalInitial[k] * CalculateShapeFunctionsLocalDerivatives(parametricCoordinates, 2, i);
                k = k + 3;
            }
            k = 1;
            for (int i = 0; i < 8; i++)
            {
                jacobianPart3[1] = jacobianPart3[1] + xNodalInitial[k] * CalculateShapeFunctionsLocalDerivatives(parametricCoordinates, 2, i);
                k = k + 3;
            }
            k = 2;
            for (int i = 0; i < 8; i++)
            {
                jacobianPart3[2] = jacobianPart3[2] + xNodalInitial[k] * CalculateShapeFunctionsLocalDerivatives(parametricCoordinates, 2, i);
                k = k + 3;
            }
            return jacobianPart3;
        }
        private double[,] CalculateBMatrix(Dictionary<string, double[]> dNlocal, double[,] jacobianMatrix,
            double[] parametricCoordinates, double[,] transformationMatrix)
        {
            double[,] BmatrixLocal = new double[6, 24];
            //-----------------------------------------------------------------
            double[] a = new double[] { 1d, 1d, 0d };
            double[] b = new double[] { -1d, 1d, 0d };
            double[] c = new double[] { -1d, -1d, 0d };
            double[] d = new double[] { 1d, -1d, 0d };
            //-----------------------------------------------------------------
            double[] e = new double[] { 1d, 0d, 1d };
            double[] f = new double[] { -1d, 0d, 1d };
            double[] g = new double[] { -1d, 0d, -1d };
            double[] h = new double[] { 1d, 0d, -1d };
            //-----------------------------------------------------------------
            double[] j = new double[] { 0d, 1d, 1d };
            double[] k = new double[] { 0d, -1d, 1d };
            double[] l = new double[] { 0d, -1d, -1d };
            double[] m = new double[] { 0d, 1d, -1d };
            for (int i = 0; i < 8; i++)
            {
                BmatrixLocal[0, i * 3] = dNlocal["ksi"][i] * jacobianMatrix[0, 0];
                BmatrixLocal[0, i * 3 + 1] = dNlocal["ksi"][i] * jacobianMatrix[0, 1];
                BmatrixLocal[0, i * 3 + 2] = dNlocal["ksi"][i] * jacobianMatrix[0, 2];

                BmatrixLocal[1, i * 3] = dNlocal["ihta"][i] * jacobianMatrix[1, 0];
                BmatrixLocal[1, i * 3 + 1] = dNlocal["ihta"][i] * jacobianMatrix[1, 1];
                BmatrixLocal[1, i * 3 + 2] = dNlocal["ihta"][i] * jacobianMatrix[1, 2];

                //εζζ Assumed strain interpolation from points A, B, C, D
                BmatrixLocal[2, i * 3] = (1.0 / 4.0) * (1.0 + parametricCoordinates[0]) * (1.0 + parametricCoordinates[1]) *
                                CalculateShapeFunctionsLocalDerivatives(a, 2, i) * CalculateJacobianPart3(a)[0] +
                                (1.0 / 4.0) * (1.0 - parametricCoordinates[0]) * (1.0 + parametricCoordinates[1]) *
                                CalculateShapeFunctionsLocalDerivatives(b, 2, i) * CalculateJacobianPart3(b)[0] +
                                (1.0 / 4.0) * (1.0 - parametricCoordinates[0]) * (1.0 - parametricCoordinates[1]) *
                                CalculateShapeFunctionsLocalDerivatives(c, 2, i) * CalculateJacobianPart3(c)[0] +
                                (1.0 / 4.0) * (1.0 + parametricCoordinates[0]) * (1.0 - parametricCoordinates[1]) *
                                CalculateShapeFunctionsLocalDerivatives(d, 2, i) * CalculateJacobianPart3(d)[0];
                BmatrixLocal[2, i * 3 + 1] = (1.0 / 4.0) * (1.0 + parametricCoordinates[0]) * (1.0 + parametricCoordinates[1]) *
                                CalculateShapeFunctionsLocalDerivatives(a, 2, i) * CalculateJacobianPart3(a)[1] +
                                (1.0 / 4.0) * (1.0 - parametricCoordinates[0]) * (1.0 + parametricCoordinates[1]) *
                                CalculateShapeFunctionsLocalDerivatives(b, 2, i) * CalculateJacobianPart3(b)[1] +
                                (1.0 / 4.0) * (1.0 - parametricCoordinates[0]) * (1.0 - parametricCoordinates[1]) *
                                CalculateShapeFunctionsLocalDerivatives(c, 2, i) * CalculateJacobianPart3(c)[1] +
                                (1.0 / 4.0) * (1.0 + parametricCoordinates[0]) * (1.0 - parametricCoordinates[1]) *
                                CalculateShapeFunctionsLocalDerivatives(d, 2, i) * CalculateJacobianPart3(d)[1];
                BmatrixLocal[2, i * 3 + 2] = (1.0 / 4.0) * (1.0 + parametricCoordinates[0]) * (1.0 + parametricCoordinates[1]) *
                                CalculateShapeFunctionsLocalDerivatives(a, 2, i) * CalculateJacobianPart3(a)[2] +
                                (1.0 / 4.0) * (1.0 - parametricCoordinates[0]) * (1.0 + parametricCoordinates[1]) *
                                CalculateShapeFunctionsLocalDerivatives(b, 2, i) * CalculateJacobianPart3(b)[2] +
                                (1.0 / 4.0) * (1.0 - parametricCoordinates[0]) * (1.0 - parametricCoordinates[1]) *
                                CalculateShapeFunctionsLocalDerivatives(c, 2, i) * CalculateJacobianPart3(c)[2] +
                                (1.0 / 4.0) * (1.0 + parametricCoordinates[0]) * (1.0 - parametricCoordinates[1]) *
                                CalculateShapeFunctionsLocalDerivatives(d, 2, i) * CalculateJacobianPart3(d)[2];
                //-----------------------------------------------------------------
                BmatrixLocal[3, i * 3] = dNlocal["ihta"][i] * jacobianMatrix[0, 0] +
                    dNlocal["ksi"][i] * jacobianMatrix[1, 0];
                BmatrixLocal[3, i * 3 + 1] = dNlocal["ihta"][i] * jacobianMatrix[0, 1] +
                    dNlocal["ksi"][i] * jacobianMatrix[1, 1];
                BmatrixLocal[3, i * 3 + 2] = dNlocal["ihta"][i] * jacobianMatrix[0, 2] +
                    dNlocal["ksi"][i] * jacobianMatrix[1, 2];
                //εηζ Assumed strain interpolation from points E,F,G,H
                BmatrixLocal[4, i * 3] = (1.0 / 4.0) * (1.0 + parametricCoordinates[0]) * (1.0 + parametricCoordinates[2]) *
                                (CalculateShapeFunctionsLocalDerivatives(e, 2, i) * CalculateJacobianPart2(e)[0] +
                                CalculateShapeFunctionsLocalDerivatives(e, 1, i) * CalculateJacobianPart3(e)[0]) +
                                (1.0 / 4.0) * (1.0 - parametricCoordinates[0]) * (1.0 + parametricCoordinates[2]) *
                                (CalculateShapeFunctionsLocalDerivatives(f, 2, i) * CalculateJacobianPart2(f)[0] +
                                CalculateShapeFunctionsLocalDerivatives(f, 1, i) * CalculateJacobianPart3(f)[0]) +
                                (1.0 / 4.0) * (1.0 - parametricCoordinates[0]) * (1.0 - parametricCoordinates[2]) *
                                (CalculateShapeFunctionsLocalDerivatives(g, 2, i) * CalculateJacobianPart2(g)[0] +
                                CalculateShapeFunctionsLocalDerivatives(g, 1, i) * CalculateJacobianPart3(g)[0]) +
                                (1.0 / 4.0) * (1.0 + parametricCoordinates[0]) * (1.0 - parametricCoordinates[2]) *
                                (CalculateShapeFunctionsLocalDerivatives(h, 2, i) * CalculateJacobianPart2(h)[0] +
                                CalculateShapeFunctionsLocalDerivatives(h, 1, i) * CalculateJacobianPart3(h)[0]);

                BmatrixLocal[4, i * 3 + 1] = (1.0 / 4.0) * (1.0 + parametricCoordinates[0]) * (1.0 + parametricCoordinates[2]) *
                                (CalculateShapeFunctionsLocalDerivatives(e, 2, i) * CalculateJacobianPart2(e)[1] +
                                CalculateShapeFunctionsLocalDerivatives(e, 1, i) * CalculateJacobianPart3(e)[1]) +
                                (1.0 / 4.0) * (1.0 - parametricCoordinates[0]) * (1.0 + parametricCoordinates[2]) *
                                (CalculateShapeFunctionsLocalDerivatives(f, 2, i) * CalculateJacobianPart2(f)[1] +
                                CalculateShapeFunctionsLocalDerivatives(f, 1, i) * CalculateJacobianPart3(f)[1]) +
                                (1.0 / 4.0) * (1.0 - parametricCoordinates[0]) * (1.0 - parametricCoordinates[2]) *
                                (CalculateShapeFunctionsLocalDerivatives(g, 2, i) * CalculateJacobianPart2(g)[1] +
                                CalculateShapeFunctionsLocalDerivatives(g, 1, i) * CalculateJacobianPart3(g)[1]) +
                                (1.0 / 4.0) * (1.0 + parametricCoordinates[0]) * (1.0 - parametricCoordinates[2]) *
                                (CalculateShapeFunctionsLocalDerivatives(h, 2, i) * CalculateJacobianPart2(h)[1] +
                                CalculateShapeFunctionsLocalDerivatives(h, 1, i) * CalculateJacobianPart3(h)[1]);

                BmatrixLocal[4, i * 3 + 2] = (1.0 / 4.0) * (1.0 + parametricCoordinates[0]) * (1.0 + parametricCoordinates[2]) *
                                (CalculateShapeFunctionsLocalDerivatives(e, 2, i) * CalculateJacobianPart2(e)[2] +
                                CalculateShapeFunctionsLocalDerivatives(e, 1, i) * CalculateJacobianPart3(e)[2]) +
                                (1.0 / 4.0) * (1.0 - parametricCoordinates[0]) * (1.0 + parametricCoordinates[2]) *
                                (CalculateShapeFunctionsLocalDerivatives(f, 2, i) * CalculateJacobianPart2(f)[2] +
                                CalculateShapeFunctionsLocalDerivatives(f, 1, i) * CalculateJacobianPart3(f)[2]) +
                                (1.0 / 4.0) * (1.0 - parametricCoordinates[0]) * (1.0 - parametricCoordinates[2]) *
                                (CalculateShapeFunctionsLocalDerivatives(g, 2, i) * CalculateJacobianPart2(g)[2] +
                                CalculateShapeFunctionsLocalDerivatives(g, 1, i) * CalculateJacobianPart3(g)[2]) +
                                (1.0 / 4.0) * (1.0 + parametricCoordinates[0]) * (1.0 - parametricCoordinates[2]) *
                                (CalculateShapeFunctionsLocalDerivatives(h, 2, i) * CalculateJacobianPart2(h)[2] +
                                CalculateShapeFunctionsLocalDerivatives(h, 1, i) * CalculateJacobianPart3(h)[2]);
                //-----------------------------------------------------------------
                //εζξ Assumed strain interpolation from points J, K, L, M
                BmatrixLocal[5, i * 3] = (1.0 / 4.0) * (1.0 + parametricCoordinates[1]) * (1.0 + parametricCoordinates[2]) *
                                (CalculateShapeFunctionsLocalDerivatives(j, 0, i) * CalculateJacobianPart3(j)[0] +
                                CalculateShapeFunctionsLocalDerivatives(j, 2, i) * CalculateJacobianPart1(j)[0]) +
                                (1.0 / 4.0) * (1.0 - parametricCoordinates[1]) * (1.0 + parametricCoordinates[2]) *
                                (CalculateShapeFunctionsLocalDerivatives(k, 0, i) * CalculateJacobianPart3(k)[0] +
                                CalculateShapeFunctionsLocalDerivatives(k, 2, i) * CalculateJacobianPart1(k)[0]) +
                                (1.0 / 4.0) * (1.0 - parametricCoordinates[1]) * (1.0 - parametricCoordinates[2]) *
                                (CalculateShapeFunctionsLocalDerivatives(l, 0, i) * CalculateJacobianPart3(l)[0] +
                                CalculateShapeFunctionsLocalDerivatives(l, 2, i) * CalculateJacobianPart1(l)[0]) +
                                (1.0 / 4.0) * (1.0 + parametricCoordinates[1]) * (1.0 - parametricCoordinates[2]) *
                                (CalculateShapeFunctionsLocalDerivatives(m, 0, i) * CalculateJacobianPart3(m)[0] +
                                CalculateShapeFunctionsLocalDerivatives(m, 2, i) * CalculateJacobianPart1(m)[0]);

                BmatrixLocal[5, i * 3 + 1] = (1.0 / 4.0) * (1.0 + parametricCoordinates[1]) * (1.0 + parametricCoordinates[2]) *
                                (CalculateShapeFunctionsLocalDerivatives(j, 0, i) * CalculateJacobianPart3(j)[1] +
                                CalculateShapeFunctionsLocalDerivatives(j, 2, i) * CalculateJacobianPart1(j)[1]) +
                                (1.0 / 4.0) * (1.0 - parametricCoordinates[1]) * (1.0 + parametricCoordinates[2]) *
                                (CalculateShapeFunctionsLocalDerivatives(k, 0, i) * CalculateJacobianPart3(k)[1] +
                                CalculateShapeFunctionsLocalDerivatives(k, 2, i) * CalculateJacobianPart1(k)[1]) +
                                (1.0 / 4.0) * (1.0 - parametricCoordinates[1]) * (1.0 - parametricCoordinates[2]) *
                                (CalculateShapeFunctionsLocalDerivatives(l, 0, i) * CalculateJacobianPart3(l)[1] +
                                CalculateShapeFunctionsLocalDerivatives(l, 2, i) * CalculateJacobianPart1(l)[1]) +
                                (1.0 / 4.0) * (1.0 + parametricCoordinates[1]) * (1.0 - parametricCoordinates[2]) *
                                (CalculateShapeFunctionsLocalDerivatives(m, 0, i) * CalculateJacobianPart3(m)[1] +
                                CalculateShapeFunctionsLocalDerivatives(m, 2, i) * CalculateJacobianPart1(m)[1]);

                BmatrixLocal[5, i * 3 + 2] = (1.0 / 4.0) * (1.0 + parametricCoordinates[1]) * (1.0 + parametricCoordinates[2]) *
                                (CalculateShapeFunctionsLocalDerivatives(j, 0, i) * CalculateJacobianPart3(j)[2] +
                                CalculateShapeFunctionsLocalDerivatives(j, 2, i) * CalculateJacobianPart1(j)[2]) +
                                (1.0 / 4.0) * (1.0 - parametricCoordinates[1]) * (1.0 + parametricCoordinates[2]) *
                                (CalculateShapeFunctionsLocalDerivatives(k, 0, i) * CalculateJacobianPart3(k)[2] +
                                CalculateShapeFunctionsLocalDerivatives(k, 2, i) * CalculateJacobianPart1(k)[2]) +
                                (1.0 / 4.0) * (1.0 - parametricCoordinates[1]) * (1.0 - parametricCoordinates[2]) *
                                (CalculateShapeFunctionsLocalDerivatives(l, 0, i) * CalculateJacobianPart3(l)[2] +
                                CalculateShapeFunctionsLocalDerivatives(l, 2, i) * CalculateJacobianPart1(l)[2]) +
                                (1.0 / 4.0) * (1.0 + parametricCoordinates[1]) * (1.0 - parametricCoordinates[2]) *
                                (CalculateShapeFunctionsLocalDerivatives(m, 0, i) * CalculateJacobianPart3(m)[2] +
                                CalculateShapeFunctionsLocalDerivatives(m, 2, i) * CalculateJacobianPart1(m)[2]);
                //-----------------------------------------------------------------
            }
            var BmatrixGlobal = MatrixOperations.MatrixProduct(transformationMatrix, BmatrixLocal);
            return BmatrixGlobal;
        }
        private double[] CalculateStrainsVector(double[,] Bmatrix)
        {
            double[] strains = VectorOperations.MatrixVectorProduct(Bmatrix, DisplacementVector);
            return strains;
        }
        private double[] CalculateIncrementStrainsVector(double[,] Bmatrix)
        {
            double[] deltaU = VectorOperations.VectorVectorSubtraction(DisplacementVector, Properties.DisplacementVectorPreviousStep);
            double[] strainsIncrement = VectorOperations.MatrixVectorProduct(Bmatrix, deltaU);
            return strainsIncrement;
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
        private Tuple<double[], double[]> LobattoPoints(int i, int j, int k)
        {
            double[] gaussPoints = new double[] { -1.0, 1.0 };
            double[] gaussWeights = new double[] { 1.0, 1.0 };

            double[] vectorWithPoints = new double[] { gaussPoints[i], gaussPoints[j], gaussPoints[k] };
            double[] vectorWithWeights = new double[] { gaussWeights[i], gaussWeights[j], gaussWeights[k] };
            return new Tuple<double[], double[]>(vectorWithPoints, vectorWithWeights);
        }
        private double[,] CreateEnhancedStrainsInterpolationMatrix(double[] ksi)
        {
            double[,] M = new double[6, 7];
            M[0, 0] = ksi[0];
            M[1, 1] = ksi[1];
            M[2, 2] = ksi[2];
            M[2, 5] = ksi[0] * ksi[2];
            M[2, 6] = ksi[1] * ksi[2];
            M[3, 3] = ksi[0];
            M[3, 4] = ksi[1];
            return M;
        }
        private double[,] CalculateEnhancedStrainMatrixGamma(double[,] transformationMat0, double[,] M,
            double detJ0, double detJ)
        {
            double scalar = detJ0 / detJ;
            double[,] gamma = MatrixOperations.ScalarMatrixProductNew(scalar, MatrixOperations.MatrixProduct(transformationMat0, M));
            return gamma;
        }
        private double[,] TransformationMatrixTransposed(double[,] jacobianMatrix)
        {
            double[,] TransposedT = new double[6, 6];
            TransposedT[0, 0] = jacobianMatrix[0, 0] * jacobianMatrix[0, 0];
            TransposedT[0, 1] = jacobianMatrix[0, 1] * jacobianMatrix[0, 1];
            TransposedT[0, 2] = jacobianMatrix[0, 2] * jacobianMatrix[0, 2];
            TransposedT[0, 3] = jacobianMatrix[0, 0] * jacobianMatrix[0, 1];
            TransposedT[0, 4] = jacobianMatrix[0, 0] * jacobianMatrix[0, 2];
            TransposedT[0, 5] = jacobianMatrix[0, 1] * jacobianMatrix[0, 2];

            TransposedT[1, 0] = jacobianMatrix[1, 0] * jacobianMatrix[1, 0];
            TransposedT[1, 1] = jacobianMatrix[1, 1] * jacobianMatrix[1, 1];
            TransposedT[1, 2] = jacobianMatrix[1, 2] * jacobianMatrix[1, 2];
            TransposedT[1, 3] = jacobianMatrix[1, 0] * jacobianMatrix[1, 1];
            TransposedT[1, 4] = jacobianMatrix[1, 0] * jacobianMatrix[1, 2];
            TransposedT[1, 5] = jacobianMatrix[1, 1] * jacobianMatrix[1, 2];

            TransposedT[2, 0] = jacobianMatrix[2, 0] * jacobianMatrix[2, 0];
            TransposedT[2, 1] = jacobianMatrix[2, 1] * jacobianMatrix[2, 1];
            TransposedT[2, 2] = jacobianMatrix[2, 2] * jacobianMatrix[2, 2];
            TransposedT[2, 3] = jacobianMatrix[2, 0] * jacobianMatrix[2, 1];
            TransposedT[2, 4] = jacobianMatrix[2, 0] * jacobianMatrix[2, 2];
            TransposedT[2, 5] = jacobianMatrix[2, 1] * jacobianMatrix[2, 2];

            TransposedT[3, 0] = 2 * jacobianMatrix[0, 0] * jacobianMatrix[1, 0];
            TransposedT[3, 1] = 2 * jacobianMatrix[0, 1] * jacobianMatrix[1, 1];
            TransposedT[3, 2] = 2 * jacobianMatrix[0, 2] * jacobianMatrix[1, 2];

            TransposedT[3, 3] = jacobianMatrix[0, 0] * jacobianMatrix[1, 1] +
                                jacobianMatrix[0, 1] * jacobianMatrix[1, 0];

            TransposedT[3, 4] = jacobianMatrix[0, 0] * jacobianMatrix[1, 2] +
                                jacobianMatrix[0, 2] * jacobianMatrix[1, 0];

            TransposedT[3, 5] = jacobianMatrix[0, 1] * jacobianMatrix[1, 2] +
                                jacobianMatrix[0, 2] * jacobianMatrix[1, 1];

            TransposedT[4, 0] = 2 * jacobianMatrix[0, 0] * jacobianMatrix[2, 0];
            TransposedT[4, 1] = 2 * jacobianMatrix[0, 1] * jacobianMatrix[2, 1];
            TransposedT[4, 2] = 2 * jacobianMatrix[0, 2] * jacobianMatrix[2, 2];

            TransposedT[4, 3] = jacobianMatrix[0, 0] * jacobianMatrix[2, 1] +
                                jacobianMatrix[0, 1] * jacobianMatrix[2, 0];

            TransposedT[4, 4] = jacobianMatrix[0, 0] * jacobianMatrix[2, 2] +
                                jacobianMatrix[0, 2] * jacobianMatrix[2, 0];

            TransposedT[4, 5] = jacobianMatrix[0, 1] * jacobianMatrix[2, 2] +
                                jacobianMatrix[0, 2] * jacobianMatrix[2, 1];

            TransposedT[5, 0] = 2 * jacobianMatrix[1, 0] * jacobianMatrix[2, 0];
            TransposedT[5, 1] = 2 * jacobianMatrix[1, 1] * jacobianMatrix[2, 1];
            TransposedT[5, 2] = 2 * jacobianMatrix[1, 2] * jacobianMatrix[2, 2];

            TransposedT[5, 3] = jacobianMatrix[1, 0] * jacobianMatrix[2, 1] +
                                jacobianMatrix[1, 1] * jacobianMatrix[2, 0];

            TransposedT[5, 4] = jacobianMatrix[1, 0] * jacobianMatrix[2, 2] +
                                jacobianMatrix[1, 2] * jacobianMatrix[2, 0];

            TransposedT[5, 5] = jacobianMatrix[1, 1] * jacobianMatrix[2, 2] +
                                jacobianMatrix[1, 2] * jacobianMatrix[2, 1];
            return TransposedT;
        }
        public double[,] CreateGlobalStiffnessMatrix()
        {
            double[,] ansK = new double[24, 24];
            double[,] LeMatrix = new double[7, 24];
            double[,] DeMatrix = new double[7, 7];
            double[,] E = CalculateStressStrainMatrix(Properties.YoungMod, Properties.PoissonRatio);
            double[] centerNaturalCoordinates = new double[3];
            Dictionary<string, double[]> localdN0 = CalculateShapeFunctionsLocalDerivatives(centerNaturalCoordinates);
            var J0 = CalculateJacobian(localdN0);
            var JInverse0 = CalculateInverseJacobian(J0);
            double detJ0 = JInverse0.Item2;
            double[,] invJacobian0 = JInverse0.Item1;
            var transformationMat0 = CalculateTransformationMatrix(invJacobian0);
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
                        var invJacobian = CalculateInverseJacobian(J);
                        var inverseJacobianMatrix = invJacobian.Item1;
                        double detJ = invJacobian.Item2;
                        var transformationMatrix = CalculateTransformationMatrix(inverseJacobianMatrix);

                        double[,] B = CalculateBMatrix(localdN, J, gP, transformationMatrix);
                        ansK = MatrixOperations.MatrixAddition(ansK, MatrixOperations.ScalarMatrixProductNew(detJ * gW[0] * gW[1] * gW[2],
                            MatrixOperations.MatrixProduct(MatrixOperations.Transpose(B), MatrixOperations.MatrixProduct(E, B))));
                        double[,] Gamma = CalculateEnhancedStrainMatrixGamma(transformationMat0, CreateEnhancedStrainsInterpolationMatrix(gP),
                        detJ0, detJ);
                        LeMatrix = MatrixOperations.MatrixAddition(LeMatrix, MatrixOperations.ScalarMatrixProductNew(detJ * gW[0] * gW[1] * gW[2],
                            MatrixOperations.MatrixProduct(
                            MatrixOperations.Transpose(Gamma), MatrixOperations.MatrixProduct(E, B))));
                        DeMatrix = MatrixOperations.MatrixAddition(DeMatrix, MatrixOperations.ScalarMatrixProductNew(detJ * gW[0] * gW[1] * gW[2],
                            MatrixOperations.MatrixProduct(
                            MatrixOperations.Transpose(Gamma), MatrixOperations.MatrixProduct(E, Gamma))));
                    }
                }
            }
            double[,] DeMatrixInv = MatrixOperations.CalculateInverseMatrixGaussJordanMethod(DeMatrix);
            //DeMatrixInv[0, 0] = 1.0 / DeMatrix[0, 0];
            double[,] tangentMatrix =
                MatrixOperations.MatrixSubtraction(
                ansK,
                MatrixOperations.MatrixProduct(
                MatrixOperations.MatrixProduct(MatrixOperations.Transpose(LeMatrix), DeMatrixInv), LeMatrix));
            Properties.EASLMatrix = LeMatrix;
            Properties.EASDMatrix = DeMatrixInv;
            return tangentMatrix;
        }
        public double[,] CreateMassMatrix()
        {
            double[,] M = new double[24, 24];
            double[] nodalX = UpdateNodalCoordinates(DisplacementVector);
            for (int i = 0; i < 2; i++)
            {
                for (int j = 0; j < 2; j++)
                {
                    for (int k = 0; k < 2; k++)
                    {
                        double[] lP = LobattoPoints(i, j, k).Item1;
                        double[] lW = LobattoPoints(i, j, k).Item2;
                        Dictionary<string, double[]> localdN = CalculateShapeFunctionsLocalDerivatives(lP);
                        double[,] J = CalculateUpdatedJacobianMatrix(nodalX, localdN);
                        //double[,] invJ = CalculateInverseJacobian(J).Item1;
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
            double[] VectorRe = new double[24];
            double[] PeVector = new double[7];
            double[,] matrixLe = Properties.EASLMatrix;
            double[,] matrixDeInv = Properties.EASDMatrix;
            double[,] E = CalculateStressStrainMatrix(Properties.YoungMod, Properties.PoissonRatio);
            double[] centerNaturalCoordinates = new double[3];
            Dictionary<string, double[]> localdN0 = CalculateShapeFunctionsLocalDerivatives(centerNaturalCoordinates);
            var J0 = CalculateJacobian(localdN0);
            var JInverse0 = CalculateInverseJacobian(J0);
            double detJ0 = JInverse0.Item2;
            double[,] invJacobian0 = JInverse0.Item1;
            var transformationMat0 = CalculateTransformationMatrix(invJacobian0);
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
                        var invJacobian = CalculateInverseJacobian(J);
                        var inverseJacobianMatrix = invJacobian.Item1;
                        double detJ = invJacobian.Item2;
                        var transformationMatrix = CalculateTransformationMatrix(inverseJacobianMatrix);
                        double[,] B = CalculateBMatrix(localdN, J, gP, transformationMatrix);
                        double[,] Gamma = CalculateEnhancedStrainMatrixGamma(transformationMat0, CreateEnhancedStrainsInterpolationMatrix(gP),
                        detJ0, detJ);
                        double[] alphaVector = Properties.EASVector;
                        double[] EnhStrainVector = VectorOperations.MatrixVectorProduct(Gamma, alphaVector);
                        double[] modifiedStrainVector = VectorOperations.VectorVectorAddition(CalculateStrainsVector(B), EnhStrainVector);
                        double[] stressVectorModified = CalculateStressVector(E, modifiedStrainVector);
                        VectorRe = VectorOperations.VectorVectorAddition(VectorRe, VectorOperations.VectorScalarProductNew(
                            VectorOperations.MatrixVectorProduct(MatrixOperations.Transpose(B), stressVectorModified), detJ * gW[0] * gW[1] * gW[2]));
                        PeVector = VectorOperations.VectorVectorAddition(PeVector,
                           VectorOperations.VectorScalarProductNew(
                           VectorOperations.MatrixVectorProduct(MatrixOperations.Transpose(Gamma), stressVectorModified),
                           detJ * gW[0] * gW[1] * gW[2]));
                    }
                }
            }
            double[] fInt = VectorOperations.VectorVectorSubtraction(VectorRe,
            VectorOperations.MatrixVectorProduct(MatrixOperations.MatrixProduct(MatrixOperations.Transpose(matrixLe), matrixDeInv), PeVector));
            Properties.EASFEnhancedVector = PeVector;
            return fInt;
        }
    }
}