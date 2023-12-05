using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace BumperCollitionSimulation
{
    class IsoparamShell18 : IElement
    {
        public Dictionary<int, INode> Nodes { get; }
        public Dictionary<int, INode> NodesMS { get; set; }

        public IElementProperties Properties { get; set; }
        public Dictionary<int, bool[]> ElementFreedomSignature { get; } = new Dictionary<int, bool[]>();
        public List<int> ElementFreedomList { get; set; }
        public double[] DisplacementVector { get; set; }
        public double[] AccelerationVector { get; set; }
        public double poisson { get; set; }
        public Dictionary<int, double> DegenerateShellNodalThicknesses { get; set; }
        public Dictionary<int, double[]> DegenerateShellNodalNormalUnitVectors { get; set; }
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
        public IsoparamShell18(IElementProperties properties, Dictionary<int, INode> nodes)
        {
            Properties = properties;
            this.Nodes = nodes;
            Dictionary<int, INode> newNodes = new Dictionary<int, INode>();
            for (int i = nodes.Keys.Min(); i <= nodes.Keys.Max() / 2; i++)
            {
                double x = (nodes.Single(n => n.Key == i).Value.XCoordinate +
                    nodes.Single(n => n.Key == i + nodes.Keys.Max() / 2).Value.XCoordinate) / 2.0;
                double y = (nodes.Single(n => n.Key == i).Value.YCoordinate +
                    nodes.Single(n => n.Key == i + nodes.Keys.Max() / 2).Value.YCoordinate) / 2.0;
                double z = (nodes.Single(n => n.Key == i).Value.ZCoordinate +
                    nodes.Single(n => n.Key == i + nodes.Keys.Max() / 2).Value.ZCoordinate) / 2.0;
                newNodes[i] = new Node(x, y, z);
            }
            this.NodesMS = newNodes;
            ElementFreedomSignature[1] = new bool[] { true, true, true, true, true, false };
            ElementFreedomSignature[2] = new bool[] { true, true, true, true, true, false };
            ElementFreedomSignature[3] = new bool[] { true, true, true, true, true, false };
            ElementFreedomSignature[4] = new bool[] { true, true, true, true, true, false };
            ElementFreedomSignature[5] = new bool[] { true, true, true, true, true, false };
            ElementFreedomSignature[6] = new bool[] { true, true, true, true, true, false };
            ElementFreedomSignature[7] = new bool[] { true, true, true, true, true, false };
            ElementFreedomSignature[8] = new bool[] { true, true, true, true, true, false };
            ElementFreedomSignature[9] = new bool[] { true, true, true, true, true, false };
            double[] pos1 = new double[] { nodes[1].XCoordinate, nodes[1].YCoordinate, nodes[1].ZCoordinate };
            double[] pos2 = new double[] { nodes[2].XCoordinate, nodes[2].YCoordinate, nodes[2].ZCoordinate };
            double[] pos3 = new double[] { nodes[3].XCoordinate, nodes[3].YCoordinate, nodes[3].ZCoordinate };
            double[] pos4 = new double[] { nodes[4].XCoordinate, nodes[4].YCoordinate, nodes[4].ZCoordinate };
            double[] pos5 = new double[] { nodes[5].XCoordinate, nodes[5].YCoordinate, nodes[5].ZCoordinate };
            double[] pos6 = new double[] { nodes[6].XCoordinate, nodes[6].YCoordinate, nodes[6].ZCoordinate };
            double[] pos7 = new double[] { nodes[7].XCoordinate, nodes[7].YCoordinate, nodes[7].ZCoordinate };
            double[] pos8 = new double[] { nodes[8].XCoordinate, nodes[8].YCoordinate, nodes[8].ZCoordinate };
            double[] pos9 = new double[] { nodes[9].XCoordinate, nodes[9].YCoordinate, nodes[9].ZCoordinate };
            double[] pos10 = new double[] { nodes[10].XCoordinate, nodes[10].YCoordinate, nodes[10].ZCoordinate };
            double[] pos11 = new double[] { nodes[11].XCoordinate, nodes[11].YCoordinate, nodes[11].ZCoordinate };
            double[] pos12 = new double[] { nodes[12].XCoordinate, nodes[12].YCoordinate, nodes[12].ZCoordinate };
            double[] pos13 = new double[] { nodes[13].XCoordinate, nodes[13].YCoordinate, nodes[13].ZCoordinate };
            double[] pos14 = new double[] { nodes[14].XCoordinate, nodes[14].YCoordinate, nodes[14].ZCoordinate };
            double[] pos15 = new double[] { nodes[15].XCoordinate, nodes[15].YCoordinate, nodes[15].ZCoordinate };
            double[] pos16 = new double[] { nodes[16].XCoordinate, nodes[16].YCoordinate, nodes[16].ZCoordinate };
            double[] pos17 = new double[] { nodes[17].XCoordinate, nodes[17].YCoordinate, nodes[17].ZCoordinate };
            double[] pos18 = new double[] { nodes[18].XCoordinate, nodes[18].YCoordinate, nodes[18].ZCoordinate };
            Dictionary<int, double> tNodal = new Dictionary<int, double>();
            tNodal.Add(1, VectorOperations.VectorNorm2(VectorOperations.VectorVectorSubtraction(pos10, pos1)));
            tNodal.Add(2, VectorOperations.VectorNorm2(VectorOperations.VectorVectorSubtraction(pos11, pos2)));
            tNodal.Add(3, VectorOperations.VectorNorm2(VectorOperations.VectorVectorSubtraction(pos12, pos3)));
            tNodal.Add(4, VectorOperations.VectorNorm2(VectorOperations.VectorVectorSubtraction(pos13, pos4)));
            tNodal.Add(5, VectorOperations.VectorNorm2(VectorOperations.VectorVectorSubtraction(pos14, pos5)));
            tNodal.Add(6, VectorOperations.VectorNorm2(VectorOperations.VectorVectorSubtraction(pos15, pos6)));
            tNodal.Add(7, VectorOperations.VectorNorm2(VectorOperations.VectorVectorSubtraction(pos16, pos7)));
            tNodal.Add(8, VectorOperations.VectorNorm2(VectorOperations.VectorVectorSubtraction(pos17, pos8)));
            tNodal.Add(9, VectorOperations.VectorNorm2(VectorOperations.VectorVectorSubtraction(pos18, pos9)));
            this.Properties.DegenerateShellNodalThicknesses = tNodal;
            Dictionary<int, double[]> unitNormalVectorsNodal = new Dictionary<int, double[]>();
            unitNormalVectorsNodal.Add(1, VectorOperations.VectorScalarProduct(
                VectorOperations.VectorVectorSubtraction(pos10, pos1), 1.0 / tNodal[1]));
            unitNormalVectorsNodal.Add(2, VectorOperations.VectorScalarProduct(
                VectorOperations.VectorVectorSubtraction(pos11, pos2), 1.0 / tNodal[2]));
            unitNormalVectorsNodal.Add(3, VectorOperations.VectorScalarProduct(
                VectorOperations.VectorVectorSubtraction(pos12, pos3), 1.0 / tNodal[3]));
            unitNormalVectorsNodal.Add(4, VectorOperations.VectorScalarProduct(
                VectorOperations.VectorVectorSubtraction(pos13, pos4), 1.0 / tNodal[4]));
            unitNormalVectorsNodal.Add(5, VectorOperations.VectorScalarProduct(
                VectorOperations.VectorVectorSubtraction(pos14, pos5), 1.0 / tNodal[5]));
            unitNormalVectorsNodal.Add(6, VectorOperations.VectorScalarProduct(
                VectorOperations.VectorVectorSubtraction(pos15, pos6), 1.0 / tNodal[6]));
            unitNormalVectorsNodal.Add(7, VectorOperations.VectorScalarProduct(
                VectorOperations.VectorVectorSubtraction(pos16, pos7), 1.0 / tNodal[7]));
            unitNormalVectorsNodal.Add(8, VectorOperations.VectorScalarProduct(
                VectorOperations.VectorVectorSubtraction(pos17, pos8), 1.0 / tNodal[8]));
            unitNormalVectorsNodal.Add(9, VectorOperations.VectorScalarProduct(
                VectorOperations.VectorVectorSubtraction(pos18, pos9), 1.0 / tNodal[9]));
            this.Properties.DegenerateShellNodalNormalUnitVectors = unitNormalVectorsNodal;
            DisplacementVector = new double[45];
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
        public List<double[]> GetStressVector()
        {
            throw new Exception("Needs to be added.");
        }
        public List<double[]> GetStrainVector()
        {
            throw new Exception("Needs to be added.");
        }
        public List<double[]> GetGaussPointsInPhysicalSpace()
        {
            throw new Exception("Needs to be added.");
        }
        public List<double[]> GetStressFromElementsNodes()
        {
            throw new Exception("Needs to be added.");
        }
        public List<double[]> GetStrainFromElementsNodes()
        {
            throw new Exception("Needs to be added.");
        }
        public List<double[]> GetStressFromElements(List<double[]> parametricCoordinatesVector)
        {
            throw new Exception("Needs to be added.");
        }
        public List<double[]> GetphysicalCoordinatesFromElements(List<double[]> parametricCoordinatesVector)
        {
            throw new Exception("Needs to be added.");
        }
        private Dictionary<int, Dictionary<int, double[]>> CalculateNodalUnitVectors(double[] msNodalCoordinates,
            Dictionary<int, double> nodalThicknesses)
        {
            Dictionary<int, Dictionary<int, double[]>> nodalUnitVectors =
                new Dictionary<int, Dictionary<int, double[]>>();
            double deltaKsi = 0.01;
            double deltaIhta = 0.01;
            List<double[]> nodalParamCoordinates = new List<double[]>();
            nodalParamCoordinates.Add(new double[]
            {-1.0,-1.0 });
            nodalParamCoordinates.Add(new double[]
            {0.0,-1.0 });
            nodalParamCoordinates.Add(new double[]
            {1.0,-1.0 });
            nodalParamCoordinates.Add(new double[]
            {1.0, 0.0 });
            nodalParamCoordinates.Add(new double[]
            {1.0, 1.0 });
            nodalParamCoordinates.Add(new double[]
            {0.0, 1.0 });
            nodalParamCoordinates.Add(new double[]
            {-1.0, 1.0 });
            nodalParamCoordinates.Add(new double[]
            {-1.0, 0.0 });
            nodalParamCoordinates.Add(new double[]
            {0.0, 0.0 });
            double deltaX1 = 0;
            double deltaY1 = 0;
            double deltaZ1 = 0;
            double deltaX2 = 0;
            double deltaY2 = 0;
            double deltaZ2 = 0;
            for (int i = 0; i< nodalParamCoordinates.Count; i++)
            {
                double[] naturalCoordVec = nodalParamCoordinates[i];
                var derivatives = CalculateShapeFunctionsLocalDerivatives(naturalCoordVec);
                for(int j = 0; j < 9; j++)
                {
                    deltaX1 += derivatives.Single(d => d.Key == "ksi").Value[j] *
                        msNodalCoordinates[j * 5] * deltaKsi;
                    deltaY1 += derivatives.Single(d => d.Key == "ksi").Value[j] *
                        msNodalCoordinates[j * 5 + 1] * deltaKsi;
                    deltaZ1 += derivatives.Single(d => d.Key == "ksi").Value[j] *
                        msNodalCoordinates[j * 5 + 2] * deltaKsi;
                    deltaX2 += derivatives.Single(d => d.Key == "ihta").Value[j] *
                        msNodalCoordinates[j * 5] * deltaIhta;
                    deltaY2 += derivatives.Single(d => d.Key == "ihta").Value[j] *
                        msNodalCoordinates[j * 5 + 1] * deltaIhta;
                    deltaZ2 += derivatives.Single(d => d.Key == "ihta").Value[j] *
                        msNodalCoordinates[j * 5 + 2] * deltaIhta;
                }

                double[] e1 = new double[] { deltaX1, deltaY1, deltaZ1 };
                double[] e2 = new double[] { deltaX2, deltaY2, deltaZ2 };
                double[] e3 = VectorOperations.VectorCrossProduct(e1, e2);
                double[] v1 = VectorOperations.VectorScalarProduct(e1,
                            1.0 / VectorOperations.VectorNorm2(e1));
                double[] v3 = VectorOperations.VectorScalarProduct(e3,
                            1.0 / VectorOperations.VectorNorm2(e3));
                double[] v = VectorOperations.VectorCrossProduct(v3, v1);
                double[] v2 = VectorOperations.VectorScalarProduct(v,
                            1.0 / VectorOperations.VectorNorm2(v));
                Dictionary<int, double[]> vectors = new Dictionary<int, double[]>();
                vectors.Add(1, v1);
                vectors.Add(2, v2);
                vectors.Add(3, v3);
                nodalUnitVectors.Add(i + 1, vectors);
            }
            return nodalUnitVectors;
        }
        private double[] UpdateNodalCoordinates(double[] displacementVector)
        {
            double[] updatedMSNodalCoor = new double[45];
            for (int i = 1; i <= 9; i++)
            {
                updatedMSNodalCoor[5 * i - 5] = NodesMS[i].XCoordinate + displacementVector[5 * i - 5];
                updatedMSNodalCoor[5 * i - 4] = NodesMS[i].YCoordinate + displacementVector[5 * i - 4];
                updatedMSNodalCoor[5 * i - 3] = NodesMS[i].ZCoordinate + displacementVector[5 * i - 3];
                updatedMSNodalCoor[5 * i - 2] = NodesMS[i].RX + displacementVector[5 * i - 2];
                updatedMSNodalCoor[5 * i - 1] = NodesMS[i].RY + displacementVector[5 * i - 1];
            }
            return updatedMSNodalCoor;
        }
        private void UpdateNodalGeometry(double[] updatedMSNodalCoor)
        {
            Dictionary<int, double> nodalT = new Dictionary<int, double>();
            Dictionary<int, double[]> nodalUnitVectors = new Dictionary<int, double[]>();
            List<double[]> parametriCoordinatesMs = new List<double[]>();
            parametriCoordinatesMs.Add(new double[] { -1.0, -1.0 });
            parametriCoordinatesMs.Add(new double[] { 0.0, -1.0 });
            parametriCoordinatesMs.Add(new double[] { 1.0, -1.0 });
            parametriCoordinatesMs.Add(new double[] { 1.0, 0.0 });
            parametriCoordinatesMs.Add(new double[] { 1.0, 1.0 });
            parametriCoordinatesMs.Add(new double[] { 0.0, 1.0 });
            parametriCoordinatesMs.Add(new double[] { -1.0, 1.0 });
            parametriCoordinatesMs.Add(new double[] { -1.0, 0.0 });
            parametriCoordinatesMs.Add(new double[] { 0.0, 0.0 });
            for(int i = 0; i< parametriCoordinatesMs.Count; i++)
            {
                double thickness = new double();
                double[] unitNormal = new double[3];
                double xUpper = new double();
                double yUpper = new double();
                double zUpper = new double();
                double xLower = new double();
                double yLower = new double();
                double zLower = new double();
                var N = CalculateShapeFunctions(parametriCoordinatesMs[i][0], parametriCoordinatesMs[i][1]);
                for(int j = 0; j < 9; j++)
                {
                    xUpper += N[j + 1] * (updatedMSNodalCoor[j * 5] +
                        Properties.DegenerateShellNodalNormalUnitVectors[j + 1][0]
                        * Properties.DegenerateShellNodalThicknesses[j + 1] / 2.0);
                    xLower += N[j + 1] * (updatedMSNodalCoor[j * 5] -
                        Properties.DegenerateShellNodalNormalUnitVectors[j + 1][0]
                        * Properties.DegenerateShellNodalThicknesses[j + 1] / 2.0);
                    yUpper += N[j + 1] * (updatedMSNodalCoor[j * 5 + 1] +
                        Properties.DegenerateShellNodalNormalUnitVectors[j + 1][1]
                        * Properties.DegenerateShellNodalThicknesses[j + 1] / 2.0);
                    yLower += N[j + 1] * (updatedMSNodalCoor[j * 5 + 1] -
                        Properties.DegenerateShellNodalNormalUnitVectors[j + 1][1]
                        * Properties.DegenerateShellNodalThicknesses[j + 1] / 2.0);
                    zUpper += N[j + 1] * (updatedMSNodalCoor[j * 5 + 2] +
                        Properties.DegenerateShellNodalNormalUnitVectors[j + 1][2]
                        * Properties.DegenerateShellNodalThicknesses[j + 1] / 2.0);
                    zLower += N[j + 1] * (updatedMSNodalCoor[j * 5 + 2] -
                        Properties.DegenerateShellNodalNormalUnitVectors[j + 1][2]
                        * Properties.DegenerateShellNodalThicknesses[j + 1] / 2.0);
                }
                thickness = VectorOperations.VectorNorm2(
                    VectorOperations.VectorVectorSubtraction(new double[] { xUpper, yUpper, zUpper },
                    new double[] { xLower, yLower, zLower }));

                nodalT.Add(i + 1, thickness);
                nodalUnitVectors.Add(i + 1, VectorOperations.VectorScalarProduct(
                VectorOperations.VectorVectorSubtraction(new double[] { xUpper, yUpper, zUpper },
                    new double[] { xLower, yLower, zLower }), 1.0 / thickness));
            }
            for (int i = 0; i < parametriCoordinatesMs.Count; i++)
            {
                Properties.DegenerateShellNodalThicknesses[i + 1] = nodalT[i + 1];
                Properties.DegenerateShellNodalNormalUnitVectors[i + 1] = nodalUnitVectors[i + 1];
            }
        }
        public Dictionary<int, INode> NodesAtFinalState()
        {
            Dictionary<int, INode> finalNodes = new Dictionary<int, INode>();
            //finalNodes[1] = new Node(Nodes[1].XCoordinate + DisplacementVector[0], Nodes[1].YCoordinate + DisplacementVector[1]);
            //finalNodes[2] = new Node(Nodes[2].XCoordinate + DisplacementVector[2], Nodes[2].YCoordinate + DisplacementVector[3]);
            //finalNodes[3] = new Node(Nodes[3].XCoordinate + DisplacementVector[4], Nodes[3].YCoordinate + DisplacementVector[5]);
            //finalNodes[4] = new Node(Nodes[4].XCoordinate + DisplacementVector[6], Nodes[4].YCoordinate + DisplacementVector[7]);
            //finalNodes[5] = new Node(Nodes[5].XCoordinate + DisplacementVector[8], Nodes[5].YCoordinate + DisplacementVector[9]);
            //finalNodes[6] = new Node(Nodes[6].XCoordinate + DisplacementVector[10], Nodes[6].YCoordinate + DisplacementVector[11]);
            //finalNodes[7] = new Node(Nodes[7].XCoordinate + DisplacementVector[12], Nodes[7].YCoordinate + DisplacementVector[13]);
            //finalNodes[8] = new Node(Nodes[8].XCoordinate + DisplacementVector[14], Nodes[8].YCoordinate + DisplacementVector[15]);
            return finalNodes;
        }

        private Dictionary<int, double> CalculateShapeFunctions(double ksi, double ihta)
        {
            Dictionary<int, double> shapeFunctions = new Dictionary<int, double>();
            double N1 = 1.0 / 4.0 * ksi * ihta * (1 - ksi) * (1 - ihta); shapeFunctions.Add(1, N1);
            double N2 = -1.0 / 2.0 * ihta * (1 + ksi) * (1 - ksi) * (1 - ihta); shapeFunctions.Add(2, N2);
            double N3 = -1.0 / 4.0 * ksi * ihta * (1 + ksi) * (1 - ihta); shapeFunctions.Add(3, N3);
            double N4 = 1.0 / 2.0 * ksi * (1 + ksi) * (1 + ihta) * (1 - ihta); shapeFunctions.Add(4, N4);
            double N5 = 1.0 / 4.0 * ksi * ihta * (1 + ksi) * (1 + ihta); shapeFunctions.Add(5, N5);
            double N6 = 1.0 / 2.0 * ihta * (1 - ksi) * (1 + ksi) * (1 + ihta); shapeFunctions.Add(6, N6);
            double N7 = -1.0 / 4.0 * ksi * ihta * (1 - ksi) * (1 + ihta); shapeFunctions.Add(7, N7);
            double N8 = -1.0 / 2.0 * ksi * (1 - ksi) * (1 + ihta) * (1 - ihta); shapeFunctions.Add(8, N8);
            double N9 = (1 - Math.Pow(ksi, 2)) * (1 + ihta) * (1 - ihta); shapeFunctions.Add(9, N9);
            return shapeFunctions;
        }

        private Tuple<double[,], Dictionary<int, double>> CalculateShapeFunctionMatrix(double ksi, double ihta)
        {
            Dictionary<int, double> shapeFunctions = CalculateShapeFunctions(ksi, ihta);
            double[,] N = new double[,]
            {
                {shapeFunctions[1], 0, 0,
                shapeFunctions[2], 0, 0,
                shapeFunctions[3], 0, 0,
                shapeFunctions[4], 0, 0,
                shapeFunctions[5], 0, 0,
                shapeFunctions[6], 0, 0,
                shapeFunctions[7], 0, 0,
                shapeFunctions[8], 0, 0,
                shapeFunctions[9], 0, 0 },
                {0, shapeFunctions[1], 0,
                0, shapeFunctions[2], 0,
                0, shapeFunctions[3], 0,
                0, shapeFunctions[4], 0,
                0, shapeFunctions[5], 0,
                0, shapeFunctions[6], 0,
                0, shapeFunctions[7], 0,
                0, shapeFunctions[8], 0,
                0, shapeFunctions[9], 0, },
                {0, 0, shapeFunctions[1],
                0, 0, shapeFunctions[2],
                0, 0, shapeFunctions[3],
                0, 0, shapeFunctions[4],
                0, 0, shapeFunctions[5],
                0, 0, shapeFunctions[6],
                0, 0, shapeFunctions[7],
                0, 0, shapeFunctions[8],
                0, 0, shapeFunctions[9]}
            };
            return new Tuple<double[,], Dictionary<int, double>> (N, shapeFunctions);
        }

        private Dictionary<string, double[]> CalculateShapeFunctionsLocalDerivatives(double[] naturalCoordinates)
        {
            double ksi = naturalCoordinates[0];
            double ihta = naturalCoordinates[1];

            double[] dN_ksi = new double[]
            {
                (1.0/4.0 * ihta * (1-ihta)*(1 - 2*ksi)),
                ((1 - ihta)*ksi*ihta),
                (-1.0/4.0 * ihta * (1 - ihta) * (1 + 2*ksi )),
                (1.0/2.0*(1 + 2*ksi) * (1 - Math.Pow(ihta,2))),
                (1.0/4.0*ihta*(1 + ihta)*(2*ksi + 1)),
                (-(ihta+1)*ksi*ihta),
                (-1.0/4.0*ihta*(1+ihta)*(1 - 2*ksi)),
                (-1.0/2.0*(1 - 2*ksi) * (1 - Math.Pow(ihta,2))),
                (-2*ksi * (1 - Math.Pow(ihta,2))),
            };
            double[] dN_ihta = new double[]
            {
                (1.0/4.0 * ksi * (1 - ksi)*(1 - 2*ihta)),
                (-1.0/2.0*(1 - Math.Pow(ksi,2))*(1 - 2*ihta)),
                (-1.0/4.0 * ksi * (1 + ksi)*(1 - 2*ihta)),
                (-(ksi + 1) * ihta * ksi),
                (1.0/4.0 * ksi * (1 + ksi) * (1 + 2*ihta)),
                (1.0/2.0*(1 - Math.Pow(ksi,2)) * (1 + 2*ihta)),
                (-1.0/4.0 * ksi * (1 - ksi)*(1 + 2*ihta)),
                ((1 - ksi)*ksi*ihta),
                (-2*ihta * (1 - Math.Pow(ksi,2))),
            };
            Dictionary<string, double[]> dN = new Dictionary<string, double[]>();
            dN.Add("ksi", dN_ksi);
            dN.Add("ihta", dN_ihta);
            return dN;
        }
        private double[,] CalculateJacobian(double zita, Dictionary<string, double[]> dN,
            Dictionary<int, double> N, double[] xUpdated,
            Dictionary<int, Dictionary<int, double[]>> nodesDirectionalCosines)
        {
            double[,] jacobianMatrix = new double[3, 3];
            int k = 0;
            for (int i = 0; i < 9; i++)
            {
                jacobianMatrix[0, 0] = jacobianMatrix[0, 0] + (xUpdated[k] + zita *
                    Properties.DegenerateShellNodalThicknesses.Single(t=>t.Key == i + 1).Value / 2.0 *
                    nodesDirectionalCosines.Single(c=>c.Key == i + 1).Value[3][0])
                    * dN["ksi"][i];
                k += 5;
            }
            k = 1;
            for (int i = 0; i < 9; i++)
            {
                jacobianMatrix[0, 1] = jacobianMatrix[0, 1] + (xUpdated[k] + zita *
                    Properties.DegenerateShellNodalThicknesses.Single(t => t.Key == i + 1).Value / 2.0 *
                    nodesDirectionalCosines.Single(c => c.Key == i + 1).Value[3][1])
                    * dN["ksi"][i];
                k += 5;
            }
            k = 2;
            for (int i = 0; i < 9; i++)
            {
                jacobianMatrix[0, 2] = jacobianMatrix[0, 2] + (xUpdated[k] + zita *
                    Properties.DegenerateShellNodalThicknesses.Single(t => t.Key == i + 1).Value / 2.0 *
                    nodesDirectionalCosines.Single(c => c.Key == i + 1).Value[3][2])
                    * dN["ksi"][i];
                k += 5;
            }

            k = 0;
            for (int i = 0; i < 9; i++)
            {
                jacobianMatrix[1, 0] = jacobianMatrix[1, 0] + (xUpdated[k] + zita *
                    Properties.DegenerateShellNodalThicknesses.Single(t => t.Key == i + 1).Value / 2.0 *
                    nodesDirectionalCosines.Single(c => c.Key == i + 1).Value[3][0]) * dN["ihta"][i];
                k += 5;
            }
            k = 1;
            for (int i = 0; i < 9; i++)
            {
                jacobianMatrix[1, 1] = jacobianMatrix[1, 1] + (xUpdated[k] + zita *
                    Properties.DegenerateShellNodalThicknesses.Single(t => t.Key == i + 1).Value / 2.0 *
                    nodesDirectionalCosines.Single(c => c.Key == i + 1).Value[3][1]) * dN["ihta"][i];
                k += 5;
            }
            k = 2;
            for (int i = 0; i < 9; i++)
            {
                jacobianMatrix[1, 2] = jacobianMatrix[1, 2] + (xUpdated[k] + zita *
                    Properties.DegenerateShellNodalThicknesses.Single(t => t.Key == i + 1).Value / 2.0 *
                    nodesDirectionalCosines.Single(c => c.Key == i + 1).Value[3][2]) * dN["ihta"][i];
                k += 5;
            }

            k = 0;
            for (int i = 0; i < 9; i++)
            {
                jacobianMatrix[2, 0] = jacobianMatrix[2, 0] +
                    Properties.DegenerateShellNodalThicknesses.Single(t => t.Key == i + 1).Value / 2.0 *
                    nodesDirectionalCosines.Single(c => c.Key == i + 1).Value[3][0] * N[i + 1];
                k += 5;
            }
            k = 1;
            for (int i = 0; i < 9; i++)
            {
                jacobianMatrix[2, 1] = jacobianMatrix[2, 1] +
                    Properties.DegenerateShellNodalThicknesses.Single(t => t.Key == i + 1).Value / 2.0 *
                    nodesDirectionalCosines.Single(c => c.Key == i + 1).Value[3][1] * N[i + 1];
                k += 5;
            }
            k = 2;
            for (int i = 0; i < 9; i++)
            {
                jacobianMatrix[2, 2] = jacobianMatrix[2, 2] +
                    Properties.DegenerateShellNodalThicknesses.Single(t => t.Key == i + 1).Value / 2.0 *
                    nodesDirectionalCosines.Single(c => c.Key == i + 1).Value[3][2] * N[i + 1];
                k += 5;
            }
            return jacobianMatrix;
        }

        private Tuple<double[,], double, double[,]> CalculateInverseJacobian(double[,] jacobianMatrix)
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
            double[,] B1 = new double[6, 9];
            B1[0, 0] = jacobianInverseMatrix[0, 0];
            B1[0, 1] = jacobianInverseMatrix[0, 1];
            B1[0, 2] = jacobianInverseMatrix[0, 2];
            B1[1, 3] = jacobianInverseMatrix[1, 0];
            B1[1, 4] = jacobianInverseMatrix[1, 1];
            B1[1, 5] = jacobianInverseMatrix[1, 2];
            B1[2, 6] = jacobianInverseMatrix[2, 0];
            B1[2, 7] = jacobianInverseMatrix[2, 1];
            B1[2, 8] = jacobianInverseMatrix[2, 2];

            B1[3, 0] = jacobianInverseMatrix[1, 0];
            B1[3, 1] = jacobianInverseMatrix[1, 1];
            B1[3, 2] = jacobianInverseMatrix[1, 2];
            B1[3, 3] = jacobianInverseMatrix[0, 0];
            B1[3, 4] = jacobianInverseMatrix[0, 1];
            B1[3, 5] = jacobianInverseMatrix[0, 2];

            B1[4, 3] = jacobianInverseMatrix[2, 0];
            B1[4, 4] = jacobianInverseMatrix[2, 1];
            B1[4, 5] = jacobianInverseMatrix[2, 2];
            B1[4, 6] = jacobianInverseMatrix[1, 0];
            B1[4, 7] = jacobianInverseMatrix[1, 1];
            B1[4, 8] = jacobianInverseMatrix[1, 2];

            B1[5, 0] = jacobianInverseMatrix[2, 0];
            B1[5, 1] = jacobianInverseMatrix[2, 1];
            B1[5, 2] = jacobianInverseMatrix[2, 2];
            B1[5, 6] = jacobianInverseMatrix[0, 0];
            B1[5, 7] = jacobianInverseMatrix[0, 1];
            B1[5, 8] = jacobianInverseMatrix[0, 2];
            return new Tuple<double[,], double, double[,]>(jacobianInverseMatrix, detj, B1);
        }
        private double[] CalculateStrainsVector(double[,] Bmatrix)
        {
            double[] strains = VectorOperations.MatrixVectorProduct(Bmatrix, DisplacementVector);
            return strains;
        }

        private double[,] CalculateBMatrix(double zita, double[,] B1matrix, Dictionary<string, double[]> localdN,
            Dictionary<int, double> N, Dictionary<int, Dictionary<int, double[]>> nodesDirectionalCosines)
        {
            double[,] Bmatrix = new double[6, 45];
            double[,] B2matrix = new double[9, 45];

            for (int j = 0; j < 9; j++)
            {
                B2matrix[0, j * 5] = localdN["ksi"][j];
                B2matrix[0, j * 5 + 3] = -localdN["ksi"][j] * nodesDirectionalCosines[j + 1][2][0] *
                    zita * Properties.DegenerateShellNodalThicknesses[j + 1] / 2.0;
                B2matrix[0, j * 5 + 4] = localdN["ksi"][j] * nodesDirectionalCosines[j + 1][1][0] *
                    zita * Properties.DegenerateShellNodalThicknesses[j + 1] / 2.0;

                B2matrix[1, j * 5] = localdN["ihta"][j];
                B2matrix[1, j * 5 + 3] = -localdN["ihta"][j] * nodesDirectionalCosines[j + 1][2][0] *
                    zita * Properties.DegenerateShellNodalThicknesses[j + 1] / 2.0;
                B2matrix[1, j * 5 + 4] = localdN["ihta"][j] * nodesDirectionalCosines[j + 1][1][0] *
                    zita * Properties.DegenerateShellNodalThicknesses[j + 1] / 2.0;

                B2matrix[2, j * 5 + 3] = -N[j + 1] * nodesDirectionalCosines[j + 1][2][0] *
                    Properties.DegenerateShellNodalThicknesses[j + 1] / 2.0;
                B2matrix[2, j * 5 + 4] = N[j + 1] * nodesDirectionalCosines[j + 1][1][0] *
                    Properties.DegenerateShellNodalThicknesses[j + 1] / 2.0;
                //
                //
                B2matrix[3, j * 5 + 1] = localdN["ksi"][j];
                B2matrix[3, j * 5 + 3] = -localdN["ksi"][j] * nodesDirectionalCosines[j + 1][2][1] *
                    zita * Properties.DegenerateShellNodalThicknesses[j + 1] / 2.0;
                B2matrix[3, j * 5 + 4] = localdN["ksi"][j] * nodesDirectionalCosines[j + 1][1][1] *
                    zita * Properties.DegenerateShellNodalThicknesses[j + 1] / 2.0;

                B2matrix[4, j * 5 + 1] = localdN["ihta"][j];
                B2matrix[4, j * 5 + 3] = -localdN["ihta"][j] * nodesDirectionalCosines[j + 1][2][1] *
                    zita * Properties.DegenerateShellNodalThicknesses[j + 1] / 2.0;
                B2matrix[4, j * 5 + 4] = localdN["ihta"][j] * nodesDirectionalCosines[j + 1][1][1] *
                    zita * Properties.DegenerateShellNodalThicknesses[j + 1] / 2.0;

                B2matrix[5, j * 5 + 3] = -N[j + 1] * nodesDirectionalCosines[j + 1][2][1] *
                    Properties.DegenerateShellNodalThicknesses[j + 1] / 2.0;
                B2matrix[5, j * 5 + 4] = N[j + 1] * nodesDirectionalCosines[j + 1][1][1] *
                    Properties.DegenerateShellNodalThicknesses[j + 1] / 2.0;
                //
                //
                B2matrix[6, j * 5 + 2] = localdN["ksi"][j];
                B2matrix[6, j * 5 + 3] = -localdN["ksi"][j] * nodesDirectionalCosines[j + 1][2][2] *
                    zita * Properties.DegenerateShellNodalThicknesses[j + 1] / 2.0;
                B2matrix[6, j * 5 + 4] = localdN["ksi"][j] * nodesDirectionalCosines[j + 1][1][2] *
                    zita * Properties.DegenerateShellNodalThicknesses[j + 1] / 2.0;

                B2matrix[7, j * 5 + 2] = localdN["ihta"][j];
                B2matrix[7, j * 5 + 3] = -localdN["ihta"][j] * nodesDirectionalCosines[j + 1][2][2] *
                    zita * Properties.DegenerateShellNodalThicknesses[j + 1] / 2.0;
                B2matrix[7, j * 5 + 4] = localdN["ihta"][j] * nodesDirectionalCosines[j + 1][1][2] *
                    zita * Properties.DegenerateShellNodalThicknesses[j + 1] / 2.0;

                B2matrix[8, j * 5 + 3] = -N[j + 1] * nodesDirectionalCosines[j + 1][2][2] *
                    Properties.DegenerateShellNodalThicknesses[j + 1] / 2.0;
                B2matrix[8, j * 5 + 4] = N[j + 1] * nodesDirectionalCosines[j + 1][1][2] *
                    Properties.DegenerateShellNodalThicknesses[j + 1] / 2.0;
            }
            Bmatrix = MatrixOperations.MatrixProduct(B1matrix, B2matrix);
            return Bmatrix;
        }

        private double[,] CalculateStressStrainMatrix(double E, double v)
        {
            double[,] Ematrix = new double[6, 6];
            //double cs = 1.0; //for ortho shear stress distribution along thickness
            double cs = 5.0 / 6.0;// parabolic distribution
            double Ehat = E / (1.0 - Math.Pow(v, 2));
            Ematrix[0, 0] = Ehat;
            Ematrix[0, 1] = Ehat * v;
            Ematrix[1, 0] = Ehat * v;
            Ematrix[1, 1] = Ehat;
            Ematrix[3, 3] = Ehat * (1.0 / 2.0) * (1.0 - v);
            Ematrix[4, 4] = Ehat * (1.0 / 2.0) * (1.0 - v) * cs;
            Ematrix[5, 5] = Ehat * (1.0 / 2.0) * (1.0 - v) * cs;
            return Ematrix;
        }

        private double[] CalculateStressVector(double[,] E, double[] strain)
        {
            double[] stressVector = VectorOperations.MatrixVectorProduct(E, strain);
            return stressVector;
        }

        private Tuple<double[], double[]> GaussPoints(int i, int j, int k)
        {
            double[] gaussPoints = new double[] { -0.77459, 0.0, 0.77459 };
            double[] gaussWeights = new double[] { 0.55555, 0.88888, 0.55555 };

            double[] vectorWithPoints = new double[] { gaussPoints[i], gaussPoints[j], gaussPoints[k] };
            double[] vectorWithWeights = new double[] { gaussWeights[i], gaussWeights[j], gaussWeights[k] };
            return new Tuple<double[], double[]>(vectorWithPoints, vectorWithWeights);
        }
        private Tuple<double[], double[]> LobattoPoints(int i, int j)
        {
            double[] lobattoPoints = new double[] { -1.0, 0.0, 1.0 };
            double[] lobattoWeights = new double[] { 1.0 / 3.0, 4.0 / 3.0, 1.0 /3.0 };

            double[] vectorWithPoints = new double[] { lobattoPoints[i], lobattoPoints[j]};
            double[] vectorWithWeights = new double[] { lobattoWeights[i], lobattoWeights[j]};
            return new Tuple<double[], double[]>(vectorWithPoints, vectorWithWeights);
        }
        private Dictionary<int, double[]> unitVectors(Dictionary<int, double> N,
            Dictionary<int, Dictionary<int, double[]>> nodesDirectionalCosines)
        {
            Dictionary<int, double[]> result = new Dictionary<int, double[]>();
            double v1X = new double();
            double v2X = new double();
            double v3X = new double();
            double v1Y = new double();
            double v2Y = new double();
            double v3Y = new double();
            double v1Z = new double();
            double v2Z = new double();
            double v3Z = new double();
            for (int i = 1; i <= 9; i++)
            {
                v1X += N[i] * nodesDirectionalCosines[i][1][0];
                v1Y += N[i] * nodesDirectionalCosines[i][1][1];
                v1Z += N[i] * nodesDirectionalCosines[i][1][2];
                v2X += N[i] * nodesDirectionalCosines[i][2][0];
                v2Y += N[i] * nodesDirectionalCosines[i][2][1];
                v2Z += N[i] * nodesDirectionalCosines[i][2][2];
                v3X += N[i] * nodesDirectionalCosines[i][3][0];
                v3Y += N[i] * nodesDirectionalCosines[i][3][1];
                v3Z += N[i] * nodesDirectionalCosines[i][3][2];
            }
            result.Add(1, new double[] { v1X, v1Y, v1Z });
            result.Add(2, new double[] { v2X, v2Y, v2Z });
            result.Add(3, new double[] { v3X, v3Y, v3Z });
            return result;
        }
        private double[,] StressStrainTransformationMatrix(Dictionary<int, double[]> unitVectors)
        {
            double[,] trMat = new double[6, 6];

            trMat[0, 0] = unitVectors[1][0] * unitVectors[1][0];
            trMat[0, 1] = unitVectors[1][1] * unitVectors[1][1];
            trMat[0, 2] = unitVectors[1][2] * unitVectors[1][2];

            trMat[1, 0] = unitVectors[2][0] * unitVectors[2][0];
            trMat[1, 1] = unitVectors[2][1] * unitVectors[2][1];
            trMat[1, 2] = unitVectors[2][2] * unitVectors[2][2];

            trMat[2, 0] = unitVectors[3][0] * unitVectors[3][0];
            trMat[2, 1] = unitVectors[3][1] * unitVectors[3][1];
            trMat[2, 2] = unitVectors[3][2] * unitVectors[3][2];
            //
            trMat[0, 3] = unitVectors[1][0] * unitVectors[1][1];
            trMat[0, 4] = unitVectors[1][1] * unitVectors[1][2];
            trMat[0, 5] = unitVectors[1][2] * unitVectors[1][0];

            trMat[1, 3] = unitVectors[2][0] * unitVectors[2][1];
            trMat[1, 4] = unitVectors[2][1] * unitVectors[2][2];
            trMat[1, 5] = unitVectors[2][2] * unitVectors[2][0];

            trMat[2, 3] = unitVectors[3][0] * unitVectors[3][1];
            trMat[2, 4] = unitVectors[3][1] * unitVectors[3][2];
            trMat[2, 5] = unitVectors[3][2] * unitVectors[3][0];
            //
            trMat[3, 0] = 2.0 * unitVectors[1][0] * unitVectors[2][0];
            trMat[4, 0] = 2.0 * unitVectors[2][0] * unitVectors[3][0];
            trMat[5, 0] = 2.0 * unitVectors[3][0] * unitVectors[1][0];

            trMat[3, 1] = 2.0 * unitVectors[1][1] * unitVectors[2][1];
            trMat[4, 1] = 2.0 * unitVectors[2][1] * unitVectors[3][1];
            trMat[5, 1] = 2.0 * unitVectors[3][1] * unitVectors[1][1];

            trMat[3, 2] = 2.0 * unitVectors[1][2] * unitVectors[2][2];
            trMat[4, 2] = 2.0 * unitVectors[2][2] * unitVectors[3][2];
            trMat[5, 2] = 2.0 * unitVectors[3][2] * unitVectors[1][2];
            //
            trMat[3, 3] = unitVectors[1][0] * unitVectors[2][1] + unitVectors[2][0] * unitVectors[1][1];
            trMat[4, 3] = unitVectors[2][0] * unitVectors[3][1] + unitVectors[3][0] * unitVectors[2][1];
            trMat[5, 3] = unitVectors[3][0] * unitVectors[1][1] + unitVectors[1][0] * unitVectors[3][1];

            trMat[3, 4] = unitVectors[1][1] * unitVectors[2][2] + unitVectors[2][1] * unitVectors[1][2];
            trMat[4, 4] = unitVectors[2][1] * unitVectors[3][2] + unitVectors[3][1] * unitVectors[2][2];
            trMat[5, 4] = unitVectors[3][1] * unitVectors[1][2] + unitVectors[1][1] * unitVectors[3][2];

            trMat[3, 5] = unitVectors[1][2] * unitVectors[2][0] + unitVectors[2][2] * unitVectors[1][0];
            trMat[4, 5] = unitVectors[2][2] * unitVectors[3][0] + unitVectors[3][2] * unitVectors[2][0];
            trMat[5, 5] = unitVectors[3][2] * unitVectors[1][0] + unitVectors[1][2] * unitVectors[3][0];
            return trMat;
        }
        public double[,] CreateGlobalStiffnessMatrix()
        {
            //var F = nodalForces();
            double[,] K = new double[45, 45];
            double[,] ElocalNodal = CalculateStressStrainMatrix(Properties.YoungMod, Properties.PoissonRatio);
            double[] xUpdated = UpdateNodalCoordinates(DisplacementVector);
            UpdateNodalGeometry(xUpdated);
            Dictionary<int, Dictionary<int, double[]>> nodesDirectionalCosines = CalculateNodalUnitVectors(xUpdated, Properties.DegenerateShellNodalThicknesses);
            for (int i = 0; i <= 2; i++)
            {
                for (int j = 0; j <= 2; j++)
                {
                    for (int k = 0; k <= 2; k++)
                    {
                        double[] gP = GaussPoints(i, j, k).Item1;
                        double[] gW = GaussPoints(i, j, k).Item2;
                        Dictionary<string, double[]> localdN = CalculateShapeFunctionsLocalDerivatives(gP);
                        var N = CalculateShapeFunctions(gP[0], gP[1]);
                        double[,] J = CalculateJacobian(gP[2], localdN, N, xUpdated, nodesDirectionalCosines);
                        //double[,] invJ = CalculateInverseJacobian(J).Item1;
                        double detJ = CalculateInverseJacobian(J).Item2;
                        double[,] B1 = CalculateInverseJacobian(J).Item3;
                        double[,] B = CalculateBMatrix(gP[2], B1, localdN, N, nodesDirectionalCosines);
                        double[,] tMatrix = StressStrainTransformationMatrix(unitVectors(N, nodesDirectionalCosines));
                        double[,] E = MatrixOperations.MatrixProduct(MatrixOperations.Transpose(tMatrix),
                            MatrixOperations.MatrixProduct(ElocalNodal, tMatrix));
                        K = MatrixOperations.MatrixAddition(K, MatrixOperations.ScalarMatrixProductNew(detJ * gW[0] * gW[1] * gW[2],
                            MatrixOperations.MatrixProduct(MatrixOperations.Transpose(B), MatrixOperations.MatrixProduct(E, B))));

                    }
                }
            }
            return K;
        }

        public double[,] CreateMassMatrix()
        {
            double[,] M = new double[45, 45];
            double[,] Mtr = new double[27, 27];
            double[] xUpdated = UpdateNodalCoordinates(DisplacementVector);
            UpdateNodalGeometry(xUpdated);
            double[,] noadalThicknessesMatrix = new double[,] {
                {Properties.DegenerateShellNodalThicknesses[1],
                Properties.DegenerateShellNodalThicknesses[2],
                Properties.DegenerateShellNodalThicknesses[3]},
            {Properties.DegenerateShellNodalThicknesses[8],
                Properties.DegenerateShellNodalThicknesses[9],
                Properties.DegenerateShellNodalThicknesses[4]},
            {Properties.DegenerateShellNodalThicknesses[7],
                Properties.DegenerateShellNodalThicknesses[6],
                Properties.DegenerateShellNodalThicknesses[5]} };
            for(int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    var lobatto = LobattoPoints(i, j);
                    double[] lobattoPoints = lobatto.Item1;
                    double[] lobattoWeights = lobatto.Item2;
                    var shapeFunctions = CalculateShapeFunctionMatrix(lobattoPoints[0], lobattoPoints[1]);
                    double[,] N = shapeFunctions.Item1;
                    Dictionary<int, double> dictionaryN = shapeFunctions.Item2;
                    Dictionary<string, double[]> localdN = CalculateShapeFunctionsLocalDerivatives(lobattoPoints);
                    Dictionary<int, Dictionary<int, double[]>> nodesDirectionalCosines = CalculateNodalUnitVectors(xUpdated, Properties.DegenerateShellNodalThicknesses);
                    double[,] J = CalculateJacobian(0.0, localdN, dictionaryN,
                        xUpdated, nodesDirectionalCosines);
                    double detJ = CalculateInverseJacobian(J).Item2;
                    Mtr = MatrixOperations.MatrixAddition(Mtr, MatrixOperations.ScalarMatrixProductNew(detJ * lobattoWeights[0] * lobattoWeights[1] * noadalThicknessesMatrix[i,j],
                        MatrixOperations.MatrixProduct(MatrixOperations.Transpose(N), MatrixOperations.ScalarMatrixProductNew(Properties.Density, N))));
                }
            }
            Dictionary<int, double> I = new Dictionary<int, double>();
            for(int i = 0; i < 9; i++)
            {
                double inertia = Mtr[i * 3, 3 * i] * Math.Pow(Properties.DegenerateShellNodalThicknesses[i + 1],
                    2.0) / 4.0;
                I.Add(i + 1, inertia);
            }
            for(int i = 0; i < 9; i++)
            {
                M[5 * i, 5 * i] = Mtr[3 * i, 3 * i];
                M[5 * i + 1, 5 * i + 1] = Mtr[3 * i + 1, 3 * i + 1];
                M[5 * i + 2, 5 * i + 2] = Mtr[3 * i + 2, 3 * i + 2];
                M[5 * i + 3, 5 * i + 3] = I[i + 1];
                M[5 * i + 4, 5 * i + 4] = I[i + 1];
            }
            //double alpha = 0.10;
            //double thickness = Properties.DegenerateShellNodalThicknesses.Max(t => t.Value);
            //double[,] M = MatrixOperations.CreateDiagonalMatrix(45, Properties.Density
            //    * thickness * alpha * alpha);
            return M;
        }

        public double[,] CreateDampingMatrix()
        {
            return new double[45, 45];
        }

        public double[] CreateInternalGlobalForcesVector()
        {
            double[] fInt = new double[45];
            double[,] ElocalNodal = CalculateStressStrainMatrix(Properties.YoungMod, Properties.PoissonRatio);
            double[] xUpdated = UpdateNodalCoordinates(DisplacementVector);
            UpdateNodalGeometry(xUpdated);
            Dictionary<int, Dictionary<int, double[]>> nodesDirectionalCosines = CalculateNodalUnitVectors(xUpdated, Properties.DegenerateShellNodalThicknesses);
            for (int i = 0; i <= 2; i++)
            {
                for (int j = 0; j <= 2; j++)
                {
                    for (int k = 0; k <= 2; k++)
                    {
                        double[] gP = GaussPoints(i, j, k).Item1;
                        double[] gW = GaussPoints(i, j, k).Item2;
                        Dictionary<string, double[]> localdN = CalculateShapeFunctionsLocalDerivatives(gP);
                        var N = CalculateShapeFunctions(gP[0], gP[1]);
                        double[,] J = CalculateJacobian(gP[2], localdN, N, xUpdated, nodesDirectionalCosines);
                        //double[,] invJ = CalculateInverseJacobian(J).Item1;
                        double detJ = CalculateInverseJacobian(J).Item2;
                        double[,] B1 = CalculateInverseJacobian(J).Item3;
                        double[,] B = CalculateBMatrix(gP[2], B1, localdN, N, nodesDirectionalCosines);
                        double[] strainVector = CalculateStrainsVector(B);
                        double[,] tMatrix = StressStrainTransformationMatrix(unitVectors(N, nodesDirectionalCosines));
                        double[,] E = MatrixOperations.MatrixProduct(MatrixOperations.Transpose(tMatrix),
                            MatrixOperations.MatrixProduct(ElocalNodal, tMatrix));
                        double[] stressVector = CalculateStressVector(E, strainVector);
                        fInt = VectorOperations.VectorVectorAddition(fInt, VectorOperations.VectorScalarProductNew(
                            VectorOperations.MatrixVectorProduct(MatrixOperations.Transpose(B), stressVector), detJ * gW[0] * gW[1] * gW[2]));
                    }

                }
            }
            return fInt;
        }
        //private double[] nodalForces()
        //{
        //    double[] nF = new double[27];
        //    double[] extSurfLoad = new double[] { 0.0, 0.0, -1.0 };
        //    double[] xUpdated = UpdateNodalCoordinates(DisplacementVector);
        //    UpdateNodalGeometry(xUpdated);
        //    Dictionary<int, Dictionary<int, double[]>> nodesDirectionalCosines = CalculateNodalUnitVectors(xUpdated, Properties.DegenerateShellNodalThicknesses);
        //    int k = 1;
        //    for (int i = 0; i <= 2; i++)
        //    {
        //        for (int j = 0; j <= 2; j++)
        //        {
        //            double[] gP = GaussPoints(i, j, k).Item1;
        //            double[] gW = GaussPoints(i, j, k).Item2;
        //            Dictionary<string, double[]> localdN = CalculateShapeFunctionsLocalDerivatives(gP);
        //            var matrixN = CalculateShapeFunctionMatrix(gP[0], gP[1]).Item1;
        //            var N = CalculateShapeFunctionMatrix(gP[0], gP[1]).Item2;
        //            double[,] J = CalculateJacobian(1.0, localdN, N, xUpdated, nodesDirectionalCosines);
        //            //double[,] invJ = CalculateInverseJacobian(J).Item1;
        //            double detJ = CalculateInverseJacobian(J).Item2;
        //            nF = VectorOperations.VectorVectorAddition(nF, VectorOperations.VectorScalarProduct(
        //                VectorOperations.MatrixVectorProduct(
        //            MatrixOperations.Transpose(matrixN), extSurfLoad), detJ * gW[0] * gW[1]));
        //        }
        //    }
        //    return nF;
        //}
        //private double[,] CalculateJacobian2(Dictionary<string, double[]> dN)
        //{
        //    double[,] jacobianMatrix = new double[2, 2];

        //    double[] xUpdated = UpdateNodalCoordinates(DisplacementVector);

        //    int k = 0;
        //    for (int i = 0; i < 9; i++)
        //    {
        //        jacobianMatrix[0, 0] = jacobianMatrix[0, 0] + xUpdated[k] * dN["ksi"][i];
        //        k = k + 2;
        //    }
        //    k = 1;
        //    for (int i = 0; i < 9; i++)
        //    {
        //        jacobianMatrix[0, 1] = jacobianMatrix[0, 1] + xUpdated[k] * dN["ksi"][i];
        //        k = k + 2;
        //    }

        //    k = 0;
        //    for (int i = 0; i < 9; i++)
        //    {
        //        jacobianMatrix[1, 0] = jacobianMatrix[1, 0] + xUpdated[k] * dN["ihta"][i];
        //        k = k + 2;
        //    }
        //    k = 1;
        //    for (int i = 0; i < 9; i++)
        //    {
        //        jacobianMatrix[1, 1] = jacobianMatrix[1, 1] + xUpdated[k] * dN["ihta"][i];
        //        k = k + 2;
        //    }
        //    return jacobianMatrix;
        //}

        //private Tuple<double[,], double> CalculateInverseJacobian2(double[,] jacobianMatrix)
        //{
        //    double[,] jacobianInverseMatrix = new double[2, 2];

        //    double detj = jacobianMatrix[0, 0] * jacobianMatrix[1, 1] - jacobianMatrix[0, 1] * jacobianMatrix[1, 0];

        //    jacobianInverseMatrix[0, 0] = jacobianMatrix[1, 1] / detj;
        //    jacobianInverseMatrix[0, 1] = -jacobianMatrix[0, 1] / detj;
        //    jacobianInverseMatrix[1, 0] = -jacobianMatrix[1, 0] / detj;
        //    jacobianInverseMatrix[1, 1] = jacobianMatrix[0, 0] / detj;

        //    return new Tuple<double[,], double>(jacobianInverseMatrix, detj);
        //}
    }
}

