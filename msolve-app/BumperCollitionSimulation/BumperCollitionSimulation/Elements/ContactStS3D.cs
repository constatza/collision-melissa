using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace BumperCollitionSimulation
{
    class ContactStS3D : IElement
    {
        public Dictionary<int, INode> Nodes { get; }
        public IElementProperties Properties { get; set; }
        public Dictionary<int, bool[]> ElementFreedomSignature { get; } = new Dictionary<int, bool[]>();
        public List<int> ElementFreedomList { get; set; }
        public double[] DisplacementVector { get; set; }
        public double[] AccelerationVector { get; set; }
        private double PenaltyFactor { get; set; }
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
        public List<double[]> GetStressFromElements(List<double[]> parametricCoordinatesVector)
        {
            List<double[]> StessVectorsList = new List<double[]>();
            StessVectorsList.Add(new double[] { 0.0, 0.0 });
            //double[,] E = CalculateStressStrainMatrix(Properties.YoungMod, Properties.PoissonRatio);
            //int count = parametricCoordinatesVector.Count;           
            //for (int i = 0; i < count; i++)
            //{
            //    double[] nodalParamCoord = parametricCoordinatesVector[i];
            //    Dictionary<string, double[]> localdN = CalculateShapeFunctionsLocalDerivatives(nodalParamCoord);
            //    double[,] J = CalculateJacobian(localdN);
            //    double[,] invJ = CalculateInverseJacobian(J).Item1;
            //    Dictionary<int, double[]> globaldN = CalculateShapeFunctionsGlobalDerivatives(localdN, invJ);
            //    double[,] B = CalculateBMatrix(globaldN);
            //    double[] strainVector = CalculateStrainsVector(B);
            //    double[] stressVector = CalculateStressVector(E, strainVector);
            //    StessVectorsList.Add(stressVector);
            //}
            return StessVectorsList;
        }
        public List<double[]> GetphysicalCoordinatesFromElements(List<double[]> parametricCoordinatesVector)
        {
            List<double[]> PositionVectorsList = new List<double[]>();
            //double[] xUpdated = UpdateNodalCoordinates(DisplacementVector);
            PositionVectorsList.Add(new double[] { 0.0, 0.0 });
            //int count = parametricCoordinatesVector.Count;
            //for (int i = 0; i < count; i++)
            //{
            //    double[] parametricCoordinatesVec = parametricCoordinatesVector[i];
            //    double[] positionVector = VectorOperations.MatrixVectorProduct(CalculateShapeFunctionMatrix(parametricCoordinatesVec[0], parametricCoordinatesVec[1], parametricCoordinatesVec[2]), xUpdated);
            //    PositionVectorsList.Add(positionVector);
            //}
            return PositionVectorsList;
        }
        public ContactStS3D(IElementProperties properties, Dictionary<int, INode> nodes)
        {
            Properties = properties;
            this.Nodes = nodes;
            int amountOfNodes = (Properties.SlaveSegmentPolynomialDegree + 1) * (Properties.SlaveSegmentPolynomialDegree + 1) +
                (Properties.MasterSegmentPolynomialDegree + 1) * (Properties.MasterSegmentPolynomialDegree + 1);
            for (int i = 1; i <= amountOfNodes; i++)
            {
                ElementFreedomSignature[i] = new bool[] { true, true, true, false, false, false };

            }
            DisplacementVector = new double[3 * amountOfNodes];
            PenaltyFactor = properties.YoungMod * properties.PenaltyFactorRatio;
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
        public Dictionary<int, INode> NodesAtFinalState()
        {
            Dictionary<int, INode> finalNodes = new Dictionary<int, INode>();
            int amountOfNodes = (Properties.SlaveSegmentPolynomialDegree + 1) * (Properties.SlaveSegmentPolynomialDegree + 1) +
                (Properties.MasterSegmentPolynomialDegree + 1) * (Properties.MasterSegmentPolynomialDegree + 1);
            int dofCount = new int();
            for (int i = 1; i <= amountOfNodes; i++)
            {
                finalNodes[i] = new Node(Nodes[i].XCoordinate + DisplacementVector[dofCount], Nodes[i].YCoordinate + DisplacementVector[dofCount + 1],
                    Nodes[i].ZCoordinate + DisplacementVector[dofCount + 2]);
                dofCount += 3;
            }
            return finalNodes;
        }
        private double[] NodalXUpdated()
        {
            int amountOfNodes = (Properties.SlaveSegmentPolynomialDegree + 1) * (Properties.SlaveSegmentPolynomialDegree + 1) +
                (Properties.MasterSegmentPolynomialDegree + 1) * (Properties.MasterSegmentPolynomialDegree + 1);
            double[] x = new double[3 * amountOfNodes];
            for (int i = 1; i <= amountOfNodes; i++)
            {
                int count = (i - 1) * 3;
                x[count] = Nodes[i].XCoordinate + DisplacementVector[count];
                x[count + 1] = Nodes[i].YCoordinate + DisplacementVector[count + 1];
                x[count + 2] = Nodes[i].ZCoordinate + DisplacementVector[count + 2];
            }
            return x;
        }
        public List<double[]> GetStressVector()
        {
            List<double[]> l = new List<double[]>();
            l.Add(new double[] { 0.0, 0.0, 0.0 });
            return l;
        }
        public List<double[]> GetStrainVector()
        {
            List<double[]> l = new List<double[]>();
            l.Add(new double[] { 0.0, 0.0, 0.0 });
            return l;
        }

        public List<double[]> GetGaussPointsInPhysicalSpace()
        {
            List<double[]> l = new List<double[]>();
            l.Add(new double[] { 0.0, 0.0 });
            return l;
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
        public double ClosestPointProjection()
        {
            throw new Exception("Alternative method <Project> has been used for higher order elements");
        }

        private Tuple<double[,], Tuple<double[,], double[,], double[,], double[,], double[,]>, Tuple<double[,], double[,], double[,], double[,]>> CalculatePositionMatrix(double ksi1, double ksi2, double ksi3, double ksi4)
        {
            if (Properties.SlaveSegmentPolynomialDegree == 1 && Properties.MasterSegmentPolynomialDegree == 1)
            {
                double N1 = 1.0 / 4.0 * (1.0 - ksi1) * (1.0 - ksi2);
                double N2 = 1.0 / 4.0 * (1.0 + ksi1) * (1.0 - ksi2);
                double N3 = 1.0 / 4.0 * (1.0 + ksi1) * (1.0 + ksi2);
                double N4 = 1.0 / 4.0 * (1.0 - ksi1) * (1.0 + ksi2);
                double N5 = 1.0 / 4.0 * (1.0 - ksi3) * (1.0 - ksi4);
                double N6 = 1.0 / 4.0 * (1.0 + ksi3) * (1.0 - ksi4);
                double N7 = 1.0 / 4.0 * (1.0 + ksi3) * (1.0 + ksi4);
                double N8 = 1.0 / 4.0 * (1.0 - ksi3) * (1.0 + ksi4);

                double dN11 = -1.0 / 4.0 * (1.0 - ksi2);
                double dN21 = 1.0 / 4.0 * (1.0 - ksi2);
                double dN31 = 1.0 / 4.0 * (1.0 + ksi2);
                double dN41 = -1.0 / 4.0 * (1.0 + ksi2);

                double dN12 = -1.0 / 4.0 * (1.0 - ksi1);
                double dN22 = -1.0 / 4.0 * (1.0 + ksi1);
                double dN32 = 1.0 / 4.0 * (1.0 + ksi1);
                double dN42 = 1.0 / 4.0 * (1.0 - ksi1);

                double dN53 = -1.0 / 4.0 * (1.0 - ksi4);
                double dN63 = 1.0 / 4.0 * (1.0 - ksi4);
                double dN73 = 1.0 / 4.0 * (1.0 + ksi4);
                double dN83 = -1.0 / 4.0 * (1.0 + ksi4);

                double dN54 = -1.0 / 4.0 * (1.0 - ksi3);
                double dN64 = -1.0 / 4.0 * (1.0 + ksi3);
                double dN74 = 1.0 / 4.0 * (1.0 + ksi3);
                double dN84 = 1.0 / 4.0 * (1.0 - ksi3);

                double dN112 = 1.0 / 4.0;
                double dN212 = -1.0 / 4.0;
                double dN312 = 1.0 / 4.0;
                double dN412 = -1.0 / 4.0;

                double[,] aMatrix = new double[,]
                    {
                    { -N1 ,0.0, 0.0 ,-N2 ,0.0 ,0.0 ,-N3, 0.0, 0.0, -N4, 0.0, 0.0, N5 ,0.0, 0.0 , N6, 0.0 ,0.0 , N7, 0.0, 0.0, N8, 0.0, 0.0 },
                    { 0.0, -N1 ,0.0, 0.0 ,-N2 ,0.0 ,0.0 ,-N3, 0.0, 0.0, -N4, 0.0, 0.0, N5 ,0.0, 0.0 , N6, 0.0 ,0.0 , N7, 0.0, 0.0, N8, 0.0 },
                    { 0.0, 0.0, -N1 ,0.0, 0.0 ,-N2 ,0.0 ,0.0 ,-N3, 0.0, 0.0, -N4, 0.0, 0.0, N5 ,0.0, 0.0 , N6, 0.0 ,0.0 , N7, 0.0, 0.0, N8 }
                    };

                double[,] da1Matrix = new double[,]
                    {
                    { -dN11 ,0.0, 0.0 ,-dN21 ,0.0 ,0.0 ,-dN31, 0.0, 0.0, -dN41, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                    { 0.0, -dN11 ,0.0, 0.0 ,-dN21 ,0.0 ,0.0 ,-dN31, 0.0, 0.0, -dN41, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0 },
                    { 0.0, 0.0, -dN11 ,0.0, 0.0 ,-dN21 ,0.0 ,0.0 ,-dN31, 0.0, 0.0, -dN41, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0 }
                    };

                double[,] da11Matrix = new double[,]
                    {
                    { 0.0 ,0.0, 0.0 ,0.0 ,0.0 ,0.0 ,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                    { 0.0, 0.0 ,0.0, 0.0 ,0.0 ,0.0 ,0.0 ,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0 },
                    { 0.0, 0.0, 0.0 ,0.0, 0.0 ,0.0 ,0.0 ,0.0 ,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0 }
                    };

                double[,] da12Matrix = new double[,]
                    {
                    { -dN112 ,0.0, 0.0 ,-dN212 ,0.0 ,0.0 ,-dN312, 0.0, 0.0, -dN412, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                    { 0.0, -dN112 ,0.0, 0.0 ,-dN212 ,0.0 ,0.0 ,-dN312, 0.0, 0.0, -dN412, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0 },
                    { 0.0, 0.0, -dN112 ,0.0, 0.0 ,-dN212 ,0.0 ,0.0 ,-dN312, 0.0, 0.0, -dN412, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0 }
                    };

                double[,] da2Matrix = new double[,]
                    {
                    { -dN12 ,0.0, 0.0 ,-dN22 ,0.0 ,0.0 ,-dN32, 0.0, 0.0, -dN42, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                    { 0.0, -dN12 ,0.0, 0.0 ,-dN22 ,0.0 ,0.0 ,-dN32, 0.0, 0.0, -dN42, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0 },
                    { 0.0, 0.0, -dN12 ,0.0, 0.0 ,-dN22 ,0.0 ,0.0 ,-dN32, 0.0, 0.0, -dN42, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0 }
                    };

                double[,] da22Matrix = new double[,]
                    {
                    { 0.0 ,0.0, 0.0 ,0.0 ,0.0 ,0.0 ,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                    { 0.0, 0.0 ,0.0, 0.0 ,0.0 ,0.0 ,0.0 ,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0 },
                    { 0.0, 0.0, 0.0 ,0.0, 0.0 ,0.0 ,0.0 ,0.0 ,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0 }
                    };

                double[,] da3Matrix = new double[,]
                    {
                    { 0.0, 0.0, 0.0 ,0.0 ,0.0 ,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, dN53 ,0.0, 0.0 , dN63, 0.0 ,0.0 , dN73, 0.0, 0.0, dN83, 0.0, 0.0 },
                    { 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, dN53 ,0.0, 0.0 , dN63, 0.0 ,0.0 , dN73, 0.0, 0.0, dN83, 0.0 },
                    { 0.0, 0.0, 0.0 ,0.0, 0.0 ,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, dN53 ,0.0, 0.0 , dN63, 0.0 ,0.0 , dN73, 0.0, 0.0, dN83 }
                    };

                double[,] da33Matrix = new double[,]
                    {
                    { 0.0 ,0.0, 0.0 ,0.0 ,0.0 ,0.0 ,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                    { 0.0, 0.0 ,0.0, 0.0 ,0.0 ,0.0 ,0.0 ,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0 },
                    { 0.0, 0.0, 0.0 ,0.0, 0.0 ,0.0 ,0.0 ,0.0 ,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0 }
                    };

                double[,] da4Matrix = new double[,]
                    {
                    { 0.0, 0.0, 0.0 ,0.0 ,0.0 ,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, dN54 ,0.0, 0.0 , dN64, 0.0 ,0.0 , dN74, 0.0, 0.0, dN84, 0.0, 0.0 },
                    { 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, dN54 ,0.0, 0.0 , dN64, 0.0 ,0.0 , dN74, 0.0, 0.0, dN84, 0.0 },
                    { 0.0, 0.0, 0.0 ,0.0, 0.0 ,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, dN54 ,0.0, 0.0 , dN64, 0.0 ,0.0 , dN74, 0.0, 0.0, dN84 }
                    };

                double[,] da44Matrix = new double[,]
                    {
                    { 0.0 ,0.0, 0.0 ,0.0 ,0.0 ,0.0 ,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                    { 0.0, 0.0 ,0.0, 0.0 ,0.0 ,0.0 ,0.0 ,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0 },
                    { 0.0, 0.0, 0.0 ,0.0, 0.0 ,0.0 ,0.0 ,0.0 ,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0 }
                    };

                var T1 = new Tuple<double[,], double[,], double[,], double[,], double[,]>(da1Matrix, da11Matrix, da2Matrix, da22Matrix, da12Matrix);
                var T2 = new Tuple<double[,], double[,], double[,], double[,]>(da3Matrix, da33Matrix, da4Matrix, da44Matrix);
                return new Tuple<double[,], Tuple<double[,], double[,], double[,], double[,], double[,]>, Tuple<double[,], double[,], double[,], double[,]>>(aMatrix, T1, T2);
            }
            else if (Properties.MasterSegmentPolynomialDegree == 1 && Properties.SlaveSegmentPolynomialDegree == 2)
            {
                double N1 = 1.0 / 4.0 * (1.0 - ksi1) * (1.0 - ksi2);
                double N2 = 1.0 / 4.0 * (1.0 + ksi1) * (1.0 - ksi2);
                double N3 = 1.0 / 4.0 * (1.0 + ksi1) * (1.0 + ksi2);
                double N4 = 1.0 / 4.0 * (1.0 - ksi1) * (1.0 + ksi2);
                double N5 = 1.0 / 4.0 * ksi3 * ksi4 * (1 - ksi3) * (1 - ksi4);
                double N6 = -1.0 / 2.0 * ksi4 * (1 + ksi3) * (1 - ksi3) * (1 - ksi4);
                double N7 = -1.0 / 4.0 * ksi3 * ksi4 * (1 + ksi3) * (1 - ksi4);
                double N8 = 1.0 / 2.0 * ksi3 * (1 + ksi3) * (1 + ksi4) * (1 - ksi4);
                double N9 = 1.0 / 4.0 * ksi3 * ksi4 * (1 + ksi3) * (1 + ksi4);
                double N10 = 1.0 / 2.0 * ksi4 * (1 - ksi3) * (1 + ksi3) * (1 + ksi4);
                double N11 = -1.0 / 4.0 * ksi3 * ksi4 * (1 - ksi3) * (1 + ksi4);
                double N12 = -1.0 / 2.0 * ksi3 * (1 - ksi3) * (1 + ksi4) * (1 - ksi4);
                double N13 = (1 - Math.Pow(ksi3, 2)) * (1 + ksi4) * (1 - ksi4);

                double dN11 = -1.0 / 4.0 * (1.0 - ksi2);
                double dN21 = 1.0 / 4.0 * (1.0 - ksi2);
                double dN31 = 1.0 / 4.0 * (1.0 + ksi2);
                double dN41 = -1.0 / 4.0 * (1.0 + ksi2);

                double dN12 = -1.0 / 4.0 * (1.0 - ksi1);
                double dN22 = -1.0 / 4.0 * (1.0 + ksi1);
                double dN32 = 1.0 / 4.0 * (1.0 + ksi1);
                double dN42 = 1.0 / 4.0 * (1.0 - ksi1);

                double dN112 = 1.0 / 4.0;
                double dN212 = -1.0 / 4.0;
                double dN312 = 1.0 / 4.0;
                double dN412 = -1.0 / 4.0;

                double dN53 = 1.0 / 4.0 * ksi4 * (1 - ksi4) * (1 - 2 * ksi3);
                double dN63 = (1 - ksi4) * ksi3 * ksi4;
                double dN73 = -1.0 / 4.0 * ksi4 * (1 - ksi4) * (1 + 2 * ksi3);
                double dN83 = 1.0 / 2.0 * (1 + 2 * ksi3) * (1 - Math.Pow(ksi4, 2));
                double dN93 = 1.0 / 4.0 * ksi4 * (1 + ksi4) * (2 * ksi3 + 1);
                double dN103 = -(ksi4 + 1) * ksi3 * ksi4;
                double dN113 = -1.0 / 4.0 * ksi4 * (1 + ksi4) * (1 - 2 * ksi3);
                double dN123 = -1.0 / 2.0 * (1 - 2 * ksi3) * (1 - Math.Pow(ksi4, 2));
                double dN133 = -2 * ksi3 * (1 - Math.Pow(ksi4, 2));

                double dN533 = -1.0 / 2.0 * ksi4 * (1 - ksi4);
                double dN633 = (1 - ksi4) * ksi4;
                double dN733 = -1.0 / 2.0 * ksi4 * (1 - ksi4);
                double dN833 = 1 - Math.Pow(ksi4, 2);
                double dN933 = 1.0 / 2.0 * ksi4 * (1 + ksi4);
                double dN1033 = -(ksi4 + 1) * ksi4;
                double dN1133 = 1.0 / 2.0 * ksi4 * (1 + ksi4);
                double dN1233 = 1 - Math.Pow(ksi4, 2);
                double dN1333 = -2 * (1 - Math.Pow(ksi4, 2));

                double dN54 = 1.0 / 4.0 * ksi3 * (1 - ksi3) * (1 - 2 * ksi4);
                double dN64 = -1.0 / 2.0 * (1 - Math.Pow(ksi3, 2)) * (1 - 2 * ksi4);
                double dN74 = -1.0 / 4.0 * ksi3 * (1 + ksi3) * (1 - 2 * ksi4);
                double dN84 = -(ksi3 + 1) * ksi4 * ksi3;
                double dN94 = 1.0 / 4.0 * ksi3 * (1 + ksi3) * (1 + 2 * ksi4);
                double dN104 = 1.0 / 2.0 * (1 - Math.Pow(ksi3, 2)) * (1 + 2 * ksi4);
                double dN114 = -1.0 / 4.0 * ksi3 * (1 - ksi3) * (1 + 2 * ksi4);
                double dN124 = (1 - ksi3) * ksi3 * ksi4;
                double dN134 = -2 * ksi4 * (1 - Math.Pow(ksi3, 2));

                double dN544 = -1.0 / 2.0 * ksi3 * (1 - ksi3);
                double dN644 = 1 - Math.Pow(ksi3, 2);
                double dN744 = 1.0 / 2.0 * ksi3 * (1 + ksi3);
                double dN844 = -(ksi3 + 1) * ksi3;
                double dN944 = 1.0 / 2.0 * ksi3 * (1 + ksi3);
                double dN1044 = 1 - Math.Pow(ksi3, 2);
                double dN1144 = -1.0 / 2.0 * ksi3 * (1 - ksi3);
                double dN1244 = (1 - ksi3) * ksi3;
                double dN1344 = -2 * (1 - Math.Pow(ksi3, 2));

                double[,] aMatrix = new double[,]
                    {
                    { -N1 ,0.0, 0.0 ,-N2 ,0.0 ,0.0 ,-N3, 0.0, 0.0, -N4, 0.0, 0.0, N5 ,0.0, 0.0 , N6, 0.0 ,0.0 , N7, 0.0, 0.0, N8, 0.0, 0.0, N9, 0.0, 0.0, N10, 0.0, 0.0, N11, 0.0, 0.0, N12, 0.0, 0.0, N13, 0.0, 0.0 },
                    { 0.0, -N1 ,0.0, 0.0 ,-N2 ,0.0 ,0.0 ,-N3, 0.0, 0.0, -N4, 0.0, 0.0, N5 ,0.0, 0.0 , N6, 0.0 ,0.0 , N7, 0.0, 0.0, N8, 0.0,0.0, N9, 0.0, 0.0, N10, 0.0, 0.0, N11, 0.0, 0.0, N12, 0.0, 0.0, N13, 0.0 },
                    { 0.0, 0.0, -N1 ,0.0, 0.0 ,-N2 ,0.0 ,0.0 ,-N3, 0.0, 0.0, -N4, 0.0, 0.0, N5 ,0.0, 0.0 , N6, 0.0 ,0.0 , N7, 0.0, 0.0, N8, 0.0,0.0, N9, 0.0, 0.0, N10, 0.0, 0.0, N11, 0.0, 0.0, N12, 0.0, 0.0, N13 }
                    };

                double[,] da1Matrix = new double[,]
                    {
                    { -dN11 ,0.0, 0.0 ,-dN21 ,0.0 ,0.0 ,-dN31, 0.0, 0.0, -dN41, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                    { 0.0, -dN11 ,0.0, 0.0 ,-dN21 ,0.0 ,0.0 ,-dN31, 0.0, 0.0, -dN41, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                    { 0.0, 0.0, -dN11 ,0.0, 0.0 ,-dN21 ,0.0 ,0.0 ,-dN31, 0.0, 0.0, -dN41, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
                    };

                double[,] da11Matrix = new double[,]
                    {
                    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
                    };

                double[,] da12Matrix = new double[,]
                    {
                    { -dN112 ,0.0, 0.0 ,-dN212 ,0.0 ,0.0 ,-dN312, 0.0, 0.0, -dN412, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                    { 0.0, -dN112 ,0.0, 0.0 ,-dN212 ,0.0 ,0.0 ,-dN312, 0.0, 0.0, -dN412, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                    { 0.0, 0.0, -dN112 ,0.0, 0.0 ,-dN212 ,0.0 ,0.0 ,-dN312, 0.0, 0.0, -dN412, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
                     };
                double[,] da2Matrix = new double[,]
                    {
                    { -dN12 ,0.0, 0.0 ,-dN22 ,0.0 ,0.0 ,-dN32, 0.0, 0.0, -dN42, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                    { 0.0, -dN12 ,0.0, 0.0 ,-dN22 ,0.0 ,0.0 ,-dN32, 0.0, 0.0, -dN42, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                    { 0.0, 0.0, -dN12 ,0.0, 0.0 ,-dN22 ,0.0 ,0.0 ,-dN32, 0.0, 0.0, -dN42, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
                    };

                double[,] da22Matrix = new double[,]
                    {
                    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
                    };

                double[,] da3Matrix = new double[,]
                    {
                    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, dN53 ,0.0, 0.0 , dN63, 0.0 ,0.0 , dN73, 0.0, 0.0, dN83, 0.0, 0.0, dN93, 0.0, 0.0, dN103, 0.0, 0.0, dN113, 0.0, 0.0, dN123, 0.0, 0.0, dN133, 0.0, 0.0 },
                    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, dN53 ,0.0, 0.0 , dN63, 0.0 ,0.0 , dN73, 0.0, 0.0, dN83, 0.0,0.0, dN93, 0.0, 0.0, dN103, 0.0, 0.0, dN113, 0.0, 0.0, dN123, 0.0, 0.0, dN133, 0.0 },
                    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, dN53 ,0.0, 0.0 , dN63, 0.0 ,0.0 , dN73, 0.0, 0.0, dN83, 0.0,0.0, dN93, 0.0, 0.0, dN103, 0.0, 0.0, dN113, 0.0, 0.0, dN123, 0.0, 0.0, dN133 }
                    };

                double[,] da33Matrix = new double[,]
                    {
                    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, dN533 ,0.0, 0.0 , dN633, 0.0 ,0.0 , dN733, 0.0, 0.0, dN833, 0.0, 0.0, dN933, 0.0, 0.0, dN1033, 0.0, 0.0, dN1133, 0.0, 0.0, dN1233, 0.0, 0.0, dN1333, 0.0, 0.0 },
                    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, dN533 ,0.0, 0.0 , dN633, 0.0 ,0.0 , dN733, 0.0, 0.0, dN833, 0.0,0.0, dN933, 0.0, 0.0, dN1033, 0.0, 0.0, dN1133, 0.0, 0.0, dN1233, 0.0, 0.0, dN1333, 0.0 },
                    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, dN533 ,0.0, 0.0 , dN633, 0.0 ,0.0 , dN733, 0.0, 0.0, dN833, 0.0,0.0, dN933, 0.0, 0.0, dN1033, 0.0, 0.0, dN1133, 0.0, 0.0, dN1233, 0.0, 0.0, dN1333 }
                    };

                double[,] da4Matrix = new double[,]
                    {
                    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, dN54 ,0.0, 0.0 , dN64, 0.0 ,0.0 , dN74, 0.0, 0.0, dN84, 0.0, 0.0, dN94, 0.0, 0.0, dN104, 0.0, 0.0, dN114, 0.0, 0.0, dN124, 0.0, 0.0, dN134, 0.0, 0.0 },
                    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, dN54 ,0.0, 0.0 , dN64, 0.0 ,0.0 , dN74, 0.0, 0.0, dN84, 0.0,0.0, dN94, 0.0, 0.0, dN104, 0.0, 0.0, dN114, 0.0, 0.0, dN124, 0.0, 0.0, dN134, 0.0 },
                    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, dN54 ,0.0, 0.0 , dN64, 0.0 ,0.0 , dN74, 0.0, 0.0, dN84, 0.0,0.0, dN94, 0.0, 0.0, dN104, 0.0, 0.0, dN114, 0.0, 0.0, dN124, 0.0, 0.0, dN134 }
                    };

                double[,] da44Matrix = new double[,]
                    {
                    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, dN544 ,0.0, 0.0 , dN644, 0.0 ,0.0 , dN744, 0.0, 0.0, dN844, 0.0, 0.0, dN944, 0.0, 0.0, dN1044, 0.0, 0.0, dN1144, 0.0, 0.0, dN1244, 0.0, 0.0, dN1344, 0.0, 0.0 },
                    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, dN544 ,0.0, 0.0 , dN644, 0.0 ,0.0 , dN744, 0.0, 0.0, dN844, 0.0,0.0, dN944, 0.0, 0.0, dN1044, 0.0, 0.0, dN1144, 0.0, 0.0, dN1244, 0.0, 0.0, dN1344, 0.0 },
                    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, dN544 ,0.0, 0.0 , dN644, 0.0 ,0.0 , dN744, 0.0, 0.0, dN844, 0.0,0.0, dN944, 0.0, 0.0, dN1044, 0.0, 0.0, dN1144, 0.0, 0.0, dN1244, 0.0, 0.0, dN1344 }
                    };

                var T1 = new Tuple<double[,], double[,], double[,], double[,], double[,]>(da1Matrix, da11Matrix, da2Matrix, da22Matrix, da12Matrix);
                var T2 = new Tuple<double[,], double[,], double[,], double[,]>(da3Matrix, da33Matrix, da4Matrix, da44Matrix);
                return new Tuple<double[,], Tuple<double[,], double[,], double[,], double[,], double[,]>, Tuple<double[,], double[,], double[,], double[,]>>(aMatrix, T1, T2);
            }
            else if (Properties.MasterSegmentPolynomialDegree == 2 && Properties.SlaveSegmentPolynomialDegree == 1)
            {
                double N1 = 1.0 / 4.0 * ksi1 * ksi2 * (1 - ksi1) * (1 - ksi2);
                double N2 = -1.0 / 2.0 * ksi2 * (1 + ksi1) * (1 - ksi1) * (1 - ksi2);
                double N3 = -1.0 / 4.0 * ksi1 * ksi2 * (1 + ksi1) * (1 - ksi2);
                double N4 = 1.0 / 2.0 * ksi1 * (1 + ksi1) * (1 + ksi2) * (1 - ksi2);
                double N5 = 1.0 / 4.0 * ksi1 * ksi2 * (1 + ksi1) * (1 + ksi2);
                double N6 = 1.0 / 2.0 * ksi2 * (1 - ksi1) * (1 + ksi1) * (1 + ksi2);
                double N7 = -1.0 / 4.0 * ksi1 * ksi2 * (1 - ksi1) * (1 + ksi2);
                double N8 = -1.0 / 2.0 * ksi1 * (1 - ksi1) * (1 + ksi2) * (1 - ksi2);
                double N9 = (1 - Math.Pow(ksi1, 2)) * (1 + ksi2) * (1 - ksi2);
                double N10 = 1.0 / 4.0 * (1.0 - ksi3) * (1.0 - ksi4);
                double N11 = 1.0 / 4.0 * (1.0 + ksi3) * (1.0 - ksi4);
                double N12 = 1.0 / 4.0 * (1.0 + ksi3) * (1.0 + ksi4);
                double N13 = 1.0 / 4.0 * (1.0 - ksi3) * (1.0 + ksi4);

                double dN11 = 1.0 / 4.0 * ksi2 * (1 - ksi2) * (1 - 2 * ksi1);
                double dN21 = (1 - ksi2) * ksi1 * ksi2;
                double dN31 = -1.0 / 4.0 * ksi2 * (1 - ksi2) * (1 + 2 * ksi1);
                double dN41 = 1.0 / 2.0 * (1 + 2 * ksi1) * (1 - Math.Pow(ksi2, 2));
                double dN51 = 1.0 / 4.0 * ksi2 * (1 + ksi2) * (2 * ksi1 + 1);
                double dN61 = -(ksi2 + 1) * ksi1 * ksi2;
                double dN71 = -1.0 / 4.0 * ksi2 * (1 + ksi2) * (1 - 2 * ksi1);
                double dN81 = -1.0 / 2.0 * (1 - 2 * ksi1) * (1 - Math.Pow(ksi2, 2));
                double dN91 = -2 * ksi1 * (1 - Math.Pow(ksi2, 2));

                double dN111 = -1.0 / 2.0 * ksi2 * (1 - ksi2);
                double dN211 = (1 - ksi2) * ksi2;
                double dN311 = -1.0 / 2.0 * ksi2 * (1 - ksi2);
                double dN411 = 1 - Math.Pow(ksi2, 2);
                double dN511 = 1.0 / 2.0 * ksi2 * (1 + ksi2);
                double dN611 = -(ksi2 + 1) * ksi2;
                double dN711 = 1.0 / 2.0 * ksi2 * (1 + ksi2);
                double dN811 = 1 - Math.Pow(ksi2, 2);
                double dN911 = -2 * (1 - Math.Pow(ksi2, 2));

                double dN12 = 1.0 / 4.0 * ksi1 * (1 - ksi1) * (1 - 2 * ksi2);
                double dN22 = -1.0 / 2.0 * (1 - Math.Pow(ksi1, 2)) * (1 - 2 * ksi2);
                double dN32 = -1.0 / 4.0 * ksi1 * (1 + ksi1) * (1 - 2 * ksi2);
                double dN42 = -(ksi1 + 1) * ksi2 * ksi1;
                double dN52 = 1.0 / 4.0 * ksi1 * (1 + ksi1) * (1 + 2 * ksi2);
                double dN62 = 1.0 / 2.0 * (1 - Math.Pow(ksi1, 2)) * (1 + 2 * ksi2);
                double dN72 = -1.0 / 4.0 * ksi1 * (1 - ksi1) * (1 + 2 * ksi2);
                double dN82 = (1 - ksi1) * ksi1 * ksi2;
                double dN92 = -2 * ksi2 * (1 - Math.Pow(ksi1, 2));

                double dN122 = -1.0 / 2.0 * ksi1 * (1 - ksi1);
                double dN222 = 1 - Math.Pow(ksi1, 2);
                double dN322 = 1.0 / 2.0 * ksi1 * (1 + ksi1);
                double dN422 = -(ksi1 + 1) * ksi1;
                double dN522 = 1.0 / 2.0 * ksi1 * (1 + ksi1);
                double dN622 = 1 - Math.Pow(ksi1, 2);
                double dN722 = -1.0 / 2.0 * ksi1 * (1 - ksi1);
                double dN822 = (1 - ksi1) * ksi1;
                double dN922 = -2 * (1 - Math.Pow(ksi1, 2));

                double dN112 = 1.0 / 4.0 * (1 - 2 * ksi2) * (1 - 2 * ksi1);
                double dN212 = (1 - 2 * ksi2) * ksi1;
                double dN312 = -1.0 / 4.0 * (1 - 2 * ksi2) * (1 + 2 * ksi1);
                double dN412 = -(1 + 2 * ksi1) * ksi2;
                double dN512 = 1.0 / 4.0 * (1 + 2 * ksi2) * (2 * ksi1 + 1);
                double dN612 = -(1 + 2 * ksi2) * ksi1;
                double dN712 = -1.0 / 4.0 * (1 + 2 * ksi2) * (1 - 2 * ksi1);
                double dN812 = (1 - 2 * ksi1) * ksi2;
                double dN912 = 4 * ksi1 * ksi2;

                double dN103 = -1.0 / 4.0 * (1.0 - ksi4);
                double dN113 = 1.0 / 4.0 * (1.0 - ksi4);
                double dN123 = 1.0 / 4.0 * (1.0 + ksi4);
                double dN133 = -1.0 / 4.0 * (1.0 + ksi4);

                double dN104 = -1.0 / 4.0 * (1.0 - ksi3);
                double dN114 = -1.0 / 4.0 * (1.0 + ksi3);
                double dN124 = 1.0 / 4.0 * (1.0 + ksi3);
                double dN134 = 1.0 / 4.0 * (1.0 - ksi3);

                double[,] aMatrix = new double[,]
                    {
                    { -N1 ,0.0, 0.0 ,-N2 ,0.0 ,0.0 ,-N3, 0.0, 0.0, -N4, 0.0, 0.0, -N5 ,0.0, 0.0 , -N6, 0.0 ,0.0 , -N7, 0.0, 0.0, -N8, 0.0, 0.0, -N9, 0.0, 0.0, N10, 0.0, 0.0, N11, 0.0, 0.0, N12, 0.0, 0.0, N13, 0.0, 0.0 },
                    { 0.0, -N1 ,0.0, 0.0 ,-N2 ,0.0 ,0.0 ,-N3, 0.0, 0.0, -N4, 0.0, 0.0, -N5 ,0.0, 0.0 , -N6, 0.0 ,0.0 , -N7, 0.0, 0.0, -N8, 0.0,0.0, -N9, 0.0, 0.0, N10, 0.0, 0.0, N11, 0.0, 0.0, N12, 0.0, 0.0, N13, 0.0 },
                    { 0.0, 0.0, -N1 ,0.0, 0.0 ,-N2 ,0.0 ,0.0 ,-N3, 0.0, 0.0, -N4, 0.0, 0.0, -N5 ,0.0, 0.0 , -N6, 0.0 ,0.0 , -N7, 0.0, 0.0, -N8, 0.0, 0.0, -N9, 0.0, 0.0, N10, 0.0, 0.0, N11, 0.0, 0.0, N12, 0.0, 0.0, N13 }
                    };

                double[,] da1Matrix = new double[,]
                    {
                    { -dN11 ,0.0, 0.0 ,-dN21 ,0.0 ,0.0 ,-dN31, 0.0, 0.0, -dN41, 0.0, 0.0, -dN51 ,0.0, 0.0 , -dN61, 0.0 ,0.0 , -dN71, 0.0, 0.0, -dN81, 0.0, 0.0, -dN91, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                    { 0.0, -dN11 ,0.0, 0.0 ,-dN21 ,0.0 ,0.0 ,-dN31, 0.0, 0.0, -dN41, 0.0, 0.0, -dN51 ,0.0, 0.0 , -dN61, 0.0 ,0.0 , -dN71, 0.0, 0.0, -dN81, 0.0,0.0, -dN91, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                    { 0.0, 0.0, -dN11 ,0.0, 0.0 ,-dN21 ,0.0 ,0.0 ,-dN31, 0.0, 0.0, -dN41, 0.0, 0.0, -dN51 ,0.0, 0.0 , -dN61, 0.0 ,0.0 , -dN71, 0.0, 0.0, -dN81, 0.0, 0.0, -dN91, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
                    };

                double[,] da11Matrix = new double[,]
                    {
                    { -dN111 ,0.0, 0.0 ,-dN211 ,0.0 ,0.0 ,-dN311, 0.0, 0.0, -dN411, 0.0, 0.0, -dN511 ,0.0, 0.0 , -dN611, 0.0 ,0.0 , -dN711, 0.0, 0.0, -dN811, 0.0, 0.0, -dN911, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                    { 0.0, -dN111 ,0.0, 0.0 ,-dN211 ,0.0 ,0.0 ,-dN311, 0.0, 0.0, -dN411, 0.0, 0.0, -dN511 ,0.0, 0.0 , -dN611, 0.0 ,0.0 , -dN711, 0.0, 0.0, -dN811, 0.0,0.0, -dN911, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                    { 0.0, 0.0, -dN111 ,0.0, 0.0 ,-dN211 ,0.0 ,0.0 ,-dN311, 0.0, 0.0, -dN411, 0.0, 0.0, -dN511 ,0.0, 0.0 , -dN611, 0.0 ,0.0 , -dN711, 0.0, 0.0, -dN811, 0.0, 0.0, -dN911, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
                    };

                double[,] da12Matrix = new double[,]
                    {
                    { -dN112 ,0.0, 0.0 ,-dN212 ,0.0 ,0.0 ,-dN312, 0.0, 0.0, -dN412, 0.0, 0.0, -dN512 ,0.0, 0.0 , -dN612, 0.0 ,0.0 , -dN712, 0.0, 0.0, -dN812, 0.0, 0.0, -dN912, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                    { 0.0, -dN112 ,0.0, 0.0 ,-dN212 ,0.0 ,0.0 ,-dN312, 0.0, 0.0, -dN412, 0.0, 0.0, -dN512 ,0.0, 0.0 , -dN612, 0.0 ,0.0 , -dN712, 0.0, 0.0, -dN812, 0.0,0.0, -dN912, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                    { 0.0, 0.0, -dN112 ,0.0, 0.0 ,-dN212 ,0.0 ,0.0 ,-dN312, 0.0, 0.0, -dN412, 0.0, 0.0, -dN512 ,0.0, 0.0 , -dN612, 0.0 ,0.0 , -dN712, 0.0, 0.0, -dN812, 0.0, 0.0, -dN912, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
                    };

                double[,] da2Matrix = new double[,]
                    {
                    { -dN12 ,0.0, 0.0 ,-dN22 ,0.0 ,0.0 ,-dN32, 0.0, 0.0, -dN42, 0.0, 0.0, -dN52 ,0.0, 0.0 , -dN62, 0.0 ,0.0 , -dN72, 0.0, 0.0, -dN82, 0.0, 0.0, -dN92, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                    { 0.0, -dN12 ,0.0, 0.0 ,-dN22 ,0.0 ,0.0 ,-dN32, 0.0, 0.0, -dN42, 0.0, 0.0, -dN52 ,0.0, 0.0 , -dN62, 0.0 ,0.0 , -dN72, 0.0, 0.0, -dN82, 0.0,0.0, -dN92, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                    { 0.0, 0.0, -dN12 ,0.0, 0.0 ,-dN22 ,0.0 ,0.0 ,-dN32, 0.0, 0.0, -dN42, 0.0, 0.0, -dN52 ,0.0, 0.0 , -dN62, 0.0 ,0.0 , -dN72, 0.0, 0.0, -dN82, 0.0, 0.0, -dN92, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
                    };

                double[,] da22Matrix = new double[,]
                    {
                    { -dN122 ,0.0, 0.0 ,-dN222 ,0.0 ,0.0 ,-dN322, 0.0, 0.0, -dN422, 0.0, 0.0, -dN522 ,0.0, 0.0 , -dN622, 0.0 ,0.0 , -dN722, 0.0, 0.0, -dN822, 0.0, 0.0, -dN922, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                    { 0.0, -dN122 ,0.0, 0.0 ,-dN222 ,0.0 ,0.0 ,-dN322, 0.0, 0.0, -dN422, 0.0, 0.0, -dN522 ,0.0, 0.0 , -dN622, 0.0 ,0.0 , -dN722, 0.0, 0.0, -dN822, 0.0,0.0, -dN922, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                    { 0.0, 0.0, -dN122 ,0.0, 0.0 ,-dN222 ,0.0 ,0.0 ,-dN322, 0.0, 0.0, -dN422, 0.0, 0.0, -dN522 ,0.0, 0.0 , -dN622, 0.0 ,0.0 , -dN722, 0.0, 0.0, -dN822, 0.0, 0.0, -dN922, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
                    };

                double[,] da3Matrix = new double[,]
                    {
                    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, dN103, 0.0, 0.0, dN113, 0.0, 0.0, dN123, 0.0, 0.0, dN133, 0.0, 0.0 },
                    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, dN103, 0.0, 0.0, dN113, 0.0, 0.0, dN123, 0.0, 0.0, dN133, 0.0 },
                    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, dN103, 0.0, 0.0, dN113, 0.0, 0.0, dN123, 0.0, 0.0, dN133 }
                    };

                double[,] da33Matrix = new double[,]
                    {
                    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
                    };
                double[,] da4Matrix = new double[,]
                    {
                    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, dN104, 0.0, 0.0, dN114, 0.0, 0.0, dN124, 0.0, 0.0, dN134, 0.0, 0.0 },
                    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, dN104, 0.0, 0.0, dN114, 0.0, 0.0, dN124, 0.0, 0.0, dN134, 0.0 },
                    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, dN104, 0.0, 0.0, dN114, 0.0, 0.0, dN124, 0.0, 0.0, dN134 }
                    };

                double[,] da44Matrix = new double[,]
                    {
                    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
                    };

                var T1 = new Tuple<double[,], double[,], double[,], double[,], double[,]>(da1Matrix, da11Matrix, da2Matrix, da22Matrix, da12Matrix);
                var T2 = new Tuple<double[,], double[,], double[,], double[,]>(da3Matrix, da33Matrix, da4Matrix, da44Matrix);
                return new Tuple<double[,], Tuple<double[,], double[,], double[,], double[,], double[,]>, Tuple<double[,], double[,], double[,], double[,]>>(aMatrix, T1, T2);
            }
            else
            {
                double N1 = 1.0 / 4.0 * ksi1 * ksi2 * (1 - ksi1) * (1 - ksi2);
                double N2 = -1.0 / 2.0 * ksi2 * (1 + ksi1) * (1 - ksi1) * (1 - ksi2);
                double N3 = -1.0 / 4.0 * ksi1 * ksi2 * (1 + ksi1) * (1 - ksi2);
                double N4 = 1.0 / 2.0 * ksi1 * (1 + ksi1) * (1 + ksi2) * (1 - ksi2);
                double N5 = 1.0 / 4.0 * ksi1 * ksi2 * (1 + ksi1) * (1 + ksi2);
                double N6 = 1.0 / 2.0 * ksi2 * (1 - ksi1) * (1 + ksi1) * (1 + ksi2);
                double N7 = -1.0 / 4.0 * ksi1 * ksi2 * (1 - ksi1) * (1 + ksi2);
                double N8 = -1.0 / 2.0 * ksi1 * (1 - ksi1) * (1 + ksi2) * (1 - ksi2);
                double N9 = (1 - Math.Pow(ksi1, 2)) * (1 + ksi2) * (1 - ksi2);
                double N10 = 1.0 / 4.0 * ksi3 * ksi4 * (1 - ksi3) * (1 - ksi4);
                double N11 = -1.0 / 2.0 * ksi4 * (1 + ksi3) * (1 - ksi3) * (1 - ksi4);
                double N12 = -1.0 / 4.0 * ksi3 * ksi4 * (1 + ksi3) * (1 - ksi4);
                double N13 = 1.0 / 2.0 * ksi3 * (1 + ksi3) * (1 + ksi4) * (1 - ksi4);
                double N14 = 1.0 / 4.0 * ksi3 * ksi4 * (1 + ksi3) * (1 + ksi4);
                double N15 = 1.0 / 2.0 * ksi4 * (1 - ksi3) * (1 + ksi3) * (1 + ksi4);
                double N16 = -1.0 / 4.0 * ksi3 * ksi4 * (1 - ksi3) * (1 + ksi4);
                double N17 = -1.0 / 2.0 * ksi3 * (1 - ksi3) * (1 + ksi4) * (1 - ksi4);
                double N18 = (1 - Math.Pow(ksi3, 2)) * (1 + ksi4) * (1 - ksi4);

                double dN11 = 1.0 / 4.0 * ksi2 * (1 - ksi2) * (1 - 2 * ksi1);
                double dN21 = (1 - ksi2) * ksi1 * ksi2;
                double dN31 = -1.0 / 4.0 * ksi2 * (1 - ksi2) * (1 + 2 * ksi1);
                double dN41 = 1.0 / 2.0 * (1 + 2 * ksi1) * (1 - Math.Pow(ksi2, 2));
                double dN51 = 1.0 / 4.0 * ksi2 * (1 + ksi2) * (2 * ksi1 + 1);
                double dN61 = -(ksi2 + 1) * ksi1 * ksi2;
                double dN71 = -1.0 / 4.0 * ksi2 * (1 + ksi2) * (1 - 2 * ksi1);
                double dN81 = -1.0 / 2.0 * (1 - 2 * ksi1) * (1 - Math.Pow(ksi2, 2));
                double dN91 = -2 * ksi1 * (1 - Math.Pow(ksi2, 2));

                double dN111 = -1.0 / 2.0 * ksi2 * (1 - ksi2);
                double dN211 = (1 - ksi2) * ksi2;
                double dN311 = -1.0 / 2.0 * ksi2 * (1 - ksi2);
                double dN411 = 1 - Math.Pow(ksi2, 2);
                double dN511 = 1.0 / 2.0 * ksi2 * (1 + ksi2);
                double dN611 = -(ksi2 + 1) * ksi2;
                double dN711 = 1.0 / 2.0 * ksi2 * (1 + ksi2);
                double dN811 = 1 - Math.Pow(ksi2, 2);
                double dN911 = -2 * (1 - Math.Pow(ksi2, 2));

                double dN112 = 1.0 / 4.0 * (1 - 2 * ksi2) * (1 - 2 * ksi1);
                double dN212 = (1 - 2 * ksi2) * ksi1;
                double dN312 = -1.0 / 4.0 * (1 - 2 * ksi2) * (1 + 2 * ksi1);
                double dN412 = -(1 + 2 * ksi1) * ksi2;
                double dN512 = 1.0 / 4.0 * (1 + 2 * ksi2) * (2 * ksi1 + 1);
                double dN612 = -(1 + 2 * ksi2) * ksi1;
                double dN712 = -1.0 / 4.0 * (1 + 2 * ksi2) * (1 - 2 * ksi1);
                double dN812 = (1 - 2 * ksi1) * ksi2;
                double dN912 = 4 * ksi1 * ksi2;

                double dN12 = 1.0 / 4.0 * ksi1 * (1 - ksi1) * (1 - 2 * ksi2);
                double dN22 = -1.0 / 2.0 * (1 - Math.Pow(ksi1, 2)) * (1 - 2 * ksi2);
                double dN32 = -1.0 / 4.0 * ksi1 * (1 + ksi1) * (1 - 2 * ksi2);
                double dN42 = -(ksi1 + 1) * ksi2 * ksi1;
                double dN52 = 1.0 / 4.0 * ksi1 * (1 + ksi1) * (1 + 2 * ksi2);
                double dN62 = 1.0 / 2.0 * (1 - Math.Pow(ksi1, 2)) * (1 + 2 * ksi2);
                double dN72 = -1.0 / 4.0 * ksi1 * (1 - ksi1) * (1 + 2 * ksi2);
                double dN82 = (1 - ksi1) * ksi1 * ksi2;
                double dN92 = -2 * ksi2 * (1 - Math.Pow(ksi1, 2));

                double dN122 = -1.0 / 2.0 * ksi1 * (1 - ksi1);
                double dN222 = 1 - Math.Pow(ksi1, 2);
                double dN322 = 1.0 / 2.0 * ksi1 * (1 + ksi1);
                double dN422 = -(ksi1 + 1) * ksi1;
                double dN522 = 1.0 / 2.0 * ksi1 * (1 + ksi1);
                double dN622 = 1 - Math.Pow(ksi1, 2);
                double dN722 = -1.0 / 2.0 * ksi1 * (1 - ksi1);
                double dN822 = (1 - ksi1) * ksi1;
                double dN922 = -2 * (1 - Math.Pow(ksi1, 2));

                double dN103 = 1.0 / 4.0 * ksi4 * (1 - ksi4) * (1 - 2 * ksi3);
                double dN113 = (1 - ksi4) * ksi3 * ksi4;
                double dN123 = -1.0 / 4.0 * ksi4 * (1 - ksi4) * (1 + 2 * ksi3);
                double dN133 = 1.0 / 2.0 * (1 + 2 * ksi3) * (1 - Math.Pow(ksi4, 2));
                double dN143 = 1.0 / 4.0 * ksi4 * (1 + ksi4) * (2 * ksi3 + 1);
                double dN153 = -(ksi4 + 1) * ksi3 * ksi4;
                double dN163 = -1.0 / 4.0 * ksi4 * (1 + ksi4) * (1 - 2 * ksi3);
                double dN173 = -1.0 / 2.0 * (1 - 2 * ksi3) * (1 - Math.Pow(ksi4, 2));
                double dN183 = -2 * ksi3 * (1 - Math.Pow(ksi4, 2));

                double dN1033 = -1.0 / 2.0 * ksi4 * (1 - ksi4);
                double dN1133 = (1 - ksi4) * ksi4;
                double dN1233 = -1.0 / 2.0 * ksi4 * (1 - ksi4);
                double dN1333 = 1 - Math.Pow(ksi4, 2);
                double dN1433 = 1.0 / 2.0 * ksi4 * (1 + ksi4);
                double dN1533 = -(ksi4 + 1) * ksi4;
                double dN1633 = 1.0 / 2.0 * ksi4 * (1 + ksi4);
                double dN1733 = 1 - Math.Pow(ksi4, 2);
                double dN1833 = -2 * (1 - Math.Pow(ksi4, 2));

                double dN104 = 1.0 / 4.0 * ksi3 * (1 - ksi3) * (1 - 2 * ksi4);
                double dN114 = -1.0 / 2.0 * (1 - Math.Pow(ksi3, 2)) * (1 - 2 * ksi4);
                double dN124 = -1.0 / 4.0 * ksi3 * (1 + ksi3) * (1 - 2 * ksi4);
                double dN134 = -(ksi3 + 1) * ksi4 * ksi3;
                double dN144 = 1.0 / 4.0 * ksi3 * (1 + ksi3) * (1 + 2 * ksi4);
                double dN154 = 1.0 / 2.0 * (1 - Math.Pow(ksi3, 2)) * (1 + 2 * ksi4);
                double dN164 = -1.0 / 4.0 * ksi3 * (1 - ksi3) * (1 + 2 * ksi4);
                double dN174 = (1 - ksi3) * ksi3 * ksi4;
                double dN184 = -2 * ksi4 * (1 - Math.Pow(ksi3, 2));

                double dN1044 = -1.0 / 2.0 * ksi3 * (1 - ksi3);
                double dN1144 = 1 - Math.Pow(ksi3, 2);
                double dN1244 = 1.0 / 2.0 * ksi3 * (1 + ksi3);
                double dN1344 = -(ksi3 + 1) * ksi3;
                double dN1444 = 1.0 / 2.0 * ksi3 * (1 + ksi3);
                double dN1544 = 1 - Math.Pow(ksi3, 2);
                double dN1644 = -1.0 / 2.0 * ksi3 * (1 - ksi3);
                double dN1744 = (1 - ksi3) * ksi3;
                double dN1844 = -2 * (1 - Math.Pow(ksi3, 2));

                double[,] aMatrix = new double[,]
                    {
                    { -N1 ,0.0, 0.0 ,-N2 ,0.0 ,0.0 ,-N3, 0.0, 0.0, -N4, 0.0, 0.0, -N5 ,0.0, 0.0 , -N6, 0.0 ,0.0 , -N7, 0.0, 0.0, -N8, 0.0, 0.0, -N9, 0.0, 0.0, N10, 0.0, 0.0, N11, 0.0, 0.0, N12, 0.0, 0.0, N13, 0.0, 0.0, N14, 0.0, 0.0, N15, 0.0, 0.0, N16, 0.0, 0.0, N17, 0.0, 0.0, N18, 0.0, 0.0 },
                    { 0.0, -N1 ,0.0, 0.0 ,-N2 ,0.0 ,0.0 ,-N3, 0.0, 0.0, -N4, 0.0, 0.0, -N5 ,0.0, 0.0 , -N6, 0.0 ,0.0 , -N7, 0.0, 0.0, -N8, 0.0,0.0, -N9, 0.0, 0.0, N10, 0.0, 0.0, N11, 0.0, 0.0, N12, 0.0, 0.0, N13, 0.0, 0.0, N14, 0.0, 0.0, N15, 0.0, 0.0, N16, 0.0, 0.0, N17, 0.0, 0.0, N18, 0.0 },
                    { 0.0, 0.0, -N1 ,0.0, 0.0 ,-N2 ,0.0 ,0.0 ,-N3, 0.0, 0.0, -N4, 0.0, 0.0, -N5 ,0.0, 0.0 , -N6, 0.0 ,0.0 , -N7, 0.0, 0.0, -N8, 0.0, 0.0, -N9, 0.0, 0.0, N10, 0.0, 0.0, N11, 0.0, 0.0, N12, 0.0, 0.0, N13, 0.0, 0.0, N14, 0.0, 0.0, N15, 0.0, 0.0, N16, 0.0, 0.0, N17, 0.0, 0.0, N18 }
                    };

                double[,] da1Matrix = new double[,]
                    {
                    { -dN11 ,0.0, 0.0 ,-dN21 ,0.0 ,0.0 ,-dN31, 0.0, 0.0, -dN41, 0.0, 0.0, -dN51 ,0.0, 0.0 , -dN61, 0.0 ,0.0 , -dN71, 0.0, 0.0, -dN81, 0.0, 0.0, -dN91, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                    { 0.0, -dN11 ,0.0, 0.0 ,-dN21 ,0.0 ,0.0 ,-dN31, 0.0, 0.0, -dN41, 0.0, 0.0, -dN51 ,0.0, 0.0 , -dN61, 0.0 ,0.0 , -dN71, 0.0, 0.0, -dN81, 0.0,0.0, -dN91, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                    { 0.0, 0.0, -dN11 ,0.0, 0.0 ,-dN21 ,0.0 ,0.0 ,-dN31, 0.0, 0.0, -dN41, 0.0, 0.0, -dN51 ,0.0, 0.0 , -dN61, 0.0 ,0.0 , -dN71, 0.0, 0.0, -dN81, 0.0, 0.0, -dN91, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
                    };

                double[,] da11Matrix = new double[,]
                    {
                    { -dN111 ,0.0, 0.0 ,-dN211 ,0.0 ,0.0 ,-dN311, 0.0, 0.0, -dN411, 0.0, 0.0, -dN511 ,0.0, 0.0 , -dN611, 0.0 ,0.0 , -dN711, 0.0, 0.0, -dN811, 0.0, 0.0, -dN911, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                    { 0.0, -dN111 ,0.0, 0.0 ,-dN211 ,0.0 ,0.0 ,-dN311, 0.0, 0.0, -dN411, 0.0, 0.0, -dN511 ,0.0, 0.0 , -dN611, 0.0 ,0.0 , -dN711, 0.0, 0.0, -dN811, 0.0,0.0, -dN911, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0  },
                    { 0.0, 0.0, -dN111 ,0.0, 0.0 ,-dN211 ,0.0 ,0.0 ,-dN311, 0.0, 0.0, -dN411, 0.0, 0.0, -dN511 ,0.0, 0.0 , -dN611, 0.0 ,0.0 , -dN711, 0.0, 0.0, -dN811, 0.0, 0.0, -dN911, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
                    };

                double[,] da12Matrix = new double[,]
                    {
                    { -dN112 ,0.0, 0.0 ,-dN212 ,0.0 ,0.0 ,-dN312, 0.0, 0.0, -dN412, 0.0, 0.0, -dN512 ,0.0, 0.0 , -dN612, 0.0 ,0.0 , -dN712, 0.0, 0.0, -dN812, 0.0, 0.0, -dN912, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                    { 0.0, -dN112 ,0.0, 0.0 ,-dN212 ,0.0 ,0.0 ,-dN312, 0.0, 0.0, -dN412, 0.0, 0.0, -dN512 ,0.0, 0.0 , -dN612, 0.0 ,0.0 , -dN712, 0.0, 0.0, -dN812, 0.0,0.0, -dN912, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0  },
                    { 0.0, 0.0, -dN112 ,0.0, 0.0 ,-dN212 ,0.0 ,0.0 ,-dN312, 0.0, 0.0, -dN412, 0.0, 0.0, -dN512 ,0.0, 0.0 , -dN612, 0.0 ,0.0 , -dN712, 0.0, 0.0, -dN812, 0.0, 0.0, -dN912, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
                    };

                double[,] da2Matrix = new double[,]
                    {
                    { -dN12 ,0.0, 0.0 ,-dN22 ,0.0 ,0.0 ,-dN32, 0.0, 0.0, -dN42, 0.0, 0.0, -dN52 ,0.0, 0.0 , -dN62, 0.0 ,0.0 , -dN72, 0.0, 0.0, -dN82, 0.0, 0.0, -dN92, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                    { 0.0, -dN12 ,0.0, 0.0 ,-dN22 ,0.0 ,0.0 ,-dN32, 0.0, 0.0, -dN42, 0.0, 0.0, -dN52 ,0.0, 0.0 , -dN62, 0.0 ,0.0 , -dN72, 0.0, 0.0, -dN82, 0.0,0.0, -dN92, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                    { 0.0, 0.0, -dN12 ,0.0, 0.0 ,-dN22 ,0.0 ,0.0 ,-dN32, 0.0, 0.0, -dN42, 0.0, 0.0, -dN52 ,0.0, 0.0 , -dN62, 0.0 ,0.0 , -dN72, 0.0, 0.0, -dN82, 0.0, 0.0, -dN92, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
                    };

                double[,] da22Matrix = new double[,]
                    {
                    { -dN122 ,0.0, 0.0 ,-dN222 ,0.0 ,0.0 ,-dN322, 0.0, 0.0, -dN422, 0.0, 0.0, -dN522 ,0.0, 0.0 , -dN622, 0.0 ,0.0 , -dN722, 0.0, 0.0, -dN822, 0.0, 0.0, -dN922,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                    { 0.0, -dN122 ,0.0, 0.0 ,-dN222 ,0.0 ,0.0 ,-dN322, 0.0, 0.0, -dN422, 0.0, 0.0, -dN522 ,0.0, 0.0 , -dN622, 0.0 ,0.0 , -dN722, 0.0, 0.0, -dN822, 0.0,0.0, -dN922, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                    { 0.0, 0.0, -dN122 ,0.0, 0.0 ,-dN222 ,0.0 ,0.0 ,-dN322, 0.0, 0.0, -dN422, 0.0, 0.0, -dN522 ,0.0, 0.0 , -dN622, 0.0 ,0.0 , -dN722, 0.0, 0.0, -dN822, 0.0, 0.0, -dN922, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
                    };

                double[,] da3Matrix = new double[,]
                    {
                    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, dN103, 0.0, 0.0, dN113, 0.0, 0.0, dN123, 0.0, 0.0, dN133, 0.0, 0.0, dN143, 0.0, 0.0, dN153, 0.0, 0.0, dN163, 0.0, 0.0, dN173, 0.0, 0.0, dN183, 0.0, 0.0 },
                    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, dN103, 0.0, 0.0, dN113, 0.0, 0.0, dN123, 0.0, 0.0, dN133, 0.0, 0.0, dN143, 0.0, 0.0, dN153, 0.0, 0.0, dN163, 0.0, 0.0, dN173, 0.0, 0.0, dN183, 0.0 },
                    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, dN103, 0.0, 0.0, dN113, 0.0, 0.0, dN123, 0.0, 0.0, dN133, 0.0, 0.0, dN143, 0.0, 0.0, dN153, 0.0, 0.0, dN163, 0.0, 0.0, dN173, 0.0, 0.0, dN183 }
                    };

                double[,] da33Matrix = new double[,]
                    {
                    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, dN1033, 0.0, 0.0, dN1133, 0.0, 0.0, dN1233, 0.0, 0.0, dN1333, 0.0, 0.0, dN1433, 0.0, 0.0, dN1533, 0.0, 0.0, dN1633, 0.0, 0.0, dN1733, 0.0, 0.0, dN1833, 0.0, 0.0 },
                    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, dN1033, 0.0, 0.0, dN1133, 0.0, 0.0, dN1233, 0.0, 0.0, dN1333, 0.0, 0.0, dN1433, 0.0, 0.0, dN1533, 0.0, 0.0, dN1633, 0.0, 0.0, dN1733, 0.0, 0.0, dN1833, 0.0 },
                    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, dN1033, 0.0, 0.0, dN1133, 0.0, 0.0, dN1233, 0.0, 0.0, dN1333, 0.0, 0.0, dN1433, 0.0, 0.0, dN1533, 0.0, 0.0, dN1633, 0.0, 0.0, dN1733, 0.0, 0.0, dN1833 }
                    };

                double[,] da4Matrix = new double[,]
                    {
                    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, dN104, 0.0, 0.0, dN114, 0.0, 0.0, dN124, 0.0, 0.0, dN134, 0.0, 0.0, dN144, 0.0, 0.0, dN154, 0.0, 0.0, dN164, 0.0, 0.0, dN174, 0.0, 0.0, dN184, 0.0, 0.0 },
                    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, dN104, 0.0, 0.0, dN114, 0.0, 0.0, dN124, 0.0, 0.0, dN134, 0.0, 0.0, dN144, 0.0, 0.0, dN154, 0.0, 0.0, dN164, 0.0, 0.0, dN174, 0.0, 0.0, dN184, 0.0 },
                    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, dN104, 0.0, 0.0, dN114, 0.0, 0.0, dN124, 0.0, 0.0, dN134, 0.0, 0.0, dN144, 0.0, 0.0, dN154, 0.0, 0.0, dN164, 0.0, 0.0, dN174, 0.0, 0.0, dN184 }
                    };

                double[,] da44Matrix = new double[,]
                    {
                    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, dN1044, 0.0, 0.0, dN1144, 0.0, 0.0, dN1244, 0.0, 0.0, dN1344, 0.0, 0.0, dN1444, 0.0, 0.0, dN1544, 0.0, 0.0, dN1644, 0.0, 0.0, dN1744, 0.0, 0.0, dN1844, 0.0, 0.0 },
                    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, dN1044, 0.0, 0.0, dN1144, 0.0, 0.0, dN1244, 0.0, 0.0, dN1344, 0.0, 0.0, dN1444, 0.0, 0.0, dN1544, 0.0, 0.0, dN1644, 0.0, 0.0, dN1744, 0.0, 0.0, dN1844, 0.0 },
                    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, dN1044, 0.0, 0.0, dN1144, 0.0, 0.0, dN1244, 0.0, 0.0, dN1344, 0.0, 0.0, dN1444, 0.0, 0.0, dN1544, 0.0, 0.0, dN1644, 0.0, 0.0, dN1744, 0.0, 0.0, dN1844 }
                    };

                var T1 = new Tuple<double[,], double[,], double[,], double[,], double[,]>(da1Matrix, da11Matrix, da2Matrix, da22Matrix, da12Matrix);
                var T2 = new Tuple<double[,], double[,], double[,], double[,]>(da3Matrix, da33Matrix, da4Matrix, da44Matrix);
                return new Tuple<double[,], Tuple<double[,], double[,], double[,], double[,], double[,]>, Tuple<double[,], double[,], double[,], double[,]>>(aMatrix, T1, T2);
            }
        }

        private Tuple<Tuple<double[], double[], double[], double[], double[]>, double[,], double, double[,], double[], double[,], double[,]> MasterSurfaceGeometry(double[,] da1Matrix, double[,] da2Matrix, double[,] da11Matrix,
            double[,] da12Matrix, double[,] da22Matrix)
        {
            double[] xupd = VectorOperations.VectorScalarProduct(NodalXUpdated(), -1);
            double[] surfaceVector1 = VectorOperations.MatrixVectorProduct(da1Matrix, xupd);
            double[] surfaceVector2 = VectorOperations.MatrixVectorProduct(da2Matrix, xupd);
            double[] surfaceVectorDerivative11 = VectorOperations.MatrixVectorProduct(da11Matrix, xupd);
            double[] surfaceVectorDerivative12 = VectorOperations.MatrixVectorProduct(da12Matrix, xupd);
            double[] surfaceVectorDerivative22 = VectorOperations.MatrixVectorProduct(da22Matrix, xupd);
            var surafaceVectors = new Tuple<double[], double[], double[], double[], double[]>(surfaceVector1, surfaceVector2,
                surfaceVectorDerivative11, surfaceVectorDerivative12, surfaceVectorDerivative22);
            double[,] m = new double[,]
                    {
                    { VectorOperations.VectorDotProduct(surfaceVector1, surfaceVector1), VectorOperations.VectorDotProduct(surfaceVector1, surfaceVector2) },
                    { VectorOperations.VectorDotProduct(surfaceVector2, surfaceVector1), VectorOperations.VectorDotProduct(surfaceVector2, surfaceVector2) }
                    };
            double detm = m[0, 0] * m[1, 1] - m[0, 1] * m[1, 0];
            double[,] mInv = new double[,]
                    {
                    { m[1,1]/detm, - m[0,1]/detm },
                    { -m[1,0]/detm, m[0,0]/detm }
                    };
            double[] normalUnitVec = VectorOperations.VectorScalarProductNew(VectorOperations.VectorCrossProduct(surfaceVector1, surfaceVector2), 1.0 / Math.Pow(detm, 0.5));
            double[,] curvatureTensor = new double[,]
                    {
                    { 0.0, 0.0 },
                    { 0.0, 0.0 }
                    };
            double[,] curvatureContravariantTensor = new double[,]
                    {
                    { 0.0, 0.0 },
                    { 0.0, 0.0 }
                    };
            if (Properties.MasterSegmentPolynomialDegree >= 2)
            {
                curvatureTensor[0, 0] = VectorOperations.VectorDotProduct(surfaceVectorDerivative11, normalUnitVec);
                curvatureTensor[0, 1] = VectorOperations.VectorDotProduct(surfaceVectorDerivative12, normalUnitVec);
                curvatureTensor[1, 0] = VectorOperations.VectorDotProduct(surfaceVectorDerivative12, normalUnitVec);
                curvatureTensor[1, 1] = VectorOperations.VectorDotProduct(surfaceVectorDerivative22, normalUnitVec);

                curvatureContravariantTensor[0, 0] = curvatureTensor[0, 0] * mInv[0, 0] * mInv[0, 0] +
                    curvatureTensor[0, 1] * mInv[0, 0] * mInv[1, 0] +
                    curvatureTensor[1, 0] * mInv[0, 1] * mInv[0, 0] +
                    curvatureTensor[1, 1] * mInv[0, 1] * mInv[1, 0];
                curvatureContravariantTensor[0, 1] = curvatureTensor[0, 0] * mInv[0, 0] * mInv[0, 1] +
                    curvatureTensor[0, 1] * mInv[0, 0] * mInv[1, 1] +
                    curvatureTensor[1, 0] * mInv[0, 1] * mInv[0, 1] +
                    curvatureTensor[1, 1] * mInv[0, 1] * mInv[1, 1];
                curvatureContravariantTensor[1, 0] = curvatureTensor[0, 0] * mInv[1, 0] * mInv[0, 0] +
                    curvatureTensor[0, 1] * mInv[1, 0] * mInv[1, 0] +
                    curvatureTensor[1, 0] * mInv[1, 1] * mInv[0, 0] +
                    curvatureTensor[1, 1] * mInv[1, 1] * mInv[1, 0];
                curvatureContravariantTensor[1, 1] = curvatureTensor[0, 0] * mInv[1, 0] * mInv[0, 1] +
                    curvatureTensor[0, 1] * mInv[1, 0] * mInv[1, 1] +
                    curvatureTensor[1, 0] * mInv[1, 1] * mInv[0, 1] +
                    curvatureTensor[1, 1] * mInv[1, 1] * mInv[1, 1];
            }
            return new Tuple<Tuple<double[], double[], double[], double[], double[]>, double[,], double, double[,], double[], double[,], double[,]>(surafaceVectors, m, detm, mInv, normalUnitVec, curvatureTensor, curvatureContravariantTensor);
        }
        private double SlaveSurfaceMetricDet(double[,] da3Matrix, double[,] da4Matrix)
        {
            double[] xupd = NodalXUpdated();
            double[] surfaceVector1 = VectorOperations.MatrixVectorProduct(da3Matrix, xupd);
            double[] surfaceVector2 = VectorOperations.MatrixVectorProduct(da4Matrix, xupd);
            double[,] m = new double[,]
                    {
                    { VectorOperations.VectorDotProduct(surfaceVector1, surfaceVector1), VectorOperations.VectorDotProduct(surfaceVector1, surfaceVector2) },
                    { VectorOperations.VectorDotProduct(surfaceVector2, surfaceVector1), VectorOperations.VectorDotProduct(surfaceVector2, surfaceVector2) }
                    };
            double detm = m[0, 0] * m[1, 1] - m[0, 1] * m[1, 0];
            return detm;
        }
        private double CalculatePenetration(double[,] aMatrix, double[] n)
        {
            double[,] AT = MatrixOperations.Transpose(aMatrix);
            double[] AT_n = VectorOperations.MatrixVectorProduct(AT, n);
            double[] xupd = NodalXUpdated();
            //double[] xupd =VectorOperations.VectorScalarProduct(NodalXUpdated(),-1d);
            double normalGap = VectorOperations.VectorDotProduct(xupd, AT_n);
            return normalGap;
        }
        private double[] CalculateDeltaKsi(double[] masterSlaveRelativeVector, double[,] aMatrix, double[,] metricTensor,
            double[] surfaceVector1, double[] surfaceVector2, double[] surfaceVectorDerivative11, double[] surfaceVectorDerivative12,
            double[] surfaceVectorDerivative22, double detm, double[] xupd)
        {
            double[] deltaKsi = new double[2];
            if (Properties.MasterSegmentPolynomialDegree == 1)
            {
                double e = VectorOperations.VectorDotProduct(surfaceVectorDerivative12, VectorOperations.MatrixVectorProduct(aMatrix, xupd));
                double[] f = new double[] { VectorOperations.VectorDotProduct(surfaceVector1, VectorOperations.MatrixVectorProduct(aMatrix, xupd)),
                VectorOperations.VectorDotProduct(surfaceVector2, VectorOperations.MatrixVectorProduct(aMatrix, xupd))};
                double[,] matrix = new double[,]
                {
                    { metricTensor[1,1], e - metricTensor[0,1] },
                    { e - metricTensor[1,0], metricTensor[0,0] }
                };
                double scalar = 1.0 / (detm - Math.Pow(e, 2) + 2 * e * metricTensor[0, 1]);
                deltaKsi = VectorOperations.VectorScalarProductNew(VectorOperations.MatrixVectorProduct(matrix, f), scalar);
            }
            else
            {
                double detDDF = (VectorOperations.VectorDotProduct(surfaceVector1, surfaceVector1) - VectorOperations.VectorDotProduct(surfaceVectorDerivative11, masterSlaveRelativeVector)) *
                    (VectorOperations.VectorDotProduct(surfaceVector2, surfaceVector2) - VectorOperations.VectorDotProduct(surfaceVectorDerivative22, masterSlaveRelativeVector)) -
                    (VectorOperations.VectorDotProduct(surfaceVector1, surfaceVector2) - VectorOperations.VectorDotProduct(surfaceVectorDerivative12, masterSlaveRelativeVector)) *
                    (VectorOperations.VectorDotProduct(surfaceVector2, surfaceVector1) - VectorOperations.VectorDotProduct(surfaceVectorDerivative12, masterSlaveRelativeVector));
                double scalar = -1.0 / detDDF;
                double[,] matrix = new double[,]
                {
                    { VectorOperations.VectorDotProduct(surfaceVector2, surfaceVector2) - VectorOperations.VectorDotProduct(surfaceVectorDerivative22, masterSlaveRelativeVector),
                      VectorOperations.VectorDotProduct(surfaceVectorDerivative12, masterSlaveRelativeVector) - VectorOperations.VectorDotProduct(surfaceVector1, surfaceVector2)},
                    { VectorOperations.VectorDotProduct(surfaceVectorDerivative12, masterSlaveRelativeVector) - VectorOperations.VectorDotProduct(surfaceVector2, surfaceVector1),
                     VectorOperations.VectorDotProduct(surfaceVector1, surfaceVector1) - VectorOperations.VectorDotProduct(surfaceVectorDerivative11, masterSlaveRelativeVector)}
                };
                double[] vector = new double[] { VectorOperations.VectorDotProduct(surfaceVector1, masterSlaveRelativeVector), VectorOperations.VectorDotProduct(surfaceVector2, masterSlaveRelativeVector) };
                vector = VectorOperations.VectorScalarProductNew(vector, -1.0);
                deltaKsi = VectorOperations.VectorScalarProductNew(VectorOperations.MatrixVectorProduct(matrix, vector), scalar);
            }
            return deltaKsi;
        }
        private double[] Project(double ksi1Initial, double ksi2Initial, double ksi3, double ksi4)
        {
            int maxIterations = 1000;
            double tol = Math.Pow(10.0, -6.0);
            double[] deltaKsi = new double[2];
            double[] ksi = new double[] { ksi1Initial, ksi2Initial };
            double[] xUpdated = NodalXUpdated();
            for (int i = 1; i <= maxIterations; i++)
            {
                Tuple<double[,], Tuple<double[,], double[,], double[,], double[,], double[,]>, Tuple<double[,], double[,], double[,], double[,]>> aMatrices = CalculatePositionMatrix(ksi[0], ksi[1], ksi3, ksi4);
                double[] masterSlaveRelativeVector = VectorOperations.MatrixVectorProduct(aMatrices.Item1, xUpdated);
                Tuple<Tuple<double[], double[], double[], double[], double[]>, double[,], double, double[,], double[], double[,], double[,]> masterSurfaceGeometry = MasterSurfaceGeometry(aMatrices.Item2.Item1, aMatrices.Item2.Item3, aMatrices.Item2.Item2, aMatrices.Item2.Item5, aMatrices.Item2.Item4);

                deltaKsi = CalculateDeltaKsi(masterSlaveRelativeVector, aMatrices.Item1, masterSurfaceGeometry.Item2,
                    masterSurfaceGeometry.Item1.Item1, masterSurfaceGeometry.Item1.Item2, masterSurfaceGeometry.Item1.Item3,
                    masterSurfaceGeometry.Item1.Item4, masterSurfaceGeometry.Item1.Item5, masterSurfaceGeometry.Item3,
                    xUpdated);
                ksi[0] += deltaKsi[0];
                ksi[1] += deltaKsi[1];
                if (VectorOperations.VectorNorm2(deltaKsi) <= tol)
                {
                    break;
                }
            }
            if (VectorOperations.VectorNorm2(deltaKsi) > tol)
            {
                throw new Exception("CPP not found in current iterations");
            }
            else
            {
                return ksi;
            }
        }
        private Tuple<double, double> GaussPoints(int i)
        {
            int iP = Properties.IntegrationPoints;
            double[] gaussPoints = new double[iP];
            double[] gaussWeights = new double[iP];
            if (iP == 1)
            {
                gaussPoints[0] = 0.0;
                gaussWeights[0] = 2.0;
            }
            else if (iP == 2)
            {
                gaussPoints[0] = -1.0 / Math.Sqrt(3);
                gaussPoints[1] = 1.0 / Math.Sqrt(3);
                gaussWeights[0] = 1.0;
                gaussWeights[1] = 1.0;
            }
            else if (iP == 3)
            {
                gaussPoints[0] = -0.77459;
                gaussPoints[1] = 0.0;
                gaussPoints[2] = 0.77459;
                gaussWeights[0] = 0.55555;
                gaussWeights[1] = 0.88888;
                gaussWeights[2] = 0.55555;
            }
            else if (iP == 4)
            {
                gaussPoints[0] = -0.86113;
                gaussPoints[1] = -0.33998;
                gaussPoints[2] = 0.33998;
                gaussPoints[3] = 0.86113;
                gaussWeights[0] = 0.34785;
                gaussWeights[1] = 0.65214;
                gaussWeights[2] = 0.65214;
                gaussWeights[3] = 0.34785;
            }
            else if (iP == 5)
            {
                gaussPoints[0] = -0.90617;
                gaussPoints[1] = -0.53846;
                gaussPoints[2] = 0.0;
                gaussPoints[3] = 0.53846;
                gaussPoints[4] = 0.90617;
                gaussWeights[0] = 0.23692;
                gaussWeights[1] = 0.47862;
                gaussWeights[2] = 0.56888;
                gaussWeights[3] = 0.47862;
                gaussWeights[4] = 0.23692;
            }
            else if (iP == 6)
            {
                gaussPoints[0] = -0.93246;
                gaussPoints[1] = -0.66120;
                gaussPoints[2] = -0.23861;
                gaussPoints[3] = 0.23861;
                gaussPoints[4] = 0.66120;
                gaussPoints[5] = 0.93246;
                gaussWeights[0] = 0.17132;
                gaussWeights[1] = 0.36076;
                gaussWeights[2] = 0.46791;
                gaussWeights[3] = 0.46791;
                gaussWeights[4] = 0.36076;
                gaussWeights[5] = 0.17132;
            }
            else if (iP == 7)
            {
                gaussPoints[0] = -0.94910;
                gaussPoints[1] = -0.74153;
                gaussPoints[2] = -0.40584;
                gaussPoints[3] = 0.0;
                gaussPoints[4] = 0.40584;
                gaussPoints[5] = 0.74153;
                gaussPoints[6] = 0.94910;
                gaussWeights[0] = 0.12948;
                gaussWeights[1] = 0.27970;
                gaussWeights[2] = 0.38183;
                gaussWeights[3] = 0.41795;
                gaussWeights[4] = 0.38183;
                gaussWeights[5] = 0.27970;
                gaussWeights[6] = 0.12948;
            }
            else if (iP == 8)
            {
                gaussPoints[0] = -0.96028;
                gaussPoints[1] = -0.79666;
                gaussPoints[2] = -0.52553;
                gaussPoints[3] = -0.18343;
                gaussPoints[4] = 0.18343;
                gaussPoints[5] = 0.52553;
                gaussPoints[6] = 0.79666;
                gaussPoints[7] = 0.96028;
                gaussWeights[0] = 0.10122;
                gaussWeights[1] = 0.22238;
                gaussWeights[2] = 0.31370;
                gaussWeights[3] = 0.36268;
                gaussWeights[4] = 0.36268;
                gaussWeights[5] = 0.31370;
                gaussWeights[6] = 0.22238;
                gaussWeights[7] = 0.10122;
            }
            else if (iP == 9)
            {
                gaussPoints[0] = -0.96816;
                gaussPoints[1] = -0.83603;
                gaussPoints[2] = -0.61337;
                gaussPoints[3] = -0.32425;
                gaussPoints[4] = 0.0;
                gaussPoints[5] = 0.32425;
                gaussPoints[6] = 0.61337;
                gaussPoints[7] = 0.83603;
                gaussPoints[8] = 0.96816;
                gaussWeights[0] = 0.08127;
                gaussWeights[1] = 0.18064;
                gaussWeights[2] = 0.26061;
                gaussWeights[3] = 0.31234;
                gaussWeights[4] = 0.33023;
                gaussWeights[5] = 0.31234;
                gaussWeights[6] = 0.26061;
                gaussWeights[7] = 0.18064;
                gaussWeights[8] = 0.08127;
            }
            else if (iP == 10)
            {
                gaussPoints[0] = -0.97390;
                gaussPoints[1] = -0.86506;
                gaussPoints[2] = -0.67940;
                gaussPoints[3] = -0.43339;
                gaussPoints[4] = -0.14887;
                gaussPoints[5] = 0.14887;
                gaussPoints[6] = 0.43339;
                gaussPoints[7] = 0.67940;
                gaussPoints[8] = 0.86506;
                gaussPoints[9] = 0.97390;
                gaussWeights[0] = 0.06667;
                gaussWeights[1] = 0.14945;
                gaussWeights[2] = 0.21908;
                gaussWeights[3] = 0.26926;
                gaussWeights[4] = 0.29552;
                gaussWeights[5] = 0.29552;
                gaussWeights[6] = 0.26926;
                gaussWeights[7] = 0.21908;
                gaussWeights[8] = 0.14945;
                gaussWeights[9] = 0.06667;
            }
            else if (iP == 11)
            {
                gaussPoints[0] = -0.97823;
                gaussPoints[1] = -0.88706;
                gaussPoints[2] = -0.73015;
                gaussPoints[3] = -0.51909;
                gaussPoints[4] = -0.26954;
                gaussPoints[5] = 0.0;
                gaussPoints[6] = 0.26954;
                gaussPoints[7] = 0.51909;
                gaussPoints[8] = 0.73015;
                gaussPoints[9] = 0.88706;
                gaussPoints[10] = 0.97823;
                gaussWeights[0] = 0.05567;
                gaussWeights[1] = 0.12558;
                gaussWeights[2] = 0.18629;
                gaussWeights[3] = 0.23319;
                gaussWeights[4] = 0.26280;
                gaussWeights[5] = 0.27292;
                gaussWeights[6] = 0.26280;
                gaussWeights[7] = 0.23319;
                gaussWeights[8] = 0.18629;
                gaussWeights[9] = 0.12558;
                gaussWeights[10] = 0.05567;
            }
            else if (iP == 12)
            {
                gaussPoints[0] = -0.98156;
                gaussPoints[1] = -0.90411;
                gaussPoints[2] = -0.76990;
                gaussPoints[3] = -0.58732;
                gaussPoints[4] = -0.36783;
                gaussPoints[5] = -0.12523;
                gaussPoints[6] = 0.12523;
                gaussPoints[7] = 0.36783;
                gaussPoints[8] = 0.58732;
                gaussPoints[9] = 0.76990;
                gaussPoints[10] = 0.90411;
                gaussPoints[11] = 0.98156;
                gaussWeights[0] = 0.047175;
                gaussWeights[1] = 0.106939;
                gaussWeights[2] = 0.16008;
                gaussWeights[3] = 0.20317;
                gaussWeights[4] = 0.23349;
                gaussWeights[5] = 0.249147;
                gaussWeights[6] = 0.249147;
                gaussWeights[7] = 0.23349;
                gaussWeights[8] = 0.20317;
                gaussWeights[9] = 0.16008;
                gaussWeights[10] = 0.106939;
                gaussWeights[11] = 0.047175;
            }
            else if (iP == 13)
            {
                gaussPoints[0] = -0.98418;
                gaussPoints[1] = -0.917598;
                gaussPoints[2] = -0.801578;
                gaussPoints[3] = -0.642349;
                gaussPoints[4] = -0.448493;
                gaussPoints[5] = -0.230458;
                gaussPoints[6] = 0.0;
                gaussPoints[7] = 0.230458;
                gaussPoints[8] = 0.448493;
                gaussPoints[9] = 0.642349;
                gaussPoints[10] = 0.801578;
                gaussPoints[11] = 0.917598;
                gaussPoints[12] = 0.98418;
                gaussWeights[0] = 0.040484;
                gaussWeights[1] = 0.0921215;
                gaussWeights[2] = 0.138874;
                gaussWeights[3] = 0.178146;
                gaussWeights[4] = 0.207816;
                gaussWeights[5] = 0.226283;
                gaussWeights[6] = 0.232552;
                gaussWeights[7] = 0.226283;
                gaussWeights[8] = 0.207816;
                gaussWeights[9] = 0.178146;
                gaussWeights[10] = 0.138874;
                gaussWeights[11] = 0.0921215;
                gaussWeights[12] = 0.040484;
            }
            else if (iP == 14)
            {
                gaussPoints[0] = -0.986284;
                gaussPoints[1] = -0.928435;
                gaussPoints[2] = -0.827201;
                gaussPoints[3] = -0.687293;
                gaussPoints[4] = -0.515249;
                gaussPoints[5] = -0.319112;
                gaussPoints[6] = -0.108055;
                gaussPoints[7] = 0.108055;
                gaussPoints[8] = 0.319112;
                gaussPoints[9] = 0.515249;
                gaussPoints[10] = 0.687293;
                gaussPoints[11] = 0.827201;
                gaussPoints[12] = 0.928435;
                gaussPoints[13] = 0.986284;
                gaussWeights[0] = 0.0351195;
                gaussWeights[1] = 0.0801581;
                gaussWeights[2] = 0.121519;
                gaussWeights[3] = 0.157203;
                gaussWeights[4] = 0.185538;
                gaussWeights[5] = 0.205198;
                gaussWeights[6] = 0.215264;
                gaussWeights[7] = 0.215264;
                gaussWeights[8] = 0.205198;
                gaussWeights[9] = 0.185538;
                gaussWeights[10] = 0.157203;
                gaussWeights[11] = 0.121519;
                gaussWeights[12] = 0.0801581;
                gaussWeights[13] = 0.0351195;
            }
            else
            {
                gaussPoints[0] = -0.987993;
                gaussPoints[1] = -0.937273;
                gaussPoints[2] = -0.848207;
                gaussPoints[3] = -0.724418;
                gaussPoints[4] = -0.570972;
                gaussPoints[5] = -0.394151;
                gaussPoints[6] = -0.201194;
                gaussPoints[7] = 0.0;
                gaussPoints[8] = 0.201194;
                gaussPoints[9] = 0.394151;
                gaussPoints[10] = 0.570972;
                gaussPoints[11] = 0.724418;
                gaussPoints[12] = 0.848207;
                gaussPoints[13] = 0.937273;
                gaussPoints[14] = 0.987993;
                gaussWeights[0] = 0.0307532;
                gaussWeights[1] = 0.070366;
                gaussWeights[2] = 0.107159;
                gaussWeights[3] = 0.139571;
                gaussWeights[4] = 0.166269;
                gaussWeights[5] = 0.186161;
                gaussWeights[6] = 0.198431;
                gaussWeights[7] = 0.202578;
                gaussWeights[8] = 0.198431;
                gaussWeights[9] = 0.186161;
                gaussWeights[10] = 0.166269;
                gaussWeights[11] = 0.139571;
                gaussWeights[12] = 0.107159;
                gaussWeights[13] = 0.070366;
                gaussWeights[14] = 0.0307532;
            }
            double GaussPoint = gaussPoints[i];
            double Weight = gaussWeights[i];
            return new Tuple<double, double>(GaussPoint, Weight);
        }
        private double[,] CalculateMainStiffnessPart(double[,] A, double[] n)
        {
            double[,] nxn = VectorOperations.VectorVectorTensorProduct(n, n);
            double[,] nxn_A = MatrixOperations.MatrixProduct(nxn, A);
            double[,] AT_nxn_A = MatrixOperations.MatrixProduct(MatrixOperations.Transpose(A), nxn_A);
            double[,] mainStiffnessMatrix = MatrixOperations.ScalarMatrixProductNew(PenaltyFactor, AT_nxn_A);
            return mainStiffnessMatrix;
        }

        private double[,] CalculateRotationalStiffnessPart(double[,] A, double[,] dA1, double[,] dA2, double[] n, double ksi3, double[,] mInv, double[] surfaceVector1, double[] surfaceVector2)
        {
            double[,] rotationalPart;
            double coef1 = PenaltyFactor * ksi3 * mInv[0, 0];
            double[,] rotationalPart1;
            double[,] n_x_surfaceVector1 = VectorOperations.VectorVectorTensorProduct(n, surfaceVector1);
            double[,] surfaceVector1_x_n = VectorOperations.VectorVectorTensorProduct(surfaceVector1, n);
            double[,] firstTerm = MatrixOperations.MatrixProduct(
                                                                    MatrixOperations.Transpose(dA1),
                                                                    MatrixOperations.MatrixProduct(n_x_surfaceVector1, A)
                                                                    );
            double[,] secondTerm = MatrixOperations.MatrixProduct(
                                                                    MatrixOperations.Transpose(A),
                                                                    MatrixOperations.MatrixProduct(surfaceVector1_x_n, dA1)
                                                                    );
            rotationalPart1 = MatrixOperations.ScalarMatrixProductNew(
                                                                        coef1,
                                                                        MatrixOperations.MatrixAddition(firstTerm, secondTerm)
                                                                        );
            double coef2 = PenaltyFactor * ksi3 * mInv[1, 0];
            double[,] rotationalPart2;
            double[,] n_x_surfaceVector2 = VectorOperations.VectorVectorTensorProduct(n, surfaceVector2);
            double[,] firstTerm2 = MatrixOperations.MatrixProduct(
                                                                    MatrixOperations.Transpose(dA1),
                                                                    MatrixOperations.MatrixProduct(n_x_surfaceVector2, A)
                                                                    );
            double[,] secondTerm2 = MatrixOperations.MatrixProduct(
                                                                    MatrixOperations.Transpose(A),
                                                                    MatrixOperations.MatrixProduct(surfaceVector1_x_n, dA2)
                                                                    );
            rotationalPart2 = MatrixOperations.ScalarMatrixProductNew(
                                                                        coef2,
                                                                        MatrixOperations.MatrixAddition(firstTerm2, secondTerm2)
                                                                        );
            double coef3 = PenaltyFactor * ksi3 * mInv[0, 1];
            double[,] rotationalPart3;
            double[,] surfaceVector2_x_n = VectorOperations.VectorVectorTensorProduct(surfaceVector2, n);
            double[,] firstTerm3 = MatrixOperations.MatrixProduct(
                                                                    MatrixOperations.Transpose(dA2),
                                                                    MatrixOperations.MatrixProduct(n_x_surfaceVector1, A)
                                                                    );
            double[,] secondTerm3 = MatrixOperations.MatrixProduct(
                                                                    MatrixOperations.Transpose(A),
                                                                    MatrixOperations.MatrixProduct(surfaceVector2_x_n, dA1)
                                                                    );
            rotationalPart3 = MatrixOperations.ScalarMatrixProductNew(
                                                                        coef3,
                                                                        MatrixOperations.MatrixAddition(firstTerm3, secondTerm3)
                                                                        );
            double coef4 = PenaltyFactor * ksi3 * mInv[1, 1];
            double[,] rotationalPart4;
            double[,] firstTerm4 = MatrixOperations.MatrixProduct(
                                                                    MatrixOperations.Transpose(dA2),
                                                                    MatrixOperations.MatrixProduct(n_x_surfaceVector2, A)
                                                                    );
            double[,] secondTerm4 = MatrixOperations.MatrixProduct(
                                                                    MatrixOperations.Transpose(A),
                                                                    MatrixOperations.MatrixProduct(surfaceVector2_x_n, dA2)
                                                                    );
            rotationalPart4 = MatrixOperations.ScalarMatrixProductNew(
                                                                        coef4,
                                                                        MatrixOperations.MatrixAddition(firstTerm4, secondTerm4)
                                                                        );
            rotationalPart = MatrixOperations.MatrixAddition(MatrixOperations.MatrixAddition(MatrixOperations.MatrixAddition(rotationalPart1, rotationalPart2), rotationalPart3), rotationalPart4);
            return rotationalPart;
        }
        private double[,] CalculateCurvatureStiffnessPart(double[,] A, double ksi3, double[,] h, double[] surfaceVector1, double[] surfaceVector2)
        {
            double[,] curvaturePart;
            double coef1 = PenaltyFactor * ksi3 * h[0, 0];
            double coef2 = PenaltyFactor * ksi3 * h[1, 0];
            double coef3 = PenaltyFactor * ksi3 * h[0, 1];
            double coef4 = PenaltyFactor * ksi3 * h[1, 1];

            double[,] surfaceVector1_x_surfaceVector1 = VectorOperations.VectorVectorTensorProduct(surfaceVector1, surfaceVector1);
            double[,] surfaceVector1_x_surfaceVector2 = VectorOperations.VectorVectorTensorProduct(surfaceVector1, surfaceVector2);
            double[,] surfaceVector2_x_surfaceVector1 = VectorOperations.VectorVectorTensorProduct(surfaceVector2, surfaceVector1);
            double[,] surfaceVector2_x_surfaceVector2 = VectorOperations.VectorVectorTensorProduct(surfaceVector2, surfaceVector2);

            double[,] Matrix1 = MatrixOperations.MatrixProduct(
                                                                    MatrixOperations.Transpose(A),
                                                                    MatrixOperations.MatrixProduct(surfaceVector1_x_surfaceVector1, A)
                                                                    );
            double[,] firstTerm = MatrixOperations.ScalarMatrixProductNew(
                                                                        coef1,
                                                                        Matrix1
                                                                        );
            double[,] Matrix2 = MatrixOperations.MatrixProduct(
                                                        MatrixOperations.Transpose(A),
                                                        MatrixOperations.MatrixProduct(surfaceVector2_x_surfaceVector1, A)
                                                        );
            double[,] secondTerm = MatrixOperations.ScalarMatrixProductNew(
                                                                        coef2,
                                                                        Matrix2
                                                                        );
            double[,] Matrix3 = MatrixOperations.MatrixProduct(
                                            MatrixOperations.Transpose(A),
                                            MatrixOperations.MatrixProduct(surfaceVector1_x_surfaceVector2, A)
                                            );
            double[,] thirdTerm = MatrixOperations.ScalarMatrixProductNew(
                                                                        coef3,
                                                                        Matrix3
                                                                        );
            double[,] Matrix4 = MatrixOperations.MatrixProduct(
                                MatrixOperations.Transpose(A),
                                MatrixOperations.MatrixProduct(surfaceVector2_x_surfaceVector2, A)
                                );
            double[,] fourthTerm = MatrixOperations.ScalarMatrixProductNew(
                                                                        coef4,
                                                                        Matrix4
                                                                        );
            curvaturePart = MatrixOperations.MatrixAddition(
                MatrixOperations.MatrixAddition(MatrixOperations.MatrixAddition(firstTerm, secondTerm), thirdTerm), fourthTerm);
            return curvaturePart;
        }
        public double[,] CreateGlobalStiffnessMatrix()
        {
            int nodesNumber = (Properties.MasterSegmentPolynomialDegree + 1) * (Properties.MasterSegmentPolynomialDegree + 1) +
                (Properties.SlaveSegmentPolynomialDegree + 1) * (Properties.SlaveSegmentPolynomialDegree + 1);
            double[,] globalStifnessMatrix = new double[3 * nodesNumber, 3 * nodesNumber];
            for (int i = 0; i < Properties.IntegrationPoints; i++)
            {
                double ihta1 = GaussPoints(i).Item1;
                double gW1 = GaussPoints(i).Item2;
                for (int j = 0; j < Properties.IntegrationPoints; j++)
                {
                    double ihta2 = GaussPoints(j).Item1;
                    double gW2 = GaussPoints(j).Item2;
                    double[] ksi = Project(0.0, 0.0, ihta1, ihta2);
                    if (Math.Abs(ksi[0]) <= 1.05 && Math.Abs(ksi[1]) <= 1.05)
                    {
                        Tuple<double[,], Tuple<double[,], double[,], double[,], double[,], double[,]>, Tuple<double[,], double[,], double[,], double[,]>> positionMatrices = CalculatePositionMatrix(ksi[0], ksi[1], ihta1, ihta2);
                        double[,] aMatrix = positionMatrices.Item1;
                        double[,] da1Matrix = positionMatrices.Item2.Item1;
                        double[,] da2Matrix = positionMatrices.Item2.Item3;
                        double[,] da11Matrix = positionMatrices.Item2.Item2;
                        double[,] da12Matrix = positionMatrices.Item2.Item5;
                        double[,] da22Matrix = positionMatrices.Item2.Item4;
                        Tuple<Tuple<double[], double[], double[], double[], double[]>, double[,], double, double[,], double[], double[,], double[,]> masterSurfaceCharacteristics = MasterSurfaceGeometry(da1Matrix, da2Matrix, da11Matrix,
                        da12Matrix, da22Matrix);
                        double[,] m = masterSurfaceCharacteristics.Item4;
                        double[] dRho1 = masterSurfaceCharacteristics.Item1.Item1;
                        double[] dRho2 = masterSurfaceCharacteristics.Item1.Item2;
                        double[] n = masterSurfaceCharacteristics.Item5;
                        double[,] h = masterSurfaceCharacteristics.Item7;
                        double ksi3 = CalculatePenetration(aMatrix, n);
                        if (ksi3 <= 0)
                        {
                            double slaveMetricTensorDet = SlaveSurfaceMetricDet(positionMatrices.Item3.Item1,
                                positionMatrices.Item3.Item3);
                            double[,] mainPart = CalculateMainStiffnessPart(aMatrix, n);
                            double[,] rotationalPart = CalculateRotationalStiffnessPart(aMatrix, da1Matrix, da2Matrix, n, ksi3,
                                m, dRho1, dRho2);
                            double[,] curvaturePart = new double[3 * nodesNumber, 3 * nodesNumber];
                            if (Properties.MasterSegmentPolynomialDegree != 1)
                            {
                                curvaturePart = CalculateCurvatureStiffnessPart(aMatrix, ksi3,
                                h, dRho1, dRho2);
                            }
                            double scalar = Math.Pow(slaveMetricTensorDet, 0.5) * gW1 * gW2;
                            double[,] StifnessMatrix = MatrixOperations.ScalarMatrixProductNew(scalar, MatrixOperations.MatrixAddition(MatrixOperations.MatrixAddition(mainPart, rotationalPart),
                                curvaturePart));
                            globalStifnessMatrix = MatrixOperations.MatrixAddition(globalStifnessMatrix, StifnessMatrix);
                        }
                    }
                }
            }
            return globalStifnessMatrix;
        }

        public double[] CreateInternalGlobalForcesVector()
        {
            int nodesNumber = (Properties.MasterSegmentPolynomialDegree + 1) * (Properties.MasterSegmentPolynomialDegree + 1) +
                (Properties.SlaveSegmentPolynomialDegree + 1) * (Properties.SlaveSegmentPolynomialDegree + 1);
            double[] internalGlobalForcesVector = new double[3 * nodesNumber];
            for (int i = 0; i < Properties.IntegrationPoints; i++)
            {
                double ihta1 = GaussPoints(i).Item1;
                double gW1 = GaussPoints(i).Item2;
                for (int j = 0; j < Properties.IntegrationPoints; j++)
                {
                    double ihta2 = GaussPoints(j).Item1;
                    double gW2 = GaussPoints(j).Item2;
                    double[] ksi = Project(0.0, 0.0, ihta1, ihta2);
                    if (Math.Abs(ksi[0]) <= 1.05 && Math.Abs(ksi[1]) <= 1.05)
                    {
                        Tuple<double[,], Tuple<double[,], double[,], double[,], double[,], double[,]>, Tuple<double[,], double[,], double[,], double[,]>> positionMatrices = CalculatePositionMatrix(ksi[0], ksi[1], ihta1, ihta2);
                        double[,] aMatrix = positionMatrices.Item1;
                        double[,] da1Matrix = positionMatrices.Item2.Item1;
                        double[,] da2Matrix = positionMatrices.Item2.Item3;
                        double[,] da11Matrix = positionMatrices.Item2.Item2;
                        double[,] da12Matrix = positionMatrices.Item2.Item5;
                        double[,] da22Matrix = positionMatrices.Item2.Item4;
                        Tuple<Tuple<double[], double[], double[], double[], double[]>, double[,], double, double[,], double[], double[,], double[,]> masterSurfaceCharacteristics = MasterSurfaceGeometry(da1Matrix, da2Matrix, da11Matrix,
                        da12Matrix, da22Matrix);
                        double[] n = masterSurfaceCharacteristics.Item5;
                        double ksi3 = CalculatePenetration(aMatrix, n);
                        if (ksi3 <= 0)
                        {
                            double slaveMetricTensor = SlaveSurfaceMetricDet(positionMatrices.Item3.Item1,
                                positionMatrices.Item3.Item3);
                            double scalar = Math.Pow(slaveMetricTensor, 0.5) * gW1 * gW2;
                            double[,] AT = MatrixOperations.Transpose(aMatrix);
                            double[] AT_n = VectorOperations.MatrixVectorProduct(AT, n);
                            double[] internalLocalForcesVector = VectorOperations.VectorScalarProductNew(VectorOperations.VectorScalarProductNew(AT_n, PenaltyFactor * ksi3), scalar);
                            internalGlobalForcesVector = VectorOperations.VectorVectorAddition(internalGlobalForcesVector, internalLocalForcesVector);
                        }
                    }
                }
            }
            return internalGlobalForcesVector;
        }

        public double[,] CreateMassMatrix()
        {
            int nodesNumber = (Properties.MasterSegmentPolynomialDegree + 1) * (Properties.MasterSegmentPolynomialDegree + 1) +
                (Properties.SlaveSegmentPolynomialDegree + 1) * (Properties.SlaveSegmentPolynomialDegree + 1);
            return new double[3 * nodesNumber, 3 * nodesNumber];
        }

        public double[,] CreateDampingMatrix()
        {
            int nodesNumber = (Properties.MasterSegmentPolynomialDegree + 1) * (Properties.MasterSegmentPolynomialDegree + 1) +
                (Properties.SlaveSegmentPolynomialDegree + 1) * (Properties.SlaveSegmentPolynomialDegree + 1);
            return new double[3 * nodesNumber, 3 * nodesNumber];
        }
    }
}

