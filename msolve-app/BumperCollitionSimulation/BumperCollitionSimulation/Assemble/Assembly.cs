using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BumperCollitionSimulation
{
    public class Assembly : IAssembly
    {
        public Dictionary<int, IElementProperties> ElementsProperties { get; set; }
        public Dictionary<int, INode> Nodes { get; set; }
        public Dictionary<int, Dictionary<int, int>> ElementsConnectivity { get; set; }
        public Dictionary<int, IElement> ElementsAssembly { get; set; }
        private int totalDOF;
        public Dictionary<int, bool[]> NodeFreedomAllocationList { get; set; }
        public bool ActivateBoundaryConditions { get; set; }
        public int[] BoundedDOFsVector { get; set; }
        public int[] ContactDOFsVector { get; set; }
        public int[] NoContactDOFsVector { get; set; }
        public bool ActivateParallelCalculations { get; set; }

        public Assembly()
        {
            ActivateParallelCalculations = false;
        }
        //private int[] boundedDOFsVector;

        //public int[] BoundedDOFsVector
        //{
        //    get
        //    {
        //        return boundedDOFsVector;
        //    }

        //    set
        //    {
        //        boundedDOFsVector = value;

        //    }
        //}

        private Dictionary<int, INode> AssignElementNodes(Dictionary<int, int> elementConnectivity)
        {
            Dictionary<int, INode> elementNodes = new Dictionary<int, INode>();
            for (int i = 1; i <= elementConnectivity.Count; i++)
            {
                int node = elementConnectivity[i];
                elementNodes[i] = Nodes[node];
            }
            return elementNodes;
        }

        private Dictionary<int, int> CreateNodeFreedomMapList()
        {
            Dictionary<int, int> nodeFMT = new Dictionary<int, int>();
            int baselineCounter = 0;
            for (int node = 1; node <= NodeFreedomAllocationList.Count; node++)
            {
                nodeFMT[node] = baselineCounter;
                bool[] nodeActiveDofs = NodeFreedomAllocationList[node];
                int nodeActiveDofsCount = nodeActiveDofs.Count(c => c == true);
                baselineCounter = baselineCounter + nodeActiveDofsCount;
            }
            totalDOF = baselineCounter;
            return nodeFMT;
        }

        private List<int> CreateElementFreedomList(Dictionary<int, int> singleElementConnectivity, Dictionary<int, int> nodeFMT,
            Dictionary<int, bool[]> elementFreedomSignature, string type, int degenerateElementsNodesCount)
        {
            List<int> globalDOFs = new List<int>();
            if (!(type== "IsoparamShell18"))
            {
                foreach (KeyValuePair<int, int> node in singleElementConnectivity)
                {
                    int localNode = node.Key;
                    int globalNode = node.Value;
                    int countActiveDOFs = elementFreedomSignature[localNode].Count(c => c == true);
                    for (int i = 0; i < countActiveDOFs; i++)
                    {
                        //
                        if (globalNode > degenerateElementsNodesCount / 2)
                        {
                            globalDOFs.Add(nodeFMT[globalNode - degenerateElementsNodesCount / 2] + i);//Needs fixing. Another method must be added to work for degenerate elements!!!

                        }
                        else
                        {
                            globalDOFs.Add(nodeFMT[globalNode] + i);//Needs fixing. Another method must be added to work for degenerate elements!!!

                        }
                    }

                }
            }
            else
            {
                for (int i =1; i<=9;i++)
                {
                    int localNode = i;
                    int globalNode = singleElementConnectivity[i];
                    int countActiveDOFs = elementFreedomSignature[localNode].Count(c => c == true);
                    for (int j = 0; j < countActiveDOFs; j++)
                    {
                        globalDOFs.Add(nodeFMT[globalNode] + j);
                    }

                }
            }
            return globalDOFs;
        }

        public void CreateElementsAssembly()
        {
            ElementsAssembly = new Dictionary<int, IElement>();
            Dictionary<int, int> nodefmt = CreateNodeFreedomMapList();
            var degenerateElementsNodesCount = DegenerateElementsCounts().Item2;
            for (int elem = 1; elem <= ElementsConnectivity.Count; elem++)
            {
                Dictionary<int, INode> elementNodes = AssignElementNodes(ElementsConnectivity[elem]);

                switch (ElementsProperties[elem].ElementType)
                { 
                    case "Hex8":
                        ElementsAssembly[elem] = new Hex8(ElementsProperties[elem], elementNodes);
                        break; 
                    case "ContactStS3D":
                        ElementsAssembly[elem] = new ContactStS3D(ElementsProperties[elem], elementNodes);
                        break;
                    case "IsoparamShell18":
                        ElementsAssembly[elem] = new IsoparamShell18(ElementsProperties[elem], elementNodes);
                        break;
                    case "ANSSolidShell8EAS":
                        ElementsAssembly[elem] = new ANSSolidShell8EAS(ElementsProperties[elem], elementNodes);
                        break;
                }
                Dictionary<int, bool[]> efs = ElementsAssembly[elem].ElementFreedomSignature;
                Dictionary<int, int> elemConnectivity = ElementsConnectivity[elem];
                List<int> eft = CreateElementFreedomList(elemConnectivity, nodefmt, efs, ElementsAssembly[elem].Properties.ElementType, degenerateElementsNodesCount);
                ElementsAssembly[elem].ElementFreedomList = eft;
            }
        }
        private Tuple<int, int> DegenerateElementsCounts()
        {
            int countElements = 0;
            int countNodes = 0;
            List<INode> degenElNodes = new List<INode>();
            for (int elem = 1; elem <= ElementsConnectivity.Count; elem++)
            {
                if (ElementsProperties[elem].ElementType == "IsoparamShell18")
                {
                    countElements += 1;
                    Dictionary<int, INode> elementNodes = AssignElementNodes(ElementsConnectivity[elem]);
                    for (int i = elementNodes.Keys.Min(); i <= elementNodes.Keys.Max(); i++)
                    {
                        if (!degenElNodes.Contains(elementNodes[i]))
                        {
                            degenElNodes.Add(elementNodes[i]);
                        }
                    }
                }
            }
                countNodes = degenElNodes.Count();
            return new Tuple<int, int>(countElements, countNodes);
        }

        public int CountElementsOfSameType(Type elementType)
        {
            int counter = 0;
            
            foreach (var item in ElementsAssembly)
            {
                Type kati = item.GetType();
                if (item.Value.GetType() == elementType)
                {
                    counter = counter + 1;
                }
            }
            return counter;
        }

        public void UpdateDisplacements(double[] totalDisplacementVector)
        {
            double[] fullTotalDisplacementVector = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(totalDisplacementVector, BoundedDOFsVector);
            for (int element = 1; element <= ElementsConnectivity.Count; element++)
            {
                int elementDofs = ElementsAssembly[element].ElementFreedomList.Count;
                double[] elementDisplacementVector = new double[elementDofs];
                for (int i = 0; i < elementDofs; i++)
                {
                    int localRow = i;
                    int globalRow = ElementsAssembly[element].ElementFreedomList[i];
                    elementDisplacementVector[localRow] = fullTotalDisplacementVector[globalRow];
                }
                ElementsAssembly[element].DisplacementVector = elementDisplacementVector;
            }
        }
        public void UpdateElementsIncrementalDisplacements(double[] deltaU)
        {
            double[] fullTotalDisplacementVector = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(deltaU, BoundedDOFsVector);
            for (int element = 1; element <= ElementsConnectivity.Count; element++)
            {
                    int elementDofs = ElementsAssembly[element].ElementFreedomList.Count;
                    double[] elementDisplacementVector = new double[elementDofs];
                    for (int i = 0; i < elementDofs; i++)
                    {
                        int localRow = i;
                        int globalRow = ElementsAssembly[element].ElementFreedomList[i];
                        elementDisplacementVector[localRow] = fullTotalDisplacementVector[globalRow];
                    }
                    if (ElementsAssembly[element].Properties.ElementType == "ContactStS3Df" ||
                    ElementsAssembly[element].Properties.ElementType == "Shell2DQuadratic4")
                    {
                        ElementsAssembly[element].UpdateIncrementalDisplacements(elementDisplacementVector);
                    }   
            }
        }
        public void MMCPCGUpdateDisplacements(double[] totalDisplacementVector)
        {
            double[] fullTotalDisplacementVector = BoundaryConditionsImposition.MMCPCGCreateFullVectorFromReducedVector(totalDisplacementVector, BoundedDOFsVector, ContactDOFsVector, NoContactDOFsVector);
            for (int element = 1; element <= ElementsConnectivity.Count; element++)
            {
                int elementDofs = ElementsAssembly[element].ElementFreedomList.Count;
                double[] elementDisplacementVector = new double[elementDofs];
                for (int i = 0; i < elementDofs; i++)
                {
                    int localRow = i;
                    int globalRow = ElementsAssembly[element].ElementFreedomList[i];
                    elementDisplacementVector[localRow] = fullTotalDisplacementVector[globalRow];
                }
                ElementsAssembly[element].DisplacementVector = elementDisplacementVector;
            }
        }
        public void MMCPCGUpdateElementsIncrementalDisplacements(double[] deltaU)
        {
            double[] fullTotalDisplacementVector = BoundaryConditionsImposition.MMCPCGCreateFullVectorFromReducedVector(deltaU, BoundedDOFsVector, ContactDOFsVector, NoContactDOFsVector);
            for (int element = 1; element <= ElementsConnectivity.Count; element++)
            {
                int elementDofs = ElementsAssembly[element].ElementFreedomList.Count;
                double[] elementDisplacementVector = new double[elementDofs];
                for (int i = 0; i < elementDofs; i++)
                {
                    int localRow = i;
                    int globalRow = ElementsAssembly[element].ElementFreedomList[i];
                    elementDisplacementVector[localRow] = fullTotalDisplacementVector[globalRow];
                }
                if (ElementsAssembly[element].Properties.ElementType == "ContactStS3Df")
                {
                    ElementsAssembly[element].UpdateIncrementalDisplacements(elementDisplacementVector);
                }
            }
        }
        public double[] MMCPCGRearrange(double[] reducedVector)
        {
            double[] RearrangedVector = BoundaryConditionsImposition.MMCPCGRearrangeReducedVector(reducedVector, BoundedDOFsVector, ContactDOFsVector, NoContactDOFsVector);
            return RearrangedVector;
        }
        public void InitializeContactTangentialProperties()
        {
            for (int element = 1; element <= ElementsConnectivity.Count; element++)
            {
                if (ElementsAssembly[element].Properties.ElementType == "ContactStS2Df" || 
                    ElementsAssembly[element].Properties.ElementType == "ContactStS3Df" ||
                    ElementsAssembly[element].Properties.ElementType == "ContactNtS3Df" ||
                    ElementsAssembly[element].Properties.ElementType == "Shell2DQuadratic4")
                {
                    ElementsAssembly[element].InitializeTangentialProperties();

                }
            }
        }
        public void UpdateContactTangentialProperties()
        {
            for (int element = 1; element <= ElementsConnectivity.Count; element++)
            {
                if 
                    (
                    ElementsAssembly[element].Properties.ElementType == "ContactStS2Df" ||
                    ElementsAssembly[element].Properties.ElementType == "ContactStS3Df" ||
                    ElementsAssembly[element].Properties.ElementType == "ContactNtS3Df"
                    )
                {
                    ElementsAssembly[element].UpdateTangentialProperties();

                }
            }
        }
        public void InitializeContactSurfaceVectors()
        {
            for (int element = 1; element <= ElementsConnectivity.Count; element++)
            {
                if (ElementsAssembly[element].Properties.ElementType == "ContactStS3Df" ||
                    ElementsAssembly[element].Properties.ElementType == "ContactNtS3Df")
                {
                    ElementsAssembly[element].InitializeContactSurfaceGeometry();

                }
            }
        }
        public void UpdateContactSurfaceVectors()
        {
            for (int element = 1; element <= ElementsConnectivity.Count; element++)
            {
                if (ElementsAssembly[element].Properties.ElementType == "ContactStS3Df" ||
                    ElementsAssembly[element].Properties.ElementType == "ContactNtS3Df")
                {
                    ElementsAssembly[element].UpdateContactSurfaceGeometry();

                }
            }
        }
        public void UpdateAccelerations(double[] totalAccelerationsVector)
        {
            double[] fullTotalAccelerationsVector = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(totalAccelerationsVector, BoundedDOFsVector);
            for (int element = 1; element <= ElementsConnectivity.Count; element++)
            {
                int elementDofs = ElementsAssembly[element].ElementFreedomList.Count;
                double[] elementAccelerationsVector = new double[elementDofs];
                for (int i = 0; i < elementDofs; i++)
                {
                    int localRow = i;
                    int globalRow = ElementsAssembly[element].ElementFreedomList[i];
                    elementAccelerationsVector[localRow] = fullTotalAccelerationsVector[globalRow];
                }
                ElementsAssembly[element].AccelerationVector = elementAccelerationsVector;
            }
        }

        public double[,] CreateTotalStiffnessMatrix()
        {
            double[,] totalStiffnessMatrix;
            if (ActivateParallelCalculations)
            {
                totalStiffnessMatrix = CalculateTotalStiffnessMatrixParallel();
                return totalStiffnessMatrix;
            }
            else
            {
                totalStiffnessMatrix = CalculateTotalStiffnessMatrix();
                return totalStiffnessMatrix;
            }
        }
        private double[,] CalculateTotalStiffnessMatrix()
        {
            double[,] totalStiffnessMatrix = new double[totalDOF, totalDOF];
            //List<int> falseStiffness = new List<int>();
            for (int element = 1; element <= ElementsConnectivity.Count; element++)
            {
                int elementDofs = ElementsAssembly[element].ElementFreedomList.Count;
                double[,] elementStiffnessMatrix = ElementsAssembly[element].CreateGlobalStiffnessMatrix();
                //if (!MatrixOperations.CheckIfSymmetric(elementStiffnessMatrix))
                //{
                //    falseStiffness.Add(element);
                //}
                for (int i = 0; i < elementDofs; i++)
                {
                    int localRow = i;
                    int globalRow = ElementsAssembly[element].ElementFreedomList[i];
                    for (int j = 0; j < elementDofs; j++)
                    {
                        int localColumn = j;
                        int globalColumn = ElementsAssembly[element].ElementFreedomList[j];
                        totalStiffnessMatrix[globalRow, globalColumn] = totalStiffnessMatrix[globalRow, globalColumn] + elementStiffnessMatrix[localRow, localColumn];
                    }
                }
            }
            //var falseStiff = falseStiffness.ToArray();
            //VectorOperations.PrintIntVectorToFile(falseStiff, @"C:\Users\Public\Documents\" + "nonSymmetricStiffnessElements.dat");

            if (ActivateBoundaryConditions)
            {
                double[,] reducedStiffnessMatrix = BoundaryConditionsImposition.ReducedTotalStiff(totalStiffnessMatrix, BoundedDOFsVector);
                return reducedStiffnessMatrix;
            }
            else
            {
                return totalStiffnessMatrix;
            }
        }

        private double[,] CalculateTotalStiffnessMatrixParallel()
        {
            double[,] totalStiffnessMatrix = new double[totalDOF, totalDOF];
            Parallel.For(1, ElementsConnectivity.Count + 1, element =>
           {
               int elementDofs = ElementsAssembly[element].ElementFreedomList.Count;
               double[,] elementStiffnessMatrix = ElementsAssembly[element].CreateGlobalStiffnessMatrix();
                for (int i = 0; i < elementDofs; i++)
               {
                   int localRow = i;
                   int globalRow = ElementsAssembly[element].ElementFreedomList[i];
                   for (int j = 0; j < elementDofs; j++)
                   {
                       int localColumn = j;
                       int globalColumn = ElementsAssembly[element].ElementFreedomList[j];
                       lock(totalStiffnessMatrix)
                       {
                           totalStiffnessMatrix[globalRow, globalColumn] = totalStiffnessMatrix[globalRow, globalColumn] + elementStiffnessMatrix[localRow, localColumn];
                       }
                   }
               }
           });
            if (ActivateBoundaryConditions)
            {
                double[,] reducedStiffnessMatrix = BoundaryConditionsImposition.ReducedTotalStiff(totalStiffnessMatrix, BoundedDOFsVector);
                return reducedStiffnessMatrix;
            }
            else
            {
                return totalStiffnessMatrix;
            }
        }
        public double[,] CreateStiffnessMatrixLinearPart()
        {
            double[,] totalStiffnessMatrix = new double[totalDOF, totalDOF];
            for (int element = 1; element <= ElementsConnectivity.Count; element++)
            {
                int elementDofs = ElementsAssembly[element].ElementFreedomList.Count;
                double[,] elementStiffnessMatrix = ElementsAssembly[element].CreateGlobalStiffnessMatrix();

                for (int i = 0; i < elementDofs; i++)
                {
                    int localRow = i;
                    int globalRow = ElementsAssembly[element].ElementFreedomList[i];
                    for (int j = 0; j < elementDofs; j++)
                    {
                        int localColumn = j;
                        int globalColumn = ElementsAssembly[element].ElementFreedomList[j];
                        totalStiffnessMatrix[globalRow, globalColumn] = totalStiffnessMatrix[globalRow, globalColumn] + elementStiffnessMatrix[localRow, localColumn];
                    }
                }
            }
            return totalStiffnessMatrix;
        }

        public double[,] UpdateStifnessMatrix(double[,] stiffnessMatrixLinearPart)
        {
            double[,] totalStiffnessMatrix = new double[totalDOF, totalDOF];
            for (int element = 1; element <= ElementsConnectivity.Count; element++)
            {
                if (ElementsAssembly[element].Properties.ElementType == "ContactNtN2D" ||
                    ElementsAssembly[element].Properties.ElementType == "ContactNtN2Df" ||
                    ElementsAssembly[element].Properties.ElementType == "ContactNtS2D" ||
                    ElementsAssembly[element].Properties.ElementType == "ContactNtS2Df" ||
                    ElementsAssembly[element].Properties.ElementType == "ContactNtS3D" ||
                    ElementsAssembly[element].Properties.ElementType == "ContactNtS3Df" ||
                    ElementsAssembly[element].Properties.ElementType == "ContactStS2D" ||
                    ElementsAssembly[element].Properties.ElementType == "ContactStS2Df" ||
                    ElementsAssembly[element].Properties.ElementType == "ContactStS3D" ||
                    ElementsAssembly[element].Properties.ElementType == "ContactStS3Df")
                {
                    int elementDofs = ElementsAssembly[element].ElementFreedomList.Count;
                    double[,] elementStiffnessMatrix = ElementsAssembly[element].CreateGlobalStiffnessMatrix();

                    for (int i = 0; i < elementDofs; i++)
                    {
                        int localRow = i;
                        int globalRow = ElementsAssembly[element].ElementFreedomList[i];
                        for (int j = 0; j < elementDofs; j++)
                        {
                            int localColumn = j;
                            int globalColumn = ElementsAssembly[element].ElementFreedomList[j];
                            totalStiffnessMatrix[globalRow, globalColumn] = totalStiffnessMatrix[globalRow, globalColumn] + elementStiffnessMatrix[localRow, localColumn];
                        }
                    }
                }
            }
            totalStiffnessMatrix = MatrixOperations.MatrixAddition(totalStiffnessMatrix, stiffnessMatrixLinearPart);
            if (ActivateBoundaryConditions)
            {
                double[,] reducedStiffnessMatrix = BoundaryConditionsImposition.ReducedTotalStiff(totalStiffnessMatrix, BoundedDOFsVector);
                return reducedStiffnessMatrix;
            }
            else
            {
                return totalStiffnessMatrix;
            }
        }

        public Tuple<double[,], double[,], double[,], double[,]> MMCPCGCreateTotalStiffnessMatrix()
        {
            int[] contactDof = ContactDOFsVector;
            int[] noContactDof = NoContactDOFsVector;
            double[,] totalStiffnessMatrix = new double[totalDOF, totalDOF];
            for (int element = 1; element <= ElementsConnectivity.Count; element++)
            {
                int elementDofs = ElementsAssembly[element].ElementFreedomList.Count;
                double[,] elementStiffnessMatrix = ElementsAssembly[element].CreateGlobalStiffnessMatrix();

                for (int i = 0; i < elementDofs; i++)
                {
                    int localRow = i;
                    int globalRow = ElementsAssembly[element].ElementFreedomList[i];
                    for (int j = 0; j < elementDofs; j++)
                    {
                        int localColumn = j;
                        int globalColumn = ElementsAssembly[element].ElementFreedomList[j];
                        totalStiffnessMatrix[globalRow, globalColumn] = totalStiffnessMatrix[globalRow, globalColumn] + elementStiffnessMatrix[localRow, localColumn];
                    }
                }
            }

            Tuple<double[,], double[,], double[,], double[,]> reducedStiffnessMatrices = BoundaryConditionsImposition.MMCPCGReducedTotalStiffMatrices(totalStiffnessMatrix, BoundedDOFsVector,
                contactDof, noContactDof);
            return reducedStiffnessMatrices;
        }

        public void SeperateContactDoF()
        {
            List<int> contactDofList = new List<int>();
            List<int> noContactDofList = new List<int>();
            for (int element = 1; element <= ElementsConnectivity.Count; element++)
            {
                if(ElementsAssembly[element].Properties.ElementType == "ContactNtN2D" ||
                    ElementsAssembly[element].Properties.ElementType == "ContactNtN2Df" ||
                    ElementsAssembly[element].Properties.ElementType == "ContactNtS2D" ||
                    ElementsAssembly[element].Properties.ElementType == "ContactNtS2Df" ||
                    ElementsAssembly[element].Properties.ElementType == "ContactNtS3D" ||
                    ElementsAssembly[element].Properties.ElementType == "ContactNtS3Df" ||
                    ElementsAssembly[element].Properties.ElementType == "ContactStS2D" ||
                    ElementsAssembly[element].Properties.ElementType == "ContactStS2Df" ||
                    ElementsAssembly[element].Properties.ElementType == "ContactStS3D" ||
                    ElementsAssembly[element].Properties.ElementType == "ContactStS3Df")
                {
                    var l = ElementsAssembly[element].ElementFreedomList;
                    foreach(var elemDof in l)
                    {
                        if (!contactDofList.Contains(elemDof + 1))
                        {
                            contactDofList.Add(elemDof + 1);
                        }
                    }
                }
                //else
                //{
                //    var l = ElementsAssembly[element].ElementFreedomList;
                //    foreach (var elemDof in l)
                //    {
                //        if (!noContactDofList.Contains(elemDof))
                //        {
                //            noContactDofList.Add(elemDof);
                //        }
                //    }
                //}
            }
            for (int element = 1; element <= ElementsConnectivity.Count; element++)
            {
                var l = ElementsAssembly[element].ElementFreedomList;
                foreach (var elemDof in l)
                {
                    if (!contactDofList.Contains(elemDof + 1) && !noContactDofList.Contains(elemDof + 1))
                    {
                        noContactDofList.Add(elemDof + 1);
                    }
                }
            }
            contactDofList = contactDofList.OrderBy(a => a).ToList();
            noContactDofList = noContactDofList.OrderBy(a => a).ToList();
            int[] contactDof = contactDofList.ToArray<int>();
            int[] noContactDof = noContactDofList.ToArray<int>();
            ContactDOFsVector = contactDof;
            NoContactDOFsVector = noContactDof;
        }

        public double[,] CreateTotalMassMatrix()
        {
            if (ActivateParallelCalculations)
            {
                var totalMassMatrix = CalculateTotalMassMatrixParallel();
                return totalMassMatrix;
            }
            else
            {
                var totalMassMatrix = CalculateTotalMassMatrix();
                return totalMassMatrix;
            }
        }
        private double[,] CalculateTotalMassMatrix()
        {
            double[,] totalMassMatrix = new double[totalDOF, totalDOF];
            for (int element = 1; element <= ElementsConnectivity.Count; element++)
            {
                int elementDofs = ElementsAssembly[element].ElementFreedomList.Count;
                double[,] elementMassMatrix = ElementsAssembly[element].CreateMassMatrix();
                for (int i = 0; i < elementDofs; i++)
                {
                    int localRow = i;
                    int globalRow = ElementsAssembly[element].ElementFreedomList[i];
                    for (int j = 0; j < elementDofs; j++)
                    {
                        int localColumn = j;
                        int globalColumn = ElementsAssembly[element].ElementFreedomList[j];
                        totalMassMatrix[globalRow, globalColumn] = totalMassMatrix[globalRow, globalColumn] + elementMassMatrix[localRow, localColumn];
                    }
                }
            }
            if (ActivateBoundaryConditions)
            {
                double[,] reducedMassMatrix = BoundaryConditionsImposition.ReducedTotalStiff(totalMassMatrix, BoundedDOFsVector);
                return reducedMassMatrix;
            }
            else
            {
                return totalMassMatrix;
            }
        }
        private double[,] CalculateTotalMassMatrixParallel()
        {
            double[,] totalMassMatrix = new double[totalDOF, totalDOF];
            Parallel.For(1, ElementsConnectivity.Count + 1, element =>
            {
                int elementDofs = ElementsAssembly[element].ElementFreedomList.Count;
                double[,] elementMassMatrix = ElementsAssembly[element].CreateMassMatrix();
                for (int i = 0; i < elementDofs; i++)
                {
                    int localRow = i;
                    int globalRow = ElementsAssembly[element].ElementFreedomList[i];
                    for (int j = 0; j < elementDofs; j++)
                    {
                        int localColumn = j;
                        int globalColumn = ElementsAssembly[element].ElementFreedomList[j];
                        lock (totalMassMatrix)
                        {
                            totalMassMatrix[globalRow, globalColumn] = totalMassMatrix[globalRow, globalColumn] + elementMassMatrix[localRow, localColumn];
                        }
                    }
                }
            });
            if (ActivateBoundaryConditions)
            {
                double[,] reducedMassMatrix = BoundaryConditionsImposition.ReducedTotalStiff(totalMassMatrix, BoundedDOFsVector);
                return reducedMassMatrix;
            }
            else
            {
                return totalMassMatrix;
            }

        }
        public double[,] CreateTotalDampingMatrix()
        {
            if (ActivateParallelCalculations)
            {
                var totalDampingMatrix = CalculateTotalDampingMatrixParallel();
                return totalDampingMatrix;
            }
            else
            {
                var totalDampingMatrix = CalculateTotalDampingMatrix();
                return totalDampingMatrix;
            }
        }
        private double[,] CalculateTotalDampingMatrix()
        {
            double[,] totalDampingMatrix = new double[totalDOF, totalDOF];
            for (int element = 1; element <= ElementsConnectivity.Count; element++)
            {
                int elementDofs = ElementsAssembly[element].ElementFreedomList.Count;
                double[,] elementDampingMatrix = ElementsAssembly[element].CreateDampingMatrix();

                for (int i = 0; i < elementDofs; i++)
                {
                    int localRow = i;
                    int globalRow = ElementsAssembly[element].ElementFreedomList[i];
                    for (int j = 0; j < elementDofs; j++)
                    {
                        int localColumn = j;
                        int globalColumn = ElementsAssembly[element].ElementFreedomList[j];
                        totalDampingMatrix[globalRow, globalColumn] = totalDampingMatrix[globalRow, globalColumn] + elementDampingMatrix[localRow, localColumn];
                    }
                }
            }

            if (ActivateBoundaryConditions)
            {
                double[,] reducedDampingMatrix = BoundaryConditionsImposition.ReducedTotalStiff(totalDampingMatrix, BoundedDOFsVector);
                return reducedDampingMatrix;
            }
            else
            {
                return totalDampingMatrix;
            }

        }
        private double[,] CalculateTotalDampingMatrixParallel()
        {
            double[,] totalDampingMatrix = new double[totalDOF, totalDOF];
            Parallel.For(1, ElementsConnectivity.Count + 1, element =>
            {
                int elementDofs = ElementsAssembly[element].ElementFreedomList.Count;
                double[,] elementDampingMatrix = ElementsAssembly[element].CreateDampingMatrix();

                for (int i = 0; i < elementDofs; i++)
                {
                    int localRow = i;
                    int globalRow = ElementsAssembly[element].ElementFreedomList[i];
                    for (int j = 0; j < elementDofs; j++)
                    {
                        int localColumn = j;
                        int globalColumn = ElementsAssembly[element].ElementFreedomList[j];
                        lock (totalDampingMatrix)
                        {
                            totalDampingMatrix[globalRow, globalColumn] = totalDampingMatrix[globalRow, globalColumn] + elementDampingMatrix[localRow, localColumn];
                        }
                    }
                }
            });

            if (ActivateBoundaryConditions)
            {
                double[,] reducedDampingMatrix = BoundaryConditionsImposition.ReducedTotalStiff(totalDampingMatrix, BoundedDOFsVector);
                return reducedDampingMatrix;
            }
            else
            {
                return totalDampingMatrix;
            }
        }
        public double[] CreateTotalInternalForcesVector()
        {
            if (ActivateParallelCalculations)
            {
                var internalForcesTotalVector = CalculateTotalInternalForcesVectorParallel();
                return internalForcesTotalVector;
            }
            else
            {
                var internalForcesTotalVector = CalculateTotalInternalForcesVector();
                return internalForcesTotalVector;
            }
        }
        private double[] CalculateTotalInternalForcesVector()
        {
            double[] internalForcesTotalVector = new double[totalDOF];
            //List<int> Elements = new List<int>();
            for (int element = 1; element <= ElementsConnectivity.Count; element++)
            {
                int elementDofs = ElementsAssembly[element].ElementFreedomList.Count;
                double[] elementInternalGlobalForcesVector = ElementsAssembly[element].CreateInternalGlobalForcesVector();
                //if (!VectorOperations.CheckForNonZeroElements(elementInternalGlobalForcesVector))
                //{
                //    Elements.Add(element);
                //}
                for (int i = 0; i < elementDofs; i++)
                {
                    int localLine = i;
                    int globalLine = ElementsAssembly[element].ElementFreedomList[i];
                    internalForcesTotalVector[globalLine] = internalForcesTotalVector[globalLine] + elementInternalGlobalForcesVector[localLine];
                }
            }
            //var falseForceElmnts = Elements.ToArray();
            //VectorOperations.PrintIntVectorToFile(falseForceElmnts, @"C:\Users\Public\Documents\" + "ElementsMass.dat");
            if (ActivateBoundaryConditions)
            {
                double[] reducedInternalForcesVector = BoundaryConditionsImposition.ReducedVector(internalForcesTotalVector, BoundedDOFsVector);
                return reducedInternalForcesVector;
            }
            return internalForcesTotalVector;
        }
        private double[] CalculateTotalInternalForcesVectorParallel()
        {
            double[] internalForcesTotalVector = new double[totalDOF];
            Parallel.For(1, ElementsConnectivity.Count + 1, element =>
            {
                int elementDofs = ElementsAssembly[element].ElementFreedomList.Count;
                double[] elementInternalGlobalForcesVector = ElementsAssembly[element].CreateInternalGlobalForcesVector();
                for (int i = 0; i < elementDofs; i++)
                {
                    int localLine = i;
                    int globalLine = ElementsAssembly[element].ElementFreedomList[i];
                    lock (internalForcesTotalVector)
                    {
                        internalForcesTotalVector[globalLine] = internalForcesTotalVector[globalLine] + elementInternalGlobalForcesVector[localLine];
                    }
                }
            });
            if (ActivateBoundaryConditions)
            {
                double[] reducedInternalForcesVector = BoundaryConditionsImposition.ReducedVector(internalForcesTotalVector, BoundedDOFsVector);
                return reducedInternalForcesVector;
            }
            return internalForcesTotalVector;
        }
        public double[] MMCPCGCreateTotalInternalForcesVectors()
        {
            int[] contactDof = ContactDOFsVector;
            int[] noContactDof = NoContactDOFsVector;
            double[] internalForcesTotalVector = new double[totalDOF];
            for (int element = 1; element <= ElementsConnectivity.Count; element++)
            {
                int elementDofs = ElementsAssembly[element].ElementFreedomList.Count;
                double[] elementInternalGlobalForcesVector = ElementsAssembly[element].CreateInternalGlobalForcesVector();
                for (int i = 0; i < elementDofs; i++)
                {
                    int localLine = i;
                    int globalLine = ElementsAssembly[element].ElementFreedomList[i];
                    internalForcesTotalVector[globalLine] = internalForcesTotalVector[globalLine] + elementInternalGlobalForcesVector[localLine];
                }
            }
            double[] reducedInternalForcesVector = BoundaryConditionsImposition.MMCPCGReducedVector(internalForcesTotalVector, BoundedDOFsVector,
                contactDof, noContactDof).Item3;
            return reducedInternalForcesVector;
        }
        public Tuple<double[], double[]> MMCPCGSeperateReducedExternalForcesVectors(double[] extForces)
        {
            int[] contactDof = ContactDOFsVector;
            int[] noContactDof = NoContactDOFsVector;
            var reducedExternalForces = BoundaryConditionsImposition.MMCPCGSeperateReducedVectorMatrices(extForces, BoundedDOFsVector,
            contactDof, noContactDof);
            return new Tuple<double[], double[]>(reducedExternalForces.Item1, reducedExternalForces.Item2);
        }

        public static Dictionary<int, INode> CalculateFinalNodalCoordinates(Dictionary<int, INode> nodesList, double[] diplacements)
        {
            if (nodesList.Count*2 != diplacements.Length)
            {
                throw new Exception("Nodes list does not match with displacements vector! Currently works only for 2 degreed of freedom per node");
            }
            Dictionary<int, INode> finalNodesList = new Dictionary<int, INode>();
            for (int i = 1; i <= nodesList.Count; i++)
            {
                double x = nodesList[i].XCoordinate + diplacements[2 * i - 2];
                double y = nodesList[i].YCoordinate + diplacements[2 * i - 1];
                INode finalNode = new Node(x, y);
                finalNodesList.Add(i, finalNode); 
            }
            return finalNodesList;
        }

        public Dictionary<int, double[]> GetElementsInternalForces(double[] totalInternalForcesVector)
        {
            UpdateDisplacements(totalInternalForcesVector);
            Dictionary<int, double[]> elementsInternalForces = new Dictionary<int, double[]>();
            for (int element = 1; element <= ElementsConnectivity.Count; element++)
            {
                double[] elementInternalGlobalForcesVector = ElementsAssembly[element].CreateInternalGlobalForcesVector();
                elementsInternalForces.Add(element, elementInternalGlobalForcesVector);

            }
            return elementsInternalForces;
        }

        public List<string> GetElementsType()
        {
            List<string> elementTypes = new List<string>();
            for (int element = 1; element <= ElementsConnectivity.Count; element++)
            {
                elementTypes.Add(ElementsAssembly[element].GetType().ToString());
            }
            return elementTypes;
        }

        public static Tuple<double[], double[]> NodalCoordinatesToVectors(Dictionary<int, INode> nodesList)
        {
            int totalNodes = nodesList.Count;
            double[] xCoorVector = new double[totalNodes];
            double[] yCoorVector = new double[totalNodes];
            for (int i = 1; i <= totalNodes; i++)
            {
                xCoorVector[i - 1] = nodesList[i].XCoordinate;
                yCoorVector[i - 1] = nodesList[i].YCoordinate;
            }
            return new Tuple<double[], double[]>(xCoorVector, yCoorVector);
        }
        public Dictionary<int, List<double[]>> GetElementsStresses(double[] totalDisplacementVector)
        {
            UpdateDisplacements(totalDisplacementVector);
            Dictionary<int, List<double[]>> elementsStresses = new Dictionary<int, List<double[]>>();
            for (int element = 1; element <= ElementsConnectivity.Count; element++)
            {
                List<double[]> elementStress = ElementsAssembly[element].GetStressVector();
                elementsStresses.Add(element, elementStress);
            }
            return elementsStresses;
        }
        public Dictionary<int, List<double[]>> GetElementsGaussPoints(double[] totalDisplacementVector)
        {
            UpdateDisplacements(totalDisplacementVector);
            Dictionary<int, List<double[]>> elementsGaussPoints = new Dictionary<int, List<double[]>>();
            for (int element = 1; element <= ElementsConnectivity.Count; element++)
            {
                List<double[]> elementStress = ElementsAssembly[element].GetGaussPointsInPhysicalSpace();
                elementsGaussPoints.Add(element, elementStress);
            }
            return elementsGaussPoints;
        }
        public Dictionary<int, List<double[]>> GetElementsStains(double[] totalDisplacementVector)
        {
            UpdateDisplacements(totalDisplacementVector);
            Dictionary<int, List<double[]>> elementsStains = new Dictionary<int, List<double[]>>();
            for (int element = 1; element <= ElementsConnectivity.Count; element++)
            {
                List<double[]> elementStress = ElementsAssembly[element].GetStrainVector();
                elementsStains.Add(element, elementStress);
            }
            return elementsStains;
        }
        public Dictionary<int, List<double[]>> GetElementsNodesStains(double[] totalDisplacementVector)
        {
            UpdateDisplacements(totalDisplacementVector);
            Dictionary<int, List<double[]>> elementsStains = new Dictionary<int, List<double[]>>();
            for (int element = 1; element <= ElementsConnectivity.Count; element++)
            {
                List<double[]> elementStress = ElementsAssembly[element].GetStrainFromElementsNodes();
                elementsStains.Add(element, elementStress);
            }
            return elementsStains;
        }
        public Dictionary<int, List<double[]>> GetElementsNodesStresses(double[] totalDisplacementVector)
        {
            UpdateDisplacements(totalDisplacementVector);
            Dictionary<int, List<double[]>> elementsStresses = new Dictionary<int, List<double[]>>();
            for (int element = 1; element <= ElementsConnectivity.Count; element++)
            {
                List<double[]> elementStress = ElementsAssembly[element].GetStressFromElementsNodes();
                elementsStresses.Add(element, elementStress);
            }
            return elementsStresses;
        }
        public double[] MMCPGCreateReducedFromFullVector(double[] fullVector) 
        {
            int[] contactDof = ContactDOFsVector;
            int[] noContactDof = NoContactDOFsVector;
            double[] reducedInternalForcesVector = BoundaryConditionsImposition.MMCPCGReducedVector(fullVector, BoundedDOFsVector, contactDof, noContactDof).Item3;
            return reducedInternalForcesVector;
        }

        public void CalculateEASMatrices()
        {
            for (int element = 1; element <= ElementsConnectivity.Count; element++)
            {
                if (ElementsAssembly[element].Properties.ElementType == "ANSSolidShell8EAS4NL" ||
                    ElementsAssembly[element].Properties.ElementType == "ANSSolidShell8LEAS1RI" ||
                    ElementsAssembly[element].Properties.ElementType == "ANSSolidShell8LEAS1FI" ||
                    ElementsAssembly[element].Properties.ElementType == "ANSSolidShell8LEAS7")
                {
                    ElementsAssembly[element].CalculateElementEASMatrices();
                }
            }
        }
        public void InitializeEASParameters()
        {
            for (int element = 1; element <= ElementsConnectivity.Count; element++)
            {
                if (ElementsAssembly[element].Properties.ElementType == "ANSSolidShell8EAS4NL" ||
                    ElementsAssembly[element].Properties.ElementType == "ANSSolidShell8LEAS1RI" ||
                    ElementsAssembly[element].Properties.ElementType == "ANSSolidShell8LEAS1FI" ||
                    ElementsAssembly[element].Properties.ElementType == "ANSSolidShell8LEAS7")
                {
                    ElementsAssembly[element].InitializeElementEASParameters();
                }
            }
        }
        public void UpdateEASParameters(double[] deltaU)
        {
            double[] fullTotalDisplacementVector = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(deltaU,
                BoundedDOFsVector);
            for (int element = 1; element <= ElementsConnectivity.Count; element++)
            {
                if (ElementsAssembly[element].Properties.ElementType == "ANSSolidShell8EAS4NL" ||
                    ElementsAssembly[element].Properties.ElementType == "ANSSolidShell8LEAS1RI" ||
                    ElementsAssembly[element].Properties.ElementType == "ANSSolidShell8LEAS1FI" ||
                    ElementsAssembly[element].Properties.ElementType == "ANSSolidShell8LEAS7")
                {
                    int elementDofs = ElementsAssembly[element].ElementFreedomList.Count;
                    double[] elementDisplacementVector = new double[elementDofs];
                    for (int i = 0; i < elementDofs; i++)
                    {
                        int localRow = i;
                        int globalRow = ElementsAssembly[element].ElementFreedomList[i];
                        elementDisplacementVector[localRow] = fullTotalDisplacementVector[globalRow];
                    }
                    ElementsAssembly[element].UpdateElementEASParameters(elementDisplacementVector);
                }
            }
        }
        public void StoreFinalStepDisplacementVector(double[] solutionVector)
        {
            double[] fullTotalDisplacementVector = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(solutionVector,
                BoundedDOFsVector);
            for (int element = 1; element <= ElementsConnectivity.Count; element++)
            {
                if (ElementsAssembly[element].Properties.ElementType == "ANSSolidShell8EAS4NL" ||
                    ElementsAssembly[element].Properties.ElementType == "ANSSolidShell8LEAS1RI" ||
                    ElementsAssembly[element].Properties.ElementType == "ANSSolidShell8LEAS1FI" ||
                    ElementsAssembly[element].Properties.ElementType == "ANSSolidShell8LEAS7")
                {
                    int elementDofs = ElementsAssembly[element].ElementFreedomList.Count;
                    double[] elementDisplacementVector = new double[elementDofs];
                    for (int i = 0; i < elementDofs; i++)
                    {
                        int localRow = i;
                        int globalRow = ElementsAssembly[element].ElementFreedomList[i];
                        elementDisplacementVector[localRow] = fullTotalDisplacementVector[globalRow];
                    }
                    ElementsAssembly[element].StoreElementFinalStepDisplacementVector(elementDisplacementVector);
                }
            }
        }
    }
}
