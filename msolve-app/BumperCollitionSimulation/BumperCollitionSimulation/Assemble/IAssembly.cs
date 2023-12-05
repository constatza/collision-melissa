using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace BumperCollitionSimulation
{
    public interface IAssembly
    {
        Dictionary<int, IElementProperties> ElementsProperties { get; set; }
        Dictionary<int, INode> Nodes { get; set; }
        Dictionary<int, Dictionary<int, int>> ElementsConnectivity { get; set; }
        Dictionary<int, bool[]> NodeFreedomAllocationList { get; set; }
        void CreateElementsAssembly();
        double[,] CreateTotalStiffnessMatrix();
        double[,] CreateStiffnessMatrixLinearPart();
        double[,] UpdateStifnessMatrix(double[,] stiffnessMatrixLinearPart);
        bool ActivateBoundaryConditions { get; set; }
        int[] BoundedDOFsVector { get; set; }
        int[] ContactDOFsVector { get; set; }
        int[] NoContactDOFsVector { get; set; }

        void UpdateDisplacements(double[] totalDisplacementVector);
        void UpdateAccelerations(double[] totalAccelerationsVector);
        double[] CreateTotalInternalForcesVector();
        double[,] CreateTotalMassMatrix();
        double[,] CreateTotalDampingMatrix();
        Dictionary<int, double[]> GetElementsInternalForces(double[] totalInternalForcesVector);
        Dictionary<int, List<double[]>> GetElementsStresses(double[] totalDisplacementVector);
        Dictionary<int, List<double[]>> GetElementsStains(double[] totalDisplacementVector);
        Dictionary<int, List<double[]>> GetElementsNodesStresses(double[] totalDisplacementVector);
        Dictionary<int, List<double[]>> GetElementsNodesStains(double[] totalDisplacementVector);
        Dictionary<int, List<double[]>> GetElementsGaussPoints(double[] totalDisplacementVector);
        void InitializeContactTangentialProperties();
        void UpdateContactTangentialProperties();
        void InitializeContactSurfaceVectors();
        void UpdateContactSurfaceVectors();
        void UpdateElementsIncrementalDisplacements(double[] deltaU);
        List<string> GetElementsType();
        Dictionary<int, IElement> ElementsAssembly { get; set; }
        int CountElementsOfSameType(Type elementType);
        Tuple<double[,], double[,], double[,], double[,]> MMCPCGCreateTotalStiffnessMatrix();
        double[] MMCPCGCreateTotalInternalForcesVectors();
        void SeperateContactDoF();
        Tuple<double[], double[]> MMCPCGSeperateReducedExternalForcesVectors(double[] extForces);
        void MMCPCGUpdateDisplacements(double[] totalDisplacementVector);
        void MMCPCGUpdateElementsIncrementalDisplacements(double[] deltaU);
        double[] MMCPCGRearrange(double[] vector);
        double[] MMCPGCreateReducedFromFullVector(double[] fullVector);
        void CalculateEASMatrices();
        void InitializeEASParameters();
        void UpdateEASParameters(double[] solutionVector);
        void StoreFinalStepDisplacementVector(double[] solutionVector);
        bool ActivateParallelCalculations { get; set; }


        //void UpdateValues(double[] totalDisplacementVector);
        //double[,] CreateTotalStiffnessMatrix();
        //double[,] CreateTotalMassMatrix();
        //double[] CreateTotalInternalForcesVector();
        //int[] BoundedDOFsVector
        //{
        //    get;
        //    set;
        //}
    }
}