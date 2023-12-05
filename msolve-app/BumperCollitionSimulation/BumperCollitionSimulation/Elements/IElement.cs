using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace BumperCollitionSimulation
{
    public interface IElement
    {
        Dictionary<int, INode> Nodes { get; }
        IElementProperties Properties { get; set; }
        Dictionary<int, bool[]> ElementFreedomSignature { get; }
        List<int> ElementFreedomList { get; set; }
        //double CalculateElementLength();
        //double CalculateElementSinus();
        //double CalculateElementCosinus();
        //double[,] CreateLocalStiffnessMatrix();
        double[,] CreateGlobalStiffnessMatrix();
        double[,] CreateMassMatrix();
        double[] DisplacementVector { get; set; }
        //Dictionary<int, double> AllIntegrationPointsStickingPoints { get; set; }
        //Dictionary<int, double> AllIntegrationPointsTangentialTractions { get; set; }
        void InitializeTangentialProperties();
        void UpdateTangentialProperties();
        void InitializeContactSurfaceGeometry();
        void UpdateContactSurfaceGeometry();
        void UpdateIncrementalDisplacements(double[] deltaU);

        double[] AccelerationVector { get; set; }
        double[] CreateInternalGlobalForcesVector();
        double[,] CreateDampingMatrix();
        Dictionary<int, INode> NodesAtFinalState();
        double ClosestPointProjection();
        List<double[]> GetStressVector();
        List<double[]> GetStrainVector();
        List<double[]> GetGaussPointsInPhysicalSpace();
        List<double[]> GetStressFromElementsNodes();
        List<double[]> GetStrainFromElementsNodes();
        List<double[]> GetStressFromElements(List<double[]> parametricCoordinatesVector);
        List<double[]> GetphysicalCoordinatesFromElements(List<double[]> parametricCoordinatesVector);
        void CalculateElementEASMatrices();
        void InitializeElementEASParameters();
        void UpdateElementEASParameters(double[] solutionVector);
        void StoreElementFinalStepDisplacementVector(double[] solutionVector);


    }
}

