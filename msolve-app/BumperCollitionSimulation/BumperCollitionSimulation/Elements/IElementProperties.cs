using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace BumperCollitionSimulation
{
    public interface IElementProperties
    {
        double YoungMod { get; set; }
        double SectionArea { get; set; }
        double MomentOfInertia { get; set; }
        string ElementType { get; set; }
        double Density { get; set; }
        double Thickness { get; set; }
        double ThermalConductivity { get; set; }
        double ContactForceValue { get; set; }
        double ContactThermalConductivity { get; set; }
        double SurfaceRoughness { get; set; }
        double YieldStrength { get; set; }
        double Dx1 { get; set; }
        double Dx2 { get; set; }
        double Dx { get; set; }
        double A { get; set; }
        double B { get; set; }
        double PoissonRatio { get; set; }
        int MasterSegmentPolynomialDegree { get; set; }
        int SlaveSegmentPolynomialDegree { get; set; }
        int IntegrationPoints { get; set; }
        double PenaltyFactorRatio { get; set; }
        double TangentialPenaltyFactorRatio { get; set; }
        double StickingCoefficient { get; set; }
        double SlidingCoefficient { get; set; }
        Dictionary<int, double> AllIntegrationPointsStickingPoints { get; set; }
        Dictionary<int, double> AllIntegrationPointsTangentialTractions { get; set; }
        Dictionary<int, double> AllIntegrationPointsStickingPoints2 { get; set; }
        Dictionary<int, double> AllIntegrationPointsTangentialTractions2 { get; set; }
        Dictionary<int, double[]> AllIntegrationPointsSurfaceVectors1 { get; set; }
        Dictionary<int, double[]> AllIntegrationPointsSurfaceVectors2 { get; set; }
        Dictionary<int, double> StoredTrialTangentialTractions1 { get; set; }
        Dictionary<int, double> StoredTrialTangentialTractions2 { get; set; }
        Dictionary<int, double> StoredIntegrationPointsStickingPoints { get; set; }
        Dictionary<int, double> StoredIntegrationPointsStickingPoints2 { get; set; }
        List<double> DU { get; set; }
        List<double> Xprev { get; set; }

        Dictionary<int, double> DegenerateShellNodalThicknesses { get; set; }
        Dictionary<int, double[]> DegenerateShellNodalNormalUnitVectors { get; set; }
        double[] EASVector { get; set; }
        double[] EASFEnhancedVector { get; set; }
        double[,] EASLMatrix { get; set; }
        double[,] EASDMatrix { get; set; }
        double[] DisplacementVectorPreviousStep { get; set; }
        double[] DisplacementVectorPreviousIncrement { get; set; }
        double[] HourglassInternalForceVector { get; set; }
        double[] HourglassInternalForceVectorPreviousConvergedSolution { get; set; }

        //Dictionary<int, double> TangentMatrixInitStickPointUpdSlideCase1 { get; set; }
        //Dictionary<int, double> TangentMatrixinitStickPointUpdSlideCase2 { get; set; }
        //Dictionary<int, double> ResidualInitStickPointUpdSlideCase1 { get; set; }
        //Dictionary<int, double> ResidualInitStickPointUpdSlideCase2 { get; set; }
    }
}

