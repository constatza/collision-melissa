using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace BumperCollitionSimulation
{
    public class ElementProperties : IElementProperties
    {
        public double YoungMod { get; set; }
        public double SectionArea { get; set; }
        public double MomentOfInertia { get; set; }
        public string ElementType { get; set; }
        public double Density { get; set; }
        public double Thickness { get; set; }
        public double ThermalConductivity { get; set; }
        public double ContactForceValue { get; set; }

        public double ContactThermalConductivity { get; set; }
        public double SurfaceRoughness { get; set; }
        public double YieldStrength { get; set; }
        public double Dx1 { get; set; }
        public double Dx2 { get; set; }
        public double Dx { get; set; }
        public double A { get; set; }
        public double B { get; set; }
        public double PoissonRatio { get; set; }
        public int MasterSegmentPolynomialDegree { get; set; }
        public int SlaveSegmentPolynomialDegree { get; set; }
        public int IntegrationPoints { get; set; }
        public double PenaltyFactorRatio { get; set; }
        public double TangentialPenaltyFactorRatio { get; set; }
        public double StickingCoefficient { get; set; }
        public double SlidingCoefficient { get; set; }
        public Dictionary<int, double> AllIntegrationPointsStickingPoints { get; set; }
        public Dictionary<int, double> AllIntegrationPointsTangentialTractions { get; set; }
        public Dictionary<int, double> AllIntegrationPointsStickingPoints2 { get; set; }
        public Dictionary<int, double> AllIntegrationPointsTangentialTractions2 { get; set; }
        public Dictionary<int, double[]> AllIntegrationPointsSurfaceVectors1 { get; set; }
        public Dictionary<int, double[]> AllIntegrationPointsSurfaceVectors2 { get; set; }
        public Dictionary<int, double> StoredTrialTangentialTractions1 { get; set; }
        public Dictionary<int, double> StoredTrialTangentialTractions2 { get; set; }
        public Dictionary<int, double> StoredIntegrationPointsStickingPoints { get; set; }
        public Dictionary<int, double> StoredIntegrationPointsStickingPoints2 { get; set; }
        public Dictionary<int, double> TangentMatrixInitStickPointUpdSlideCase1 { get; set; }
        public Dictionary<int, double> TangentMatrixinitStickPointUpdSlideCase2 { get; set; }
        public Dictionary<int, double> ResidualInitStickPointUpdSlideCase1 { get; set; }
        public Dictionary<int, double> ResidualInitStickPointUpdSlideCase2 { get; set; }
        public List<double> DU { get; set; }
        public List<double> Xprev { get; set; }
        public Dictionary<int, double> DegenerateShellNodalThicknesses { get; set; }
        public Dictionary<int, double[]> DegenerateShellNodalNormalUnitVectors { get; set; }
        public double[] EASVector { get; set; }
        public double[] DisplacementVectorPreviousStep { get; set; }

        public double[] EASFEnhancedVector { get; set; }
        public double[,] EASLMatrix { get; set; }
        public double[,] EASDMatrix { get; set; }
        public double[] DisplacementVectorPreviousIncrement { get; set; }
        public double[] HourglassInternalForceVector { get; set; }
        public double[] HourglassInternalForceVectorPreviousConvergedSolution { get; set; }
        public ElementProperties(double youngMod, double sectionArea, string elementType)
        {
            YoungMod = youngMod;
            SectionArea = sectionArea;
            ElementType = elementType;
            PoissonRatio = 0.30;
            PenaltyFactorRatio = 10.0;
        }
        /// <summary>
        /// New constructor after the inclusion of Poissons ratio in element properties
        /// </summary>
        /// <param name="youngMod"></param>
        /// <param name="poissonRatio"></param>
        /// <param name="sectionArea"></param>
        /// <param name="thickness"></param>
        /// <param name="density"></param>
        /// <param name="elementType"></param>
        public ElementProperties(double youngMod, double poissonRatio, double sectionArea, double thickness, double density, string elementType)
        {
            YoungMod = youngMod;
            SectionArea = sectionArea;
            ElementType = elementType;
            PoissonRatio = poissonRatio;
            Thickness = thickness;
            Density = density;
            PenaltyFactorRatio = 10.0;
        }


        public ElementProperties(double youngMod, string elementType)
        {
            YoungMod = youngMod;
            ElementType = elementType;
        }


        public ElementProperties(double youngMod, double sectionArea, double momentOfInertia, string elementType)
        {
            YoungMod = youngMod;
            SectionArea = sectionArea;
            MomentOfInertia = momentOfInertia;
            ElementType = elementType;
            PoissonRatio = 0.30;
            PenaltyFactorRatio = 10.0;
        }
        /// <summary>
        /// Constructor for Frictional NtS contact elements.
        /// </summary>
        public ElementProperties(double youngMod, string elementType, double normalPenaltyFactorMultiplier, double tangentialPenaltyFactorMultiplier,
            double stickingCoefficient, double slidingCoefficient)
        {
            YoungMod = youngMod;
            PenaltyFactorRatio = normalPenaltyFactorMultiplier;
            TangentialPenaltyFactorRatio = tangentialPenaltyFactorMultiplier;
            ElementType = elementType;
            StickingCoefficient = stickingCoefficient;
            SlidingCoefficient = slidingCoefficient;
        }
        /// <summary>
        /// New constructor for StS2D contact element. Properties are added in order to 
        /// include higher order elements
        /// </summary>
        public ElementProperties(double youngMod, double contactArea, string elementType, double penaltyFactorRatio,
            int integrationPoints, int slaveSegmentPolynomialDegree, int masterSegmentPolynomialDegree)
        {
            if (youngMod > 0.0)
            {
                YoungMod = youngMod;

            }
            else
            {
                throw new Exception("Young's modulus must be positive real number");

            }
            SectionArea = contactArea;
            if (penaltyFactorRatio > 0.0)
            {
                PenaltyFactorRatio = penaltyFactorRatio;

            }
            else
            {
                PenaltyFactorRatio = 10.0;//default value

            }
            ElementType = elementType;
            if (integrationPoints <= 15 && integrationPoints >= 1)
            {
                IntegrationPoints = integrationPoints;

            }
            else
            {
                throw new Exception("The amount of integration points must be between 1 & 15");
            }
            if (slaveSegmentPolynomialDegree<=2 && slaveSegmentPolynomialDegree >= 1)
            {
                SlaveSegmentPolynomialDegree = slaveSegmentPolynomialDegree;

            }
            else
            {
                throw new Exception("Only polynomial degrees <= 2");
            }
            if (masterSegmentPolynomialDegree <= 2 && masterSegmentPolynomialDegree >= 1)
            {
                MasterSegmentPolynomialDegree = masterSegmentPolynomialDegree;

            }
            else
            {
                throw new Exception("Only polynomial degrees <= 2");
            }
        }
        /// <summary>
        /// New constructor for Frictional StS2D contact element.
        /// </summary>
        public ElementProperties(double youngMod, string elementType, double penaltyFactorRatio,
            int integrationPoints, int slaveSegmentPolynomialDegree, int masterSegmentPolynomialDegree, double tangentialPenaltyFactorRatio,
            double stickingCoefficient, double slidingCoefficient)
        {
            if (youngMod > 0.0)
            {
                YoungMod = youngMod;

            }
            else
            {
                throw new Exception("Young's modulus must be positive real number");

            }
            if (penaltyFactorRatio > 0.0)
            {
                PenaltyFactorRatio = penaltyFactorRatio;

            }
            else
            {
                PenaltyFactorRatio = 10.0;//default value

            }
            if (tangentialPenaltyFactorRatio > 0.0)
            {
                TangentialPenaltyFactorRatio = tangentialPenaltyFactorRatio;

            }
            else
            {
                TangentialPenaltyFactorRatio = 10.0;//default value

            }
            ElementType = elementType;
            if (integrationPoints <= 15 && integrationPoints >= 1)
            {
                IntegrationPoints = integrationPoints;

            }
            else
            {
                throw new Exception("The amount of integration points must be between 1 & 15");
            }
            if (slaveSegmentPolynomialDegree <= 2 && slaveSegmentPolynomialDegree >= 1)
            {
                SlaveSegmentPolynomialDegree = slaveSegmentPolynomialDegree;

            }
            else
            {
                throw new Exception("Only polynomial degrees <= 2");
            }
            if (masterSegmentPolynomialDegree <= 2 && masterSegmentPolynomialDegree >= 1)
            {
                MasterSegmentPolynomialDegree = masterSegmentPolynomialDegree;

            }
            else
            {
                throw new Exception("Only polynomial degrees <= 2");
            }
            StickingCoefficient = stickingCoefficient;
            SlidingCoefficient = slidingCoefficient;
        }
        public ElementProperties()
        {

        }
    }
}
