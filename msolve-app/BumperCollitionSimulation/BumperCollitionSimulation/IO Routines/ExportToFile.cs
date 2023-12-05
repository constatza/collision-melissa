using System;
using System.Collections.Generic;
using System.Linq;
using System.IO;

namespace BumperCollitionSimulation
{
    public static class ExportToFile
    {
        public static void ExportExplicitResults(Dictionary<int, double[]> solution, Dictionary<int, double> timeAtEachStep, int dofNumber, int intervals)
        {
            string[] lines = new string[solution.Count / intervals];
            int step = 0;
            int line = 0;
            while (step < solution.Count - 1)
            {
                double[] sol = solution[step];
                lines[line] = line.ToString() + " " + timeAtEachStep[step].ToString() + " " + sol[dofNumber].ToString();
                line = line + 1;
                step = step + intervals;
                //if (step >= solution.Count-1)
                //{
                //    break;
                //}
            }
            //File.WriteAllLines(@"D:\WriteLines2.txt", lines);
        }

        
        public static void ExportMatlabInitialGeometry(IAssembly assembly)
        {
            List<string> coordinateData = new List<string>();
            List<string> connectivityData = new List<string>();
            List<string> contactConnectivityData = new List<string>();

            foreach (var node in assembly.Nodes)
            {
                coordinateData.Add(node.Key.ToString() + "\t" + node.Value.XCoordinate.ToString() + "\t" + node.Value.YCoordinate.ToString() + "\t" + node.Value.ZCoordinate.ToString());
            }

            foreach (var element in assembly.ElementsAssembly)
            {
                string line = element.Key.ToString();

                if (element.Value is ContactStS3D)
                {
                    foreach (var elementNode in assembly.ElementsConnectivity[element.Key])
                    {
                        line = line + "\t" + elementNode.Value.ToString();
                    }
                    contactConnectivityData.Add(line);
                }
                else
                {
                    foreach (var elementNode in assembly.ElementsConnectivity[element.Key])
                    {
                        line = line + "\t" + elementNode.Value.ToString();
                    }
                    connectivityData.Add(line);
                }
            }

            File.WriteAllLines("coordinateData.dat", coordinateData);
            File.WriteAllLines("connectivityData.dat", connectivityData);
            File.WriteAllLines("contactConnectivityData.dat", contactConnectivityData);
        }

        public static void ExportMatlabFinalGeometry(IAssembly assembly, double[] displacements)
        {
            List<string> finalCoordinateData = new List<string>();

            Dictionary<int, INode> nodesAtFinalState = Assembly.CalculateFinalNodalCoordinates(assembly.Nodes, displacements);

            foreach (var node in nodesAtFinalState)
            {
                finalCoordinateData.Add(node.Key.ToString() + "\t" + node.Value.XCoordinate.ToString() + "\t" + node.Value.YCoordinate.ToString());
            }
            File.WriteAllLines("finalCoordinateData.dat", finalCoordinateData);
        }
        public static void ExportUpdatedNodalCoordinates(IAssembly assembly, double[] displacements, string fileName)
        {
            List<string> finalCoordinateData = new List<string>();

            Dictionary<int, INode> nodesAtFinalState = Assembly.CalculateFinalNodalCoordinates(assembly.Nodes, displacements);

            foreach (var node in nodesAtFinalState)
            {
                finalCoordinateData.Add(node.Key.ToString() + "\t" + node.Value.XCoordinate.ToString() + "\t" + node.Value.YCoordinate.ToString());
            }
            File.WriteAllLines(fileName+".dat", finalCoordinateData);
        }
        
        public static void ExportContactForcesForAllLoadSteps(Dictionary<int, Dictionary<int, double[]>> allStepsContactForces)
        {
            int k = 0;
            foreach (KeyValuePair<int, Dictionary<int, double[]>> loadStep in allStepsContactForces)
            {
                k = k + 1;
                Dictionary<int, double[]> contactForcesForElements = loadStep.Value;
                int index = contactForcesForElements.Keys.Max();
                int componentsOfVector = contactForcesForElements[index].Length;
                string row;
                List<string> totalData = new List<string>();
                for (int i = 0; i < componentsOfVector; i++)
                {
                    row = "";
                    foreach (KeyValuePair<int, double[]> forcesVector in contactForcesForElements)
                    {
                        row = row + "\t" + forcesVector.Value[i];
                    }
                    totalData.Add(row);
                }
                File.WriteAllLines("ContactForces" + k + ".dat", totalData);
            }
        }

       
        public static void ExportConvergenceResultsToFile(List<string> convergenceResults)
        {
            File.WriteAllLines("ConvergenceResults.dat", convergenceResults);
        }
    }
}
