using System;
using System.Collections.Generic;
using System.Globalization;

namespace BumperCollitionSimulation
{
    class Program
    {        
        static void Main(string[] args)
        {
            //Read input as double[]
            double[] inputValues = new double[args.Length];
            int numberOfTimeSteps = 10;
            int numberOfDOFs = 5134 * 3;
            try
            {
                NumberFormatInfo provider = new NumberFormatInfo();
                provider.NumberDecimalSeparator = ".";
                provider.NumberGroupSeparator = ",";
                for (int i = 0; i < args.Length; i++)
                {
                    inputValues[i] = Convert.ToDouble(args[i], provider);
                }
            }
            catch (Exception ex)
            {
                Console.WriteLine(ex.Message);
                return;
            }

            // Kostas code for analysis
            double modulus = inputValues[0];
            //double modulus = 200.0 * 1e9;
            //Console.WriteLine("Attempting to read files.");
            ImportedData importedData = new ImportedData();
            importedData.ImportSlaveSurfaceConnectivity("../../SlaveSurfaceConnectivity.txt");
            importedData.ImportMasterSurfaceConnectivity("../../MasterSurfaceConnectivity.txt");
            importedData.ImportNodes("../../BumperImpactNodes.txt");
            importedData.ImportFixedNodes("../../BumperImpactFixedNodesList.txt");
            importedData.ImportElementsConnectivity("../../BumperImpactConnectivity.txt");
			
            //Console.WriteLine("Files have been imported successfully.");

            var results = TruckBumperImpact.RunDynamicExample(importedData, modulus);   
            			
			//Output to be passed to melissa
            for (int i = 0; i < numberOfTimeSteps; i++)
            {
                for (int dof = 0; dof < numberOfDOFs; dof++)
                {
					Console.WriteLine(results[i][dof]);
                }
            }
        }
    }
}
