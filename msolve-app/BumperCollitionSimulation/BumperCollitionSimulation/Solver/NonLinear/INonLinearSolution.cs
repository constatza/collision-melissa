using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace BumperCollitionSimulation
{
    public interface INonLinearSolution
    {
        double[] Solve(IAssembly assembly, ILinearSolution linearScheme, double[] forceVector);
        double[] Solve(IAssembly assembly, ILinearSolution linearScheme, Dictionary<int, double[]> forceVectors,
            Dictionary<int, int> numberOfLoadStepsForEachLoad);
        //double[] Solve(IAssembly assembly, double[] totalForceVector);

        int numberOfLoadSteps { get; set; }
        Dictionary<int, double[]> InternalForces { get; set; }
        Dictionary<int, double[]> Solutions { get; set; }
        double Tolerance { get; set; }
        int MaxIterations { get; set; }
        event EventHandler<ConvergenceValues> convergenceResult;
        List<string> LoadStepConvergence { get; set; }
    }
}