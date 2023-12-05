using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace BumperCollitionSimulation
{
    public class NonLinearSolution : INonLinearSolution
    {
        public int numberOfLoadSteps { get; set; } = 10;
        protected int[] boundaryDof;
        protected IAssembly discretization;
        protected double lambda;
        public double Tolerance { get; set; } = 1e-5;
        public int MaxIterations { get; set; } = 100;
        public bool PrintResidual { get; set; } = false;
        protected ILinearSolution linearSolver;
        public Dictionary<int, double[]> InternalForces { get; set; }
        public Dictionary<int, double[]> Solutions { get; set; }
        public event EventHandler<ConvergenceValues> convergenceResult;
        public List<string> LoadStepConvergence { get; set; }

        public virtual double[] Solve(IAssembly assembly, ILinearSolution linearScheme, double[] forceVector)
        {
            throw new Exception("LinearSolution.Solve not implemented");
        }
        public virtual double[] Solve(IAssembly assembly, ILinearSolution linearScheme, Dictionary<int, double[]> forceVectors,
            Dictionary<int, int> numberOfLoadStepsForEachLoad)
        {
            throw new Exception("LinearSolution.Solve not implemented");
        }
        //protected void OnConvergenceResult(string message)


        protected void OnConvergenceResult(ConvergenceValues convergenceValues)

        {
            if (convergenceResult != null)
            {
                convergenceResult.Invoke(this, convergenceValues);
            }
            
        }
    }
}
