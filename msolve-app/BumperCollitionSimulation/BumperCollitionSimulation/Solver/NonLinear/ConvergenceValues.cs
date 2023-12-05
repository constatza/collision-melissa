using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BumperCollitionSimulation
{
    public class ConvergenceValues
    {
        public double ResidualNorm { get; set; }
        public double Tolerance { get; set; }
        public int LoadStep { get; set; }
        public int Iteration { get; set; }
        public bool ConvergenceResult { get; set; }
    }
}
