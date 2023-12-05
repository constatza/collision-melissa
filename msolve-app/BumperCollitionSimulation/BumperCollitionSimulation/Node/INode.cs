using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace BumperCollitionSimulation
{
    public interface INode
    {
        double XCoordinate { get; set; }
        double YCoordinate { get; set; }
        double ZCoordinate { get; set; }
        bool UXDof { get; set; }
        bool UYDof { get; set; }
        bool UZDof { get; set; }
        bool RXDof { get; set; }
        bool RYDof { get; set; }
        bool RZDof { get; set; }
        double RX { get; set; }
        double RY { get; set; }
        double RZ { get; set; }
    }
}
