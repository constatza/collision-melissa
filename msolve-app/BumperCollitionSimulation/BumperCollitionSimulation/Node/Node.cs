using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace BumperCollitionSimulation
{
    class Node : INode
    {
        public double XCoordinate { get; set; }
        public double YCoordinate { get; set; }
        public double ZCoordinate { get; set; }
        public bool UXDof { get; set; }
        public bool UYDof { get; set; }
        public bool UZDof { get; set; }
        public bool RXDof { get; set; }
        public bool RYDof { get; set; }
        public bool RZDof { get; set; }
        public double RX { get; set; }
        public double RY { get; set; }
        public double RZ { get; set; }
        public Node(double xCoordinate, double yCoordinate)
        {
            XCoordinate = xCoordinate;
            YCoordinate = yCoordinate;
            UXDof = true;
            UYDof = true;
            RZDof = true;
        }

        public Node(double xCoordinate, double yCoordinate, double zCoordinate)
        {
            XCoordinate = xCoordinate;
            YCoordinate = yCoordinate;
            ZCoordinate = zCoordinate;
            UXDof = true;
            UYDof = true;
            UZDof = true;
            RXDof = true;
            RYDof = true;
            RZDof = true;
        }
        public Node(double xCoordinate, double yCoordinate, double zCoordinate, double v1Rotation, double v2Rotation)
        {
            XCoordinate = xCoordinate;
            YCoordinate = yCoordinate;
            ZCoordinate = zCoordinate;
            RX = v1Rotation;
            RY = v2Rotation;
            UXDof = true;
            UYDof = true;
            UZDof = true;
            RXDof = true;
            RYDof = true;
        }
        public Node(double xCoordinate, double yCoordinate, double zCoordinate, double xRotation, double yRotation, double zRotation)
        {
            XCoordinate = xCoordinate;
            YCoordinate = yCoordinate;
            ZCoordinate = zCoordinate;
            RX = xRotation;
            RY = yRotation;
            RZ = zRotation;
            UXDof = true;
            UYDof = true;
            UZDof = true;
            RXDof = true;
            RYDof = true;
            RZDof = true;
        }
    }
}
