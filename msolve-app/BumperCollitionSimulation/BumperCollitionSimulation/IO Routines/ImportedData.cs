using System;
using System.Collections.Generic;
using System.IO;
using System.Text;

namespace BumperCollitionSimulation
{
    public class ImportedData
    {
        public Dictionary<int, Dictionary<int, int>> elementsConnectivity;
        public Dictionary<int, Dictionary<int, int>> masterConnectivity;
        public Dictionary<int, Dictionary<int, int>> slaveConnectivity;
        public Dictionary<int, INode> nodes = new Dictionary<int, INode>();
        public List<int> loadedNodes = new List<int>();
        public List<int> fixedNodes = new List<int>();

        public ImportedData()
        {
            masterConnectivity = new Dictionary<int, Dictionary<int, int>>();
            elementsConnectivity = new Dictionary<int, Dictionary<int, int>>();
            slaveConnectivity = new Dictionary<int, Dictionary<int, int>>();
            nodes = new Dictionary<int, INode>();
            fixedNodes = new List<int>();
        }

        public void ImportSlaveSurfaceConnectivity(string path)
        {
            try
            {
                StreamReader stream = new StreamReader(path);
                string file = stream.ReadToEnd();
                List<string> lines = new List<string>(file.Split(new string[] { "\r\n", "\r", "\n" }, StringSplitOptions.RemoveEmptyEntries));
                //lines.RemoveAt(0);
                int elementIndex = 1;
                foreach (var line in lines)
                {
                    // in case of first line ...
                    string[] fields = line.Split(new string[] { "\t" }, StringSplitOptions.None);
                    //int elementIndex = int.Parse(fields[0]);

                    var element = new Dictionary<int, int>()
                        {
                            { 1, int.Parse(fields[0]) },
                            { 2, int.Parse(fields[1]) },
                            { 3, int.Parse(fields[2]) },
                            { 4, int.Parse(fields[3]) }
                        };
                    slaveConnectivity[elementIndex] = element;
                    elementIndex += 1;
                }
            }
            catch (Exception ex)
            {
                Console.WriteLine("Houston we have a problem 1");
				throw ex;
            }
        }

        public void ImportMasterSurfaceConnectivity(string path)
        {
            try
            {
                StreamReader stream = new StreamReader(path);
                string file = stream.ReadToEnd();
                List<string> lines = new List<string>(file.Split(new string[] { "\r\n", "\r", "\n" }, StringSplitOptions.RemoveEmptyEntries));
                //lines.RemoveAt(0);
                int elementIndex = 1;
                foreach (var line in lines)
                {
                    // in case of first line ...
                    string[] fields = line.Split(new string[] { "\t" }, StringSplitOptions.None);
                    //int elementIndex = int.Parse(fields[0]);

                    var element = new Dictionary<int, int>()
                        {
                            { 1, int.Parse(fields[0]) },
                            { 2, int.Parse(fields[1]) },
                            { 3, int.Parse(fields[2]) },
                            { 4, int.Parse(fields[3]) }
                        };
                    masterConnectivity[elementIndex] = element;
                    elementIndex += 1;
                }
            }
            catch (Exception ex)
            {
				Console.WriteLine("Houston we have a problem 2");
                throw ex;
            }
        }

        public void ImportElementsConnectivity(string path)
        {
            try
            {
                StreamReader stream = new StreamReader(path);
                string file = stream.ReadToEnd();
                List<string> lines = new List<string>(file.Split(new string[] { "\r\n", "\r", "\n" }, StringSplitOptions.RemoveEmptyEntries));
                //lines.RemoveAt(0);
                int elementIndex = 1;
                foreach (var line in lines)
                {
                    // in case of first line ...
                    string[] fields = line.Split(new string[] { "\t" }, StringSplitOptions.None);
                    //int elementIndex = int.Parse(fields[0]);

                    var element = new Dictionary<int, int>()
                        {
                            { 1, int.Parse(fields[0]) },
                            { 2, int.Parse(fields[1]) },
                            { 3, int.Parse(fields[2]) },
                            { 4, int.Parse(fields[3]) },
                            { 5, int.Parse(fields[4]) },
                            { 6, int.Parse(fields[5]) },
                            { 7, int.Parse(fields[6]) },
                            { 8, int.Parse(fields[7]) }
                        };
                    elementsConnectivity[elementIndex] = element;
                    elementIndex += 1;
                }
            }
            catch (Exception ex)
            {
				Console.WriteLine("Houston we have a problem 3");
                throw ex;
            }
        }

        public void ImportNodes(string path)
        {
            try
            {
                StreamReader stream = new StreamReader(path);
                string file = stream.ReadToEnd();
                List<string> lines = new List<string>(file.Split(new string[] { "\r\n", "\r", "\n" }, StringSplitOptions.RemoveEmptyEntries));
                //lines.RemoveAt(0);
                int nodeIndex = 1;
                foreach (var line in lines)
                {
                    // in case of first line ...
                    string[] fields = line.Split(new string[] { "\t" }, StringSplitOptions.None);
                    //int nodeIndex = int.Parse(fields[0]);
                    var node = new Node(double.Parse(fields[0]), double.Parse(fields[1]), double.Parse(fields[2]));
                    nodes[nodeIndex] = node;
                    nodeIndex += 1;
                }
            }
            catch (Exception ex)
            {
				Console.WriteLine("Houston we have a problem 4");
                throw ex;
            }
        }

        public void ImportLoadedNodes(string path)
        {
            try
            {
                StreamReader stream = new StreamReader(path);
                string file = stream.ReadToEnd();
                List<string> lines = new List<string>(file.Split(new string[] { "\r\n", "\r", "\n" }, StringSplitOptions.RemoveEmptyEntries));
                //lines.RemoveAt(0);
                int nodeIndex = 1;
                foreach (var line in lines)
                {
                    // in case of first line ...
                    string[] fields = line.Split(new string[] { "\t" }, StringSplitOptions.None);
                    //int nodeIndex = int.Parse(fields[0]);
                    var node = (int)(double.Parse(fields[0]));
                    loadedNodes.Add(node);
                    nodeIndex += 1;
                }
            }
            catch (Exception ex)
            {
				Console.WriteLine("Houston we have a problem 5");
                throw ex;
            }
        }
        public void ImportFixedNodes(string path)
        {
            try
            {
                StreamReader stream = new StreamReader(path);
                string file = stream.ReadToEnd();
                List<string> lines = new List<string>(file.Split(new string[] { "\r\n", "\r", "\n" }, StringSplitOptions.RemoveEmptyEntries));
                //lines.RemoveAt(0);
                int nodeIndex = 1;
                foreach (var line in lines)
                {
                    // in case of first line ...
                    string[] fields = line.Split(new string[] { "\t" }, StringSplitOptions.None);
                    //int nodeIndex = int.Parse(fields[0]);
                    var node = (int)(double.Parse(fields[0]));
                    fixedNodes.Add(node);
                    nodeIndex += 1;
                }
            }
            catch (Exception ex)
            {
                throw ex;
            }
        }
    }
}
