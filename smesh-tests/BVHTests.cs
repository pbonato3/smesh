using SMesh;
using System;
using System.IO;

namespace smesh_tests
{
    [TestClass]
    public class BVHTests
    {
    const double EPSILON = 0.00001;

        [TestMethod]
        public void PlaneTransformations()
        {
            string workingDirectory = Environment.CurrentDirectory;
            // or: Directory.GetCurrentDirectory() gives the same result

            // This will get the current PROJECT bin directory (ie ../bin/)
            string solutionDir = Directory.GetParent(Directory.GetParent(workingDirectory).Parent.FullName).Parent.FullName;

            var mesh = SMObj.ParseFile(solutionDir + "/test/sample.obj");
            var bvh = SphereBVH.BuildBVH(mesh);

            var searchPoint = new Vector3(3, -2, 14);
            var neigh = bvh.Search(searchPoint, 0);
            Console.WriteLine(Vector3String(searchPoint));
            Assert.IsTrue(neigh.Count == 12);
        }



        static string Vector2String(Vector2 vec)
        {
            string str = "(";
            str += vec.X.ToString();
            str += ", ";
            str += vec.Y.ToString();
            str += ")";
            return str;
        }

        static string Vector3String(Vector3 vec)
        {
            string str = "(";
            str += vec.X.ToString();
            str += ", ";
            str += vec.Y.ToString();
            str += ", ";
            str += vec.Z.ToString();
            str += ")";
            return str;
        }

        static void RecursiveBVHToString(SphereBVH.SphereNode currNode, int level)
        {

            string spaces = "";
            for (int i = 0; i < level; i++)
            {
                spaces += " _ ";
            }

            if (currNode.Id >= 0)
            {
                System.Console.WriteLine(spaces + currNode.Id + "\n");
                return;
            }

            System.Console.WriteLine(spaces + "{\n");

            for (int i = 0; i < currNode.Children.Count; i++)
            {
                RecursiveBVHToString(currNode.Children[i], level + 1);
            }

            System.Console.WriteLine(spaces + "}\n");

        }

    }
}