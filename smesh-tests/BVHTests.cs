using SMesh;

namespace smesh_tests
{
    [TestClass]
    public class BVHTests
    {
        const double EPSILON = 0.00001;

        [TestMethod]
        public void PlaneTransformations()
        {
            var mesh = SMObj.ParseFile("/test/sample.obj");
            var bvh = SphereBVH.BuildBVH(mesh);

            var searchPoint = new Vector3(3, -2, 14);
            var neigh = bvh.Search(searchPoint, 0);
            Console.WriteLine(Vector3String(searchPoint));
        }
    }
}