using SMesh;

namespace smesh_tests
{
    [TestClass]
    public class MathTests
    {
        const double EPSILON = 0.00001;

        [TestMethod]
        public void PlaneTransformations()
        {
            Vector3 a = new Vector3(2, 2, 2);
            Vector3 b = new Vector3(3, 6, 2);
            Vector3 c = new Vector3(7, 6, 2);

            Vector3 p = new Vector3(3, 5, 2);

            Plane plane = SMMath.PlaneFrom3Points(a, b, c);
            Matrix to3d;
            Matrix to2d;
            SMMath.PlaneTransformations(plane, a, b, out to3d, out to2d);

            SMMath.TransformTo2D(to2d, a);
            SMMath.TransformTo2D(to2d, b);
            SMMath.TransformTo2D(to2d, c);
            SMMath.TransformTo2D(to2d, p);
            Assert.IsTrue(SMMath.Vector3Distance(SMMath.TransformTo3D(to3d, SMMath.TransformTo2D(to2d, a)), a) < EPSILON);
            Assert.IsTrue(SMMath.Vector3Distance(SMMath.TransformTo3D(to3d, SMMath.TransformTo2D(to2d, b)), b) < EPSILON);
            Assert.IsTrue(SMMath.Vector3Distance(SMMath.TransformTo3D(to3d, SMMath.TransformTo2D(to2d, c)), c) < EPSILON);
            Assert.IsTrue(SMMath.Vector3Distance(SMMath.TransformTo3D(to3d, SMMath.TransformTo2D(to2d, p)), p) < EPSILON);
            Assert.IsTrue(SMMath.Vector3Distance(SMMath.TransformTo3D(to3d, new Vector2(0, 0)), a) < EPSILON);
        }


        [TestMethod]
        public void Lines3dIntersection()
        {
            var sga = new Segment3(new Vector3(-11, 0, -1), new Vector3(-9, 0, -1));
            var sgb = new Segment3(new Vector3(-10, -1, 1), new Vector3(-10, 1, 1));

            Vector3 pta;
            Vector3 ptb;
            Assert.IsTrue(SMMath.LineLineClosestPoint(sga, sgb, out pta, out ptb, EPSILON));
            Assert.IsTrue(SMMath.Vector3Distance(pta, new Vector3(-10, 0, -1)) == 0);
            Assert.IsTrue(SMMath.Vector3Distance(ptb, new Vector3(-10, 0, 1)) == 0);

            sga.A = new Vector3(2, 6, -9);
            sga.B = SMMath.Vector3Add(sga.A, new Vector3(3, 4, -4));
            sgb.A = new Vector3(-1, -2, 3);
            sgb.B = SMMath.Vector3Add(sgb.A, new Vector3(2, -6, 1));

            Assert.IsTrue(SMMath.LineLineClosestPoint(sga, sgb, out pta, out ptb, EPSILON));
            Assert.IsTrue(SMMath.Vector3Distance(pta, ptb) - 4.740201166731855 < EPSILON);
        }
    }
}