using SMesh;
using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Text;
using System.Threading.Tasks;

namespace SMesh
{
    public class SMCSG
    {
        private class FaceBuilder {
            public List<Vector3> Vertices;
            public List<Vector3> Normals;
            public List<bool> Weld;

            public List<int> Triangles;

            public List<Vector2> Vertices2D;

            // Cuts, definded as vertices indices
            public List<int> Cuts;

            public AABB BBox;
            public Plane FacePlane;

            public bool FromA;
            public bool Marked;
            public bool IsInside;

            public Matrix To2D;
            public Matrix To3D;

            public Vector3 FaceVertex(int nT, int nV) {
                return Vertices[Triangles[nT * 3 + nV]];
            }
            public Vector2 FaceVertex2D(int nT, int nV)
            {
                return Vertices2D[Triangles[nT * 3 + nV]];
            }

            public Vector2 ToLocal(Vector3 v) {
                return SMMath.TransformTo2D(To2D, v);
            }

            public Vector3 ToWorld(Vector2 v)
            {
                return SMMath.TransformTo3D(To3D, v);
            }

            public bool TestCut(int idxA, int idxB, Vector3 testNormal) {
                var ptA2d = Vertices2D[idxA];
                var ptB2d = Vertices2D[idxB];
                var dir = SMMath.Vector2Subtract(ptB2d, ptA2d);
                var rightDir = new Vector2(dir.Y, -dir.X);
                var testPt2d = SMMath.Vector2Add(ptA2d, rightDir);
                var testPt3d = ToWorld(testPt2d);
                var testDir = SMMath.Vector3Subtract(testPt3d, Vertices[idxA]);

                return SMMath.Vector3Dot(testDir, testNormal) < 0;
            }

            public bool AddCut(int idxA, int idxB, bool invert) {
                int a = idxA;
                int b = idxB;
                if (invert) {
                    a = idxB;
                    b = idxA;
                }

                for (int i = 0; i < Cuts.Count - 1; i += 2) {
                    // If we have reversed cut, remove it
                    if (Cuts[i] == b && Cuts[i + 1] == a)
                    {
                        Cuts.RemoveAt(i + 1);
                        Cuts.RemoveAt(i);
                        return false;
                    }
                    // If cut already exists, there is no need to add it again
                    else if (Cuts[i] == a && Cuts[i + 1] == b) {
                        return true;
                    }

                }
                // Simply add the cut
                Cuts.Add(a);
                Cuts.Add(b);
                return true;
            }
        }



        public void Split(ref Mesh meshA, ref Mesh meshB, double tol = 0.00001) { 
            FaceBuilder[] buildersA = new FaceBuilder[meshA.FaceCount];
            FaceBuilder[] buildersB = new FaceBuilder[meshB.FaceCount];

            // init builders for mesh A
            for (int f = 0; f < meshA.FaceCount; f++)
            {
                var a = meshA.Vertices[meshA.Indices[f * 3 + 0]];
                var b = meshA.Vertices[meshA.Indices[f * 3 + 1]];
                var c = meshA.Vertices[meshA.Indices[f * 3 + 2]];

                var na = meshA.Vertices[meshA.Indices[f * 3 + 0]];
                var nb = meshA.Vertices[meshA.Indices[f * 3 + 1]];
                var nc = meshA.Vertices[meshA.Indices[f * 3 + 2]];

                buildersA[f].Vertices = new List<Vector3>(3) { a, b, c };
                buildersA[f].Normals = new List<Vector3>(3) { na, nb, nc };
                buildersA[f].Triangles = new List<int>(3) { 0, 1, 2 };

                buildersA[f].BBox = ThreePointsAABB(a, b, c);
            }

            // init builders for mesh B
            for (int f = 0; f < meshB.FaceCount; f++)
            {
                var a = meshB.Vertices[meshB.Indices[f * 3 + 0]];
                var b = meshB.Vertices[meshB.Indices[f * 3 + 1]];
                var c = meshB.Vertices[meshB.Indices[f * 3 + 2]];

                var na = meshB.Vertices[meshB.Indices[f * 3 + 0]];
                var nb = meshB.Vertices[meshB.Indices[f * 3 + 1]];
                var nc = meshB.Vertices[meshB.Indices[f * 3 + 2]];

                buildersB[f].Vertices = new List<Vector3>(3) { a, b, c };
                buildersB[f].Normals = new List<Vector3>(3) { na, nb, nc };
                buildersB[f].Triangles = new List<int>(3) { 0, 1, 2 };

                buildersB[f].BBox = ThreePointsAABB(a, b, c);
            }

            //Build bvh for mesh B
            var bvhB = SMBVH.BuildBVH(meshB);

            // Check every face of A agains B's BVH
            for (int f = 0; f < meshA.FaceCount; f++)
            {
                List<int> collisions = new List<int>();

                SMBVH.CheckCollisions(ref collisions, bvhB, buildersA[f].BBox);

                var builderA = buildersA[f];
                for (int c = 0; c < collisions.Count; c++) {
                    var builderB = buildersB[collisions[c]];

                    // Note: in Vertices 0, 1, 2 there will always be the original triangle
                    builderA.FacePlane = SMMath.PlaneFrom3Points(builderA.Vertices[0], builderA.Vertices[1], builderA.Vertices[2]);
                    builderB.FacePlane = SMMath.PlaneFrom3Points(builderB.Vertices[0], builderB.Vertices[1], builderB.Vertices[2]);

                    SplitBuilders(ref builderA, ref builderB, tol);
                }
            }
        }


        private static void AddPointsToFaceBuilder(List<Vector3> points, ref FaceBuilder builder, double tol) {
            if (points.Count == 0) {
                return;
            }

            // Get transformation to and from local 2d triangle
            SMMath.PlaneTransformations(builder.FacePlane, builder.Vertices[0], builder.Vertices[1], out builder.To3D, out builder.To2D);

            // Convert all veritces to local 2d space
            builder.Vertices2D = new List<Vector2>(builder.Vertices.Count);
            for (int i = 0; i < builder.Vertices.Count; ++i)
            {
                builder.Vertices2D[i] = SMMath.TransformTo2D(builder.To2D, builder.Vertices[i]);
            }

            // Convert all intersection points to local 2d space
            var points2D = new List<Vector2>(points.Count);
            for (int i = 0; i < points.Count; ++i)
            {
                points2D[i] = SMMath.TransformTo2D(builder.To2D, points[i]);
            }

            // Prepare a list of triangles to add and one of triangles to delte
            var trianglesToAdd = new List<int>();
            var trianglesToDelte = new List<int>();

            // Just on point
            if (points2D.Count == 1) {
                // If it is on an existing point do nothing
                for (int i = 0; i < builder.Vertices2D.Count; ++i) {
                    if (SMMath.AreEquals(points2D[0], builder.Vertices2D[0], tol)) {
                        return;
                    }
                }

                // TODO: avoid unecessary checks
                // Check if the point is on a segment
                Vector2 res;
                for (int i = 0; i < builder.Triangles.Count / 3; ++i)
                {
                    var a = builder.Triangles[i * 3 + 0];
                    var b = builder.Triangles[i * 3 + 1];
                    var c = builder.Triangles[i * 3 + 2];
                    var va = builder.Vertices2D[a];
                    var vb = builder.Vertices2D[b];
                    var vc = builder.Vertices2D[c];

                    if (SMMath.GetPointOnSegment(points2D[0], va, vb, tol, out res)) {
                        var d = builder.Vertices.Count;
                        trianglesToDelte.Add(i);
                        builder.Vertices2D.Add(res);
                        builder.Vertices.Add(SMMath.TransformTo3D(builder.To3D, res));
                        trianglesToAdd.Append(a);
                        trianglesToAdd.Append(d);
                        trianglesToAdd.Append(c);

                        trianglesToAdd.Append(b);
                        trianglesToAdd.Append(c);
                        trianglesToAdd.Append(d);
                        continue;
                    }

                    // TODO: the code above should be repeated 3 times, avoid code duplications 
                }

                // If we added to a segment, apply changes and return
                if (trianglesToAdd.Count > 0) {
                    trianglesToDelte.Sort();
                    trianglesToDelte.Reverse();
                    foreach (var idx in trianglesToDelte) {
                        builder.Triangles.RemoveAt(idx * 3 + 2);
                        builder.Triangles.RemoveAt(idx * 3 + 1);
                        builder.Triangles.RemoveAt(idx * 3 + 0);
                    }
                    foreach (var idx in trianglesToAdd) { 
                        builder.Triangles.Add(idx);
                    }
                    return;
                }

                // Check if the point is inside one of the triangles
                for (int i = 0; i < builder.Triangles.Count / 3; ++i)
                {
                    var a = builder.Triangles[i * 3 + 0];
                    var b = builder.Triangles[i * 3 + 1];
                    var c = builder.Triangles[i * 3 + 2];
                    var va = builder.Vertices2D[a];
                    var vb = builder.Vertices2D[b];
                    var vc = builder.Vertices2D[c];

                    if (SMMath.IsPointInTriangle(points2D[0], va, vb, vc)) {
                        var d = builder.Vertices.Count;
                        trianglesToDelte.Add(i);
                        builder.Vertices2D.Add(points2D[0]);
                        builder.Vertices.Add(SMMath.TransformTo3D(builder.To3D, points2D[0]));

                        builder.Triangles.RemoveAt(i * 3 + 2);
                        builder.Triangles.RemoveAt(i * 3 + 1);
                        builder.Triangles.RemoveAt(i * 3 + 0);

                        builder.Triangles.Add(a);
                        builder.Triangles.Add(b);
                        builder.Triangles.Add(d);

                        builder.Triangles.Add(b);
                        builder.Triangles.Add(c);
                        builder.Triangles.Add(d);

                        builder.Triangles.Add(c);
                        builder.Triangles.Add(a);
                        builder.Triangles.Add(d);
                        return;
                    }
                }
            }

            // TODO: this is becoming too much complex, try to simplify a little
        }

        private static AABB ThreePointsAABB(Vector3 a, Vector3 b, Vector3 c) {
            AABB bbox = new AABB();
            bbox.Min.X = Math.Min(a.X, Math.Min(b.X, c.X));
            bbox.Min.Y = Math.Min(a.Y, Math.Min(b.Y, c.Y));
            bbox.Min.Z = Math.Min(a.Z, Math.Min(b.Z, c.Z));

            bbox.Max.X = Math.Max(a.X, Math.Max(b.X, c.X));
            bbox.Max.Y = Math.Max(a.Y, Math.Max(b.Y, c.Y));
            bbox.Max.Z = Math.Max(a.Z, Math.Max(b.Z, c.Z));

            return bbox;
        }


        private static bool AddPointToBuilder(ref FaceBuilder builder, Vector3 pt, double tol, out int newPt){
            var pt2d = SMMath.TransformTo2D(builder.To2D, pt);
            // If point doesn't snap to an existing one, it will be always added as last one.
            newPt = builder.Vertices2D.Count;

            // If vertex is colose enough to an existing one do nothing
            for (int i = 0; i < builder.Vertices2D.Count; ++i) {
                if (SMMath.Vector2Distance(builder.Vertices2D[i], pt2d) <= tol) {
                    newPt = i;
                    return true;
                }
            }

            for (int i = 0; i < builder.Triangles.Count; ++i)
            {
                for (int j = 0; j < 3; ++j)
                {
                    var closest = SMMath.SegmentClosestPoint(pt2d, new Segment2(builder.FaceVertex2D(i, j), builder.FaceVertex2D(i, (j + 1) % 3)));
                    if (SMMath.Vector2Distance(closest, pt2d) > tol)
                    {
                        continue;
                    }

                    // TODO: is it correct to add the closest?
                    builder.Vertices.Add(SMMath.TransformTo3D(builder.To3D, closest));
                    builder.Vertices2D.Add(closest);
                    int vertIdx = builder.Vertices.Count - 1;

                    bool foundOther = false;
                    // Check if one of the remaining triangles shares the splitted edge
                    for (int ii = i + 1; ii < builder.Triangles.Count; ++ii)
                    {
                        for (int jj = 0; jj < 3; ++jj)
                        {
                            if (builder.Triangles[i * 3 + j] == builder.Triangles[ii * 3 + jj] &&
                                builder.Triangles[i * 3 + (j + 1) % 3] == builder.Triangles[ii * 3 + (jj + 1) % 3]) {

                                builder.Triangles.Add(builder.Triangles[ii * 3 + (jj + 2) % 3]);
                                builder.Triangles.Add(builder.Triangles[ii * 3 + jj]);
                                builder.Triangles.Add(builder.Triangles[vertIdx]);

                                builder.Triangles.Add(builder.Triangles[ii * 3 + (jj + 2) % 3]);
                                builder.Triangles.Add(builder.Triangles[vertIdx]);
                                builder.Triangles.Add(builder.Triangles[ii * 3 + (jj + 1) % 3]);

                                builder.Triangles.RemoveAt(ii * 3);
                                builder.Triangles.RemoveAt(ii * 3);
                                builder.Triangles.RemoveAt(ii * 3);

                                foundOther = true; break;
                            }
                        }
                        if (foundOther) { break; }
                    }

                    // Add 2 new triangles and remove the old one

                    builder.Triangles.Add(builder.Triangles[i * 3 + (j + 2) % 3]);
                    builder.Triangles.Add(builder.Triangles[i * 3 + j]);
                    builder.Triangles.Add(builder.Triangles[vertIdx]);

                    builder.Triangles.Add(builder.Triangles[i * 3 + (j + 2) % 3]);
                    builder.Triangles.Add(builder.Triangles[vertIdx]);
                    builder.Triangles.Add(builder.Triangles[i * 3 + (j + 1) % 3]);

                    builder.Triangles.RemoveAt(i * 3);
                    builder.Triangles.RemoveAt(i * 3);
                    builder.Triangles.RemoveAt(i * 3);

                    return true;
                }

                // If the point is inside the triangle, add the point, remove the triangle and add 3 new triangles
                if (SMMath.IsPointInTriangle(pt2d, builder.FaceVertex2D(i, 0), builder.FaceVertex2D(i, 1), builder.FaceVertex2D(i, 2)))
                {
                    builder.Vertices.Add(pt);
                    builder.Vertices2D.Add(pt2d);
                    int vertIdx = builder.Vertices.Count - 1;

                    builder.Triangles.Add(builder.Triangles[i * 3 + 0]);
                    builder.Triangles.Add(builder.Triangles[i * 3 + 1]);
                    builder.Triangles.Add(vertIdx);

                    builder.Triangles.Add(builder.Triangles[i * 3 + 1]);
                    builder.Triangles.Add(builder.Triangles[i * 3 + 2]);
                    builder.Triangles.Add(vertIdx);

                    builder.Triangles.Add(builder.Triangles[i * 3 + 2]);
                    builder.Triangles.Add(builder.Triangles[i * 3 + 0]);
                    builder.Triangles.Add(vertIdx);

                    builder.Triangles.RemoveAt(i * 3);
                    builder.Triangles.RemoveAt(i * 3);
                    builder.Triangles.RemoveAt(i * 3);
                    return true;
                }
            }
            return false;
        }

        private static bool CutBuilderFaces(ref FaceBuilder builder, Vector3 ptA, Vector3 ptB, double tol, out int cutA, out int cutB)
        {
            cutA = -1;
            cutB = -1;

            var ptA2d = SMMath.TransformTo2D(builder.To2D, ptA);
            var ptB2d = SMMath.TransformTo2D(builder.To2D, ptB);

            var seg = new Segment2(ptA2d, ptB2d);
            var edges = new HashSet<Edge>();
            var segments = new List<Segment2>();

            for (int i = 0; i < builder.Triangles.Count; ++i)
            {
                for (int j = 0; j < 3; ++j)
                {
                    var indexA = builder.Triangles[i * 3 + j];
                    var indexB = builder.Triangles[i * 3 + (j + 1) % 3];
                    var edge = new Edge(Math.Min(indexA, indexB), Math.Max(indexA, indexB));

                    if (edges.Contains(edge)) {
                        continue;
                    }

                    edges.Add(edge);
                    segments.Add(new Segment2(builder.Vertices2D[edge.A], builder.Vertices2D[edge.B]));
                }
            }

            bool cut = false;

            int addedPt;
            Vector2 currCutA = ptB2d;
            Vector2 currCutB = ptA2d;

            for (int i = 0; i<segments.Count; ++i)
            {
                Vector2 intersection;
                if (SMMath.SegmentSegmentIntersection(seg, segments[i], out intersection))
                {
                    if (AddPointToBuilder(ref builder, SMMath.TransformTo3D(builder.To3D, intersection), tol, out addedPt)) {
                        if (SMMath.Vector2Distance(ptA2d, builder.Vertices2D[addedPt]) < SMMath.Vector2Distance(ptA2d, currCutA)) {
                            cutA = addedPt;
                            currCutA = builder.Vertices2D[cutA];
                        }
                        if (SMMath.Vector2Distance(ptB2d, builder.Vertices2D[addedPt]) < SMMath.Vector2Distance(ptB2d, currCutB))
                        {
                            cutB = addedPt;
                            currCutB = builder.Vertices2D[cutB];
                        }
                        cut = true;
                    }
                }
            }

            return cut;
        }

        private static void SplitBuilders(ref FaceBuilder builderA, ref FaceBuilder builderB, double tol) {
            // Points of intersection with the plane ( 0 to 3 )
            List<Vector3> AonB = new List<Vector3>();
            List<Vector3> BonA = new List<Vector3>();

            // If a point is on the plane there is no need to check also the two connected segment for intersection
            bool[] segA = new bool[3] { false, false, false };
            bool[] segB = new bool[3] { false, false, false };

            // Check if points of the face are on the plane of the other face
            for (int i = 0; i < 3; ++i)
            {
                if (SMMath.PointOnPlane(builderB.FacePlane, builderA.Vertices[i], tol))
                {
                    AonB.Add(builderA.Vertices[i]);
                    segA[i] = true;
                    segA[(i + 2) % 3] = true;
                }

                if (SMMath.PointOnPlane(builderA.FacePlane, builderB.Vertices[i], tol))
                {
                    BonA.Add(builderB.Vertices[i]);
                    segB[i] = true;
                    segB[(i + 2) % 3] = true;
                }
            }

            // Check for face's segments intersections with the other face's plane
            for (int i = 0; i < 3; ++i)
            {
                if (!segA[i])
                {
                    var pt = new Vector3();
                    var sg = new Segment3(builderA.Vertices[i], builderA.Vertices[(i + 1) % 3]);
                    if (SMMath.SegmentPlaneIntersection(sg, builderB.FacePlane, out pt))
                    {
                        AonB.Add(pt);
                        segA[i] = true;
                    }
                }

                if (!segB[i])
                {
                    var pt = new Vector3();
                    var sg = new Segment3(builderB.Vertices[i], builderB.Vertices[(i + 1) % 3]);
                    if (SMMath.SegmentPlaneIntersection(sg, builderA.FacePlane, out pt))
                    {
                        BonA.Add(pt);
                        segB[i] = true;
                    }
                }
            }

            if (BonA.Count > 0)
            {
                for (int p = 0; p < BonA.Count; ++p) {
                    int newPt;
                    AddPointToBuilder(ref builderA, BonA[p], tol, out newPt);
                }

                if (BonA.Count == 2)
                {
                    int cutA, cutB;
                    if (CutBuilderFaces(ref builderA, BonA[0], BonA[1], tol, out cutA, out cutB))
                    {
                        // Add cuts
                        builderA.AddCut(cutA, cutB, builderA.TestCut(cutA, cutB, builderB.FacePlane.Normal));
                    }
                    else if (BonA.Count == 3)
                    {
                        var a = builderA.ToLocal(BonA[0]);
                        var b = builderA.ToLocal(BonA[1]);
                        var c = builderA.ToLocal(BonA[2]);
                        var dirb = SMMath.Vector2Subtract(b, a);
                        var dirc = SMMath.Vector2Subtract(c, a);
                        var right = new Vector2(dirb.Y, -dirb.X);
                        bool invert = SMMath.Vector2Dot(right, dirc) < 0;

                        if (CutBuilderFaces(ref builderA, BonA[0], BonA[1], tol, out cutA, out cutB))
                        {
                            builderA.AddCut(cutA, cutB, invert);
                        }
                        if (CutBuilderFaces(ref builderA, BonA[1], BonA[2], tol, out cutA, out cutB))
                        {
                            builderA.AddCut(cutA, cutB, invert);
                        }
                        if (CutBuilderFaces(ref builderA, BonA[2], BonA[0], tol, out cutA, out cutB))
                        {
                            builderA.AddCut(cutA, cutB, invert);
                        }
                    }
                }
            }

            if (AonB.Count > 0)
            {
                for (int p = 0; p < AonB.Count; ++p)
                {
                    int newPt;
                    AddPointToBuilder(ref builderB, AonB[p], tol, out newPt);
                }

                int cutA, cutB;
                if (BonA.Count == 2)
                {
                    if (CutBuilderFaces(ref builderB, AonB[0], AonB[1], tol, out cutA, out cutB))
                    {
                        builderB.AddCut(cutA, cutB, builderB.TestCut(cutA, cutB, builderA.FacePlane.Normal));
                    }
                }
                else if (BonA.Count == 3)
                {
                    var a = builderB.ToLocal(AonB[0]);
                    var b = builderB.ToLocal(AonB[1]);
                    var c = builderB.ToLocal(AonB[2]);
                    var dirb = SMMath.Vector2Subtract(b, a);
                    var dirc = SMMath.Vector2Subtract(c, a);
                    var right = new Vector2(dirb.Y, -dirb.X);
                    bool invert = SMMath.Vector2Dot(right, dirc) < 0;

                    if (CutBuilderFaces(ref builderB, AonB[0], AonB[1], tol, out cutA, out cutB))
                    {
                        builderB.AddCut(cutA, cutB, invert);
                    }
                    if (CutBuilderFaces(ref builderB, AonB[1], AonB[2], tol, out cutA, out cutB))
                    {
                        builderB.AddCut(cutA, cutB, invert);
                    }
                    if (CutBuilderFaces(ref builderB, AonB[2], AonB[0], tol, out cutA, out cutB))
                    {
                        builderB.AddCut(cutA, cutB, invert);
                    }
                    // TODO: add cuts
                }
            }
        }

    }
}
