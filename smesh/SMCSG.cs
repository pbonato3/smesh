using SMesh;
using System;
using System.Collections.Generic;
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

            public List<int> cuts;

            public AABB BBox;
            public Plane FacePlane;

            public bool FromA;
            public bool Marked;
            public bool inside;

            public Matrix To2D;
            public Matrix To3D;
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

                    // Points of intersection with the plane ( 0 to 3 )
                    List<Vector3> AonB = new List<Vector3>();
                    List<Vector3> BonA = new List<Vector3>();

                    // If a point is on the plane there is no need to check also the two connected segment for intersection
                    bool[] segA = new bool[3] { false, false, false };
                    bool[] segB = new bool[3] { false, false, false };

                    // Check if points of the face are on the plane of the other face
                    for (int i = 0; i < 3; ++i) {
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
                            if (SMMath.SegmentPlaneIntersection(sg, builderB.FacePlane, out pt)) {
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

                    // TODO: Add points to face builders. Check against each sub face!
                    if (BonA.Count > 0)
                    {
                        for (int i = 0; i < builderA.Triangles.Count() / 3; ++i)
                        {
                            var va = builderA.Vertices[builderA.Triangles[i * 3 + 0]];
                            var vb = builderA.Vertices[builderA.Triangles[i * 3 + 1]];
                            var vc = builderA.Vertices[builderA.Triangles[i * 3 + 2]];

                        }
                    }

                    if (AonB.Count > 0)
                    {
                        for (int i = 0; i < builderB.Triangles.Count() / 3; ++i)
                        {
                            var va = builderB.Vertices[builderB.Triangles[i * 3 + 0]];
                            var vb = builderB.Vertices[builderB.Triangles[i * 3 + 1]];
                            var vc = builderB.Vertices[builderB.Triangles[i * 3 + 2]];

                        }
                    }

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



    }
}
