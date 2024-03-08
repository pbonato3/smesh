namespace SMesh
{
    public class SMCSG
    {
        public class Operation {
            public bool invertA;
            public bool invertB;
            public bool returnA;
            public bool returnB;
            public bool merge;

            public Operation(bool iA, bool iB, bool rA, bool rB, bool m) {
                invertA = iA;
                invertB = iB;
                returnA = rA;
                returnB = rB;
                merge   = m;
            }
            public static Operation SplitOperation() { 
                return new Operation(false, false, true, true, false);
            }

            public static Operation SplitAOperation()
            {
                return new Operation(false, false, true, false, false);
            }

            public static Operation SplitBOperation()
            {
                return new Operation(false, false, false, true, false);
            }

            public static Operation SubtractBfromAOperation()
            {
                return new Operation(false, true, true, true, true);
            }

            public static Operation SubtractAfromBOperation()
            {
                return new Operation(true, false, true, true, true);
            }

            public static Operation AddOperation()
            {
                return new Operation(false, false, true, true, true);
            }

            public static Operation IntersectOperation()
            {
                return new Operation(true, true, true, true, true);
            }
        }


        private class FaceBuilder {
            public List<Vector3> Vertices;
            public List<Vector3> Normals;

            public List<int> Triangles;

            // Cuts, definded as vertices indices
            public List<int> Cuts;                  // TODO: A set may improve performance when checking if cut already exists?
            // One cut flag for each vertex
            public List<bool> Cutted;

            public Sphere Sphere;
            public Plane FacePlane;

            public bool[] Marked;
            public bool[] Keep;

            public FaceBuilder() { 
                Vertices    = new List<Vector3>();
                Normals     = new List<Vector3>();
                Triangles   = new List<int>();
                Cuts        = new List<int>();
                Cutted      = new List<bool>();
                Marked      = new bool[1];
                Keep    = new bool[1];
            }

            public Vector3 FaceVertex(int nT, int nV) {
                return new Vector3(Vertices[Triangles[nT * 3 + nV]]);
            }

            public Vector3 FaceNormal(int nT, int nV)
            {
                return new Vector3(Normals[Triangles[nT * 3 + nV]]);
            }

            public bool TestCut(int idxA, int idxB, Plane testPlane) {
                var testIdx = 0;
                var maxDist = 0.0;
                var cutSign = true;
                for (int i = 0; i < 3; ++i) {
                    var sd = SMMath.PointPlaneSignedDistance(Vertices[i], testPlane);
                    var ud = Math.Abs(sd);
                    if (ud > maxDist) { 
                        testIdx = i;
                        maxDist = ud;
                        cutSign = sd > 0;
                    }
                }

                var vc = SMMath.Vector3Subtract(Vertices[idxB], Vertices[idxA]);
                var vt = SMMath.Vector3Subtract(Vertices[testIdx], Vertices[idxA]);
                var cross = SMMath.Vector3Cross(vt, vc);

                return cutSign ? SMMath.Vector3Dot(cross, FacePlane.Normal) > 0 : SMMath.Vector3Dot(cross, FacePlane.Normal) < 0;

                /*
                for (int i = Triangles.Count - 3; i >= 0; i-=3) { 
                    var a = Triangles[i + 0];
                    var b = Triangles[i + 1];
                    var c = Triangles[i + 2];

                    if (idxA == a && idxB == b) { return SMMath.PointPlaneSignedDistance(Vertices[c], testPlane) > 0; }
                    if (idxA == b && idxB == c) { return SMMath.PointPlaneSignedDistance(Vertices[a], testPlane) > 0; }
                    if (idxA == c && idxB == a) { return SMMath.PointPlaneSignedDistance(Vertices[b], testPlane) > 0; }

                    if (idxA == a && idxB == c) { return SMMath.PointPlaneSignedDistance(Vertices[b], testPlane) < 0; }
                    if (idxA == b && idxB == a) { return SMMath.PointPlaneSignedDistance(Vertices[c], testPlane) < 0; }
                    if (idxA == c && idxB == b) { return SMMath.PointPlaneSignedDistance(Vertices[a], testPlane) < 0; }
                }
                Console.WriteLine("ERROR - No face found to test cut direction.");
                return true;
                */
            }

            public bool AddCut(int idxA, int idxB, bool invert) {
                if (idxA == idxB) { 
                    return false;
                }

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

        //TODO: Use this method to force one mark in separated mesh components
        private static bool IsPointInsideMesh(Vector3 pt, SphereBVH.SphereNode tree, FaceBuilder[] builders, double tol = 0.00001) {
            if (SMMath.Vector3Distance(pt, tree.Sphere.Point) > tree.Sphere.Radius) {
                return false;
            }
            var rayOrigin = new Vector3(pt);
            rayOrigin.Z = rayOrigin.Z - tree.Sphere.Radius - 1;
            var rayDirection = SMMath.Vector3Subtract(pt, rayOrigin);
            Segment3 sg = new Segment3(rayOrigin, pt);

            List<int> collisions = tree.Search(rayOrigin, rayDirection);

            int intersectionCount = 0;
            for (int i = 0; i < collisions.Count; i++)
            {
                Vector3 intersection;
                var fb = builders[collisions[i]];
                if (SMMath.SegmentPlaneIntersection(sg, fb.FacePlane, out intersection, tol)) {
                    if (SMMath.IsPointInTriangle(intersection, fb.Vertices[0], fb.Vertices[1], fb.Vertices[2], tol))
                    {
                        intersectionCount++;
                    }
                }
            }
            return intersectionCount % 2 == 1;
        }


        // TODO: is this the right place for this 3 methods?
        public static List<int> [] GetVerticesConnectivity(Mesh mesh) {
            List<int>[] vConn = new List<int> [mesh.VertCount];

            for (int i = 0; i < mesh.FaceCount; i++) {
                for (int n = 0; n < 3; ++n) { 
                    int idx = mesh.Indices[i * 3 + n];
                    if (vConn[idx] == null) {
                        vConn[idx] = new List<int>();
                    }
                    vConn[idx].Add(i);
                }
            }

            return vConn;
        }

        // TODO: Use this method in Mark Expansion to limit the number of builders to explore
        public static HashSet<int>[] GetFacesConnectivity(Mesh mesh) {

            HashSet<int>[] fConn = new HashSet<int>[mesh.FaceCount];
            var vConn = GetVerticesConnectivity(mesh);

            for (int i = 0; i < mesh.FaceCount; i++)
            {
                for (int n = 0; n < 3; ++n)
                {
                    int idx = mesh.Indices[i * 3 + n];
                    var faces = vConn[idx];

                    for (int f = 0; f < faces.Count; f++) {
                        if (faces[f] == i) {
                            continue;
                        }

                        if (fConn[i] == null) { fConn[i] = new HashSet<int>(); }
                        fConn[i].Add(faces[f]);
                    }
                }
            }

            return fConn;
        }

        public List<HashSet<int>> GetMeshComponents(Mesh mesh) { 
            var connVerts = new List<HashSet<int>>();
            var connFaces = new List<HashSet<int>>();

            for (int i = 0; i < mesh.FaceCount; ++i)
            {
                int a = mesh.Indices[i * 3 + 0];
                int b = mesh.Indices[i * 3 + 1];
                int c = mesh.Indices[i * 3 + 2];

                List<int> matches = new List<int>();
                for (int j = 0; j < connVerts.Count; ++j)
                {
                    if (connVerts[j].Contains(a) || connVerts[j].Contains(b) || connVerts[j].Contains(c))
                    {
                        matches.Add(j);
                        connVerts[j].Add(a);
                        connVerts[j].Add(b);
                        connVerts[j].Add(c);
                        connFaces[j].Add(i);
                    }
                }

                matches.Sort();

                // A new component
                if (matches.Count == 0)
                {
                    connVerts.Add(new HashSet<int>());
                    connFaces.Add(new HashSet<int>());
                    matches.Add(connVerts.Count - 1);
                }

                // Merge components
                for (int n = matches.Count - 1; n > 1; n--)
                {
                    foreach (var v in connVerts[matches[n]])
                    {
                        connVerts[matches[0]].Add(v);
                    }
                    foreach (var v in connFaces[matches[n]])
                    {
                        connFaces[matches[0]].Add(v);
                    }
                    connVerts.RemoveAt(matches[n]);
                    connFaces.RemoveAt(matches[n]);
                }
            }

            return connFaces;
        }


        // DEPRECATED
        public static Mesh[] Split(Mesh meshA, Mesh meshB, double tol = 0.00001) { 
            FaceBuilder[] buildersA = new FaceBuilder[meshA.FaceCount];
            FaceBuilder[] buildersB = new FaceBuilder[meshB.FaceCount];

            // init builders for mesh A
            for (int f = 0; f < meshA.FaceCount; f++)
            {
                buildersA[f] = new FaceBuilder();

                var a = meshA.Vertices[meshA.Indices[f * 3 + 0]];
                var b = meshA.Vertices[meshA.Indices[f * 3 + 1]];
                var c = meshA.Vertices[meshA.Indices[f * 3 + 2]];

                buildersA[f].Vertices = new List<Vector3>(3) { a, b, c };
                buildersA[f].Triangles = new List<int>(3) { 0, 1, 2 };

                if (meshA.Normals != null && meshA.Normals.Length == meshA.Vertices.Length)
                {
                    var na = meshA.Normals[meshA.Indices[f * 3 + 0]];
                    var nb = meshA.Normals[meshA.Indices[f * 3 + 1]];
                    var nc = meshA.Normals[meshA.Indices[f * 3 + 2]];
                    buildersA[f].Normals = new List<Vector3>(3) { na, nb, nc };
                }


                buildersA[f].Sphere = ThreePointsSphere(a, b, c);
            }

            // init builders for mesh B
            for (int f = 0; f < meshB.FaceCount; f++)
            {
                buildersB[f] = new FaceBuilder();

                var a = meshB.Vertices[meshB.Indices[f * 3 + 0]];
                var b = meshB.Vertices[meshB.Indices[f * 3 + 1]];
                var c = meshB.Vertices[meshB.Indices[f * 3 + 2]];

                buildersB[f].Vertices = new List<Vector3>(3) { a, b, c };
                buildersB[f].Triangles = new List<int>(3) { 0, 1, 2 };

                if (meshB.Normals != null && meshB.Normals.Length == meshB.Vertices.Length) { 
                    var na = meshB.Normals[meshB.Indices[f * 3 + 0]];
                    var nb = meshB.Normals[meshB.Indices[f * 3 + 1]];
                    var nc = meshB.Normals[meshB.Indices[f * 3 + 2]];
                    buildersB[f].Normals = new List<Vector3>(3) { na, nb, nc };
                }

                buildersB[f].Sphere = ThreePointsSphere(a, b, c);
            }

            //Build bvh for mesh B
            var bvhB = SphereBVH.BuildBVH(meshB);

            // Check every face of A agains B's BVH
            for (int f = 0; f < meshA.FaceCount; f++)
            {
                List<int> collisions = bvhB.Search(buildersA[f].Sphere);

                var builderA = buildersA[f];
                for (int c = 0; c < collisions.Count; c++) {
                    var builderB = buildersB[collisions[c]];

                    // TODO: Do not recompute plane
                    // Note: in Vertices 0, 1, 2 there will always be the original triangle
                    builderA.FacePlane = SMMath.PlaneFrom3Points(builderA.Vertices[0], builderA.Vertices[1], builderA.Vertices[2]);
                    builderB.FacePlane = SMMath.PlaneFrom3Points(builderB.Vertices[0], builderB.Vertices[1], builderB.Vertices[2]);

                    SplitBuilders(ref builderA, ref builderB, tol);
                }
            }

            List<FaceBuilder> outbuild = new List<FaceBuilder>();
            for (int i = 0; i < buildersA.Length; ++i) {
                MarkBuilderFaces(ref buildersA[i], true);
                outbuild.Add(buildersA[i]);
            }
            for (int i = 0; i < buildersB.Length; ++i)
            {
                MarkBuilderFaces(ref buildersB[i], true);
                outbuild.Add(buildersB[i]);
            }

            ExpandMarks(buildersA);
            ExpandMarks(buildersB);

            var result = new List<Mesh>();
            result.Add(BuildIndexedMesh(outbuild.ToArray()));
            //result.Add(BuildIndexedMesh(buildersB.ToArray()));
            return result.ToArray();
        }

        public static Mesh[] RunOperation(Mesh meshA, Mesh meshB, Operation op, double tol = 0.00001)
        {
            FaceBuilder[] buildersA = new FaceBuilder[meshA.FaceCount];
            FaceBuilder[] buildersB = new FaceBuilder[meshB.FaceCount];

            // init builders for mesh A
            for (int f = 0; f < meshA.FaceCount; f++)
            {
                buildersA[f] = new FaceBuilder();

                var a = meshA.Vertices[meshA.Indices[f * 3 + 0]];
                var b = meshA.Vertices[meshA.Indices[f * 3 + 1]];
                var c = meshA.Vertices[meshA.Indices[f * 3 + 2]];

                buildersA[f].Vertices = new List<Vector3>(3) { a, b, c };
                buildersA[f].Triangles = new List<int>(3) { 0, 1, 2 };

                if (meshA.Normals != null && meshA.Normals.Length == meshA.Vertices.Length)
                {
                    var na = meshA.Normals[meshA.Indices[f * 3 + 0]];
                    var nb = meshA.Normals[meshA.Indices[f * 3 + 1]];
                    var nc = meshA.Normals[meshA.Indices[f * 3 + 2]];
                    buildersA[f].Normals = new List<Vector3>(3) { na, nb, nc };
                }


                buildersA[f].Sphere = ThreePointsSphere(a, b, c);
            }

            // init builders for mesh B
            for (int f = 0; f < meshB.FaceCount; f++)
            {
                buildersB[f] = new FaceBuilder();

                var a = meshB.Vertices[meshB.Indices[f * 3 + 0]];
                var b = meshB.Vertices[meshB.Indices[f * 3 + 1]];
                var c = meshB.Vertices[meshB.Indices[f * 3 + 2]];

                buildersB[f].Vertices = new List<Vector3>(3) { a, b, c };
                buildersB[f].Triangles = new List<int>(3) { 0, 1, 2 };

                if (meshB.Normals != null && meshB.Normals.Length == meshB.Vertices.Length)
                {
                    var na = meshB.Normals[meshB.Indices[f * 3 + 0]];
                    var nb = meshB.Normals[meshB.Indices[f * 3 + 1]];
                    var nc = meshB.Normals[meshB.Indices[f * 3 + 2]];
                    buildersB[f].Normals = new List<Vector3>(3) { na, nb, nc };
                }

                buildersB[f].Sphere = ThreePointsSphere(a, b, c);
            }

            //Build bvh for mesh B
            var bvhB = SphereBVH.BuildBVH(meshB);

            // Check every face of A agains B's BVH
            for (int f = 0; f < meshA.FaceCount; f++)
            {
                List<int> collisions = bvhB.Search(buildersA[f].Sphere);

                var builderA = buildersA[f];
                for (int c = 0; c < collisions.Count; c++)
                {
                    var builderB = buildersB[collisions[c]];

                    // TODO: Do not recompute plane
                    // Note: in Vertices 0, 1, 2 there will always be the original triangle
                    builderA.FacePlane = SMMath.PlaneFrom3Points(builderA.Vertices[0], builderA.Vertices[1], builderA.Vertices[2]);
                    builderB.FacePlane = SMMath.PlaneFrom3Points(builderB.Vertices[0], builderB.Vertices[1], builderB.Vertices[2]);

                    SplitBuilders(ref builderA, ref builderB, tol);
                }
            }


            for (int i = 0; i < buildersA.Length; ++i)
            {
                MarkBuilderFaces(ref buildersA[i], op.invertA);
            }
            for (int i = 0; i < buildersB.Length; ++i)
            {
                MarkBuilderFaces(ref buildersB[i], op.invertB);
            }

            var fConnA = GetFacesConnectivity(meshA);
            var fConnB = GetFacesConnectivity(meshB);

            if (!ExpandMarks2(buildersA, fConnA)) {
                for (int b = 0; b < buildersA.Length; ++b) {
                    if (buildersA[b].Marked[0]) { 
                        continue;
                    }
                    var keep = !(IsPointInsideMesh(buildersA[b].Vertices[0], bvhB.Root, buildersB, tol) || op.invertA);
                    var tomark = new Queue<int>();
                    tomark.Enqueue(b);

                    while (tomark.Count > 0) {
                        var index = tomark.Dequeue();
                        if (buildersA[index].Marked[0]) { continue; }
                        for (int t = 0; t < buildersA[index].Triangles.Count / 3; ++t)
                        {
                            buildersA[index].Marked[t] = true;
                            buildersA[index].Keep[t] = keep;
                        }
                        foreach (var c in fConnA[b]) {
                            tomark.Enqueue(c);
                        }
                    }

                }
            }
            if (!ExpandMarks2(buildersB, fConnB))
            {
                var bvhA = SphereBVH.BuildBVH(meshA);
                for (int b = 0; b < buildersB.Length; ++b)
                {
                    if (buildersB[b].Marked[0])
                    {
                        continue;
                    }
                    var keep = !(IsPointInsideMesh(buildersB[b].Vertices[0], bvhA.Root, buildersA, tol) || op.invertB);
                    var tomark = new Queue<int>();
                    tomark.Enqueue(b);

                    while (tomark.Count > 0)
                    {
                        var index = tomark.Dequeue();
                        if (buildersB[index].Marked[0]) { continue; }
                        for (int t = 0; t < buildersB[index].Triangles.Count / 3; ++t)
                        {
                            buildersB[index].Marked[t] = true;
                            buildersB[index].Keep[t] = keep;
                        }
                        foreach (var c in fConnB[b])
                        {
                            tomark.Enqueue(c);
                        }
                    }
                }

            }

            var result = new List<Mesh>();

            if (op.merge)
            {
                List<FaceBuilder> outbuild = new List<FaceBuilder>();
                for (int i = 0; i < buildersA.Length; ++i)
                {
                    outbuild.Add(buildersA[i]);
                }
                for (int i = 0; i < buildersB.Length; ++i)
                {
                    outbuild.Add(buildersB[i]);
                }
                result.Add(BuildIndexedMesh(outbuild.ToArray()));
            }
            else {
                if (op.returnA) { 
                    result.Add(BuildIndexedMesh(buildersA.ToArray()));
                    result.Add(BuildIndexedMesh(buildersA.ToArray(), true));
                }
                if (op.returnB) { 
                    result.Add(BuildIndexedMesh(buildersB.ToArray()));
                    result.Add(BuildIndexedMesh(buildersB.ToArray(), true));
                }
            }

            return result.ToArray();
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

        private static Sphere ThreePointsSphere(Vector3 a, Vector3 b, Vector3 c)
        {
            var sphere = new Sphere();
            sphere.Point = SMMath.Vector3Divide(SMMath.Vector3Add(SMMath.Vector3Add(a, b), c), 3);
            sphere.Radius = Math.Max(Math.Max(SMMath.Vector3Distance(a, sphere.Point), SMMath.Vector3Distance(b, sphere.Point)), SMMath.Vector3Distance(c, sphere.Point));
            return sphere;
        }


        private static bool AddPointToBuilder(ref FaceBuilder builder, Vector3 pt, double tol, out int newPt){

            // If point doesn't snap to an existing one, it will be always added as last one.
            newPt = builder.Vertices.Count;

            // If vertex is colose enough to an existing one do nothing
            for (int i = 0; i < builder.Vertices.Count; ++i) {
                if (SMMath.Vector3Distance(builder.Vertices[i], pt) <= tol) {
                    newPt = i;
                    return true;
                }
            }

            for (int i = 0; i < builder.Triangles.Count / 3; ++i)
            {
                for (int j = 0; j < 3; ++j)
                {
                    var closest = SMMath.SegmentClosestPoint(pt, new Segment3(builder.FaceVertex(i, j), builder.FaceVertex(i, (j + 1) % 3)));
                    if (SMMath.Vector3Distance(closest, pt) > tol)
                    {
                        continue;
                    }

                    // TODO: add pt or closest?
                    builder.Vertices.Add(closest);
                    if (builder.Normals != null && builder.Normals.Count > 0)
                    {
                        builder.Normals.Add(SampleNormal(builder, builder.Vertices[builder.Vertices.Count - 1]));
                    }
                    int vertIdx = builder.Vertices.Count - 1;

                    var indexA = builder.Triangles[i * 3 + j];
                    var indexB = builder.Triangles[i * 3 + (j + 1) % 3];
                    var indexC = builder.Triangles[i * 3 + (j + 2) % 3];

                    // Add 2 new triangles and remove the old one

                    builder.Triangles.Add(indexC);
                    builder.Triangles.Add(indexA);
                    builder.Triangles.Add(vertIdx);

                    builder.Triangles.Add(indexC);
                    builder.Triangles.Add(vertIdx);
                    builder.Triangles.Add(indexB);
   
                    builder.Triangles.RemoveAt(i * 3);
                    builder.Triangles.RemoveAt(i * 3);
                    builder.Triangles.RemoveAt(i * 3);

                    // Do the same thing with the opposite triangle on that segment if exists
                    for (int ii = i; ii < builder.Triangles.Count / 3 - 2; ++ii)
                    {
                        for (int jj = 0; jj < 3; ++jj)
                        {
                            if (indexB == builder.Triangles[ii * 3 + jj] && indexA == builder.Triangles[ii * 3 + (jj + 1) % 3]) {
                                var indexOpposite = builder.Triangles[ii * 3 + (jj + 2) % 3];

                                builder.Triangles.Add(indexOpposite);
                                builder.Triangles.Add(indexB);
                                builder.Triangles.Add(vertIdx);

                                builder.Triangles.Add(indexOpposite);
                                builder.Triangles.Add(vertIdx);
                                builder.Triangles.Add(indexA);

                                builder.Triangles.RemoveAt(ii * 3);
                                builder.Triangles.RemoveAt(ii * 3);
                                builder.Triangles.RemoveAt(ii * 3);

                                return true;
                            }
                        }
                    }
                    return true;
                }

                // If the point is inside the triangle, add the point, remove the triangle and add 3 new triangles
                if (SMMath.IsPointInTriangle(pt, builder.FaceVertex(i, 0), builder.FaceVertex(i, 1), builder.FaceVertex(i, 2), tol))
                {
                    builder.Vertices.Add(pt);
                    if (builder.Normals != null && builder.Normals.Count > 0) { 
                        builder.Normals.Add(SampleNormal(builder, pt));
                    }
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

        private static bool CutBuilderFaces(ref FaceBuilder builder, Vector3 ptA, Vector3 ptB, double tol, out List<int> cuts)
        {
            cuts = new List<int>();

            var seg = new Segment3(ptA, ptB);
            var edges = new HashSet<Edge>();
            var segments = new List<Segment3>();

            for (int i = 0; i < builder.Triangles.Count / 3; ++i)
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
                    segments.Add(new Segment3(builder.Vertices[edge.A], builder.Vertices[edge.B]));
                }
            }

            bool cut = false;

            int addedPt;

            for (int i = 0; i<segments.Count; ++i)
            {
                Vector3 intersection;
                if (SMMath.SegmentSegmentIntersection(seg, segments[i], out intersection, tol))
                {
                    if (AddPointToBuilder(ref builder, intersection, tol, out addedPt))
                    {
                        var newDist = SMMath.Vector3Distance(ptA, builder.Vertices[addedPt]);

                        for (int c = 0; c <= cuts.Count; ++c)
                        {
                            if (c == cuts.Count)
                            {
                                cuts.Add(addedPt);
                                break;
                            }

                            var d = SMMath.Vector3Distance(ptA, builder.Vertices[cuts[c]]);
                            if (newDist <= d)
                            {
                                cuts.Insert(c, addedPt);
                                break;
                            }
                        }
                        cut = true;
                    }
                    else {
                        Console.WriteLine("ERROR HERE!");
                        //TODO: solve numeric issues
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
                    if (SMMath.SegmentPlaneIntersection(sg, builderB.FacePlane, out pt, tol))
                    {
                        AonB.Add(pt);
                        segA[i] = true;
                    }
                }

                if (!segB[i])
                {
                    var pt = new Vector3();
                    var sg = new Segment3(builderB.Vertices[i], builderB.Vertices[(i + 1) % 3]);
                    if (SMMath.SegmentPlaneIntersection(sg, builderA.FacePlane, out pt, tol))
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

                List<int> cuts;
                if (BonA.Count == 2)
                {
                    if (CutBuilderFaces(ref builderA, BonA[0], BonA[1], tol, out cuts))
                    {
                        // Add cuts
                        for (int k = 0; k < cuts.Count - 1; ++k)
                        {
                            if (cuts[k] != cuts[k + 1])
                            {
                                builderA.AddCut(cuts[k], cuts[k + 1], builderA.TestCut(cuts[k], cuts[k + 1], builderB.FacePlane));
                            }
                        }
                    }
                }
                else if (BonA.Count == 3)
                {
                    bool invert = SMMath.Vector3Dot(builderA.FacePlane.Normal, builderB.FacePlane.Normal) < 0;

                    if (CutBuilderFaces(ref builderA, BonA[0], BonA[1], tol, out cuts))
                    {
                        for (int k = 0; k < cuts.Count - 1; ++k)
                        {
                            builderA.AddCut(cuts[k], cuts[k + 1], invert);
                        }
                    }
                    if (CutBuilderFaces(ref builderA, BonA[1], BonA[2], tol, out cuts))
                    {
                        for (int k = 0; k < cuts.Count - 1; ++k)
                        {
                            builderA.AddCut(cuts[k], cuts[k + 1], invert);
                        }
                    }
                    if (CutBuilderFaces(ref builderA, BonA[2], BonA[0], tol, out cuts))
                    {
                        for (int k = 0; k < cuts.Count - 1; ++k)
                        {
                            builderA.AddCut(cuts[k], cuts[k + 1], invert);
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

                List<int> cuts;
                if (AonB.Count == 2)
                {
                    if (CutBuilderFaces(ref builderB, AonB[0], AonB[1], tol, out cuts))
                    {
                        for (int k = 0; k < cuts.Count - 1; ++k)
                        {
                            if (cuts[k] != cuts[k + 1])
                            {
                                builderB.AddCut(cuts[k], cuts[k + 1], builderB.TestCut(cuts[k], cuts[k + 1], builderA.FacePlane));
                            }
                        }
                    }
                }
                else if (AonB.Count == 3)
                {
                    bool invert = SMMath.Vector3Dot(builderA.FacePlane.Normal, builderB.FacePlane.Normal) < 0;

                    if (CutBuilderFaces(ref builderB, AonB[0], AonB[1], tol, out cuts))
                    {
                        for (int k = 0; k < cuts.Count - 1; ++k)
                        {
                            builderB.AddCut(cuts[k], cuts[k + 1], invert);
                        }
                    }
                    if (CutBuilderFaces(ref builderB, AonB[1], AonB[2], tol, out cuts))
                    {
                        for (int k = 0; k < cuts.Count - 1; ++k)
                        {
                            builderB.AddCut(cuts[k], cuts[k + 1], invert);
                        }
                    }
                    if (CutBuilderFaces(ref builderB, AonB[2], AonB[0], tol, out cuts))
                    {
                        for (int k = 0; k < cuts.Count - 1; ++k)
                        {
                            builderB.AddCut(cuts[k], cuts[k + 1], invert);
                        }
                    }
                }
            }
        }


        private static void MarkBuilderFaces(ref FaceBuilder builder, bool invert) {
            // TODO: maybe this initialization should be better outside this method
            builder.Marked = new bool[builder.Triangles.Count / 3];
            builder.Keep = new bool[builder.Triangles.Count / 3];
            for (int i = 0; i < builder.Marked.Count(); ++i) {
                builder.Marked[i] = false;
                builder.Keep[i] = false; // TODO: initialization not necessary, only debug???
            }

            // If no cuts, no marks possible
            if (builder.Cuts.Count == 0)
            {
                return;
            }

            HashSet<(int, int)> keepEdges = new HashSet<(int, int)> ();
            HashSet<(int, int)> skipEdges = new HashSet<(int, int)> ();

            // Mark triangles using cuts
            for (int i = 0; i < builder.Triangles.Count / 3; ++i)
            {
                var a = builder.Triangles[i * 3 + 0];
                var b = builder.Triangles[i * 3 + 1];
                var c = builder.Triangles[i * 3 + 2];

                bool keep = false;
                bool skip = false;
                for (int j = 0; j < builder.Cuts.Count / 2; ++j)
                {
                    var ca = builder.Cuts[j * 2 + 0];
                    var cb = builder.Cuts[j * 2 + 1];

                    if ((ca == a && cb == b) || (ca == b && cb == c) || (ca == c && cb == a))
                    {
                        keep = true;
                    }
                    if ((ca == a && cb == c) || (ca == c && cb == b) || (ca == b && cb == a))
                    {
                        skip = true;
                    }

                    // TODO: this is for debug only
                    if (keep && skip) { 
                        continue;
                    }
                }

                if (invert) {
                    var swap = keep;
                    keep = skip;
                    skip = swap;
                }

                if (keep) {
                    builder.Marked[i] = true;
                    builder.Keep[i] = true;
                    keepEdges.Add((a, c));
                    keepEdges.Add((b, a));
                    keepEdges.Add((c, b));
                }
                if (skip) {
                    builder.Marked[i] = true;
                    builder.Keep[i] = false;
                    skipEdges.Add((a, c));
                    skipEdges.Add((b, a));
                    skipEdges.Add((c, b));
                }
            }


            // It should be always possible to expand marks inside the same builder if there is at least one marked face
            int markedCount = keepEdges.Count + skipEdges.Count;
            bool repeat = markedCount > 0 && markedCount < builder.Triangles.Count;
            while (repeat) {
                for (int i = 0; i < builder.Triangles.Count / 3; ++i) {
                    if (builder.Marked[i]) {
                        continue;
                    }
                    var a = builder.Triangles[i * 3 + 0];
                    var b = builder.Triangles[i * 3 + 1];
                    var c = builder.Triangles[i * 3 + 2];

                    if (keepEdges.Contains((a, b)) || keepEdges.Contains((b, c)) || keepEdges.Contains((c, a)))
                    {
                        builder.Marked[i] = true;
                        builder.Keep[i] = true;
                        keepEdges.Add((a, c));
                        keepEdges.Add((b, a));
                        keepEdges.Add((c, b));
                    }
                    else if (skipEdges.Contains((a, b)) || skipEdges.Contains((b, c)) || skipEdges.Contains((c, a)))
                    {
                        builder.Marked[i] = true;
                        builder.Keep[i] = false;
                        skipEdges.Add((a, c));
                        skipEdges.Add((b, a));
                        skipEdges.Add((c, b));
                    }
                }
                markedCount = keepEdges.Count + skipEdges.Count;
                repeat = markedCount > 0 && markedCount < builder.Triangles.Count;
            }
        }

        private static void ExpandMarks(FaceBuilder[] builders)
        {
            // TODO: Try to describe connectivity of the mesh and expand marks in a more efficient way

            HashSet<(Vector3, Vector3)> keepEdges = new HashSet<(Vector3, Vector3)>();
            HashSet<(Vector3, Vector3)> skipEdges = new HashSet<(Vector3, Vector3)>();

            foreach (var b in builders) {
                for (var t = 0; t < b.Triangles.Count / 3; ++t) {
                    if ( b.Marked[t] ) {
                        if (b.Keep[t])
                        {
                            keepEdges.Add((b.FaceVertex(t, 0), b.FaceVertex(t, 2)));
                            keepEdges.Add((b.FaceVertex(t, 1), b.FaceVertex(t, 0)));
                            keepEdges.Add((b.FaceVertex(t, 2), b.FaceVertex(t, 1)));
                        }
                        else {
                            skipEdges.Add((b.FaceVertex(t, 0), b.FaceVertex(t, 2)));
                            skipEdges.Add((b.FaceVertex(t, 1), b.FaceVertex(t, 0)));
                            skipEdges.Add((b.FaceVertex(t, 2), b.FaceVertex(t, 1)));
                        }
                    }
                }
            }

            var missingMarks = false;
            var somethingChanged = true;
            while (somethingChanged)
            {
                missingMarks = false;
                somethingChanged = false;
                foreach (var build in builders)
                {
                    for (var t = 0; t < build.Triangles.Count / 3; ++t)
                    {
                        if (build.Marked[t])
                        {
                            continue;
                        }


                        Vector3[] rounded = new Vector3[3]{ build.FaceVertex(t, 0), build.FaceVertex(t, 1), build.FaceVertex(t, 2)};

                        var a = build.FaceVertex(t, 0);
                        var b = build.FaceVertex(t, 1);
                        var c = build.FaceVertex(t, 2);

                        var countKeep = 0;
                        var countSkip = 0;

                        // TODO: this count is probably unusefull
                        if (keepEdges.Contains((a, b))) { countKeep++; }
                        if (keepEdges.Contains((b, c))) { countKeep++; }
                        if (keepEdges.Contains((c, a))) { countKeep++; }
                        if (skipEdges.Contains((a, b))) { countSkip++; }
                        if (skipEdges.Contains((b, c))) { countSkip++; }
                        if (skipEdges.Contains((c, a))) { countSkip++; }

                        if (countKeep > countSkip)
                        {
                            build.Marked[t] = true;
                            build.Keep[t] = true;
                            somethingChanged = true;

                            keepEdges.Add((a, c));
                            keepEdges.Add((b, a));
                            keepEdges.Add((c, b));
                            continue;
                        }
                        else if (countKeep < countSkip)
                        {
                            build.Marked[t] = true;
                            build.Keep[t] = false;
                            somethingChanged = true;

                            skipEdges.Add((a, c));
                            skipEdges.Add((b, a));
                            skipEdges.Add((c, b));
                            continue;
                        }

                        if (countKeep > 0) {
                            Console.WriteLine("MANNAGGIA");
                        }
                        missingMarks = true;
                    }
                }
            }
            Console.WriteLine(missingMarks);
        }

        private static bool ExpandMarks2(FaceBuilder[] builders, HashSet<int>[] fConn)
        {
            // We suppose that if a builder's face is marked, then all other faces are

            HashSet<int> marked = new HashSet<int>();
            HashSet<int> tomark = new HashSet<int>();

            for (int i = 0; i < builders.Length; i++) {
                if (builders[i].Marked[0]) {
                    marked.Add(i);
                    tomark.Remove(i);
                    foreach (var conn in fConn[i]) {
                        if (!marked.Contains(conn)) { 
                            tomark.Add(conn);
                        }
                    }
                }
            }

            var markCount = 0;
            while (marked.Count > markCount) {
                markCount = marked.Count;
                var newFrontier = new HashSet<int>();

                foreach (var index in tomark) {
                    var builder = builders[index];
                    var connections = fConn[index];
                    HashSet<(Vector3, Vector3)> edges = new HashSet<(Vector3, Vector3)>();
                    for (int t = 0; t < builder.Triangles.Count / 3; ++t)
                    {
                        edges.Add((builder.FaceVertex(t, 0), builder.FaceVertex(t, 1)));
                        edges.Add((builder.FaceVertex(t, 1), builder.FaceVertex(t, 2)));
                        edges.Add((builder.FaceVertex(t, 2), builder.FaceVertex(t, 0)));
                    }

                    foreach (var conn in connections) {
                        if (!marked.Contains(conn)) {
                            continue;
                        }

                        var cbuilder = builders[conn];
                        for (int t = 0; t < cbuilder.Triangles.Count / 3; ++t) {
                            if (edges.Contains((cbuilder.FaceVertex(t, 0), cbuilder.FaceVertex(t, 2))) ||
                                edges.Contains((cbuilder.FaceVertex(t, 1), cbuilder.FaceVertex(t, 0))) ||
                                edges.Contains((cbuilder.FaceVertex(t, 2), cbuilder.FaceVertex(t, 1)))) {
                                for (int k = 0; k < builder.Triangles.Count / 3; ++k) {
                                    builder.Marked[k] = true;
                                    builder.Keep[k] = cbuilder.Keep[t];
                                }
                                break;
                            }
                        }

                        if (builder.Marked[0])
                        {
                            marked.Add(index);
                            foreach (var c in fConn[index])
                            {
                                if (!marked.Contains(c))
                                {
                                    newFrontier.Add(c);
                                }
                            }
                            break;
                        }
                        else {
                            newFrontier.Add(index);
                        }

                    }
                    tomark = newFrontier;
                }
            }

            return marked.Count == builders.Length;
        }


        private static Mesh BuildUnindexedMesh(FaceBuilder[] builders)
        {
            Mesh mesh = new Mesh();
            var vertices = new List<Vector3>();
            var indices = new List<int>();

            foreach (var builder in builders) {
                for (int i = 0; i < builder.Triangles.Count / 3; ++i) {
                    for (int j = 0; j < 3; ++j)
                    {
                        indices.Add(indices.Count);
                        vertices.Add(builder.FaceVertex(i, j));
                    }
                }
            }

            mesh.VertCount = vertices.Count;
            mesh.Vertices = vertices.ToArray();
            mesh.FaceCount = indices.Count / 3;
            mesh.Indices = indices.ToArray();

            RebuildNormals(ref mesh);

            return mesh;
        }
        

        private static Mesh BuildIndexedMesh(FaceBuilder[] builders, bool negate = false)
        {
            Mesh mesh = new Mesh();
            var vertices = new List<Vector3>();
            var normals = new List<Vector3>();
            var indices = new List<int>();
            var weldMap = new Dictionary<Vector3, List<int>>(); // TODO: use a bvh?

            foreach (var builder in builders)
            {
                for (int i = 0; i < builder.Triangles.Count / 3; ++i)
                {
                    if (builder.Keep[i] != negate) {
                        continue;
                    }
                    for (int j = 0; j < 3; ++j)
                    {
                        if (weldMap.ContainsKey(builder.FaceVertex(i, j)))
                        {
                            bool weld = false;
                            Vector3 normal = builder.FaceNormal(i, j);
                            for (int n = 0; n < weldMap[builder.FaceVertex(i, j)].Count; ++n)
                            {
                                int index = weldMap[builder.FaceVertex(i, j)][n];
                                if (SMMath.AreEquals(normal, normals[index], 0)) {
                                    indices.Add(index);
                                    weld = true;
                                    break;
                                }
                            }
                            if (!weld) {
                                indices.Add(vertices.Count);
                                weldMap[builder.FaceVertex(i, j)].Add(vertices.Count);
                                vertices.Add(builder.FaceVertex(i, j));
                                normals.Add(builder.FaceNormal(i, j));
                            }
                        }
                        else {
                            indices.Add(vertices.Count);
                            weldMap.Add(builder.FaceVertex(i, j), new List<int> { vertices.Count });
                            vertices.Add(builder.FaceVertex(i, j));
                            normals.Add(builder.FaceNormal(i, j));
                        }
                    }
                }
            }

            mesh.VertCount = vertices.Count;
            mesh.Vertices = vertices.ToArray();
            mesh.Normals = normals.ToArray();
            mesh.FaceCount = indices.Count / 3;
            mesh.Indices = indices.ToArray();

            return mesh;
        }

        public static void RebuildNormals(ref Mesh mesh) {
            mesh.Normals = new Vector3[mesh.Vertices.Length];
            var count = new int[mesh.Vertices.Length];
            for (int i = 0; i < count.Length; i++) { 
                count[i] = 0;
                mesh.Normals[i] = new Vector3(0, 0, 0);
            }

            for (int i = 0; i < mesh.FaceCount; ++i) { 
                var a = mesh.Indices[i * 3 + 0];
                var b = mesh.Indices[i * 3 + 1];
                var c = mesh.Indices[i * 3 + 2];

                var va = mesh.Vertices[a];
                var vb = mesh.Vertices[b];
                var vc = mesh.Vertices[c];

                var vab = SMMath.Vector3Subtract(vb, va);
                var vac = SMMath.Vector3Subtract(vc, va);
                var norm = SMMath.Vector3Normal(SMMath.Vector3Cross(vac, vab));

                mesh.Normals[a] = SMMath.Vector3Add(mesh.Normals[a], norm);
                mesh.Normals[b] = SMMath.Vector3Add(mesh.Normals[b], norm);
                mesh.Normals[c] = SMMath.Vector3Add(mesh.Normals[c], norm);

                count[a]++;
                count[b]++;
                count[c]++;
            }

            for (int i = 0; i < count.Length; i++)
            {
                mesh.Normals[i] = SMMath.Vector3Normal(SMMath.Vector3Scale(mesh.Normals[i], 1.0 / count[i]));
            }
        }

        private static Vector3 SampleNormal(FaceBuilder builder, Vector3 pt) {
            var uvw = SMMath.GetBarycentric(pt, builder.Vertices[0], builder.Vertices[1], builder.Vertices[2]);
            Vector3 norm = SMMath.Vector3Add(
                SMMath.Vector3Add(
                    SMMath.Vector3Scale(builder.Normals[0], uvw.X),
                    SMMath.Vector3Scale(builder.Normals[1], uvw.Y)
                    ), SMMath.Vector3Scale(builder.Normals[2], uvw.Z));
            return SMMath.Vector3Normal(norm);
        }
    }
}
