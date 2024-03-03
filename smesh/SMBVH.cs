using System.Reflection.Metadata.Ecma335;

namespace SMesh
{

    public class SphereBVH {
        public class SphereNode {
            public int Id;
            public Sphere Sphere;
            public List<SphereNode> Children;

            public SphereNode()
            {
                Id = -1;
                Sphere = new Sphere();
                Children = new List<SphereNode>();
            }

            public SphereNode(Sphere sphere, int id = -1) {
                Id = id;
                Sphere = sphere;
                Children = new List<SphereNode>();
            }

            public void SetVolume(List<SphereNode> spheres) {
                //TODO: find a better algorithm for this
                if (spheres.Count < 1) {
                    return;
                }

                Sphere.Point = spheres[0].Sphere.Point;
                Sphere.Radius = spheres[0].Sphere.Radius;
                for (int i = 1; i < spheres.Count; i++) {
                    Sphere.Union(spheres[i].Sphere);
                }
            }

            public List<int> Search(Vector3 point, double range)
            {
                var results = new List<int>();
                SearchRecursion(this, ref results, new Sphere(point, range));
                return results;
            }

            private void SearchRecursion(SphereNode node, ref List<int> results, Sphere search)
            {
                if (!search.Intersects(node.Sphere))
                {
                    return;
                }

                if (node.Children.Count > 0)
                {
                    for (int i = 0; i < node.Children.Count; ++i)
                    {
                        SearchRecursion(node.Children[i], ref results, search);
                    }
                }
                else
                {
                    results.Add(node.Id);
                }
            }

            public List<int> Search(Vector3 rayOrigin, Vector3 rayDirection)
            {
                var results = new List<int>();
                SearchRecursion(this, ref results, rayOrigin, rayDirection);
                return results;
            }

            private void SearchRecursion(SphereNode node, ref List<int> results, Vector3 rayOrigin, Vector3 rayDirection)
            {
                if (!node.Sphere.Intersects(rayOrigin, rayDirection))
                {
                    return;
                }

                if (node.Children.Count > 0)
                {
                    for (int i = 0; i < node.Children.Count; ++i)
                    {
                        SearchRecursion(node.Children[i], ref results, rayOrigin, rayDirection);
                    }
                }
                else
                {
                    results.Add(node.Id);
                }
            }
        }

        public enum Sorting { 
            NO_SORT,
            X_SORT,
            Y_SORT,
            Z_SORT
        }
        public SphereNode Root;

        public SphereBVH() {
            Root = new SphereNode();
        }

        public void BuildBVH(List<SphereNode> nodes, int maxLeafCount = 8, int breadth = 8, Sorting sorting = Sorting.X_SORT) {
            Root = BuildBVHRecursion(nodes, maxLeafCount, breadth, sorting);
        }

        private SphereNode BuildBVHRecursion(List<SphereNode> nodes, int max, int breadth, Sorting sorting) {
            if (nodes.Count == 1) {
                return nodes[0];
            }

            var node = new SphereNode();
            node.SetVolume(nodes);

            if (nodes.Count <= max) {
                for (int i = 0; i < nodes.Count; ++i) {
                    node.Children.Add(nodes[i]);
                }
                return node;
            }

            var nextSort = sorting;
            switch (sorting)
            {
                case Sorting.X_SORT:
                    nodes.Sort((a, b) => { return a.Sphere.Point.X.CompareTo(b.Sphere.Point.X); });
                    nextSort = Sorting.Y_SORT;
                    break;
                case Sorting.Y_SORT:
                    nodes.Sort((a, b) => { return a.Sphere.Point.Y.CompareTo(b.Sphere.Point.Y); });
                    nextSort = Sorting.Z_SORT;
                    break;
                case Sorting.Z_SORT:
                    nodes.Sort((a, b) => { return a.Sphere.Point.Z.CompareTo(b.Sphere.Point.Z); });
                    nextSort = Sorting.X_SORT;
                    break;
                default:
                    break;
            }

            int nnodes = nodes.Count / breadth;
            int remain = nodes.Count % breadth;

            int c = 0;
            for (int i = 0; i < breadth; ++i) {
                var maxIndex = c + nnodes + (remain > 0 ? 1 : 0);
                remain--;
                var toAdd = new List<SphereNode>();
                for (; c < maxIndex; c++) {
                    toAdd.Add(nodes[c]);
                }
                if (toAdd.Count > 0) { 
                    node.Children.Add(BuildBVHRecursion(toAdd, max, breadth, nextSort));
                }
            }

            return node;
        }

        public void BuildBVH(List<Vector3> points, double range){
            var spheres = new List<SphereNode>();
            for (int i = 0; i < points.Count; ++i) {
                spheres.Add(new SphereNode(new Sphere(points[i], range), i));
            }
            BuildBVH(spheres);
        }

        public static SphereBVH BuildBVH(Mesh mesh)
        {
            var bvh = new SphereBVH();
            var spheres = new List<SphereNode>();
            for (int i = 0; i < mesh.FaceCount; ++i) {
                var a = mesh.Vertices[mesh.Indices[i * 3 + 0]];
                var b = mesh.Vertices[mesh.Indices[i * 3 + 1]];
                var c = mesh.Vertices[mesh.Indices[i * 3 + 2]];

                spheres.Add(new SphereNode(FaceSphere(a, b, c), i));
            }
            bvh.BuildBVH(spheres);
            return bvh;
        }

        public static Sphere FaceSphere(Vector3 a, Vector3 b, Vector3 c) {
            var d = SMMath.Vector3Divide(SMMath.Vector3Add(SMMath.Vector3Add(a, b), c), 3);
            var r = Math.Max(Math.Max(SMMath.Vector3Distance(a, d), SMMath.Vector3Distance(b, d)), SMMath.Vector3Distance(c, d));
            return new Sphere(d, r);
        }

        public List<int> Search(Sphere sphere)
        {
            return Search(sphere.Point, sphere.Radius);
        }

        public List<int> Search(Vector3 point, double range) {
            return Root.Search(point, range);
        }

        // TODO: Ray intersection
    }

    public static class SMBVH
    {
        public static BVHNode BuildBVH(Mesh mesh)
        {
            BVHNode root = new BVHNode();

            AABB[] boxes = new AABB[mesh.FaceCount];
            int[]  toAdd = new int[mesh.FaceCount];
            for (int i = 0; i < mesh.FaceCount; ++i) { 
                boxes[i] = FaceAABB(mesh, i);
                toAdd[i] = i;
            }

            RecursiveBVH(ref root, boxes, toAdd);

            return root;
        }

        public static void FindNeighbours(ref List<Vector3> neighbours, BVHNode node, Vector3 test, double dist) {
            if (!SMMath.IsInside(test, node.BBox))
            {
                return;
            }

            var pt = SMMath.AABBCenter(node.BBox);
            if (node.Index >= 0 && SMMath.Vector3Distance(test, pt) < dist)
            {
                neighbours.Add(pt);
                return;
            }

            FindNeighbours(ref neighbours, node.Right, test, dist);
            FindNeighbours(ref neighbours, node.Left, test, dist);
        }

        public static void CheckCollisions(ref List<int> res, BVHNode node, AABB box) {
            if (!SMMath.Collision(node.BBox, box)) {
                return;
            }

            if (node.Index >= 0) { 
                res.Add(node.Index);
                return;
            }

            CheckCollisions(ref res, node.Right, box);
            CheckCollisions(ref res, node.Left, box);
        }

        public static void RecursiveBVH(ref BVHNode currNode, AABB[] boxes, int[] toAdd) {
            if (toAdd.Count() == 1) {
                // Leaf situation
                currNode.BBox = boxes[toAdd[0]];
                currNode.Index = toAdd[0];
                return;
            }

            currNode.BBox = MergeAABBs(boxes, toAdd);
            var curr = new AABB(currNode.BBox);
            var toAddLeft = new List<int>();
            var toAddRight = new List<int>();

            // Split to add list in left and right
            if (toAdd.Count() == 2)
            {
                // Optimized to split 2 nodes
                toAddLeft.Add(toAdd[0]);
                toAddRight.Add(toAdd[1]);
            }
            else
            {
                // More thant two nodes: keep halving the bbox on bigger side untill there is a division
                while (toAddLeft.Count() == 0 || toAddRight.Count() == 0)
                {
                    toAddLeft = new List<int>();
                    toAddRight = new List<int>();

                    var half = HalveAABB(curr, true);

                    for (int i = 0; i < toAdd.Count(); ++i)
                    {
                        var box = boxes[toAdd[i]];
                        var boxCenter = SMMath.AABBCenter(box);

                        if (SMMath.IsInside(boxCenter, half))
                        {
                            toAddLeft.Add(toAdd[i]);
                        }
                        else
                        {
                            toAddRight.Add(toAdd[i]);
                        }
                    }

                    if (toAddLeft.Count() == 0)
                    {
                        curr = HalveAABB(curr, false);
                    }
                    else {
                        curr = half;
                    }
                }
            }

            currNode.Left = new BVHNode();
            currNode.Right = new BVHNode();

            RecursiveBVH(ref currNode.Left, boxes, toAddLeft.ToArray());
            RecursiveBVH(ref currNode.Right, boxes, toAddRight.ToArray());
        }

        private static AABB HalveAABB(AABB bbox, bool isLeft = true) {
            var half = new AABB(bbox);

            var center = SMMath.AABBCenter(bbox);
            var maxDim = SMMath.MaxDim(bbox);
            if (isLeft)
            {
                if (maxDim == 0)
                {
                    half.Max.X = center.X;
                }
                else if (maxDim == 1)
                {
                    half.Max.Y = center.Y;
                }
                else
                {
                    half.Max.Z = center.Z;
                }
            }
            else {
                if (maxDim == 0)
                {
                    half.Min.X = center.X;
                }
                else if (maxDim == 1)
                {
                    half.Min.Y = center.Y;
                }
                else
                {
                    half.Min.Z = center.Z;
                }
            }

            return half;
        }

        private static AABB FaceAABB(Mesh mesh, int f)
        {
            var a = mesh.Vertices[mesh.Indices[f * 3 + 0]];
            var b = mesh.Vertices[mesh.Indices[f * 3 + 1]];
            var c = mesh.Vertices[mesh.Indices[f * 3 + 2]];

            AABB bbox = new AABB();

            bbox.Min.X = Math.Min(a.X, Math.Min(b.X, c.X));
            bbox.Min.Y = Math.Min(a.Y, Math.Min(b.Y, c.Y));
            bbox.Min.Z = Math.Min(a.Z, Math.Min(b.Z, c.Z));

            bbox.Max.X = Math.Max(a.X, Math.Max(b.X, c.X));
            bbox.Max.Y = Math.Max(a.Y, Math.Max(b.Y, c.Y));
            bbox.Max.Z = Math.Max(a.Z, Math.Max(b.Z, c.Z));

            return bbox;
        }

        private static AABB MergeAABBs(AABB[] boxes, int[] toMerge) { 
            AABB bbox = SMMath.ZeroAABB();

            for (int i = 0; i < toMerge.Length; i++) {
                var box = boxes[toMerge[i]];

                bbox.Min.X = Math.Min(bbox.Min.X, box.Min.X);
                bbox.Min.Y = Math.Min(bbox.Min.Y, box.Min.Y);
                bbox.Min.Z = Math.Min(bbox.Min.Z, box.Min.Z);

                bbox.Max.X = Math.Max(bbox.Max.X, box.Max.X);
                bbox.Max.Y = Math.Max(bbox.Max.Y, box.Max.Y);
                bbox.Max.Z = Math.Max(bbox.Max.Z, box.Max.Z);
            }

            return bbox;
        }
    }

    public class BVHNode
    {
        public int Index;
        public AABB BBox;

        public BVHNode? Left;
        public BVHNode? Right;

        public BVHNode() {
            Index = -1;
            BBox  = new AABB();
            Left  = null;
            Right = null;
        }
    }


}
