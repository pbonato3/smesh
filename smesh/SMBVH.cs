using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Xml.Schema;

namespace SMesh
{
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
