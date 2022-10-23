using System.ComponentModel;
using System.Runtime.CompilerServices;

namespace SMesh
{
    public static class SMMath
    {
        // Vector2 Methods
        public static Vector2 Vector2Add(Vector2 a, Vector2 b)
        {
            return new Vector2(a.X + b.X, a.Y + b.Y);
        }

        public static Vector2 Vector2Subtract(Vector2 a, Vector2 b)
        {
            return new Vector2(a.X - b.X, a.Y - b.Y);
        }

        public static Vector2 Vector2Scale(Vector2 v, double s)
        {
            return new Vector2(v.X * s, v.Y * s);
        }

        public static double Vector2Length(Vector2 v)
        {
            return Math.Sqrt(v.X * v.X + v.Y * v.Y);
        }

        public static double Vector2Dot(Vector2 a, Vector2 b)
        {
            return (a.X * b.X + a.Y * b.Y);
        }

        public static Vector2 Vector2Normal(Vector2 v) {
            var invMag = 1.0 / Math.Sqrt(v.X * v.X + v.Y * v.Y);
            return new Vector2(v.X * invMag, v.Y * invMag);
        }

        public static double Vector2Distance(Vector2 a, Vector2 b) {
            return Vector2Length(Vector2Subtract(b, a));
        }

        public static bool AreEquals(Vector2 a, Vector2 b, double tol) { 
            return Vector2Distance(a, b) <= tol;
        }

        // Vector3 Methods
        public static Vector3 Vector3Add(Vector3 a, Vector3 b)
        {
            return new Vector3(a.X + b.X, a.Y + b.Y, a.Z + b.Z);
        }

        public static Vector3 Vector3Subtract(Vector3 a, Vector3 b)
        {
            return new Vector3(a.X - b.X, a.Y - b.Y, a.Z - b.Z);
        }

        public static Vector3 Vector3Scale(Vector3 v, double s)
        {
            return new Vector3(v.X * s, v.Y * s, v.Z * s);
        }

        public static double Vector3Length(Vector3 v)
        {
            return Math.Sqrt(v.X * v.X + v.Y * v.Y + v.Z * v.Z);
        }

        public static double Vector3Dot(Vector3 a, Vector3 b)
        {
            return (a.X * b.X + a.Y * b.Y + a.Z * b.Z);
        }

        public static Vector3 Vector3Cross(Vector3 a, Vector3 b)
        {
            return new Vector3(a.Y * b.Z - a.Z * b.Y, a.Z * b.X - a.X * b.Z, a.X * b.Y - a.Y * b.X );
        }

        public static Vector3 Vector3Normal(Vector3 v)
        {
            var invMag = 1.0 / Math.Sqrt(v.X * v.X + v.Y * v.Y + v.Z * v.Z);
            return new Vector3(v.X * invMag, v.Y * invMag, v.Z * invMag);
        }
        public static double Vector3Distance(Vector3 a, Vector3 b)
        {
            return Vector3Length(Vector3Subtract(b, a));
        }

        // Matrix Methods

        public static Matrix IdentityMatrix() { 
            Matrix mat = new Matrix();

            mat.M00 = 1; mat.M01 = 0; mat.M02 = 0; mat.M03 = 0;
            mat.M10 = 0; mat.M11 = 1; mat.M12 = 0; mat.M13 = 0;
            mat.M20 = 0; mat.M21 = 0; mat.M22 = 1; mat.M23 = 0;
            mat.M30 = 0; mat.M31 = 0; mat.M32 = 0; mat.M33 = 1;

            return mat;
        }

        public static Vector2 Transform(Matrix m, Vector2 v)
        {
            var o = new Vector2();

            o.X = v.X * m.M00 + v.Y * m.M01 + m.M03;
            o.Y = v.X * m.M10 + v.Y * m.M11 + m.M13;

            return o;
        }

        public static Vector3 Transform(Matrix m, Vector3 v) {
            var o = new Vector3();

            o.X = v.X * m.M00 + v.Y * m.M01 + v.Z * m.M02 + m.M03;
            o.Y = v.X * m.M10 + v.Y * m.M11 + v.Z * m.M12 + m.M13;
            o.Z = v.X * m.M20 + v.Y * m.M21 + v.Z * m.M22 + m.M23;

            return o;
        }


        public static Vector3 TransformTo3D(Matrix m, Vector2 vec)
        {
            Vector3 v = new Vector3(vec.X, vec.Y, 0);
            var o = new Vector3();

            o.X = v.X * m.M00 + v.Y * m.M01 + v.Z * m.M02 + m.M03;
            o.Y = v.X * m.M10 + v.Y * m.M11 + v.Z * m.M12 + m.M13;
            o.Z = v.X * m.M20 + v.Y * m.M21 + v.Z * m.M22 + m.M23;

            return o;
        }

        public static Vector2 TransformTo2D(Matrix m, Vector3 v)
        {
            var o = new Vector2();

            o.X = v.X * m.M00 + v.Y * m.M01 + v.Z * m.M02 + m.M03;
            o.Y = v.X * m.M10 + v.Y * m.M11 + v.Z * m.M12 + m.M13;

            return o;
        }

        public static void PlaneTransformations(Plane plane, Vector3 origin, Vector3 xPoint, out Matrix toWorld, out Matrix toLocal) {
            toWorld = IdentityMatrix();
            toLocal = IdentityMatrix();

            var xdir = Vector3Normal(Vector3Subtract(xPoint, origin));
            var ydir = Vector3Normal(Vector3Cross(xdir, plane.Normal));

            toWorld.M00 = xdir.X;
            toWorld.M01 = xdir.Y;
            toWorld.M02 = xdir.Z;

            toWorld.M10 = ydir.X;
            toWorld.M11 = ydir.Y;
            toWorld.M12 = ydir.Z;

            toWorld.M20 = plane.Normal.X;
            toWorld.M21 = plane.Normal.Y;
            toWorld.M22 = plane.Normal.Z;

            toWorld.M30 = origin.X;
            toWorld.M31 = origin.Y;
            toWorld.M32 = origin.Z;

            toLocal = InvertMatrix(toWorld);
        }

        public static Matrix InvertMatrix(Matrix mat)
        {
            Matrix result = new Matrix();

            // Cache the matrix values (speed optimization)
            double a00 = mat.M00, a01 = mat.M01, a02 = mat.M02, a03 = mat.M03;
            double a10 = mat.M10, a11 = mat.M11, a12 = mat.M12, a13 = mat.M13;
            double a20 = mat.M20, a21 = mat.M21, a22 = mat.M22, a23 = mat.M23;
            double a30 = mat.M30, a31 = mat.M31, a32 = mat.M32, a33 = mat.M33;

            double b00 = a00 * a11 - a01 * a10;
            double b01 = a00 * a12 - a02 * a10;
            double b02 = a00 * a13 - a03 * a10;
            double b03 = a01 * a12 - a02 * a11;
            double b04 = a01 * a13 - a03 * a11;
            double b05 = a02 * a13 - a03 * a12;
            double b06 = a20 * a31 - a21 * a30;
            double b07 = a20 * a32 - a22 * a30;
            double b08 = a20 * a33 - a23 * a30;
            double b09 = a21 * a32 - a22 * a31;
            double b10 = a21 * a33 - a23 * a31;
            double b11 = a22 * a33 - a23 * a32;

            // Calculate the invert determinant (inlined to avoid double-caching)
            double invDet = 1.0f / (b00 * b11 - b01 * b10 + b02 * b09 + b03 * b08 - b04 * b07 + b05 * b06);

            result.M00 = (a11 * b11 - a12 * b10 + a13 * b09) * invDet;
            result.M10 = (-a01 * b11 + a02 * b10 - a03 * b09) * invDet;
            result.M20 = (a31 * b05 - a32 * b04 + a33 * b03) * invDet;
            result.M30 = (-a21 * b05 + a22 * b04 - a23 * b03) * invDet;
            result.M01 = (-a10 * b11 + a12 * b08 - a13 * b07) * invDet;
            result.M11 = (a00 * b11 - a02 * b08 + a03 * b07) * invDet;
            result.M21 = (-a30 * b05 + a32 * b02 - a33 * b01) * invDet;
            result.M31 = (a20 * b05 - a22 * b02 + a23 * b01) * invDet;
            result.M02 = (a10 * b10 - a11 * b08 + a13 * b06) * invDet;
            result.M12 = (-a00 * b10 + a01 * b08 - a03 * b06) * invDet;
            result.M22 = (a30 * b04 - a31 * b02 + a33 * b00) * invDet;
            result.M32 = (-a20 * b04 + a21 * b02 - a23 * b00) * invDet;
            result.M03 = (-a10 * b09 + a11 * b07 - a12 * b06) * invDet;
            result.M13 = (a00 * b09 - a01 * b07 + a02 * b06) * invDet;
            result.M23 = (-a30 * b03 + a31 * b01 - a32 * b00) * invDet;
            result.M33 = (a20 * b03 - a21 * b01 + a22 * b00) * invDet;

            return result;
        }


        // Segment Methods

        private static bool EvaluateOnSegment(Segment2 s, Vector2 pt ) {
            var minX = Math.Min(s.A.X, s.B.X);
            var maxX = Math.Max(s.A.X, s.B.X);

            var minY = Math.Max(s.A.Y, s.B.Y);
            var maxY = Math.Max(s.A.Y, s.B.Y);

            return pt.X >= minX && pt.X <= maxX && pt.Y >= minY && pt.Y <= maxY;
        }

        private static bool EvaluateOnSegment(Segment3 s, Vector3 pt)
        {
            var minX = Math.Min(s.A.X, s.B.X);
            var maxX = Math.Max(s.A.X, s.B.X);

            var minY = Math.Max(s.A.Y, s.B.Y);
            var maxY = Math.Max(s.A.Y, s.B.Y);

            var minZ = Math.Max(s.A.Z, s.B.Z);
            var maxZ = Math.Max(s.A.Z, s.B.Z);

            return pt.X >= minX && pt.X <= maxX && pt.Y >= minY && pt.Y <= maxY && pt.Z >= minZ && pt.Z <= maxZ;
        }

        private static double triSign(Vector2 p1, Vector2 p2, Vector2 p3)
        {
            return (p1.X - p3.X) * (p2.Y - p3.Y) - (p2.X - p3.X) * (p1.Y - p3.Y);
        }

        public static bool IsPointInTriangle(Vector2 pt, Vector2 v1, Vector2 v2, Vector2 v3)
        {
            double d1, d2, d3;
            bool has_neg, has_pos;

            d1 = triSign(pt, v1, v2);
            d2 = triSign(pt, v2, v3);
            d3 = triSign(pt, v3, v1);

            has_neg = (d1 < 0) || (d2 < 0) || (d3 < 0);
            has_pos = (d1 > 0) || (d2 > 0) || (d3 > 0);

            return !(has_neg && has_pos);
        }

        public static bool GetPointOnSegment(Vector2 pt, Vector2 va, Vector2 vb, double tol, out Vector2 res) {
            if (AreEquals(pt, va, tol)) {
                res = va;
                return true;
            }
            if (AreEquals(pt, vb, tol))
            {
                res = vb;
                return true;
            }
            res = SegmentClosestPoint(pt, new Segment2(va, vb));
            return AreEquals(pt, res, tol);
        }

        public static Vector2 LineClosestPoint(Vector2 pt, Segment2 seg) {
            var AB = Vector2Subtract(seg.B, seg.A);
            var AP = Vector2Subtract(pt, seg.A);
            double t = Vector2Dot(AP, AB) / Vector2Dot(AB, AB);
            return Vector2Add(seg.A, Vector2Scale(AB, t));
        }

        public static Vector3 LineClosestPoint(Vector3 pt, Segment3 seg)
        {
            var AB = Vector3Subtract(seg.B, seg.A);
            var AP = Vector3Subtract(pt, seg.A);
            double t = Vector3Dot(AP, AB) / Vector3Dot(AB, AB);
            return Vector3Add(seg.A, Vector3Scale(AB, t));
        }

        public static Vector2 SegmentClosestPoint(Vector2 pt, Segment2 seg)
        {
            var AB = Vector2Subtract(seg.B, seg.A);
            var AP = Vector2Subtract(pt, seg.A);
            double t = Vector2Dot(AP, AB) / Vector2Dot(AB, AB);
            if (t <= 0) { return seg.A; }
            if (t >= 1) { return seg.B; }
            return Vector2Add(seg.A, Vector2Scale(AB, t));
        }

        public static Vector3 SegmentClosestPoint(Vector3 pt, Segment3 seg)
        {
            var AB = Vector3Subtract(seg.B, seg.A);
            var AP = Vector3Subtract(pt, seg.A);
            double t = Vector3Dot(AP, AB) / Vector3Dot(AB, AB);
            if (t <= 0) { return seg.A; }
            if (t >= 1) { return seg.B; }
            return Vector3Add(seg.A, Vector3Scale(AB, t));
        }

        public static double DistancePointLine(Vector2 pt, Segment2 seg) {
            var cpt = LineClosestPoint(pt, seg);
            return Vector2Distance(cpt, pt);
        }

        public static double DistancePointLine(Vector3 pt, Segment3 seg)
        {
            var cpt = LineClosestPoint(pt, seg);
            return Vector3Distance(cpt, pt);
        }

        public static double DistancePointSegment(Vector2 pt, Segment2 seg)
        {
            var cpt = SegmentClosestPoint(pt, seg);
            return Vector2Distance(cpt, pt);
        }

        public static double DistancePointSegment(Vector3 pt, Segment3 seg)
        {
            var cpt = SegmentClosestPoint(pt, seg);
            return Vector3Distance(cpt, pt);
        }

        public static bool AreParallel(Segment2 l0, Segment2 l1, double tol) { 
            var AB0 = Vector2Normal(Vector2Subtract(l0.B, l0.A));
            var AB1 = Vector2Normal(Vector2Subtract(l1.B, l1.A));

            return 1 - Math.Abs(Vector2Dot(AB0, AB1)) < tol;
        }

        public static bool LineLineIntersection(Segment2 l0, Segment2 l1, out Vector2 pt) {
            pt = new Vector2();

            var M0 = (l0.B.Y - l0.A.Y)/(l0.B.X - l0.A.Y);
            var M1 = (l1.B.Y - l1.A.Y)/(l1.B.X - l1.A.Y);

            if (M0 == M1) { return false; }

            pt.X = (l1.A.Y - l0.A.Y + M0 * l0.A.X - M1 * l1.A.X) / (M0 - M1);
            pt.Y = M0 * (pt.X - l0.A.X) + l0.A.Y;

            return true;
        }

        public static bool LineSegmentIntersection(Segment2 l0, Segment2 l1, out Vector2 pt)
        {
            var result = LineLineIntersection(l0, l1, out pt);
            if (!result) { return result; }

            return EvaluateOnSegment(l1, pt);
        }

        public static bool SegmentSegmentIntersection(Segment2 l0, Segment2 l1, out Vector2 pt)
        {
            var result = LineSegmentIntersection(l0, l1, out pt);
            if (!result) { return result; }

            return EvaluateOnSegment(l0, pt);
        }



        // Plane Methods

        public static Plane PlaneFrom3Points(Vector3 pa, Vector3 pb, Vector3 pc) { 
            Plane plane = new Plane();
            var aaaa = Vector3Cross(Vector3Subtract(pa, pb), Vector3Subtract(pc, pb));
            plane.Normal = Vector3Normal(Vector3Cross(Vector3Subtract(pa, pb), Vector3Subtract(pc, pb)));
            plane.Dist = Vector3Dot(pb, plane.Normal);
            return plane;
        }

        public static bool PointOnPlane(Plane plane, Vector3 pt, double tol) {
            var d = Vector3Dot(pt, plane.Normal);
            return Math.Abs(d - plane.Dist) < tol;
        }

        public static bool IsLineParallelToPlane(Segment3 s, Plane p, double tol) {
            var ab = Vector3Subtract(s.B, s.A);
            return Math.Abs(Vector3Dot(p.Normal, ab)) < tol;
        }

        public static bool LinePlaneIntersection(Segment3 s, Plane p, out Vector3 pt) { 
            pt = new Vector3();
            var ab = Vector3Subtract(s.B, s.A);
            var dotAB = Vector3Dot(p.Normal, ab);

            if (dotAB == 0) { return false; }

            var ptOnPlane = Vector3Scale(p.Normal, p.Dist);
            var ap = Vector3Subtract(ptOnPlane, s.A);
            var dotAP = Vector3Dot(p.Normal, ap);

            pt = new Vector3(Vector3Add(s.A, Vector3Scale(ab, dotAP / dotAB)));
            return true;
        }

        public static bool SegmentPlaneIntersection(Segment3 s, Plane p, out Vector3 pt)
        {
            if (!LinePlaneIntersection(s, p, out pt)) {
                return false;
            }

            return (EvaluateOnSegment(s, pt));
        }


        // AABB Methods

        public static bool IsInside(Vector3 p, AABB bbox) {
            return p.X >= bbox.Min.X && p.X <= bbox.Max.X &&
                   p.Y >= bbox.Min.Y && p.Y <= bbox.Max.Y &&
                   p.Z >= bbox.Min.Z && p.Z <= bbox.Max.Z;
            //ti amo
        }

        public static bool Collision(AABB a, AABB b) {
            return a.Max.X >= b.Min.X && a.Min.X <= b.Max.X &&
                   a.Max.Y >= b.Min.Y && a.Min.Y <= b.Max.Y &&
                   a.Max.Z >= b.Min.Z && a.Min.Z <= b.Max.Z;
        }

        public static AABB ZeroAABB() { 
            AABB bbox = new AABB();
            bbox.Min = new Vector3(Double.MaxValue, Double.MaxValue, Double.MaxValue);
            bbox.Max = new Vector3(Double.MinValue, Double.MinValue, Double.MinValue);
            return bbox;
        }

        public static bool IsValid(AABB bbox) {
            return bbox.Min.X <= bbox.Max.X && bbox.Min.Y <= bbox.Max.Y && bbox.Min.Z <= bbox.Max.Z;
        }

        public static Vector3 AABBCenter(AABB bbox) {
            return Vector3Scale(Vector3Add(bbox.Min, bbox.Max), 0.5);
        }

        public static AABB Expand(AABB bbox, Vector3 pt) { 
            var result = new AABB();

            result.Min.X = Math.Min(bbox.Min.X, pt.X);
            result.Min.Y = Math.Min(bbox.Min.Y, pt.Y);
            result.Min.Z = Math.Min(bbox.Min.Z, pt.Z);

            result.Max.X = Math.Max(bbox.Max.X, pt.X);
            result.Max.Y = Math.Max(bbox.Max.Y, pt.Y);
            result.Max.Z = Math.Max(bbox.Max.Z, pt.Z);

            return result;
        }

        public static int MaxDim(AABB bbox)
        {
            var diagonal = Vector3Subtract(bbox.Max, bbox.Min);
            if (diagonal.X >= diagonal.Y && diagonal.X >= diagonal.Z) {
                return 0;
            }
            if (diagonal.Y >= diagonal.Z) {
                return 1;
            }
            return 2;
        }

        // Edge Methods

        public static Edge OrderedEdge(int a, int b) {
            return new Edge(Math.Min(a, b), Math.Max(a, b));
        }
    }
}