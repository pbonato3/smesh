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


        // AABB Methods

        public static bool IsInside(Vector3 p, AABB bbox) {
            return p.X >= bbox.Min.X && p.X <= bbox.Max.X &&
                   p.Y >= bbox.Min.Y && p.Y <= bbox.Max.Y &&
                   p.Z >= bbox.Min.Z && p.Z <= bbox.Max.Z;
            //ti amo
        }

        public static AABB ZeroAABB() { 
            AABB bbox = new AABB();
            bbox.Min = new Vector3(Double.MaxValue, Double.MaxValue, Double.MaxValue);
            bbox.Max = new Vector3(Double.MinValue, Double.MinValue, Double.MinValue);
            return bbox;
        }

        public static Vector3 AABBCenter(AABB bbox) {
            return Vector3Scale(Vector3Add(bbox.Min, bbox.Max), 0.5);
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
    }
}