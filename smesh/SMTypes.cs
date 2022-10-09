namespace SMesh
{

    public struct Vector2
    {
        public double X, Y;

        public Vector2()
        {
            this.X = 0;
            this.Y = 0;
        }

        public Vector2(double x, double y)
        {
            this.X = x;
            this.Y = y;
        }

        public Vector2(Vector2 v)
        {
            this.X = v.X;
            this.Y = v.Y;
        }
    }

    public struct Vector3
    {
        public double X, Y, Z;

        public Vector3()
        {
            this.X = 0;
            this.Y = 0;
            this.Z = 0;
        }

        public Vector3(double x, double y, double z)
        {
            this.X = x;
            this.Y = y;
            this.Z = z;
        }

        public Vector3(Vector3 v)
        {
            this.X = v.X;
            this.Y = v.Y;
            this.Z = v.Z;
        }
    }

    public struct Mesh
    {
        public int VertCount;
        public int FaceCount;

        public int[]? Indices;
        public Vector3[]? Vertices;
        public Vector3[]? Normals;
        public Vector2[]? TextureCoord;
    };

    public struct AABB
    {
        public Vector3 Min, Max;

        public AABB()
        {
            Min = new Vector3();
            Max = new Vector3();
        }

        public AABB(AABB b)
        {
            this.Min = new Vector3(b.Min);
            this.Max = new Vector3(b.Max);
        }

    };

}
