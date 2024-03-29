﻿namespace SMesh
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

    public struct Matrix {
        public double M00, M01, M02, M03;
        public double M10, M11, M12, M13;
        public double M20, M21, M22, M23;
        public double M30, M31, M32, M33;
    }

    public struct Segment2
    {
        public Vector2 A, B;

        public Segment2()
        {
            A = new Vector2();
            B = new Vector2();
        }

        public Segment2(Vector2 a, Vector2 b)
        {
            A = a;
            B = b;
        }
    }

    public struct Segment3 {
        public Vector3 A, B;

        public Segment3() { 
            A = new Vector3();
            B = new Vector3();
        }

        public Segment3(Vector3 a, Vector3 b) {
            A = a;
            B = b;
        }
    }

    public struct Plane {
        public Vector3 Normal;
        public double Dist;

        public Plane() {
            Normal = new Vector3();
            Dist = 0;
        }

        public Plane(Vector3 norm, double d) {
            Normal = norm;
            Dist = d;
        }
    }

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

    public class Sphere
    {
        public Vector3 Point;
        public double Radius;

        public Sphere()
        {
            Point = new Vector3();
            Radius = 0;
        }

        public Sphere(Vector3 point, double radius)
        {
            Point = point;
            Radius = radius;
        }

        public bool Intersects(Sphere sphere)
        {
            return SMMath.Vector3Distance(Point, sphere.Point) < Radius + sphere.Radius;
        }

        public bool Intersects(Vector3 rayOrigin, Vector3 rayDirection)
        {
            var OP = SMMath.Vector3Subtract(Point, rayOrigin);
            var d = SMMath.Vector3Dot(OP, rayDirection);
            var rp = SMMath.Vector3Add(SMMath.Vector3Scale(SMMath.Vector3Normal(rayDirection), d), rayOrigin);
            return SMMath.Vector3Distance(Point, rp) <= Radius;
        }

        public void Union(Sphere sphere)
        {
            var dist = SMMath.Vector3Distance(Point, sphere.Point);
            if (Radius >= dist + sphere.Radius)
            {
                return;
            }
            if (sphere.Radius >= dist + Radius)
            {
                Point = sphere.Point;
                Radius = sphere.Radius;
                return;
            }

            var dir = SMMath.Vector3Normal(SMMath.Vector3Subtract(sphere.Point, Point));
            var r = (dist + Radius + sphere.Radius) * 0.5;
            var m = r - Radius;

            var movement = SMMath.Vector3Scale(dir, m);

            Point = SMMath.Vector3Add(Point, movement);
            Radius = r;
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

    public struct Edge
    {
        public int A, B;

        public Edge(int a, int b) {
            A = a;
            B = b;
        }
    }
}
