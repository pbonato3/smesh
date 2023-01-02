﻿using SMesh;



// TEST TRANSFORM
//Vector3 a = new Vector3(2, 6, 0);
//Vector3 b = new Vector3(3, 6, 2);
//Vector3 c = new Vector3(7, 6,-2);

//Vector3 p = new Vector3(3, 5, 0);

Vector3 a = new Vector3(2, 2, 2);
Vector3 b = new Vector3(3, 6, 2);
Vector3 c = new Vector3(7, 6, 2);

Vector3 p = new Vector3(3, 5, 2);

Plane plane = SMMath.PlaneFrom3Points(a, b, c);
Matrix to3d = new Matrix();
Matrix to2d = new Matrix();
SMMath.PlaneTransformations(plane, a, b, out to3d, out to2d);


Console.WriteLine("A: " + Vector2String(SMMath.TransformTo2D(to2d, a)));
Console.WriteLine("B: " + Vector2String(SMMath.TransformTo2D(to2d, b)));
Console.WriteLine("C: " + Vector2String(SMMath.TransformTo2D(to2d, c)));
Console.WriteLine("P: " + Vector2String(SMMath.TransformTo2D(to2d, p)));

Console.WriteLine("A: " + Vector3String(SMMath.TransformTo3D(to3d, SMMath.TransformTo2D(to2d, a))));
Console.WriteLine("B: " + Vector3String(SMMath.TransformTo3D(to3d, SMMath.TransformTo2D(to2d, b))));
Console.WriteLine("C: " + Vector3String(SMMath.TransformTo3D(to3d, SMMath.TransformTo2D(to2d, c))));
Console.WriteLine("P: " + Vector3String(SMMath.TransformTo3D(to3d, SMMath.TransformTo2D(to2d, p))));

Console.WriteLine("O: " + Vector3String(SMMath.TransformTo3D(to3d, new Vector2(0, 0))));


//var meshA = SMObj.ParseFile("C:/Users/paolo/source/repos/smesh/test/sample.obj");
var meshA = new Mesh();

meshA.VertCount = 3;
meshA.Vertices = new Vector3[3] { new Vector3(-10, -10, 0), new Vector3(10, -10, 0), new Vector3(-10, 10, 0) };
meshA.FaceCount = 1;
meshA.Indices = new int[3] { 0, 1, 2 };


var meshB = new Mesh();

meshB.VertCount = 3;
//meshB.Vertices = new Vector3[3] { new Vector3(1.5, -10, -20), new Vector3(1.5, 10, 20), new Vector3(1.5, -10, 20) };
meshB.Vertices = new Vector3[3] { new Vector3(0, -5, -10), new Vector3(0, -5, 10), new Vector3(0, 15, 0) };
meshB.FaceCount = 1;
meshB.Indices = new int[3] { 0, 1, 2 };

var splitted = SMCSG.Split(meshA, meshB);
SMObj.WriteFile("C:/Users/paolo/source/repos/smesh/test/out.obj", splitted);



//SMObj.WriteFile("C:/Users/paolo/source/repos/smesh/test/out.obj", mesh);

//var root = SMBVH.BuildBVH(mesh);
//RecursiveBVHToString(root, 0);


/* TEST TRANSFORM
Vector3 a = new Vector3(2, 6, 0);
Vector3 b = new Vector3(3, 6, 2);
Vector3 c = new Vector3(7, 6,-2);

Vector3 p = new Vector3(3, 5, 0);

Plane plane = SMMath.PlaneFrom3Points(a, b, c);
Matrix to3d = new Matrix();
Matrix to2d = new Matrix();
SMMath.PlaneTransformations(plane, a, b, out to3d, out to2d);

Console.WriteLine("B: " + Vector3String(SMMath.Vector3Transform(to2d, b)));
Console.WriteLine("C: " + Vector3String(SMMath.Vector3Transform(to2d, c)));
Console.WriteLine("P: " + Vector3String(SMMath.Vector3Transform(to2d, p)));
*/


static string Vector2String(Vector2 vec)
{
    string str = "(";
    str += vec.X.ToString();
    str += ", ";
    str += vec.Y.ToString();
    str += ")";
    return str;
}

static string Vector3String(Vector3 vec) {
    string str = "(";
    str += vec.X.ToString();
    str += ", ";
    str += vec.Y.ToString();
    str += ", ";
    str += vec.Z.ToString();
    str += ")";
    return str;
}

static void RecursiveBVHToString(BVHNode currNode, int level)
{

    string spaces = "";
    for (int i = 0; i < level; i++)
    {
        spaces += "  ";
    }

    if (currNode.Index >= 0)
    {
        System.Console.WriteLine(spaces + currNode.Index + "\n");
        return;
    }

    System.Console.WriteLine(spaces + "{\n");

    RecursiveBVHToString(currNode.Left, level + 1);
    RecursiveBVHToString(currNode.Right, level + 1);

    System.Console.WriteLine(spaces + "}\n");

}
