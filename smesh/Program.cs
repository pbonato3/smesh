using SMesh;
using System.Runtime.CompilerServices;



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


var meshA = SMObj.ParseFile("C:/Users/paolo/source/repos/smesh/test/sample.obj");
//var meshA = new Mesh();

//meshA.VertCount = 4;
//meshA.Vertices = new Vector3[4] { new Vector3(-10, -10, 0), new Vector3(10, -10, 0), new Vector3(-10, 10, 0), new Vector3(10, 10, 0) };
//meshA.FaceCount = 2;
//meshA.Indices = new int[6] { 0, 1, 2, 1, 3, 2 };


meshA = new Mesh();
meshA.VertCount = 4;
meshA.Vertices = new Vector3[4] { new Vector3(-2, -2, 0), new Vector3(2, 2, 0), new Vector3(-2, 2, 0), new Vector3( 2, -2, 0) };
meshA.FaceCount = 2;
meshA.Indices = new int[6] { 0, 1, 2, 0, 3, 1 };

meshA = SMObj.ParseFile("C:/Users/paolo/source/repos/smesh/test/terrain.obj");
for (int i = 0; i < meshA.VertCount; i++)
{
    meshA.Vertices[i] = SMMath.Vector3Scale(meshA.Vertices[i], 100);   // TODO: automatic scale
}

//var meshB = new Mesh();
//meshB.VertCount = 4;
//meshB.Vertices = new Vector3[4] { new Vector3(1.5, -20, -20), new Vector3(1.5, 20, 20), new Vector3(1.5, -20, 20), new Vector3(1.5, 20, -20) };
//meshB.FaceCount = 2;
//meshB.Indices = new int[6] { 0, 1, 2, 0, 3, 1 };


var meshB = SMObj.ParseFile("C:/Users/paolo/source/repos/smesh/test/armadillo.obj");
var move = new Vector3(0.2, 0.2, 0.2);
for (int i = 0; i < meshB.VertCount; i++) {
    //meshB.Vertices[i] = SMMath.Vector3Add(meshB.Vertices[i], move);
    meshB.Vertices[i] = SMMath.Vector3Scale(meshB.Vertices[i], 100);
}

//var meshB = SMObj.ParseFile("C:/Users/paolo/source/repos/smesh/test/cow.obj");


if (meshA.Normals == null || meshA.Normals.Length != meshA.Vertices.Length) {
    SMCSG.RebuildNormals(ref meshA);
}
if (meshB.Normals == null || meshB.Normals.Length != meshB.Vertices.Length)
{
    SMCSG.RebuildNormals(ref meshB);
}

var splitted = SMCSG.Split(meshB, meshA);
for (int i = 0; i < splitted.Length; ++i) { 
    SMObj.WriteFile("C:/Users/paolo/source/repos/smesh/test/out_" + i + ".obj", splitted[i]);
}



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
