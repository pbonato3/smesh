using SMesh;

var mesh = SMObj.ParseFile("C:/Users/paolo/source/repos/smesh/test/sample.obj");
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
