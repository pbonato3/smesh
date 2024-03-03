using SMesh;


//var mesh = SMObj.ParseFile("C:/Users/paolo/source/repos/smesh/test/sample.obj");

//var root = SphereBVH.BuildBVH(mesh);
//RecursiveBVHToString(root.Root, 0);

//var meshA = SMObj.ParseFile("C:/Users/paolo/source/repos/smesh/test/sample.obj");
//var meshA = new Mesh();

//meshA.VertCount = 4;
//meshA.Vertices = new Vector3[4] { new Vector3(-10, -10, 0), new Vector3(10, -10, 0), new Vector3(-10, 10, 0), new Vector3(10, 10, 0) };
//meshA.FaceCount = 2;
//meshA.Indices = new int[6] { 0, 1, 2, 1, 3, 2 };


//meshA = new Mesh();
//meshA.VertCount = 4;
//meshA.Vertices = new Vector3[4] { new Vector3(-2, 0, -2), new Vector3(2, 0, 2), new Vector3(-2, 0, 2), new Vector3( 2, 0, -2) };
//meshA.FaceCount = 2;
//meshA.Indices = new int[6] { 0, 1, 2, 0, 3, 1 };

var meshA = SMObj.ParseFile("C:/Users/paolo/source/repos/smesh/test/cow.obj");
var moveA = new Vector3(0, 0, 1);
/*
for (int i = 0; i < meshA.VertCount; i++)
{
    meshA.Vertices[i] = SMMath.Vector3Add(meshA.Vertices[i], moveA);
    meshA.Vertices[i] = SMMath.Vector3Scale(meshA.Vertices[i], 10);   // TODO: automatic scale
}
*/

//var meshB = new Mesh();
//meshB.VertCount = 4;
//meshB.Vertices = new Vector3[4] { new Vector3(1.5, -20, -20), new Vector3(1.5, 20, 20), new Vector3(1.5, -20, 20), new Vector3(1.5, 20, -20) };
//meshB.FaceCount = 2;
//meshB.Indices = new int[6] { 0, 1, 2, 0, 3, 1 };


var meshB = SMObj.ParseFile("C:/Users/paolo/source/repos/smesh/test/bunny.obj");
var moveB = new Vector3(0.5, 0.5, 0.5);
for (int i = 0; i < meshB.VertCount; i++) {
    meshB.Vertices[i] = SMMath.Vector3Add(meshB.Vertices[i], moveB);
    //meshB.Vertices[i] = SMMath.Vector3Scale(meshB.Vertices[i], 10);
}

//var meshB = SMObj.ParseFile("C:/Users/paolo/source/repos/smesh/test/cow.obj");


if (meshA.Normals == null || meshA.Normals.Length != meshA.Vertices.Length) {
    SMCSG.RebuildNormals(ref meshA);
}
if (meshB.Normals == null || meshB.Normals.Length != meshB.Vertices.Length)
{
    SMCSG.RebuildNormals(ref meshB);
}

var result = SMCSG.RunOperation(meshA, meshB, SMCSG.Operation.AddOperation());
for (int i = 0; i < result.Length; ++i) { 
    SMObj.WriteFile("C:/Users/paolo/source/repos/smesh/test/out_" + i + ".obj", result[i]);
}



//SMObj.WriteFile("C:/Users/paolo/source/repos/smesh/test/out.obj", mesh);






