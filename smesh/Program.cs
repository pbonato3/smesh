using SMesh;

var mesh = SMObj.ParseFile("C:/Users/paolo/source/repos/smesh/test/sample.obj");
//SMObj.WriteFile("C:/Users/paolo/source/repos/smesh/test/out.obj", mesh);

var root = SMBVH.BuildBVH(mesh);
//RecursiveBVHToString(root, 0);


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
