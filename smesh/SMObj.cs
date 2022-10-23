using System;
using System.Collections.Generic;
using System.Globalization;

namespace SMesh
{
    public static class SMObj
    {
        static CultureInfo IC = CultureInfo.InvariantCulture;

        public static Mesh ParseFile(string path)
        {
            List<int> faces = new List<int>(); ;
            List<Vector3> vertices = new List<Vector3>();
            List<Vector3> normals  = new List<Vector3>();
            List<Vector2> textures = new List<Vector2>();


            string[] lines = System.IO.File.ReadAllLines(path);
            foreach (string line in lines)
            {
                var chunks = line.Split(' ');
                var cleanedChunks = new List<string>();
                for (int i = 0; i < chunks.Length; i++) {
                    if (chunks[i].Length > 0) { 
                        cleanedChunks.Add(chunks[i]);
                    }
                }
                chunks = cleanedChunks.ToArray();
                if (chunks.Length == 0) {
                    continue;
                }

                if (chunks[0] == "v")
                {
                    vertices.Add(new Vector3(double.Parse(chunks[1], IC), double.Parse(chunks[2], IC), double.Parse(chunks[3], IC)));
                    continue;
                }
                if (chunks[0] == "f")
                {
                    var ca = Int32.Parse(chunks[1].Split('/')[0]) - 1;
                    var cb = Int32.Parse(chunks[2].Split('/')[0]) - 1;
                    var cc = Int32.Parse(chunks[3].Split('/')[0]) - 1;

                    faces.Add(ca);
                    faces.Add(cb);
                    faces.Add(cc);

                    if (chunks.Length > 4) {
                        var cd = Int32.Parse(chunks[4].Split('/')[0]) - 1;

                        faces.Add(cc);
                        faces.Add(cd);
                        faces.Add(ca);
                    }
                    continue;
                }
                if (chunks[0] == "vn")
                {
                    normals.Add(new Vector3(double.Parse(chunks[1], IC), double.Parse(chunks[2], IC), double.Parse(chunks[3], IC)));
                    continue;
                }
                if (chunks[0] == "vt")
                {
                    textures.Add(new Vector2(double.Parse(chunks[1], IC), double.Parse(chunks[2], IC)));
                    continue;
                }

            }

            Mesh mesh = new Mesh();

            mesh.VertCount = vertices.Count;
            mesh.FaceCount = faces.Count / 3;

            mesh.Vertices =     vertices.Count() > 0 ? vertices.ToArray() : null;
            mesh.Indices =      faces.Count() > 0    ? faces.ToArray()    : null;
            mesh.Normals =      normals.Count() > 0  ? normals.ToArray()  : null;
            mesh.TextureCoord = textures.Count() > 0 ? textures.ToArray() : null;

            return mesh;
        }

        public static void WriteFile(string path, Mesh mesh)
        {
            var lineCount = mesh.VertCount + mesh.FaceCount;
            if (mesh.TextureCoord != null) {
                lineCount += mesh.VertCount;
            }
            if (mesh.Normals != null)
            {
                lineCount += mesh.VertCount;
            }

            string[] lines = new string[lineCount];

            var counter = 0;
            if (mesh.Vertices != null)
            {
                for (int i = 0; i < mesh.VertCount; i++)
                {
                    lines[counter] = $"v {mesh.Vertices[i].X.ToString(IC)} {mesh.Vertices[i].Y.ToString(IC)} {mesh.Vertices[i].Z.ToString(IC)} {1.0.ToString(IC)}";
                    counter++;
                }
            }

            if (mesh.TextureCoord != null)
            {
                for (int i = 0; i < mesh.TextureCoord.Count(); i++)
                {
                    lines[counter] = $"vt {mesh.TextureCoord[i].X.ToString(IC)} {mesh.TextureCoord[i].Y.ToString(IC)}";
                    counter++;
                }
            }

            if (mesh.Normals != null)
            {
                for (int i = 0; i < mesh.Normals.Count(); i++)
                {
                    lines[counter] = $"vn {mesh.Normals[i].X.ToString(IC)} {mesh.Normals[i].Y.ToString(IC)} {mesh.Normals[i].Z.ToString(IC)}";
                    counter++;
                }
            }

            if (mesh.Indices != null)
            {
                for (int i = 0; i < mesh.FaceCount; i++)
                {
                    var a = mesh.Indices[i * 3 + 0] + 1;
                    var b = mesh.Indices[i * 3 + 1] + 1;
                    var c = mesh.Indices[i * 3 + 2] + 1;

                    lines[counter] = $"f {a}//{a} {b}//{b} {c}//{c}";
                    counter++;
                }
            }

            System.IO.File.WriteAllLines(path, lines);
        }
    }
}
