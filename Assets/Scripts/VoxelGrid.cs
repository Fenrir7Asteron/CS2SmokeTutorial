using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using UnityEngine;
using UnityEngine.Serialization;
using UnityEngine.Windows;

public struct VoxelData
{
    public Vector4 position;
    public Color color;
    public Vector3Int flags;  // x - isVisible, y - is occupied by static mesh, z - flood fill remaining power
}

public struct VoxelGridParameters
{
    public Vector4 VoxelSize;
    public int CubeCount;
    public Vector3Int VoxelCountsPerAxis;
    public Vector4 VoxelGridExtents;
    public Vector4 VoxelGridCenter;
}

public struct SmokeParameters
{
    public Vector4 smokeSpawnPosition;
    public Vector4 smokeRadius;
    public Vector2 spawnTime;
    public Vector2 density; // x - smoke density, y - shadow density
    public Vector2 coefficients; // x - absorption coefficient, y - scattering coefficient
}

public class VoxelGrid : MonoBehaviour
{
    private const string SmokeVoxelPositionsKernel = "ComputeSmokeVoxelPositions";
    private const string SmokeVoxelColorsKernel = "ComputeSmokeVoxelColors";
    private const string FloodFillKernel = "ComputeFloodFill";
    
    public ComputeShader computeShader;
    public Vector3 zoneExtents;
    public float voxelSize;
    public int floodDistance = 16;
    [SerializeField] private float smokeDensity = 1.0f;
    [SerializeField] private float shadowDensity = 1.0f;
    [SerializeField, Range(0.0f, 1.0f)] private float absorptionCoefficient = 1.0f;
    public float ScatteringCoefficient { get => 1.0f - absorptionCoefficient; }
    [SerializeField] private bool showDebugVoxels;
    public Mesh debugMesh;
    public Material debugMaterial;

    private VoxelGridParameters _voxelGridParameters;


    private VoxelData[] _data;

    private static readonly int VoxelsPropertyID = Shader.PropertyToID("_voxels");
    private static readonly int VoxelGridParametersPropertyID = Shader.PropertyToID("_voxelGridParameters");
    private static readonly int SmokeParametersPropertyID = Shader.PropertyToID("_smokeParameters");
    private static readonly int VoxelSizePropertyID = Shader.PropertyToID("_voxelSize");
    private ComputeBuffer _voxelGridConstantBuffer;
    private ComputeBuffer _voxelsBuffer;
    private SmokeParameters _smokeParameters;
    private static readonly string VoxelizedSceneJsonPath = Application.streamingAssetsPath + "/VoxelizedSceneOccupation.json";

    public void Start()
    {
        List<int> occupiedVoxels;
        
        if (File.Exists(VoxelizedSceneJsonPath))
        {
            string json = Encoding.ASCII.GetString(File.ReadAllBytes(VoxelizedSceneJsonPath));
            occupiedVoxels = JsonUtility.FromJson<ListInt>(json).List;
            Debug.Log($"Occupied voxel count: {occupiedVoxels.Count}");
        }
        else
        {
            occupiedVoxels = new List<int>();
        }
        
        InitVoxelGrid(occupiedVoxels);
    }

    public void OnDestroy()
    {
        _voxelGridConstantBuffer.Dispose();
    }

    public void Update()
    {
        RecalculateVoxelsPositionsAndVisibility();

        if (showDebugVoxels)
        {
            Bounds bounds = new Bounds(_voxelGridParameters.VoxelGridCenter,
                _voxelGridParameters.VoxelGridExtents * 2.0f);
            
            debugMaterial.SetBuffer(VoxelsPropertyID, _voxelsBuffer);
            debugMaterial.SetFloat(VoxelSizePropertyID, voxelSize);
            
            Graphics.DrawMeshInstancedProcedural(debugMesh, 0, debugMaterial, bounds, _data.Length);
        }
    }

    private void InitVoxelGrid(in List<int> occupiedVoxels)
    {
        Vector3 zoneSizes = zoneExtents * 2 / voxelSize;

        _voxelGridParameters.VoxelSize = Vector4.one * voxelSize;
        _voxelGridParameters.VoxelCountsPerAxis = new Vector3Int((int) zoneSizes.x, (int) zoneSizes.y, (int) zoneSizes.z);
        _voxelGridParameters.CubeCount = _voxelGridParameters.VoxelCountsPerAxis.x * _voxelGridParameters.VoxelCountsPerAxis.y * _voxelGridParameters.VoxelCountsPerAxis.z;
        _voxelGridParameters.VoxelGridExtents = new Vector4(zoneExtents.x, zoneExtents.y, zoneExtents.z);
        Vector3 voxelGridCenter = transform.position;
        _voxelGridParameters.VoxelGridCenter = new Vector4(voxelGridCenter.x, voxelGridCenter.y, voxelGridCenter.z);

        Debug.Log($"Zone sizes integer: {_voxelGridParameters.VoxelCountsPerAxis}, voxel count: {_voxelGridParameters.CubeCount}");
        _data = new VoxelData[_voxelGridParameters.CubeCount];

        for (int voxelIdx = 0; voxelIdx < _voxelGridParameters.CubeCount; ++voxelIdx)
        {
            int isOccupied = occupiedVoxels.Contains(voxelIdx) ? 1 : 0;
            
            VoxelData voxelDataData = new VoxelData()
            {
                position = new Vector4(),
                color = Color.magenta,
                flags = new Vector3Int(0, isOccupied, 0),
            };
        
            _data[voxelIdx] = voxelDataData;
        }


        _voxelsBuffer = new ComputeBuffer(_data.Length,
            System.Runtime.InteropServices.Marshal.SizeOf(typeof(VoxelData)));
        _voxelsBuffer.SetData(_data);

        _voxelGridConstantBuffer = new ComputeBuffer(1, 
            System.Runtime.InteropServices.Marshal.SizeOf(typeof(VoxelGridParameters)),
            ComputeBufferType.Constant);
        _voxelGridConstantBuffer.SetData(new VoxelGridParameters[] {_voxelGridParameters});

        RecalculateVoxelsPositionsAndVisibility(true);
        RandomizeVoxelColors();
        
        _smokeParameters = new SmokeParameters()
        {
            spawnTime = new Vector2(-1, -1),
            density = new Vector2(smokeDensity, shadowDensity),
            coefficients = new Vector2(absorptionCoefficient, ScatteringCoefficient),
        };
    }

    public void SpawnSmoke(Vector3 spawnPosition, Vector3 smokeRadius, float smokeSpawnDuration)
    {
        CleanAllVoxels();
        Vector3 voxelGridCenter = transform.position;
        int voxelIdx = CalculateVoxelIdxByWorldPosition(spawnPosition, voxelGridCenter, zoneExtents);

        _data[voxelIdx].flags.z = floodDistance;
        
        _voxelsBuffer.SetData(_data);

        _smokeParameters.smokeSpawnPosition = new Vector4(spawnPosition.x, spawnPosition.y, spawnPosition.z, 0.0f);
        _smokeParameters.smokeRadius = new Vector4(smokeRadius.x, smokeRadius.y, smokeRadius.z, 0.0f);
        _smokeParameters.spawnTime = new Vector2(smokeSpawnDuration, Time.time);
        _smokeParameters.density = new Vector2(smokeDensity, shadowDensity);
        _smokeParameters.coefficients = new Vector2(absorptionCoefficient, ScatteringCoefficient);
        
        RecalculateVoxelsPositionsAndVisibility();
    }

    private int CalculateVoxelIdxByWorldPosition(in Vector3 worldPosition, in Vector3 voxelGridCenter, in Vector3 extents)
    {
        Vector3 zoneSize = extents * 2;
        Vector3 voxelPosition = worldPosition - voxelGridCenter + extents;
        
        voxelPosition.x /= zoneSize.x;
        voxelPosition.y /= zoneSize.y;
        voxelPosition.z /= zoneSize.z;
        
        voxelPosition = Vector3.Scale(voxelPosition, _voxelGridParameters.VoxelCountsPerAxis);
        
        Vector3Int voxelPositionInt = new Vector3Int((int) voxelPosition.x, (int) voxelPosition.y, (int) voxelPosition.z);
        voxelPositionInt.x = Math.Clamp(voxelPositionInt.x, 0, _voxelGridParameters.VoxelCountsPerAxis.x - 1);
        voxelPositionInt.y = Math.Clamp(voxelPositionInt.y, 0, _voxelGridParameters.VoxelCountsPerAxis.y - 1);
        voxelPositionInt.z = Math.Clamp(voxelPositionInt.z, 0, _voxelGridParameters.VoxelCountsPerAxis.z - 1);
        
        int voxelIdx = voxelPositionInt.x
                   + voxelPositionInt.y * _voxelGridParameters.VoxelCountsPerAxis.x
                   + voxelPositionInt.z * _voxelGridParameters.VoxelCountsPerAxis.x *
                   _voxelGridParameters.VoxelCountsPerAxis.y;
        
        Debug.Log($"World position: {worldPosition}, voxel position: {voxelPositionInt}, voxel idx: {voxelIdx}");
        
        return voxelIdx;
    }
    
    private Vector3 VoxelIdxToVoxelPosition(int voxelIdx, Vector2Int voxelCounts, Vector3 voxelZoneCenter,
        Vector3 voxelZoneExtents, Vector3 voxelSize)
    {
        float xPos = voxelIdx % voxelCounts.x;
        float yPos = (voxelIdx / voxelCounts.x) % voxelCounts.y;
        float zPos = voxelIdx / (voxelCounts.x * voxelCounts.y);
        return Vector3.Scale(new Vector3(xPos, yPos, zPos), voxelSize) 
            + voxelSize * 0.5f 
            + voxelZoneCenter - voxelZoneExtents;
    }

    private void RecalculateVoxelsPositionsAndVisibility(bool force = false)
    {
        if (!IsSmokeSpawned() && !force) // smoke was not created yet
        {
            return;
        }
        
        // Compute Flood fill
        int kernelIndex = computeShader.FindKernel(FloodFillKernel);

        computeShader.SetBuffer(kernelIndex, VoxelsPropertyID, _voxelsBuffer);
        computeShader.SetConstantBuffer(VoxelGridParametersPropertyID, _voxelGridConstantBuffer, 0,
            System.Runtime.InteropServices.Marshal.SizeOf(typeof(VoxelGridParameters)));

        computeShader.Dispatch(kernelIndex,
            _voxelGridParameters.CubeCount / 8 / 8,
            _voxelGridParameters.CubeCount / 8 / 8,
            1);
        
        // Compute positions
        kernelIndex = computeShader.FindKernel(SmokeVoxelPositionsKernel);

        computeShader.SetBuffer(kernelIndex, VoxelsPropertyID, _voxelsBuffer);
        computeShader.SetConstantBuffer(VoxelGridParametersPropertyID, _voxelGridConstantBuffer, 0,
            System.Runtime.InteropServices.Marshal.SizeOf(typeof(VoxelGridParameters)));
        
        ComputeBuffer smokeBuffer = new ComputeBuffer(1,
            System.Runtime.InteropServices.Marshal.SizeOf(typeof(SmokeParameters)),
            ComputeBufferType.Structured);
        
        smokeBuffer.SetData(new SmokeParameters[] {_smokeParameters});
        computeShader.SetBuffer(kernelIndex, SmokeParametersPropertyID, smokeBuffer);

        computeShader.Dispatch(kernelIndex,
            _voxelGridParameters.CubeCount / 8 / 8,
            _voxelGridParameters.CubeCount / 8 / 8,
            1);
        
        smokeBuffer.Dispose();
    }

    private bool IsSmokeSpawned()
    {
        return _smokeParameters.spawnTime.x >= 0;
    }

    private void RandomizeVoxelColors()
    {
        int kernelIndex = computeShader.FindKernel(SmokeVoxelColorsKernel);
        
        computeShader.SetBuffer(kernelIndex, VoxelsPropertyID, _voxelsBuffer);
        computeShader.SetConstantBuffer(VoxelGridParametersPropertyID, _voxelGridConstantBuffer, 0, 
            System.Runtime.InteropServices.Marshal.SizeOf(typeof(VoxelGridParameters)));

        computeShader.Dispatch(kernelIndex,
            _voxelGridParameters.CubeCount / 8 / 8, 
            _voxelGridParameters.CubeCount / 8 / 8, 
            1);

        _voxelsBuffer.GetData(_data);
    }

    private void CleanAllVoxels()
    {
        for (int i = 0; i < _data.Length; i++)
        {
            _data[i].flags = new Vector3Int(0, _data[i].flags.y, 0);
        }
    }

    private void OnGUI()
    {
        if (GUI.Button(new Rect(0, 0, 100, 50), "Bake scene as voxels"))
        {
            VoxelizeSceneGeometry();
        }
    }

    private void VoxelizeSceneGeometry()
    {
        MeshFilter[] sceneMeshes = FindObjectsOfType<MeshFilter>()
            .Where(meshFilter => meshFilter.gameObject.layer == LayerMask.NameToLayer("SceneGeometry")
                                 && meshFilter.gameObject.activeInHierarchy)
            .ToArray();

        ListInt occupiedVoxels = new ListInt(FindOccupiedVoxels(sceneMeshes));
        Debug.Log($"Occupied voxel count: {occupiedVoxels.List.Count}");
        
        foreach (int voxelIdx in occupiedVoxels.List)
        {
            _data[voxelIdx].flags.y = 1;
        }
        
        _voxelsBuffer.SetData(_data);

        string json = JsonUtility.ToJson(occupiedVoxels);
        File.WriteAllBytes(VoxelizedSceneJsonPath, Encoding.ASCII.GetBytes(json));
    }

    private List<int> FindOccupiedVoxels(MeshFilter[] sceneMeshes)
    {
        var intersectingVoxels = new List<int>();
        
        float voxelExtent = voxelSize * 0.5f;

        for (int voxelIdx = 0; voxelIdx < _data.Length; ++voxelIdx)
        {
            List<Vector3> voxelVertices = new List<Vector3>();
            Vector3 voxelPosition = VoxelIdxToVoxelPosition(
                voxelIdx, 
                new Vector2Int(_voxelGridParameters.VoxelCountsPerAxis.x, _voxelGridParameters.VoxelCountsPerAxis.y),
                _voxelGridParameters.VoxelGridCenter,
                _voxelGridParameters.VoxelGridExtents, 
                _voxelGridParameters.VoxelSize);
                        
            voxelVertices.Add(voxelPosition + (Vector3.right + Vector3.up + Vector3.forward) * voxelExtent);
            voxelVertices.Add(voxelPosition + (Vector3.right + Vector3.up + Vector3.back) * voxelExtent);
            voxelVertices.Add(voxelPosition + (Vector3.right + Vector3.down + Vector3.forward) * voxelExtent);
            voxelVertices.Add(voxelPosition + (Vector3.right + Vector3.down + Vector3.back) * voxelExtent);
            voxelVertices.Add(voxelPosition + (Vector3.left + Vector3.up + Vector3.forward) * voxelExtent);
            voxelVertices.Add(voxelPosition + (Vector3.left + Vector3.up + Vector3.back) * voxelExtent);
            voxelVertices.Add(voxelPosition + (Vector3.left + Vector3.down + Vector3.forward) * voxelExtent);
            voxelVertices.Add(voxelPosition + (Vector3.left + Vector3.down + Vector3.back) * voxelExtent);
            
            List<PointWithDir> voxelFaceNormals = GetVoxelFaceNormals(voxelPosition);

            foreach (MeshFilter meshFilter in sceneMeshes)
            {
                bool voxelIntersectsMeshBoundingBox = meshFilter.gameObject.GetComponent<Collider>().bounds
                    .Intersects(new Bounds(voxelPosition, voxelExtent * 2.0f * Vector3.one));
                
                if (!voxelIntersectsMeshBoundingBox)
                {
                    continue;
                }
                
                //Debug.Log($"{meshFilter.gameObject.name}, {voxelIdx}");

                int[] triangleIndices = meshFilter.sharedMesh.GetTriangles(0);

                List<Vector3> vertices = new List<Vector3>();
                meshFilter.sharedMesh.GetVertices(vertices);
                Matrix4x4 localToWorldMatrix = meshFilter.transform.localToWorldMatrix;


                for (int idx = 0; idx < triangleIndices.Length; idx += 3)
                {
                    Vector3 v1 = localToWorldMatrix.MultiplyPoint3x4(vertices[triangleIndices[idx]] * 0.999f);
                    Vector3 v2 = localToWorldMatrix.MultiplyPoint3x4(vertices[triangleIndices[idx + 1]] * 0.999f);
                    Vector3 v3 = localToWorldMatrix.MultiplyPoint3x4(vertices[triangleIndices[idx + 2]] * 0.999f);

                    List<Vector3> triangleVertices = new List<Vector3>();
                    triangleVertices.Add(v1);
                    triangleVertices.Add(v2);
                    triangleVertices.Add(v3);

                    List<PointWithDir> intersectionLines = new List<PointWithDir>();
                    PointWithDir triangleFaceNormal = GetSceneTriangleFaceNormals(v1, v2, v3);

                    for (int triangleVertexIdx = 0; triangleVertexIdx < triangleVertices.Count; ++triangleVertexIdx)
                    {
                        triangleVertices[triangleVertexIdx] -= triangleFaceNormal.Dir * 0.001f;
                    }

                    intersectionLines.Add(triangleFaceNormal);
                    intersectionLines.AddRange(voxelFaceNormals);
                    intersectionLines.AddRange(GetEdgeCrossNormals(triangleVertices, triangleFaceNormal,
                        voxelFaceNormals));
                    
                    bool intersects = true;
                    foreach (PointWithDir intersectionLine in intersectionLines)
                    {
                        if (!ObjectsIntersectOnProjectionLine(voxelVertices, triangleVertices, intersectionLine))
                        {
                            intersects = false;
                            break;
                        }
                    }

                    if (intersects)
                    {
                        intersectingVoxels.Add(voxelIdx);
                        //Debug.Log($"INTERSECTS {voxelIdx}");
                        break;
                    }
                }
            }
        }
        
        return intersectingVoxels;
    }

    private bool ObjectsIntersectOnProjectionLine(in List<Vector3> vertices1, in List<Vector3> vertices2,
        in PointWithDir intersectionLine)
    {
        //Debug.Log("Voxel");
        GetProjectionDistances(vertices1, intersectionLine, out float minDistance1, out float maxDistance1, false);
        //Debug.Log("Mesh");
        GetProjectionDistances(vertices2, intersectionLine, out float minDistance2, out float maxDistance2, false);
        
        return maxDistance1 > minDistance2 && maxDistance2 > minDistance1;
    }

    private void GetProjectionDistances(in List<Vector3> vertices, in PointWithDir intersectionLine,
        out float minDistance, out float maxDistance, bool drawDebug = false)
    {
        minDistance = Single.PositiveInfinity;
        maxDistance = Single.NegativeInfinity;
        foreach (Vector3 vertex in vertices)
        {
            Vector3 toVertex = (vertex - intersectionLine.Point);
            Vector3 projection = Vector3.Project(toVertex, intersectionLine.Dir);
            float distance = projection.sqrMagnitude * Math.Sign(Vector3.Dot(projection, intersectionLine.Dir));
            minDistance = Math.Min(minDistance, distance);
            maxDistance = Math.Max(maxDistance, distance);

            if (drawDebug)
            {
                Debug.DrawLine(intersectionLine.Point, vertex, Color.cyan, Single.PositiveInfinity);
               // Debug.DrawLine(intersectionLine.Point, intersectionLine.Point + intersectionLine.Dir, Color.red, Single.PositiveInfinity);
                Debug.DrawLine(intersectionLine.Point, intersectionLine.Point + projection, Color.green, Single.PositiveInfinity);
                Debug.DrawLine(intersectionLine.Point, intersectionLine.Point + intersectionLine.Dir, Color.yellow, Single.PositiveInfinity);
                //Debug.DrawLine(intersectionLine.Point, intersectionLine.Point + intersectionLine.Dir * distance, Color.red, Single.PositiveInfinity);
            }
            
        }

        if (drawDebug)
        {
            Debug.Log($"Min: {minDistance}, Max: {maxDistance}");
        }
    }

    private List<PointWithDir> GetVoxelFaceNormals(in Vector3 voxelPosition)
    {
        List<PointWithDir> voxelFaceNormals = new List<PointWithDir>();
        
        PointWithDir voxelFaceNormal1 = GetVoxelFaceNormal(Vector3.right, voxelPosition);
        PointWithDir voxelFaceNormal2 = GetVoxelFaceNormal(Vector3.forward, voxelPosition);
        PointWithDir voxelFaceNormal3 = GetVoxelFaceNormal(Vector3.up, voxelPosition);
        
        voxelFaceNormals.Add(voxelFaceNormal1);
        voxelFaceNormals.Add(voxelFaceNormal2);
        voxelFaceNormals.Add(voxelFaceNormal3);

        return voxelFaceNormals;
    }

    private PointWithDir GetSceneTriangleFaceNormals(in Vector3 v1, in Vector3 v2, in Vector3 v3)
    {
        Vector3 triangleCenter = (v1 + v2 + v3) / 3;
        Vector3 triangleNormal = (Vector3.Cross(v2 - v1, v3 - v1)).normalized;

        return new PointWithDir(triangleCenter, triangleNormal);
    }

    private PointWithDir GetVoxelFaceNormal(Vector3 dir, Vector3 voxelPosition)
    {
        Vector3 voxelNormal = dir;
        Vector3 voxelNormalCenter = voxelPosition + voxelNormal * (voxelSize * 0.5f);
        return new PointWithDir(voxelNormalCenter, voxelNormal);
    }

    private IEnumerable<PointWithDir> GetEdgeCrossNormals(in List<Vector3> triangleVertices,
        in PointWithDir triangleFaceNormal, in List<PointWithDir> voxelFaceNormals)
    {
        Vector3 triangleEdgeNormal0 = Vector3.Cross((triangleVertices[1] - triangleVertices[0]).normalized, triangleFaceNormal.Dir).normalized;
        Vector3 triangleEdgeNormal1 = Vector3.Cross((triangleVertices[2] - triangleVertices[1]).normalized, triangleFaceNormal.Dir).normalized;
        Vector3 triangleEdgeNormal2 = Vector3.Cross((triangleVertices[0] - triangleVertices[2]).normalized, triangleFaceNormal.Dir).normalized;
        Vector3 triangleEdgeCenter0 = triangleVertices[0] + (triangleVertices[1] - triangleVertices[0]) * 0.5f;
        Vector3 triangleEdgeCenter1 = triangleVertices[1] + (triangleVertices[2] - triangleVertices[1]) * 0.5f;
        Vector3 triangleEdgeCenter2 = triangleVertices[2] + (triangleVertices[0] - triangleVertices[2]) * 0.5f;
        List<PointWithDir> triangleEdgeNormals = new List<PointWithDir>()
        {
            new PointWithDir(triangleEdgeCenter0, triangleEdgeNormal0),
            new PointWithDir(triangleEdgeCenter1, triangleEdgeNormal1),
            new PointWithDir(triangleEdgeCenter2, triangleEdgeNormal2),
        };
        
        // Debug.DrawLine(triangleEdgeCenter0, triangleEdgeCenter0 + triangleEdgeNormal0 * 0.2f, Color.magenta, float.PositiveInfinity);
        // Debug.DrawLine(triangleEdgeCenter1, triangleEdgeCenter1 + triangleEdgeNormal1 * 0.2f, Color.black, float.PositiveInfinity);
        // Debug.DrawLine(triangleEdgeCenter2, triangleEdgeCenter2 + triangleEdgeNormal2 * 0.2f, Color.green, float.PositiveInfinity);
        // Debug.DrawLine(triangleVertices[0], triangleVertices[0] + Vector3.up * 0.2f, Color.red, float.PositiveInfinity);
        // Debug.DrawLine(triangleVertices[1], triangleVertices[1] + Vector3.up * 0.2f, Color.red, float.PositiveInfinity);
        // Debug.DrawLine(triangleVertices[2], triangleVertices[2] + Vector3.up * 0.2f, Color.red, float.PositiveInfinity);
        // DebugDrawTriangle(triangleVertices[0], triangleVertices[1], triangleVertices[2], Color.cyan);
        // Debug.DrawLine(triangleFaceNormal.Point, triangleFaceNormal.Point + triangleFaceNormal.Dir * 0.3f, Color.yellow, float.PositiveInfinity);

        List<PointWithDir> edgeCrossNormals = new List<PointWithDir>();
        for (int i = 0; i < voxelFaceNormals.Count; ++i)
        {
            for (int j = i + 1; j < triangleEdgeNormals.Count; ++j)
            {
                Vector3 normal = Vector3.Cross(voxelFaceNormals[i].Dir, triangleEdgeNormals[j].Dir).normalized;
                edgeCrossNormals.Add(new PointWithDir(triangleEdgeNormals[j].Point, normal));
                //Debug.DrawLine(triangleEdgeNormals[j].Point, triangleEdgeNormals[j].Point + normal * 0.4f, Color.blue, float.PositiveInfinity);
            }
        }

        return edgeCrossNormals;
    }

    private void DebugDrawTriangle(Vector3 v1, Vector3 v2, Vector3 v3, Color color)
    {
        Debug.DrawLine(v1, v2, color, float.PositiveInfinity);
        Debug.DrawLine(v2, v3, color, float.PositiveInfinity);
        Debug.DrawLine(v3, v1, color, float.PositiveInfinity);
    }

    private void OnDrawGizmos()
    {
        Gizmos.color = Color.red;
        Gizmos.DrawWireCube(transform.position, zoneExtents * 2.0f);
    }

    public ComputeBuffer GetVoxelBuffer()
    {
        return _voxelsBuffer;
    }

    public VoxelGridParameters GetGridParameters()
    {
        return _voxelGridParameters;
    }

    public SmokeParameters GetSmokeParameters()
    {
        _smokeParameters.density = new Vector2(smokeDensity, shadowDensity);
        _smokeParameters.coefficients = new Vector2(absorptionCoefficient, ScatteringCoefficient);
        return _smokeParameters;
    }

    public int GetMaxFloodValue()
    {
        return floodDistance;
    }
}

[Serializable]
public class ListInt
{
    public List<int> List;

    public ListInt(List<int> list)
    {
        List = list;
    }
}

public class PointWithDir
{
    public Vector3 Point;
    public Vector3 Dir;

    public PointWithDir(Vector3 point, Vector3 dir)
    {
        Point = point;
        Dir = dir;
    }
}