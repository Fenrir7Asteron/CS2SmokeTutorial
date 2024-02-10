using UnityEngine;

public class Raymarching : MonoBehaviour 
{
    [SerializeField] private RenderTexture smokeColorTextureLowRes;
    [SerializeField] private RenderTexture smokeMaskTextureLowRes;
    [SerializeField] private RenderTexture smokeColorTextureFullRes;
    [SerializeField] private RenderTexture smokeMaskTextureFullRes;
    [SerializeField] private Material raymarchingMaterial;
    [SerializeField] private Material blendingMaterial;
    [SerializeField] private Light sunLight;
    [SerializeField] private VoxelGrid voxelGrid;
    
    private Camera _currentCamera;
    
    private static readonly int VoxelDataPropertyID = Shader.PropertyToID("_smokeVoxels");
    private static readonly int SmokeParametersPropertyID = Shader.PropertyToID("_smokeParameters");
    private static readonly int VoxelGridParametersPropertyID = Shader.PropertyToID("_voxelGridParameters");
    private static readonly int FrustumCornersESPropertyID = Shader.PropertyToID("_FrustumCornersES");
    private static readonly int CameraInvViewMatrixPropertyID = Shader.PropertyToID("_CameraInvViewMatrix");
    private static readonly int CameraWsPropertyID = Shader.PropertyToID("_CameraWS");
    private static readonly int LightDirPropertyID = Shader.PropertyToID("_LightDir");
    private static readonly int LightColorPropertyID = Shader.PropertyToID("_LightColor");
    private static readonly int MouseYPositionPropertyID = Shader.PropertyToID("_MouseYPosition");
    private static readonly int MaxFloodValuePropertyID = Shader.PropertyToID("_MaxFloodValue");
    private static readonly int _DraineCoefficientsPropertyID = Shader.PropertyToID("_DraineCoefficients");
    private static readonly int _WaterDropletDiameterPropertyID = Shader.PropertyToID("_WaterDropletDiameter");

    public Camera CurrentCamera 
    { 
        get 
        {
            if (_currentCamera != null)
            {
                return _currentCamera; 
            }

#if UNITY_EDITOR
            _currentCamera = Camera.main;
            //_currentCamera = GetSceneEditorCamera();
#else
            _currentCamera = Camera.main;
#endif

            return _currentCamera;
        }
        set => _currentCamera = value; 
    }

    [ImageEffectOpaque]
    void OnRenderImage(RenderTexture source, RenderTexture destination)
    {
        if (!raymarchingMaterial || !CurrentCamera)
        {
            Graphics.Blit(source, destination); // do nothing
            return;
        }

        ComputeBuffer voxelBuffer = voxelGrid.GetVoxelBuffer();

        if (voxelBuffer == null)
        {
            Graphics.Blit(source, destination); // do nothing
            return;
        }
            

        // pass frustum rays to shader
        raymarchingMaterial.SetMatrix(FrustumCornersESPropertyID, GetFrustumCorners(CurrentCamera));
        raymarchingMaterial.SetMatrix(CameraInvViewMatrixPropertyID, CurrentCamera.cameraToWorldMatrix);
        raymarchingMaterial.SetVector(CameraWsPropertyID, CurrentCamera.transform.position);
        raymarchingMaterial.SetVector(LightDirPropertyID, sunLight ? sunLight.transform.forward : Vector3.down);
        raymarchingMaterial.SetVector(LightColorPropertyID, sunLight ? sunLight.color : Vector3.one);
        raymarchingMaterial.SetFloat(MouseYPositionPropertyID, Mathf.Clamp01(Input.mousePosition.y / Screen.height));
        raymarchingMaterial.SetInt(MaxFloodValuePropertyID, voxelGrid.GetMaxFloodValue());
        raymarchingMaterial.SetVector(_DraineCoefficientsPropertyID, new Vector2(voxelGrid.draineG, voxelGrid.draineAlpha));
        raymarchingMaterial.SetFloat(_WaterDropletDiameterPropertyID, voxelGrid.waterDropletDiameter);

        raymarchingMaterial.SetBuffer(VoxelDataPropertyID, voxelBuffer);

        GraphicsBuffer smokeParametersBuffer = new GraphicsBuffer(GraphicsBuffer.Target.Structured, 1,
            System.Runtime.InteropServices.Marshal.SizeOf(typeof(SmokeParameters)));
        smokeParametersBuffer.SetData(new []{voxelGrid.GetSmokeParameters()});
        raymarchingMaterial.SetBuffer(SmokeParametersPropertyID, smokeParametersBuffer);

        VoxelGridParameters voxelGridParameters = voxelGrid.GetGridParameters();
        GraphicsBuffer voxelGridConstantBuffer = new GraphicsBuffer(GraphicsBuffer.Target.Constant, 1,
            System.Runtime.InteropServices.Marshal.SizeOf(typeof(VoxelGridParameters)));
        voxelGridConstantBuffer.SetData(new VoxelGridParameters[] {voxelGridParameters});
        
        raymarchingMaterial.SetConstantBuffer(VoxelGridParametersPropertyID, voxelGridConstantBuffer, 0,
            System.Runtime.InteropServices.Marshal.SizeOf(typeof(VoxelGridParameters)));

        CustomGraphicsBlit(source, new []{smokeColorTextureLowRes, smokeMaskTextureLowRes}, raymarchingMaterial, 0); // use given effect shader as image effect
        
        Graphics.Blit(smokeColorTextureLowRes, smokeColorTextureFullRes);
        Graphics.Blit(smokeMaskTextureLowRes, smokeMaskTextureFullRes);
        
        blendingMaterial.SetTexture("_MainTex", source);
        blendingMaterial.SetTexture("_SecondTex", smokeColorTextureFullRes);
        blendingMaterial.SetTexture("_MaskTex", smokeMaskTextureFullRes);
        Graphics.Blit(source, destination, blendingMaterial);
        
        smokeParametersBuffer.Dispose();
        voxelGridConstantBuffer.Dispose();
    }

    /// \brief Stores the normalized rays representing the camera frustum in a 4x4 matrix.  Each row is a vector.
    /// 
    /// The following rays are stored in each row (in eyespace, not worldspace):
    /// Top Left corner:     row=0
    /// Top Right corner:    row=1
    /// Bottom Right corner: row=2
    /// Bottom Left corner:  row=3
    private Matrix4x4 GetFrustumCorners(Camera cam)
    {
        float camFov = cam.fieldOfView;
        float camAspect = cam.aspect;

        Matrix4x4 frustumCorners = Matrix4x4.identity;

        float fovWHalf = camFov * 0.5f;

        float tan_fov = Mathf.Tan(fovWHalf * Mathf.Deg2Rad);

        Vector3 toRight = Vector3.right * tan_fov * camAspect;
        Vector3 toTop = Vector3.up * tan_fov;

        Vector3 topLeft = (-Vector3.forward - toRight + toTop);
        Vector3 topRight = (-Vector3.forward + toRight + toTop);
        Vector3 bottomRight = (-Vector3.forward + toRight - toTop);
        Vector3 bottomLeft = (-Vector3.forward - toRight - toTop);

        frustumCorners.SetRow(0, topLeft);
        frustumCorners.SetRow(1, topRight);
        frustumCorners.SetRow(2, bottomRight);
        frustumCorners.SetRow(3, bottomLeft);

        return frustumCorners;
    }
    
    /// \brief Custom version of Graphics.Blit that encodes frustum corner indices into the input vertices.
    /// 
    /// In a shader you can expect the following frustum cornder index information to get passed to the z coordinate:
    /// Top Left vertex:     z=0, u=0, v=0
    /// Top Right vertex:    z=1, u=1, v=0
    /// Bottom Right vertex: z=2, u=1, v=1
    /// Bottom Left vertex:  z=3, u=1, v=0
    /// 
    /// \warning You may need to account for flipped UVs on DirectX machines due to differing UV semantics
    ///          between OpenGL and DirectX.  Use the shader define UNITY_UV_STARTS_AT_TOP to account for this.
    static void CustomGraphicsBlit(RenderTexture source, RenderTexture[] destinations, Material fxMaterial, int passNr)
    {
        RenderBuffer[] rb = new RenderBuffer[destinations.Length];
        
        // Set targets needs the color buffers so make a array from
        // each textures buffer.
        for(int i = 0; i < destinations.Length; i++)
            rb[i] = destinations[i].colorBuffer;
        
        //Set the targets to render into.
        //Will use the depth buffer of the
        //first render texture provided.
        Graphics.SetRenderTarget(rb, destinations[0].depthBuffer);

        GL.Clear(true, true, Color.clear);

        fxMaterial.SetTexture("_MainTex", source);

        GL.PushMatrix();
        GL.LoadOrtho(); // Note: z value of vertices don't make a difference because we are using ortho projection

        fxMaterial.SetPass(passNr);

        GL.Begin(GL.QUADS);

        // Here, GL.MultitexCoord2(0, x, y) assigns the value (x, y) to the TEXCOORD0 slot in the shader.
        // GL.Vertex3(x,y,z) queues up a vertex at position (x, y, z) to be drawn.  Note that we are storing
        // our own custom frustum information in the z coordinate.
        GL.MultiTexCoord2(0, 0.0f, 0.0f);
        GL.Vertex3(0.0f, 0.0f, 3.0f); // BL

        GL.MultiTexCoord2(0, 1.0f, 0.0f);
        GL.Vertex3(1.0f, 0.0f, 2.0f); // BR

        GL.MultiTexCoord2(0, 1.0f, 1.0f);
        GL.Vertex3(1.0f, 1.0f, 1.0f); // TR

        GL.MultiTexCoord2(0, 0.0f, 1.0f);
        GL.Vertex3(0.0f, 1.0f, 0.0f); // TL
    
        GL.End();
        GL.PopMatrix();
    }

#if UNITY_EDITOR
    private Camera GetSceneEditorCamera()
    {
        return UnityEditor.EditorWindow.GetWindow<UnityEditor.SceneView>().camera;
    }
#endif
}