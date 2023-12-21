// Upgrade NOTE: replaced 'mul(UNITY_MATRIX_MVP,*)' with 'UnityObjectToClipPos(*)'

Shader "Hidden/RaymarchingSmoke"
{
    Properties
    {
        _MainTex ("Texture", 2D) = "white" {}
    }
    SubShader
    {
        // No culling or depth
        Cull Off ZWrite Off ZTest Always

        Pass
        {
            CGPROGRAM
            #pragma vertex vert
            #pragma fragment frag
            #include "UnityCG.cginc"

            // Provided by our script
            uniform float4x4 _FrustumCornersES;
            uniform sampler2D _MainTex;
            uniform float4 _MainTex_TexelSize;
            uniform float4x4 _CameraInvViewMatrix;
            uniform float3 _CameraWS;
            uniform float3 _LightDir;
            UNITY_DECLARE_DEPTH_TEXTURE(_CameraDepthTexture);

            struct Cube
            {
                float4 position;
                float4 color;
                int3 flags; // x - isVisible, y - is occupied by static mesh, z - flood fill remaining power
            };

            CBUFFER_START(_voxelGridParameters)
                float4 _voxelSize;
                int4 _voxelCounts; // x stores total voxelCount, yzw store voxel counts for each xyz axis 
                float4 _voxelZoneExtents; // xyz store zone extents
                float4 _voxelZoneCenter; // xyz store zone center
            CBUFFER_END

            StructuredBuffer<Cube> _smokeVoxels;

            // Input to vertex shader
            struct appdata
            {
                // Remember, the z value here contains the index of _FrustumCornersES to use
                float4 vertex : POSITION;
                float2 uv : TEXCOORD0;
            };

            // Output of vertex shader / input to fragment shader
            struct v2f
            {
                float4 pos : SV_POSITION;
                float2 uv : TEXCOORD0;
                float3 ray : TEXCOORD1;
            };

            float insideBox3D(float3 v, float3 bottomLeft, float3 topRight) {
                float3 s = step(bottomLeft, v) - step(topRight, v);
                return s.x * s.y * s.z; 
            }
            
            int CalculateVoxelIdxByWorldPosition(float3 worldPosition, float3 voxelGridCenter, float3 extents,
                int3 voxelCounts)
            {
                const float3 zoneSize = extents * 2;
                float3 voxelPosition = worldPosition - voxelGridCenter + extents;
                
                voxelPosition = voxelPosition / zoneSize * voxelCounts;
                
                int3 voxelPositionInt = voxelPosition;

                const float voxelIdx = voxelPositionInt.x
                           + voxelPositionInt.y * voxelCounts.x
                           + voxelPositionInt.z * voxelCounts.x * voxelCounts.y;

                const float insideVoxelGrid = insideBox3D(worldPosition, voxelGridCenter - extents,
                    voxelGridCenter + extents);
                
                return (voxelIdx + 1) * insideVoxelGrid - 1;
            }

            // This is the distance field function.  The distance field represents the closest distance to the surface
            // of any object we put in the scene.  If the given point (point p) is inside of an object, we return a
            // negative answer.
            float SampleSmokeDensity(float3 worldPosition)
            {
                const int voxelIdx =
                    CalculateVoxelIdxByWorldPosition(worldPosition,
                        _voxelZoneCenter.xyz, _voxelZoneExtents.xyz, _voxelCounts.yzw);

                const float isValidVoxelIdx = step(0.0, (float) voxelIdx);
                return _smokeVoxels[voxelIdx].flags.x * isValidVoxelIdx;
            }

            float CalculateTransmittance(float density)
            {
                return exp(-density);
            }
            
            // Raymarch along given ray
            float Raymarch(float3 rayOrigin, float3 rayDirection, float sceneDepth) {
                const int maxStepCount = 128;
                const float stepSize = 0.25f;
                float totalDensity = 0.0f;
                
                for (int i = 0; i < maxStepCount; ++i) {
                    const float3 samplePosition = rayOrigin + rayDirection * stepSize * i;

                    const float depthTest = step(length(samplePosition - rayOrigin), sceneDepth);
                    
                    const float density = SampleSmokeDensity(samplePosition) * depthTest;   
                    totalDensity += density;

                    if (depthTest < 1.0)
                    {
                        break;
                    }
                }

                totalDensity *= stepSize;
                return CalculateTransmittance(totalDensity);
            }

            v2f vert (appdata v)
            {
                v2f o;
                
                // Index passed via custom blit function in RaymarchGeneric.cs
                half index = v.vertex.z;
                v.vertex.z = 0.1;
                
                o.pos = UnityObjectToClipPos(v.vertex);
                o.uv = v.uv.xy;
                
                #if UNITY_UV_STARTS_AT_TOP
                if (_MainTex_TexelSize.y < 0)
                    o.uv.y = 1 - o.uv.y;
                #endif

                // Get the eyespace view ray (normalized)
                o.ray = _FrustumCornersES[(int)index].xyz;

                // Transform the ray from eyespace to worldspace
                // Note: _CameraInvViewMatrix was provided by the script
                o.ray = mul(_CameraInvViewMatrix, o.ray);
                return o;
            }

            fixed4 frag (v2f i) : SV_Target
            {
                const float3 rayDirection = normalize(i.ray.xyz);
                const float3 rayOrigin = _CameraWS;

                const float sceneDepth = LinearEyeDepth(SAMPLE_DEPTH_TEXTURE(_CameraDepthTexture, i.uv));
                //const float sceneDepth = SAMPLE_DEPTH_TEXTURE(_CameraDepthTexture, i.uv) * 50;
                //return float4(sceneDepth, 0.0, 0.0, 1.0);

                const float transmittance = Raymarch(rayOrigin, rayDirection, sceneDepth);

                // Returns final color using alpha blending
                const float3 originalColor = tex2D(_MainTex,i.uv); // Color of the scene before this shader was run
                return float4(originalColor * transmittance + float3(1.0, 1.0, 1.0) * (1.0 - transmittance), 1.0);
            }
            ENDCG
        }
    }
}