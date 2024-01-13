// Upgrade NOTE: replaced 'mul(UNITY_MATRIX_MVP,*)' with 'UnityObjectToClipPos(*)'

Shader "Hidden/RaymarchingSmoke"
{
    Properties
    {
        _MainTex ("Texture", 2D) = "white" {}
        _SmokeAlbedo ("Smoke Albedo", Color) = (1.0, 1.0, 1.0, 1.0)
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
            #include "WorleyNoise.cginc"
            #include "Interpolation.cginc"
            #include "VoxelUtils.cginc"

            // Raymarching constants
            #define MAX_STEP_COUNT 256
            #define STEP_SIZE 0.05

            // Noise constants
            #define WORLEY_NOISE_SLICES_COUNT 128
            #define WORLEY_NOISE_FREQUENCY 4.0

            #define PI 3.141592

            // Provided by our script
            uniform float4x4 _FrustumCornersES;
            uniform sampler2D _MainTex;
            uniform float4 _SmokeAlbedo;
            uniform float4 _MainTex_TexelSize;
            uniform float4x4 _CameraInvViewMatrix;
            uniform float3 _CameraWS;
            uniform float3 _LightDir;
            uniform float3 _LightColor;
            uniform float _MouseYPosition;
            uniform int _MaxFloodValue;
            UNITY_DECLARE_DEPTH_TEXTURE(_CameraDepthTexture);

            struct Cube
            {
                float4 position;
                float4 color;
                int3 flags; // x - isVisible, y - is occupied by static mesh, z - flood fill remaining power
            };

            struct SmokeParameters
            {
                float4 spawnPosition; // xyz is spawn position
                float4 radius; // xyz are ellipsoid radii
                float2 spawnTime; // x - spawn duration, y - spawn start timestamp 
            };

            CBUFFER_START(_voxelGridParameters)
                float4 _voxelSize;
                int4 _voxelCounts; // x stores total voxelCount, yzw store voxel counts for each xyz axis 
                float4 _voxelZoneExtents; // xyz store zone extents
                float4 _voxelZoneCenter; // xyz store zone center
            CBUFFER_END

            StructuredBuffer<Cube> _smokeVoxels;
            StructuredBuffer<SmokeParameters> _smokeParameters;

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
            
            int CalculateVoxelIdxByWorldPosition(float3 worldPosition, float3 voxelGridCenter, float3 extents,
                const int3 voxelCounts)
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

            static int3 offsets[8] = {
                int3(-1, -1, -1),
                int3(1, -1, -1),
                int3(-1, -1, 1),
                int3(1, -1, 1),
                int3(-1, 1, -1),
                int3(1, 1, -1),
                int3(-1, 1, 1),
                int3(1, 1, 1),
            };

            // This is the distance field function.  The distance field represents the closest distance to the surface
            // of any object we put in the scene.  If the given point (point p) is inside of an object, we return a
            // negative answer.
            float SampleSmokeDensity(float3 worldPosition)
            {
                const int voxelIdx =
                    CalculateVoxelIdxByWorldPosition(worldPosition,
                        _voxelZoneCenter.xyz, _voxelZoneExtents.xyz, _voxelCounts.yzw);

                const float isValidVoxelIdx = step(0.0, (float) voxelIdx);
                
                const Cube v0 = _smokeVoxels[voxelIdx];

                const SmokeParameters smokeParameters = _smokeParameters[0];
                const float distanceEuclidian = length((worldPosition - smokeParameters.spawnPosition.xyz) / smokeParameters.radius.xyz);
                const float distanceFlood = 1.0 - v0.flags.z / (float) _MaxFloodValue;
                const float distance = max(distanceEuclidian, distanceFlood);
                const float falloff = smoothstep(0.05, 1.0, distance);

                // const int3 voxelIdx3D = VoxelIdxToInt3(voxelIdx, _voxelCounts.yz);
                // int densities[8];
                // for (int i = 0; i < 8; ++i)
                // {
                //     const int3 neighbourIdx3D = voxelIdx3D + offsets[i];
                //     const float isInsideVoxelGrid = insideBox3D(neighbourIdx3D, float3(0.0, 0.0, 0.0), _voxelCounts.yzw);
// 
                //     const uint neighbourVoxelIdx = VoxelIdx3DToVoxelIdx(neighbourIdx3D, _voxelCounts.yz);
// 
                //     densities[i] = _smokeVoxels[neighbourVoxelIdx].flags.x * isInsideVoxelGrid;
                // }

                const float density = v0.flags.x;

                //const int3 referenceIdx3D = voxelIdx3D + offsets[0];
                //const uint referenceVoxelIdx = VoxelIdx3DToVoxelIdx(referenceIdx3D, _voxelCounts.yz);
                //const float3 voxelPosition = VoxelIdxToVoxelPosition(referenceVoxelIdx, _voxelCounts.yz, _voxelZoneCenter.xyz, _voxelZoneExtents.xyz, _voxelSize.xxx);
                //const float3 normalizedOffsetFromVoxelPosition = (worldPosition - voxelPosition) / (_voxelSize * 2.0).xxx;
                //const float density = interpolate3D(densities[0], densities[1], densities[2], densities[3],
                //    densities[4], densities[5], densities[6], densities[7], normalizedOffsetFromVoxelPosition);

                return density
                    * isValidVoxelIdx 
                    * (1.0 - falloff)
                    ;
            }

            float CalculateTransmittance(float density)
            {
                return exp(-density);
            }

            // Raymarch along given ray and return transmittance
            float RaymarchCalculateTransmittance(float3 rayOrigin, float3 rayDirection) {
                const int maxStepCount = 32;
                const float stepSize = STEP_SIZE;
                float totalDensity = 0.0f;
                
                for (int i = 0; i < maxStepCount; ++i) {
                    const float3 samplePosition = rayOrigin + rayDirection * stepSize * i;

                    const float density = SampleSmokeDensity(samplePosition);   
                    totalDensity += density;
                }

                totalDensity *= stepSize;
                return CalculateTransmittance(totalDensity);
            }

            // Raymarch along given ray and return smoke color
            float3 RaymarchCalculatePixelColor(const float3 rayOrigin, const float3 rayDirection, const float sceneDepth, const float3 originalColor) {
                const int maxStepCount = MAX_STEP_COUNT;
                const float stepSize = STEP_SIZE;
                float totalDensity = 0.0f;
                float3 color3 = originalColor;
                const float3 volumeAlbedo = _SmokeAlbedo.xyz;
                float baseTransmittance = 1.0;
                float previousBaseTransmittance = 1.0;
                const float3 ambientColor = unity_AmbientGround;
                
                int sampleIdx = 0;

                const SmokeParameters smokeParameters = _smokeParameters[0];
                float3 minAABB = max(smokeParameters.spawnPosition.xyz - smokeParameters.radius.xyz, _voxelZoneCenter.xyz - _voxelZoneExtents.xyz);
                float3 maxAABB = min(smokeParameters.spawnPosition.xyz + smokeParameters.radius.xyz, _voxelZoneCenter.xyz + _voxelZoneExtents.xyz);

                const float3 boxIntersection = RayToBoxIntersection(rayOrigin, rayDirection, minAABB, maxAABB);

                const float isIntersectingVoxelGrid = boxIntersection.x;
                const float maxDistance = boxIntersection.y;
                const float minDistance = boxIntersection.z;

                const float3 rayOriginInsideGrid = rayOrigin + rayDirection * minDistance * step(0.0, minDistance);

                if (boxIntersection.x < 1.0)
                {
                    return originalColor;
                }

                for (int sampleIdx = 0; sampleIdx < maxStepCount; ++sampleIdx) {
                    const float3 samplePosition = rayOriginInsideGrid + rayDirection * stepSize * sampleIdx;
                    float distance = length(samplePosition - rayOrigin);

                    float4 worley = float4(0.0, 0.0, 0.0, 0.0);

                    const float slices = WORLEY_NOISE_SLICES_COUNT; // number of layers of the 3d texture
                    const float freq = WORLEY_NOISE_FREQUENCY;
                    
                    float3 n = normalize(samplePosition - smokeParameters.spawnPosition);
                    float u = atan2(n.x, n.z) / (2*PI) + 0.5;
                    float v = n.y * 0.5 + 0.5;
                    float2 uv = float2(u, v);

                    float pfbm= lerp(1.0, perlinfbm(float3(uv, floor(_MouseYPosition*slices)/slices), 4.0, 7), 0.5);
                    pfbm = abs(pfbm * 2.0 - 1.0); // billowy perlin noise
                    
                    worley.g += worleyFbm(float3(uv, floor(_MouseYPosition*slices)/slices), freq);
                    worley.b += worleyFbm(float3(uv, floor(_MouseYPosition*slices)/slices), freq*2.0);
                    worley.a += worleyFbm(float3(uv, floor(_MouseYPosition*slices)/slices), freq*4.0);
                    worley.r += remap(pfbm, 0.0, 1.0, worley.g, 1.0); // perlin-worley

                    distance += worley.r * 2.5;

                    const float depthTest = step(distance, sceneDepth);

                    if (distance > maxDistance || depthTest < 1.0)
                    {
                        break;
                    }

                    const float density = SampleSmokeDensity(samplePosition) * depthTest;   
                    totalDensity += density * stepSize;

                    baseTransmittance *= CalculateTransmittance(totalDensity);
                    const float stepTransmittance = previousBaseTransmittance - baseTransmittance;

                    //const float sunTransmittance = RaymarchCalculateTransmittance(samplePosition, -_LightDir);
                    //color3 += sunTransmittance * _LightColor * volumeAlbedo * stepTransmittance;
                    color3 += stepTransmittance * volumeAlbedo * ambientColor;

                    previousBaseTransmittance = baseTransmittance;
                }

                return color3;
            }

            v2f vert (appdata v)
            {
                v2f o;
                
                // Index passed via custom blit function in RaymarchGeneric.cs
                half index = v.vertex.z;
                v.vertex.z = 0.0;
                
                o.pos = UnityObjectToClipPos(v.vertex);
                o.uv = v.uv.xy;
                
                #if UNITY_UV_STARTS_AT_TOP
                if (_MainTex_TexelSize.y < 0)
                    o.uv.y = 1 - o.uv.y;
                #endif

                // Get the eyespace view ray (normalized)
                o.ray = _FrustumCornersES[(int)index].xyz;
                o.ray /= abs(o.ray.z);

                // Transform the ray from eyespace to worldspace
                // Note: _CameraInvViewMatrix was provided by the script
                o.ray = mul(_CameraInvViewMatrix, o.ray);
                
                return o;
            }

            fixed4 frag (v2f i) : SV_Target
            {
                float3 rayDirection = normalize(i.ray.xyz);
                const float3 rayOrigin = _CameraWS;

                //return float4(rayDirection, 1.0);

                

                float sceneDepth = LinearEyeDepth(SAMPLE_DEPTH_TEXTURE(_CameraDepthTexture, i.uv));
                sceneDepth *= length(i.ray);

                const float3 originalColor = tex2D(_MainTex, i.uv); // Color of the scene before this shader was run

                const float3 raymarchedColor = RaymarchCalculatePixelColor(rayOrigin, rayDirection, sceneDepth, originalColor);
                return float4(raymarchedColor, 1.0);
            }
            ENDCG
        }
    }
}