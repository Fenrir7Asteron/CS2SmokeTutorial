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
            #include "AutoLight.cginc"

            // Raymarching constants
            #define MAX_STEP_COUNT 512
            #define MAX_SHADOW_STEP_COUNT 32
            #define STEP_SIZE 0.04
            #define SHADOW_STEP_SIZE 0.04
            #define ALPHA_THRESHOLD 0.01

            // Noise constants
            #define WORLEY_NOISE_SLICES_COUNT 256
            #define WORLEY_NOISE_FREQUENCY 0.7
            #define NOISE_AMPLITUDE 2

            #define PI 3.141592

            // Provided by our script
            uniform float4x4 _FrustumCornersES;
            uniform sampler2D _MainTex;
            uniform float4 _SmokeAlbedo;
            uniform float4 _MainTex_TexelSize;
            uniform float4x4 _CameraInvViewMatrix;
            uniform float3 _CameraWS;
            uniform float3 _LightDir;
            uniform float4 _LightColor0;
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
                float2 density; // x - smoke density, y - shadow density
                float2 coefficients; // x - absorption coefficient, y - scattering coefficient
            };

            CBUFFER_START(_voxelGridParameters)
                int4 _voxelCounts; // x stores total voxelCount, yzw store voxel counts for each xyz axis 
                float4 _voxelZoneExtents; // xyz store zone extents
                float4 _voxelZoneCenter; // xyz store zone center
                float _voxelSize;
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

            struct fragOutput
            {
                float4 col0 : COLOR0;
                float4 col1 : COLOR1;
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
                int3(0, 0, 0),
                int3(1, 0, 0),
                int3(0, 0, 1),
                int3(1, 0, 1),
                int3(0, 1, 0),
                int3(1, 1, 0),
                int3(0, 1, 1),
                int3(1, 1, 1),
            };

            float CalculateWorleyNoise(const float3 samplePosition, const float3 spawnPosition, const float3 smokeRadius)
            {
                const float timeOffset = _Time.x;
                float3 n = (samplePosition - spawnPosition) / smokeRadius;
                
                // float u = atan2(n.x, n.z) / (2*PI) + 0.5;
                // float v = n.y * 0.5 + 0.5;
                // float2 uv = float2(u, v + timeOffset);
                
                const float slices = WORLEY_NOISE_SLICES_COUNT; // number of layers of the 3d texture
                const float freq = WORLEY_NOISE_FREQUENCY;
                
                float3 noiseUVW = float3(n * slices) / slices
                    - float3(timeOffset * 1.5, timeOffset, 0.0)
                ;
                float pfbm= lerp(1.0, perlinfbm(noiseUVW, 4.0, 7), 0.5);
                pfbm = abs(pfbm * 2.0 - 1.0); // billowy perlin noise
                //return pfbm;
                
                float worleyNoise = worleyFbm(noiseUVW, freq);
                worleyNoise = remap(pfbm, 0.0, 1.0, worleyNoise, 1.0); // perlin-worley

                worleyNoise = saturate(worleyNoise) * NOISE_AMPLITUDE;

                return NOISE_AMPLITUDE - worleyNoise;
            }

            float ease_smoke_progress(float x) {
                if (x < 0.5)
                {
                    return 2.0 * pow(x, 2);   
                }

                return 1.0 - (0.025 / (0.2 * pow(x, 2)));
            }

            // This is the distance field function.  The distance field represents the closest distance to the surface
            // of any object we put in the scene.  If the given point (point p) is inside of an object, we return a
            // negative answer.
            float SampleSmokeDensity(in const float3 worldPosition, in const float3 spawnPosition,
                in float distanceFromSmokeOrigin, in const float3 smokeRadius)
            {
                const int voxelIdx =
                    CalculateVoxelIdxByWorldPosition(worldPosition,
                        _voxelZoneCenter.xyz, _voxelZoneExtents.xyz, _voxelCounts.yzw);

                const float isValidVoxelIdx = step(0.0, (float) voxelIdx);
                
                const Cube v0 = _smokeVoxels[voxelIdx];
                const float distanceFlood = 1.0 - v0.flags.z / (float) _MaxFloodValue;
                const float distance = max(distanceFromSmokeOrigin, distanceFlood);
                const float noise = CalculateWorleyNoise(worldPosition, spawnPosition, smokeRadius);
                //const float falloff = smoothstep(0.01, ease_smoke_progress(smokeProgress), distance + noise);
                const float falloff = max(0.7, smoothstep(0.01, 1.0, (distance + noise)));
                
                const int3 voxelIdx3D = VoxelIdxToInt3(voxelIdx, _voxelCounts.yz);
                float densities[8];
                const float3 bottomLeftCorner = float3(0.0, 0.0, 0.0);
                for (int i = 0; i < 8; ++i)
                {
                    const int3 neighbourIdx3D = voxelIdx3D + offsets[i];
                    const float isInsideVoxelGrid = insideBox3D(neighbourIdx3D, bottomLeftCorner, _voxelCounts.yzw);
                
                    const uint neighbourVoxelIdx = VoxelIdx3DToVoxelIdx(neighbourIdx3D, _voxelCounts.yz);
                
                    densities[i] = _smokeVoxels[neighbourVoxelIdx].flags.x * isInsideVoxelGrid;
                }
                
                const float3 voxelPosition = v0.position;
                const float3 normalizedOffsetFromVoxelPosition = (worldPosition - voxelPosition + _voxelSize.xxx * 0.5) / (_voxelSize);
                float density = saturate(interpolate3D(densities[0], densities[1], densities[2], densities[3],
                    densities[4], densities[5], densities[6], densities[7], normalizedOffsetFromVoxelPosition));

                //density = v0.flags.x;

                return density
                    * isValidVoxelIdx 
                    * (1.0 - falloff)
                    ;
            }

            float CalculateTransmittance(float density)
            {
                return exp(-density);
            }

            float RayleighScattering(float cosThetaSquared)
            {
                return (3.0 * PI / 16.0) * (1 + cosThetaSquared);
            }

            // Raymarch along given ray and return smoke color
            float4 RaymarchCalculatePixelColor(const float3 rayOrigin, const float3 rayDirection, const float sceneDepth,
                out float alpha) {
                float4 resColor = float4(0.0, 0.0, 0.0, 0.0);
                const float3 volumeAlbedo = _SmokeAlbedo.xyz;
                
                const SmokeParameters smokeParameters = _smokeParameters[0];
                const float3 minAABB = max(smokeParameters.spawnPosition.xyz - smokeParameters.radius.xyz, _voxelZoneCenter.xyz - _voxelZoneExtents.xyz);
                const float3 maxAABB = min(smokeParameters.spawnPosition.xyz + smokeParameters.radius.xyz, _voxelZoneCenter.xyz + _voxelZoneExtents.xyz);

                const float3 boxIntersection = RayToBoxIntersection(rayOrigin, rayDirection, minAABB, maxAABB);

                const float isIntersectingVoxelGrid = boxIntersection.x;
                const float maxDistance = boxIntersection.y;
                const float minDistance = boxIntersection.z;

                const float3 rayOriginInsideGrid = rayOrigin + rayDirection * minDistance * step(0.0, minDistance);

                if (isIntersectingVoxelGrid < 1.0)
                {
                    alpha = 1.0;
                    return resColor;
                }

                float3 samplePosition = rayOriginInsideGrid;
                float distanceTraveled = length(samplePosition - rayOrigin);
                const float3 sampleOffset = rayDirection * STEP_SIZE;
                const float3 shadowSampleOffset = -_LightDir * SHADOW_STEP_SIZE;
                const float densityMultiplier = smokeParameters.density.x;
                const float shadowDensityMultiplier = smokeParameters.density.y;
                const float absorptionCoefficient = smokeParameters.coefficients.x;
                const float scatteringCoefficient = smokeParameters.coefficients.y;
                const float extinctionCoefficient = absorptionCoefficient + scatteringCoefficient;
                

                const float3 smokeSpawnPosition = smokeParameters.spawnPosition.xyz;
                const float2 spawnTime = smokeParameters.spawnTime.x;
                //const float smokeProgress = (clamp((_Time.y - spawnTime.y) / spawnTime.x, 0.01, 1.0));
                const float3 smokeRadius = smokeParameters.radius.xyz;

                const float rayToSunAngle = dot(rayDirection, _LightDir);
                const float rayToSunAngleSquared = rayToSunAngle * rayToSunAngle;
                const float phase = RayleighScattering(rayToSunAngleSquared);

                float thickness = 0.0f;
                float totalDensity = 0.0f;

                for (int sampleIdx = 0; sampleIdx < MAX_STEP_COUNT; ++sampleIdx) {
                    const float distanceFromSmokeOrigin =
                        saturate(length((samplePosition - smokeSpawnPosition) / smokeRadius));
                    
                    const float depthTest = step(distanceTraveled, sceneDepth);

                    if (distanceTraveled > maxDistance || depthTest < 1.0)
                    {
                        break;
                    }
                    
                    float sampleDensity = SampleSmokeDensity(samplePosition, smokeSpawnPosition, distanceFromSmokeOrigin, smokeRadius)
                        * depthTest;

                    thickness += sampleDensity * STEP_SIZE;
                    totalDensity += sampleDensity * densityMultiplier;

                    alpha = CalculateTransmittance(thickness * totalDensity * extinctionCoefficient);

                    if (sampleDensity > 0.001)
                    {
                        float totalShadowDensity = 0.0f;
                        float3 shadowSamplePosition = samplePosition;
                        
                        for (int shadowSampleIdx = 0; shadowSampleIdx < MAX_SHADOW_STEP_COUNT; ++shadowSampleIdx) {
                            const float shadowSmokeDensity = SampleSmokeDensity(shadowSamplePosition,
                                smokeSpawnPosition, distanceFromSmokeOrigin, smokeRadius);
                            
                            totalShadowDensity += shadowSmokeDensity;
                            shadowSamplePosition += shadowSampleOffset;
                        }

                        sampleDensity = (sampleDensity * densityMultiplier * STEP_SIZE);
                        const float shadowTransmittance = CalculateTransmittance(totalShadowDensity * SHADOW_STEP_SIZE * shadowDensityMultiplier * extinctionCoefficient);

                        float4 color = _LightColor0 * shadowTransmittance * alpha * phase * scatteringCoefficient * sampleDensity;
                        resColor += color;
                    }

                    samplePosition += sampleOffset;
                    distanceTraveled += STEP_SIZE;

                    if (alpha < ALPHA_THRESHOLD)
                    {
                        alpha = 0.0;
                        break;   
                    }
                }

                // if (distanceTraveled > sceneDepth)
                // {
                //     
                // }

                return resColor;
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

            fragOutput frag (v2f i) : SV_Target
            {
                fragOutput output;
                
                float3 rayDirection = normalize(i.ray.xyz);
                const float3 rayOrigin = _CameraWS;

                float sceneDepth = LinearEyeDepth(SAMPLE_DEPTH_TEXTURE(_CameraDepthTexture, i.uv));
                sceneDepth *= length(i.ray);

                //const float3 originalColor = tex2D(_MainTex, i.uv); // Color of the scene before this shader was run

                float alpha = 1.0;
                const float4 raymarchedColor = RaymarchCalculatePixelColor(rayOrigin, rayDirection, sceneDepth, alpha);
                //float3 color = originalColor * (1.0 - raymarchedColor.a) + raymarchedColor.rgb;
                output.col0 = raymarchedColor;
                output.col1 = float4(alpha, alpha, alpha, alpha);
                return output;
            }
            ENDCG
        }
    }
}