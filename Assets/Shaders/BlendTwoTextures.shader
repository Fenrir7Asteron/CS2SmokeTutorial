// Upgrade NOTE: replaced 'mul(UNITY_MATRIX_MVP,*)' with 'UnityObjectToClipPos(*)'

Shader "Hidden/RaymarchingSmoke"
{
    Properties
    {
        _MainTex ("Texture", 2D) = "white" {}
        _SecondTex ("Second Texture", 2D) = "white" {}
        _MaskTex ("Mask Texture", 2D) = "white" {}
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
            
            // Provided by our script
            uniform sampler2D _MainTex;
            uniform sampler2D _SecondTex;
            uniform sampler2D _MaskTex;

            
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
            };
            
            v2f vert (appdata v)
            {
                v2f o;
                
                o.pos = UnityObjectToClipPos(v.vertex);
                o.uv = v.uv.xy;
                
                return o;
            }

            fixed4 frag (v2f i) : SV_Target
            {
                const float4 color1 = tex2D(_MainTex, i.uv); // Color of the scene before this shader was run                
                const float4 color2 = tex2D(_SecondTex, i.uv); // Color of the scene before this shader was run
                const float alpha = 1.0 -  saturate(tex2D(_MaskTex, i.uv).r);
                
                // float x = 1.0 - alpha;
                // float u = x * x * (3.0 - 2.0 * x ); // cubic curve
                const float3 colorf = color1.rgb * alpha + color2.rgb;
                const float alphaf = alpha + color2.a * (1.0 - alpha);
                return float4(colorf, alphaf);
            }
            ENDCG
        }
    }
}