Shader "Instanced/Simple" {
  //show values to edit in inspector
  Properties{
  }

  SubShader{
    //the material is completely non-transparent and is rendered at the same time as the other opaque geometry
    Tags{ "RenderType"="Opaque" "Queue"="Geometry"}

    Pass{
      CGPROGRAM
// Upgrade NOTE: excluded shader from DX11, OpenGL ES 2.0 because it uses unsized arrays
//#pragma exclude_renderers d3d11 gles
      //allow instancing
      #pragma multi_compile_instancing
      #pragma target 4.5
      #pragma instancing_options procedural:ConfigureProcedural

      //shader functions
      #pragma vertex vert
      #pragma fragment frag

      //use unity shader library
      #include "UnityCG.cginc"

      //per vertex data that comes from the model/parameters
      struct appdata{
        float4 vertex : POSITION;
        UNITY_VERTEX_INPUT_INSTANCE_ID
      };

      //per vertex data that gets passed from the vertex to the fragment function
      struct v2f{
        float4 position : SV_POSITION;
        float4 color : TEXCOORD0;
        UNITY_VERTEX_INPUT_INSTANCE_ID
      };

      struct Cube
      {
          float4 position;
          float4 color;
          int3 flags; // x - isVisible, y - is occupied by static mesh, z - flood fill remaining power
      };

      StructuredBuffer<Cube> _voxels;
      float _voxelSize;

      void ConfigureProcedural()
      {
        float3 position = _voxels[unity_InstanceID].position.xyz;
        float isVisible = float(_voxels[unity_InstanceID].flags.x);

        unity_ObjectToWorld = 0.0;
				unity_ObjectToWorld._m03_m13_m23_m33 = float4(position, 1.0);
				unity_ObjectToWorld._m00_m11_m22 = _voxelSize * step(0.1, isVisible);
      }

      v2f vert(appdata v){
        v2f o;

        //setup instance id
        UNITY_SETUP_INSTANCE_ID(v);
        UNITY_TRANSFER_INSTANCE_ID(v, o);

        //calculate the position in clip space to render the object
        o.position = UnityObjectToClipPos(v.vertex);
        o.color = _voxels[unity_InstanceID].color;
        
        return o;
      }

      fixed4 frag(v2f i) : SV_TARGET{
          //setup instance id
          UNITY_SETUP_INSTANCE_ID(i);
        return i.color;
      }

      ENDCG
    }
  }
}