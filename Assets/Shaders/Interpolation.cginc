float interpolate1D(float v1, float v2, float x){
    return abs(v1*(1-x) + v2*x);
}

float interpolate2D(float v1, float v2, float v3, float v4, float2 p){

    float s = interpolate1D(v1, v2, p.x);
    float t = interpolate1D(v3, v4, p.x);
    return interpolate1D(s, t, p.y);
}

float interpolate3D(float v1, float v2, float v3, float v4,
    float v5, float v6, float v7, float v8, float3 p)
{
    float s = interpolate2D(v1, v2, v3, v4, p.xz);
    float t = interpolate2D(v5, v6, v7, v8, p.xz);
    return interpolate1D(s, t, p.y);
}