int3 VoxelIdxToInt3(uint voxelIdx, const int2 voxelCounts)
{
    const int x = voxelIdx % voxelCounts.x;
    const int y = (voxelIdx / voxelCounts.x) % voxelCounts.y;
    const int z = voxelIdx / (voxelCounts.x * voxelCounts.y);
    return int3(x, y, z);
}

uint VoxelIdx3DToVoxelIdx(const int3 voxelIndices, const int2 voxelCounts)
{
    return voxelIndices.x + voxelIndices.y * voxelCounts.x + voxelIndices.z * voxelCounts.x * voxelCounts.y;
}

float insideBox3D(float3 v, float3 bottomLeft, float3 topRight) {
    float3 s = step(bottomLeft, v) - step(topRight, v);
    return s.x * s.y * s.z; 
}

float3 VoxelIdxToVoxelPosition(const uint voxelIdx, const int2 voxelCounts, const float3 voxelZoneCenter,
    const float3 voxelZoneExtents, const float3 voxelSize)
{
    const float xPos = voxelIdx % voxelCounts.x;
    const float yPos = (voxelIdx / voxelCounts.x) % voxelCounts.y;
    const float zPos = voxelIdx / (voxelCounts.x * voxelCounts.y);
    const float3 voxelSizeVector = voxelSize;
    return float3(xPos, yPos, zPos) * voxelSizeVector + voxelSizeVector * 0.5
        + voxelZoneCenter - voxelZoneExtents;
}

// Output x - isIntersected, y - maxDistance, z - minDistance
float3 RayToBoxIntersection(float3 ray_origin, float3 ray_dir, float3 minpos, float3 maxpos) {
  float3 inverse_dir = 1.0 / ray_dir;
  float3 tbot = inverse_dir * (minpos - ray_origin);
  float3 ttop = inverse_dir * (maxpos - ray_origin);
  float3 tmin = min(ttop, tbot);
  float3 tmax = max(ttop, tbot);
  float2 traverse = max(tmin.xx, tmin.yz);
  float traverselow = max(traverse.x, traverse.y);
  traverse = min(tmax.xx, tmax.yz);
  float traversehi = min(traverse.x, traverse.y);
  return float3(float(traversehi > max(traverselow, 0.0)), traversehi, traverselow);
}