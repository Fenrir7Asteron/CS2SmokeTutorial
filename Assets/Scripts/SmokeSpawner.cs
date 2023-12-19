using UnityEngine;

public class SmokeSpawner : MonoBehaviour
{
    [Header("Smoke parameters")]
    [SerializeField] private Vector3 smokeEllipsoidRadii;
    [SerializeField] private float smokeSpawnDuration;
    
    [Header("Scene references")]
    [SerializeField] private Camera mainCamera;
    [SerializeField] private VoxelGrid voxelGrid;

    public void Update()
    {
        if (Input.GetMouseButtonDown(0))
        {
            Ray screenRay = mainCamera.ScreenPointToRay(new Vector3(Input.mousePosition.x, Input.mousePosition.y, 1));
            if (Physics.Raycast(screenRay.origin, screenRay. direction * 1000.0f, out RaycastHit hitInfo))
            {
                voxelGrid.SpawnSmoke(hitInfo.point, smokeEllipsoidRadii, smokeSpawnDuration);
            }
        }
    }
}