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
            Transform cameraTransform = mainCamera.transform;
            Vector3 origin = cameraTransform.position;
            Vector3 direction = cameraTransform.forward;
            Debug.Log(origin + " " + direction);
            if (Physics.Raycast(origin, direction, out RaycastHit hitInfo, 1000.0f, LayerMask.GetMask("Default")))
            {
                Debug.Log(hitInfo.point);
                voxelGrid.SpawnSmoke(hitInfo.point, smokeEllipsoidRadii, smokeSpawnDuration);
            }
        }
    }
}