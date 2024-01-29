using UnityEngine;
using UnityEngine.Serialization;

public class BlendImages : MonoBehaviour
{
    [SerializeField] private RenderTexture smokeColorTexture;
    [SerializeField] private RenderTexture smokeMaskTexture;
    [SerializeField] private Material blendingMaterial;
    void OnRenderImage(RenderTexture source, RenderTexture destination)
    {
        //Graphics.Blit(smokeTexture, destination);
        blendingMaterial.SetTexture("_SecondTex", smokeColorTexture);
        blendingMaterial.SetTexture("_MaskTex", smokeMaskTexture);
        Graphics.Blit(source, destination, blendingMaterial);
    }
}