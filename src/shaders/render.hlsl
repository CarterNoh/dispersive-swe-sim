// --- RENDER CONSTANTS ---
cbuffer RenderParams : register(b0) {
    matrix viewProjection; // Camera matrix
    float gridSizeX;
    float gridSizeY;
    float cellSize;
};

// --- INPUT TEXTURE (SRV) ---
Texture2D<float> HeightMap : register(t0);
Texture2D<float> DispXMap  : register(t1);
Texture2D<float> DispYMap  : register(t2);
SamplerState PointSampler  : register(s0);

struct Vertex {
    float2 position : POSITION; // The flat X, Z grid coordinates
};

struct Pixel {
    float4 position : SV_POSITION; // 3D screen position
    float3 worldPos : TEXCOORD0;   // 3D world coordinate
    float  height   : TEXCOORD1;   // Pass height to pixel shader for coloring
};

// --- VERTEX SHADER ---
Pixel VSMain(Vertex input) {
    Pixel output;
    
    // Convert flat vertex position to texture UV coordinates (0.0 to 1.0)
    float2 uv = input.position / float2(gridSizeX - 1, gridSizeY - 1);
    float h = HeightMap.SampleLevel(PointSampler, uv, 0);
    float dx = DispXMap.SampleLevel(PointSampler, uv, 0);
    float dy = DispYMap.SampleLevel(PointSampler, uv, 0);
    float scale = 1.f;
    // float3 wPos = float3(input.position.x, h * scale, input.position.y);
    float3 wPos = float3(input.position.x + dx, h * scale, input.position.y + dy);
    
    // Project to the screen
    output.position = mul(float4(wPos, 1.0f), viewProjection);
    output.worldPos = wPos;
    output.height = h;
    
    return output;
}

// --- PIXEL SHADER ---
float4 PSMain(Pixel input) : SV_TARGET {
    // Simple water color gradient based on height
    // float3 deepWater = float3(0.0f, 0.2f, 0.5f);
    // float3 shallowWater = float3(0.0f, 0.6f, 0.8f);
    // float3 foam = float3(1.0f, 1.0f, 1.0f);
    
    // float3 color = lerp(deepWater, shallowWater, saturate(input.height / 5.0f));
    // if (input.height > 8.0f) color = lerp(color, foam, saturate((input.height - 8.0f) / 2.0f));

    // Normalize the height from 0.0 to 1.0
    // Adjust hMax depending on how high your water level / waves get!
    float hMin = -5.f;
    float hMax = 5.f; // Adjust based on your wave heights!
    float t = saturate((input.height - hMin) / (hMax - hMin));

    // "viterbi" color scheme
    float3 v0 = float3(0.267f, 0.004f, 0.329f); // Dark Purple
    float3 v1 = float3(0.231f, 0.322f, 0.545f); // Blue
    float3 v2 = float3(0.129f, 0.569f, 0.549f); // Teal
    float3 v3 = float3(0.365f, 0.788f, 0.388f); // Green
    float3 v4 = float3(0.992f, 0.906f, 0.145f); // Yellow

    float3 baseColor;
    if      (t < 0.25f) baseColor = lerp(v0, v1, t / 0.25f);
    else if (t < 0.50f) baseColor = lerp(v1, v2, (t - 0.25f) / 0.25f);
    else if (t < 0.75f) baseColor = lerp(v2, v3, (t - 0.50f) / 0.25f);
    else                baseColor = lerp(v3, v4, (t - 0.75f) / 0.25f);

    // ddx and ddy tell us how the 3D position changes from this pixel to the next
    float3 dx = ddx(input.worldPos);
    float3 dy = ddy(input.worldPos);
    
    // Cross product gives us the vector pointing straight out of the triangle
    // (If the lighting looks "inside out", swap dx and dy here)
    float3 normal = normalize(cross(dx, dy));

    // Apply Lighting
    float3 lightDir = normalize(float3(-1.0f, -1.5f, 1.0f)); 
    float diffuse = saturate(dot(normal, -lightDir));
    float ambient = 0.4f; 
    float3 finalColor = baseColor * (ambient + (diffuse * 0.7f));

    return float4(finalColor, 1.0f);
}