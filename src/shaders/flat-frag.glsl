#version 300 es
precision highp float;

uniform vec3 u_Eye, u_Ref, u_Up;
uniform vec2 u_Dimensions;
uniform float u_Time;

in vec2 fs_Pos;
out vec4 out_Col;

// Concept: Bouncing balls with various shaders (lambert, dielectric, normals)

struct Ray {
  vec3 origin;
  vec3 direction;
};

vec3 pointOnRay(Ray r, float t) {
  return r.origin + t * r.direction;
}

// SDF for a sphere centered at the origin
float sphereSDF(vec3 p, float r) {
  return length(p) - r;
}

// Compute a surface normal
  // Source: http://jamie-wong.com/2016/07/15/ray-marching-signed-distance-functions/
vec3 surfaceNormal(vec3 p, float r) {
  float e = 0.001;
  vec3 n;
  n.x = sphereSDF(vec3(p.x + e, p.y, p.z), r) - sphereSDF(vec3(p.x - e, p.y, p.z), r);
  n.y = sphereSDF(vec3(p.x, p.y + e, p.z), r) - sphereSDF(vec3(p.x, p.y - e, p.z), r);
  n.z = sphereSDF(vec3(p.x, p.y, p.z + e), r) - sphereSDF(vec3(p.x, p.y, p.z - e), r);
  return normalize(n);
}

// Sphere marching algorithm
  // Source: http://jamie-wong.com/2016/07/15/ray-marching-signed-distance-functions/
vec4 raymarch(Ray r, const float start, const int maxIterations) {
  float depth = start;
  for (int i = 0; i < maxIterations; ++i) {
    vec3 p = pointOnRay(r, depth);
    float rad = 1.2;

    // Euclidean distance to the sphere
    float toSphere = sphereSDF(p, rad);

    // We're inside the sphere, so the ray hit it; return color of sphere
    if (abs(toSphere) <= 0.01) {
      // Lambertian shading
      p = pointOnRay(r, depth + toSphere);
      vec3 nHat = surfaceNormal(p, rad);
      vec3 lHat = normalize(vec3(1.0, 1.0, -1.0)); //normalize(vec3(5.0, 5.0, 5.0) - p);
      float intensity = 1.0;
      vec4 color = vec4(1.0, 0.0, 0.0, 1.0);
      return vec4(vec3(dot(lHat, nHat)), 1.0)  * color * intensity;
    }

    // We've yet to hit the sphere; continue marching
    depth += toSphere;
  }

  // Background color
  return vec4(0.5 * (r.direction + vec3(1.0, 1.0, 1.0)), 1.0);
}

void main() {
  // From 560's ray casting slides
  vec3 eyeToRef = u_Ref - u_Eye;
  float len = length(eyeToRef);
  vec3 uLook = normalize(eyeToRef);
  vec3 uRight = cross(uLook, u_Up);

  float vFOV = 90.0;
  float tanAlpha = tan(vFOV * 0.5);
  float aspect = u_Dimensions.x / u_Dimensions.y;
  vec3 V = u_Up * len * tanAlpha;
  vec3 H = uRight * len * aspect * tanAlpha;
  
  vec3 p = u_Ref + fs_Pos.x * H + fs_Pos.y * V;
  vec3 dir = normalize(p - u_Eye);

  // SDF Coloring
  Ray r = Ray(u_Eye, dir);
  out_Col = raymarch(r, 0.001, 20);

  // out_Col = vec4(0.5 * (dir + vec3(1.0, 1.0, 1.0)), 1.0);
  //out_Col = vec4(0.5 * (fs_Pos + vec2(1.0)), 0.5 * (sin(u_Time * 3.14159 * 0.01) + 1.0), 1.0);
}
