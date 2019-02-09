#version 300 es
precision highp float;

uniform vec3 u_Eye, u_Ref, u_Up;
uniform vec2 u_Dimensions;
uniform float u_Time;

in vec2 fs_Pos;
out vec4 out_Col;

// Concept: Bouncing balls with various shaders (lambert, dielectric, normals)
  // Requirements:
    // BVH
    // Intersection, subtraction, smooth blend -- combine all in one object
    // Easing function x 1 -- use in background shading
    // Noise -- use in background shading
    // dat.GUI x 2 (Number of balls, coeff of restitution)

struct Ball {
  vec4 baseColor;
  vec3 center;
  float size;
  int reflectionModel;
};

struct Ray {
  vec3 origin;
  vec3 direction;
};

float hash1D(float x) {
  return fract(sin(x * 1324623.13274867) * 356211.532);
}

// From Mariano's github
float hash3D(vec3 x) {
	float i = dot(x, vec3(123.4031, 46.5244876, 91.106168));
	return fract(sin(i * 7.13) * 268573.103291);
}

// Randomize the shading model for each ball
int randomShadingModel(vec3 n) {
  float r = hash3D(n);
  if (r < 0.33) {
    return 0;
  }
  else if (r < 0.66) {
    return 1;
  }
  else {
    return 2;
  }
}

// March t units along ray r
vec3 pointOnRay(Ray r, float t) {
  return r.origin + t * r.direction;
}

// SDF for a sphere centered at c
  // Source: http://www.michaelwalczyk.com/blog/2017/5/25/ray-marching
float sphereSDF(vec3 p, vec3 c, float r) {
  return length(p - c) - r;
}

float floorSDF(vec3 p, vec4 n) {
  return dot(p, n.xyz) + n.w;
}

// Compute log base <base> of <arg>
float logBase(float arg, float base) {
  return log(arg) / log(base);
}

// Explicit function for a bounce; xi must be in [0, 1]
  // Source: https://physics.stackexchange.com/questions/245791/explicit-function-for-bouncing-ball
float bounce(float t, float xi, float v0, float y0) {
  // Slow down time to make bouncing last longer
  t /= 24.0;

  // Assume y(0) = y0, yDot(0) = v0, g = 1
  float logBaseXi = logBase(t * (xi - 1.0) / (2.0 * v0) + 1.0, xi);
  float k = floor(logBaseXi);
  float a = t - 2.0 * v0 * (pow(xi, k) - 1.0) / (xi - 1.0);
  float y = (v0 * pow(xi, k) * a - 0.5 * a * a + y0);

  // Let the balls settle near their starting position
  if (y - y0 < 0.001) {
    return y0;
  }
  return y;
}

// Compute a surface normal
  // Source: http://jamie-wong.com/2016/07/15/ray-marching-signed-distance-functions/
vec3 surfaceNormal(vec3 p, vec3 c, float r, vec4 floorNorm, int sdfToCall) {
  float e = 0.001;
  vec3 n;
  if (sdfToCall == 0) {
    n.x = sphereSDF(vec3(p.x + e, p.y, p.z), c, r) - sphereSDF(vec3(p.x - e, p.y, p.z), c, r);
    n.y = sphereSDF(vec3(p.x, p.y + e, p.z), c, r) - sphereSDF(vec3(p.x, p.y - e, p.z), c, r);
    n.z = sphereSDF(vec3(p.x, p.y, p.z + e), c, r) - sphereSDF(vec3(p.x, p.y, p.z - e), c, r);
  }
  else if (sdfToCall == 1) {
    n.x = floorSDF(vec3(p.x + e, p.y, p.z), floorNorm) - floorSDF(vec3(p.x - e, p.y, p.z), floorNorm);
    n.y = floorSDF(vec3(p.x, p.y + e, p.z), floorNorm) - floorSDF(vec3(p.x, p.y - e, p.z), floorNorm);
    n.z = floorSDF(vec3(p.x, p.y, p.z + e), floorNorm) - floorSDF(vec3(p.x, p.y, p.z - e), floorNorm);
  }
  return normalize(n);
}

// From IQ; used to display multiple SDFs in one scene
float opUnion( float d1, float d2 ) {  
  return min(d1, d2);
}

// Sphere marching algorithm
  // Source: http://jamie-wong.com/2016/07/15/ray-marching-signed-distance-functions/
vec4 raymarch(Ray r, Ball balls[20], const float start, const int maxIterations, float t) {
  float depth = start;

  for (int i = 0; i < maxIterations; ++i) {
    vec3 p = pointOnRay(r, depth);

    // Find closest SDF
    float minDist = 100000000.0;
    int idx = 0;
    float shiftedY = 0.0;
    float shiftedX = 0.0;
    for (int j = 0; j < balls.length(); ++j) {
      Ball b = balls[j];

      // Offset time so the bounces are jittered
      float tOffset = hash1D(float(j)) * 100.0;

      // Compute bounce height and horizontal move
      b.center.y = bounce(t + tOffset, 0.99, 4.0, b.size);
      // Reset x using mod to make the balls appear infinite
      b.center.x = mod((t + tOffset) / 60.0, 4.0);

      float currDist = opUnion(minDist, sphereSDF(p, b.center, b.size));
      if (currDist < minDist) {
        minDist = currDist;
        idx = j;

        shiftedX = b.center.x;
        shiftedY = b.center.y;
      }
    }

    // To be passed to the surfaceNormal() function; identifies which SDF to call
    int sdfID = 0;

    Ball b = balls[idx];

    // Check distance to floor separately
    vec4 nFloor = normalize(vec4(0.0, 1.0, 0.0, 0.0));
    float toFloor = floorSDF(p, nFloor);
    if (toFloor < minDist) {
      shiftedX = 0.0;
      shiftedY = 0.0;
      minDist = toFloor;
      sdfID = 1;
    }

    // Move the ball as calculated above
    b.center.x = shiftedX;
    b.center.y = shiftedY;

    // Euclidean distance to the shape
    float toShape = minDist;

    // We're inside the shape, so the ray hit it; return color of shape + shading
    if (abs(toShape) <= 0.01) {
      // Lambert (applied to all shapes)
      vec3 nHat = surfaceNormal(p, b.center, b.size, nFloor, sdfID);
      vec3 lHat = normalize(vec3(1.0, 1.0, -1.0)); // Light vector
      float intensity = 2.0;
      vec4 color = b.baseColor; // Makes a Voronoi diagram on the floor plane!
      float dProd = max(dot(lHat, nHat), 0.0);
      color = vec4(vec3(dProd), 1.0) * color * intensity;

      if (b.reflectionModel == 1) {
        // Blinn-Phong
          // Source: https://en.wikipedia.org/wiki/Blinn%E2%80%93Phong_shading_model
        vec3 h = normalize(lHat - r.direction);
        float angle = max(dot(h, nHat), 0.0);
        float spec = pow(angle, 25.0);
        vec4 specColor = vec4(1.0);
        color += specColor * vec4(vec3(spec), 1.0) * intensity;
      }
      else if (b.reflectionModel == 2) {
        // Surface normals
        color = vec4(0.5 * (nHat + vec3(1.0)), 1.0);
        // Diversify the floor colors
        if (sdfID == 1) {
          color *= vec4(normalize(float(idx) * 2.0 * vec3(b.center.x, 1.0, b.center.z)), 1.0);
        }
      }
      return color;
    }

    // We've yet to hit the sphere; continue marching
    depth += toShape;
  }

  // We miss the object entirely; return the background color
  return vec4(0.5 * (r.direction + vec3(1.0, 1.0, 1.0)), 1.0);
}

void main() {
  // Basic ray casting (from 560 slides)
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

  // Create bouncing balls
  const int nBalls = 20;
  Ball ballArr[nBalls];
  for (float i = 0.0; i < float(nBalls); i += 1.0) {
    highp int idx = int(i);
    ballArr[idx] = Ball(vec4(hash3D(vec3(38.324 * i)) + 0.1, hash3D(vec3(10.4924 * i)) + 0.08, hash3D(vec3(102.521 * i)) + 0.1, 1.0),
                      vec3(hash3D(vec3(2.538 * i)) + 10.0 * hash3D(vec3(12.53892538 * i)), 2.0 * hash3D(vec3(235.202 * i)), 4.0 * hash3D(vec3(59.3423 * i)) - 2.0),
                      hash3D(vec3(12.53892538 * i * i + 3452.31)) + 0.2,
                      randomShadingModel(vec3(95.3829 * i)));
    ballArr[idx].center += vec3(-2.0, 0.0, -2.0) * 4.0;
    ballArr[idx].size /= 3.0;
  }

  // SDF Coloring and motion
  Ray r = Ray(u_Eye, dir);
  out_Col = raymarch(r, ballArr, 0.001, 150, u_Time);
}
