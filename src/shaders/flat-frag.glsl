#version 300 es
precision highp float;

uniform vec3 u_Eye, u_Ref, u_Up;
uniform vec2 u_Dimensions;
uniform float u_Time;
uniform float u_Balls;
uniform float u_Bounce;

in vec2 fs_Pos;
out vec4 out_Col;

// Properties of each ball defined in the scene
struct Ball {
  vec4 baseColor;
  vec3 center;
  float size;
  int reflectionModel;
};

// Properties of each ray cast into the scene
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
    // Lambert
    return 0;
  }
  else if (r < 0.58) {
    // Blinn-Phong
    return 1;
  }
  else {
    // Normals
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

// SDF for a plane; from IQ
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

// 3D noise
float noise(vec3 p) {
  vec3 bCorner = floor(p);
  vec3 inCell = fract(p);

  float bLL = hash3D(bCorner);
  float bUL = hash3D(bCorner + vec3(0.0, 0.0, 1.0));
  float bLR = hash3D(bCorner + vec3(0.0, 1.0, 0.0));
  float bUR = hash3D(bCorner + vec3(0.0, 1.0, 1.0));
  float b0 = mix(bLL, bUL, inCell.z);
  float b1 = mix(bLR, bUR, inCell.z);
  float b = mix(b0, b1, inCell.y);

  vec3 fCorner = bCorner + vec3(1.0, 0.0, 0.0);
  float fLL = hash3D(fCorner);
  float fUL = hash3D(fCorner + vec3(0.0, 0.0, 1.0));
  float fLR = hash3D(fCorner + vec3(0.0, 1.0, 0.0));
  float fUR = hash3D(fCorner + vec3(0.0, 1.0, 1.0));
  float f0 = mix(fLL, fUL, inCell.z);
  float f1 = mix(fLR, fUR, inCell.z);
  float f = mix(f0, f1, inCell.y);

  return mix(b, f, inCell.x);
}

// 5-octave FBM
float fbm(vec3 q) {
  float acc = 0.0;
  float freqScale = 2.0;
  float invScale = 1.0 / freqScale;
  float freq = 1.0;
  float amp = 1.0;

  for (int i = 0; i < 5; ++i) {
    freq *= freqScale;
    amp *= invScale;
    acc += noise(q * freq) * amp;
  }
  return acc;
}

//////// From IQ ////////
float pattern(in vec3 p) {
  vec3 q = vec3(fbm(p + vec3(0.0)),
                fbm(p + vec3(5.2, 1.3, 2.8)),
                fbm(p + vec3(1.2, 3.4, 1.2)));

  return fbm(p + 4.0 * q);
}

vec3 palette(in float t, in vec3 a, in vec3 b, in vec3 c, in vec3 d) {
    return a + b * cos(6.28318 * (c * t + d));
}

// Used to display multiple SDFs in one scene
float opUnion(float d1, float d2) {  
  return min(d1, d2);
}

float opSmoothUnion(float d1, float d2, float k) {
  float h = clamp(0.5 + 0.5 * (d2 - d1) / k, 0.0, 1.0);
  return mix(d2, d1, h) - k * h * (1.0 - h);
}

float opSubtraction(float d1, float d2) {
  return max(-d1,d2);
}

float opIntersection(float d1, float d2) {
  return max(d1,d2);
}
////////////////////////

// A more interesting cloudy background
vec4 computeBackgroundColor(Ray r, float time) {
  float t = pattern(r.direction * time);
  vec3 a = vec3(-0.992, -0.212, 0.128);
  vec3 b = vec3(1.328, 0.588, 0.318);
  vec3 c = vec3(0.958, 0.478, 0.638);
  vec3 d = vec3(-0.480, 0.908, 0.688);
  return vec4(palette(t, a, b, c, d), 1.0) * vec4(0.6, 0.2, 0.7, 0.5);
}

// Sphere marching algorithm
  // Source: http://jamie-wong.com/2016/07/15/ray-marching-signed-distance-functions/
vec4 raymarch(Ray r, Ball balls[25], const float start, const int maxIterations, float t) {
  float depth = start;

  for (int i = 0; i < maxIterations; ++i) {
    vec3 p = pointOnRay(r, depth);

    // Find closest SDF
    float minDist = 100000000.0;
    int idx = 0;
    float shiftedY = 0.0;
    float shiftedX = 0.0;
    for (int j = 0; j < 25; ++j) {
      Ball b = balls[j];

      // Offset time so the bounces are jittered
      float tOffset = hash1D(float(j)) * 100.0;

      // Compute bounce height and horizontal move
      b.center.y = bounce(t + tOffset, u_Bounce, u_Balls, -b.size * 1.5);
      // Reset x using mod to make the balls appear infinite
      b.center.x = mod((t + tOffset) / 20.0, 6.25);

      float currDist = opUnion(minDist, sphereSDF(p, b.center, b.size));
      currDist = opIntersection(opSubtraction(minDist, currDist), currDist);
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
    minDist = opSmoothUnion(minDist, 2.0 * toFloor, b.size * 3.0);
    if (toFloor < minDist) {
      shiftedX = 0.0;
      shiftedY = 0.0;
      // Make the balls appear to squish into the floor
      minDist = opUnion(toFloor, minDist);
      // minDist = opSmoothUnion(minDist, 2.0 * toFloor, 1.0);
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
        vec4 specColor = vec4(1.0); // Light color
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
      if (sdfID == 1) {
        vec4 tempColor = computeBackgroundColor(r, t * 0.0001);
        float lum = 0.2126 * tempColor.x + 0.7152 * tempColor.y + 0.0722 * tempColor.z;
        color = vec4(vec3(lum), 1.0);
      }
      return color;
    }

    // We've yet to hit anything; continue marching
    depth += toShape;
  }

  // We miss the object entirely; return the background color
  return computeBackgroundColor(r, t * 0.0003);
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
  Ball ballArr[25];
  for (float i = 0.0; i < 25.0; i += 1.0) {
    // if (i > u_Balls) { break; }
    highp int idx = int(i);
    ballArr[idx] = Ball(vec4(hash3D(vec3(38.324 * i)) + 0.1, hash3D(vec3(10.4924 * i)) + 0.08, hash3D(vec3(102.521 * i)) + 0.1, 1.0),
                        vec3(hash3D(vec3(2.538 * i)) + 10.0 * hash3D(vec3(12.53892538 * i)), 2.0 * hash3D(vec3(235.202 * i)), 8.0 * hash3D(vec3(59.3423 * i)) - 4.0),
                        hash3D(vec3(12.53892538 * i * i + 3452.31)) + 0.2,
                        randomShadingModel(vec3(95.3829 * i)));
    ballArr[idx].center += vec3(-2.0, 0.0, -2.0) * 4.0;
    ballArr[idx].size /= 3.0;
  }

  // SDF Coloring and motion
  Ray r = Ray(u_Eye, dir);
  out_Col = raymarch(r, ballArr, 0.001, 150, u_Time);
}
