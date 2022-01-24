#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include <stdio.h>
#include <math.h>

// Vector math helpers
//

/**
 * Computes distance between two points in 3D
 * 
 * @param vec1 vector 1
 * @param vec2 vector 2
 * @return result distance
 */
float pointDistance3D(float vec1[3], float vec2[3]){
    return sqrt(pow(vec2[0]-vec1[0], 2)+pow(vec2[1]-vec1[1], 2)+pow(vec2[2]-vec1[2], 2));
}

/**
 * Copies vector 1 into vector 2
 * 
 * @param vec2 vector to be copied into
 * @param vec1 vector to copy from
 * @return result vector
 */
float * vectorCpy(float vec2[3], float vec1[3]){
    vec2[0] = vec1[0];
    vec2[1] = vec1[1];
    vec2[2] = vec1[2];
    return vec2;
}

/**
 * Computes vector from adding vec2 to vec1
 * 
 * @param vec3 output vector
 * @param vec1 left operand
 * @param vec2 right operand
 * @return result vector
 */
float * vectorAdd(float vec3[3], float vec1[3], float vec2[3]){
    vec3[0] = vec1[0] + vec2[0];
    vec3[1] = vec1[1] + vec2[1];
    vec3[2] = vec1[2] + vec2[2];
    return vec3;
}

/**
 * Computes vector from subbing vec2 from vec1
 * 
 * @param vec3 output vector
 * @param vec1 left operand
 * @param vec2 right operand
 * @return result vector
 */
float * vectorSub(float vec3[3], float vec1[3], float vec2[3]){
    vec3[0] = vec1[0] - vec2[0];
    vec3[1] = vec1[1] - vec2[1];
    vec3[2] = vec1[2] - vec2[2];
    return vec3;
}

/**
 * Normalizes vector with 3 dimensions
 * 
 * @param vec pixel location vector
 * @return normalized vector
 */
float * normalize(float vec[3]){
    float magnitude = sqrt(pow(vec[0],2.0) + pow(vec[1], 2.0) + pow(vec[2],2.0));
    vec[0] = vec[0] / magnitude;
    vec[1] = vec[1] / magnitude;
    vec[2] = vec[2] / magnitude;
    return vec;
}

/**
 * Scale vector by multiplying with scalar
 * 
 * @param vec vector
 * @param sc scalar
 * @return scaled vector
 */
float * scale(float vec[3], float sc){
    vec[0] = vec[0] * sc;
    vec[1] = vec[1] * sc;
    vec[2] = vec[2] * sc;
    return vec;
}

/**
 * fidget vector in direction by given scale
 * 
 * @param vec vector to fidget
 * @param d direction to fidget in
 * @param sc scale to fidget by
 * @return moved vector
 */
float * fidget(float vec[3], float d[3], float sc){
    vec[0] += sc*d[0];
    vec[1] += sc*d[1];
    vec[2] += sc*d[2];
    return vec;
}

/**
 * Computes new vector from ray equation
 * 
 * @param vec3 output vector
 * @param e start vector
 * @param d direction vector
 * @param t time value
 * @return resulting vector
 */
float * rayEquation(float vec3[3], float e[3], float d[3], float t){
    return vectorAdd(vec3, e, scale(d,t));
}

/**
 * Computes dot product of two vectors
 * 
 * @param vec1 vector 1
 * @param vec2 vector 2
 * @return dot product of input vectors
 */
float dotProduct(float vec1[3], float vec2[3]){
    return vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2];
}

/**
 * Clamp float between two extremes
 * 
 * @param num number to clamp
 * @param min min num
 * @param max max num
 * @return clamped float
 */
float clamp(float num, float min, float max){
    if (num < min) {num = min;}
    if (num > max) {num = max;}
    return num;
}

/**
 * Get diffuse
 * 
 * @param norm normal vector
 * @param vec direction vector
 * @return new color value
 */
float getDiffuse(float norm[3], float vec[3]){
    return clamp(dotProduct(normalize(norm), normalize(vec)), 0.2, 1.0);
}

// Materials
//

/**
 * Ray struct
 */
typedef struct {
     float color[3]; //RGB color values
     int reflective; //reflective bool
} material;

// Materials pre-defined
material refl = { .color = {0,1,0}, .reflective = 1 };
material blue = { .color = {0,0,1}, .reflective = 0 };
material red =  { .color = {1,0,0}, .reflective = 0 };
material white= { .color = {1,1,1}, .reflective = 0 };
material green= { .color = {0,1,0}, .reflective = 0 };
material gray= { .color = {0.5,0.5,0.5}, .reflective = 0 };

// Ray struct
//

/**
 * Ray struct
 */
typedef struct {
    float pos[3]; //point of origin
    float dir[3]; //direction
} ray;

/**
 * Creates ray struct
 * 
 * @param ori origin vector
 * @param loc goal vector
 * @return ray object
 */
ray getRay(float ori[3], float loc[3]){
    ray ray;

    //set vector position
    ray.pos[0] = ori[0];
    ray.pos[1] = ori[1];
    ray.pos[2] = ori[2];

    //calculate vector direction
    vectorSub(ray.dir, loc, ray.pos);

    return ray;
}

/**
 * Creates ray struct starting at origin
 * 
 * @param loc goal vector
 * @return ray object
 */
ray getRayO(float loc[3]){
    ray ray;

    //set vector position
    ray.pos[0] = 0.0;
    ray.pos[1] = 0.0;
    ray.pos[2] = 0.0;

    //calculate vector direction
    vectorSub(ray.dir, loc, ray.pos);

    return ray;
}

// Ball struct
//

/**
 * Sphere struct
 */
typedef struct {
    float cen[3]; //center point
    float rad; //radius
    material mat; //material
} ball;

/**
 * Creates ball struct
 * 
 * @param cen vector of center of ball
 * @param rad radius of ball
 * @param mat material of ball
 * @return ball object
 */
ball getBall(float cen[3], float rad, material mat){
    ball sphere;

    //set sphere center vector
    sphere.cen[0] = cen[0];
    sphere.cen[1] = cen[1];
    sphere.cen[2] = cen[2];

    //set sphere radius
    sphere.rad = rad;

    //set sphere material
    sphere.mat = mat;

    return sphere;
}

// Triangle struct
//

/**
 * Triangle struct
 */
typedef struct {
    float v[3][3]; //verticies
    material mat; //material
} triangle;

/**
 * Creates triangle struct
 * 
 * @param v   vectors of triangle points
 * @param mat material of triangle
 * @return ball object
 */
triangle getTriangle(float v[3][3], material mat){
    triangle tri;

    //set triangle points
    for(int i = 0; i < 3; i++){
        tri.v[i][0] = v[i][0];
        tri.v[i][1] = v[i][1];
        tri.v[i][2] = v[i][2];
    }

    //set triangle material
    tri.mat = mat;

    return tri;
}

// Rayhit struct
//

/**
 * Rayhit struct
 */
typedef struct {
    float t; //time ray hits object surface
    material mat; //material of object hit
    float p[3]; //point object hit at
    float norm[3]; //norm of objects surface
} rayhit;

// Helper functions
//

/**
 * Get color of pixel given material and diffuse
 * 
 * @param col color array
 * @param dif diffuse
 * @return new color values
 */
float * getPixColor(float col[3], float dif){
    col[0] = col[0] * dif * 255;
    col[1] = col[1] * dif * 255;
    col[2] = col[2] * dif * 255;
    return col;
}

/**
 * Sets rgb values in image at given index
 * 
 * @param image image array
 * @param color color array
 * @param index index
 */
void setColorAtIndex(unsigned char *image, float color[3], int index){
    image[index] = color[0];
    image[index+1] = color[1];
    image[index+2] = color[2];
}

// Intersections
//

/**
 * Computes all values contained in raystruct for intersect with sphere
 * 
 * @param ray ray struct
 * @param sphere ball struct
 * @return ray hit struct
 */
rayhit sphereIntersect(ray ray, ball sphere){
    float res1, res2;
    float newVec[3]; 
    vectorSub(newVec, ray.pos, sphere.cen);
    rayhit res;
    
    //calculate discriminate
    float disc = pow( dotProduct(ray.dir, newVec), 2) - dotProduct(ray.dir, ray.dir) * ( dotProduct( newVec, newVec ) - pow( sphere.rad, 2 ));
    //if discriminate is negative, return arbitrary negative
    if (disc < 0){
        res.t = -1;
        return res;
    }

    //check both solutions (could be the same)
    res1 = ( -dotProduct( ray.dir, newVec ) + sqrt( disc ) ) / dotProduct( ray.dir, ray.dir );
    res2 = ( -dotProduct( ray.dir, newVec ) - sqrt( disc ) ) / dotProduct( ray.dir, ray.dir );
    
    //if either res is negative, return arbitrary negative
    if (res1 < 0 || res2 < 0) {
        res.t = -1;
        return res;
    }

    //set lowest of the two results
    if (res1 < res2) {
        res.t = res1;
    } else {
        res.t = res2;
    }

    //set material
    res.mat = sphere.mat;

    //set intersction
    rayEquation(res.p, ray.pos, ray.dir, res.t);

    //set norm
    normalize(vectorSub(res.norm, res.p, sphere.cen));

    return res;
}

/**
 * Computes all values contained in raystruct for intersect with triangle
 * 
 * @param ray ray struct
 * @param tri triangle struct
 * @return ray hit struct
 */
rayhit triangleIntersect(ray ray, triangle tri, float closest){
    rayhit res;

    //Triangle example format
    //          x  y  z  
    //          0  1  2
    //   a 0 { -8,-2,-20 }
    //   b 1 {  8,-2,-20 }
    //   c 2 {  8,10,-20 }

    //Set vars needed for math
    float a = tri.v[0][0] - tri.v[1][0];
    float b = tri.v[0][1] - tri.v[1][1]; 
    float c = tri.v[0][2] - tri.v[1][2]; 
    float d = tri.v[0][0] - tri.v[2][0]; 
    float e = tri.v[0][1] - tri.v[2][1]; 
    float f = tri.v[0][2] - tri.v[2][2]; 
    float g = ray.dir[0];
    float h = ray.dir[1];
    float i = ray.dir[2];
    float j = tri.v[0][0] - ray.pos[0];
    float k = tri.v[0][1] - ray.pos[1];
    float l = tri.v[0][2] - ray.pos[2];
    float m = a*(e*i-h*f)+b*(g*f-d*i)+c*(d*h-e*g);

    float t = -(f*(a*k-j*b)+e*(j*c-a*l)+d*(b*l-k*c))/m;
    if (t < 0 || t > closest) {
        res.t = -1;
        return res;
    }

    float gamma = (i*(a*k-j*b)+h*(j*c-a*l)+g*(b*l-k*c))/m;
    if (gamma < 0 || gamma > 1) {
        res.t = -1;
        return res;
    }

    float beta = (j*(e*i-h*f)+k*(g*f-d*i)+l*(d*h-e*g))/m; 
    if (beta < 0 || beta > 1-gamma) {
        res.t = -1;
        return res;
    }

    //set t
    res.t = t;

    //set material
    res.mat = tri.mat;

    //set intersction
    rayEquation(res.p, ray.pos, ray.dir, res.t);

    //set norm
    float v1[3];
    vectorSub(v1, tri.v[1], tri.v[0]);

    float v2[3];
    vectorSub(v2, tri.v[2], tri.v[0]);

    res.norm[0] = (v1[1]*v2[2]) - (v1[2]*v2[1]);
    res.norm[1] = (v1[2]*v2[0]) - (v1[0]*v2[2]);
    res.norm[2] = (v1[0]*v2[1]) - (v1[1]*v2[0]);

    normalize(res.norm);

    return res;
}

// Main
//

/**
 * Main function
 * 
 * @return int
 */
int main(){
    //vars to edit
    int width = 255; //height of image
    int height = width; //width of image
    int outputNum = 60; //number of images to output for zoom gif
    float zoomSpeed = 0.3; //speed to zoom in by
    float light[3] = {3.0, 5.0, -15.0}; //light pos
    float camera[3] = {0.0, 0.0, 0.0}; //camera pos
    float missColor[3] = {0.0, 0.0, 0.0}; //color to print on missed ray

    //misc vars
    int chan = 3;
    int arrLen = width*height*chan;
    unsigned char* image = (unsigned char*) malloc(arrLen*sizeof(int));

    int index;
    float colors[4][3];

    //var for image plane and camera
    float pixWidth = 2.0/width;
    float pixTopLeft[3] = {-1 + (pixWidth/2.0), 1 - (pixWidth/2.0), -2.0};

    //structs
    ray rayOrig;
    ball spheres[100];

    //create balls and walls
    spheres[0] = (ball) { .cen = { 0,0,-16 },   .rad = 2, .mat = refl };
    spheres[1] = (ball) { .cen = { 3,-1,-14 },  .rad = 1, .mat = refl };
    spheres[2] = (ball) { .cen = { -3,1,-14 }, .rad = 1, .mat = red };
    spheres[3] = (ball) { .cen = { -3,-1,-14 }, .rad = 1, .mat = green };
    spheres[4] = (ball) { .cen = { -3,2.5,-14 }, .rad = 0.5, .mat = refl };
    spheres[5] = (ball) { .cen = { 0, -1.5, -11 }, .rad = 0.5, .mat = refl };
    spheres[6] = (ball) { .cen = { 6, 0, -18 }, .rad = 2, .mat = gray };
    int numSpheres = 7;

    triangle triangles[100];
    // back wall
    triangles[0] = (triangle) { .v = { { -8,-2,-20 }, {8,-2,-20}, {8,10,-20} }, .mat = blue };
    triangles[1] = (triangle) { .v = { { -8,-2,-20 }, {8,10,-20}, {-8,10,-20} }, .mat = blue };
    // floor
    triangles[2] = (triangle) { .v = { { -8,-2,-20 }, {8,-2,-10}, {8,-2,-20}}, .mat = white };
    triangles[3] = (triangle) { .v = { { -8,-2,-20 }, {-8,-2,-10}, {8,-2,-10}}, .mat = white };
    // right red triangle
    triangles[4] = (triangle) { .v = { { 8,-2,-20 }, {8,-2,-10}, {8,10,-20}}, .mat = red };
    //left gray triangle
    triangles[5] = (triangle) { .v = { { -8,-2,-20 }, {-8,-2,-10}, {-8,10,-20}}, .mat = gray };
    // overhang
    triangles[6] = (triangle) { .v = { { -8,-2,-10 }, {8,-2,-10}, {8,-4,-10} }, .mat = blue };
    triangles[7] = (triangle) { .v = { { -8,-2,-10 }, {8,-4,-10}, {-8,-4,-10} }, .mat = blue };
    int numTriangles = 8;

    //main loop, loops each pixel in image
    for (int n = 0; n <= outputNum; n++){
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                //index of pixel in image array
                index = (x + (width * y)) * 3;

                for (int pp = 0; pp < 4; pp++) {
                    //Prep Work
                    //
                    int bounces = 0;
                    int shadowCheck = 0;
                    int shadowHit = 0;
                    int shooting = 1;
                    int hit = 0;
                    int distanceToLight = __INT_MAX__;
                    ray rayLight = getRayO(light); //avoid cppcheck error
                    rayhit keeper, check;

                    //location of pixel on image plane
                    float vec[3];
                    vec[0] = pixTopLeft[0] + pixWidth * x;
                    vec[1] = pixTopLeft[1] - pixWidth * y;
                    vec[2] = pixTopLeft[2];
                    switch (pp)
                    {
                    case 0:
                        vec[0] -= pixWidth / 4;
                        vec[1] += pixWidth / 4;
                        break;
                    case 1:
                        vec[0] += pixWidth / 4;
                        vec[1] += pixWidth / 4;
                        break;
                    case 2:
                        vec[0] -= pixWidth / 4;
                        vec[1] -= pixWidth / 4;
                        break;
                    case 3:
                        vec[0] += pixWidth / 4;
                        vec[1] -= pixWidth / 4;
                        break;
                    
                    default:
                        break;
                    }

                    //get ray struct for pixel
                    rayOrig = getRay(camera, normalize(vec));

                    //Shoot ray and find surface hit if any
                    //
                    while (shooting) {
                        int closest = __INT_MAX__;
                        
                        //loop through spheres
                        for(int i = 0; i < numSpheres; i++){
                            if (shadowCheck) {
                                //determine if hit shadow
                                check = sphereIntersect(rayLight, spheres[i]);
                                if (check.t > 0 && check.t < distanceToLight) {
                                    shadowHit = 1;
                                }
                            } else {
                                //determine first sphere hit
                                check = sphereIntersect(rayOrig, spheres[i]);
                                if (check.t > 0 && check.t < closest) {
                                    hit = 1;
                                    keeper = check;
                                    closest = check.t;
                                }
                            }
                        }

                        //loop through triangles
                        for(int i = 0; i < numTriangles; i++){
                            if (shadowCheck) {
                                //determine if hit shadow
                                check = triangleIntersect(rayLight, triangles[i], __INT_MAX__);
                                if (check.t > 0 && check.t < distanceToLight) {
                                    shadowHit = 1;
                                }
                            } else {
                                //determine first triangle hit
                                check = triangleIntersect(rayOrig, triangles[i], closest);
                                if (check.t > 0 && check.t < closest) {
                                    hit = 1;
                                    keeper = check;
                                    closest = check.t;
                                }
                            }
                        }

                        //if no hit or complete shadow check, exit
                        if (shadowCheck || !hit) {
                            shooting = 0;
                        }  
                        //if hit non reflective surface, check for shadows
                        if ((hit == 1 && keeper.mat.reflective == 0) || bounces >= 10) {
                            shadowCheck = 1;
                            rayLight = getRay(keeper.p, light);
                            fidget(rayLight.pos, rayLight.dir, 0.0001);
                            distanceToLight = pointDistance3D(rayLight.pos, light);
                        }
                        //otherwise, keep shooting 
                        else 
                        {
                            hit = 0;
                            bounces += 1;
                            vectorCpy(rayOrig.pos, keeper.p);
                            scale(keeper.norm, -2*(dotProduct(rayOrig.dir, keeper.norm)));
                            vectorAdd(rayOrig.dir, rayOrig.dir, keeper.norm);
                            //fidget(rayOrig.pos, rayOrig.dir, 0.0001);
                        }
                    }

                    //color pixel based on hit
                    //
                    if (hit) {
                        float diffuse;
                        if (!shadowHit) {
                            diffuse = getDiffuse(keeper.norm, rayLight.dir);
                        } else {
                            diffuse = 0.2;
                        }
                        vectorCpy(colors[pp], getPixColor(keeper.mat.color, diffuse));
                    } else {
                        vectorCpy(colors[pp], missColor);
                    }
                }

                //color pixel based colors gathered
                //
                float average[3];
                average[0] = (colors[0][0] + colors[1][0] + colors[2][0] + colors[3][0]) / 4;
                average[1] = (colors[0][1] + colors[1][1] + colors[2][1] + colors[3][1]) / 4;
                average[2] = (colors[0][2] + colors[1][2] + colors[2][2] + colors[3][2]) / 4;
                setColorAtIndex(image, average, index);
            }
        }

        //write image to file 
        char fileName[100];
        sprintf(fileName, "%03d.png", n);
        stbi_write_png(fileName, width, height, chan, image, width*chan);

        //move camera
        if (n < outputNum*0.667) {
            camera[2] = camera[2]+zoomSpeed;
        } else {
            camera[2] = camera[2]-(zoomSpeed*2);
        }
    }
    
    //free image memory
    free(image);

    return 0;
}