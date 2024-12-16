#pragma once
#include "config.h"
#include "render/debugrender.h"
#include <iostream>
#include "render/model.h"

namespace Game
{

struct Plane
{
    glm::vec3 position = glm::vec3(0);
    glm::vec3 normal = glm::vec3(0.0f, 1.0f, 0.0f);
    Plane() = default;
    Plane(glm::vec3 position, glm::vec3 normal) : position(position), normal(glm::normalize(normal)) {}

    void DrawDebugPlane(float size) {
        // Check if the normal is almost equal to up
        glm::vec3 up(0, 1, 0);
        if (glm::dot(normal, up) > 0.999f) {
            up = glm::vec3(1, 0, 0);
        }

        // Compute right and forward directions relative to the plane's normal
        glm::vec3 right = glm::normalize(glm::cross(normal, up));
        glm::vec3 forward = glm::normalize(glm::cross(right, normal));

        //// Four corners of the plane rectangle
        glm::vec3 topLeft = position - right * size / 2.0f - forward * size / 2.0f;
        glm::vec3 topRight = position + right * size / 2.0f - forward * size / 2.0f;
        glm::vec3 bottomLeft = position - right * size / 2.0f + forward * size / 2.0f;
        glm::vec3 bottomRight = position + right * size / 2.0f + forward * size / 2.0f;

        //// Draw lines for the plane rectangle
        Debug::DrawLine(topLeft, topRight, 1.0f, glm::vec4(0, 1, 0, 1), glm::vec4(0, 1, 0, 1), Debug::RenderMode::Normal);
        Debug::DrawLine(topRight, bottomRight, 1.0f, glm::vec4(0, 1, 0, 1), glm::vec4(0, 1, 0, 1), Debug::RenderMode::Normal);
        Debug::DrawLine(bottomRight, bottomLeft, 1.0f, glm::vec4(0, 1, 0, 1), glm::vec4(0, 1, 0, 1), Debug::RenderMode::Normal);
        Debug::DrawLine(bottomLeft, topLeft, 1.0f, glm::vec4(0, 1, 0, 1), glm::vec4(0, 1, 0, 1), Debug::RenderMode::Normal);
        Debug::DrawLine(position - right * size / 2.0f, position + right * size / 2.0f, 1.0f, glm::vec4(0, 1, 0, 1), glm::vec4(0, 1, 0, 1), Debug::RenderMode::Normal);
        Debug::DrawLine(position - forward * size / 2.0f, position + forward * size / 2.0f, 1.0f, glm::vec4(0, 1, 0, 1), glm::vec4(0, 1, 0, 1), Debug::RenderMode::Normal);
        Debug::DrawLine(position, position + normal * 1.0f, 1.0f, glm::vec4(0, 0, 1, 1), glm::vec4(0, 0, 1, 1), Debug::RenderMode::Normal);
    }
};

struct Ray
{
    glm::vec3 origin;
    glm::vec3 direction;

    struct RaycastHit
    {
        bool hit;
        glm::vec3 point;
        glm::vec3 normal;
        float distance;
        RaycastHit() : hit(false), point(glm::vec3(0)), normal(glm::vec3(0)), distance(0.0f) {}
        RaycastHit(const glm::vec3& point, const glm::vec3& normal, float distance) : hit(true), point(point), normal(normal), distance(distance) {}
    };

    struct IntersectionTriangle {
        glm::vec3 a, b, c;
        IntersectionTriangle() = default;
        IntersectionTriangle(const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& v3) : a(v1), b(v2), c(v3) {}
        void DrawDebugTriangle() {
            Debug::DrawLine(a, b, 1.0f, glm::vec4(0, 1, 0, 1), glm::vec4(0, 1, 0, 1), Debug::RenderMode::Normal);
            Debug::DrawLine(b, c, 1.0f, glm::vec4(0, 1, 0, 1), glm::vec4(0, 1, 0, 1), Debug::RenderMode::Normal);
            Debug::DrawLine(c, a, 1.0f, glm::vec4(0, 1, 0, 1), glm::vec4(0, 1, 0, 1), Debug::RenderMode::Normal);
        }
    };

    Ray(glm::vec3 origin, glm::vec3 direction) : origin(origin), direction(glm::normalize(direction)) {}

    void DrawDebugRay(RaycastHit hitresult,float lenght) {
        //Draw the ray
        Debug::DrawLine(origin, origin + direction * lenght, 1.0f, glm::vec4(1, 0, 0, 1), glm::vec4(1, 0, 0, 1), Debug::RenderMode::Normal);
        //Draw hit point or not
        if (hitresult.hit) {
            glm::vec3 scaleFactors = glm::vec3(.05f, .05f, .05f);
            glm::mat4 scaleMatrix = glm::scale(glm::mat4(1.0f), scaleFactors);
            glm::mat4 translationMatrix = glm::translate(glm::mat4(1.0f), hitresult.point);
            glm::mat4 combinedMatrix = translationMatrix * scaleMatrix;
            Debug::DrawBox(combinedMatrix, glm::vec4(1, 0, 0, 1), Debug::Normal);
        }
    }

    RaycastHit IntersectPlane(Plane& p) {
        // Calculate the denominator of the intersection formula (N . D)
        float denom = glm::dot(p.normal, direction);

        // If the denominator is close to zero, the ray is parallel to the plane
        if (abs(denom) > 1e-6) { // Use an epsilon value to prevent division by zero
            // Calculate the nominator part of the t formula (N . (Q - O))
            float t = glm::dot(p.normal, p.position - origin) / denom;

            // If t is negative, the intersection point is behind the ray origin
            if (t >= 0) {
                // Calculate the actual intersection point
                glm::vec3 point = origin + direction * t;

                // Return a RaycastHit with the intersection information
                return RaycastHit(point, p.normal, t);
            }
        }

        // If the ray does not intersect or intersects behind the origin,
        // return a RaycastHit with hit set to false.
        return RaycastHit();
    }

    RaycastHit IntersectTriangle(IntersectionTriangle& tri) { 
        // Define the plane in which the triangle lies
        Plane triPlane;
        triPlane.position = tri.a; // Use vertex a as a point on the plane
        glm::vec3 edge1 = tri.b - tri.a;
        glm::vec3 edge2 = tri.c - tri.a;
        triPlane.normal = glm::normalize(glm::cross(edge1, edge2)); // Plane's normal is perpendicular to the triangle

        // Use IntersectPlane to find the intersection point with the plane
        RaycastHit planeHit = IntersectPlane(triPlane);

        // If there's no hit, return no hit
        if (!planeHit.hit) return RaycastHit();

        // Perform the inside-outside test
        glm::vec3 pt = planeHit.point; // Intersection point
        glm::vec3 C; // Vector perpendicular to triangle's plane

        // Edge 0
        glm::vec3 edge0 = tri.b - tri.a;
        glm::vec3 vp0 = pt - tri.a;
        C = glm::cross(edge0, vp0);
        if (glm::dot(triPlane.normal, C) < 0) return RaycastHit(); // P is on the right side

        // Edge 1
        glm::vec3 edge1b = tri.c - tri.b;
        glm::vec3 vp1 = pt - tri.b;
        C = glm::cross(edge1b, vp1);
        if (glm::dot(triPlane.normal, C) < 0)  return RaycastHit(); // P is on the right side

        // Edge 2
        glm::vec3 edge2b = tri.a - tri.c;
        glm::vec3 vp2 = pt - tri.c;
        C = glm::cross(edge2b, vp2);
        if (glm::dot(triPlane.normal, C) < 0) return RaycastHit(); // P is on the right side

        // If the point passes all edge tests, it lies within the triangle
        return planeHit;
    }

    bool IntersectRayAABB(Render::Model::AABB aabb, glm::mat4 modelTransform) {
        // Extract translation component from the model transform
        glm::vec3 modelPosition = glm::vec3(modelTransform[3]);

        // Transform the AABB using the model transform matrix
        /*glm::vec3 aabbMin = glm::vec3(modelTransform * glm::vec4(aabb.min, 1.0f));
        glm::vec3 aabbMax = glm::vec3(modelTransform * glm::vec4(aabb.max, 1.0f));*/
        glm::mat3 rotationMatrix = glm::mat3(modelTransform);
        glm::vec3 aabbExtents = aabb.max - aabb.min;
        glm::vec3 rotatedExtents = glm::abs(rotationMatrix[0]) * aabbExtents.x +
            glm::abs(rotationMatrix[1]) * aabbExtents.y +
            glm::abs(rotationMatrix[2]) * aabbExtents.z;

        glm::vec3 aabbMin = modelPosition - rotatedExtents * 0.5f;
        glm::vec3 aabbMax = modelPosition + rotatedExtents * 0.5f;

        glm::vec3 invDir = 1.0f / direction;

        glm::vec3 t1 = (aabbMin - origin) * invDir;
        glm::vec3 t2 = (aabbMax - origin) * invDir;

        glm::vec3 tmin = glm::min(t1, t2);
        glm::vec3 tmax = glm::max(t1, t2);

        float tMin = glm::max(tmin.x, glm::max(tmin.y, tmin.z));
        float tMax = glm::min(tmax.x, glm::min(tmax.y, tmax.z));

        // Draw the AABB using Debug::DrawLine
        Debug::DrawLine(aabbMin, glm::vec3(aabbMax.x, aabbMin.y, aabbMin.z), 1.0f, glm::vec4(1, 1, 1, 1), glm::vec4(1, 1, 1, 1), Debug::RenderMode::Normal);
        Debug::DrawLine(glm::vec3(aabbMax.x, aabbMin.y, aabbMin.z), glm::vec3(aabbMax.x, aabbMax.y, aabbMin.z), 1.0f, glm::vec4(1, 1, 1, 1), glm::vec4(1, 1, 1, 1), Debug::RenderMode::Normal);
        Debug::DrawLine(glm::vec3(aabbMax.x, aabbMax.y, aabbMin.z), glm::vec3(aabbMin.x, aabbMax.y, aabbMin.z), 1.0f, glm::vec4(1, 1, 1, 1), glm::vec4(1, 1, 1, 1), Debug::RenderMode::Normal);
        Debug::DrawLine(glm::vec3(aabbMin.x, aabbMax.y, aabbMin.z), aabbMin, 1.0f, glm::vec4(1, 1, 1, 1), glm::vec4(1, 1, 1, 1), Debug::RenderMode::Normal);

        Debug::DrawLine(glm::vec3(aabbMin.x, aabbMin.y, aabbMax.z), glm::vec3(aabbMax.x, aabbMin.y, aabbMax.z), 1.0f, glm::vec4(1, 1, 1, 1), glm::vec4(1, 1, 1, 1), Debug::RenderMode::Normal);
        Debug::DrawLine(glm::vec3(aabbMax.x, aabbMin.y, aabbMax.z), glm::vec3(aabbMax.x, aabbMax.y, aabbMax.z), 1.0f, glm::vec4(1, 1, 1, 1), glm::vec4(1, 1, 1, 1), Debug::RenderMode::Normal);
        Debug::DrawLine(glm::vec3(aabbMax.x, aabbMax.y, aabbMax.z), glm::vec3(aabbMin.x, aabbMax.y, aabbMax.z), 1.0f, glm::vec4(1, 1, 1, 1), glm::vec4(1, 1, 1, 1), Debug::RenderMode::Normal);
        Debug::DrawLine(glm::vec3(aabbMin.x, aabbMax.y, aabbMax.z), glm::vec3(aabbMin.x, aabbMin.y, aabbMax.z), 1.0f, glm::vec4(1, 1, 1, 1), glm::vec4(1, 1, 1, 1), Debug::RenderMode::Normal);

        Debug::DrawLine(aabbMin, glm::vec3(aabbMin.x, aabbMin.y, aabbMax.z), 1.0f, glm::vec4(1, 1, 1, 1), glm::vec4(1, 1, 1, 1), Debug::RenderMode::Normal);
        Debug::DrawLine(glm::vec3(aabbMax.x, aabbMin.y, aabbMin.z), glm::vec3(aabbMax.x, aabbMin.y, aabbMax.z), 1.0f, glm::vec4(1, 1, 1, 1), glm::vec4(1, 1, 1, 1), Debug::RenderMode::Normal);
        Debug::DrawLine(glm::vec3(aabbMax.x, aabbMax.y, aabbMin.z), glm::vec3(aabbMax.x, aabbMax.y, aabbMax.z), 1.0f, glm::vec4(1, 1, 1, 1), glm::vec4(1, 1, 1, 1), Debug::RenderMode::Normal);
        Debug::DrawLine(glm::vec3(aabbMin.x, aabbMax.y, aabbMax.z), glm::vec3(aabbMin.x, aabbMax.y, aabbMin.z), 1.0f, glm::vec4(1, 1, 1, 1), glm::vec4(1, 1, 1, 1), Debug::RenderMode::Normal);

        return tMax >= tMin;
    }

    RaycastHit IntersectModel(const Render::Model& model, glm::mat4 modelTransform) {
        RaycastHit closestHit;
        closestHit.distance = std::numeric_limits<float>::max(); // Set to maximum possible distance

        /*
        // Compute or retrieve the AABB of the model
        AABB modelAABB = ComputeOrRetrieveModelAABB(model);

        // Transform AABB with the same model transformation
        modelAABB.Transform(modelTransform);
        */

        // Check for intersection with the AABB first
        if (!IntersectRayAABB(model.aabb, modelTransform)) {
            return RaycastHit(); // No intersection with the model's AABB, return no hit
        }

        for (const auto& mesh : model.meshes) {
            for (const auto& primitive : mesh.primitives) {
                for (size_t i = 0; i < primitive.numIndices / 3; ++i) {
                    IntersectionTriangle tri;
                    tri.a = glm::vec3(modelTransform * glm::vec4(primitive.vertices[primitive.indices[3 * i + 0]], 1.0f));
                    tri.b = glm::vec3(modelTransform * glm::vec4(primitive.vertices[primitive.indices[3 * i + 1]], 1.0f));
                    tri.c = glm::vec3(modelTransform * glm::vec4(primitive.vertices[primitive.indices[3 * i + 2]], 1.0f));
                    RaycastHit hit = IntersectTriangle(tri);
                    if (hit.hit && hit.distance < closestHit.distance) {
                        closestHit = hit; // Update closest hit if this hit is closer
                    }
                }
                glBindVertexArray(0);
            }
        }

        if (closestHit.distance < std::numeric_limits<float>::max()) {
            return closestHit; // Return the closest intersection
        }
        else {
            return RaycastHit(); // Return no hit if no intersections were found
        }
    }

};

}