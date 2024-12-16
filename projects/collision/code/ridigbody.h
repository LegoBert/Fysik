#pragma once
#include "config.h"
#include "render/model.h"
#include "render/debugrender.h"

using namespace glm;
using namespace std;

class Collider {
public:
    std::vector<glm::vec3> localVertices;       // Original collider vertices in local space
    std::vector<glm::vec3> worldVertices;       // Transformed collider vertices in world space

    Collider(Render::ModelId modelId, glm::mat4 transform) {
        for (const auto& mesh : Render::GetModel(modelId).meshes) {
            for (const auto& primitive : mesh.primitives) {
                for (const auto& vertex : primitive.vertices) {
                    localVertices.push_back(vertex);
                    worldVertices.push_back(glm::vec3(transform * glm::vec4(vertex, 1.0f)));
                }
            }
        }
    }

    void UpdateCollider(glm::mat4 transform) {
        worldVertices.clear();
        for (const auto& vertex : localVertices) {
            worldVertices.push_back(glm::vec3(transform * glm::vec4(vertex, 1.0f)));
        }
    }

    bool DetectCollision(Collider& other) {
        return GJK(worldVertices, other.worldVertices);
    }

    // GJK
    glm::vec3 FindFurthestPoint(std::vector<glm::vec3>& vertices, glm::vec3 direction) {
        float maxDistance = -FLT_MAX;
        glm::vec3 furthestPoint(0);

        for (const auto& vertex : vertices) {
            float distance = glm::dot(vertex, direction);
            if (distance > maxDistance) {
                maxDistance = distance;
                furthestPoint = vertex;
            }
        }

        return furthestPoint;
    }

    glm::vec3 Support(std::vector<glm::vec3>& shapeA, std::vector<glm::vec3>& shapeB, glm::vec3& direction) {
        return FindFurthestPoint(shapeA, direction) - FindFurthestPoint(shapeB, -direction);
    }

    bool DoSimplex(std::vector<glm::vec3>& simplex, glm::vec3& direction) {
        const float epsilon = 1e-6f; // tolerance to prevent infinite loops due to floating point inaccuracies

        // Line case
        if (simplex.size() == 2) {
            glm::vec3 A = simplex[1];   // Most recent point
            glm::vec3 B = simplex[0];   // Previous point

            glm::vec3 AB = B - A;
            glm::vec3 AO = -A;

            // Check if origin is on the AB line segment
            if (glm::dot(AB, AO) > 0) {
                direction = glm::cross(glm::cross(AB, AO), AB);
            }
            else {
                direction = AO;
            }

            return false;
        }

        // Triangle case
        if (simplex.size() == 3) {
            glm::vec3 A = simplex[2];   // Most recent point
            glm::vec3 B = simplex[1];
            glm::vec3 C = simplex[0];

            glm::vec3 AB = B - A;
            glm::vec3 AC = C - A;
            glm::vec3 AO = -A;

            glm::vec3 ABC = glm::cross(AB, AC);

            // Check which region the origin lies in
            if (glm::dot(glm::cross(ABC, AC), AO) > 0) {
                if (glm::dot(AC, AO) > 0) {
                    simplex = { C, A };
                    direction = glm::cross(glm::cross(AC, AO), AC);
                }
                else {
                    if (glm::dot(AB, AO) > 0) {
                        simplex = { B, A };
                        direction = glm::cross(glm::cross(AB, AO), AB);
                    }
                    else {
                        simplex = { A };
                        direction = AO;
                    }
                }
            }
            else {
                if (glm::dot(glm::cross(AB, ABC), AO) > 0) {
                    if (glm::dot(AB, AO) > 0) {
                        simplex = { B, A };
                        direction = glm::cross(glm::cross(AB, AO), AB);
                    }
                    else {
                        simplex = { A };
                        direction = AO;
                    }
                }
                else {
                    if (glm::dot(ABC, AO) > 0) {
                        direction = ABC;
                    }
                    else {
                        direction = -ABC;
                    }
                }
            }

            return false;
        }

        // Tetrahedron case
        if (simplex.size() == 4) {
            glm::vec3 A = simplex[3];   // Most recent point
            glm::vec3 B = simplex[2];
            glm::vec3 C = simplex[1];
            glm::vec3 D = simplex[0];

            glm::vec3 AB = B - A;
            glm::vec3 AC = C - A;
            glm::vec3 AD = D - A;
            glm::vec3 AO = -A;

            glm::vec3 ABC = glm::cross(AB, AC);
            glm::vec3 ACD = glm::cross(AC, AD);
            glm::vec3 ADB = glm::cross(AD, AB);

            // Check which face the origin is on
            if (glm::dot(ABC, AO) > 0) {
                simplex = { C, B, A };
                direction = ABC;
            }
            else if (glm::dot(ACD, AO) > 0) {
                simplex = { D, C, A };
                direction = ACD;
            }
            else if (glm::dot(ADB, AO) > 0) {
                simplex = { B, D, A };
                direction = ADB;
            }
            else {
                return true;  // The origin is inside the tetrahedron
            }

            if (glm::length(direction) < epsilon) return true;  // Collision if direction is too small
            return false;
        }

        return false;
    }

    bool GJK(std::vector<glm::vec3>& shapeA, std::vector<glm::vec3>& shapeB) {
        auto direction = glm::vec3(1, 0, 0); // random
        std::vector<glm::vec3> simplex = { Support(shapeA, shapeB, direction) };
        direction = -simplex[0];

        while (true) {
            // Get a new point in the direction of the origin
            glm::vec3 A = Support(shapeA, shapeB, direction);

            // If the new point isn't past the origin in the given direction, there's no collision
            if (glm::dot(A, direction) <= 0) {
                return false;   // No collision
            }

            // Add the new point to the simplex
            simplex.push_back(A);

            // Check if the simplex contains the origin
            if (DoSimplex(simplex, direction)) {
                Debug::DrawLine(simplex[0], simplex[1], 1.0f, glm::vec4(1, 1, 1, 1), glm::vec4(1, 1, 1, 1), Debug::RenderMode::Normal);
                Debug::DrawLine(simplex[1], simplex[2], 1.0f, glm::vec4(1, 1, 1, 1), glm::vec4(1, 1, 1, 1), Debug::RenderMode::Normal);
                Debug::DrawLine(simplex[2], simplex[0], 1.0f, glm::vec4(1, 1, 1, 1), glm::vec4(1, 1, 1, 1), Debug::RenderMode::Normal);
                Debug::DrawLine(simplex[0], simplex[3], 1.0f, glm::vec4(1, 1, 1, 1), glm::vec4(1, 1, 1, 1), Debug::RenderMode::Normal);
                Debug::DrawLine(simplex[1], simplex[3], 1.0f, glm::vec4(1, 1, 1, 1), glm::vec4(1, 1, 1, 1), Debug::RenderMode::Normal);
                Debug::DrawLine(simplex[2], simplex[3], 1.0f, glm::vec4(1, 1, 1, 1), glm::vec4(1, 1, 1, 1), Debug::RenderMode::Normal);

                vec3 epaResult = EPA(simplex, shapeA, shapeB).Normal;
                printf("Collision Normal: (%f, %f, %f)\n", epaResult.x, epaResult.y, epaResult.z);
                Debug::DrawLine(vec3(0), epaResult, 1.0f, glm::vec4(1, 0.5, 1, 1), glm::vec4(1, 0.5, 1, 1), Debug::RenderMode::Normal);

                return true;    // Collision detected
            }
        }
    }

    // EPA
    std::pair<std::vector<vec4>, size_t> GetFaceNormals(const std::vector<vec3>& polytope, const std::vector<size_t>& faces)
    {
        std::vector<vec4> normals;
        size_t minTriangle = 0;
        float  minDistance = FLT_MAX;

        for (size_t i = 0; i < faces.size(); i += 3) {
            vec3 a = polytope[faces[i]];
            vec3 b = polytope[faces[i + 1]];
            vec3 c = polytope[faces[i + 2]];

            vec3 normal = glm::normalize(cross(b - a, c - a));
            float distance = dot(normal, a);

            if (distance < 0) {
                normal *= -1;
                distance *= -1;
            }

            normals.emplace_back(normal, distance);

            if (distance < minDistance) {
                minTriangle = i / 3;
                minDistance = distance;
            }
        }

        return { normals, minTriangle };
    }

    void AddIfUniqueEdge(std::vector<std::pair<size_t, size_t>>& edges, const std::vector<size_t>& faces, size_t a, size_t b)
    {
        auto reverse = std::find(                       //      0--<--3
            edges.begin(),                              //     / \ B /   A: 2-0
            edges.end(),                                //    / A \ /    B: 0-2
            std::make_pair(faces[b], faces[a]) //   1-->--2
        );

        if (reverse != edges.end()) {
            edges.erase(reverse);
        }

        else {
            edges.emplace_back(faces[a], faces[b]);
        }
    }

    struct CollisionPoints {
        glm::vec3 Normal;
        double PenetrationDepth;
        bool HasCollision;
    };

    CollisionPoints EPA(std::vector<glm::vec3> simplex, std::vector<glm::vec3>& colliderA, std::vector<glm::vec3>& colliderB) {
        std::vector<vec3> polytope(simplex.begin(), simplex.end());
        std::vector<size_t> faces = {
            0, 1, 2,
            0, 3, 1,
            0, 2, 3,
            1, 3, 2
        };

        // list: vec4(normal, distance), index: min distance
        auto [normals, minFace] = GetFaceNormals(polytope, faces);

        vec3  minNormal;
        float minDistance = FLT_MAX;

        while (minDistance == FLT_MAX) {
            minNormal = normals[minFace];
            minDistance = normals[minFace].w;

            vec3 support = Support(colliderA, colliderB, minNormal);
            float sDistance = dot(minNormal, support);

            if (abs(sDistance - minDistance) > 0.001f) {
                minDistance = FLT_MAX;
                std::vector<std::pair<size_t, size_t>> uniqueEdges;

                for (size_t i = 0; i < normals.size(); i++) {
                    if (glm::dot(vec3(normals[i]), support - polytope[faces[i * 3]]) > 0) {
                        size_t f = i * 3;

                        AddIfUniqueEdge(uniqueEdges, faces, f, f + 1);
                        AddIfUniqueEdge(uniqueEdges, faces, f + 1, f + 2);
                        AddIfUniqueEdge(uniqueEdges, faces, f + 2, f);

                        faces[f + 2] = faces.back(); faces.pop_back();
                        faces[f + 1] = faces.back(); faces.pop_back();
                        faces[f] = faces.back(); faces.pop_back();

                        normals[i] = normals.back(); // pop-erase
                        normals.pop_back();

                        i--;
                    }
                }
                std::vector<size_t> newFaces;
                for (auto [edgeIndex1, edgeIndex2] : uniqueEdges) {
                    newFaces.push_back(edgeIndex1);
                    newFaces.push_back(edgeIndex2);
                    newFaces.push_back(polytope.size());
                }

                polytope.push_back(support);

                auto [newNormals, newMinFace] = GetFaceNormals(polytope, newFaces);
                float oldMinDistance = FLT_MAX;
                for (size_t i = 0; i < normals.size(); i++) {
                    if (normals[i].w < oldMinDistance) {
                        oldMinDistance = normals[i].w;
                        minFace = i;
                    }
                }

                if (newNormals[newMinFace].w < oldMinDistance) {
                    minFace = newMinFace + normals.size();
                }

                faces.insert(faces.end(), newFaces.begin(), newFaces.end());
                normals.insert(normals.end(), newNormals.begin(), newNormals.end());
            }
        }

        CollisionPoints points;
        points.Normal = minNormal;
        points.PenetrationDepth = minDistance + 0.001f;
        points.HasCollision = true;

        return points;
    }

};

class RidigBody
{
public:
    // Physical state variables
    glm::vec3 position;          // Center of mass position
    glm::vec3 linearVelocity;    // Linear velocity
    glm::vec3 linearAcceleration;// Linear acceleration
    glm::quat orientation;       // Rotation represented as quaternion
    glm::vec3 angularVelocity;   // Angular velocity (rotational)
    glm::vec3 angularAcceleration; // Angular acceleration

    // Net force and torque
    glm::vec3 netForce;          // Accumulated force
    glm::vec3 netTorque;         // Accumulated torque

    // Mass and inertia properties
    float mass;
    float momentOfInertia;

    // Coefficient of restitution for collisions
    float restitutionCoefficient;

    glm::mat4 transform;
    
public:
    RidigBody(float mass, float restitutionCoefficient, glm::vec3 pos);
    void ApplyForceAtCenter(glm::vec3 force);
    void ApplyForceAtPoint(const glm::vec3& force, const glm::vec3& pointOfContact);
    void Update(float deltaTime);
    //bool GJK(const std::vector<glm::vec3>& shapeA, const std::vector<glm::vec3>& shapeB);
};