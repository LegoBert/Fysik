#include "ridigbody.h"

RidigBody::RidigBody(float mass, float restitutionCoefficient, glm::vec3 pos) {
    this->position = pos;
    this->linearVelocity = glm::vec3(0.0f);
    this->linearAcceleration = glm::vec3(0.0f);
    this->orientation = glm::quat(1.0f, 0.0f, 0.0f, 0.0f);
    this->angularVelocity = glm::vec3(0.0f);
    this->angularAcceleration = glm::vec3(0.0f);

    this->netForce = glm::vec3(0.0f);
    this->netTorque = glm::vec3(0.0f);

    this->mass = mass;
    this->momentOfInertia = 1.0f;
    this->restitutionCoefficient = restitutionCoefficient;

    this->transform = glm::translate(glm::mat4(1.0f), position) * glm::mat4_cast(orientation);
}

// Method to apply a force at the center of mass
void RidigBody::ApplyForceAtCenter(glm::vec3 force)
{
    netForce += force;
}

// Method to apply a force at a specific point (generating torque if not at the center)
void RidigBody::ApplyForceAtPoint(const glm::vec3& force, const glm::vec3& pointOfContact) {
    // Force always contributes to linear motion
    netForce += force;

    // If force is applied away from the center of mass, it generates torque
    glm::vec3 r = pointOfContact - position;  // Vector from center of mass to the point
    glm::vec3 torque = glm::cross(r, force);  // Torque is r x F (cross product)
    netTorque += torque;  // Accumulate the torque
}

void RidigBody::Update(float deltaTime) {
    // Update linear motion
    linearAcceleration = netForce / mass;  // a = F/m
    //linearAcceleration += glm::vec3(0, -9.807f, 0); // gravity
    linearVelocity += linearAcceleration * deltaTime;
    position += linearVelocity * deltaTime;
    
    // Update angular motion
    angularAcceleration = netTorque / momentOfInertia; // Need inertia tensor!
    angularVelocity += angularAcceleration * deltaTime;
    //printf("angularVelocity: %f, %f, %f \n", angularVelocity.x, angularVelocity.y, angularVelocity.z);

    // Update orientation using quaternion integration
    glm::quat angularVelocityQuat(0, angularVelocity);
    orientation += 0.5f * angularVelocityQuat * orientation * deltaTime; // Update orientation
    orientation = glm::normalize(orientation);  // Normalize to avoid drift
    
    // Debug draw angular velocity
    Debug::DrawLine(position, position + glm::normalize(angularVelocity) * 5.0f, 1.0f, glm::vec4(1, 1, 0, 1), glm::vec4(1, 1, 0, 1), Debug::RenderMode::Normal);
    Debug::DrawLine(position, position - glm::normalize(angularVelocity) * 5.0f, 1.0f, glm::vec4(1, 1, 0, 1), glm::vec4(1, 1, 0, 1), Debug::RenderMode::Normal);

    // Update transform
    transform = glm::translate(glm::mat4(1.0f), position) * glm::mat4_cast(orientation);

    // Reset net forces and torques for the next frame
    netForce = glm::vec3(0.0f);
    netTorque = glm::vec3(0.0f);
}

// Find the furthest point in a certain direction
//glm::vec3 Support(const std::vector<glm::vec3>& shape, const glm::vec3& direction) {
//    float maxDot = -std::numeric_limits<float>::infinity();
//    glm::vec3 furthestPoint;
//
//    for (const auto& point : shape) {
//        float dotProduct = glm::dot(point, direction);
//        if (dotProduct > maxDot) {
//            maxDot = dotProduct;
//            furthestPoint = point;
//        }
//    }
//
//    return furthestPoint;
//}
//
//glm::vec3 MinkowskiDiffrence(const std::vector<glm::vec3>& shapeA, const std::vector<glm::vec3>& shapeB, const glm::vec3& direction) {
//    glm::vec3 pointA = Support(shapeA, direction);
//    glm::vec3 pointB = Support(shapeB, -direction);
//    return pointA - pointB;
//}
//
//bool DoSimplex(std::vector<glm::vec3>& simplex, glm::vec3& direction) {
//    const float epsilon = 1e-6f; // tolerance to prevent infinite loops due to floating point inaccuracies
//
//    // Line case
//    if (simplex.size() == 2) {
//        glm::vec3 A = simplex[1];   // Most recent point
//        glm::vec3 B = simplex[0];   // Previous point
//
//        glm::vec3 AB = B - A;
//        glm::vec3 AO = -A;
//
//        // Check if origin is on the AB line segment
//        if (glm::dot(AB, AO) > 0) {
//            direction = glm::cross(glm::cross(AB, AO), AB);
//        }
//        else {
//            direction = AO;
//        }
//
//        return false;
//    }
//
//    // Triangle case
//    if (simplex.size() == 3) {
//        glm::vec3 A = simplex[2];   // Most recent point
//        glm::vec3 B = simplex[1];
//        glm::vec3 C = simplex[0];
//
//        glm::vec3 AB = B - A;
//        glm::vec3 AC = C - A;
//        glm::vec3 AO = -A;
//
//        glm::vec3 ABC = glm::cross(AB, AC);
//
//        // Check which region the origin lies in
//        if (glm::dot(glm::cross(ABC, AC), AO) > 0) {
//            if (glm::dot(AC, AO) > 0) {
//                simplex = { C, A };
//                direction = glm::cross(glm::cross(AC, AO), AC);
//            }
//            else {
//                if (glm::dot(AB, AO) > 0) {
//                    simplex = { B, A };
//                    direction = glm::cross(glm::cross(AB, AO), AB);
//                }
//                else {
//                    simplex = { A };
//                    direction = AO;
//                }
//            }
//        }
//        else {
//            if (glm::dot(glm::cross(AB, ABC), AO) > 0) {
//                if (glm::dot(AB, AO) > 0) {
//                    simplex = { B, A };
//                    direction = glm::cross(glm::cross(AB, AO), AB);
//                }
//                else {
//                    simplex = { A };
//                    direction = AO;
//                }
//            }
//            else {
//                if (glm::dot(ABC, AO) > 0) {
//                    direction = ABC;
//                }
//                else {
//                    direction = -ABC;
//                }
//            }
//        }
//
//        return false;
//    }
//
//    // Tetrahedron case
//    if (simplex.size() == 4) {
//        glm::vec3 A = simplex[3];   // Most recent point
//        glm::vec3 B = simplex[2];
//        glm::vec3 C = simplex[1];
//        glm::vec3 D = simplex[0];
//
//        glm::vec3 AB = B - A;
//        glm::vec3 AC = C - A;
//        glm::vec3 AD = D - A;
//        glm::vec3 AO = -A;
//
//        glm::vec3 ABC = glm::cross(AB, AC);
//        glm::vec3 ACD = glm::cross(AC, AD);
//        glm::vec3 ADB = glm::cross(AD, AB);
//
//        // Check which face the origin is on
//        if (glm::dot(ABC, AO) > 0) {
//            simplex = { C, B, A };
//            direction = ABC;
//        }
//        else if (glm::dot(ACD, AO) > 0) {
//            simplex = { D, C, A };
//            direction = ACD;
//        }
//        else if (glm::dot(ADB, AO) > 0) {
//            simplex = { B, D, A };
//            direction = ADB;
//        }
//        else {
//            return true;  // The origin is inside the tetrahedron
//        }
//
//        if (glm::length(direction) < epsilon) return true;  // Collision if direction is too small
//        return false;
//    }
//
//    return false;
//}
//
//// EPA
//
//std::pair<std::vector<vec4>, size_t> GetFaceNormals(
//    const std::vector<vec3>& polytope,
//    const std::vector<size_t>& faces)
//{
//    std::vector<vec4> normals;
//    size_t minTriangle = 0;
//    float  minDistance = FLT_MAX;
//
//    for (size_t i = 0; i < faces.size(); i += 3) {
//        vec3 a = polytope[faces[i]];
//        vec3 b = polytope[faces[i + 1]];
//        vec3 c = polytope[faces[i + 2]];
//
//        vec3 normal = glm::normalize(cross(b - a, c - a));
//        float distance = dot(normal, a);
//
//        if (distance < 0) {
//            normal *= -1;
//            distance *= -1;
//        }
//
//        normals.emplace_back(normal, distance);
//
//        if (distance < minDistance) {
//            minTriangle = i / 3;
//            minDistance = distance;
//        }
//    }
//
//    return { normals, minTriangle };
//}
//
//void AddIfUniqueEdge(
//    std::vector<std::pair<size_t, size_t>>& edges,
//    const std::vector<size_t>& faces,
//    size_t a,
//    size_t b)
//{
//    auto reverse = std::find(                       //      0--<--3
//        edges.begin(),                              //     / \ B /   A: 2-0
//        edges.end(),                                //    / A \ /    B: 0-2
//        std::make_pair(faces[b], faces[a]) //   1-->--2
//    );
//
//    if (reverse != edges.end()) {
//        edges.erase(reverse);
//    }
//
//    else {
//        edges.emplace_back(faces[a], faces[b]);
//    }
//}
//
//struct CollisionPoints {
//    glm::vec3 Normal;
//    double PenetrationDepth;
//    bool HasCollision;
//};
//
//CollisionPoints EPA(const std::vector<glm::vec3> simplex, const std::vector<glm::vec3> colliderA, const std::vector<glm::vec3> colliderB) {
//    std::vector<vec3> polytope(simplex.begin(), simplex.end());
//    std::vector<size_t> faces = {
//        0, 1, 2,
//        0, 3, 1,
//        0, 2, 3,
//        1, 3, 2
//    };
//
//    // list: vec4(normal, distance), index: min distance
//    auto [normals, minFace] = GetFaceNormals(polytope, faces);
//
//    vec3  minNormal;
//    float minDistance = FLT_MAX;
//
//    while (minDistance == FLT_MAX) {
//        minNormal = normals[minFace];
//        minDistance = normals[minFace].w;
//
//        vec3 support = MinkowskiDiffrence(colliderA, colliderB, minNormal);
//        float sDistance = dot(minNormal, support);
//
//        if (abs(sDistance - minDistance) > 0.001f) {
//            minDistance = FLT_MAX;
//            std::vector<std::pair<size_t, size_t>> uniqueEdges;
//
//            for (size_t i = 0; i < normals.size(); i++) {
//                if (glm::dot(vec3(normals[i]), support - polytope[faces[i * 3]]) > 0) {
//                    size_t f = i * 3;
//
//                    AddIfUniqueEdge(uniqueEdges, faces, f, f + 1);
//                    AddIfUniqueEdge(uniqueEdges, faces, f + 1, f + 2);
//                    AddIfUniqueEdge(uniqueEdges, faces, f + 2, f);
//
//                    faces[f + 2] = faces.back(); faces.pop_back();
//                    faces[f + 1] = faces.back(); faces.pop_back();
//                    faces[f] = faces.back(); faces.pop_back();
//
//                    normals[i] = normals.back(); // pop-erase
//                    normals.pop_back();
//
//                    i--;
//                }
//            }
//            std::vector<size_t> newFaces;
//            for (auto [edgeIndex1, edgeIndex2] : uniqueEdges) {
//                newFaces.push_back(edgeIndex1);
//                newFaces.push_back(edgeIndex2);
//                newFaces.push_back(polytope.size());
//            }
//
//            polytope.push_back(support);
//
//            auto [newNormals, newMinFace] = GetFaceNormals(polytope, newFaces);
//            float oldMinDistance = FLT_MAX;
//            for (size_t i = 0; i < normals.size(); i++) {
//                if (normals[i].w < oldMinDistance) {
//                    oldMinDistance = normals[i].w;
//                    minFace = i;
//                }
//            }
//
//            if (newNormals[newMinFace].w < oldMinDistance) {
//                minFace = newMinFace + normals.size();
//            }
//
//            faces.insert(faces.end(), newFaces.begin(), newFaces.end());
//            normals.insert(normals.end(), newNormals.begin(), newNormals.end());
//        }
//    }
//
//    CollisionPoints points;
//    points.Normal = minNormal;
//    points.PenetrationDepth = minDistance + 0.001f;
//    points.HasCollision = true;
//
//    return points;
//}
//
//bool RidigBody::GJK(const std::vector<glm::vec3>& shapeA, const std::vector<glm::vec3>& shapeB) {
//    auto direction = glm::vec3(1, 0, 0); // random
//    std::vector<glm::vec3> simplex = { MinkowskiDiffrence(shapeA, shapeB, direction) };
//    direction = -simplex[0];
//
//    while (true) {
//        // Get a new point in the direction of the origin
//        glm::vec3 A = MinkowskiDiffrence(shapeA, shapeB, direction);
//
//        // If the new point isn't past the origin in the given direction, there's no collision
//        if (glm::dot(A, direction) <= 0) {
//            printf("false\n");
//            return false;   // No collision
//        }
//
//        // Add the new point to the simplex
//        simplex.push_back(A);
//
//        // Check if the simplex contains the origin
//        if (DoSimplex(simplex, direction)) {
//            //printf("true\n");
//            Debug::DrawLine(simplex[0], simplex[1], 1.0f, glm::vec4(1, 1, 1, 1), glm::vec4(1, 1, 1, 1), Debug::RenderMode::Normal);
//            Debug::DrawLine(simplex[1], simplex[2], 1.0f, glm::vec4(1, 1, 1, 1), glm::vec4(1, 1, 1, 1), Debug::RenderMode::Normal);
//            Debug::DrawLine(simplex[2], simplex[0], 1.0f, glm::vec4(1, 1, 1, 1), glm::vec4(1, 1, 1, 1), Debug::RenderMode::Normal);
//            Debug::DrawLine(simplex[0], simplex[3], 1.0f, glm::vec4(1, 1, 1, 1), glm::vec4(1, 1, 1, 1), Debug::RenderMode::Normal);
//            Debug::DrawLine(simplex[1], simplex[3], 1.0f, glm::vec4(1, 1, 1, 1), glm::vec4(1, 1, 1, 1), Debug::RenderMode::Normal);
//            Debug::DrawLine(simplex[2], simplex[3], 1.0f, glm::vec4(1, 1, 1, 1), glm::vec4(1, 1, 1, 1), Debug::RenderMode::Normal);
//            vec3 epaResult = EPA(simplex, shapeA, shapeB).Normal;
//            // Print the collision normal
//            printf("Collision Normal: (%f, %f, %f)\n", epaResult.x, epaResult.y, epaResult.z);
//            Debug::DrawLine(vec3(0), epaResult, 1.0f, glm::vec4(1, 0.5, 1, 1), glm::vec4(1, 0.5, 1, 1), Debug::RenderMode::Normal);
//
//            return true;    // Collision detected
//        }
//    }
//}