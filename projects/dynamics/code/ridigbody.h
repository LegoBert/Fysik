#pragma once
#include "config.h"
#include "render/model.h"
#include "render/debugrender.h"

using namespace glm;
using namespace std;

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
    /*glm::vec3 minAABB;
    glm::vec3 maxAABB;*/
public:
    RidigBody(float mass, float restitutionCoefficient);
    void ApplyForceAtCenter(glm::vec3 force);
    void ApplyForceAtPoint(const glm::vec3& force, const glm::vec3& pointOfContact);
    void Update(float deltaTime);
};