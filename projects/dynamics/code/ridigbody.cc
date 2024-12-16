#include"ridigbody.h"

RidigBody::RidigBody(float mass, float restitutionCoefficient) {
    this->position = glm::vec3(0.0f);
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
    linearVelocity += linearAcceleration * deltaTime;
    position += linearVelocity * deltaTime;
    
    // Update angular motion
    angularAcceleration = netTorque / momentOfInertia;
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