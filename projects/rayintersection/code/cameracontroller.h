#pragma once
#include "config.h"
#include "raycast.h"
#include <vector>
using namespace glm;

namespace Game
{

struct CameraController
{
    glm::vec3 position = glm::vec3(0);
    glm::vec3 forward = glm::vec3(0.0f, 0.0f, 1.0f);
    glm::vec3 right = glm::cross(forward, vec3(0.0f, 1.0f, 0.0f));

    float sensitivity = 0.2f;
    float speed = 10.0f;
    float yaw = 90;
    float pitch = 0;

    std::vector<Ray> cameraRays;

    void Update(float dt);
    void SpawnRay(glm::vec2 screenPos);
};

}