#include "cameracontroller.h"
#include "config.h"
#include "render/cameramanager.h"
#include "render/input/inputserver.h"
#include <iostream>

using namespace Input;
using namespace glm;
using namespace Render;

namespace Game
{
    void CameraController::Update(float dt)
    {
		Mouse* mouse = GetDefaultMouse();
		Keyboard* kbd = GetDefaultKeyboard();
		Camera* cam = CameraManager::GetCamera(CAMERA_MAIN);

		// Update camera rotation
		if(mouse->held[1]){
			glm::vec2 mouseDelta = mouse->delta;

			yaw += mouseDelta.x * sensitivity;
			pitch -= mouseDelta.y * sensitivity;
			pitch = glm::clamp(pitch, -89.99f, 89.99f);

			forward.x = cos(glm::radians(yaw)) * cos(glm::radians(pitch));
			forward.y = sin(glm::radians(pitch));
			forward.z = sin(glm::radians(yaw)) * cos(glm::radians(pitch));
			right = glm::cross(forward, glm::vec3(0, 1, 0));
		}
		if (mouse->pressed[0]) {
			SpawnRay(mouse->position);
		}

		// Update camera position
		if (kbd->held[Key::W])
		{
			this->position += forward * speed * dt;
		}
		if (kbd->held[Key::A])
		{
			this->position -= right * speed * dt;
		}
		if (kbd->held[Key::S])
		{
			this->position -= forward * speed * dt;
		}
		if (kbd->held[Key::D])
		{
			this->position += right * speed * dt;
		}
		
		// update camera view transform
		cam->view = glm::lookAt(this->position, this->position + forward, glm::vec3(0, 1, 0));
    }

	void CameraController::SpawnRay(glm::vec2 screenPos)
	{
		Camera* cam = CameraManager::GetCamera(CAMERA_MAIN);

		glm::vec2 ndcPos((2.0f * (screenPos.x / 1600)) - 1.0f, 1.0f - (2.0f * (screenPos.y / 900)));
		// Debug print the NDC position.
		//std::cout << "Screen Position: " << screenPos.x << ", " << screenPos.y << std::endl;
		//std::cout << "NDC Position: " << ndcPos.x << ", " << ndcPos.y << std::endl;

		glm::vec4 clipCoords = glm::vec4(ndcPos.x, ndcPos.y, -1.0f, 1.0f);
		glm::vec4 eyeCoords = cam->invProjection * clipCoords;
		eyeCoords.z = -1.0f;
		eyeCoords.w = 0.0f;
		glm::vec3 dir = glm::normalize(cam->invView * eyeCoords);

		// Add the ray to the camera's rays
		this->cameraRays.push_back(Ray(this->position, dir));
		if (this->cameraRays.size() > 5) {
			this->cameraRays.erase(this->cameraRays.begin());
		}
	}
}