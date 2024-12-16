//------------------------------------------------------------------------------
// spacegameapp.cc
// (C) 2022 Individual contributors, see AUTHORS file
//------------------------------------------------------------------------------
#include "config.h"
#include "spacegameapp.h"
#include <cstring>
#include "imgui.h"
#include "render/renderdevice.h"
#include "render/shaderresource.h"
#include <vector>
#include "render/textureresource.h"
#include "render/model.h"
#include "render/cameramanager.h"
#include "render/lightserver.h"
#include "render/debugrender.h"
#include "core/random.h"
#include "render/input/inputserver.h"
#include "core/cvar.h"
#include "render/physics.h"
#include <chrono>

#include "cameracontroller.h"
#include "raycast.h"
#include "ridigbody.h"

using namespace Display;
using namespace Render;

namespace Game
{

//------------------------------------------------------------------------------

SpaceGameApp::SpaceGameApp()
{

}

//------------------------------------------------------------------------------

SpaceGameApp::~SpaceGameApp()
{

}

//------------------------------------------------------------------------------

bool
SpaceGameApp::Open()
{
	App::Open();
	this->window = new Display::Window;
    this->window->SetSize(1600, 900);

    if (this->window->Open())
	{
		// set clear color to gray
		glClearColor(0.1f, 0.1f, 0.1f, 1.0f);

        RenderDevice::Init();

		// set ui rendering function
		this->window->SetUiRender([this]()
		{
			this->RenderUI();
		});
        
        return true;
	}
	return false;
}

//------------------------------------------------------------------------------

void
SpaceGameApp::Run()
{
    int w;
    int h;
    this->window->GetSize(w, h);
    glm::mat4 projection = glm::perspective(glm::radians(90.0f), float(w) / float(h), 0.01f, 1000.f);
    Camera* cam = CameraManager::GetCamera(CAMERA_MAIN);
    CameraController cameraController;
    cam->projection = projection;
    cam->view = glm::mat4(1);

    Plane plane(glm::vec3(0,0,0), glm::vec3(0.0f, 1.0f, 0.0f));

    // load all resources
    ModelId quad = LoadQuad(glm::vec4(0,1,0,1));
    ModelId asteroid = LoadModel("assets/space/Asteroid_1.glb");
    RidigBody rb(20, 20, glm::vec3(0));
    RidigBody rb2(20, 20, glm::vec3(5,0,0));
    Collider col(asteroid, rb.transform);
    Collider col2(asteroid, rb2.transform);
    
    // Setup skybox
    std::vector<const char*> skybox
    {
        "assets/space/bg.png",
        "assets/space/bg.png",
        "assets/space/bg.png",
        "assets/space/bg.png",
        "assets/space/bg.png",
        "assets/space/bg.png"
    };
    TextureResourceId skyboxId = TextureResource::LoadCubemap("skybox", skybox, true);
    RenderDevice::SetSkybox(skyboxId);
    
    Input::Keyboard* kbd = Input::GetDefaultKeyboard();

    std::clock_t c_start = std::clock();
    double dt = 0.01667f;

    // game loop
    while (this->window->IsOpen())
	{
        auto timeStart = std::chrono::steady_clock::now();
		glClear(GL_DEPTH_BUFFER_BIT);
		glEnable(GL_DEPTH_TEST);
		glEnable(GL_CULL_FACE);
		glCullFace(GL_BACK);
        
        this->window->Update();

        if (kbd->pressed[Input::Key::Code::End])
        {
            ShaderResource::ReloadShaders();
        }
        cameraController.Update(dt);
        plane.DrawDebugPlane(10);

        rb.Update(dt);
        col.UpdateCollider(rb.transform);
        rb2.Update(dt);

        RenderDevice::Draw(asteroid, rb.transform);
        RenderDevice::Draw(asteroid, rb2.transform);

        //std::vector<glm::vec3> shapeA;
        //std::vector<glm::vec3> shapeB;

        //for (const auto& mesh : Render::GetModel(asteroid).meshes) {
        //    for (const auto& primitive : mesh.primitives) {
        //        // Loop through the vertices (not indices)
        //        for (const auto& vertex : primitive.vertices) {
        //            // Transform the vertex by the first rigid body's transform
        //            glm::vec3 transformedVertexA = glm::vec3(rb.transform * glm::vec4(vertex, 1.0f));
        //            shapeA.push_back(transformedVertexA);
        //            // Transform the vertex by the second rigid body's transform
        //            glm::vec3 transformedVertexB = glm::vec3(rb2.transform * glm::vec4(vertex, 1.0f));
        //            shapeB.push_back(transformedVertexB);
        //        }
        //    }
        //}

        //rb.GJK(shapeA, shapeB);
        col.DetectCollision(col2);

       
        // Draw and update rays 
        for (Ray& r : cameraController.cameraRays)
        {
            Ray::RaycastHit hitresult = r.IntersectModel(Render::GetModel(asteroid), rb.transform);
            if (hitresult.hit)
            {
                rb.ApplyForceAtPoint(normalize(r.direction) * 10.0f, hitresult.point);
            }
            r.DrawDebugRay(hitresult,20);
        }
        cameraController.cameraRays.clear();

        // Execute the entire rendering pipeline
        RenderDevice::Render(this->window, dt);

		// transfer new frame to window
		this->window->SwapBuffers();

        auto timeEnd = std::chrono::steady_clock::now();
        dt = std::min(0.04, std::chrono::duration<double>(timeEnd - timeStart).count());

        if (kbd->pressed[Input::Key::Code::Escape])
            this->Exit();
	}
}

//------------------------------------------------------------------------------

void
SpaceGameApp::Exit()
{
    this->window->Close();
}

//------------------------------------------------------------------------------

void
SpaceGameApp::RenderUI()
{
	if (this->window->IsOpen())
	{
        ImGui::Begin("Debug");
        Core::CVar* r_draw_light_spheres = Core::CVarGet("r_draw_light_spheres");
        int drawLightSpheres = Core::CVarReadInt(r_draw_light_spheres);
        if (ImGui::Checkbox("Draw Light Spheres", (bool*)&drawLightSpheres))
            Core::CVarWriteInt(r_draw_light_spheres, drawLightSpheres);
        
        Core::CVar* r_draw_light_sphere_id = Core::CVarGet("r_draw_light_sphere_id");
        int lightSphereId = Core::CVarReadInt(r_draw_light_sphere_id);
        if (ImGui::InputInt("LightSphereId", (int*)&lightSphereId))
            Core::CVarWriteInt(r_draw_light_sphere_id, lightSphereId);
        
        ImGui::End();

        Debug::DispatchDebugTextDrawing();
	}
}

} // namespace Game