#pragma once
//------------------------------------------------------------------------------
/**
    Render::Model

    (C) 2022 Individual contributors, see AUTHORS file
*/
//------------------------------------------------------------------------------
#include "GL/glew.h"
#include <string>
#include <vector>
#include "renderdevice.h"
#include "resourceid.h"
#include "textureresource.h"

namespace Render
{

struct Model
{   
    struct VertexAttribute
    {
        GLuint slot = 0;
        GLint components = 0;
        GLenum type = GL_NONE;
        GLsizei stride = 0;
        GLsizei offset = 0;
        GLboolean normalized = GL_FALSE;
    };

    struct Material
    {
        enum
        {
            TEXTURE_BASECOLOR,
            TEXTURE_NORMAL,
            TEXTURE_METALLICROUGHNESS,
            TEXTURE_EMISSIVE,
            TEXTURE_OCCLUSION,
            NUM_TEXTURES
        };

        glm::vec4 baseColorFactor = glm::vec4(1.0f);
        glm::vec4 emissiveFactor = glm::vec4(1.0f);
        float metallicFactor = 1.0f;
        float roughnessFactor = 1.0f;
        TextureResourceId textures[NUM_TEXTURES] = {
            InvalidResourceId,
            InvalidResourceId,
            InvalidResourceId,
            InvalidResourceId,
            InvalidResourceId
        };

        enum class AlphaMode : uint8_t
        {
            Opaque,
            Mask,
            Blend
        };

        float alphaCutoff{ 0.5f };
        AlphaMode alphaMode{ AlphaMode::Opaque };
        bool doubleSided{ false };

        static Material CreateCustomMaterial(glm::vec4& rgba)
        {
            Material material;
            material.baseColorFactor = rgba;//glm::vec4(1.0f, 0.0f, 0.0f, 1.0f); // RGBA
            material.metallicFactor = 0.0f;
            material.roughnessFactor = 0.5f; // Adjust this value as needed.
            material.alphaMode = Material::AlphaMode::Opaque;
            material.textures[Model::Material::TEXTURE_BASECOLOR] = TextureResource::GetWhiteTexture();
            material.textures[Model::Material::TEXTURE_NORMAL] = TextureResource::GetDefaultNormalTexture();
            material.textures[Model::Material::TEXTURE_METALLICROUGHNESS] = TextureResource::GetDefaultMetallicRoughnessTexture();
            material.textures[Model::Material::TEXTURE_EMISSIVE] = TextureResource::GetBlackTexture();
            material.textures[Model::Material::TEXTURE_OCCLUSION] = TextureResource::GetBlackTexture();
            return material;
        }
    };

    struct Mesh
    {
        struct Primitive
        {
            GLuint vao;
            GLuint vbo;
            GLuint ebo;
            GLuint numIndices;
            GLuint offset = 0;
            GLenum indexType;
            Material material;
            std::vector<glm::vec3> vertices;
            std::vector<GLuint> indices;
        };

        std::vector<Primitive> primitives;
        std::vector<uint16_t> opaquePrimitives; // contains ids of all primitives to be rendered opaque or mask mode
        std::vector<uint16_t> blendPrimitives; // contains ids of all primitives to be rendered with blend mode

        void AddQuad(glm::vec4& rgba) {
            Primitive quadPrimitive;
            quadPrimitive.material = Material::CreateCustomMaterial(rgba);

            // Define vertices (only positions, for simplicity)
            quadPrimitive.vertices = {
                glm::vec3(-0.5f, -0.5f, 0.0f),
                glm::vec3(0.5f, -0.5f, 0.0f),
                glm::vec3(0.5f, 0.5f, 0.0f),
                glm::vec3(-0.5f, 0.5f, 0.0f)
            };

            // Define indices
            quadPrimitive.indices = { 0, 1, 2, 0, 2, 3 };

            quadPrimitive.numIndices = quadPrimitive.indices.size();
            quadPrimitive.indexType = GL_UNSIGNED_INT;

            // Generate and bind the VAO
            glGenVertexArrays(1, &quadPrimitive.vao);
            glBindVertexArray(quadPrimitive.vao);

            // Generate and bind the VBO
            glGenBuffers(1, &quadPrimitive.vbo);
            glBindBuffer(GL_ARRAY_BUFFER, quadPrimitive.vbo);
            glBufferData(GL_ARRAY_BUFFER, quadPrimitive.vertices.size() * sizeof(glm::vec3), quadPrimitive.vertices.data(), GL_STATIC_DRAW);

            // Generate and bind the EBO
            glGenBuffers(1, &quadPrimitive.ebo);
            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, quadPrimitive.ebo);
            glBufferData(GL_ELEMENT_ARRAY_BUFFER, quadPrimitive.indices.size() * sizeof(GLuint), quadPrimitive.indices.data(), GL_STATIC_DRAW);

            // Set the vertex attribute pointers
            // Position attribute
            glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), (GLvoid*)0);
            glEnableVertexAttribArray(0);

            // Unbind VAO
            glBindVertexArray(0);

            // Store the Primitive in the mesh
            primitives.push_back(quadPrimitive);
            opaquePrimitives.push_back(static_cast<uint16_t>(primitives.size() - 1));
        }
    };

    struct AABB
    {
        glm::vec3 max;
        glm::vec3 min;
    };

    AABB aabb;
    std::vector<Mesh> meshes;
    std::vector<GLuint> buffers;
    uint refcount;
    //std::vector<TextureResourceId> textures;
};

ModelId LoadModel(std::string name);

void UnloadModel(ModelId);

bool const IsModelValid(ModelId);

Model const& GetModel(ModelId id);

ModelId LoadQuad(glm::vec4 rgba);

} // namespace Render
