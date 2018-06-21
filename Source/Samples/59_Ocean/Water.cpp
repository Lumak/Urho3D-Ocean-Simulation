//=============================================================================
// Copyright (c) 2016 Lumak
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
//
//=============================================================================

#include <Urho3D/Urho3D.h>

#include <Urho3D/Graphics/Camera.h>
#include <Urho3D/Core/CoreEvents.h>
#include <Urho3D/Engine/Engine.h>
#include <Urho3D/IO/File.h>
#include <Urho3D/IO/FileSystem.h>
#include <Urho3D/UI/Font.h>
#include <Urho3D/Graphics/Graphics.h>
#include <Urho3D/Input/Input.h>
#include <Urho3D/Graphics/Light.h>
#include <Urho3D/Graphics/Material.h>
#include <Urho3D/Graphics/Model.h>
#include <Urho3D/Graphics/Octree.h>
#include <Urho3D/Scene/Scene.h>
#include <Urho3D/Graphics/Skybox.h>
#include <Urho3D/Graphics/StaticModel.h>
#include <Urho3D/Graphics/Terrain.h>
#include <Urho3D/UI/Text.h>
#include <Urho3D/Graphics/Texture2D.h>
#include <Urho3D/UI/UI.h>
#include <Urho3D/Graphics/Zone.h>
#include <SDL/SDL_log.h>
#include <Urho3D/Graphics/StaticModel.h>
#include <Urho3D/Graphics/DebugRenderer.h>

#include "Water.h"
#include "Ocean.h"

#include <Urho3D/DebugNew.h>

//=============================================================================
//=============================================================================
URHO3D_DEFINE_APPLICATION_MAIN(Water)

//=============================================================================
//=============================================================================
Water::Water(Context* context) 
    : Sample(context)
    , fpsCounter_(0)
    , m_dbgShow(false)
{
    // register objects
    DStaticModel::RegisterObject(context);
    Ocean::RegisterObject(context);
}

Water::~Water()
{
}

//=============================================================================
//=============================================================================
void Water::Setup()
{
    // Modify engine startup parameters
    engineParameters_["WindowTitle"]  = GetTypeName();
    engineParameters_["LogName"]      = GetSubsystem<FileSystem>()->GetAppPreferencesDir("urho3d", "logs") + GetTypeName() + ".log";
    engineParameters_["FullScreen"]   = false;
    engineParameters_["Headless"]     = false;
    engineParameters_["WindowWidth"]  = 1280; 
    engineParameters_["WindowHeight"] = 720;
}

//=============================================================================
//=============================================================================
void Water::Start()
{
    // Execute base class startup
    Sample::Start();

    // Create the scene content
    CreateScene();

    // Create the UI content
    CreateInstructions();

    // Setup the viewport for displaying the scene
    SetupViewport();

    // Hook up to the frame update event
    SubscribeToEvents();
}

//=============================================================================
//=============================================================================
void Water::CreateScene()
{
    ResourceCache* cache = GetSubsystem<ResourceCache>();

    scene_ = new Scene(context_);

    m_pDbgRenderer = scene_->CreateComponent<DebugRenderer>();

    // Create octree, use default volume (-1000, -1000, -1000) to (1000, 1000, 1000)
    scene_->CreateComponent<Octree>();

    // Create a Zone component for ambient lighting & fog control
    Node* zoneNode = scene_->CreateChild("Zone");
    Zone* zone = zoneNode->CreateComponent<Zone>();
    zone->SetBoundingBox(BoundingBox(-2000.0f, 2000.0f));
    zone->SetAmbientColor(Color(0.3f, 0.3f, 0.3f));
    zone->SetFogColor(Color(0.1f, 0.2f, 0.3f));
    zone->SetFogStart(800.0f);
    zone->SetFogEnd(1000.0f);

    // Create a directional light to the world. Enable cascaded shadows on it
    Node* lightNode = scene_->CreateChild("DirectionalLight");
    lightNode->SetDirection(Vector3(0.6f, -1.0f, 0.8f));
    Light* light = lightNode->CreateComponent<Light>();
    light->SetLightType(LIGHT_DIRECTIONAL);
    light->SetCastShadows(true);
    light->SetShadowBias(BiasParameters(0.00025f, 0.5f));
    light->SetShadowCascade(CascadeParameters(10.0f, 50.0f, 200.0f, 0.0f, 0.8f));
    light->SetSpecularIntensity(0.1f);
    // Apply slightly overbright lighting to match the skybox
    light->SetColor(Color(1.0f, 1.0f, 1.0f));

    // Create skybox. The Skybox component is used like StaticModel, but it will be always located at the camera, giving the
    // illusion of the box planes being far away. Use just the ordinary Box model and a suitable material, whose shader will
    // generate the necessary 3D texture coordinates for cube mapping
    Node* skyNode = scene_->CreateChild("Sky");
    skyNode->SetScale(500.0f); // The scale actually does not matter
    Skybox* skybox = skyNode->CreateComponent<Skybox>();
    skybox->SetModel(cache->GetResource<Model>("Models/Box.mdl"));
    skybox->SetMaterial(cache->GetResource<Material>("Materials/Skybox.xml"));

    // Create heightmap terrain
    Node* terrainNode = scene_->CreateChild("Terrain");
    terrainNode->SetPosition(Vector3(0.0f, -4.0f, 0.0f));
    Terrain* terrain = terrainNode->CreateComponent<Terrain>();
    terrain->SetPatchSize(32);
    terrain->SetSpacing(Vector3(2.2f, 0.4f, 2.2f)); // Spacing between vertices and vertical resolution of the height map
    //terrain->SetSmoothing(true);
    terrain->SetHeightMap( cache->GetResource<Image>("Ocean/HeightMap-001.png") );
    terrain->SetMaterial(cache->GetResource<Material>("Ocean/TerrainGeneric.xml"));

    waterNode_ = scene_->CreateChild("Water");
    waterNode_->SetScale(Vector3(2048.0f, 1.0f, 2048.0f));
    waterNode_->SetPosition(Vector3(0.0f, 5.0f, 0.0f));

    cameraNode_ = scene_->CreateChild("camNode");
    Camera* camera = cameraNode_->CreateComponent<Camera>();
    camera->SetFarClip(4000.0f);

    // Set an initial position for the camera scene node above the ground
    cameraNode_->SetPosition(Vector3(0.0f, 67.0f, -200.0f));

    oceanNode_ = scene_->CreateChild( "Ocean" );
    oceanNode_->SetPosition(Vector3(0.0f, 10.0f, 0.0f));
    oceanNode_->SetScale(Vector3(2.0f, 2.0f, 2.0f));

    // create and start
    m_pOcean = oceanNode_->CreateComponent<Ocean>();
    m_pOcean->InitOcean();

    m_pStaticModelOcean = oceanNode_->CreateComponent<DStaticModel>();
    m_pStaticModelOcean->SetModel( m_pOcean->GetOceanModel() );
    m_pStaticModelOcean->SetMaterial(cache->GetResource<Material>("Ocean/MatOcean.xml"));
    m_pStaticModelOcean->SetViewMask(0x80000000);
}

//=============================================================================
//=============================================================================
void Water::CreateInstructions()
{
    Graphics* graphics = GetSubsystem<Graphics>();
    ResourceCache* cache = GetSubsystem<ResourceCache>();
    UI* ui = GetSubsystem<UI>();
    
    // Construct new Text object, set string to display and font to use
    Text* instructionText = ui->GetRoot()->CreateChild<Text>();
    instructionText->SetText("Use WASD keys and mouse/touch to move");
    instructionText->SetFont(cache->GetResource<Font>("Fonts/Anonymous Pro.ttf"), 15);
    instructionText->SetTextAlignment(HA_CENTER);
    
    // Position the text relative to the screen center
    instructionText->SetHorizontalAlignment(HA_CENTER);
    instructionText->SetVerticalAlignment(VA_CENTER);
    instructionText->SetPosition(0, ui->GetRoot()->GetHeight() / 4);

    // fps text
    fpsText_  = ui->GetRoot()->CreateChild<Text>();
    fpsText_->SetFont(cache->GetResource<Font>("Fonts/Anonymous Pro.ttf"), 15);
    fpsText_->SetPosition(graphics->GetWidth() - 100, 5);
    fpsText_->SetColor(Color::BLACK);
    fpsText_->SetText("60");
}

//=============================================================================
//=============================================================================
void Water::SetupViewport()
{
    Graphics* graphics = GetSubsystem<Graphics>();
    Renderer* renderer = GetSubsystem<Renderer>();
    ResourceCache* cache = GetSubsystem<ResourceCache>();

    // Set up a viewport to the Renderer subsystem so that the 3D scene can be seen
    SharedPtr<Viewport> viewport(new Viewport(context_, scene_, cameraNode_->GetComponent<Camera>()));
    renderer->SetViewport(0, viewport);
#ifdef ADD_WATER_REFLECTION
    // Create a mathematical plane to represent the water in calculations
    waterPlane_ = Plane(waterNode_->GetWorldRotation() * Vector3(0.0f, 1.0f, 0.0f), waterNode_->GetWorldPosition());
    // Create a downward biased plane for reflection view clipping. Biasing is necessary to avoid too aggressive clipping
    waterClipPlane_ = Plane(waterNode_->GetWorldRotation() * Vector3(0.0f, 1.0f, 0.0f), waterNode_->GetWorldPosition() -
        Vector3(0.0f, 0.1f, 0.0f));

    // Create camera for water reflection
    // It will have the same farclip and position as the main viewport camera, but uses a reflection plane to modify
    // its position when rendering
    reflectionCameraNode_ = cameraNode_->CreateChild();
    Camera* reflectionCamera = reflectionCameraNode_->CreateComponent<Camera>();
    reflectionCamera->SetFarClip(750.0);
    reflectionCamera->SetViewMask(0x7fffffff); // Hide objects with only bit 31 in the viewmask (the water plane)
    reflectionCamera->SetAutoAspectRatio(false);
    reflectionCamera->SetUseReflection(true);
    reflectionCamera->SetReflectionPlane(waterPlane_);
    reflectionCamera->SetUseClipping(true); // Enable clipping of geometry behind water plane
    reflectionCamera->SetClipPlane(waterClipPlane_);
    // The water reflection texture is rectangular. Set reflection camera aspect ratio to match
    reflectionCamera->SetAspectRatio((float)graphics->GetWidth() / (float)graphics->GetHeight());
    // View override flags could be used to optimize reflection rendering. For example disable shadows
    //reflectionCamera->SetViewOverrideFlags(VO_DISABLE_SHADOWS);

    // Create a texture and setup viewport for water reflection. Assign the reflection texture to the diffuse
    // texture unit of the water material
    int texSize = 1024;
    SharedPtr<Texture2D> renderTexture(new Texture2D(context_));
    renderTexture->SetSize(texSize, texSize, Graphics::GetRGBFormat(), TEXTURE_RENDERTARGET);
    renderTexture->SetFilterMode(FILTER_BILINEAR);
    RenderSurface* surface = renderTexture->GetRenderSurface();
    SharedPtr<Viewport> rttViewport(new Viewport(context_, scene_, reflectionCamera));
    surface->SetViewport(0, rttViewport);
    Material* waterMat = cache->GetResource<Material>("Materials/Water.xml");
    waterMat->SetTexture(TU_DIFFUSE, renderTexture);
#endif
}

//=============================================================================
//=============================================================================
void Water::SubscribeToEvents()
{
    SubscribeToEvent(E_UPDATE, URHO3D_HANDLER(Water, HandleUpdate));
}

//=============================================================================
//=============================================================================
void Water::HandleUpdate(StringHash eventType, VariantMap& eventData)
{
    using namespace Update;

    // Take the frame time step, which is stored as a float
    float timeStep = eventData[P_TIMESTEP].GetFloat();

    // Move the camera, scale movement with time step
    MoveCamera(timeStep);

    BoundingBox bbox = m_pOcean->GetBoundingBox();

    if ( (bbox.Size() - m_boundingbox.Size()).Length() > 10.0f )
    {
        m_boundingbox.Merge( bbox );
        m_pStaticModelOcean->DSetBoundingBox( m_boundingbox );
    }

    // fps text
    fpsCounter_++;

    if ( timerFps_.GetMSec(false) >= 1000 )
    {
        fpsText_->SetText(String("fps: ") + String(fpsCounter_));
        fpsCounter_ = 0;
        timerFps_.Reset();
    }

    Input* input = GetSubsystem<Input>();

    if ( input->GetKeyPress( KEY_F7 ) )
    {
        m_dbgShow = !m_dbgShow;
    }

    if ( m_dbgShow )
    {
        m_pOcean->DbgRender();
    }

}

//=============================================================================
//=============================================================================
void Water::MoveCamera(float timeStep)
{
    // Do not move if the UI has a focused element (the console)
    if (GetSubsystem<UI>()->GetFocusElement())
        return;

    Input* input = GetSubsystem<Input>();

    // Movement speed as world units per second
    const float MOVE_SPEED = 90.0f;
    // Mouse sensitivity as degrees per pixel
    const float MOUSE_SENSITIVITY = 0.1f;

    // Use this frame's mouse motion to adjust camera node yaw and pitch. Clamp the pitch between -90 and 90 degrees
    IntVector2 mouseMove = input->GetMouseMove();
    yaw_ += MOUSE_SENSITIVITY * mouseMove.x_;
    pitch_ += MOUSE_SENSITIVITY * mouseMove.y_;
    pitch_ = Clamp(pitch_, -90.0f, 90.0f);

    // Construct new orientation for the camera scene node from yaw and pitch. Roll is fixed to zero
    cameraNode_->SetRotation(Quaternion(pitch_, yaw_, 0.0f));

    // Read WASD keys and move the camera scene node to the corresponding direction if they are pressed
    if (input->GetKeyDown('W'))
        cameraNode_->Translate(Vector3::FORWARD * MOVE_SPEED * timeStep);
    if (input->GetKeyDown('S'))
        cameraNode_->Translate(Vector3::BACK * MOVE_SPEED * timeStep);
    if (input->GetKeyDown('A'))
        cameraNode_->Translate(Vector3::LEFT * MOVE_SPEED * timeStep);
    if (input->GetKeyDown('D'))
        cameraNode_->Translate(Vector3::RIGHT * MOVE_SPEED * timeStep);

#ifdef ADD_WATER_REFLECTION
    // In case resolution has changed, adjust the reflection camera aspect ratio
    Graphics* graphics = GetSubsystem<Graphics>();
    Camera* reflectionCamera = reflectionCameraNode_->GetComponent<Camera>();
    reflectionCamera->SetAspectRatio((float)graphics->GetWidth() / (float)graphics->GetHeight());
#endif
}



