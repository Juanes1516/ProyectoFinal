/*
 * Copyright 2019-2022 Diligent Graphics LLC
 * Copyright 2015-2019 Egor Yusov
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 * In no event and under no legal theory, whether in tort (including negligence),
 * contract, or otherwise, unless required by applicable law (such as deliberate
 * and grossly negligent acts) or agreed to in writing, shall any Contributor be
 * liable for any damages, including any direct, indirect, special, incidental,
 * or consequential damages of any character arising as a result of this License or
 * out of the use or inability to use the software (including but not limited to damages
 * for loss of goodwill, work stoppage, computer failure or malfunction, or any and
 * all other commercial damages or losses), even if such Contributor has been advised
 * of the possibility of such damages.
 */

#include "Tutorial01_HelloTriangle.hpp"
#include "MapHelper.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include "GraphicsTypes.h"
#include "RefCntAutoPtr.hpp"
#include "RenderDevice.h"
#include "DeviceContext.h"
#include "SwapChain.h"
#include "BasicMath.hpp"
#include "CommonDefinitions.h"
#include "imgui.h"
#include "GraphicsUtilities.h"

#ifndef PI_F
#    define PI_F 3.14159265358979323846f
#endif

namespace Diligent
{

SampleBase* CreateSample()
{
    return new Tutorial01_HelloTriangle();
}

// Vertex Shader para 3D con Normales
static const char* VSSource_Droplets_3D_Shaded = R"(
struct VSInput
{
    float3 Pos    : ATTRIB0; 
    float3 Normal : ATTRIB1; // Normal del vértice
    float3 Color  : ATTRIB2; // Color base del vértice
};

struct PSInput 
{ 
    float4 PosWVP : SV_POSITION; // Posición transformada por MVP
    float3 NormalW: NORMAL;      // Normal en espacio del mundo (o vista, según se transforme)
    float3 Color  : COLOR0;      // Color base de la gota
    // float3 PosW   : TEXCOORD0; // Posición en espacio del mundo (opcional para efectos más avanzados)
};

cbuffer VSConstants
{
    matrix g_MVP;
    // Si se necesitara una matriz específica para normales (ej. inversa transpuesta de ModelView):
    // matrix g_NormalMatrix; 
};

void main(in  VSInput VSIn,
          out PSInput PSIn) 
{
    PSIn.PosWVP = mul(float4(VSIn.Pos, 1.0), g_MVP); 
    
    // Para esferas uniformemente escaladas y no rotadas individualmente,
    // la normal local (que es la dirección desde el centro) es suficiente
    // si la iluminación se calcula en espacio del mundo y la matriz de modelo es solo traslación/escala uniforme.
    // Si hubiera rotaciones o escalas no uniformes por objeto, se necesitaría g_NormalMatrix.
    // Por ahora, asumimos que VSIn.Normal ya está correctamente orientada para el espacio del mundo.
    PSIn.NormalW = normalize(VSIn.Normal); 
    
    PSIn.Color   = VSIn.Color;
    // PSIn.PosW    = VSIn.Pos; // Se puede pasar la posición en mundo si se necesita para el PS
}
)";

// Pixel Shader para 3D con Iluminación Difusa
static const char* PSSource_Droplets_3D_Shaded = R"(
struct PSInput 
{ 
    float4 PosWVP : SV_POSITION; 
    float3 NormalW: NORMAL;      
    float3 Color  : COLOR0;      
    // float3 PosW   : TEXCOORD0;   
};

struct PSOutput
{ 
    float4 Color : SV_TARGET; 
};

void main(in  PSInput  PSIn,
          out PSOutput PSOut)
{
    float3 lightDirection = normalize(float3(0.577, 0.577, -0.577)); // Dirección de luz fija (ej. desde arriba-derecha-atrás)
    float3 ambientColor   = float3(0.2, 0.2, 0.25);                 // Luz ambiental suave
    float3 lightColor     = float3(0.8, 0.8, 0.7);                  // Color de la luz direccional (blanco cálido)

    float3 normal = normalize(PSIn.NormalW); // Asegurar que la normal interpolada esté normalizada

    float diffuseFactor = max(0.0, dot(normal, lightDirection));
    
    float3 finalColor = (ambientColor + lightColor * diffuseFactor) * PSIn.Color;
    
    PSOut.Color = float4(finalColor, 1.0); 
}
)";

Tutorial01_HelloTriangle::DimensionlessNumbers Tutorial01_HelloTriangle::CalculateDimensionlessGroups(
    double Qc,
    double Qd,
    double mu_c,
    double mu_d,
    double gamma,
    double w,
    double h)
{
    double Qc_SI = Qc * 1e-9 / 60.0;
    double Qd_SI = Qd * 1e-9 / 60.0;
    double w_SI  = w * 1e-6;
    double h_SI  = (h < 0) ? w_SI : h * 1e-6;

    double A   = w_SI * h_SI;
    double U_c = Qc_SI / A;

    double mu_c_SI  = mu_c / 1000.0;
    double mu_d_SI  = mu_d / 1000.0;
    double gamma_SI = gamma / 1000.0;

    DimensionlessNumbers nums;
    nums.Ca_c     = mu_c_SI * U_c / gamma_SI;
    nums.Phi      = Qd / Qc;
    nums.lambda_v = mu_d / mu_c;
    nums.w_SI     = w_SI;
    nums.U_c      = U_c;
    return nums;
}

void Tutorial01_HelloTriangle::PredictDropletMetrics(const DimensionlessNumbers& dimensionless, double w_um)
{
    double Ca_c      = dimensionless.Ca_c;
    double Phi       = dimensionless.Phi;
    double lam       = dimensionless.lambda_v;
    double w_microns = w_um;

    double D_base_um;
    double AR_base;
    int    regime;

    if (Ca_c < 0.01)
    {
        regime          = 0;
        double beta     = 1.0;
        double L_over_w = 1.0 + beta * Phi;
        D_base_um       = L_over_w * w_microns;
        AR_base         = 1.05;
    }
    else if (Ca_c <= 0.02)
    {
        regime   = 1;
        double K = 1.0;
        double a = 0.33, b = -0.20, c = 0.10;
        D_base_um = K * w_microns * std::pow(Phi, a) * std::pow(Ca_c, b) * std::pow(lam, c);
        AR_base   = 1.15 + 0.3 * std::pow(lam, 0.25);
    }
    else
    {
        regime    = 2;
        D_base_um = 0.8 * w_microns * std::pow(Ca_c, -0.50);
        AR_base   = 1.0 + 0.1 * std::pow(lam, 0.2);
    }

    D_base_um = std::max(D_base_um, w_microns * 0.1);
    D_base_um = std::min(D_base_um, w_microns * 2.0);
    AR_base   = std::max(AR_base, 0.5);
    AR_base   = std::min(AR_base, 3.0);

    DropletInfo newDroplet;
    newDroplet.DiameterModel = static_cast<float>(D_base_um);
    newDroplet.AspectRatio   = static_cast<float>(AR_base);
    newDroplet.Regime        = regime;
    m_Droplets.push_back(newDroplet);
}

void Tutorial01_HelloTriangle::CreateUniformBuffers()
{
    CreateUniformBuffer(m_pDevice, sizeof(float4x4), "VS Constants CB", &m_pVSConstants);
}


void Tutorial01_HelloTriangle::PrepareDropletData()
{
    m_Droplets.clear();

    double Qc_val    = static_cast<double>(m_Qc_ul_min);
    double Qd_val    = static_cast<double>(m_Qd_ul_min);
    double gamma_val = static_cast<double>(m_gamma_mN_m);
    double mu_c_val  = static_cast<double>(m_mu_c_cP);
    double mu_d_val  = static_cast<double>(m_mu_d_cP);
    double w_um_val  = static_cast<double>(m_w_um);
    double h_um_val  = static_cast<double>(m_h_um);

    const int num_droplets_to_generate = 10;

    float initial_x_world          = -1.5f;
    float current_x_position_world = initial_x_world;

    float spacing_base_model = m_w_um * 0.5f;
    m_spacing_world          = spacing_base_model * m_world_scale_factor;

    float base_velocity_world = 0.5f;

    for (int i = 0; i < num_droplets_to_generate; ++i)
    {
        double current_Qd_sim = Qd_val;

        DimensionlessNumbers d_nums = CalculateDimensionlessGroups(Qc_val, current_Qd_sim, mu_c_val, mu_d_val, gamma_val, w_um_val, h_um_val);
        PredictDropletMetrics(d_nums, w_um_val);

        if (m_Droplets.size() > static_cast<size_t>(i))
        {
            DropletInfo& current_droplet = m_Droplets[i];

            current_droplet.RadiusWorld = (current_droplet.DiameterModel / 2.0f) * m_world_scale_factor;

            current_droplet.CurrentPositionWorld = float3(current_x_position_world, 0.0f, 0.0f);

            current_droplet.VelocityX_World = base_velocity_world;
            current_droplet.Concentration   = static_cast<float>(i % 5) / 4.0f;

            if (current_droplet.Regime == 0) current_droplet.Color = float3(0.8f, 0.2f, 0.2f);
            else if (current_droplet.Regime == 1)
                current_droplet.Color = float3(0.2f, 0.8f, 0.2f);
            else
                current_droplet.Color = float3(0.2f, 0.2f, 0.8f);

            current_x_position_world += (current_droplet.DiameterModel * m_world_scale_factor) + m_spacing_world;
        }
    }

    if (!m_Droplets.empty())
    {
        m_next_wrap_around_spawn_x_world = initial_x_world - (m_Droplets[0].DiameterModel * m_world_scale_factor) - m_spacing_world;
    }
    else
    {
        m_next_wrap_around_spawn_x_world = -2.0f;
    }
}


void Tutorial01_HelloTriangle::CreateDropletGeometryAndBuffers()
{
    m_DropletVertices.clear();
    const int latitudeBands  = 10;
    const int longitudeBands = 10;

    for (const auto& droplet : m_Droplets)
    {
        float radius = droplet.RadiusWorld;
        if (radius <= 1e-5f) continue;

        float radiusY = radius / droplet.AspectRatio;
        float radiusZ = radius / droplet.AspectRatio;

        float3 baseColor         = droplet.Color;
        float  intensity         = droplet.Concentration;
        float3 finalDropletColor = baseColor * (0.3f + 0.7f * intensity); // Este color se modulará con la iluminación en el shader

        for (int latNumber = 0; latNumber < latitudeBands; ++latNumber)
        {
            for (int longNumber = 0; longNumber < longitudeBands; ++longNumber)
            {
                float theta1 = static_cast<float>(latNumber) * PI_F / static_cast<float>(latitudeBands);
                float theta2 = static_cast<float>(latNumber + 1) * PI_F / static_cast<float>(latitudeBands);

                float phi1 = static_cast<float>(longNumber) * 2.0f * PI_F / static_cast<float>(longitudeBands);
                float phi2 = static_cast<float>(longNumber + 1) * 2.0f * PI_F / static_cast<float>(longitudeBands);

                // Vértices locales (centrados en el origen de la esfera)
                float3 p0_local_pos = float3(std::cos(phi1) * std::sin(theta1) * radius,
                                             std::cos(theta1) * radiusY,
                                             std::sin(phi1) * std::sin(theta1) * radiusZ);
                float3 p1_local_pos = float3(std::cos(phi2) * std::sin(theta1) * radius,
                                             std::cos(theta1) * radiusY,
                                             std::sin(phi2) * std::sin(theta1) * radiusZ);
                float3 p2_local_pos = float3(std::cos(phi1) * std::sin(theta2) * radius,
                                             std::cos(theta2) * radiusY,
                                             std::sin(phi1) * std::sin(theta2) * radiusZ);
                float3 p3_local_pos = float3(std::cos(phi2) * std::sin(theta2) * radius,
                                             std::cos(theta2) * radiusY,
                                             std::sin(phi2) * std::sin(theta2) * radiusZ);

                // Normales: para una esfera, son las posiciones locales normalizadas
                
                float3 n0 = normalize(p0_local_pos);
                float3 n1 = normalize(p1_local_pos);
                float3 n2 = normalize(p2_local_pos);
                float3 n3 = normalize(p3_local_pos);

                // Añadir dos triángulos para el quad, con sus normales
                m_DropletVertices.push_back({p0_local_pos + droplet.CurrentPositionWorld, n0, finalDropletColor});
                m_DropletVertices.push_back({p2_local_pos + droplet.CurrentPositionWorld, n2, finalDropletColor});
                m_DropletVertices.push_back({p1_local_pos + droplet.CurrentPositionWorld, n1, finalDropletColor});

                m_DropletVertices.push_back({p1_local_pos + droplet.CurrentPositionWorld, n1, finalDropletColor});
                m_DropletVertices.push_back({p2_local_pos + droplet.CurrentPositionWorld, n2, finalDropletColor});
                m_DropletVertices.push_back({p3_local_pos + droplet.CurrentPositionWorld, n3, finalDropletColor});
            }
        }
    }

    m_TotalVertexCount = static_cast<Uint32>(m_DropletVertices.size());

    if (m_TotalVertexCount == 0)
    {
        if (m_pDropletVertexBuffer) m_pDropletVertexBuffer.Release();
        return;
    }

    if (m_pDropletVertexBuffer) m_pDropletVertexBuffer.Release();

    BufferDesc VertBuffDesc;
    VertBuffDesc.Name      = "Droplet Vertex Buffer 3D Shaded";
    VertBuffDesc.Usage     = USAGE_IMMUTABLE;
    VertBuffDesc.BindFlags = BIND_VERTEX_BUFFER;
    VertBuffDesc.Size      = m_TotalVertexCount * sizeof(DropletVertex);

    BufferData VBData;
    VBData.pData    = m_DropletVertices.data();
    VBData.DataSize = VertBuffDesc.Size;

    m_pDevice->CreateBuffer(VertBuffDesc, &VBData, &m_pDropletVertexBuffer);
    if (!m_pDropletVertexBuffer)
    {
        LOG_ERROR_MESSAGE("Failed to create droplet vertex buffer.");
    }
}

void Tutorial01_HelloTriangle::UpdateUI()
{
    ImGui::Begin("Parámetros de Simulación de Gotas");

    ImGui::Text("Flujos (uL/min):");
    if (ImGui::SliderFloat("Qc (fase continua)", &m_Qc_ul_min, 1.0f, 500.0f, "%.1f")) m_simulationParametersChanged = true;
    if (ImGui::SliderFloat("Qd (fase dispersa)", &m_Qd_ul_min, 1.0f, 100.0f, "%.1f")) m_simulationParametersChanged = true;

    ImGui::Separator();
    ImGui::Text("Propiedades de Fluidos:");
    if (ImGui::SliderFloat("Visc. continua (mu_c) [cP]", &m_mu_c_cP, 0.1f, 100.0f, "%.2f", ImGuiSliderFlags_Logarithmic)) m_simulationParametersChanged = true;
    if (ImGui::SliderFloat("Visc. dispersa (mu_d) [cP]", &m_mu_d_cP, 0.1f, 100.0f, "%.2f", ImGuiSliderFlags_Logarithmic)) m_simulationParametersChanged = true;
    if (ImGui::SliderFloat("Tensión Interf. (gamma) [mN/m]", &m_gamma_mN_m, 0.1f, 70.0f, "%.2f")) m_simulationParametersChanged = true;

    ImGui::Separator();
    ImGui::Text("Geometría del Canal (micrómetros):");
    if (ImGui::SliderFloat("Ancho (w)", &m_w_um, 10.0f, 500.0f, "%.0f")) m_simulationParametersChanged = true;
    if (ImGui::SliderFloat("Altura (h)", &m_h_um, 10.0f, 500.0f, "%.0f")) m_simulationParametersChanged = true;

    ImGui::Separator();
    ImGui::Text("Parámetros de Visualización 3D:");
    if (ImGui::SliderFloat("Escala Mundial", &m_world_scale_factor, 0.001f, 0.1f, "%.3f", ImGuiSliderFlags_Logarithmic)) m_simulationParametersChanged = true;


    if (ImGui::Button("Resetear Parámetros"))
    {
        m_Qc_ul_min                   = 100.0f;
        m_Qd_ul_min                   = 10.0f;
        m_gamma_mN_m                  = 5.0f;
        m_mu_c_cP                     = 1.0f;
        m_mu_d_cP                     = 5.0f;
        m_w_um                        = 50.0f;
        m_h_um                        = 50.0f;
        m_world_scale_factor          = 0.025f;
        m_simulationParametersChanged = true;
    }

    ImGui::Separator();
    if (!m_Droplets.empty())
    {
        DimensionlessNumbers d_nums = CalculateDimensionlessGroups(
            static_cast<double>(m_Qc_ul_min), static_cast<double>(m_Qd_ul_min),
            static_cast<double>(m_mu_c_cP), static_cast<double>(m_mu_d_cP),
            static_cast<double>(m_gamma_mN_m), static_cast<double>(m_w_um), static_cast<double>(m_h_um));
        ImGui::Text("Números Adimensionales Calculados:");
        ImGui::Text("  Ca_c: %.4e", d_nums.Ca_c);
        ImGui::Text("  Phi:  %.3f", d_nums.Phi);
        ImGui::Text("  Lambda_v: %.3f", d_nums.lambda_v);
    }
    ImGui::End();
}


void Tutorial01_HelloTriangle::Initialize(const SampleInitInfo& InitInfo)
{
    SampleBase::Initialize(InitInfo);

    m_Qc_ul_min                   = 100.0f;
    m_Qd_ul_min                   = 10.0f;
    m_gamma_mN_m                  = 5.0f;
    m_mu_c_cP                     = 1.0f;
    m_mu_d_cP                     = 5.0f;
    m_w_um                        = 50.0f;
    m_h_um                        = 50.0f;
    m_world_scale_factor          = 0.025f;
    m_simulationParametersChanged = true;

    CreateUniformBuffers();

    GraphicsPipelineStateCreateInfo PSOCreateInfo;
    PSOCreateInfo.PSODesc.Name         = "Droplet PSO 3D Shaded";
    PSOCreateInfo.PSODesc.PipelineType = PIPELINE_TYPE_GRAPHICS;

    PSOCreateInfo.GraphicsPipeline.NumRenderTargets = 1;
    PSOCreateInfo.GraphicsPipeline.RTVFormats[0]    = m_pSwapChain->GetDesc().ColorBufferFormat;
    PSOCreateInfo.GraphicsPipeline.DSVFormat        = m_pSwapChain->GetDesc().DepthBufferFormat;

    PSOCreateInfo.GraphicsPipeline.PrimitiveTopology = PRIMITIVE_TOPOLOGY_TRIANGLE_LIST;

    PSOCreateInfo.GraphicsPipeline.RasterizerDesc.CullMode      = CULL_MODE_BACK;
    PSOCreateInfo.GraphicsPipeline.DepthStencilDesc.DepthEnable = True;

    ShaderCreateInfo ShaderCI;
    ShaderCI.SourceLanguage                  = SHADER_SOURCE_LANGUAGE_HLSL;
    ShaderCI.Desc.UseCombinedTextureSamplers = true;
    ShaderCI.CompileFlags                    = SHADER_COMPILE_FLAG_PACK_MATRIX_ROW_MAJOR;


    RefCntAutoPtr<IShader> pVS;
    {
        ShaderCI.Desc.ShaderType = SHADER_TYPE_VERTEX;
        ShaderCI.EntryPoint      = "main";
        ShaderCI.Desc.Name       = "Droplet Vertex Shader 3D Shaded";
        ShaderCI.Source          = VSSource_Droplets_3D_Shaded;
        m_pDevice->CreateShader(ShaderCI, &pVS);
    }

    RefCntAutoPtr<IShader> pPS;
    {
        ShaderCI.Desc.ShaderType = SHADER_TYPE_PIXEL;
        ShaderCI.EntryPoint      = "main";
        ShaderCI.Desc.Name       = "Droplet Pixel Shader 3D Shaded";
        ShaderCI.Source          = PSSource_Droplets_3D_Shaded;
        m_pDevice->CreateShader(ShaderCI, &pPS);
    }

    
    LayoutElement LayoutElems[] =
        {
            LayoutElement{0, 0, 3, VT_FLOAT32, False, INPUT_ELEMENT_FREQUENCY_PER_VERTEX}, // Pos
            LayoutElement{1, 0, 3, VT_FLOAT32, False, INPUT_ELEMENT_FREQUENCY_PER_VERTEX}, // Normal
            LayoutElement{2, 0, 3, VT_FLOAT32, False, INPUT_ELEMENT_FREQUENCY_PER_VERTEX}  // Color
        };
    PSOCreateInfo.GraphicsPipeline.InputLayout.LayoutElements = LayoutElems;
    PSOCreateInfo.GraphicsPipeline.InputLayout.NumElements    = _countof(LayoutElems);

    PSOCreateInfo.pVS = pVS;
    PSOCreateInfo.pPS = pPS;

    PSOCreateInfo.PSODesc.ResourceLayout.DefaultVariableType = SHADER_RESOURCE_VARIABLE_TYPE_STATIC;

    m_pDevice->CreateGraphicsPipelineState(PSOCreateInfo, &m_pPSO);

    if (m_pPSO && m_pVSConstants)
    {
        m_pPSO->GetStaticVariableByName(SHADER_TYPE_VERTEX, "VSConstants")->Set(m_pVSConstants);
    }
    else
    {
        LOG_ERROR_MESSAGE("Failed to create PSO or m_pVSConstants is null, cannot bind VSConstants to PSO.");
    }

    if (m_pPSO)
    {
        m_pPSO->CreateShaderResourceBinding(&m_pSRB, true);
        if (!m_pSRB)
        {
            LOG_ERROR_MESSAGE("Failed to create SRB for Droplet PSO 3D Shaded.");
        }
    }
    else
    {
        LOG_ERROR_MESSAGE("Failed to create PSO for 3D Shaded droplets, SRB cannot be created.");
    }
}

void Tutorial01_HelloTriangle::Render()
{
    const float   ClearColor[] = {0.1f, 0.1f, 0.2f, 1.0f};
    ITextureView* pRTV         = m_pSwapChain->GetCurrentBackBufferRTV();
    ITextureView* pDSV         = m_pSwapChain->GetDepthBufferDSV();
    m_pImmediateContext->ClearRenderTarget(pRTV, ClearColor, RESOURCE_STATE_TRANSITION_MODE_TRANSITION);
    m_pImmediateContext->ClearDepthStencil(pDSV, CLEAR_DEPTH_FLAG, 1.f, 0, RESOURCE_STATE_TRANSITION_MODE_TRANSITION);

    if (!m_pPSO || !m_pDropletVertexBuffer || m_TotalVertexCount == 0 || !m_pVSConstants || !m_pSRB)
    {
        return;
    }

    {
        MapHelper<float4x4> CBConstants(m_pImmediateContext, m_pVSConstants, MAP_WRITE, MAP_FLAG_DISCARD);
        *CBConstants = m_ViewMatrix * m_ProjMatrix;
    }

    m_pImmediateContext->SetPipelineState(m_pPSO);
    m_pImmediateContext->CommitShaderResources(m_pSRB, RESOURCE_STATE_TRANSITION_MODE_TRANSITION);


    IBuffer* pBuffs[]  = {m_pDropletVertexBuffer};
    Uint64   Offsets[] = {0};
    m_pImmediateContext->SetVertexBuffers(0,
                                          1,
                                          pBuffs,
                                          Offsets,
                                          RESOURCE_STATE_TRANSITION_MODE_TRANSITION,
                                          SET_VERTEX_BUFFERS_FLAG_RESET);

    DrawAttribs drawAttrs;
    drawAttrs.NumVertices = m_TotalVertexCount;
    drawAttrs.Flags       = DRAW_FLAG_VERIFY_ALL;
    m_pImmediateContext->Draw(drawAttrs);
}

void Tutorial01_HelloTriangle::Update(double CurrTime, double ElapsedTime)
{
    SampleBase::Update(CurrTime, ElapsedTime);

    UpdateUI();

    if (m_simulationParametersChanged)
    {
        PrepareDropletData();
        CreateDropletGeometryAndBuffers();
        m_simulationParametersChanged = false;
    }

    bool needs_rebuild_for_animation = false;
    if (m_Droplets.empty())
    {
        if (m_pDropletVertexBuffer)
        {
            m_DropletVertices.clear();
            m_TotalVertexCount = 0;
            m_pDropletVertexBuffer.Release();
        }
        return;
    }

    for (auto& droplet : m_Droplets)
    {
        droplet.CurrentPositionWorld.x += droplet.VelocityX_World * static_cast<float>(ElapsedTime);

        float right_boundary_world = 2.5f;

        if (droplet.CurrentPositionWorld.x - droplet.RadiusWorld > right_boundary_world)
        {
            droplet.CurrentPositionWorld.x   = m_next_wrap_around_spawn_x_world;
            m_next_wrap_around_spawn_x_world = droplet.CurrentPositionWorld.x - droplet.RadiusWorld - m_spacing_world - droplet.RadiusWorld;
        }
        needs_rebuild_for_animation = true;
    }

    if (needs_rebuild_for_animation)
    {
        if (!m_simulationParametersChanged)
        {
            CreateDropletGeometryAndBuffers();
        }
    }

    float camera_distance    = 5.0f;
    float camera_height      = 0.75f; // Elevar un poco más la cámara
    float camera_side_offset = 1.0f;  // Mover un poco a un lado

 
    m_ViewMatrix = float4x4::Translation(-camera_side_offset, -camera_height, camera_distance);

    float NearClip = 0.01f;
    float FarClip  = 200.0f;
    m_ProjMatrix   = GetAdjustedProjectionMatrix(PI_F / 4.0f, NearClip, FarClip);
}

} // namespace Diligent
