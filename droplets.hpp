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

#pragma once

#include "SampleBase.hpp"
#include "BasicMath.hpp" // Para float2, float3, float4x4

namespace Diligent
{
// Estructura para definir los atributos de un vértice de gota 
struct DropletVertex
{
    float3 Pos;    // Posición 3D del vértice
    float3 Normal; // Normal del vértice para iluminación
    float3 Color;  // Color base del vértice
};

// Estructura para almacenar la información de cada gota
struct DropletInfo
{
    // Propiedades de la simulación empírica
    float  DiameterModel;
    float  AspectRatio;
    float3 Color;
    int    Regime;

    // Propiedades para la visualización y animación 3D
    float3 CurrentPositionWorld;
    float  RadiusWorld;
    float  VelocityX_World;
    float  Concentration;
};

class Tutorial01_HelloTriangle final : public SampleBase
{
public:
    virtual void Initialize(const SampleInitInfo& InitInfo) override final;

    virtual void        Render() override final;
    virtual void        Update(double CurrTime, double ElapsedTime) override final;
    virtual const Char* GetSampleName() const override final { return "Tutorial01: 3D Shaded Droplets"; }

private:
    RefCntAutoPtr<IPipelineState>         m_pPSO;
    RefCntAutoPtr<IBuffer>                m_pDropletVertexBuffer;
    RefCntAutoPtr<IShaderResourceBinding> m_pSRB;
    std::vector<DropletInfo>              m_Droplets;
    std::vector<DropletVertex>            m_DropletVertices;
    Uint32                                m_TotalVertexCount = 0;

    // --- Parámetros de visualización y animación ---
    float m_world_scale_factor             = 0.025f; // Ajustado para mejor visibilidad inicial
    float m_spacing_world                  = 0.0f;
    float m_next_wrap_around_spawn_x_world = -2.0f;

    // --- Parámetros de simulación controlados por UI ---
    float m_Qc_ul_min  = 100.0f;
    float m_Qd_ul_min  = 10.0f;
    float m_gamma_mN_m = 5.0f;
    float m_mu_c_cP    = 1.0f;
    float m_mu_d_cP    = 5.0f;
    float m_w_um       = 50.0f;
    float m_h_um       = 50.0f;

    bool m_simulationParametersChanged = true;

    // --- Matrices para 3D ---
    RefCntAutoPtr<IBuffer> m_pVSConstants;
    float4x4               m_ViewMatrix;
    float4x4               m_ProjMatrix;
   
    // float4x4 m_NormalMatrix;


    // --- Funciones ---
    void PrepareDropletData();
    void CreateDropletGeometryAndBuffers();
    void UpdateUI();
    void CreateUniformBuffers();

    // --- Lógica de predicción empírica ---
    struct DimensionlessNumbers
    {
        double Ca_c;
        double Phi;
        double lambda_v;
        double w_SI;
        double U_c;
    };

    DimensionlessNumbers CalculateDimensionlessGroups(double Qc, double Qd, double mu_c, double mu_d, double gamma, double w, double h = -1.0);
    void                 PredictDropletMetrics(const DimensionlessNumbers& dimensionless, double w_um);
};

} // namespace Diligent
