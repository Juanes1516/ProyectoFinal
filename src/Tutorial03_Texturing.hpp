/*
 *  Copyright 2019-2022 Diligent Graphics LLC
 *  Copyright 2015-2019 Egor Yusov
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *  In no event and under no legal theory, whether in tort (including negligence),
 *  contract, or otherwise, unless required by applicable law (such as deliberate
 *  and grossly negligent acts) or agreed to in writing, shall any Contributor be
 *  liable for any damages, including any direct, indirect, special, incidental,
 *  or consequential damages of any character arising as a result of this License or
 *  out of the use or inability to use the software (including but not limited to damages
 *  for loss of goodwill, work stoppage, computer failure or malfunction, or any and
 *  all other commercial damages or losses), even if such Contributor has been advised
 *  of the possibility of such damages.
 */

#pragma once

#include "SampleBase.hpp"
#include "BasicMath.hpp"

namespace Diligent
{

class Tutorial03_Texturing final : public SampleBase
{
public:
    virtual void Initialize(const SampleInitInfo& InitInfo) override final;

    virtual void Render() override final;
    virtual void Update(double CurrTime, double ElapsedTime) override final;

    virtual const Char* GetSampleName() const override final { return "Tutorial03: Texturing"; }

private:
    void CreatePipelineState();
    void CreateVertexBuffer();
    void CreateIndexBuffer();
    void LoadTexture();

    RefCntAutoPtr<IPipelineState>         m_pPSO;
    RefCntAutoPtr<IBuffer>                m_CubeVertexBuffer;
    RefCntAutoPtr<IBuffer>                m_CubeIndexBuffer;
    RefCntAutoPtr<IBuffer>                m_VSConstants;
    RefCntAutoPtr<ITextureView>           m_TextureSRV;
    RefCntAutoPtr<IShaderResourceBinding> m_SRB;
    float4x4                              m_WorldViewProjMatrix;

    void                   CreateSphereBuffers();
    RefCntAutoPtr<IBuffer> m_SphereVertexBuffer;
    RefCntAutoPtr<IBuffer> m_SphereIndexBuffer;
    Uint32                 m_NumSphereIndices = 0;

    // Número de partículas (ajústalo a tu gusto)
    static const Uint32 kNumParticles = 1024;

    // Buffers y vistas
    RefCntAutoPtr<IBuffer>     m_pParticleAttribsBuf;
    RefCntAutoPtr<IBufferView> m_pParticleAttribsSRV, m_pParticleAttribsUAV;
    RefCntAutoPtr<IBuffer>     m_pParticleListHeadsBuf, m_pParticleListsBuf;
    RefCntAutoPtr<IBufferView> m_pParticleListHeadsUAV, m_pParticleListsUAV;

    // Pipeline States (compute + render)
    RefCntAutoPtr<IPipelineState> m_pResetListsPSO, m_pMovePSO, m_pCollidePosPSO, m_pCollideVelPSO;
    RefCntAutoPtr<IPipelineState> m_pRenderParticlesPSO;
    RefCntAutoPtr<IBuffer>        m_pWireframeCB;

    RefCntAutoPtr<IPipelineState>         m_pWireframePSO;
    RefCntAutoPtr<IShaderResourceBinding> m_pWireframeSRB;
    float                                 m_WireframeIntensity = 0.0f;

    // Shader resource bindings
    RefCntAutoPtr<IShaderResourceBinding> m_pResetListsSRB, m_pMoveSRB, m_pCollideSRB, m_pRenderParticlesSRB;

    // Constantes globales (grid size, deltaTime…)
    struct GlobalConstants
    {
        int2   i2ParticleGridSize;
        Uint32 uiNumParticles;
        float2 f2Scale;
        float  fDeltaTime;
    } m_Constants;
    RefCntAutoPtr<IBuffer> m_pConstBuf; // para Compute
};

} // namespace Diligent
