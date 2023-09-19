///////////////////////////////////////////////////////////////////////////////
//
// File: ForcingProgrammatic.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: Programmatic forcing
//
///////////////////////////////////////////////////////////////////////////////

#include <boost/core/ignore_unused.hpp>

#include <SolverUtils/Forcing/ForcingProgrammatic.h>
#include <SolverUtils/Filters/FilterInterfaces.hpp>

namespace Nektar
{
namespace SolverUtils
{

std::string ForcingProgrammatic::className =
    GetForcingFactory().RegisterCreatorFunction(
        "Programmatic", ForcingProgrammatic::create, "Programmatic Forcing");

ForcingProgrammatic::ForcingProgrammatic(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const std::weak_ptr<EquationSystem> &pEquation)
    : Forcing(pSession, pEquation)
{
    pSession->LoadParameter("LZ", m_homogeneousLength, 0.0);
    pSession->LoadParameter("P0", m_p0, 0.001);
    pSession->LoadParameter("beta", m_beta, 1.0);
}

Array<OneD, Array<OneD, NekDouble>> &ForcingProgrammatic::UpdateForces()
{
    return m_Forcing;
}

void ForcingProgrammatic::v_InitObject(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
    const unsigned int &pNumForcingFields, const TiXmlElement *pForce)
{
    boost::ignore_unused(pForce);

    m_NumVariable = pNumForcingFields;
    int nq        = pFields[0]->GetTotPoints();

    m_Forcing = Array<OneD, Array<OneD, NekDouble>>(m_NumVariable);
    for (int i = 0; i < m_NumVariable; ++i)
    {
        m_Forcing[i] = Array<OneD, NekDouble>(nq, 0.0);
    }

}

void ForcingProgrammatic::v_Apply(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble &time)
{
    boost::ignore_unused(fields, inarray, time);

    int i, nPoints = fields[0]->GetNpoints();

    auto equ = m_equ.lock();
    ASSERTL0(equ, "Weak pointer to the equation system is expired");
    auto FluidEq = std::dynamic_pointer_cast<FluidInterface>(equ);
    //FluidEq->GetVelocity(m_Forcing);

     // Store physical values in an array
    Array<OneD, Array<OneD, NekDouble>> physfields(fields.size());
    for (i = 0; i < fields.size(); ++i)
    {
        physfields[i] = fields[i]->GetPhys();
    }

    // Calculate kinetic energy.
    int dim = 3;
    if (fields[0]->GetExpType() == MultiRegions::e2D)
    {
        dim = 2;
    }
    NekDouble E = 0.0;
    Array<OneD, NekDouble> tmp(nPoints, 0.0);
    Array<OneD, NekDouble> density;
    Array<OneD, Array<OneD, NekDouble>> u(dim);
    for (i = 0; i < dim; ++i)
    {
        u[i] = Array<OneD, NekDouble>(nPoints);
    }
    FluidEq->GetVelocity(physfields, u);

    bool m_homogeneous = false;
    if (fields[0]->GetExpType() == MultiRegions::e3DH1D)
    {
        m_homogeneous = true;
    }

    for (i = 0; i < dim; ++i)
    {
        if (m_homogeneous && fields[i]->GetWaveSpace())
        {
            fields[i]->HomogeneousBwdTrans(u[i], u[i]);
        }

        Vmath::Vvtvp(nPoints, u[i], 1, u[i], 1, tmp, 1, tmp, 1);
    }

    if (!FluidEq->HasConstantDensity())
    {
        density = Array<OneD, NekDouble>(nPoints);
        FluidEq->GetDensity(physfields, density);
        Vmath::Vmul(nPoints, density, 1, tmp, 1, tmp, 1);
    }

    if (m_homogeneous)
    {
        Array<OneD, NekDouble> tmp2(nPoints, 0.0);
        fields[0]->HomogeneousFwdTrans(tmp, tmp2);
        E = fields[0]->GetPlane(0)->Integral(tmp2) * m_homogeneousLength;
    }
    else
    {
        E = fields[0]->Integral(tmp);
    }

    E /= 2.0;
    
    for (int i = 0; i < m_NumVariable; i++)
    {
        Vmath::Svtvp(outarray[i].size(), m_p0*m_beta/E, u[i], 1, m_Forcing[i], 1, m_Forcing[i], 1);
    }


    for (int i = 0; i < m_NumVariable; i++)
    {
        Vmath::Vadd(outarray[i].size(), outarray[i], 1, m_Forcing[i], 1,
                    outarray[i], 1);
    }
}

} // namespace SolverUtils
} // namespace Nektar
