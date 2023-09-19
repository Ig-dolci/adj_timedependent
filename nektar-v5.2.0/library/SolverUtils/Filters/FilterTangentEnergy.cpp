///////////////////////////////////////////////////////////////////////////////
//
// File FilterTangentEnergy.cpp
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
// Description: Output kinetic energy and enstrophy.
//
///////////////////////////////////////////////////////////////////////////////

#include <iomanip>

#include <boost/core/ignore_unused.hpp>

#include <SolverUtils/Filters/FilterTangentEnergy.h>
#include <SolverUtils/Filters/FilterInterfaces.hpp>

using namespace std;

namespace Nektar
{
namespace SolverUtils
{
std::string FilterTangentEnergy::className =
    SolverUtils::GetFilterFactory().RegisterCreatorFunction(
        "TangentEnergy", FilterTangentEnergy::create);

FilterTangentEnergy::FilterTangentEnergy(const LibUtilities::SessionReaderSharedPtr &pSession,
                           const std::weak_ptr<EquationSystem> &pEquation,
                           const ParamMap &pParams)
    : Filter(pSession, pEquation), m_index(-1), m_homogeneous(false), m_planes()
{
    std::string outName;

    // OutputFile
    auto it = pParams.find("OutputFile");
    if (it == pParams.end())
    {
        outName = m_session->GetSessionName();
    }
    else
    {
        ASSERTL0(it->second.length() > 0, "Missing parameter 'OutputFile'.");
        outName = it->second;
    }
    outName += ".eny";

    m_comm = pSession->GetComm();
    if (m_comm->GetRank() == 0)
    {
        m_outFile.open(outName.c_str());
        ASSERTL0(m_outFile.good(), "Unable to open: '" + outName + "'");
        m_outFile.setf(ios::scientific, ios::floatfield);
        m_outFile << "# Time                Kinetic energy        " << endl
                  << "# ---------------------------------------------"
                  << "--------------" << endl;
    }
    pSession->LoadParameter("LZ", m_homogeneousLength, 0.0);

    // OutputFrequency
    it = pParams.find("OutputFrequency");
    ASSERTL0(it != pParams.end(), "Missing parameter 'OutputFrequency'.");
    LibUtilities::Equation equ(m_session->GetInterpreter(), it->second);
    m_outputFrequency = round(equ.Evaluate());
}

FilterTangentEnergy::~FilterTangentEnergy()
{
}

void FilterTangentEnergy::v_Initialise(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    m_index = -1;
    MultiRegions::ExpListSharedPtr areaField;

    ASSERTL0(pFields[0]->GetExpType() != MultiRegions::e1D,
             "1D expansion not supported for energy filter");

    //ASSERTL0(pFields[0]->GetExpType() != MultiRegions::e2D,
    //         "2D expansion not supported for energy filter");

    ASSERTL0(pFields[0]->GetExpType() != MultiRegions::e3DH2D,
             "Homogeneous 2D expansion not supported for energy filter");

    if (pFields[0]->GetExpType() == MultiRegions::e3DH1D)
    {
        m_homogeneous = true;
    }

    // Calculate area/volume of domain.
    if (m_homogeneous)
    {
        m_planes  = pFields[0]->GetZIDs();
        areaField = pFields[0]->GetPlane(0);
    }
    else
    {
        areaField = pFields[0];
    }

    Array<OneD, NekDouble> inarray(areaField->GetNpoints(), 1.0);
    m_area = areaField->Integral(inarray);

    if (m_homogeneous)
    {
        m_area *= m_homogeneousLength;
    }

    // Output values at initial time.
    v_Update(pFields, time);
}

void FilterTangentEnergy::v_Update(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    int i, nPoints = pFields[0]->GetNpoints();

    m_index++;

    if (m_index % m_outputFrequency > 0)
    {
        return;
    }


    // Lock equation system pointer
    auto equ = m_equ.lock();
    ASSERTL0(equ, "Weak pointer expired");

    auto fluidEqu = std::dynamic_pointer_cast<FluidInterface>(equ);
    ASSERTL0(fluidEqu, "Energy filter is incompatible with this solver.");

    // Store physical values in an array
    Array<OneD, Array<OneD, NekDouble>> physfields(pFields.size());
    for (i = 0; i < pFields.size(); ++i)
    {
        physfields[i] = pFields[i]->GetPhys();
    }

    // Calculate kinetic energy.
    int dim = 3;
    if (pFields[0]->GetExpType() == MultiRegions::e2D)
    {
        dim = 2;
    }
    NekDouble E = 0.0;
    Array<OneD, NekDouble> tmp(nPoints, 0.0);
    Array<OneD, NekDouble> density;
    Array<OneD, Array<OneD, NekDouble>> u(dim);
    Array<OneD, Array<OneD, NekDouble>> v(dim);
    for (i = 0; i < dim; ++i)
    {
        u[i] = Array<OneD, NekDouble>(nPoints);
        v[i] = Array<OneD, NekDouble>(nPoints);
    }
    fluidEqu->GetVelocity(physfields, v);
    fluidEqu->GetBaseFlow2(u);

    for (i = 0; i < dim; ++i)
    {
        if (m_homogeneous && pFields[i]->GetWaveSpace())
        {
            pFields[i]->HomogeneousBwdTrans(u[i], u[i]);
            pFields[i]->HomogeneousBwdTrans(v[i], v[i]);
        }

        Vmath::Vvtvp(nPoints, u[i], 1, v[i], 1, tmp, 1, tmp, 1);
    }

    if (!fluidEqu->HasConstantDensity())
    {
        density = Array<OneD, NekDouble>(nPoints);
        fluidEqu->GetDensity(physfields, density);
        Vmath::Vmul(nPoints, density, 1, tmp, 1, tmp, 1);
    }

    if (m_homogeneous)
    {
        Array<OneD, NekDouble> tmp2(nPoints, 0.0);
        pFields[0]->HomogeneousFwdTrans(tmp, tmp2);
        E = pFields[0]->GetPlane(0)->Integral(tmp2) * m_homogeneousLength;
    }
    else
    {
        E = pFields[0]->Integral(tmp);
    }


    if (m_comm->GetRank() == 0)
    {
        m_outFile << setw(17) << setprecision(8) << time << setw(22)
                  << setprecision(11) << E << endl;
    }
}

void FilterTangentEnergy::v_Finalise(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    boost::ignore_unused(pFields, time);
    m_outFile.close();
}

bool FilterTangentEnergy::v_IsTimeDependent()
{
    return true;
}

} // namespace SolverUtils
} // namespace Nektar
