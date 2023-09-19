///////////////////////////////////////////////////////////////////////////////
//
// File FilterAdjointGrad.cpp
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
// Description: Output values of aerodynamic forces during time-stepping.
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <iomanip>
#include <boost/core/ignore_unused.hpp>
#include <boost/algorithm/string.hpp>
#include <LocalRegions/Expansion1D.h>
#include <LocalRegions/Expansion2D.h>
#include <LocalRegions/Expansion3D.h>
#include <MultiRegions/ExpList.h>     
#include <MultiRegions/ExpList3DHomogeneous1D.h>
#include <SolverUtils/Filters/FilterAdjointGrad.h>
#include <SolverUtils/Filters/FilterInterfaces.hpp>
#include <LibUtilities/BasicUtils/ParseUtils.h>

using namespace std;

namespace Nektar
{
namespace SolverUtils
{
std::string FilterAdjointGrad::className =
        GetFilterFactory().RegisterCreatorFunction(
                "AdjointGrad", FilterAdjointGrad::create);

/**
 *
 */
FilterAdjointGrad::FilterAdjointGrad(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const std::weak_ptr<EquationSystem>      &pEquation,
    const ParamMap &pParams) :
    Filter(pSession, pEquation)
{
    // OutputFile
    auto it = pParams.find("OutputFile");
    if (it == pParams.end())
    {
        m_outputFile = m_session->GetSessionName();
    }
    else
    {
        ASSERTL0(it->second.length() > 0, "Missing parameter 'OutputFile'.");
        m_outputFile = it->second;
    }
    if (!(m_outputFile.length() >= 4
          && m_outputFile.substr(m_outputFile.length() - 4) == ".fce"))
    {
        m_outputFile += ".fce";
    }

    // OutputFrequency
    it = pParams.find("OutputFrequency");
    if (it == pParams.end())
    {
        m_outputFrequency = 1;
    }
    else
    {
        LibUtilities::Equation equ(
            m_session->GetInterpreter(), it->second);
        m_outputFrequency = round(equ.Evaluate());
    }

    // Time after which we need to calculate the forces
    it = pParams.find("StartTime");
    if (it == pParams.end())
    {
        m_startTime = 0;
    }
    else
    {
        LibUtilities::Equation equ(
            m_session->GetInterpreter(), it->second);
        m_startTime = equ.Evaluate();
    }

    // For 3DH1D, OutputAllPlanes or just average forces? just average forces

    m_session->MatchSolverInfo("Homogeneous", "1D",
                               m_isHomogeneous1D, false);
    m_outputAllPlanes = false;

    // Boundary (to calculate forces on)
    it = pParams.find("Boundary");
    ASSERTL0(it != pParams.end(),   "Missing parameter 'Boundary'.");
    ASSERTL0(it->second.length() > 0, "Empty parameter 'Boundary'.");
    m_BoundaryString = it->second;

    //
    // Directions (to project forces)
    //

    // Allocate m_directions
    m_directions = Array<OneD, Array<OneD, NekDouble> > (3);
    //Initialise directions to default values (ex, ey, ez)
    for (int i = 0; i < 3; ++i)
    {
        m_directions[i] = Array<OneD, NekDouble>(3, 0.0);
        m_directions[i][i] = 1.0;
    }
    std::stringstream       directionStream;
    std::string             directionString;
    //Override with input from xml file (if defined)
    for (int i = 0; i < 3; ++i)
    {
        std::stringstream tmp;
        tmp << i+1;
        std::string dir = "Direction" + tmp.str();
        it = pParams.find(dir);
        if ( it != pParams.end() )
        {
            ASSERTL0(!(it->second.empty()),
                     "Missing parameter '"+dir+"'.");
            directionStream.str(it->second);
            // Guarantee the stream is in its start position
            //      before extracting
            directionStream.clear();
            // normalisation factor
            NekDouble norm = 0.0;
            for (int j = 0; j < 3; j++)
            {
                directionStream >> directionString;
                if (!directionString.empty())
                {
                    LibUtilities::Equation equ(
                        m_session->GetInterpreter(), directionString);
                    m_directions[i][j] = equ.Evaluate();
                    norm += m_directions[i][j]*m_directions[i][j];
                }
            }
            // Normalise direction
            for( int j = 0; j < 3; j++)
            {
                m_directions[i][j] /= sqrt(norm);
            }
        }
    }


    // Case definition
    it = pParams.find("Case");
    ASSERTL0(it != pParams.end(),   "Missing parameter 'Case'.");
    ASSERTL0(it->second.length() > 0, "Empty parameter 'Case'.");
    m_case = it->second;

}


/**
 *
 */
FilterAdjointGrad::~FilterAdjointGrad()
{

}

/**
 *
 */
void FilterAdjointGrad::v_Initialise(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    // Load mapping
    m_mapping = GlobalMapping::Mapping::Load(m_session, pFields);


    // Parse the boundary regions into a list.
    std::string::size_type FirstInd =
                            m_BoundaryString.find_first_of('[') + 1;
    std::string::size_type LastInd =
                            m_BoundaryString.find_last_of(']') - 1;

    ASSERTL0(FirstInd <= LastInd,
            (std::string("Error reading boundary region definition:") +
             m_BoundaryString).c_str());

    std::string IndString =
            m_BoundaryString.substr(FirstInd, LastInd - FirstInd + 1);
    bool parseGood = ParseUtils::GenerateSeqVector(IndString,
                                                   m_boundaryRegionsIdList);
    ASSERTL0(parseGood && !m_boundaryRegionsIdList.empty(),
             (std::string("Unable to read boundary regions index "
              "range for FilterAdjointGrad: ") + IndString).c_str());

    // determine what boundary regions need to be considered
    int cnt;
    unsigned int numBoundaryRegions =
                        pFields[0]->GetBndConditions().size();
    m_boundaryRegionIsInList.insert(m_boundaryRegionIsInList.end(),
                                    numBoundaryRegions, 0);

    SpatialDomains::BoundaryConditions bcs(m_session,
                                            pFields[0]->GetGraph());
    const SpatialDomains::BoundaryRegionCollection &bregions =
                                            bcs.GetBoundaryRegions();

    cnt = 0;
    for (auto &it : bregions)
    {
        if ( std::find(m_boundaryRegionsIdList.begin(),
                       m_boundaryRegionsIdList.end(), it.first) !=
                m_boundaryRegionsIdList.end() )
        {
            m_boundaryRegionIsInList[cnt] = 1;
        }
        cnt++;
    }


    // Create map for element and edge/face of each boundary expansion
    if(m_isHomogeneous1D)
    {
        pFields[0]->GetPlane(0)->GetBoundaryToElmtMap
                                        (m_BCtoElmtID,m_BCtoTraceID);
    }
    else
    {
        pFields[0]->GetBoundaryToElmtMap(m_BCtoElmtID,m_BCtoTraceID);
    }

    // Define number of planes  to calculate the forces
    //     in the Homogeneous direction ( if m_outputAllPlanes is false,
    //      consider only first plane in wave space)
    // If flow has no Homogeneous direction, use 1 to make code general
    if(m_isHomogeneous1D &&(m_outputAllPlanes || m_mapping->IsDefined()))
    {
        m_nPlanes = pFields[0]->GetHomogeneousBasis()->
                                            GetZ().size();
    }
    else
    {
        m_nPlanes = 1;
    }


    // Create map for Planes ID for Homogeneous direction
    //    If flow has no Homogeneous direction, create trivial map
    int j;
    m_planesID = Array<OneD, int> (m_nPlanes,-1);
    if(m_isHomogeneous1D)
    {
        Array<OneD, const unsigned int> IDs = pFields[0]->GetZIDs();
        //Loop through all planes
        for(cnt = 0; cnt < m_nPlanes; cnt++)
        {
            for(j = 0; j < IDs.size(); ++j)
            {
                //Look for current plane ID in this process
                if(IDs[j] == cnt)
                {
                    break;
                }
            }
            // Assign ID to planesID
            // If plane is not found in this process, leave it with -1
            if(j != IDs.size())
            {
                m_planesID[cnt] = j;
            }
        }
    }
    else
    {
        m_planesID[0] = 0;
    }


    // Case definition
    GetCaseInfo();
    m_adjgrad = Array<OneD, NekDouble>(m_nParam, 0.0);
    m_adjgrad2 = Array<OneD, NekDouble>(m_nParam, 0.0);
    m_aux = 1;
    m_Ftp = Array<OneD, NekDouble>(m_nParam, 0.0);
    m_Ftpp = Array<OneD, NekDouble>(m_nParam, 0.0);
    

    LibUtilities::CommSharedPtr vComm = pFields[0]->GetComm();

    // Write header
    //int expdim = pFields[0]->GetGraph()->GetMeshDimension();
    if (vComm->GetRank() == 0)
    {
        // Open output stream
        bool adaptive;
        m_session->MatchSolverInfo("Driver", "Adaptive",
                                    adaptive, false);
        if (adaptive)
        {
            m_outputStream.open(m_outputFile.c_str(), ofstream::app);
        }
        else
        {
            m_outputStream.open(m_outputFile.c_str());
        }
        m_outputStream << "# Adjoint gradient" << endl;
        m_outputStream << "# Case: " << m_case << endl;
        m_outputStream << "# Boundary regions: " << IndString.c_str() << endl;
        m_outputStream << "#";
        m_outputStream.width(7);
        m_outputStream << "Time";
        for( int i = 1; i <= m_nParam; i++ )
        {
            //m_outputStream.width(8);
            //m_outputStream <<  "grad" << i << "-press";
            //m_outputStream.width(9);
            //m_outputStream <<  "grad" << i << "-visc";
            m_outputStream.width(11);
            m_outputStream <<  "grad" << i;
            m_outputStream.width(16);
            m_outputStream <<  "1st_int_grad" << i;
            m_outputStream.width(15);
            m_outputStream <<  "2nd_int_grad" << i;
        }
        if (m_session->DefinesSolverInfo("HomoStrip"))
        {
            ASSERTL0(m_outputAllPlanes==false,
                    "Output forces on all planes not compatible with HomoStrips");
            m_outputStream.width(10);
            m_outputStream << "Strip";
        }

        m_outputStream << endl;
    }

    m_lastTime = -1;
    m_index = 0;
    
    v_Update(pFields, time);
}

/**
 *
 */
void FilterAdjointGrad::v_Update(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    // Only output every m_outputFrequency.
    if ((m_index++) % m_outputFrequency  || (time < m_startTime))
    {
        return;
    }
    // Calculate the forces
    for (int i = 0; i < m_nParam; i++) 
    {
        m_Ftpp[i] = m_Ftp[i];
        m_Ftp[i] = m_Ft[i];

    }
    CalculateForces(pFields, time);

    for (int i = 0; i < m_nParam; i++) 
    {
        m_adjgrad[i] += (m_Ft[i]+m_Ftp[i])/2 * m_session->GetParameter("TimeStep") * m_outputFrequency;
        if (m_aux > 0)
        {
            m_adjgrad2[i] += (m_Ft[i]+4*m_Ftp[i]+m_Ftpp[i])/3 * m_session->GetParameter("TimeStep") * m_outputFrequency;
        }

    }
    m_aux = m_aux*(-1);



    // Communicators to exchange results
    LibUtilities::CommSharedPtr vComm = pFields[0]->GetComm();

    //Write Results
    if (vComm->GetRank() == 0)
    {

        // Output average (or total) force
        m_outputStream.width(8);
        m_outputStream << setprecision(6) << time;
        for( int i = 0; i < m_nParam; i++)
        {
            //m_outputStream.width(15);
            //m_outputStream << setprecision(8) << m_Fp[i];
            //m_outputStream.width(15);
            //m_outputStream << setprecision(8) << m_Fv[i];
            m_outputStream.width(15);
            m_outputStream << setprecision(8) << m_Ft[i];
            m_outputStream.width(15);
            m_outputStream << setprecision(8) << m_adjgrad[i];
            m_outputStream.width(15);
            if(m_aux < 0)
            {
                m_outputStream << setprecision(8) << m_adjgrad2[i];
            }
            else
            {
                m_outputStream << setprecision(8) << 0.0;
            }
        }

        if( m_session->DefinesSolverInfo("HomoStrip"))
        {
            // The result we already wrote is for strip 0
            m_outputStream.width(10);
            m_outputStream << 0;
            m_outputStream << endl;

            // Now get result from other strips
            int nstrips;
            m_session->LoadParameter("Strip_Z", nstrips);
            for(int i = 1; i<nstrips; i++)
            {
                //vComm->GetColumnComm()->Recv(i, m_Fp);
                //vComm->GetColumnComm()->Recv(i, m_Fv);
                vComm->GetColumnComm()->Recv(i, m_Ft);

                m_outputStream.width(8);
                m_outputStream << setprecision(6) << time;
                for( int j = 0; j < m_nParam; j++)
                {
                    //m_outputStream.width(15);
                    //m_outputStream << setprecision(8) << m_Fp[j];
                    //m_outputStream.width(15);
                    //m_outputStream << setprecision(8) << m_Fv[j];
                    m_outputStream.width(15);
                    m_outputStream << setprecision(8) << m_Ft[j];
                }
            m_outputStream.width(10);
            m_outputStream << i;
            m_outputStream << endl;
            }
        }
        else
        {
            m_outputStream << endl;
        }
    }
    else
    {
        // In the homostrips case, we need to send result to root
        if (m_session->DefinesSolverInfo("HomoStrip") &&
                (vComm->GetRowComm()->GetRank() == 0) )
        {
                //vComm->GetColumnComm()->Send(0, m_Fp);
                //vComm->GetColumnComm()->Send(0, m_Fv);
                vComm->GetColumnComm()->Send(0, m_Ft);
        }
    }

}


/**
 *
 */
void FilterAdjointGrad::v_Finalise(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    boost::ignore_unused(time);

    if (pFields[0]->GetComm()->GetRank() == 0)
    {
        m_outputStream.close();
    }
}


/**
 *
 */
bool FilterAdjointGrad::v_IsTimeDependent()
{
    return true;
}

/**
 *     This function calculates the forces
 */
void FilterAdjointGrad::CalculateForces(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time)
{
    // Store time so we can check if result is up-to-date
    m_lastTime = time;

    std::cout << "bla1" << endl;

    if (m_mapping->IsDefined())
    {
        CalculateForcesMapping( pFields, time);
        return;
    }

    std::cout << "bla2" << endl;

    // Lock equation system weak pointer
    auto equ = m_equ.lock();
    ASSERTL0(equ, "Weak pointer expired");

    auto fluidEqu = std::dynamic_pointer_cast<FluidInterface>(equ);
    ASSERTL0(fluidEqu, "Aero forces filter is incompatible with this solver.");

    int i, j, k, n, l, cnt, elmtid, nq, offset, boundary, plane;
    // Get number of quadrature points and dimensions
    int physTot = pFields[0]->GetNpoints();
    int expdim  = pFields[0]->GetGraph()->GetMeshDimension();
    int nVel    = expdim;
    if (m_isHomogeneous1D)
    {
        nVel = nVel + 1;
    }

    LocalRegions::ExpansionSharedPtr elmt;

    // Fields used to calculate forces (a single plane for 3DH1D)
    Array<OneD, MultiRegions::ExpListSharedPtr> fields(pFields.size());

    // Arrays of variables in field
    Array<OneD, Array<OneD, NekDouble>> physfields(pFields.size());
    Array<OneD, Array<OneD, NekDouble>> velocity(nVel);
    Array<OneD, NekDouble> pressure;

    // Arrays of variables in the element
    Array<OneD, Array<OneD, NekDouble>> velElmt(expdim);
    Array<OneD, NekDouble> pElmt(physTot);

    // Velocity gradient
    Array<OneD, Array<OneD, NekDouble>> grad(expdim * expdim);
    Array<OneD, NekDouble> div;

    // Values at the boundary
    Array<OneD, NekDouble> Pb;
    Array<OneD, Array<OneD, NekDouble>> gradb(expdim * expdim);

    // Communicators to exchange results
    LibUtilities::CommSharedPtr vComm   = pFields[0]->GetComm();
    LibUtilities::CommSharedPtr rowComm = vComm->GetRowComm();
    LibUtilities::CommSharedPtr colComm =
        m_session->DefinesSolverInfo("HomoStrip")
            ? vComm->GetColumnComm()->GetColumnComm()
            : vComm->GetColumnComm();

    // Arrays with forces in each plane
    Array<OneD, NekDouble> Fpplane(m_nPlanes, 0.0);
    Array<OneD, NekDouble> Fvplane(m_nPlanes, 0.0);

    // Forces per element length in a boundary
    Array<OneD, Array<OneD, NekDouble>> fp(expdim);
    Array<OneD, Array<OneD, NekDouble>> fv(expdim);

    // Adjoint gradient per element length in a boundary
    Array<OneD, NekDouble> gradp;
    Array<OneD, NekDouble> gradv;

    // Get viscosity
    NekDouble mu;
    if (m_session->DefinesParameter("Kinvis"))
    {
        NekDouble rho = (m_session->DefinesParameter("rho"))
                            ? (m_session->GetParameter("rho"))
                            : 1;
        mu            = rho * m_session->GetParameter("Kinvis");
    }
    else
    {
        mu = m_session->GetParameter("mu");
    }
    NekDouble lambda = -2.0 / 3.0;

    // Perform BwdTrans: when we only want the mean force in a 3DH1D
    //     we work in wavespace, otherwise we use physical space
    for (i = 0; i < pFields.size(); ++i)
    {
        if (m_isHomogeneous1D && m_outputAllPlanes)
        {
            pFields[i]->SetWaveSpace(false);
        }
        pFields[i]->BwdTrans(pFields[i]->GetCoeffs(), pFields[i]->UpdatePhys());
        pFields[i]->SetPhysState(true);
    }

    // Define boundary expansions
    Array<OneD, const SpatialDomains::BoundaryConditionShPtr> BndConds;
    Array<OneD, MultiRegions::ExpListSharedPtr> BndExp;
    if (m_isHomogeneous1D)
    {
        BndConds = pFields[0]->GetPlane(0)->GetBndConditions();
        BndExp   = pFields[0]->GetPlane(0)->GetBndCondExpansions();
    }
    else
    {
        BndConds = pFields[0]->GetBndConditions();
        BndExp   = pFields[0]->GetBndCondExpansions();
    }

    m_Fp = Array<OneD, NekDouble>(m_nParam, 0.0);
    m_Fv = Array<OneD, NekDouble>(m_nParam, 0.0);
    m_Ft = Array<OneD, NekDouble>(m_nParam, 0.0);

    for (int param = 0; param < m_nParam; param++)
    {

        for (l = 0; l < m_nPlanes; l++)
        {
            Fpplane[l] = 0.0;
            Fvplane[l] = 0.0;
        }

        // For Homogeneous, calculate force on each 2D plane
        // Otherwise, m_nPlanes = 1, and loop only runs once
        for (plane = 0; plane < m_nPlanes; plane++)
        {
            // Check if plane is in this proc
            if (m_planesID[plane] != -1)
            {
                // For Homogeneous, consider the 2D expansion
                //      on the current plane
                if (m_isHomogeneous1D)
                {
                    for (n = 0; n < pFields.size(); n++)
                    {
                        fields[n] = pFields[n]->GetPlane(m_planesID[plane]);
                    }
                }
                else
                {
                    for (n = 0; n < pFields.size(); n++)
                    {
                        fields[n] = pFields[n];
                    }
                }

                // Get velocity and pressure values
                for (n = 0; n < physfields.size(); ++n)
                {
                    physfields[n] = fields[n]->GetPhys();
                }
                for (n = 0; n < nVel; ++n)
                {
                    velocity[n] =
                        Array<OneD, NekDouble>(fields[n]->GetTotPoints());
                }
                pressure = Array<OneD, NekDouble>(fields[0]->GetTotPoints());
                fluidEqu->GetVelocity(physfields, velocity);
                fluidEqu->GetPressure(physfields, pressure);

                // Loop all the Boundary Regions
                for (cnt = n = 0; n < BndConds.size(); n++)
                {
                    if (m_boundaryRegionIsInList[n] == 1)
                    {
                        for (i = 0; i < BndExp[n]->GetExpSize(); ++i, cnt++)
                        {
                            elmtid = m_BCtoElmtID[cnt];
                            elmt   = fields[0]->GetExp(elmtid);
                            nq     = elmt->GetTotPoints();
                            offset = fields[0]->GetPhys_Offset(elmtid);

                            // Extract  fields on this element
                            for (j = 0; j < expdim; j++)
                            {
                                velElmt[j] = velocity[j] + offset;
                            }
                            pElmt = pressure + offset;

                            // Compute the velocity gradients
                            div = Array<OneD, NekDouble>(nq, 0.0);
                            for (j = 0; j < expdim; j++)
                            {
                                for (k = 0; k < expdim; k++)
                                {
                                    grad[j * expdim + k] =
                                        Array<OneD, NekDouble>(nq, 0.0);
                                    elmt->PhysDeriv(k, velElmt[j],
                                                    grad[j * expdim + k]);

                                    if (j == k)
                                    {
                                        Vmath::Vadd(nq, grad[j * expdim + k], 1,
                                                    div, 1, div, 1);
                                    }
                                }
                            }
                            // Scale div by lambda (for compressible flows)
                            Vmath::Smul(nq, lambda, div, 1, div, 1);

                            // identify boundary of element
                            boundary = m_BCtoTraceID[cnt];

                            // Dimension specific part for obtaining values
                            //   at boundary and normal vector
                            Array<OneD, Array<OneD, NekDouble>> normals =
                                elmt->GetTraceNormal(boundary);

                            // Get expansion on boundary
                            LocalRegions::ExpansionSharedPtr bc =
                                BndExp[n]->GetExp(i);

                            // Get number of points on the boundary
                            int nbc = bc->GetTotPoints();

                            // Extract values at boundary
                            Pb = Array<OneD, NekDouble>(nbc, 0.0);
                            elmt->GetTracePhysVals(boundary, bc, pElmt, Pb);

                            for (j = 0; j < expdim * expdim; ++j)
                            {
                                gradb[j] = Array<OneD, NekDouble>(nbc, 0.0);
                                elmt->GetTracePhysVals(boundary, bc, grad[j],
                                                       gradb[j]);
                            }

                            // Calculate forces per unit length

                            // Pressure component: fp[j] = p*n[j]
                            for (j = 0; j < expdim; j++)
                            {
                                fp[j] = Array<OneD, NekDouble>(nbc, 0.0);
                                Vmath::Vmul(nbc, Pb, 1, normals[j], 1, fp[j],
                                            1);
                            }

                            // Viscous component:
                            //     fv[j] = -mu*{(grad[k,j]+grad[j,k]) *n[k]}
                            for (j = 0; j < expdim; j++)
                            {
                                fv[j] = Array<OneD, NekDouble>(nbc, 0.0);
                                for (k = 0; k < expdim; k++)
                                {
                                    //Vmath::Vvtvp(nbc, gradb[k * expdim + j], 1,
                                    //             normals[k], 1, fv[j], 1, fv[j],
                                    //             1);
                                    Vmath::Vvtvp(nbc, gradb[j * expdim + k], 1,
                                                 normals[k], 1, fv[j], 1, fv[j],
                                                 1);
                                }
                                if (!fluidEqu->HasConstantDensity())
                                {
                                    // Add gradient term
                                    Vmath::Vvtvp(nbc, div, 1, normals[j], 1,
                                                 fv[j], 1, fv[j], 1);
                                }
                                Vmath::Smul(nbc, -mu, fv[j], 1, fv[j], 1);
                            }

                            // Multiply adjoint force by the derivative of the
                            // inlet velocity wrt the parameter
                            gradp = Array<OneD, NekDouble>(nbc, 0.0);
                            gradv = Array<OneD, NekDouble>(nbc, 0.0);

                            Array<OneD, NekDouble> x0(nbc, 0.0);
                            Array<OneD, NekDouble> x1(nbc, 0.0);
                            Array<OneD, NekDouble> x2 (nbc, 0.0);
                            Array<OneD, NekDouble> caseduda (nbc);



                            bc->GetCoords(x0, x1, x2);

                            for (k = 0; k < expdim; k++)
                            {

                                LibUtilities::Equation caseduda_equ(m_session->GetInterpreter(), m_caseduda_exp[param][k]);

                                caseduda_equ.Evaluate(x0, x1, x2, time, caseduda);
                                
                                Vmath::Vvtvp(nbc, fv[k], 1, caseduda, 1,
                                             gradv, 1, gradv, 1);
                                Vmath::Vvtvp(nbc, fp[k], 1, caseduda, 1,
                                             gradp, 1, gradp, 1);
                            }

                            // Integrate to obtain force
                            Fpplane[plane] += BndExp[n]->GetExp(i)->Integral(gradp);
                            Fvplane[plane] += BndExp[n]->GetExp(i)->Integral(gradv);
                        }
                    }
                    else
                    {
                        cnt += BndExp[n]->GetExpSize();
                    }
                }
            }
        }

        // Combine contributions from different processes
        //    this is split between row and col comm because of
        //      homostrips case, which only keeps its own strip
        //for (i = 0; i < expdim; i++)
        //{
        //    rowComm->AllReduce(Fpplane[i], LibUtilities::ReduceSum);
        //    colComm->AllReduce(Fpplane[i], LibUtilities::ReduceSum);

        //    rowComm->AllReduce(Fvplane[i], LibUtilities::ReduceSum);
        //    colComm->AllReduce(Fvplane[i], LibUtilities::ReduceSum);
        //}

        m_Fp[param] = (-1)*Vmath::Vsum(m_nPlanes, Fpplane, 1) / m_nPlanes;
        m_Fv[param] = (-1)*Vmath::Vsum(m_nPlanes, Fvplane, 1) / m_nPlanes;
        m_Ft[param] = m_Fp[param] + m_Fv[param];
    }

    // Put results back to wavespace, if necessary
    if (m_isHomogeneous1D && m_outputAllPlanes)
    {
        for (i = 0; i < pFields.size(); ++i)
        {
            pFields[i]->SetWaveSpace(true);
            pFields[i]->HomogeneousFwdTrans(pFields[i]->GetPhys(),
                                            pFields[i]->UpdatePhys());
        }
    }
}

/**
 *     This function calculates the forces when we have a mapping
 *         defining a coordinate system transformation
 */
void FilterAdjointGrad::CalculateForcesMapping(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time)
{
    boost::ignore_unused(time);

    int i, j, k, n, l, cnt, elmtid, offset, boundary, plane;
    // Get number of quadrature points and dimensions
    int physTot = pFields[0]->GetNpoints();
    int expdim = pFields[0]->GetGraph()->GetMeshDimension();
    int nVel = expdim;
    if( m_isHomogeneous1D )
    {
        nVel = nVel + 1;
    }

    LocalRegions::ExpansionSharedPtr elmt;

    // Pressure stress tensor
    //    (global, in a plane, in element and boundary)
    Array<OneD, MultiRegions::ExpListSharedPtr>  P      ( nVel*nVel);
    Array<OneD, MultiRegions::ExpListSharedPtr>  PPlane ( nVel*nVel);
    Array<OneD, Array<OneD, NekDouble> >         PElmt  ( nVel*nVel);
    Array<OneD, Array<OneD, NekDouble> >         PBnd  ( nVel*nVel);
    // Velocity gradient
    Array<OneD, MultiRegions::ExpListSharedPtr>  grad     ( nVel*nVel);
    Array<OneD, MultiRegions::ExpListSharedPtr>  gradPlane( nVel*nVel);
    Array<OneD, Array<OneD, NekDouble> >         gradElmt ( nVel*nVel);
    Array<OneD, Array<OneD, NekDouble> >         gradBnd  ( nVel*nVel);

    // Transformation to cartesian system
    Array<OneD, MultiRegions::ExpListSharedPtr>  C     ( nVel*nVel);
    Array<OneD, MultiRegions::ExpListSharedPtr>  CPlane( nVel*nVel);
    Array<OneD, Array<OneD, NekDouble> >         CElmt ( nVel*nVel);
    Array<OneD, Array<OneD, NekDouble> >         CBnd  ( nVel*nVel);

    // Jacobian
    MultiRegions::ExpListSharedPtr  Jac;
    MultiRegions::ExpListSharedPtr  JacPlane;
    Array<OneD, NekDouble>          JacElmt;
    Array<OneD, NekDouble>          JacBnd;

    // Communicators to exchange results
    LibUtilities::CommSharedPtr vComm = pFields[0]->GetComm();
    LibUtilities::CommSharedPtr rowComm = vComm->GetRowComm();
    LibUtilities::CommSharedPtr colComm =
                            m_session->DefinesSolverInfo("HomoStrip") ?
                                vComm->GetColumnComm()->GetColumnComm():
                                vComm->GetColumnComm();

    // Arrays with forces in each plane
    Array<OneD, NekDouble> Fpplane(m_nPlanes, 0.0);
    Array<OneD, NekDouble> Fvplane(m_nPlanes, 0.0);

    // Forces per element length in a boundary
    Array<OneD, Array<OneD, NekDouble> >       fp( expdim );
    Array<OneD, Array<OneD, NekDouble> >       fv( expdim );

    // Adjoint gradient per element length in a boundary
    Array<OneD, NekDouble> gradp;
    Array<OneD, NekDouble> gradv;

    // Get viscosity
    NekDouble rho = (m_session->DefinesParameter("rho"))
            ? (m_session->GetParameter("rho"))
            : 1;
    NekDouble mu = rho*m_session->GetParameter("Kinvis");

    // Perform BwdTrans: for case with mapping, we can never work
    //                   in wavespace
    for(i = 0; i < pFields.size(); ++i)
    {
        if (m_isHomogeneous1D)
        {
            pFields[i]->SetWaveSpace(false);
        }
        pFields[i]->BwdTrans(pFields[i]->GetCoeffs(),
                             pFields[i]->UpdatePhys());
        pFields[i]->SetPhysState(true);
    }

    // Define boundary expansions
    Array<OneD, const SpatialDomains::BoundaryConditionShPtr > BndConds;
    Array<OneD, MultiRegions::ExpListSharedPtr>  BndExp;
    if(m_isHomogeneous1D)
    {
        BndConds = pFields[0]->GetPlane(0)->GetBndConditions();
        BndExp   = pFields[0]->GetPlane(0)->GetBndCondExpansions();
    }
    else
    {
        BndConds = pFields[0]->GetBndConditions();
        BndExp   = pFields[0]->GetBndCondExpansions();
    }

    //
    // Calculate pressure stress tensor, velocity gradient
    //      and get informations about the mapping

    // Initialise variables
    switch (expdim)
    {
        case 2:
        {
            if (m_isHomogeneous1D)
            {
                MultiRegions::ExpList3DHomogeneous1DSharedPtr Exp3DH1;
                Exp3DH1 = std::dynamic_pointer_cast
                                <MultiRegions::ExpList3DHomogeneous1D>
                                                    (pFields[0]);
                for(i = 0; i < nVel*nVel; i++)
                {
                    grad[i] = MemoryManager<MultiRegions::ExpList3DHomogeneous1D>::
                                AllocateSharedPtr(*Exp3DH1);

                    P[i] = MemoryManager<MultiRegions::ExpList3DHomogeneous1D>::
                                AllocateSharedPtr(*Exp3DH1);

                    C[i] = MemoryManager<MultiRegions::ExpList3DHomogeneous1D>::
                                AllocateSharedPtr(*Exp3DH1);
                }
                Jac = MemoryManager<MultiRegions::ExpList3DHomogeneous1D>::
                                AllocateSharedPtr(*Exp3DH1);
            }
            else
            {
                MultiRegions::ExpListSharedPtr Exp2D;
                Exp2D = std::dynamic_pointer_cast
                                <MultiRegions::ExpList>
                                                    (pFields[0]);
                for(i = 0; i < nVel*nVel; i++)
                {
                    grad[i] = MemoryManager<MultiRegions::ExpList>::
                                AllocateSharedPtr(*Exp2D);

                    P[i] = MemoryManager<MultiRegions::ExpList>::
                                AllocateSharedPtr(*Exp2D);

                    C[i] = MemoryManager<MultiRegions::ExpList>::
                                AllocateSharedPtr(*Exp2D);
                }
                Jac = MemoryManager<MultiRegions::ExpList>::
                                AllocateSharedPtr(*Exp2D);
            }
            break;
        }
        case 3:
        {
            MultiRegions::ExpListSharedPtr Exp3D;
            Exp3D = std::dynamic_pointer_cast
                            <MultiRegions::ExpList>
                                                (pFields[0]);
            for(i = 0; i < nVel*nVel; i++)
            {
                grad[i] = MemoryManager<MultiRegions::ExpList>::
                            AllocateSharedPtr(*Exp3D);

                P[i] = MemoryManager<MultiRegions::ExpList>::
                            AllocateSharedPtr(*Exp3D);

                C[i] = MemoryManager<MultiRegions::ExpList>::
                            AllocateSharedPtr(*Exp3D);
            }
            Jac = MemoryManager<MultiRegions::ExpList>::
                            AllocateSharedPtr(*Exp3D);

            break;
        }
        default:
            ASSERTL0(false,"Expansion dimension not supported by FilterAeroForces");
            break;
    }


    // Get g^ij
    Array<OneD, Array<OneD, NekDouble> >        tmp( nVel*nVel );
    m_mapping->GetInvMetricTensor(tmp);

    // Calculate P^ij = g^ij*p
    for (i = 0; i < nVel*nVel; i++)
    {
        Vmath::Vmul(physTot, tmp[i], 1,
                            pFields[nVel]->GetPhys(), 1,
                            P[i]->UpdatePhys(), 1);
    }

    // Calculate covariant derivatives of U = u^i_,k
    Array<OneD, Array<OneD, NekDouble> >        wk( nVel );
    for (i=0; i<nVel; i++)
    {
        wk[i] = Array<OneD, NekDouble>(physTot, 0.0);
        Vmath::Vcopy(physTot, pFields[i]->GetPhys(), 1,
                            wk[i], 1);
    }
    m_mapping->ApplyChristoffelContravar(wk, tmp);
    for (i=0; i< nVel; i++)
    {
        for (k=0; k< nVel; k++)
        {
            pFields[0]->PhysDeriv(MultiRegions::DirCartesianMap[k],
                                   wk[i], grad[i*nVel+k]->UpdatePhys());

            Vmath::Vadd(physTot,tmp[i*nVel+k],1,
                                grad[i*nVel+k]->UpdatePhys(), 1,
                                grad[i*nVel+k]->UpdatePhys(), 1);
        }
    }
    // Raise index to obtain Grad^ij = g^jk u^i_,k
    for (i=0; i< nVel; i++)
    {
        for (k=0; k< nVel; k++)
        {
            Vmath::Vcopy(physTot, grad[i*nVel+k]->GetPhys(), 1,
                                  wk[k], 1);
        }
        m_mapping->RaiseIndex(wk, wk);
        for (j=0; j<nVel; j++)
        {
            Vmath::Vcopy(physTot, wk[j], 1,
                                  grad[i*nVel+j]->UpdatePhys(), 1);
        }
    }

    // Get Jacobian
    m_mapping->GetJacobian( Jac->UpdatePhys());

    std::cout << "bla4" << endl;

    // Get transformation to Cartesian system
    for (i=0; i< nVel; i++)
    {
        // Zero wk storage
        for (k=0; k< nVel; k++)
        {
            wk[k] = Array<OneD, NekDouble>(physTot, 0.0);
        }
        // Make wk[i] = 1
        wk[i] = Array<OneD, NekDouble>(physTot, 1.0);
        // Transform wk to Cartesian
        m_mapping->ContravarToCartesian(wk,wk);
        // Copy result to a column in C
        for (k=0; k< nVel; k++)
        {
            Vmath::Vcopy(physTot, wk[k], 1,
                                  C[k*nVel+i]->UpdatePhys(), 1);
        }
    }

    //
    // Calculate forces
    //

    m_Fp = Array<OneD, NekDouble>(m_nParam, 0.0);
    m_Fv = Array<OneD, NekDouble>(m_nParam, 0.0);
    m_Ft = Array<OneD, NekDouble>(m_nParam, 0.0);

    for (int param = 0; param < m_nParam; param++)
    {

        for (l = 0; l < m_nPlanes; l++)
        {
            Fpplane[l] = 0.0;
            Fvplane[l] = 0.0;
        }


    std::cout << "blaA" << endl;


    // For Homogeneous, calculate force on each 2D plane
    // Otherwise, m_nPlanes = 1, and loop only runs once
    for(plane = 0; plane < m_nPlanes; plane++ )
    {
        // Check if plane is in this proc
        if( m_planesID[plane] != -1 )
        {
            // For Homogeneous, consider the 2D expansion
            //      on the current plane
            if(m_isHomogeneous1D)
            {
                for(n = 0; n < nVel*nVel; n++)
                {
                    PPlane[n]    = P[n]->GetPlane(m_planesID[plane]);
                    gradPlane[n] = grad[n]->GetPlane(m_planesID[plane]);
                    CPlane[n]    = C[n]->GetPlane(m_planesID[plane]);
                }
                JacPlane = Jac->GetPlane(m_planesID[plane]);
            }
            else
            {
                for(n = 0; n < nVel*nVel; n++)
                {
                    PPlane[n]    = P[n];
                    gradPlane[n] = grad[n];
                    CPlane[n] = C[n];
                }
                JacPlane = Jac;
            }


            std::cout << "blaB" << endl;

            //Loop all the Boundary Regions
            for( cnt = n = 0; n < BndConds.size(); n++)
            {
                if(m_boundaryRegionIsInList[n] == 1)
                {
                    for (i=0; i < BndExp[n]->GetExpSize(); ++i,cnt++)
                    {
                        elmtid = m_BCtoElmtID[cnt];
                        elmt   = PPlane[0]->GetExp(elmtid);
                        offset = PPlane[0]->GetPhys_Offset(elmtid);

                        // Extract  fields on this element
                        for( j=0; j<nVel*nVel; j++)
                        {
                            PElmt[j]    = PPlane[j]->GetPhys()
                                        + offset;
                            gradElmt[j] = gradPlane[j]->GetPhys()
                                        + offset;
                            CElmt[j]    = CPlane[j]->GetPhys()
                                        + offset;
                        }
                        JacElmt = JacPlane->GetPhys() + offset;

                        // identify boundary of element
                        boundary = m_BCtoTraceID[cnt];

                        // Dimension specific part for obtaining values
                        //   at boundary and normal vector
                        Array<OneD, Array<OneD, NekDouble> > normals;
                        // Get normals
                        normals = elmt->GetTraceNormal(boundary);

                        // Get expansion on boundary
                        LocalRegions::ExpansionSharedPtr bc = 
                            BndExp[n]->GetExp(i);
                        
                        // Get number of points on the boundary
                        int nbc = bc->GetTotPoints();

                        // Extract values at boundary
                        for(j = 0; j < nVel*nVel; ++j)
                        {
                            gradBnd[j] = Array<OneD, NekDouble>(nbc,0.0);
                            elmt->GetTracePhysVals(boundary,bc,gradElmt[j],gradBnd[j]);
                            
                            PBnd[j] = Array<OneD, NekDouble>(nbc,0.0);
                            elmt->GetTracePhysVals(boundary,bc,PElmt[j],PBnd[j]);

                            CBnd[j] = Array<OneD, NekDouble>(nbc,0.0);
                            elmt->GetTracePhysVals(boundary,bc,CElmt[j],CBnd[j]);
                        }

                        JacBnd = Array<OneD, NekDouble>(nbc,0.0);
                        elmt->GetTracePhysVals(boundary,bc,JacElmt,JacBnd);


                        std::cout << "blaC" << endl;


                        // Calculate forces per unit length

                        // Pressure component: fp[j] = P[j,k]*n[k]
                        for ( j = 0; j < nVel; j++)
                        {
                            fp[j] = Array<OneD, NekDouble> (nbc,0.0);
                            // Normals only has expdim elements
                            for ( k = 0; k < expdim; k++)
                            {
                                Vmath::Vvtvp (nbc, PBnd[ j*nVel + k], 1,
                                                   normals[k], 1,
                                                   fp[j], 1,
                                                   fp[j], 1);
                            }
                        }

                        // Viscous component:
                        //     fv[j] = -mu*{(grad[k,j]+grad[j,k]) *n[k]}
                        for ( j = 0; j < nVel; j++ )
                        {
                            fv[j] = Array<OneD, NekDouble> (nbc,0.0);
                            for ( k = 0; k < expdim; k++ )
                            {
                                Vmath::Vvtvp (nbc,gradBnd[k*nVel+j],1,
                                                   normals[k], 1,
                                                   fv[j], 1,
                                                   fv[j], 1);
                                Vmath::Vvtvp (nbc,gradBnd[j*nVel+k],1,
                                                   normals[k], 1,
                                                   fv[j], 1,
                                                   fv[j], 1);
                            }
                            Vmath::Smul(nbc, -mu, fv[j], 1, fv[j], 1);
                        }

                        // Multiply by Jacobian
                        for ( k = 0; k < nVel; k++ )
                        {
                            Vmath::Vmul(nbc, JacBnd, 1, fp[k], 1,
                                                        fp[k], 1);
                            Vmath::Vmul(nbc, JacBnd, 1, fv[k], 1,
                                                        fv[k], 1);
                        }

                        // Convert to cartesian system
                        for ( k = 0; k < nVel; k++ )
                        {
                            wk[k] = Array<OneD, NekDouble>(physTot,0.0);
                            for ( j = 0; j < nVel; j++ )
                            {
                                Vmath::Vvtvp(nbc, CBnd[k*nVel+j], 1,
                                                    fp[j], 1,
                                                    wk[k], 1,
                                                    wk[k], 1);
                            }
                        }
                        for ( k = 0; k < nVel; k++ )
                        {
                            Vmath::Vcopy(nbc, wk[k], 1, fp[k], 1);
                        }

                        for ( k = 0; k < nVel; k++ )
                        {
                            wk[k] = Array<OneD, NekDouble>(physTot,0.0);
                            for ( j = 0; j < nVel; j++ )
                            {
                                Vmath::Vvtvp(nbc, CBnd[k*nVel+j], 1,
                                                    fv[j], 1,
                                                    wk[k], 1,
                                                    wk[k], 1);
                            }
                        }
                        for ( k = 0; k < nVel; k++ )
                        {
                            Vmath::Vcopy(nbc, wk[k], 1, fv[k], 1);
                        }


                        // Multiply adjoint force by the derivative of the
                            // inlet velocity wrt the parameter
                            gradp = Array<OneD, NekDouble>(nbc, 0.0);
                            gradv = Array<OneD, NekDouble>(nbc, 0.0);

                            Array<OneD, NekDouble> x0(nbc, 0.0);
                            Array<OneD, NekDouble> x1(nbc, 0.0);
                            Array<OneD, NekDouble> x2 (nbc, 0.0);
                            Array<OneD, NekDouble> caseduda (nbc);


                            std::cout << "blaX" << endl;



                            bc->GetCoords(x0, x1, x2);

                            for (k = 0; k < expdim; k++)
                            {
                                std::cout << "blaX1" << endl;

                                LibUtilities::Equation caseduda_equ(m_session->GetInterpreter(), m_caseduda_exp[param][k]);

                                std::cout << "blaX2" << endl;

                                caseduda_equ.Evaluate(x0, x1, x2, time, caseduda);


                                std::cout << "blaY" << endl;
                                
                                Vmath::Vvtvp(nbc, fv[k], 1, caseduda, 1,
                                             gradv, 1, gradv, 1);
                                Vmath::Vvtvp(nbc, fp[k], 1, caseduda, 1,
                                             gradp, 1, gradp, 1);
                            }

                            std::cout << "blaZ" << endl;

                        
                            // Integrate to obtain force
                            Fpplane[plane] += BndExp[n]->GetExp(i)->Integral(gradp);
                            Fvplane[plane] += BndExp[n]->GetExp(i)->Integral(gradv);
                    }
                }
                else
                {
                    cnt += BndExp[n]->GetExpSize();
                }
            }
        }
    }

    std::cout << "bla5" << endl;

    // Combine contributions from different processes
        //    this is split between row and col comm because of
        //      homostrips case, which only keeps its own strip
        //for (i = 0; i < expdim; i++)
        //{
        //    rowComm->AllReduce(Fpplane[i], LibUtilities::ReduceSum);
        //    colComm->AllReduce(Fpplane[i], LibUtilities::ReduceSum);

        //    rowComm->AllReduce(Fvplane[i], LibUtilities::ReduceSum);
        //    colComm->AllReduce(Fvplane[i], LibUtilities::ReduceSum);
        //}

        m_Fp[param] = (-1)*Vmath::Vsum(m_nPlanes, Fpplane, 1) / m_nPlanes;
        m_Fv[param] = (-1)*Vmath::Vsum(m_nPlanes, Fvplane, 1) / m_nPlanes;
        m_Ft[param] = m_Fp[param] + m_Fv[param];
    }

    // Put results back to wavespace, if necessary
    if( m_isHomogeneous1D)
    {
        for (i = 0; i < pFields.size(); ++i)
        {
            pFields[i]->SetWaveSpace(true);
            pFields[i]->HomogeneousFwdTrans(pFields[i]->GetPhys(),
                                            pFields[i]->UpdatePhys());
        }
    }
}


void FilterAdjointGrad::GetCaseInfo() {


    if (m_case == "snake")
    {
        m_nParam = 1;
        m_caseduda_exp.resize(m_nParam);
        for (int i = 0; i < m_nParam; i++) {
            m_caseduda_exp[i].resize(3);
        }
        m_caseduda_exp[0][0] = "1.0";
        m_caseduda_exp[0][1] = "0.0";
        m_caseduda_exp[0][2] = "0.0";
    }
    else if (m_case == "cylsnake")
    {
        m_nParam = 1;
        m_caseduda_exp.resize(m_nParam);
        for (int i = 0; i < m_nParam; i++) {
            m_caseduda_exp[i].resize(3);
        }
        m_caseduda_exp[0][0] = "1.5*(1.0-y*y/h/h)";
        m_caseduda_exp[0][1] = "0.0";
        m_caseduda_exp[0][2] = "0.0";
    }
    else if (m_case == "aeroforces")
    {
        m_nParam = 1;
        m_caseduda_exp.resize(m_nParam);
        for (int i = 0; i < m_nParam; i++) {
            m_caseduda_exp[i].resize(3);
        }
        m_caseduda_exp[0][0] = "1.0";
        m_caseduda_exp[0][1] = "0.0";
        m_caseduda_exp[0][2] = "0.0";
    }
    else if (m_case == "periodic")
    {
        m_nParam = 2;
        m_caseduda_exp.resize(m_nParam);
        for (int i = 0; i < m_nParam; i++) {
            m_caseduda_exp[i].resize(3);
        }
        //m_caseduda_exp[0][0] = "cos(freq*(adjperiod-t))*freq*(adjperiod)";
        m_caseduda_exp[0][0] = "sin(freq*(adjperiod-t))";
        m_caseduda_exp[0][1] = "0.0";
        m_caseduda_exp[0][2] = "0.0";
        m_caseduda_exp[1][0] = "cos(freq*(adjperiod-t))*freq*(adjperiod-t)";
        m_caseduda_exp[1][1] = "0.0";
        m_caseduda_exp[1][2] = "0.0";
    }
    else if (m_case == "cavity")
    {
        m_nParam = 2;
        m_caseduda_exp.resize(m_nParam);
        for (int i = 0; i < m_nParam; i++) {
            m_caseduda_exp[i].resize(3);
        }
        //m_caseduda_exp[0][0] = "cos(freq*(adjperiod-t))*freq*(adjperiod)";
        m_caseduda_exp[0][0] = "(1-(2*x-1)*(2*x-1))*(1-(2*x-1)*(2*x-1))*sin(freq*(adjperiod-t))";
        m_caseduda_exp[0][1] = "0.0";
        m_caseduda_exp[0][2] = "0.0";
        m_caseduda_exp[1][0] = "(1-(2*x-1)*(2*x-1))*(1-(2*x-1)*(2*x-1))*cos(freq*(adjperiod-t))*freq*(adjperiod-t)";
        m_caseduda_exp[1][1] = "0.0";
        m_caseduda_exp[1][2] = "0.0";
    }
    else if (m_case == "newuperiodic")
    {
        m_nParam = 2;
        m_caseduda_exp.resize(m_nParam);
        for (int i = 0; i < m_nParam; i++) {
            m_caseduda_exp[i].resize(3);
        }
        m_caseduda_exp[0][0] = "-sin(freq*(adjperiod-t))/freq";
        m_caseduda_exp[0][1] = "0.0";
        m_caseduda_exp[0][2] = "0.0";
        m_caseduda_exp[1][0] = "-cos(freq*(adjperiod-t))*(adjperiod-t) + sin(freq*(adjperiod-t))/freq";
        m_caseduda_exp[1][1] = "0.0";
        m_caseduda_exp[1][2] = "0.0";
    }
    else if (m_case == "snakeperiodic")
    {
        m_nParam = 2;
        m_caseduda_exp.resize(m_nParam);
        for (int i = 0; i < m_nParam; i++) {
            m_caseduda_exp[i].resize(3);
        }
        m_caseduda_exp[0][0] = "-1*((cos(1/2*sqrt(2)*sqrt(B*Re*freq)*y)*cosh(1/2*sqrt(2)*sqrt(B*Re*freq)*y)*sin(1/2*sqrt(2)*sqrt(B*Re*freq)*h)*sinh(1/2*sqrt(2)*sqrt(B*Re*freq)*h) - cos(1/2*sqrt(2)*sqrt(B*Re*freq)*h)*cosh(1/2*sqrt(2)*sqrt(B*Re*freq)*h)*sin(1/2*sqrt(2)*sqrt(B*Re*freq)*y)*sinh(1/2*sqrt(2)*sqrt(B*Re*freq)*y))*cos(B*(adjperiod-t)*freq) + (cos(1/2*sqrt(2)*sqrt(B*Re*freq)*h)^2*cosh(1/2*sqrt(2)*sqrt(B*Re*freq)*h)^2 - cos(1/2*sqrt(2)*sqrt(B*Re*freq)*h)*cos(1/2*sqrt(2)*sqrt(B*Re*freq)*y)*cosh(1/2*sqrt(2)*sqrt(B*Re*freq)*h)*cosh(1/2*sqrt(2)*sqrt(B*Re*freq)*y) + sin(1/2*sqrt(2)*sqrt(B*Re*freq)*h)^2*sinh(1/2*sqrt(2)*sqrt(B*Re*freq)*h)^2 - sin(1/2*sqrt(2)*sqrt(B*Re*freq)*h)*sin(1/2*sqrt(2)*sqrt(B*Re*freq)*y)*sinh(1/2*sqrt(2)*sqrt(B*Re*freq)*h)*sinh(1/2*sqrt(2)*sqrt(B*Re*freq)*y))*sin(B*(adjperiod-t)*freq))/((cos(1/2*sqrt(2)*sqrt(B*Re*freq)*h)^2*cosh(1/2*sqrt(2)*sqrt(B*Re*freq)*h)^2 + sin(1/2*sqrt(2)*sqrt(B*Re*freq)*h)^2*sinh(1/2*sqrt(2)*sqrt(B*Re*freq)*h)^2)*B*freq)" ;
        m_caseduda_exp[0][1] = "0.0";
        m_caseduda_exp[0][2] = "0.0";
        m_caseduda_exp[1][0] = "A*(-1/4*(4*(cos(1/2*sqrt(2)*sqrt(B*Re*freq)*h)^2*cosh(1/2*sqrt(2)*sqrt(B*Re*freq)*h)^2 - cos(1/2*sqrt(2)*sqrt(B*Re*freq)*h)*cos(1/2*sqrt(2)*sqrt(B*Re*freq)*y)*cosh(1/2*sqrt(2)*sqrt(B*Re*freq)*h)*cosh(1/2*sqrt(2)*sqrt(B*Re*freq)*y) + sin(1/2*sqrt(2)*sqrt(B*Re*freq)*h)^2*sinh(1/2*sqrt(2)*sqrt(B*Re*freq)*h)^2 - sin(1/2*sqrt(2)*sqrt(B*Re*freq)*h)*sin(1/2*sqrt(2)*sqrt(B*Re*freq)*y)*sinh(1/2*sqrt(2)*sqrt(B*Re*freq)*h)*sinh(1/2*sqrt(2)*sqrt(B*Re*freq)*y))*(adjperiod-t)*freq*cos(B*(adjperiod-t)*freq) - 4*(cos(1/2*sqrt(2)*sqrt(B*Re*freq)*y)*cosh(1/2*sqrt(2)*sqrt(B*Re*freq)*y)*sin(1/2*sqrt(2)*sqrt(B*Re*freq)*h)*sinh(1/2*sqrt(2)*sqrt(B*Re*freq)*h) - cos(1/2*sqrt(2)*sqrt(B*Re*freq)*h)*cosh(1/2*sqrt(2)*sqrt(B*Re*freq)*h)*sin(1/2*sqrt(2)*sqrt(B*Re*freq)*y)*sinh(1/2*sqrt(2)*sqrt(B*Re*freq)*y))*(adjperiod-t)*freq*sin(B*(adjperiod-t)*freq) + (sqrt(2)*Re*h*freq*cos(1/2*sqrt(2)*sqrt(B*Re*freq)*y)*cosh(1/2*sqrt(2)*sqrt(B*Re*freq)*h)*cosh(1/2*sqrt(2)*sqrt(B*Re*freq)*y)*sin(1/2*sqrt(2)*sqrt(B*Re*freq)*h)/sqrt(B*Re*freq) - sqrt(2)*Re*freq*y*cos(1/2*sqrt(2)*sqrt(B*Re*freq)*h)*cosh(1/2*sqrt(2)*sqrt(B*Re*freq)*h)*cosh(1/2*sqrt(2)*sqrt(B*Re*freq)*y)*sin(1/2*sqrt(2)*sqrt(B*Re*freq)*y)/sqrt(B*Re*freq) + sqrt(2)*Re*h*freq*cos(1/2*sqrt(2)*sqrt(B*Re*freq)*h)*cos(1/2*sqrt(2)*sqrt(B*Re*freq)*y)*cosh(1/2*sqrt(2)*sqrt(B*Re*freq)*y)*sinh(1/2*sqrt(2)*sqrt(B*Re*freq)*h)/sqrt(B*Re*freq) - sqrt(2)*Re*freq*y*cosh(1/2*sqrt(2)*sqrt(B*Re*freq)*y)*sin(1/2*sqrt(2)*sqrt(B*Re*freq)*h)*sin(1/2*sqrt(2)*sqrt(B*Re*freq)*y)*sinh(1/2*sqrt(2)*sqrt(B*Re*freq)*h)/sqrt(B*Re*freq) - sqrt(2)*Re*freq*y*cos(1/2*sqrt(2)*sqrt(B*Re*freq)*h)*cos(1/2*sqrt(2)*sqrt(B*Re*freq)*y)*cosh(1/2*sqrt(2)*sqrt(B*Re*freq)*h)*sinh(1/2*sqrt(2)*sqrt(B*Re*freq)*y)/sqrt(B*Re*freq) + sqrt(2)*Re*h*freq*cosh(1/2*sqrt(2)*sqrt(B*Re*freq)*h)*sin(1/2*sqrt(2)*sqrt(B*Re*freq)*h)*sin(1/2*sqrt(2)*sqrt(B*Re*freq)*y)*sinh(1/2*sqrt(2)*sqrt(B*Re*freq)*y)/sqrt(B*Re*freq) + sqrt(2)*Re*freq*y*cos(1/2*sqrt(2)*sqrt(B*Re*freq)*y)*sin(1/2*sqrt(2)*sqrt(B*Re*freq)*h)*sinh(1/2*sqrt(2)*sqrt(B*Re*freq)*h)*sinh(1/2*sqrt(2)*sqrt(B*Re*freq)*y)/sqrt(B*Re*freq) - sqrt(2)*Re*h*freq*cos(1/2*sqrt(2)*sqrt(B*Re*freq)*h)*sin(1/2*sqrt(2)*sqrt(B*Re*freq)*y)*sinh(1/2*sqrt(2)*sqrt(B*Re*freq)*h)*sinh(1/2*sqrt(2)*sqrt(B*Re*freq)*y)/sqrt(B*Re*freq))*cos(B*(adjperiod-t)*freq) - (2*sqrt(2)*Re*h*freq*cos(1/2*sqrt(2)*sqrt(B*Re*freq)*h)*cosh(1/2*sqrt(2)*sqrt(B*Re*freq)*h)^2*sin(1/2*sqrt(2)*sqrt(B*Re*freq)*h)/sqrt(B*Re*freq) - sqrt(2)*Re*h*freq*cos(1/2*sqrt(2)*sqrt(B*Re*freq)*y)*cosh(1/2*sqrt(2)*sqrt(B*Re*freq)*h)*cosh(1/2*sqrt(2)*sqrt(B*Re*freq)*y)*sin(1/2*sqrt(2)*sqrt(B*Re*freq)*h)/sqrt(B*Re*freq) - sqrt(2)*Re*freq*y*cos(1/2*sqrt(2)*sqrt(B*Re*freq)*h)*cosh(1/2*sqrt(2)*sqrt(B*Re*freq)*h)*cosh(1/2*sqrt(2)*sqrt(B*Re*freq)*y)*sin(1/2*sqrt(2)*sqrt(B*Re*freq)*y)/sqrt(B*Re*freq) - 2*sqrt(2)*Re*h*freq*cos(1/2*sqrt(2)*sqrt(B*Re*freq)*h)^2*cosh(1/2*sqrt(2)*sqrt(B*Re*freq)*h)*sinh(1/2*sqrt(2)*sqrt(B*Re*freq)*h)/sqrt(B*Re*freq) + sqrt(2)*Re*h*freq*cos(1/2*sqrt(2)*sqrt(B*Re*freq)*h)*cos(1/2*sqrt(2)*sqrt(B*Re*freq)*y)*cosh(1/2*sqrt(2)*sqrt(B*Re*freq)*y)*sinh(1/2*sqrt(2)*sqrt(B*Re*freq)*h)/sqrt(B*Re*freq) - 2*sqrt(2)*Re*h*freq*cosh(1/2*sqrt(2)*sqrt(B*Re*freq)*h)*sin(1/2*sqrt(2)*sqrt(B*Re*freq)*h)^2*sinh(1/2*sqrt(2)*sqrt(B*Re*freq)*h)/sqrt(B*Re*freq) + sqrt(2)*Re*freq*y*cosh(1/2*sqrt(2)*sqrt(B*Re*freq)*y)*sin(1/2*sqrt(2)*sqrt(B*Re*freq)*h)*sin(1/2*sqrt(2)*sqrt(B*Re*freq)*y)*sinh(1/2*sqrt(2)*sqrt(B*Re*freq)*h)/sqrt(B*Re*freq) - 2*sqrt(2)*Re*h*freq*cos(1/2*sqrt(2)*sqrt(B*Re*freq)*h)*sin(1/2*sqrt(2)*sqrt(B*Re*freq)*h)*sinh(1/2*sqrt(2)*sqrt(B*Re*freq)*h)^2/sqrt(B*Re*freq) + sqrt(2)*Re*freq*y*cos(1/2*sqrt(2)*sqrt(B*Re*freq)*h)*cos(1/2*sqrt(2)*sqrt(B*Re*freq)*y)*cosh(1/2*sqrt(2)*sqrt(B*Re*freq)*h)*sinh(1/2*sqrt(2)*sqrt(B*Re*freq)*y)/sqrt(B*Re*freq) + sqrt(2)*Re*h*freq*cosh(1/2*sqrt(2)*sqrt(B*Re*freq)*h)*sin(1/2*sqrt(2)*sqrt(B*Re*freq)*h)*sin(1/2*sqrt(2)*sqrt(B*Re*freq)*y)*sinh(1/2*sqrt(2)*sqrt(B*Re*freq)*y)/sqrt(B*Re*freq) + sqrt(2)*Re*freq*y*cos(1/2*sqrt(2)*sqrt(B*Re*freq)*y)*sin(1/2*sqrt(2)*sqrt(B*Re*freq)*h)*sinh(1/2*sqrt(2)*sqrt(B*Re*freq)*h)*sinh(1/2*sqrt(2)*sqrt(B*Re*freq)*y)/sqrt(B*Re*freq) + sqrt(2)*Re*h*freq*cos(1/2*sqrt(2)*sqrt(B*Re*freq)*h)*sin(1/2*sqrt(2)*sqrt(B*Re*freq)*y)*sinh(1/2*sqrt(2)*sqrt(B*Re*freq)*h)*sinh(1/2*sqrt(2)*sqrt(B*Re*freq)*y)/sqrt(B*Re*freq))*sin(B*(adjperiod-t)*freq))/((cos(1/2*sqrt(2)*sqrt(B*Re*freq)*h)^2*cosh(1/2*sqrt(2)*sqrt(B*Re*freq)*h)^2 + sin(1/2*sqrt(2)*sqrt(B*Re*freq)*h)^2*sinh(1/2*sqrt(2)*sqrt(B*Re*freq)*h)^2)*B*freq) - 1/2*(sqrt(2)*Re*h*freq*cos(1/2*sqrt(2)*sqrt(B*Re*freq)*h)*cosh(1/2*sqrt(2)*sqrt(B*Re*freq)*h)^2*sin(1/2*sqrt(2)*sqrt(B*Re*freq)*h)/sqrt(B*Re*freq) - sqrt(2)*Re*h*freq*cos(1/2*sqrt(2)*sqrt(B*Re*freq)*h)^2*cosh(1/2*sqrt(2)*sqrt(B*Re*freq)*h)*sinh(1/2*sqrt(2)*sqrt(B*Re*freq)*h)/sqrt(B*Re*freq) - sqrt(2)*Re*h*freq*cosh(1/2*sqrt(2)*sqrt(B*Re*freq)*h)*sin(1/2*sqrt(2)*sqrt(B*Re*freq)*h)^2*sinh(1/2*sqrt(2)*sqrt(B*Re*freq)*h)/sqrt(B*Re*freq) - sqrt(2)*Re*h*freq*cos(1/2*sqrt(2)*sqrt(B*Re*freq)*h)*sin(1/2*sqrt(2)*sqrt(B*Re*freq)*h)*sinh(1/2*sqrt(2)*sqrt(B*Re*freq)*h)^2/sqrt(B*Re*freq))*((cos(1/2*sqrt(2)*sqrt(B*Re*freq)*y)*cosh(1/2*sqrt(2)*sqrt(B*Re*freq)*y)*sin(1/2*sqrt(2)*sqrt(B*Re*freq)*h)*sinh(1/2*sqrt(2)*sqrt(B*Re*freq)*h) - cos(1/2*sqrt(2)*sqrt(B*Re*freq)*h)*cosh(1/2*sqrt(2)*sqrt(B*Re*freq)*h)*sin(1/2*sqrt(2)*sqrt(B*Re*freq)*y)*sinh(1/2*sqrt(2)*sqrt(B*Re*freq)*y))*cos(B*(adjperiod-t)*freq) + (cos(1/2*sqrt(2)*sqrt(B*Re*freq)*h)^2*cosh(1/2*sqrt(2)*sqrt(B*Re*freq)*h)^2 - cos(1/2*sqrt(2)*sqrt(B*Re*freq)*h)*cos(1/2*sqrt(2)*sqrt(B*Re*freq)*y)*cosh(1/2*sqrt(2)*sqrt(B*Re*freq)*h)*cosh(1/2*sqrt(2)*sqrt(B*Re*freq)*y) + sin(1/2*sqrt(2)*sqrt(B*Re*freq)*h)^2*sinh(1/2*sqrt(2)*sqrt(B*Re*freq)*h)^2 - sin(1/2*sqrt(2)*sqrt(B*Re*freq)*h)*sin(1/2*sqrt(2)*sqrt(B*Re*freq)*y)*sinh(1/2*sqrt(2)*sqrt(B*Re*freq)*h)*sinh(1/2*sqrt(2)*sqrt(B*Re*freq)*y))*sin(B*(adjperiod-t)*freq))/((cos(1/2*sqrt(2)*sqrt(B*Re*freq)*h)^2*cosh(1/2*sqrt(2)*sqrt(B*Re*freq)*h)^2 + sin(1/2*sqrt(2)*sqrt(B*Re*freq)*h)^2*sinh(1/2*sqrt(2)*sqrt(B*Re*freq)*h)^2)^2*B*freq) + ((cos(1/2*sqrt(2)*sqrt(B*Re*freq)*y)*cosh(1/2*sqrt(2)*sqrt(B*Re*freq)*y)*sin(1/2*sqrt(2)*sqrt(B*Re*freq)*h)*sinh(1/2*sqrt(2)*sqrt(B*Re*freq)*h) - cos(1/2*sqrt(2)*sqrt(B*Re*freq)*h)*cosh(1/2*sqrt(2)*sqrt(B*Re*freq)*h)*sin(1/2*sqrt(2)*sqrt(B*Re*freq)*y)*sinh(1/2*sqrt(2)*sqrt(B*Re*freq)*y))*cos(B*(adjperiod-t)*freq) + (cos(1/2*sqrt(2)*sqrt(B*Re*freq)*h)^2*cosh(1/2*sqrt(2)*sqrt(B*Re*freq)*h)^2 - cos(1/2*sqrt(2)*sqrt(B*Re*freq)*h)*cos(1/2*sqrt(2)*sqrt(B*Re*freq)*y)*cosh(1/2*sqrt(2)*sqrt(B*Re*freq)*h)*cosh(1/2*sqrt(2)*sqrt(B*Re*freq)*y) + sin(1/2*sqrt(2)*sqrt(B*Re*freq)*h)^2*sinh(1/2*sqrt(2)*sqrt(B*Re*freq)*h)^2 - sin(1/2*sqrt(2)*sqrt(B*Re*freq)*h)*sin(1/2*sqrt(2)*sqrt(B*Re*freq)*y)*sinh(1/2*sqrt(2)*sqrt(B*Re*freq)*h)*sinh(1/2*sqrt(2)*sqrt(B*Re*freq)*y))*sin(B*(adjperiod-t)*freq))/((cos(1/2*sqrt(2)*sqrt(B*Re*freq)*h)^2*cosh(1/2*sqrt(2)*sqrt(B*Re*freq)*h)^2 + sin(1/2*sqrt(2)*sqrt(B*Re*freq)*h)^2*sinh(1/2*sqrt(2)*sqrt(B*Re*freq)*h)^2)*B^2*freq))";
        m_caseduda_exp[1][1] = "0.0";
        m_caseduda_exp[1][2] = "0.0";
    }
    else
    {    
        ASSERTL0(1==0,   "Given 'Case' parameter does not match any of the implemented cases.");
    
    }
}





}
}
