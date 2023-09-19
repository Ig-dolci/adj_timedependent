///////////////////////////////////////////////////////////////////////////////
//
// File FilterAdjointGrad.h
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
// Description: Outputs values at specific points during time-stepping.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_FILTERS_FILTERFORCES_H
#define NEKTAR_SOLVERUTILS_FILTERS_FILTERFORCES_H

#include <SolverUtils/Filters/Filter.h>
#include <LocalRegions/Expansion3D.h>
#include <GlobalMapping/Mapping.h>

namespace Nektar
{
namespace SolverUtils
{
class FilterAdjointGrad : public Filter
{
public:
    friend class MemoryManager<FilterAdjointGrad>;

    /// Creates an instance of this class
    static FilterSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const std::weak_ptr<EquationSystem>      &pEquation,
        const std::map<std::string, std::string>   &pParams)
    {
        FilterSharedPtr p = MemoryManager<FilterAdjointGrad>::
                                AllocateSharedPtr(pSession, pEquation, pParams);
        return p;
    }

    ///Name of the class
    static std::string className;

    SOLVER_UTILS_EXPORT FilterAdjointGrad(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const std::weak_ptr<EquationSystem>      &pEquation,
        const std::map<std::string, std::string>   &pParams);

    SOLVER_UTILS_EXPORT virtual ~FilterAdjointGrad();


protected:
    virtual void v_Initialise(
            const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
            const NekDouble &time);
    virtual void v_Update(
            const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
            const NekDouble &time);
    virtual void v_Finalise(
            const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
            const NekDouble &time);
    virtual bool v_IsTimeDependent();

private:
    /// ID's of boundary regions where we want the forces
    std::vector<unsigned int>       m_boundaryRegionsIdList;
    /// Determines if a given Boundary Region is in
    /// m_boundaryRegionsIdList
    std::vector<bool>               m_boundaryRegionIsInList;
    unsigned int                    m_index;
    unsigned int                    m_outputFrequency;
    /// if using a homogeneous1D expansion, determine if should output
    ///     all planes or just the average
    bool                            m_outputAllPlanes;
    bool                            m_isHomogeneous1D;
    std::string                     m_outputFile;
    std::ofstream                   m_outputStream;
    LibUtilities::BasisSharedPtr    m_homogeneousBasis;
    std::string                     m_BoundaryString;
    Array<OneD, int>                m_BCtoElmtID;
    Array<OneD, int>                m_BCtoTraceID;
    /// number of planes for homogeneous1D expansion
    int                             m_nPlanes;
    Array<OneD, int>                m_planesID;
    // Time when we start calculating the forces
    NekDouble                       m_startTime;
    // Directions on which the forces will be projected
    Array<OneD, Array<OneD, NekDouble> >    m_directions;
    // Arrays storing the forces calculated at current timestep
    Array<OneD, NekDouble>          m_Fp;
    Array<OneD, NekDouble>          m_Fv;
    Array<OneD, NekDouble>          m_Ft;
    // Array storing the forces calculated at previous timesteps
    Array<OneD, NekDouble>          m_Ftp;
    Array<OneD, NekDouble>          m_Ftpp;
    
    // Adjoint gradient
    std::string                     m_case;
    int                             m_nParam;
    std::vector<std::vector<std::string>>    m_caseduda_exp;
    // First order time integration adjoint gradient (trapezoidal rule)
    Array<OneD, NekDouble>          m_adjgrad;
    // Second order time integration adjoint gradient (Simpson's rule)
    Array<OneD, NekDouble>          m_adjgrad2;
    int                             m_aux;
    

    NekDouble                       m_lastTime;
    GlobalMapping::MappingSharedPtr m_mapping;

    void CalculateForces(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time);

    void CalculateForcesMapping(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time);

    void GetCaseInfo();


};

typedef std::shared_ptr<FilterAdjointGrad>  FilterAdjointGradSharedPtr;
}
}

#endif /* NEKTAR_SOLVERUTILS_FILTERS_FILTERCHECKPOINT_H */
