///////////////////////////////////////////////////////////////////////////////
//
// File: StandardExtrapolate.cpp
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
// Description: Abstract base class for StandardExtrapolate.
//
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/EquationSystems/StandardExtrapolate.h>
#include <LibUtilities/Communication/Comm.h>

namespace Nektar
{
    NekDouble StandardExtrapolate::DuDt_Coeffs[3][4] = {
        { 1.0,  -1., 0.0, 0.0},
        { 2.5, -4.0, 1.5, 0.0},
        { 13./3., -9.5, 7.0, -11.0/6.0}};

    /**
     * Registers the class with the Factory.
     */
    std::string StandardExtrapolate::className = GetExtrapolateFactory().RegisterCreatorFunction(
        "Standard",
        StandardExtrapolate::create,
        "Standard");

    StandardExtrapolate::StandardExtrapolate(
        const LibUtilities::SessionReaderSharedPtr pSession,
        Array<OneD, MultiRegions::ExpListSharedPtr> pFields,
        MultiRegions::ExpListSharedPtr pPressure,
        const Array<OneD, int> pVel,
        const SolverUtils::AdvectionSharedPtr advObject)
        : Extrapolate(pSession,pFields,pPressure,pVel,advObject)
    {
    }

    StandardExtrapolate::~StandardExtrapolate()
    {
    }


    /** 
     * Function to extrapolate the new pressure boundary condition.
     * Based on the velocity field and on the advection term.
     * Acceleration term is also computed.
     * This routine is a general one for 2d and 3D application and it can be called
     * directly from velocity correction scheme. Specialisation on dimensionality is
     * redirected to the CalcNeumannPressureBCs method.
     */
    void StandardExtrapolate::v_EvaluatePressureBCs(
        const Array<OneD, const Array<OneD, NekDouble> > &fields,
        const Array<OneD, const Array<OneD, NekDouble> >  &N,
        NekDouble kinvis,
        const NekDouble time)
    {
        m_pressureCalls++;

        if(m_HBCnumber>0)
        {
            // Calculate non-linear and viscous BCs at current level
            // and put in m_pressureHBCs[0]
            CalcNeumannPressureBCs(fields,N,kinvis);
            
            // Extrapolate to n+1
            ExtrapolateArray(m_pressureHBCs);
            
            // Add (phi,Du/Dt) term to m_presureHBC
            AddDuDt(time);

            // Copy m_pressureHBCs to m_PbndExp
            CopyPressureHBCsToPbndExp();            
        }

        CalcOutflowBCs(fields, kinvis);
    }

    
    /** 
     * 
     */
    void StandardExtrapolate::v_SubSteppingTimeIntegration(
        const LibUtilities::TimeIntegrationSchemeSharedPtr & IntegrationScheme )
    {
        if ( IntegrationScheme->GetName() == "IMEX" ||
             IntegrationScheme->GetName() == "IMEXGear" )
        {
            m_intSteps = IntegrationScheme->GetOrder();
        }
        else
        {
            NEKERROR(ErrorUtil::efatal, "Integration method not suitable: "
                     "Options include IMEXGear or IMEXOrder{1,2,3,4}");
        }
    }

    /** 
     * 
     */
    void StandardExtrapolate::v_SubStepSetPressureBCs(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray, 
        NekDouble Aii_DT,
        NekDouble kinvis)
    {
    }


    /** 
     * 
     */
    void StandardExtrapolate::v_SubStepAdvance(
        int nstep,
        NekDouble time )
    {
    }


    /** 
     * 
     */
    void StandardExtrapolate::v_SubStepSaveFields(
        int nstep)
    {
    }

    /** 
     * 
     */
    void StandardExtrapolate::v_MountHOPBCs(
        int HBCdata, 
        NekDouble kinvis, 
        Array<OneD, NekDouble> &Q, 
        Array<OneD, const NekDouble> &Advection)
    {
        Vmath::Svtvp(HBCdata,-kinvis,Q,1,Advection,1,Q,1);
    }

    /**
     *    At the start, the newest value is stored in array[nlevels-1]
     *        and the previous values in the first positions
     *    At the end, the acceleration from BDF is stored in array[nlevels-1]
     *        and the storage has been updated to included the new value
     */
    void StandardExtrapolate::v_AccelerationBDF(
            Array<OneD, Array<OneD, NekDouble> > &array,
            const NekDouble time)
    {
        int nlevels = array.size();
        int nPts    = array[0].size();
        int n, cnt, i, j, totpoints;
        
        //std::cout << "bla " << nPts << endl;

        if (nPts)
        {
            // Update array
            RollOver(array);

            Array<OneD, NekDouble> accelerationTerm(nPts, 0.0);

            Array<OneD, NekDouble> IProdVnTmp;
            Array<OneD, Array<OneD, NekDouble>> velbc(m_bnd_dim);
            Array<OneD, Array<OneD, MultiRegions::ExpListSharedPtr>> VelBndExp(
                m_bnd_dim);

            for (i = 0; i < m_bnd_dim; ++i)
            {
                VelBndExp[i] = m_fields[m_velocity[i]]->GetBndCondExpansions();
            }

            // Calculate acceleration using Backward Differentiation Formula
            if (m_pressureCalls > 2)
            {
                int acc_order = min(m_pressureCalls - 2, m_intSteps);
                Vmath::Smul(nPts, DuDt_Coeffs[acc_order - 1][0], array[0], 1,
                            accelerationTerm, 1);

                for (int i = 0; i < acc_order; i++)
                {
                    Vmath::Svtvp(nPts, DuDt_Coeffs[acc_order - 1][i + 1],
                                 array[i + 1], 1, accelerationTerm, 1,
                                 accelerationTerm, 1);
                }

                //std::cout << "BDF ";
                //std::cout << accelerationTerm[nPts/2] << endl;
            }

            // Check if there is any analytical expression for acceleration

            for (n = cnt = 0; n < m_PBndConds.size(); ++n)
            {

                // High order boundary condition;
                if (m_hbcType[n] == eHBCNeumann)
                {
                    auto collectionIter = m_accelerationConditions.find(n);
                    ASSERTL1(collectionIter != m_accelerationConditions.end(),
                             "Unable to locate collection " +
                                 boost::lexical_cast<string>(n));

                    if (collectionIter != m_accelerationConditions.end())
                    {

                        for (i = 0; i < m_bnd_dim; ++i)
                        {
                            totpoints = VelBndExp[i][n]->GetTotPoints();

                            Array<OneD, NekDouble> x0(totpoints, 0.0);
                            Array<OneD, NekDouble> x1(totpoints, 0.0);
                            Array<OneD, NekDouble> x2(totpoints, 0.0);

                            VelBndExp[i][n]->GetCoords(x0, x1, x2);

                            for (int j = 0; j < totpoints; ++j)
                            {

                            //std::cout << i << " x " <<  x0[j] << endl;
                            //std::cout << i << " y " <<  x1[j] << endl;
                            }

                            velbc[i] = Array<OneD, NekDouble>(totpoints, 0.0);

                            // VelBndExp[i][n]->BwdTrans(VelBndExp[i][n]->GetCoeffs(),
                            // velbc[i]);

                            const AccelerationMapShPtr accMap =
                                (*collectionIter).second;
                            auto conditionMapIter = accMap->find(m_session->GetVariable(i));
                            ASSERTL1(conditionMapIter != accMap->end(),
                                     "Unable to locate condition map.");

                            AccelerationShPtr accPtr =
                                std::static_pointer_cast<AccelerationCondition>(
                                    (*conditionMapIter).second);

                            if (accPtr->m_expr != "")
                            {
                                LibUtilities::Equation condition =
                                    std::static_pointer_cast<
                                        AccelerationCondition>(
                                        (*conditionMapIter).second)
                                        ->m_accelerationCondition;

                                condition.Evaluate(x0, x1, x2, time, velbc[i]);
                            }

                            // analyticaccel[i]->Evaluate(x0, x1, x2, time,
                            // velbc[i]);

                            Vmath::Smul(totpoints, m_timestep, velbc[i], 1,
                                        velbc[i], 1);
                            
                            //std::cout << i << " " <<  velbc[i][totpoints/2] << endl;

                        }

                        IProdVnTmp = accelerationTerm + cnt;

                        m_PBndExp[n]->NormVectorIProductWRTBase(velbc,
                                                                IProdVnTmp);
                    }

                    cnt += m_PBndExp[n]->GetNcoeffs();
                }
                else if (m_hbcType[n] == eConvectiveOBC ||
                         (m_hbcType[n] == eHBCNeumann))
                {
                    // skip over convective OBC
                    cnt += m_PBndExp[n]->GetNcoeffs();
                }
            }

            //std::cout << "final ";
            //std::cout << accelerationTerm[nPts/2] << endl;

            array[nlevels - 1] = accelerationTerm;
            
            
        }
    }
}
