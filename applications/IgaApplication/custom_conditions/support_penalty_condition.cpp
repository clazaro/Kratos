//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ricky Aristio
//                   Tobias Teschemacher
//

// System includes

// External includes

// Project includes
#include "custom_conditions/support_penalty_condition.h"


namespace Kratos
{
    void SupportPenaltyCondition::CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
    )
    {
        KRATOS_TRY
        const double penalty = GetProperties()[PENALTY_FACTOR];

        const auto& r_geometry = GetGeometry();
        const SizeType number_of_nodes = r_geometry.size();
        const SizeType mat_size = r_geometry.WorkingSpaceDimension() * number_of_nodes;

        // Memory allocation
        if (CalculateStiffnessMatrixFlag) {
            if (rLeftHandSideMatrix.size1() != mat_size) {
                rLeftHandSideMatrix.resize(mat_size, mat_size, false);
            }
            noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size);
        }
        if (CalculateResidualVectorFlag) {
            if (rRightHandSideVector.size() != mat_size) {
                rRightHandSideVector.resize(mat_size, false);
            }
            rRightHandSideVector = ZeroVector(mat_size);
        }


        if (Has(DISPLACEMENT))
        {
            // Integration
            const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints();

            // initial determinant of jacobian 
            Vector determinat_jacobian_vector_initial(integration_points.size());
            DeterminantOfJacobianInitial(r_geometry, determinat_jacobian_vector_initial);

            for (IndexType point_number = 0; point_number < integration_points.size(); point_number++)
            {
                const Matrix& N = r_geometry.ShapeFunctionsValues();

                //FOR DISPLACEMENTS
                Matrix H = ZeroMatrix(3, mat_size);
                for (IndexType i = 0; i < number_of_nodes; i++)
                {
                    IndexType index = 3 * i;
                    //if (Is(IgaFlags::FIX_DISPLACEMENT_X))
                        H(0, index) = N(point_number, i);
                    //if (Is(IgaFlags::FIX_DISPLACEMENT_Y))
                        H(1, index + 1) = N(point_number, i);
                    //if (Is(IgaFlags::FIX_DISPLACEMENT_Z))
                        H(2, index + 2) = N(point_number, i);
                }

                // Differential area
                const double integration_weight = integration_points[point_number].Weight();
                const double determinat_jacobian = determinat_jacobian_vector_initial[point_number];

                // Assembly
                if (CalculateStiffnessMatrixFlag) {
                    noalias(rLeftHandSideMatrix) += prod(trans(H), H)
                        * penalty* integration_weight * determinat_jacobian;
                }
                if (CalculateResidualVectorFlag) {

                    const array_1d<double, 3>& displacement = this->GetValue(DISPLACEMENT);

                    Vector u(mat_size);
                    for (IndexType i = 0; i < number_of_nodes; i++)
                    {
                        const array_1d<double, 3> disp = r_geometry[i].FastGetSolutionStepValue(DISPLACEMENT);
                        IndexType index = 3 * i;
                        u[index]     = (disp[0] - displacement[0]);
                        u[index + 1] = (disp[1] - displacement[1]);
                        u[index + 2] = (disp[2] - displacement[2]);
                    }

                    noalias(rRightHandSideVector) -= prod(prod(trans(H), H), u)
                        * penalty * integration_weight * determinat_jacobian;
                }
            }
        }
        KRATOS_CATCH("")
    }

    void SupportPenaltyCondition::DeterminantOfJacobianInitial(
        const GeometryType& rGeometry,
        Vector& rDeterminantOfJacobian)
    {
        const IndexType nb_integration_points = rGeometry.IntegrationPointsNumber();
        if (rDeterminantOfJacobian.size() != nb_integration_points) {
            rDeterminantOfJacobian.resize(nb_integration_points, false);
        }

        const SizeType working_space_dimension = rGeometry.WorkingSpaceDimension();
        const SizeType local_space_dimension = rGeometry.LocalSpaceDimension();
        const SizeType nb_nodes = rGeometry.PointsNumber();

        Matrix J = ZeroMatrix(working_space_dimension, local_space_dimension);
        for (IndexType pnt = 0; pnt < nb_integration_points; pnt++)
        {
            const Matrix& r_DN_De = rGeometry.ShapeFunctionsLocalGradients()[pnt];
            J.clear();
            for (IndexType i = 0; i < nb_nodes; ++i) {
                const array_1d<double, 3>& r_coordinates = rGeometry[i].GetInitialPosition();
                for (IndexType k = 0; k < working_space_dimension; ++k) {
                    const double value = r_coordinates[k];
                    for (IndexType m = 0; m < local_space_dimension; ++m) {
                        J(k, m) += value * r_DN_De(i, m);
                    }
                }
            }

            //Compute the tangent and  the normal to the boundary vector
            array_1d<double, 3> local_tangent;
            GetGeometry().Calculate(LOCAL_TANGENT, local_tangent);

            array_1d<double, 3> a_1 = column(J, 0);
            array_1d<double, 3> a_2 = column(J, 1);

            rDeterminantOfJacobian[pnt] = norm_2(a_1 * local_tangent[0] + a_2 * local_tangent[1]);
        }
    }

    int SupportPenaltyCondition::Check(const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;
        KRATOS_ERROR_IF_NOT(GetProperties().Has(PENALTY_FACTOR)) << "No penalty factor (PENALTY_FACTOR) defined in property of SupportPenaltyPointDiscreteCondition" << std::endl;
        //KRATOS_WATCH(GetGeometry())
        return 0;
        KRATOS_CATCH("");
    }   

    void SupportPenaltyCondition::EquationIdVector(
        EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo
    )
    {
        const auto& r_geometry = GetGeometry();
        const SizeType number_of_nodes = r_geometry.size();

        if (rResult.size() != 3 * number_of_nodes)
            rResult.resize(3 * number_of_nodes, false);

        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const IndexType index = i * 3;
            const auto& r_node = r_geometry[i];
            rResult[index]     = r_node.GetDof(DISPLACEMENT_X).EquationId();
            rResult[index + 1] = r_node.GetDof(DISPLACEMENT_Y).EquationId();
            rResult[index + 2] = r_node.GetDof(DISPLACEMENT_Z).EquationId();
        }
    }

    void SupportPenaltyCondition::GetDofList(
        DofsVectorType& rElementalDofList,
        ProcessInfo& rCurrentProcessInfo
    )
    {
        const auto& r_geometry = GetGeometry();
        const SizeType number_of_nodes = r_geometry.size();

        rElementalDofList.resize(0);
        rElementalDofList.reserve(3 * number_of_nodes);

        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const auto& r_node = r_geometry[i];
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Y));
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Z));
        }
    };

} // Namespace Kratos