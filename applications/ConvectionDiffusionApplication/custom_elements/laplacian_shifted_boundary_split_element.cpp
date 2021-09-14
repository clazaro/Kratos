// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _|
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla, Franziska Wahl
//

// System includes

// External includes

// Project includes
#include "includes/checks.h"
#include "includes/define.h"
#include "includes/kratos_flags.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "includes/convection_diffusion_settings.h"
#include "utilities/geometry_utilities.h"

// Application includes
#include "custom_elements/laplacian_shifted_boundary_split_element.h"

#include "modified_shape_functions/triangle_2d_3_modified_shape_functions.h"
#include "modified_shape_functions/tetrahedra_3d_4_modified_shape_functions.h"

namespace Kratos
{

template<std::size_t TDim>
LaplacianShiftedBoundarySplitElement<TDim>::LaplacianShiftedBoundarySplitElement(
    IndexType NewId,
    GeometryType::Pointer pGeometry)
    : LaplacianElement(
        NewId,
        pGeometry)
{
}

template<std::size_t TDim>
LaplacianShiftedBoundarySplitElement<TDim>::LaplacianShiftedBoundarySplitElement(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties)
    : LaplacianElement(
        NewId,
        pGeometry,
        pProperties)
{
}

template<std::size_t TDim>
Element::Pointer LaplacianShiftedBoundarySplitElement<TDim>::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<LaplacianShiftedBoundarySplitElement>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

template<std::size_t TDim>
Element::Pointer LaplacianShiftedBoundarySplitElement<TDim>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<LaplacianShiftedBoundarySplitElement>(NewId, pGeom, pProperties);
}

template<std::size_t TDim>
LaplacianShiftedBoundarySplitElement<TDim>::~LaplacianShiftedBoundarySplitElement()
{
}

template<std::size_t TDim>
void LaplacianShiftedBoundarySplitElement<TDim>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Check if the element belongs to the intersected ones
    // Note that the BOUNDARY flag is assumed to be set in the elements which are cut by the embedded geometry
    if (Is(BOUNDARY)) {

        const auto& r_geom = GetGeometry();
        const std::size_t num_points = r_geom.size();

        //TODO only as big as number of Gauss points ?!?
        //resizing and resetting the LHS
        if(rLeftHandSideMatrix.size1() != num_points)
            rLeftHandSideMatrix.resize(num_points,num_points,false);
        noalias(rLeftHandSideMatrix) = ZeroMatrix(num_points,num_points);

        //resizing and resetting the RHS
        if(rRightHandSideVector.size() != num_points)
            rRightHandSideVector.resize(num_points,false);
        noalias(rRightHandSideVector) = ZeroVector(num_points);

        // Set up the distances vector
        //SetNodalDistancesVector(r_geom, mNodalDistances);

        //Get nodal distances and set splitting and shape functions
        InitializeGeometryData(r_geom);

        //TODO integrate over positive (??) side of the element
        /*
        // Iterate over the positive side volume integration points
        const std::size_t number_of_positive_gauss_points = data.PositiveSideWeights.size();
        for (std::size_t g = 0; g < number_of_positive_gauss_points; ++g) {
            const std::size_t gauss_pt_index = g;
            this->UpdateIntegrationPointData(data, gauss_pt_index, data.PositiveSideWeights[g], row(data.PositiveSideN, g), data.PositiveSideDNDX[g]);
            this->AddTimeIntegratedSystem(data, rLeftHandSideMatrix, rRightHandSideVector);
        }

        // Iterate over the negative side volume integration points
        const std::size_t number_of_negative_gauss_points = data.NegativeSideWeights.size();
        for (std::size_t g = 0; g < number_of_negative_gauss_points; ++g) {
            const std::size_t gauss_pt_index = g + number_of_positive_gauss_points;
            this->UpdateIntegrationPointData(data, gauss_pt_index, data.NegativeSideWeights[g], row(data.NegativeSideN, g), data.NegativeSideDNDX[g]);
            this->AddTimeIntegratedSystem(data, rLeftHandSideMatrix, rRightHandSideVector);
        }*/

    } else {
        // Add base Laplacian contribution
        BaseType::CalculateLocalSystem(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);
    }

    KRATOS_CATCH("")
}

template<std::size_t TDim>
void LaplacianShiftedBoundarySplitElement<TDim>::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    VectorType temp(0);
    CalculateLocalSystem(rLeftHandSideMatrix, temp, rCurrentProcessInfo);
}

template<std::size_t TDim>
void LaplacianShiftedBoundarySplitElement<TDim>::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    MatrixType temp(0,0);
    CalculateLocalSystem(temp, rRightHandSideVector, rCurrentProcessInfo);
}

template<std::size_t TDim>
int LaplacianShiftedBoundarySplitElement<TDim>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    // Base Laplacian element check
    return BaseType::Check(rCurrentProcessInfo);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Protected functions
///////////////////////////////////////////////////////////////////////////////////////////////////

template<std::size_t TDim>
void LaplacianShiftedBoundarySplitElement<TDim>::InitializeGeometryData(
    const GeometryType& rGeometry)
{
    // Get nodal distances
    if (mNodalDistances.size() != NumNodes) {
        mNodalDistances.resize(NumNodes);
    }
    for (std::size_t i = 0; i < NumNodes; ++i) {
        mNodalDistances[i] = rGeometry[i].FastGetSolutionStepValue(DISTANCE);
    }

    // Number and indices of positive and negative distance function values
    mNumPositiveNodes = 0;
    mNumNegativeNodes = 0;
    mPositiveIndices.clear();
    mNegativeIndices.clear();

    for (std::size_t i = 0; i < NumNodes; ++i){
        if (mNodalDistances[i] > 0.0){
            mNumPositiveNodes++;
            mPositiveIndices.push_back(i);
        } else {
            mNumNegativeNodes++;
            mNegativeIndices.push_back(i);
        }
    }

    // Get shape function calculator
    ModifiedShapeFunctions::Pointer p_calculator =
        ShiftedBoundarySplitInternals::GetContinuousShapeFunctionCalculator<TDim, NumNodes>(
            *this,
            mNodalDistances);

    // Positive side volume
    p_calculator->ComputePositiveSideShapeFunctionsAndGradientsValues(
        mPositiveSideN,
        mPositiveSideDNDX,
        mPositiveSideWeights,
        GeometryData::GI_GAUSS_2);

    // Negative side volume
    /*p_calculator->ComputeNegativeSideShapeFunctionsAndGradientsValues(
        mNegativeSideN,
        mNegativeSideDNDX,
        mNegativeSideWeights,
        GeometryData::GI_GAUSS_2);*/
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Private functions
///////////////////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////////////////////
// Helper functions for template specialization
///////////////////////////////////////////////////////////////////////////////////////////////////

namespace ShiftedBoundarySplitInternals {

template <>
ModifiedShapeFunctions::Pointer GetContinuousShapeFunctionCalculator<2, 3>(
    const Element& rElement,
    const Vector& rNodalDistances)
{
    return ModifiedShapeFunctions::Pointer(new Triangle2D3ModifiedShapeFunctions(rElement.pGetGeometry(), rNodalDistances));
}

template <>
ModifiedShapeFunctions::Pointer GetContinuousShapeFunctionCalculator<3, 4>(
    const Element& rElement,
    const Vector& rNodalDistances)
{
    return ModifiedShapeFunctions::Pointer(new Tetrahedra3D4ModifiedShapeFunctions(rElement.pGetGeometry(), rNodalDistances));
}

}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation

template class LaplacianShiftedBoundarySplitElement<2>;
template class LaplacianShiftedBoundarySplitElement<3>;

///////////////////////////////////////////////////////////////////////////////////////////////////

} // Namespace Kratos
