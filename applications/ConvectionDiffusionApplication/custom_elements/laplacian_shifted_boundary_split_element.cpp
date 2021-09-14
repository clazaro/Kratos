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

template<std::size_t TTDim>
LaplacianShiftedBoundarySplitElement<TTDim>::LaplacianShiftedBoundarySplitElement(
    IndexType NewId,
    GeometryType::Pointer pGeometry)
    : LaplacianElement(
        NewId,
        pGeometry)
{
}

template<std::size_t TTDim>
LaplacianShiftedBoundarySplitElement<TTDim>::LaplacianShiftedBoundarySplitElement(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties)
    : LaplacianElement(
        NewId,
        pGeometry,
        pProperties)
{
}

template<std::size_t TTDim>
Element::Pointer LaplacianShiftedBoundarySplitElement<TTDim>::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<LaplacianShiftedBoundarySplitElement>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

template<std::size_t TTDim>
Element::Pointer LaplacianShiftedBoundarySplitElement<TTDim>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<LaplacianShiftedBoundarySplitElement>(NewId, pGeom, pProperties);
}

template<std::size_t TTDim>
LaplacianShiftedBoundarySplitElement<TTDim>::~LaplacianShiftedBoundarySplitElement()
{
}

template<std::size_t TTDim>
void LaplacianShiftedBoundarySplitElement<TTDim>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Check if the element belongs to the intersected ones
    // Note that the BOUNDARY flag is assumed to be set in the elements which are cut by the embedded geometry
    if (Is(BOUNDARY)) {

        const auto& r_geom = GetGeometry();

        //Get nodal distances and set splitting and shape functions
        InitializeGeometryData(r_geom);

        //Calculate local system for the positive side of the element
        AddPositiveSideToSystem(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, r_geom);

    } else {
        // Add base Laplacian contribution
        BaseType::CalculateLocalSystem(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);
    }

    KRATOS_CATCH("")
}

template<std::size_t TTDim>
void LaplacianShiftedBoundarySplitElement<TTDim>::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    VectorType temp(0);
    CalculateLocalSystem(rLeftHandSideMatrix, temp, rCurrentProcessInfo);
}

template<std::size_t TTDim>
void LaplacianShiftedBoundarySplitElement<TTDim>::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    MatrixType temp(0,0);
    CalculateLocalSystem(temp, rRightHandSideVector, rCurrentProcessInfo);
}

template<std::size_t TTDim>
int LaplacianShiftedBoundarySplitElement<TTDim>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    // Base Laplacian element check
    return BaseType::Check(rCurrentProcessInfo);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Protected functions
///////////////////////////////////////////////////////////////////////////////////////////////////

template<std::size_t TTDim>
void LaplacianShiftedBoundarySplitElement<TTDim>::InitializeGeometryData(
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
        ShiftedBoundarySplitInternals::GetContinuousShapeFunctionCalculator<TTDim, NumNodes>(
            *this,
            mNodalDistances);

    // Positive side volume
    p_calculator->ComputePositiveSideShapeFunctionsAndGradientsValues(
        mPositiveSideN,
        mPositiveSideDNDX,
        mPositiveSideWeights,
        this->GetIntegrationMethod()); //GeometryData::GI_GAUSS_2

    // Negative side volume
    /*p_calculator->ComputeNegativeSideShapeFunctionsAndGradientsValues(
        mNegativeSideN,
        mNegativeSideDNDX,
        mNegativeSideWeights,
        this->GetIntegrationMethod());*/
}

template<std::size_t TTDim>
void LaplacianShiftedBoundarySplitElement<TTDim>::AddPositiveSideToSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    const GeometryType& rGeometry)
{
    auto& r_settings = *rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS];

    const Variable<double>& r_unknown_var = r_settings.GetUnknownVariable();
    const Variable<double>& r_diffusivity_var = r_settings.GetDiffusionVariable();
    const Variable<double>& r_volume_source_var = r_settings.GetVolumeSourceVariable();

    //resizing and resetting the LHS
    if(rLeftHandSideMatrix.size1() != NumNodes)
        rLeftHandSideMatrix.resize(NumNodes,NumNodes,false);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(NumNodes,NumNodes);

    //resizing and resetting the RHS
    if(rRightHandSideVector.size() != NumNodes)
        rRightHandSideVector.resize(NumNodes,false);
    noalias(rRightHandSideVector) = ZeroVector(NumNodes);

    // Get heat flux, conductivity and temp (RHS = ExtForces - K*temp) nodal vectors
    Vector heat_flux_local(NumNodes);
    Vector nodal_conductivity(NumNodes);
    Vector temp(NumNodes);
    for(std::size_t n = 0; n < NumNodes; ++n) {
        heat_flux_local[n] = rGeometry[n].FastGetSolutionStepValue(r_volume_source_var);
        nodal_conductivity[n] = rGeometry[n].FastGetSolutionStepValue(r_diffusivity_var);
        temp[n] = rGeometry[n].GetSolutionStepValue(r_unknown_var);
    }

    // Iterate over the positive side volume integration points 
    // = number of integration points * number of subdivisions on positive side of element
    const std::size_t number_of_positive_gauss_points = mPositiveSideWeights.size();
    for (std::size_t g = 0; g < number_of_positive_gauss_points; ++g) {

        const auto& N = row(mPositiveSideN, g); 
        const auto& DN_DX = mPositiveSideDNDX[g];
        const double IntToReferenceWeight = mPositiveSideWeights[g]; 

        //Calculate the local conductivity
        const double conductivity_gauss = inner_prod(N, nodal_conductivity);

        noalias(rLeftHandSideMatrix) += IntToReferenceWeight * conductivity_gauss * prod(DN_DX, trans(DN_DX)); 

        // Calculate the local RHS
        const double qgauss = inner_prod(N, heat_flux_local);

        noalias(rRightHandSideVector) += IntToReferenceWeight * qgauss * N;
    }
    
    //RHS -= K*temp
    noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,temp);  

    // Iterate over the negative side volume integration points
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
