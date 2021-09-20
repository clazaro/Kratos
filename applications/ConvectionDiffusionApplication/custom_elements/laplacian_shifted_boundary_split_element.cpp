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
#include "utilities/element_size_calculator.h"

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

    SplitElementData data;
    data.Initialize(*this);

    // Check if the element belongs to the intersected ones
    // Note that the BOUNDARY flag is assumed to be set in the elements which are cut by the embedded geometry
    if (data.IsSplit()) {

        // Get nodal distances and set splitting and shape functions
        InitializeGeometryData(data);

        // Resizing and resetting the LHS
        if(rLeftHandSideMatrix.size1() != NumNodes)
            rLeftHandSideMatrix.resize(NumNodes,NumNodes,false);
        noalias(rLeftHandSideMatrix) = ZeroMatrix(NumNodes,NumNodes);

        // Resizing and resetting the RHS
        if(rRightHandSideVector.size() != NumNodes)
            rRightHandSideVector.resize(NumNodes,false);
        noalias(rRightHandSideVector) = ZeroVector(NumNodes);

        // Calculate and add local system for the positive side of the element
        AddPositiveElementSide(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, data);
        
        // Calculate and add boundary terms
        AddPositiveBoundaryTerms(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, data);

        //Calculate and add Nitsche terms for weak imposition of boundary condition
        //AddNormalPenaltyContribution()
        //AddNormalSymmetricCounterpartContribution()
        //AddTangentialPenaltyContribution()
        //AddTangentialSymmetricCounterpartContribution()

    } else {
        // Add base Laplacian contribution (standard Galerkin)
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
    SplitElementData& rData)
{
    // Get shape function calculator
    ModifiedShapeFunctions::Pointer p_calculator =
        LaplacianSplitInternals::GetContinuousShapeFunctionCalculator<TTDim, NumNodes>(
            *this,
            rData.NodalDistances);

    // Positive side volume
    p_calculator->ComputePositiveSideShapeFunctionsAndGradientsValues(
        rData.PositiveSideN,
        rData.PositiveSideDNDX,
        rData.PositiveSideWeights,
        this->GetIntegrationMethod()); //GeometryData::GI_GAUSS_2

    // Negative side volume
    /*p_calculator->ComputeNegativeSideShapeFunctionsAndGradientsValues(
        rData.NegativeSideN,
        rData.NegativeSideDNDX,
        rData.NegativeSideWeights,
        this->GetIntegrationMethod());*/

    // Positive side interface
    p_calculator->ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues(
        rData.PositiveInterfaceN,
        rData.PositiveInterfaceDNDX,
        rData.PositiveInterfaceWeights,
        this->GetIntegrationMethod());

    // Positive side interface normals
    p_calculator->ComputePositiveSideInterfaceAreaNormals(
        rData.PositiveInterfaceUnitNormals,
        this->GetIntegrationMethod());

    // Normalize the normals
    // Note: we calculate h here (and we don't use the value in rData.ElementSize)
    // because rData.ElementSize might still be uninitialized: some data classes define it at the Gauss point.
    double h = ElementSizeCalculator<TTDim,NumNodes>::MinimumElementSize(this->GetGeometry());
    const double tolerance = std::pow(1e-3 * h, TTDim-1);
    this->NormalizeInterfaceNormals(rData.PositiveInterfaceUnitNormals, tolerance);
}

template<std::size_t TTDim>
void LaplacianShiftedBoundarySplitElement<TTDim>::NormalizeInterfaceNormals(
    typename SplitElementData::InterfaceNormalsType& rNormals,
    double Tolerance) const
{
    for (std::size_t i = 0; i < rNormals.size(); ++i) {
        double norm = norm_2(rNormals[i]);
        rNormals[i] /= std::max(norm,Tolerance);
    }
}

template<std::size_t TTDim>
void LaplacianShiftedBoundarySplitElement<TTDim>::AddPositiveElementSide(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    const SplitElementData& rData)
{
    const auto& r_geom = GetGeometry();
    auto& r_settings = *rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS];

    const Variable<double>& r_unknown_var = r_settings.GetUnknownVariable();
    const Variable<double>& r_diffusivity_var = r_settings.GetDiffusionVariable();
    const Variable<double>& r_volume_source_var = r_settings.GetVolumeSourceVariable();

    // Get heat flux, conductivity and temp (RHS = ExtForces - K*temp) nodal vectors
    Vector heat_flux_local(NumNodes);
    Vector nodal_conductivity(NumNodes);
    Vector temp(NumNodes);
    for(std::size_t n = 0; n < NumNodes; ++n) {
        heat_flux_local[n] = r_geom[n].FastGetSolutionStepValue(r_volume_source_var);
        nodal_conductivity[n] = r_geom[n].FastGetSolutionStepValue(r_diffusivity_var);
        temp[n] = r_geom[n].GetSolutionStepValue(r_unknown_var);
    }

    // Iterate over the positive side volume integration points 
    // = number of integration points * number of subdivisions on positive side of element
    const std::size_t number_of_positive_gauss_points = rData.PositiveSideWeights.size();
    for (std::size_t g = 0; g < number_of_positive_gauss_points; ++g) {

        const auto& N = row(rData.PositiveSideN, g); 
        const auto& DN_DX = rData.PositiveSideDNDX[g];
        const double weight = rData.PositiveSideWeights[g]; 

        //Calculate the local conductivity
        const double conductivity_gauss = inner_prod(N, nodal_conductivity);

        noalias(rLeftHandSideMatrix) += weight * conductivity_gauss * prod(DN_DX, trans(DN_DX)); 

        // Calculate the local RHS (external source)
        const double qgauss = inner_prod(N, heat_flux_local);

        noalias(rRightHandSideVector) += weight * qgauss * N;
    }
    
    //RHS -= K*temp
    noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,temp);  

    // Iterate over the negative side volume integration points
}

template<std::size_t TTDim>
void LaplacianShiftedBoundarySplitElement<TTDim>::AddPositiveBoundaryTerms(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    const SplitElementData& rData)
{
    const auto& r_geom = GetGeometry();
    auto& r_settings = *rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS];

    const Variable<double>& r_unknown_var = r_settings.GetUnknownVariable();
    //const Variable<double>& r_diffusivity_var = r_settings.GetDiffusionVariable();
    //const Variable<double>& r_volume_source_var = r_settings.GetVolumeSourceVariable();

    // Get heat flux, conductivity and temp (RHS = ExtForces - K*temp) nodal vectors
    //Vector heat_flux_local(NumNodes);
    //Vector nodal_conductivity(NumNodes);
    Vector temp(NumNodes);
    for(std::size_t n = 0; n < NumNodes; ++n) {
        //heat_flux_local[n] = r_geom[n].FastGetSolutionStepValue(r_volume_source_var);
        //nodal_conductivity[n] = r_geom[n].FastGetSolutionStepValue(r_diffusivity_var);
        temp[n] = r_geom[n].GetSolutionStepValue(r_unknown_var);
    }

    // Iterate over the positive side volume integration points 
    // = number of integration points * number of subdivisions on positive side of element
    const std::size_t number_of_positive_gauss_points = rData.PositiveInterfaceWeights.size();
    for (std::size_t g = 0; g < number_of_positive_gauss_points; ++g) {

        const auto& N = row(rData.PositiveInterfaceN, g); 
        const auto& DN_DX = rData.PositiveInterfaceDNDX[g];
        const double weight = rData.PositiveInterfaceWeights[g]; 
        const auto& unit_normal = rData.PositiveInterfaceUnitNormals[g];

        //Calculate the local conductivity
        //const double conductivity_gauss = inner_prod(N, nodal_conductivity);

        //const BoundedMatrix<double, NumNodes, TTDim> aux_N_projected = prod(N, trans(unit_normal));        

        //noalias(aux_LHS) += weight * prod(aux_N_projected, trans(DN_DX)); 

        // ALTERNATIVE
        for (std::size_t i = 0; i < NumNodes; ++i) {
            for (std::size_t j = 0; j < NumNodes; ++j) {
                for (std::size_t d = 0; d < TTDim; ++d) {
                    const double aux = weight * N(i) * unit_normal(d) * DN_DX(j,d);
                    rLeftHandSideMatrix(i, j) += aux;
                    rRightHandSideVector(i) -= aux * temp(j);
                }
            }
        }
    }
    
    //noalias(rLeftHandSideMatrix) += aux_LHS;
    //RHS -= K*temp
    //noalias(rRightHandSideVector) -= prod(aux_LHS,temp);  

    // Iterate over the negative side volume integration points
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Private functions
///////////////////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////////////////////
// Helper functions for template specialization
///////////////////////////////////////////////////////////////////////////////////////////////////

namespace LaplacianSplitInternals {

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

template<std::size_t TTDim>
void SplitElementData<TTDim>::Initialize(
    const Element& rElement)
{
    const auto& r_geom = rElement.GetGeometry();;

    // Get nodal distances
    if (NodalDistances.size() != NumNodes) {
        NodalDistances.resize(NumNodes);
    }
    for (std::size_t i = 0; i < NumNodes; ++i) {
        NodalDistances[i] = r_geom[i].FastGetSolutionStepValue(DISTANCE);
    }

    // Number and indices of positive and negative distance function values
    NumPositiveNodes = 0;
    NumNegativeNodes = 0;
    PositiveIndices.clear();

    for (std::size_t i = 0; i < NumNodes; ++i){
        if (NodalDistances[i] > 0.0){
            NumPositiveNodes++;
            PositiveIndices.push_back(i);
        } else {
            NumNegativeNodes++;
        }
    }
}

template<std::size_t TTDim>
bool SplitElementData<TTDim>::IsSplit()
{
    return (NumPositiveNodes > 0) && (NumNegativeNodes > 0);
}

}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation

template class LaplacianShiftedBoundarySplitElement<2>;
template class LaplacianShiftedBoundarySplitElement<3>;

///////////////////////////////////////////////////////////////////////////////////////////////////

} // Namespace Kratos
