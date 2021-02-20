//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes

// External includes

// Project includes
#include "includes/checks.h"
#include "includes/element.h"
#include "includes/properties.h"
#include "utilities/adjoint_extensions.h"
#include "utilities/geometrical_sensitivity_utility.h"
#include "utilities/indirect_scalar.h"

// Application includes
#include "fluid_dynamics_application_variables.h"

#include "custom_elements/data_containers/qs_vms/qs_vms_adjoint_element_data.h"

// Include base h
#include "fluid_adjoint_element.h"

namespace Kratos
{
template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::ThisExtensions::ThisExtensions(Element* pElement)
    : mpElement{pElement}
{
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::ThisExtensions::GetFirstDerivativesVector(
    std::size_t NodeId,
    std::vector<IndirectScalar<double>>& rVector,
    std::size_t Step)
{
    auto& r_node = mpElement->GetGeometry()[NodeId];
    rVector.resize(TDim + 1);

    const auto& dofs_list = TAdjointElementData::GetDofVariablesList();

    for (unsigned int i = 0; i < TDim; ++i) {
        rVector[i] = MakeIndirectScalar(r_node, (*dofs_list[i]).GetTimeDerivative(), Step);
    }

    rVector[TDim] = IndirectScalar<double>{}; // pressure
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::ThisExtensions::GetSecondDerivativesVector(
    std::size_t NodeId,
    std::vector<IndirectScalar<double>>& rVector,
    std::size_t Step)
{
    auto& r_node = mpElement->GetGeometry()[NodeId];
    rVector.resize(TDim + 1);

    const auto& dofs_list = TAdjointElementData::GetDofVariablesList();

    for (unsigned int i = 0; i < TDim; ++i) {
        rVector[i] = MakeIndirectScalar(
            r_node, (*dofs_list[i]).GetTimeDerivative().GetTimeDerivative(), Step);
    }

    rVector[TDim] = IndirectScalar<double>{}; // pressure
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::ThisExtensions::GetAuxiliaryVector(
    std::size_t NodeId,
    std::vector<IndirectScalar<double>>& rVector,
    std::size_t Step)
{
    auto& r_node = mpElement->GetGeometry()[NodeId];
    rVector.resize(TDim + 1);

    const auto& dofs_list = TAdjointElementData::GetDofVariablesList();

    for (unsigned int i = 0; i < TDim; ++i) {
        rVector[i] = MakeIndirectScalar(
            r_node, (*dofs_list[i]).GetTimeDerivative().GetTimeDerivative().GetTimeDerivative(),
            Step);
    }

    rVector[TDim] = IndirectScalar<double>{}; // pressure
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::ThisExtensions::GetFirstDerivativesVariables(
    std::vector<VariableData const*>& rVariables) const
{
    rVariables.resize(1);
    rVariables[0] = &ADJOINT_FLUID_VECTOR_2;
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::ThisExtensions::GetSecondDerivativesVariables(
    std::vector<VariableData const*>& rVariables) const
{
    rVariables.resize(1);
    rVariables[0] = &ADJOINT_FLUID_VECTOR_3;
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::ThisExtensions::GetAuxiliaryVariables(
    std::vector<VariableData const*>& rVariables) const
{
    rVariables.resize(1);
    rVariables[0] = &AUX_ADJOINT_FLUID_VECTOR_1;
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::FluidAdjointElement(IndexType NewId)
    : BaseType(NewId)
{
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::FluidAdjointElement(
    IndexType NewId,
    GeometryType::Pointer pGeometry)
    : BaseType(NewId, pGeometry)
{
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::FluidAdjointElement(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties)
    : BaseType(NewId, pGeometry, pProperties)
{
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::~FluidAdjointElement()
{
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
Element::Pointer FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<FluidAdjointElement>(
        NewId, Element::GetGeometry().Create(ThisNodes), pProperties);
    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
Element::Pointer FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<FluidAdjointElement>(NewId, pGeom, pProperties);
    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
Element::Pointer FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::Clone(
    IndexType NewId,
    NodesArrayType const& ThisNodes) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<FluidAdjointElement>(
        NewId, Element::GetGeometry().Create(ThisNodes), Element::pGetProperties());
    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
int FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    TAdjointElementData::Check(*this, rCurrentProcessInfo);
    return 0;
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::EquationIdVector(
    EquationIdVectorType& rElementalEquationIdList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    if (rElementalEquationIdList.size() != TElementLocalSize) {
        rElementalEquationIdList.resize(TElementLocalSize, false);
    }

    const auto& r_variables_list = TAdjointElementData::GetDofVariablesList();

    IndexType local_index = 0;
    for (IndexType i = 0; i < TNumNodes; ++i) {
        const auto& r_node = this->GetGeometry()[i];
        for (const auto p_variable : r_variables_list) {
            rElementalEquationIdList[local_index++] =
                r_node.GetDof(*p_variable).EquationId();
        }
    }
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    if (rElementalDofList.size() != TElementLocalSize) {
        rElementalDofList.resize(TElementLocalSize);
    }

    const auto& r_variables_list = TAdjointElementData::GetDofVariablesList();

    IndexType local_index = 0;
    for (IndexType i = 0; i < TNumNodes; ++i) {
        const auto& r_node = this->GetGeometry()[i];
        for (const auto p_variable : r_variables_list) {
            rElementalDofList[local_index++] = r_node.pGetDof(*p_variable);
        }
    }
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::GetValuesVector(
    VectorType& rValues,
    int Step) const
{
    if (rValues.size() != TElementLocalSize) {
        rValues.resize(TElementLocalSize, false);
    }

    const auto& r_geometry = this->GetGeometry();
    IndexType local_index = 0;
    for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
        const auto& r_node = r_geometry[i_node];
        const auto& r_velocity =
            r_node.FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_1, Step);
        for (IndexType d = 0; d < TDim; ++d) {
            rValues[local_index++] = r_velocity[d];
        }
        rValues[local_index++] =
            r_node.FastGetSolutionStepValue(ADJOINT_FLUID_SCALAR_1, Step);
    }
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::GetFirstDerivativesVector(
    VectorType& rValues,
    int Step) const
{
    if (rValues.size() != TElementLocalSize) {
        rValues.resize(TElementLocalSize, false);
    }
    rValues.clear();
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::GetSecondDerivativesVector(
    VectorType& rValues,
    int Step) const
{
    if (rValues.size() != TElementLocalSize) {
        rValues.resize(TElementLocalSize);
    }

    const auto& r_geometry = this->GetGeometry();
    IndexType local_index = 0;
    for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
        const auto& r_acceleration =
            r_geometry[i_node].FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_3, Step);
        for (IndexType d = 0; d < TDim; ++d) {
            rValues[local_index++] = r_acceleration[d];
        }
        rValues[local_index++] = 0.0; // pressure dof
    }
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    // If we are restarting, the constitutive law will be already defined
    if (mpConstitutiveLaw == nullptr) {
        const Properties& r_properties = this->GetProperties();

        KRATOS_ERROR_IF_NOT(r_properties.Has(CONSTITUTIVE_LAW))
            << "In initialization of Element " << this->Info()
            << ": No CONSTITUTIVE_LAW defined for property "
            << r_properties.Id() << "." << std::endl;

        mpConstitutiveLaw = r_properties[CONSTITUTIVE_LAW]->Clone();

        const GeometryType& r_geometry = this->GetGeometry();
        const auto& r_shape_functions =
            r_geometry.ShapeFunctionsValues(GeometryData::GI_GAUSS_1);

        mpConstitutiveLaw->InitializeMaterial(r_properties, r_geometry,
                                              row(r_shape_functions, 0));
    }

    this->SetValue(ADJOINT_EXTENSIONS, Kratos::make_shared<ThisExtensions>(this));

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Resize and initialize output
    if (rLeftHandSideMatrix.size1() != 0)
        rLeftHandSideMatrix.resize(0, 0, false);

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rLeftHandSideMatrix.size1() != TElementLocalSize || rLeftHandSideMatrix.size2() != TElementLocalSize) {
        rLeftHandSideMatrix.resize(TElementLocalSize, TElementLocalSize, false);
    }

    rLeftHandSideMatrix.clear();
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_ERROR << "FluidAdjointElement::"
                    "CalculateRightHandSide method is not implemented.";

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::CalculateLocalVelocityContribution(
    MatrixType &rDampMatrix,
    VectorType &rRightHandSideVector,
    const ProcessInfo &rCurrentProcessInfo)
{
    if (rRightHandSideVector.size() != TElementLocalSize)
        rRightHandSideVector.resize(TElementLocalSize, false);

    rRightHandSideVector.clear();

    AddFluidResidualsContributions(rRightHandSideVector, rCurrentProcessInfo);
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::CalculateFirstDerivativesLHS(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rLeftHandSideMatrix.size1() != TElementLocalSize ||
        rLeftHandSideMatrix.size2() != TElementLocalSize) {
        rLeftHandSideMatrix.resize(TElementLocalSize, TElementLocalSize, false);
    }

    rLeftHandSideMatrix.clear();

    AddFluidFirstDerivatives(rLeftHandSideMatrix, rCurrentProcessInfo);
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::CalculateSecondDerivativesLHS(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rLeftHandSideMatrix.size1() != TElementLocalSize ||
        rLeftHandSideMatrix.size2() != TElementLocalSize) {
        rLeftHandSideMatrix.resize(TElementLocalSize, TElementLocalSize, false);
    }

    rLeftHandSideMatrix.clear();

    AddFluidSecondDerivatives(rLeftHandSideMatrix, rCurrentProcessInfo);
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::CalculateMassMatrix(
    MatrixType& rMassMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    rMassMatrix.resize(0, 0);
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::CalculateDampingMatrix(
    MatrixType& rDampingMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_ERROR << "FluidAdjointElement::"
                    "CalculateDampingMatrix method is not implemented.";

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::CalculateSensitivityMatrix(
    const Variable<array_1d<double, 3>>& rSensitivityVariable,
    Matrix& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rSensitivityVariable == SHAPE_SENSITIVITY) {
        if (rOutput.size1() != TCoordsLocalSize || rOutput.size2() != TElementLocalSize) {
            rOutput.resize(TCoordsLocalSize, TElementLocalSize, false);
        }

        rOutput.clear();
        AddFluidShapeDerivatives(rOutput, rCurrentProcessInfo);
    } else {
        KRATOS_ERROR << "Sensitivity variable " << rSensitivityVariable
                     << " not supported." << std::endl;
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
std::string FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::Info() const
{
    std::stringstream buffer;
    buffer << "FluidAdjointElement #" << Element::Id();
    return buffer.str();
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "FluidAdjointElement #" << Element::Id();
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::PrintData(std::ostream& rOStream) const
{
    Element::pGetGeometry()->PrintData(rOStream);
}

///@}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::AddFluidResidualsContributions(
    VectorType& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto& integration_method = TAdjointElementData::TResidualsDerivatives::GetIntegrationMethod();

    Vector Ws;
    Matrix Ns;
    ShapeFunctionDerivativesArrayType dNdXs;
    this->CalculateGeometryData(Ws, Ns, dNdXs, integration_method);

    typename TAdjointElementData::Primal::Data element_data(*this, *mpConstitutiveLaw, rCurrentProcessInfo);
    typename TAdjointElementData::Primal::ResidualsContributions residual_contributions(element_data);

    VectorF residual = ZeroVector(TElementLocalSize);

    for (IndexType g = 0; g < Ws.size(); ++g) {
        const Vector& N = row(Ns, g);
        const Matrix& dNdX = dNdXs[g];
        const double W = Ws[g];

        element_data.CalculateGaussPointData(W, N, dNdX);
        residual_contributions.AddGaussPointResidualsContributions(residual, W, N, dNdX);
    }

    AssembleSubVectorToVector(rOutput, residual);

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::AddFluidFirstDerivatives(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto& integration_method = TAdjointElementData::TResidualsDerivatives::GetIntegrationMethod();

    Vector Ws;
    Matrix Ns;
    ShapeFunctionDerivativesArrayType dNdXs;
    this->CalculateGeometryData(Ws, Ns, dNdXs, integration_method);

    using Derivatives = typename TAdjointElementData::StateDerivatives::FirstDerivatives;

    typename Derivatives::Data     element_data(*this, *mpConstitutiveLaw, rCurrentProcessInfo);
    typename Derivatives::Velocity velocity_derivative(element_data);
    typename Derivatives::Pressure pressure_derivative(element_data);

    BoundedVector<double, TElementLocalSize> residual;
    BoundedMatrix<double, TNumNodes, TDim> dNdXDerivative = ZeroMatrix(TNumNodes, TDim);

    for (IndexType g = 0; g < Ws.size(); ++g) {
        const double W = Ws[g];
        const Vector& N = row(Ns, g);
        const Matrix& dNdX = dNdXs[g];

        element_data.CalculateGaussPointData(W, N, dNdX);

        IndexType row = 0;
        for (IndexType c = 0; c < TNumNodes; ++c) {
            for (IndexType k = 0; k < TDim; ++k) {
                velocity_derivative.CalculateGaussPointResidualsDerivativeContributions(residual, c, k, W, N, dNdX, 0.0, 0.0, dNdXDerivative);
                AssembleSubVectorToMatrix(rLeftHandSideMatrix, row++, residual);
            }

            pressure_derivative.CalculateGaussPointResidualsDerivativeContributions(residual, c, 0, W, N, dNdX, 0.0, 0.0, dNdXDerivative);
            AssembleSubVectorToMatrix(rLeftHandSideMatrix, row++, residual);
        }
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::AddFluidSecondDerivatives(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto& integration_method = TAdjointElementData::TResidualsDerivatives::GetIntegrationMethod();

    Vector Ws;
    Matrix Ns;
    ShapeFunctionDerivativesArrayType dNdXs;
    this->CalculateGeometryData(Ws, Ns, dNdXs, integration_method);

    using Derivatives = typename TAdjointElementData::StateDerivatives::SecondDerivatives;

    typename Derivatives::Data         element_data(*this, *mpConstitutiveLaw, rCurrentProcessInfo);
    typename Derivatives::Acceleration acceleration_derivative(element_data);

    VectorF residual;

    for (IndexType g = 0; g < Ws.size(); ++g) {
        const double W = Ws[g];
        const Vector& N = row(Ns, g);
        const Matrix& dNdX = dNdXs[g];

        element_data.CalculateGaussPointData(W, N, dNdX);

        IndexType row = 0;
        for (IndexType c = 0; c < TNumNodes; ++c) {
            for (IndexType k = 0; k < TDim; ++k) {
                acceleration_derivative.CalculateGaussPointResidualsDerivativeContributions(residual, c, k, W, N, dNdX);
                AssembleSubVectorToMatrix(rLeftHandSideMatrix, row++, residual);
            }

            // skip derivatives w.r.t. pressure time derivative
            ++row;
        }
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::AddFluidShapeDerivatives(
    Matrix& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto& integration_method = TAdjointElementData::TResidualsDerivatives::GetIntegrationMethod();

    Vector Ws;
    Matrix Ns;
    ShapeFunctionDerivativesArrayType dNdXs;
    this->CalculateGeometryData(Ws, Ns, dNdXs, integration_method);

    using Derivatives = typename TAdjointElementData::SensitivityDerivatives;

    typename Derivatives::Data  element_data(*this, *mpConstitutiveLaw, rCurrentProcessInfo);
    typename Derivatives::Shape derivative(element_data);

    VectorF residual;

    for (IndexType g = 0; g < Ws.size(); ++g) {
        const Vector& N = row(Ns, g);
        const Matrix& dNdX = dNdXs[g];
        const double W = Ws[g];

        element_data.CalculateGaussPointData(W, N, dNdX);

        Geometry<Point>::JacobiansType J;
        this->GetGeometry().Jacobian(J, integration_method);
        const auto& DN_De = this->GetGeometry().ShapeFunctionsLocalGradients(integration_method);

        GeometricalSensitivityUtility::ShapeFunctionsGradientType dNdX_derivative;
        const Matrix& rJ = J[g];
        const Matrix& rDN_De = DN_De[g];
        const double inv_detJ = 1.0 / MathUtils<double>::DetMat(rJ);
        GeometricalSensitivityUtility geom_sensitivity(rJ, rDN_De);

        ShapeParameter deriv;
        IndexType row = 0;
        for (deriv.NodeIndex = 0; deriv.NodeIndex < TNumNodes; ++deriv.NodeIndex) {
            for (deriv.Direction = 0; deriv.Direction < TDim; ++deriv.Direction) {
                double detJ_derivative;
                geom_sensitivity.CalculateSensitivity(deriv, detJ_derivative, dNdX_derivative);
                const double W_derivative = detJ_derivative * inv_detJ * W;

                derivative.CalculateGaussPointResidualsDerivativeContributions(residual, deriv.NodeIndex, deriv.Direction, W, N, dNdX, W_derivative, detJ_derivative, dNdX_derivative);
                AssembleSubVectorToMatrix(rOutput, row++, residual);
            }
        }
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::CalculateGeometryData(
    Vector& rGaussWeights,
    Matrix& rNContainer,
    ShapeFunctionDerivativesArrayType& rDN_DX,
    const GeometryData::IntegrationMethod& rIntegrationMethod) const
{
    const auto& r_geometry = this->GetGeometry();
    const IndexType number_of_gauss_points =
        r_geometry.IntegrationPointsNumber(rIntegrationMethod);

    Vector DetJ;
    r_geometry.ShapeFunctionsIntegrationPointsGradients(rDN_DX, DetJ, rIntegrationMethod);

    if (rNContainer.size1() != number_of_gauss_points || rNContainer.size2() != TNumNodes) {
        rNContainer.resize(number_of_gauss_points, TNumNodes, false);
    }
    rNContainer = r_geometry.ShapeFunctionsValues(rIntegrationMethod);

    const auto& IntegrationPoints = r_geometry.IntegrationPoints(rIntegrationMethod);

    if (rGaussWeights.size() != number_of_gauss_points) {
        rGaussWeights.resize(number_of_gauss_points, false);
    }

    for (IndexType g = 0; g < number_of_gauss_points; ++g) {
        rGaussWeights[g] = DetJ[g] * IntegrationPoints[g].Weight();
    }
}

// template instantiations
template class FluidAdjointElement<2, 3, QSVMSAdjointElementData<2, 3>>;
template class FluidAdjointElement<2, 4, QSVMSAdjointElementData<2, 4>>;
template class FluidAdjointElement<3, 4, QSVMSAdjointElementData<3, 4>>;
template class FluidAdjointElement<3, 8, QSVMSAdjointElementData<3, 8>>;

} // namespace Kratos