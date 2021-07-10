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

#if !defined(KRATOS_RANS_VARIABLE_DATA_TRANSFER_PROCESS_H_INCLUDED)
#define KRATOS_RANS_VARIABLE_DATA_TRANSFER_PROCESS_H_INCLUDED

// System includes
#include <string>
#include <vector>
#include <tuple>

// External includes

// Project includes
#include "containers/model.h"
#include "containers/variable.h"
#include "includes/model_part.h"
#include "processes/process.h"

// Application incldues
#include "custom_processes/rans_formulation_process.h"

namespace Kratos
{
///@addtogroup RANSApplication
///@{

///@name Kratos Classes
///@{

/**
 * @brief Apply a specific flag for nodes and conditions
 *
 * This process apply a given flag to nodes of the modelpart.
 * Then, if preferred, applies to all conditions in the given model, which has
 * the given flag applied to all the nodes in the specific condition.
 *
 */

class KRATOS_API(RANS_APPLICATION) RansVariableDataTransferProcess : public RansFormulationProcess
{
public:
    ///@name Type Definitions
    ///@{

    using NodeType = ModelPart::NodeType;

    using CopyVariableDataListItem = std::tuple<const std::string, const bool, const std::string, const bool>;

    /// Pointer definition of RansVariableDataTransferProcess
    KRATOS_CLASS_POINTER_DEFINITION(RansVariableDataTransferProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    RansVariableDataTransferProcess(
        Model& rModel,
        Parameters rParameters);

    RansVariableDataTransferProcess(
        Model& rModel,
        const std::string& rSourceModelPartName,
        const std::string& rDestinationModelPartName,
        const std::vector<std::string>& rCopyExecutionPoints,
        const std::vector<CopyVariableDataListItem>& rCopyVariableDataList,
        const int EchoLevel);

    RansVariableDataTransferProcess(
        Model& rSourceModel,
        Model& rDestinationModel,
        const std::string& rSourceModelPartName,
        const std::string& rDestinationModelPartName,
        const std::vector<std::string>& rCopyExecutionPoints,
        const std::vector<CopyVariableDataListItem>& rCopyVariableDataList,
        const int EchoLevel);

    /// Destructor.
    ~RansVariableDataTransferProcess() override = default;

    /// Assignment operator.
    RansVariableDataTransferProcess& operator=(RansVariableDataTransferProcess const& rOther) = delete;

    /// Copy constructor.
    RansVariableDataTransferProcess(RansVariableDataTransferProcess const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    int Check() override;

    void ExecuteInitialize() override;

    void ExecuteInitializeSolutionStep() override;

    void ExecuteBeforeCouplingSolveStep() override;

    void Execute() override;

    void ExecuteAfterCouplingSolveStep() override;

    void ExecuteFinalizeSolutionStep() override;

    void ExecuteFinalize() override;

    const Parameters GetDefaultParameters() const override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override;

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override;

    ///@}

private:
    ///@name Private Classes
    ///@{

    enum ExecutionPoint
    {
        INITIALIZE = 0,
        INITIALIZE_SOLUTION_STEP = 1,
        BEFORE_COUPLING_SOLVE_STEP = 2,
        EXECUTE = 3,
        AFTER_COUPLING_SOLVE_STEP = 4,
        FINALIZE_SOLUTION_STEP = 5,
        FINALIZE = 6
    };

    template<class TVariableType>
    class CopyVariableData
    {
    public:
        ///@name Life Cycle
        ///@{

        explicit CopyVariableData(
            const TVariableType& rSourceVariable,
            const bool IsSourceVariableHistorical,
            const TVariableType& rDestinationVariable,
            const bool IsDestinationVariableHistorical);

        ///@}
        ///@name Public Operations
        ///@{

        void CopyData(const NodeType& rSourceNode, NodeType& rDestinationNode) const;

        std::string Info() const;

        ///@}

    private:
        ///@name Private Member Variables
        ///@{

        const TVariableType& mrSourceVariable;
        const bool mIsSourceVariableHistorical;
        const TVariableType& mrDestinationVariable;
        const bool mIsDestinationVariableHistorical;

        void (CopyVariableData::*mCopyMethod)(const NodeType& rSourceNode, NodeType& rDestinationNode) const;

        void CopyHistoricalToHistorical(const NodeType& rSourceNode, NodeType& rDestinationNode) const;

        void CopyHistoricalToNonHistorical(const NodeType& rSourceNode, NodeType& rDestinationNode) const;

        void CopyNonHistoricalToHistorical(const NodeType& rSourceNode, NodeType& rDestinationNode) const;

        void CopyNonHistoricalToNonHistorical(const NodeType& rSourceNode, NodeType& rDestinationNode) const;

        ///@}
    };

    ///@}
    ///@name Private Member Variables
    ///@{

    Model& mrSourceModel;
    Model& mrDestinationModel;
    int mEchoLevel;

    std::string mSourceModelPartName;
    std::string mDestinationModelPartName;

    std::vector<ExecutionPoint> mExecutionPointsList;

    std::vector<CopyVariableData<Variable<double>>> mCopyDoubleVariableDataList;
    std::vector<CopyVariableData<Variable<array_1d<double, 3>>>> mCopyArray3DVariableDataList;
    std::vector<CopyVariableData<Variable<Vector>>> mCopyVectorVariableDataList;
    std::vector<CopyVariableData<Variable<Matrix>>> mCopyMatrixVariableDataList;

    ///@}
    ///@name Private Operations
    ///@{

    void ExecuteCopy();

    void UpdateCopyExecutionPointsList(const std::vector<std::string>& rCopyExecutionPointsList);

    void UpdateCopyVariableDataList(const std::vector<CopyVariableDataListItem>& rCopyVariableDataList);

    ///@}

}; // Class RansVariableDataTransferProcess

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const RansVariableDataTransferProcess& rThis);

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_RANS_VARIABLE_DATA_TRANSFER_PROCESS_H_INCLUDED defined