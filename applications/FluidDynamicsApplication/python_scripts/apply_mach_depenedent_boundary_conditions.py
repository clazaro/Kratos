import KratosMultiphysics
from KratosMultiphysics.FluidDynamicsApplication import ApplyMachDependentBoundaryConditions

def Factory(settings, Model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyMachDependentBoundaryConditions(Model, settings["Parameters"])