#> Main script
# Add Python bindings directory to PATH
import sys, os
import math
import numpy as np
#from fuzzywuzzy import fuzz
#from scipy.optimize import minimize


# Intialise OpenCMISS
from opencmiss.iron import iron
#iron.OutputSetOn("Testing")


# Set problem parameters
growthModel = 1 # Type of growth. 1 - volumetric; 2 - stress based; 3 - strain based.
isotropic = False # True if the problem is isotropic, False if the problem is anisotropic
homogeneous = False # True if the growth rates are homogeneous, False if the problem is heterogeneous
useFibres = True # True if fibres are used for anisotropic problems
heterogeneousFibres = False # True if the fibre angles vary in space, False if not.
fixBottomRing = True # True if the bottom ring of nodes is fixed, False if not
fixTopRing = False # True if the top ring of nodes is fixed, False if not


# Tube geometry
length = 3.0 # The length of the tube
innerRadius = 2.0 # The inner radius of the tube
outerRadius = 5.0 # The outer radius of the tube


numberOfLengthElements=1 # Number of elements along the length of the tube
numberOfCircumfrentialElementsPerQuarter=1 # Number of elements in the circumfrential direction in one quarter of the tube
numberOfWallElements=1 # Number of elements through the wall of the tube

iron.DiagnosticsSetOn(iron.DiagnosticTypes.IN,[1,2,3,4,5],"Diagnostics",["DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE"])

# Hydrostatic pressure
pInit = -8.0 # The initial hydrostatic pressure

# Fibre angle
#fibreAngle = math.pi/2.0 # The fibre angle wrt the for anisotropic fibres
fibreAngle = 0.0 # The fibre angle wrt the for anisotropic fibres

# Times
startTime = 0.0 # The start time for the growth simulation
stopTime1 = 10.0 # The stop time for the growth simulation
timeIncrement = 1.0 # The time increment for the growth simulation

# Number of Gauss points used
numberOfGaussXi = 3

numberOfCircumfrentialElements = 4*numberOfCircumfrentialElementsPerQuarter
numberOfLengthNodes = numberOfLengthElements+1
numberOfCircumfrentialNodes = numberOfCircumfrentialElements
numberOfWallNodes = numberOfWallElements+1

coordinateSystemUserNumber = 1
regionUserNumber = 1
tricubicHermiteBasisUserNumber = 1
trilinearLagrangeBasisUserNumber = 2
meshUserNumber = 1
decompositionUserNumber = 1
geometricFieldUserNumber = 1
originalGeometricFieldUserNumber = 13
DELUDELNGeometricFieldUserNumber = 14
fibreFieldUserNumber = 2
dependentFieldUserNumber = 3
equationsSetUserNumber = 1
equationsSetFieldUserNumber = 5
growthCellMLUserNumber = 1
growthCellMLModelsFieldUserNumber = 6
growthCellMLStateFieldUserNumber = 7
growthCellMLParametersFieldUserNumber = 8
constitutiveCellMLUserNumber = 2
constitutiveCellMLModelsFieldUserNumber = 9
constitutiveCellMLParametersFieldUserNumber = 10
constitutiveCellMLIntermediateFieldUserNumber = 11
problemUserNumber = 1

# Get the number of computational nodes and this computational node number
numberOfComputationalNodes = iron.ComputationalNumberOfNodesGet()
computationalNodeNumber = iron.ComputationalNodeNumberGet()


# Create a 3D rectangular cartesian coordinate system
coordinateSystem = iron.CoordinateSystem()
coordinateSystem.CreateStart(coordinateSystemUserNumber)
# Set the number of dimensions to 3
coordinateSystem.DimensionSet(3)
# Finish the creation of the coordinate system
coordinateSystem.CreateFinish()

# Create a region and assign the coordinate system to the region
region = iron.Region()
region.CreateStart(regionUserNumber,iron.WorldRegion)
region.LabelSet("HeartTubeRegion")
# Set the regions coordinate system to the 3D RC coordinate system that we have created
region.coordinateSystem = coordinateSystem
# Finish the creation of the region
region.CreateFinish()

# Define basis
# Start the creation of a tricubic Hermite basis function
tricubicHermiteBasis = iron.Basis()
tricubicHermiteBasis.CreateStart(tricubicHermiteBasisUserNumber)
tricubicHermiteBasis.type = iron.BasisTypes.LAGRANGE_HERMITE_TP
tricubicHermiteBasis.numberOfXi = 3
tricubicHermiteBasis.interpolationXi = [iron.BasisInterpolationSpecifications.CUBIC_HERMITE]*3
tricubicHermiteBasis.quadratureNumberOfGaussXi = [numberOfGaussXi]*3
tricubicHermiteBasis.CreateFinish()
# Start the creation of a trilinear Hermite basis function
trilinearLagrangeBasis = iron.Basis()
trilinearLagrangeBasis.CreateStart(trilinearLagrangeBasisUserNumber)
trilinearLagrangeBasis.type = iron.BasisTypes.LAGRANGE_HERMITE_TP
trilinearLagrangeBasis.numberOfXi = 3
trilinearLagrangeBasis.interpolationXi = [iron.BasisInterpolationSpecifications.LINEAR_LAGRANGE]*3
trilinearLagrangeBasis.quadratureNumberOfGaussXi = [numberOfGaussXi]*3
trilinearLagrangeBasis.CreateFinish()

# Start the creation of a manually generated mesh in the region
numberOfNodes = numberOfCircumfrentialElements*(numberOfLengthElements+1)*(numberOfWallElements+1)
numberOfElements = numberOfCircumfrentialElements*numberOfLengthElements*numberOfWallElements

# Define nodes for the mesh
nodes = iron.Nodes()
nodes.CreateStart(region,numberOfNodes)
nodes.CreateFinish()
mesh = iron.Mesh()

# Create the mesh. The mesh will have two components - 1. tricubic Hermite elements; 2. trilinear Lagrange elements
mesh.CreateStart(meshUserNumber,region,3)
mesh.NumberOfComponentsSet(2)
mesh.NumberOfElementsSet(numberOfElements)

tricubicHermiteElements = iron.MeshElements()
tricubicHermiteElements.CreateStart(mesh,1,tricubicHermiteBasis)
trilinearLagrangeElements = iron.MeshElements()
trilinearLagrangeElements.CreateStart(mesh,2,trilinearLagrangeBasis)

elementNumber = 0
for wallElementIdx in range(1,numberOfWallElements+1):
    for lengthElementIdx in range(1,numberOfLengthElements+1):
        for circumfrentialElementIdx in range(1,numberOfCircumfrentialElements+1):
            elementNumber = elementNumber + 1
            localNode1 = circumfrentialElementIdx + (lengthElementIdx-1)*numberOfCircumfrentialNodes + \
                (wallElementIdx-1)*numberOfCircumfrentialNodes*numberOfLengthNodes
            if circumfrentialElementIdx == numberOfCircumfrentialElements:
                localNode2 = 1 + (lengthElementIdx-1)*numberOfCircumfrentialNodes + \
                    (wallElementIdx-1)*numberOfCircumfrentialNodes*numberOfLengthNodes
            else:
                localNode2 = localNode1 + 1
            localNode3 = localNode1 + numberOfCircumfrentialNodes
            localNode4 = localNode2 + numberOfCircumfrentialNodes
            localNode5 = localNode1 + numberOfCircumfrentialNodes*numberOfLengthNodes
            localNode6 = localNode2 + numberOfCircumfrentialNodes*numberOfLengthNodes
            localNode7 = localNode3 + numberOfCircumfrentialNodes*numberOfLengthNodes
            localNode8 = localNode4 + numberOfCircumfrentialNodes*numberOfLengthNodes
            localNodes = [localNode1,localNode2,localNode3,localNode4,localNode5,localNode6,localNode7,localNode8]
            tricubicHermiteElements.NodesSet(elementNumber,localNodes)
            trilinearLagrangeElements.NodesSet(elementNumber,localNodes)

tricubicHermiteElements.CreateFinish()
trilinearLagrangeElements.CreateFinish()
mesh.CreateFinish()

# Create a decomposition for the mesh
decomposition = iron.Decomposition()
decomposition.CreateStart(decompositionUserNumber,mesh)
decomposition.type = iron.DecompositionTypes.CALCULATED
decomposition.numberOfDomains = numberOfComputationalNodes
decomposition.CreateFinish()

# Create a field for the geometry
geometricField = iron.Field()
geometricField.CreateStart(geometricFieldUserNumber,region)
geometricField.MeshDecompositionSet(decomposition)
geometricField.TypeSet(iron.FieldTypes.GEOMETRIC)
geometricField.VariableLabelSet(iron.FieldVariableTypes.U,"Geometry")
geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,1)
geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,2,1)
geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,3,1)
geometricField.fieldScalingType = iron.FieldScalingTypes.ARITHMETIC_MEAN
geometricField.CreateFinish()

# defining another Geometric  field with the same set-up to be used for any new iteration of optimization 
originalGeometricField = iron.Field()
originalGeometricField.CreateStart(originalGeometricFieldUserNumber,region)
originalGeometricField.MeshDecompositionSet(decomposition)
originalGeometricField.TypeSet(iron.FieldTypes.GEOMETRIC)
originalGeometricField.VariableLabelSet(iron.FieldVariableTypes.U,"OriginalGeometry")
originalGeometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,1)
originalGeometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,2,1)
originalGeometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,3,1)
originalGeometricField.fieldScalingType = iron.FieldScalingTypes.ARITHMETIC_MEAN
originalGeometricField.CreateFinish()

# defining another Geometric  field with the same set-up to be used for any new iteration of optimization 
DELUDELNGeometricField = iron.Field()
DELUDELNGeometricField.CreateStart(DELUDELNGeometricFieldUserNumber,region)
DELUDELNGeometricField.MeshDecompositionSet(decomposition)
DELUDELNGeometricField.TypeSet(iron.FieldTypes.GEOMETRIC)
DELUDELNGeometricField.VariableLabelSet(iron.FieldVariableTypes.U,"DELUDELNGeometry")
DELUDELNGeometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,1)
DELUDELNGeometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,2,1)
DELUDELNGeometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,3,1)
DELUDELNGeometricField.fieldScalingType = iron.FieldScalingTypes.ARITHMETIC_MEAN
DELUDELNGeometricField.CreateFinish()


# Create the geometric field
for wallNodeIdx in range(1,numberOfWallNodes+1):
    for lengthNodeIdx in range(1,numberOfLengthNodes+1):
        for circumfrentialNodeIdx in range(1,numberOfCircumfrentialNodes+1):
            nodeNumber = circumfrentialNodeIdx + (lengthNodeIdx-1)*numberOfCircumfrentialNodes + \
                (wallNodeIdx-1)*numberOfCircumfrentialNodes*numberOfLengthNodes
            radius = innerRadius + (outerRadius - innerRadius)*float(wallNodeIdx-1)/float(numberOfWallNodes)
            theta = float(circumfrentialNodeIdx-1)/float(numberOfCircumfrentialNodes)*2.0*math.pi
            x = radius*math.cos(theta)
            y = radius*math.sin(theta)
            xtangent = -math.sin(theta)
            ytangent = math.cos(theta)
            xnormal = math.cos(theta)
            ynormal = math.sin(theta)
            z = float(lengthNodeIdx-1)/float(numberOfLengthElements)*length
            geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,1,x)
            geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,2,y)
            geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,3,z)
            geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,1,xtangent)
            geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,2,ytangent)
            geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,3,0.0)
            geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,1,0.0)
            geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,2,0.0)
            geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,3,1.0)
            geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3,nodeNumber,1,xnormal)
            geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3,nodeNumber,2,ynormal)
            geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3,nodeNumber,3,0.0)
            originalGeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,1,x)
            originalGeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,2,y)
            originalGeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,3,z)
            originalGeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,1,xtangent)
            originalGeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,2,ytangent)
            originalGeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,3,0.0)
            originalGeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,1,0.0)
            originalGeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,2,0.0)
            originalGeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,3,1.0)
            originalGeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3,nodeNumber,1,xnormal)
            originalGeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3,nodeNumber,2,ynormal)
            originalGeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3,nodeNumber,3,0.0)
            DELUDELNGeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,1,0.0)
            DELUDELNGeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,2,0.0)
            DELUDELNGeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,3,0.0)
            DELUDELNGeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,1,0.0)
            DELUDELNGeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,2,0.0)
            DELUDELNGeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,3,0.0)
            DELUDELNGeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,1,0.0)
            DELUDELNGeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,2,0.0)
            DELUDELNGeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,3,0.0)
            DELUDELNGeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3,nodeNumber,1,0.0)
            DELUDELNGeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3,nodeNumber,2,0.0)
            DELUDELNGeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3,nodeNumber,3,0.0)


# Update the geometric field
geometricField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
geometricField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
# Update the original geometric field
originalGeometricField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
originalGeometricField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
# Update the DELUDELN geometric field
DELUDELNGeometricField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
DELUDELNGeometricField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)


if useFibres:
    # Create a fibre field and attach it to the geometric field
    fibreField = iron.Field()
    fibreField.CreateStart(fibreFieldUserNumber,region)
    fibreField.TypeSet(iron.FieldTypes.FIBRE)
    fibreField.MeshDecompositionSet(decomposition)
    fibreField.GeometricFieldSet(geometricField)
    fibreField.VariableLabelSet(iron.FieldVariableTypes.U,"Fibre")
    fibreField.NumberOfComponentsSet(iron.FieldVariableTypes.U,3)
    fibreField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,2)   # using the mesh component number 2 for the fibers ... 
    fibreField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,2,2)   # using the mesh component number 2 for the fibers ...
    fibreField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,3,2)   # using the mesh component number 2 for the fibers ... 
    fibreField.CreateFinish()
    #Initialise the fibre field
    for wallNodeIdx in range(1,numberOfWallNodes+1):
        for lengthNodeIdx in range(1,numberOfLengthNodes+1):
            for circumfrentialNodeIdx in range(1,numberOfCircumfrentialNodes+1):
                nodeNumber = circumfrentialNodeIdx + (lengthNodeIdx-1)*numberOfCircumfrentialNodes + \
                    (wallNodeIdx-1)*numberOfCircumfrentialNodes*numberOfLengthNodes
                # Set the fibre angle
                if heterogeneousFibres == True:
            	    theta = float(circumfrentialNodeIdx-1)/float(numberOfCircumfrentialNodes)*2.0*math.pi
		    angle = fibreAngle*math.sin(theta)
                else:
		    angle = fibreAngle
                fibreField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                        1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,1,angle)
                fibreField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                        1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,2,0.0)
                fibreField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                        1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,3,0.0)
    # Update the fibre field
    fibreField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
    fibreField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)

# Create the dependent field
dependentField = iron.Field()
dependentField.CreateStart(dependentFieldUserNumber,region)
dependentField.TypeSet(iron.FieldTypes.GEOMETRIC_GENERAL)
# Set the decomposition
dependentField.MeshDecompositionSet(decomposition)
# Set the geometric field
dependentField.GeometricFieldSet(geometricField) 
dependentField.DependentTypeSet(iron.FieldDependentTypes.DEPENDENT)
# Set the field variables for displacement, traction, strain, stress and growth
dependentField.NumberOfVariablesSet(5)
dependentField.VariableTypesSet([iron.FieldVariableTypes.U,iron.FieldVariableTypes.DELUDELN,
                                 iron.FieldVariableTypes.U1,iron.FieldVariableTypes.U2,iron.FieldVariableTypes.U3])
dependentField.VariableLabelSet(iron.FieldVariableTypes.U,"Dependent")
dependentField.VariableLabelSet(iron.FieldVariableTypes.DELUDELN,"del U/del n")
dependentField.VariableLabelSet(iron.FieldVariableTypes.U1,"Strain")
dependentField.VariableLabelSet(iron.FieldVariableTypes.U2,"Stress")
dependentField.VariableLabelSet(iron.FieldVariableTypes.U3,"Growth")
dependentField.NumberOfComponentsSet(iron.FieldVariableTypes.U,4)
dependentField.NumberOfComponentsSet(iron.FieldVariableTypes.DELUDELN,4)
dependentField.NumberOfComponentsSet(iron.FieldVariableTypes.U1,6)
dependentField.NumberOfComponentsSet(iron.FieldVariableTypes.U2,6)
dependentField.NumberOfComponentsSet(iron.FieldVariableTypes.U3,3)
# Set the hydrostatic pressure to use tri-linear Lagrange elements
dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,4,2)
dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN,4,2)
# Set the strain, stress and growth variables to be Gauss point based.
dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U1,1,
                                         iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U1,2,
                                         iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U1,3,
                                         iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U1,4,
                                         iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U1,5,
                                         iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U1,6,
                                         iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U2,1,
                                         iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U2,2,
                                         iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U2,3,
                                         iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U2,4,
                                         iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U2,5,
                                         iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U2,6,
                                         iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U3,1,
                                         iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U3,2,
                                         iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U3,3,
                                         iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
# Set the field scaling
dependentField.fieldScalingType = iron.FieldScalingTypes.ARITHMETIC_MEAN
# Finish creating the field
dependentField.CreateFinish()

# Initialise dependent field from undeformed geometry
iron.Field.ParametersToFieldParametersComponentCopy(
    geometricField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,
    dependentField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1)
iron.Field.ParametersToFieldParametersComponentCopy(
    geometricField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2,
    dependentField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2)
iron.Field.ParametersToFieldParametersComponentCopy(
    geometricField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,3,
    dependentField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,3)
# Initialise the hydrostatic pressure
iron.Field.ComponentValuesInitialiseDP(dependentField,iron.FieldVariableTypes.U,
                                       iron.FieldParameterSetTypes.VALUES,4,pInit)

# Create the equations_set
equationsSetField = iron.Field()
equationsSet = iron.EquationsSet()
# Specify a finite elasticity equations set with the growth and constitutive law in CellML
equationsSetSpecification = [iron.EquationsSetClasses.ELASTICITY,
    iron.EquationsSetTypes.FINITE_ELASTICITY,
    iron.EquationsSetSubtypes.CONSTIT_AND_GROWTH_LAW_IN_CELLML]
if useFibres:
    equationsSet.CreateStart(equationsSetUserNumber,region,fibreField,
                             equationsSetSpecification,equationsSetFieldUserNumber,
                             equationsSetField)
else:
    equationsSet.CreateStart(equationsSetUserNumber,region,geometricField,
                             equationsSetSpecification,equationsSetFieldUserNumber,
                             equationsSetField)
equationsSet.CreateFinish()

# Set up the equation set dependent field
equationsSet.DependentCreateStart(dependentFieldUserNumber,dependentField)
equationsSet.DependentCreateFinish()

# Create equations
equations = iron.Equations()
equationsSet.EquationsCreateStart(equations)
# Use sparse equations
equations.sparsityType = iron.EquationsSparsityTypes.SPARSE
# Do not output any equations information
equations.outputType = iron.EquationsOutputTypes.NONE
# Finish creating the equations
equationsSet.EquationsCreateFinish()

# Set up the growth CellML model
growthCellML = iron.CellML()
growthCellML.CreateStart(growthCellMLUserNumber,region)
growthCellMLIdx = growthCellML.ModelImport("simplegrowth.cellml")
growthCellML.VariableSetAsKnown(growthCellMLIdx,"Main/fibrerate")
growthCellML.VariableSetAsKnown(growthCellMLIdx,"Main/sheetrate")
growthCellML.VariableSetAsKnown(growthCellMLIdx,"Main/normalrate")
growthCellML.CreateFinish()

# Create CellML <--> OpenCMISS field maps
growthCellML.FieldMapsCreateStart()
growthCellML.CreateCellMLToFieldMap(growthCellMLIdx,"Main/lambda1",iron.FieldParameterSetTypes.VALUES,
	dependentField,iron.FieldVariableTypes.U3,1,iron.FieldParameterSetTypes.VALUES)
growthCellML.CreateCellMLToFieldMap(growthCellMLIdx,"Main/lambda2",iron.FieldParameterSetTypes.VALUES,
	dependentField,iron.FieldVariableTypes.U3,2,iron.FieldParameterSetTypes.VALUES)
growthCellML.CreateCellMLToFieldMap(growthCellMLIdx,"Main/lambda3",iron.FieldParameterSetTypes.VALUES,
        dependentField,iron.FieldVariableTypes.U3,3,iron.FieldParameterSetTypes.VALUES)
growthCellML.FieldMapsCreateFinish()

# Create the CELL models field
growthCellMLModelsField = iron.Field()
growthCellML.ModelsFieldCreateStart(growthCellMLModelsFieldUserNumber,growthCellMLModelsField)
growthCellMLModelsField.VariableLabelSet(iron.FieldVariableTypes.U,"GrowthModelMap")
growthCellML.ModelsFieldCreateFinish()

# Create the CELL parameters field
growthCellMLParametersField = iron.Field()
growthCellML.ParametersFieldCreateStart(growthCellMLParametersFieldUserNumber,growthCellMLParametersField)
growthCellMLParametersField.VariableLabelSet(iron.FieldVariableTypes.U,"GrowthParameters")
growthCellML.ParametersFieldCreateFinish()


# Create the CELL state field
growthCellMLStateField = iron.Field()
growthCellML.StateFieldCreateStart(growthCellMLStateFieldUserNumber,growthCellMLStateField)
growthCellMLStateField.VariableLabelSet(iron.FieldVariableTypes.U,"GrowthState")
growthCellML.StateFieldCreateFinish()

# Create the CellML environment for the consitutative law
constitutiveCellML = iron.CellML()
constitutiveCellML.CreateStart(constitutiveCellMLUserNumber,region)
constitutiveCellMLIdx = constitutiveCellML.ModelImport("mooneyrivlin.cellml")
# Flag the CellML variables that OpenCMISS will supply
constitutiveCellML.VariableSetAsKnown(constitutiveCellMLIdx,"equations/E11")
constitutiveCellML.VariableSetAsKnown(constitutiveCellMLIdx,"equations/E12")
constitutiveCellML.VariableSetAsKnown(constitutiveCellMLIdx,"equations/E13")
constitutiveCellML.VariableSetAsKnown(constitutiveCellMLIdx,"equations/E22")
constitutiveCellML.VariableSetAsKnown(constitutiveCellMLIdx,"equations/E23")
constitutiveCellML.VariableSetAsKnown(constitutiveCellMLIdx,"equations/E33")
#constitutiveCellML.VariableSetAsKnown(constitutiveCellMLIdx,"equations/c1")
#constitutiveCellML.VariableSetAsKnown(constitutiveCellMLIdx,"equations/c2")
# Flag the CellML variables that OpenCMISS will obtain
constitutiveCellML.VariableSetAsWanted(constitutiveCellMLIdx,"equations/Tdev11")
constitutiveCellML.VariableSetAsWanted(constitutiveCellMLIdx,"equations/Tdev12")
constitutiveCellML.VariableSetAsWanted(constitutiveCellMLIdx,"equations/Tdev13")
constitutiveCellML.VariableSetAsWanted(constitutiveCellMLIdx,"equations/Tdev22")
constitutiveCellML.VariableSetAsWanted(constitutiveCellMLIdx,"equations/Tdev23")
constitutiveCellML.VariableSetAsWanted(constitutiveCellMLIdx,"equations/Tdev33")
constitutiveCellML.CreateFinish()

# Create CellML <--> OpenCMISS field maps
constitutiveCellML.FieldMapsCreateStart()
constitutiveCellML.CreateFieldToCellMLMap(dependentField,iron.FieldVariableTypes.U1,1,iron.FieldParameterSetTypes.VALUES,
    constitutiveCellMLIdx,"equations/E11",iron.FieldParameterSetTypes.VALUES)
constitutiveCellML.CreateFieldToCellMLMap(dependentField,iron.FieldVariableTypes.U1,2,iron.FieldParameterSetTypes.VALUES,
    constitutiveCellMLIdx,"equations/E12",iron.FieldParameterSetTypes.VALUES)
constitutiveCellML.CreateFieldToCellMLMap(dependentField,iron.FieldVariableTypes.U1,3,iron.FieldParameterSetTypes.VALUES,
    constitutiveCellMLIdx,"equations/E13",iron.FieldParameterSetTypes.VALUES)
constitutiveCellML.CreateFieldToCellMLMap(dependentField,iron.FieldVariableTypes.U1,4,iron.FieldParameterSetTypes.VALUES,
    constitutiveCellMLIdx,"equations/E22",iron.FieldParameterSetTypes.VALUES)
constitutiveCellML.CreateFieldToCellMLMap(dependentField,iron.FieldVariableTypes.U1,5,iron.FieldParameterSetTypes.VALUES,
    constitutiveCellMLIdx,"equations/E23",iron.FieldParameterSetTypes.VALUES)
constitutiveCellML.CreateFieldToCellMLMap(dependentField,iron.FieldVariableTypes.U1,6,iron.FieldParameterSetTypes.VALUES,
    constitutiveCellMLIdx,"equations/E33",iron.FieldParameterSetTypes.VALUES)
constitutiveCellML.CreateCellMLToFieldMap(constitutiveCellMLIdx,"equations/Tdev11",iron.FieldParameterSetTypes.VALUES,
    dependentField,iron.FieldVariableTypes.U2,1,iron.FieldParameterSetTypes.VALUES)
constitutiveCellML.CreateCellMLToFieldMap(constitutiveCellMLIdx,"equations/Tdev12",iron.FieldParameterSetTypes.VALUES,
    dependentField,iron.FieldVariableTypes.U2,2,iron.FieldParameterSetTypes.VALUES)
constitutiveCellML.CreateCellMLToFieldMap(constitutiveCellMLIdx,"equations/Tdev13",iron.FieldParameterSetTypes.VALUES,
    dependentField,iron.FieldVariableTypes.U2,3,iron.FieldParameterSetTypes.VALUES)
constitutiveCellML.CreateCellMLToFieldMap(constitutiveCellMLIdx,"equations/Tdev22",iron.FieldParameterSetTypes.VALUES,
    dependentField,iron.FieldVariableTypes.U2,4,iron.FieldParameterSetTypes.VALUES)
constitutiveCellML.CreateCellMLToFieldMap(constitutiveCellMLIdx,"equations/Tdev23",iron.FieldParameterSetTypes.VALUES,
    dependentField,iron.FieldVariableTypes.U2,5,iron.FieldParameterSetTypes.VALUES)
constitutiveCellML.CreateCellMLToFieldMap(constitutiveCellMLIdx,"equations/Tdev33",iron.FieldParameterSetTypes.VALUES,
    dependentField,iron.FieldVariableTypes.U2,6,iron.FieldParameterSetTypes.VALUES)
constitutiveCellML.FieldMapsCreateFinish()

# Create the CELL models field
constitutiveCellMLModelsField = iron.Field()
constitutiveCellML.ModelsFieldCreateStart(constitutiveCellMLModelsFieldUserNumber,constitutiveCellMLModelsField)
constitutiveCellMLModelsField.VariableLabelSet(iron.FieldVariableTypes.U,"ConstitutiveModelMap")
constitutiveCellML.ModelsFieldCreateFinish()

# Create the CELL parameters field
constitutiveCellMLParametersField = iron.Field()
constitutiveCellML.ParametersFieldCreateStart(constitutiveCellMLParametersFieldUserNumber,constitutiveCellMLParametersField)
constitutiveCellMLParametersField.VariableLabelSet(iron.FieldVariableTypes.U,"ConstitutiveParameters")
constitutiveCellML.ParametersFieldCreateFinish()

# Create the CELL intermediate field
constitutiveCellMLIntermediateField = iron.Field()
constitutiveCellML.IntermediateFieldCreateStart(constitutiveCellMLIntermediateFieldUserNumber,constitutiveCellMLIntermediateField)
constitutiveCellMLIntermediateField.VariableLabelSet(iron.FieldVariableTypes.U,"ConstitutiveIntermediate")
constitutiveCellML.IntermediateFieldCreateFinish()

# Define the problem
problem = iron.Problem()
problemSpecification = [iron.ProblemClasses.ELASTICITY,iron.ProblemTypes.FINITE_ELASTICITY,
        iron.ProblemSubtypes.FINITE_ELASTICITY_WITH_GROWTH_CELLML]
problem.CreateStart(problemUserNumber,problemSpecification)
problem.CreateFinish()

# Create control loops
timeLoop = iron.ControlLoop()
problem.ControlLoopCreateStart()
problem.ControlLoopGet([iron.ControlLoopIdentifiers.NODE],timeLoop)
problem.ControlLoopCreateFinish()

# Create problem solvers
odeIntegrationSolver = iron.Solver()
nonlinearSolver = iron.Solver()
linearSolver = iron.Solver()
cellMLEvaluationSolver = iron.Solver()
problem.SolversCreateStart()
problem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,odeIntegrationSolver)
problem.SolverGet([iron.ControlLoopIdentifiers.NODE],2,nonlinearSolver)
# nonlinearSolver.outputType = iron.SolverOutputTypes.MONITOR
nonlinearSolver.NewtonJacobianCalculationTypeSet(iron.JacobianCalculationTypes.FD)
nonlinearSolver.NewtonCellMLSolverGet(cellMLEvaluationSolver)
nonlinearSolver.NewtonLinearSolverGet(linearSolver)
linearSolver.linearType = iron.LinearSolverTypes.DIRECT
problem.SolversCreateFinish()

# Create nonlinear equations and add equations set to solver equations
nonlinearEquations = iron.SolverEquations()
problem.SolverEquationsCreateStart()
nonlinearSolver.SolverEquationsGet(nonlinearEquations)
nonlinearEquations.sparsityType = iron.SolverEquationsSparsityTypes.SPARSE
nonlinearEquationsSetIndex = nonlinearEquations.EquationsSetAdd(equationsSet)
problem.SolverEquationsCreateFinish()

# Create CellML equations and add growth and constitutive equations to the solvers
growthEquations = iron.CellMLEquations()
constitutiveEquations = iron.CellMLEquations()
problem.CellMLEquationsCreateStart()
odeIntegrationSolver.CellMLEquationsGet(growthEquations)
growthEquationsIndex = growthEquations.CellMLAdd(growthCellML)
cellMLEvaluationSolver.CellMLEquationsGet(constitutiveEquations)
constitutiveEquationsIndex = constitutiveEquations.CellMLAdd(constitutiveCellML)
problem.CellMLEquationsCreateFinish()

# Prescribe boundary conditions (absolute nodal parameters)
boundaryConditions = iron.BoundaryConditions()
nonlinearEquations.BoundaryConditionsCreateStart(boundaryConditions)

fixBottomRing = True
fixTopRing = False

for lengthNodeIdx in range(1,3):
    if (lengthNodeIdx == 1 and fixBottomRing) or (lengthNodeIdx == 2 and fixTopRing):
        for wallNodeIdx in range(1,numberOfWallNodes+1):
            for circumfrentialNodeIdx in range(1,numberOfCircumfrentialNodes+1):
                nodeNumber = circumfrentialNodeIdx + (lengthNodeIdx-1)*(numberOfLengthNodes-1)*numberOfCircumfrentialNodes + \
                             (wallNodeIdx-1)*numberOfCircumfrentialNodes*numberOfLengthNodes
                #print "z and all derivs ",nodeNumber,  
                # Fix z direction
                boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,
                                           1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,3,
                                           iron.BoundaryConditionsTypes.FIXED,0.0)
                # Fix S1 (circumfrential) direction derivatives
                boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,
                                           1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,1,
                                           iron.BoundaryConditionsTypes.FIXED,0.0)
                boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,
                                           1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,2,
                                           iron.BoundaryConditionsTypes.FIXED,0.0)
                boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,
                                           1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,3,
                                           iron.BoundaryConditionsTypes.FIXED,0.0)
                # Fix S2 (length) direction derivatives
                boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,
                                           1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,1,
                                           iron.BoundaryConditionsTypes.FIXED,0.0)
                boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,
                                           1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,2,
                                           iron.BoundaryConditionsTypes.FIXED,0.0)
                boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,
                                           1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,3,
                                           iron.BoundaryConditionsTypes.FIXED,0.0)
                # Fix S3 (wall) direction derivatives
                boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,
                                           1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3,nodeNumber,1,
                                           iron.BoundaryConditionsTypes.FIXED,0.0)
                boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,
                                           1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3,nodeNumber,2,
                                           iron.BoundaryConditionsTypes.FIXED,0.0)
                boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,
                                           1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3,nodeNumber,3,
                                           iron.BoundaryConditionsTypes.FIXED,0.0)
        #Set symmetry conditions on the ring to prevent rotation
        nodeNumber = 1 + (lengthNodeIdx-1)*(numberOfLengthNodes-1)*numberOfCircumfrentialNodes + \
                     (wallNodeIdx-1)*numberOfCircumfrentialNodes*numberOfLengthNodes
        #print "y direction  ",nodeNumber
        boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,
                                   1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,2,
                                   iron.BoundaryConditionsTypes.FIXED,0.0)
        nodeNumber = nodeNumber + numberOfCircumfrentialElementsPerQuarter
        #print "x direction  ",nodeNumber
        boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,
                                   1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,1,
                                   iron.BoundaryConditionsTypes.FIXED,0.0)
        nodeNumber = nodeNumber + numberOfCircumfrentialElementsPerQuarter
        #print "y direction  ",nodeNumber
        boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,
                                   1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,2,
                                   iron.BoundaryConditionsTypes.FIXED,0.0)
        nodeNumber = nodeNumber + numberOfCircumfrentialElementsPerQuarter
        #print "x direction  ",nodeNumber
        boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,
                                   1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,1,
                                   iron.BoundaryConditionsTypes.FIXED,0.0)
nonlinearEquations.BoundaryConditionsCreateFinish()

#
#
# ====================================================
# ====================================================
# ====================================================
# ====================================================  ENTERING THE SOLVE FOR THE FORWARD MODEL 
# ====================================================
# ====================================================
# ====================================================
#
#

maxFibreGrowthRate = 0.033
maxSheetGrowthRate = 0.033
maxNormalGrowthRate = 0.033

#rates 
growthRates = np.zeros((4,3))
growthRates[:,0] = 0.02
growthRates[:,1] = 0.02
growthRates[:,2] = 0.02
#print "growth tensor", growthRates
timesteps = int(growthRates.max()/maxFibreGrowthRate)+1
if timesteps == 0:
    timesteps = 1
#print "timesteps", timesteps
#print growthRates[:,:]/timesteps
def SetGrowthPatternParameters(growthPhase):
    if not homogeneous:
        for elementNumber in range (1,4+1):
            print growthRates[elementNumber-1,0]/timesteps, growthRates[elementNumber-1,1]/timesteps, growthRates[elementNumber-1,2]/timesteps
            for gaussPointNumber in range (1,27+1):
                growthCellMLParametersField.ParameterSetUpdateGaussPointDP(iron.FieldVariableTypes.U,
                    iron.FieldParameterSetTypes.VALUES,gaussPointNumber,elementNumber,1,growthRates[elementNumber-1,0]/timesteps)
                growthCellMLParametersField.ParameterSetUpdateGaussPointDP(iron.FieldVariableTypes.U,
                    iron.FieldParameterSetTypes.VALUES,gaussPointNumber,elementNumber,2,growthRates[elementNumber-1,1]/timesteps)
                growthCellMLParametersField.ParameterSetUpdateGaussPointDP(iron.FieldVariableTypes.U,
                    iron.FieldParameterSetTypes.VALUES,gaussPointNumber,elementNumber,3,growthRates[elementNumber-1,2]/timesteps)
                growthCellMLParametersField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
                growthCellMLParametersField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)


def SolveGrowthAndFittingProblems(time):
    timeLoop.TimesSet(time,time+timeIncrement,timeIncrement)
    # Solve the problem
    problem.Solve()
    # Set geometric field to current deformed geometry
    iron.Field.ParametersToFieldParametersComponentCopy(
        dependentField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,
        geometricField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1)
    iron.Field.ParametersToFieldParametersComponentCopy(
        dependentField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2,
        geometricField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2)
    iron.Field.ParametersToFieldParametersComponentCopy(
        dependentField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,3,
        geometricField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,3)
    # Reset growth state to 1.0
    growthCellMLStateField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                       iron.FieldParameterSetTypes.VALUES,1,1.0)
    growthCellMLStateField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                       iron.FieldParameterSetTypes.VALUES,2,1.0)
    growthCellMLStateField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                       iron.FieldParameterSetTypes.VALUES,3,1.0)    

# Export results
time = startTime
timeString = format(time)
filename = "HeartTubeGrowth_"+timeString
fields = iron.Fields()
fields.CreateRegion(region)
fields.NodesExport(filename,"FORTRAN")
fields.ElementsExport(filename,"FORTRAN")
fields.Finalise()

#SetGrowthPatternParameters(1)

if not homogeneous:
    for elementNumber in range (1,4+1):
        #print growthRates[elementNumber-1,0]/timesteps, growthRates[elementNumber-1,1]/timesteps , growthRates[elementNumber-1,2]/timesteps
        for gaussPointNumber in range (1,27+1):
            growthCellMLParametersField.ParameterSetUpdateGaussPointDP(iron.FieldVariableTypes.U,
                    iron.FieldParameterSetTypes.VALUES,gaussPointNumber,elementNumber,1,growthRates[elementNumber-1,0]/timesteps)
            growthCellMLParametersField.ParameterSetUpdateGaussPointDP(iron.FieldVariableTypes.U,
                    iron.FieldParameterSetTypes.VALUES,gaussPointNumber,elementNumber,2,growthRates[elementNumber-1,1]/timesteps)
            growthCellMLParametersField.ParameterSetUpdateGaussPointDP(iron.FieldVariableTypes.U,
                    iron.FieldParameterSetTypes.VALUES,gaussPointNumber,elementNumber,3,growthRates[elementNumber-1,2]/timesteps)
            growthCellMLParametersField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
            growthCellMLParametersField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)


# Loop over the time steps
while time<timesteps:
    print "enters the solver",str(time),"round" 
    #SolveGrowthAndFittingProblems(time)
    timeLoop.TimesSet(time,time+timeIncrement,timeIncrement)
    # Solve the problem
    print growthRates[:,:]/timesteps
    problem.Solve()
    # Set geometric field to current deformed geometry
    iron.Field.ParametersToFieldParametersComponentCopy(
        dependentField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,
        geometricField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1)
    iron.Field.ParametersToFieldParametersComponentCopy(
        dependentField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2,
        geometricField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2)
    iron.Field.ParametersToFieldParametersComponentCopy(
        dependentField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,3,
        geometricField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,3)
    # Reset growth state to 1.0
    growthCellMLStateField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                       iron.FieldParameterSetTypes.VALUES,1,1.0)
    growthCellMLStateField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                       iron.FieldParameterSetTypes.VALUES,2,1.0)
    growthCellMLStateField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                       iron.FieldParameterSetTypes.VALUES,3,1.0)    
    time = time + timeIncrement
    timeString = format(time)
    filename = "HeartTubeGrowth_"+timeString
    fields = iron.Fields()
    fields.CreateRegion(region)
    fields.NodesExport(filename,"FORTRAN")
    fields.ElementsExport(filename,"FORTRAN")
    time = time+timeIncrement

# update the geometric fields to be used to calculate the linelength the 
geometricField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
geometricField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)


# this function to initialize U variables 
def copyingFields(fromField, toField):
    for i in range(1,4):
        iron.Field.ParametersToFieldParametersComponentCopy(
                    fromField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,i,
                    toField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,i)
    toField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
    toField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)

# this function to initialize DELUDELN variables 
def copyingDervativeFields(fromField, toField):
    for i in range(1,4):
        iron.Field.ParametersToFieldParametersComponentCopy(
                    fromField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,i,
                    toField,iron.FieldVariableTypes.DELUDELN,iron.FieldParameterSetTypes.VALUES,i)
    toField.ParameterSetUpdateStart(iron.FieldVariableTypes.DELUDELN,iron.FieldParameterSetTypes.VALUES)
    toField.ParameterSetUpdateFinish(iron.FieldVariableTypes.DELUDELN,iron.FieldParameterSetTypes.VALUES)


## block finished 


def findObjective(localTime):
    nodesLocation = np.zeros((numberOfNodes, 3))   
    maxFibreGrowthRate = 0.033
    maxSheetGrowthRate = 0.033
    maxNormalGrowthRate = 0.033

    growthRates = np.zeros((4,3))
    growthRates[:,0] = 0.02
    growthRates[:,1] = 0.02
    growthRates[:,2] = 0.02
    #print "growth tensor", growthRates
    timesteps = int(growthRates.max()/maxFibreGrowthRate)+1
    if timesteps == 0:
        timesteps = 1
    #print "timesteps", timesteps
    #print growthRates[:,:]/timesteps

    print timesteps        
    try:
        #  ================================================================
        # 
        #          Restarting All The Fields   
        #
        #  ================================================================
        # restarting geometric filed 
        copyingFields(originalGeometricField,geometricField)
        # restarting dependent field  
        copyingFields(originalGeometricField,dependentField)
        # restarting DELUDELN variable  
        copyingDervativeFields(DELUDELNGeometricField,dependentField)
        # Initialise the hydrostatic pressure and its derivative
        iron.Field.ComponentValuesInitialiseDP(dependentField,iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES,4,pInit)
        iron.Field.ComponentValuesInitialiseDP(dependentField,iron.FieldVariableTypes.DELUDELN, iron.FieldParameterSetTypes.VALUES,4,0.0)   
        # restarting stress, strain and growth fields  
        for gaussPointNumber in range (1,27+1):
            for elementNumber in range (1,numberOfElements+1):
                for i in range (1,6+1):
                    dependentField.ParameterSetUpdateGaussPointDP(iron.FieldVariableTypes.U1,
                                iron.FieldParameterSetTypes.VALUES,gaussPointNumber,elementNumber,i,0.0)
                    dependentField.ParameterSetUpdateGaussPointDP(iron.FieldVariableTypes.U2,
                                iron.FieldParameterSetTypes.VALUES,gaussPointNumber,elementNumber,i,0.0)
                    if (i < 4):
                        dependentField.ParameterSetUpdateGaussPointDP(iron.FieldVariableTypes.U3,
                                    iron.FieldParameterSetTypes.VALUES,gaussPointNumber,elementNumber,i,0.0)
        dependentField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
        dependentField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)        
        # restarting constitutive fields  
        for gaussPointNumber in range (1,27+1):
            for elementNumber in range (1,numberOfElements+1):
                for i in range (1,6+1):
                    constitutiveCellMLParametersField.ParameterSetUpdateGaussPointDP(iron.FieldVariableTypes.U,
                                iron.FieldParameterSetTypes.VALUES,gaussPointNumber,elementNumber,i,0.0)
                    constitutiveCellMLIntermediateField.ParameterSetUpdateGaussPointDP(iron.FieldVariableTypes.U,
                                iron.FieldParameterSetTypes.VALUES,gaussPointNumber,elementNumber,i,0.0)
        constitutiveCellMLParametersField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
        constitutiveCellMLParametersField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)        
        constitutiveCellMLIntermediateField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
        constitutiveCellMLIntermediateField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)        
        # restarting growth fields 
        for gaussPointNumber in range (1,27+1):
            for elementNumber in range (1,numberOfElements+1):
                growthCellMLParametersField.ParameterSetUpdateGaussPointDP(iron.FieldVariableTypes.U,
                                iron.FieldParameterSetTypes.VALUES,gaussPointNumber,elementNumber,1,growthRates[elementNumber-1,0])
                growthCellMLParametersField.ParameterSetUpdateGaussPointDP(iron.FieldVariableTypes.U,
                                iron.FieldParameterSetTypes.VALUES,gaussPointNumber,elementNumber,2,growthRates[elementNumber-1,1])
                growthCellMLParametersField.ParameterSetUpdateGaussPointDP(iron.FieldVariableTypes.U,
                                iron.FieldParameterSetTypes.VALUES,gaussPointNumber,elementNumber,3,growthRates[elementNumber-1,2])
                growthCellMLParametersField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
                growthCellMLParametersField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)

        # 
        #   EXPORTING THE UNDEFORMED FIELDS AGAIN FOR THE NEW RUN 
        # 

        # exporting fields from a deformed state ...  
        time = localTime
        #timeString = format(str(localTime))
        timeString = format(time)
        print timesteps
        filename = "HeartTubeGrowth_"+timeString
        # Export results
        fields = iron.Fields()
        fields.CreateRegion(region)
        fields.NodesExport(filename,"FORTRAN")
        fields.ElementsExport(filename,"FORTRAN")
        fields.Finalise()
        for t in range(timesteps):
            timeLoop.TimesSet(time,time+timeIncrement,timeIncrement)
            try:
                print growthRates[:,:]
                problem.Solve()
            except: 
                pass
            #  
            #    R E S E T T I N G     T H E   F I E L D S      F O R      T I M E     S T E P S 
            #  
            # Set geometric field to current deformed geometry
            copyingFields(dependentField,geometricField)
            # Reset growth state to 1.0
            growthCellMLStateField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                               iron.FieldParameterSetTypes.VALUES,1,1.0)
            growthCellMLStateField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                               iron.FieldParameterSetTypes.VALUES,2,1.0)
            growthCellMLStateField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                               iron.FieldParameterSetTypes.VALUES,3,1.0)
            growthCellMLParametersField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
            growthCellMLParametersField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
            print "the round",time, "finished"
            time = time+timeIncrement
            timeString = format(time)
            filename = "HeartTubeGrowth_"+timeString
            # Export results
            fields = iron.Fields()
            fields.CreateRegion(region)
            fields.NodesExport(filename,"FORTRAN")
            fields.ElementsExport(filename,"FORTRAN")
            fields.Finalise()      

        geometricField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
        geometricField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)

    except:
        pass

for i in range (100):
    ts = 100*(i+1)
    findObjective(ts)

