classdef NetworkConfiguration_Class
    
    properties
        % Network Detailed Properties
        Network_Name
        PlantSpecies
        ControllerSpecies
        PlantCoefficients_Reactants
        PlantCoefficients_Products
        ControllerCoefficients_Reactants
        ControllerCoefficients_Products
        PlantPropensities
        ControllerPropensities
        PlantParameters_String
        ControllerParameters_String
        % Image
        ImageFileName
        % Plant
        PlantParameters     
        N_PlantKnobs
        PlantKnobs_Labels
        PlantKnobs_InitialValues
        PlantKnobs_InitialRanges
        % Controller
        ControlParameters
        N_ControlKnobs
        ControlKnobs_Labels
        ControlKnobs_InitialValues
        ControlKnobs_InitialRanges
        % Stoichiometry Matrix
        StoichiometryMatrix
        % Propensity Function
        PropensityFunction
        % Initial Condition
        InitialCondition
        % Names of State Variables (for each Species)
        Names_Deterministic
        Names_Stochastic
        % Output Index
        OutputIndex
        % Tuning Curves Panel
        TuningCurvesPanel_Flag = 'off'
        
    end
    
    methods
    end
end

