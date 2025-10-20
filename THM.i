pp_ini_bc = 30e6 # Pa initial value of pore pressure
T_ini_bc = 373.15 # K  initial value of temperature
sigmaV_ini_bc = 67.3e6 # Pa vertical
sigmaH_ini_bc = 42.9e6 # Pa horizontal
nu = 0.225

[Mesh]
  file = fractures.msh
[]

[Postprocessors]
  [A_inj]
    type = AreaPostprocessor
    boundary = injection
    execute_on = initial
  []
  [A_prod]
    type = AreaPostprocessor
    boundary = production
    execute_on = initial
  []
[]

[Functions]
  [ramp_all]
    type = PiecewiseLinear
    x = '0   1000  5000'
    y = '0     1e-3     1'
  []
  [ramp]
    type = PiecewiseLinear
    x = '0   5000  10000'
    y = '0     0.5     1'
  []
  [inj_flux]
    type = ParsedFunction
    symbol_names = 'A R'
    symbol_values = 'A_inj ramp_all'
    expression = '-3 / A * R' # kg/m^2/s
  []
  [prod_flux]
    type = ParsedFunction
    symbol_names = 'A R'
    symbol_values = 'A_prod ramp_all'
    expression = '3 / A * R' # kg/m^2/s
  []
  [sigmaV_flux]
    type = ParsedFunction
    symbol_names = 'R'
    symbol_values = 'ramp_all'
    expression = 'R * ${sigmaV_ini_bc}' # 最终值 = 67.3e6 Pa
  []
  [sigmaH_flux]
    type = ParsedFunction
    symbol_names = 'R'
    symbol_values = 'ramp_all'
    expression = 'R * ${sigmaH_ini_bc}' # 最终值 = 42.9e6 Pa
  []
  [sigma_zz]
    type = ParsedFunction
    symbol_names = 'sxx syy nu'
    symbol_values = '${sigmaH_ini_bc} ${sigmaV_ini_bc} ${nu}' # 0.225 是你的泊松比
    expression = '-1*nu * (sxx + syy)'
  []
[]

[GlobalParams]
  PorousFlowDictator = dictator
  displacements = 'disp_x disp_y'
[]

[PorousFlowFullySaturated]
  fluid_properties_type = PorousFlowSingleComponentFluid
  add_darcy_aux = true
  coupling_type = ThermoHydroMechanical
  gravity = '0 -9.8 0'
  porepressure = porepressure
  temperature = temperature
  fp = water
  mass_fraction_vars = tracer_concentration
  use_displaced_mesh = False
  biot_coefficient = 0.9
  eigenstrain_names = 'thermal_contribution'
  block = 'cap1 cap2 block'
  # relative_permeability_exponent = 3
  # relative_permeability_type = Corey
  # residual_saturation = 0.1
  # van_genuchten_alpha = 1E-6
  # van_genuchten_m = 0.6
[]

[Variables]
  [porepressure]
    initial_condition = ${pp_ini_bc}
  []
  [temperature]
    initial_condition = ${T_ini_bc}
  []
  [disp_x]
    initial_condition = 0
  []
  [disp_y]
    initial_condition = 0
  []
  [tracer_concentration]
  []
[]

[AuxVariables]
  [enthalpy]
  []
  [internal_energy]
  []
  [density]
  []
[]

[ICs]
  [tracer_concentration]
    type = FunctionIC
    function = '0.5*if((x)*(x) + (y+50)*(y+50) < 0.1*0.1 , 1, 0 )'
    variable = tracer_concentration
  []
  [enthalpy]
    type = PorousFlowFluidPropertyIC
    variable = enthalpy
    property = enthalpy
    porepressure = porepressure
    temperature = temperature
    fp = water
  []
  [internal_energy]
    type = PorousFlowFluidPropertyIC
    variable = internal_energy
    property = internal_energy
    porepressure = porepressure
    temperature = temperature
    fp = water
  []
  [density]
    type = PorousFlowFluidPropertyIC
    variable = density
    property = density
    porepressure = porepressure
    temperature = temperature
    fp = water
  []
[]

[AuxVariables]
  [density]
    order = CONSTANT
    family = MONOMIAL
  []
  [viscosity]
    order = CONSTANT
    family = MONOMIAL
  []
  [porosity]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[AuxKernels]
  [stress_xx]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xx
    index_i = 0
    index_j = 0
  []
  [stress_yy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_yy
    index_i = 1
    index_j = 1
  []
  [density]
    type = MaterialRealAux
    variable = density
    property = PorousFlow_fluid_phase_density_qp0
    execute_on = TIMESTEP_END
  []
  [viscosity]
    type = MaterialRealAux
    variable = viscosity
    property = PorousFlow_viscosity_qp0
    execute_on = TIMESTEP_END
  []
  [porosity]
    type = MaterialRealAux
    property = PorousFlow_porosity_qp
    variable = porosity
    execute_on = TIMESTEP_END
  []
[]

# [FluidProperties]
#   [water]
#     type = Water97FluidProperties
#   []
# []

[FluidProperties]
  [water]
    type = SimpleFluidProperties
    bulk_modulus = 2.27e9
    density0 = 970.0
    viscosity = 0.3394e-3
    cv = 4149.0
    cp = 4149.0
    porepressure_coefficient = 1
    thermal_expansion = 0.000214
  []
[]

[Materials]
  [phi_cap1]
    type = PorousFlowPorosity
    block = cap1
    porosity_zero = 0.01
    mechanical = true
    fluid = true
    thermal = true
    thermal_expansion_coeff = 1.0e-5
    solid_bulk = 1
  []
  [./strain]
    type = ComputeSmallStrain
    block = 'cap1 cap2 block'
    eigenstrain_names = 'thermal_contribution'
  [../]
  [k_cap1]
    type = PorousFlowPermeabilityConst
    block = cap1
    permeability = '1e-15 0 0   0 1e-15 0   0 0 1e-15'
  []
  [phi_blk]
    type = PorousFlowPorosity
    block = block
    porosity_zero = 0.15
    mechanical = true
    fluid = true
    thermal = true
    thermal_expansion_coeff = 1.0e-5
    solid_bulk = 1
  []
  [k_blk]
    type = PorousFlowPermeabilityKozenyCarman
    block = block
    PorousFlowDictator = dictator
    poroperm_function = kozeny_carman_phi0
    k0 = 1e-13
    phi0 = 0.15
    n = 3
    m = 2
  []
  [phi_cap2]
    type = PorousFlowPorosity
    block = cap2
    porosity_zero = 0.01
    mechanical = true
    fluid = true
    thermal = true
    thermal_expansion_coeff = 1.0e-5
    solid_bulk = 1.0 # Required but irrelevant when biot_coefficient is unity
  []
  [k_cap2]
    type = PorousFlowPermeabilityConst
    block = cap2
    permeability = '1e-15 0 0   0 1e-15 0   0 0 1e-15'
  []
  [elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    shear_modulus = 8.38e9
    poissons_ratio = 0.225
    block = 'cap1 cap2 block'
  []
  [thermal_contribution]
    type = ComputeThermalExpansionEigenstrain
    temperature = temperature
    stress_free_temperature = ${T_ini_bc}
    thermal_expansion_coeff = 1.0e-5
    eigenstrain_name = thermal_contribution
    block = 'cap1 cap2 block'
  []
  [stress]
    type = ComputeLinearElasticStress
  []
  [solid_internal_energy]
    type = PorousFlowMatrixInternalEnergy
    PorousFlowDictator = dictator
    specific_heat_capacity = 920.0
    density = 2600.0
    block = 'cap1 block cap2'
  []
  [kth_eff]
    type = PorousFlowThermalConductivityIdeal
    PorousFlowDictator = dictator
    dry_thermal_conductivity = '2.5 0 0   0 2.5 0   0 0 2.5'
    block = 'cap1 block cap2'
  []
  [density]
    type = GenericConstantMaterial
    prop_names = 'density'
    prop_values = '2600' # kg/m^3
    block = 'cap1 block cap2'
  []
#   [phi_f1]
#     type = PorousFlowPorosity
#     block = f1
#     porosity_zero = 0.25
#     mechanical = true
#     fluid = true
#     thermal = true
#     thermal_expansion_coeff = 1.2e-5
#     solid_bulk = 0.8
#   []
#   [k_f1]
#     type = PorousFlowPermeabilityConst
#     block = f1
#     permeability = '1e-11 0 0   0 1e-11 0   0 0 1e-11'
#   []
#   [elasticity_tensor_f1]
#     type = ComputeIsotropicElasticityTensor
#     shear_modulus = 2e9
#     poissons_ratio = 0.30
#     block = f1
#   []

#   [phi_f2]
#     type = PorousFlowPorosity
#     block = f2
#     porosity_zero = 0.25
#     mechanical = true
#     fluid = true
#     thermal = true
#     thermal_expansion_coeff = 1.2e-5
#     solid_bulk = 0.8
#   []
#   [k_f2]
#     type = PorousFlowPermeabilityConst
#     block = f2
#     permeability = '1e-11 0 0   0 1e-11 0   0 0 1e-11'   # 高渗透断层
#   []
#   [elasticity_tensor_f2]
#     type = ComputeIsotropicElasticityTensor
#     shear_modulus = 2e9
#     poissons_ratio = 0.30
#     block = f2
#   []
[]

[BCs]
  [m_in]
    type = PorousFlowSink
    boundary = injection
    variable = porepressure
    flux_function = inj_flux
    fluid_phase = 0
  []
  [left_T]
    type = PorousFlowEnthalpySink
    variable = temperature
    boundary = injection
    T_in = 300
    fp = wa
    flux_function = inj_flux
  []
  [m_out]
    type = PorousFlowSink
    boundary = production
    variable = porepressure
    flux_function = prod_flux
    fluid_phase = 0
    use_mobility = true
    use_relperm = true
    use_enthalpy = true
  []
  [m_out_heat]
    type = PorousFlowSink
    boundary = production
    variable = temperature
    flux_function = prod_flux
    fluid_phase = 0
    use_mobility = true
    use_relperm = true
    use_enthalpy = true
  []
  [roller_bottom_x]
    type = DirichletBC
    variable = disp_x
    boundary = 'b'
    value = 0
  []
  [roller_bottom]
    type = DirichletBC
    variable = disp_y
    boundary = 'b'
    value = 0
  []
  [injection_dispx]
    type = DirichletBC
    variable = disp_x
    boundary = 'injection'
    value = 0
  []
  [injection_dispy]
    type = DirichletBC
    variable = disp_y
    boundary = 'injection'
    value = 0
  []
  [production_dispy]
    type = DirichletBC
    variable = disp_y
    boundary = 'production'
    value = 0
  []
  [production_dispx]
    type = DirichletBC
    variable = disp_x
    boundary = 'production'
    value = 0
  []
  [traction_top_bottom]
    type = Pressure
    variable = disp_y
    boundary = 't'
    factor = ${sigmaV_ini_bc}
  []
  [traction_left]
    type = Pressure
    variable = disp_x
    boundary = 'l'
    factor = ${sigmaH_ini_bc}
  []
  [traction_right]
    type = Pressure
    variable = disp_x
    boundary = 'r'
    factor = ${sigmaH_ini_bc}
  []
  # [Porepressure_BC]
  #   type = DirichletBC
  #   variable = porepressure
  #   boundary = 't b l r'
  #   value = ${pp_ini_bc}
  # []
  # [T_BC]
  #   type = DirichletBC
  #   variable = temperature
  #   boundary = 't b l r'
  #   value = ${T_ini_bc}
  # []
  [injected_tracer]
    type = DirichletBC
    variable = tracer_concentration
    value = 0.5
    boundary = injection
  []
[]

[Preconditioning]
  [SMP]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
    petsc_options_value = 'lu mumps'
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  line_search = 'none'
  automatic_scaling = true
  nl_abs_tol = 1e-5
  nl_rel_tol = 1e-5
  l_tol = 1e-6
  l_max_its = 1000
  dt = 1
  [TimeStepper]
    type = IterationAdaptiveDT
    growth_factor = 1.3
    optimal_iterations = 12
    cutback_factor = 0.5
    dt = 1
  []
  end_time = 1e8
[]

[VectorPostprocessors]
  [point_sample1]
    type = PositionsFunctorValueSampler
    discontinuous = False
    functors = 'disp_x disp_y porepressure temperature'
    positions = 'all_elements'
    sort_by = 'x'
    execute_on = TIMESTEP_END
  []
[]

[Positions]
  [all_elements]
    type = ElementCentroidPositions
  []
[]

[Outputs]
  exodus = true
  # [CSV]
  #   type = CSV
  #   execute_on = FINAL
  # []
[]

# [Debug]
#   show_actions = true
# []
