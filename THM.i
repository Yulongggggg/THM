pp_ini_bc = 20e6 # Pa initial value of pore pressure
T_ini_bc = 373.15 # K  initial value of temperature
sigmaV_ini_bc = 67.3e6 # Pa vertical
sigmaH_ini_bc = 42.9e6 # Pa horizontal
nu = 0.225

a_frac   = 6e-4    # [m]  fracture aperture (开度)
phi_f    = 1.0     # [-]  intrinsic fracture porosity (常取 ~1)
k_f      = 3e-8    # [m^2] intrinsic fracture permeability (面内)

[Mesh]
  # file = fractures_copy.msh
  file = FN.msh
[]

[Functions]
  [ramp_injection]
    type = PiecewiseLinear
    x = '0   500  1000'
    y = '0     3e-1     1'
  []
  [ramp_production]
    type = PiecewiseLinear
    x = '1000  2000'
    y = '0     1'
  []
  [ramp_all]
    type = PiecewiseLinear
    x = '0   50  100'
    y = '0     0.5     1'
  []
  [inj_flux]
    type = ParsedFunction
    symbol_names = 'A R'
    symbol_values = 'A_prod ramp_injection'
    expression = '-0.2 /A * R' # kg/m^2/s
  []
  [prod_flux]
    type = ParsedFunction
    symbol_names = 'A R'
    symbol_values = 'A_prod ramp_production'
    expression = '0.2 / A * R' # kg/m^2/s
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
  [frac_phi_eff]
    type = ParsedFunction
    symbol_names  = 'a phi'
    symbol_values = '${a_frac} ${phi_f}'
    expression    = 'a * phi'          # = a * phi_f
  []
  [frac_k_eff]
    type = ParsedFunction
    symbol_names  = 'a k'
    symbol_values = '${a_frac} ${k_f}'
    expression    = 'a * k'            # = a * k_f
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
  # coupling_type = ThermoHydro
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
    # scaling = 1E-8
  []
  [temperature]
    initial_condition = ${T_ini_bc}
    # scaling = 1E-6
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

[Kernels]
  # -----------------------
  # 裂隙：流体（组件0 = 压力方程）
  # -----------------------
  [pp_mass_fracs]
    type = PorousFlowMassTimeDerivative
    variable = porepressure          # 变量=pp
    fluid_component = 0
    block = fractures
  []
  [pp_adv_fracs]
    # 全上风：等价 Action 的“full upwinding”分支
    type = PorousFlowFullySaturatedAdvectiveFlux
    variable = porepressure
    fluid_component = 0
    multiply_by_density = true       # 与文档一致，别关
    block = fractures
    gravity = '0 -9.8 0'
  []
  [pp_disp_fracs]
    # 若不要弥散/扩散，可注释掉这一条
    type = PorousFlowDispersiveFlux
    variable = porepressure
    fluid_component = 0
    disp_trans = 0
    disp_long  = 0
    block = fractures
    gravity = '0 -9.8 0'
  []
  # -----------------------
  # 裂隙：示踪（组件1 = massfrac0）
  # -----------------------
  [tr_mass_fracs]
    type = PorousFlowMassTimeDerivative
    variable = tracer_concentration
    fluid_component = 1
    gravity = '0 -9.8 0'
    block = fractures
  []
  [tr_adv_fracs]
    type = PorousFlowFullySaturatedAdvectiveFlux
    variable = tracer_concentration
    fluid_component = 1
    multiply_by_density = true
    gravity = '0 -9.8 0'
    block = fractures
  []
  [tr_disp_fracs]
    type = PorousFlowDispersiveFlux
    variable = tracer_concentration
    fluid_component = 1
    gravity = '0 -9.8 0'
    disp_trans = 0
    disp_long  = 0
    block = fractures
  []
  # -----------------------
  # （可选）裂隙：热方程
  # 只在你要让裂隙也解温度时启用；否则整段注释掉
  # -----------------------
  [T_mass_fracs]
    type = PorousFlowEnergyTimeDerivative
    gravity = '0 -9.8 0'
    variable = temperature
    block = fractures
  []
  [T_cond_fracs]
    type = PorousFlowHeatConduction
    gravity = '0 -9.8 0'
    variable = temperature
    block = fractures
  []
  [T_adv_fracs]
    type = PorousFlowFullySaturatedHeatAdvection
    gravity = '0 -9.8 0'
    variable = temperature
    multiply_by_density = true       # 文档规定：必须乘ρ
    block = fractures
  []
[]

[AuxVariables]
  [viscosity]
    order = CONSTANT
    family = MONOMIAL
  []
  [porosity]
    family = MONOMIAL
    order = CONSTANT
  []
  [enthalpy]
    family = MONOMIAL
    order  = CONSTANT
  []
  [internal_energy]
    family = MONOMIAL
    order  = CONSTANT
  []
  [density]
    family = MONOMIAL
    order  = CONSTANT
  []
  [heat_outflow]
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


[AuxKernels]
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

[FluidProperties]
  # [water]
  #   type = Water97FluidProperties
  # []
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
    # mechanical = true  may not converge, be careful before using it
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
    # mechanical = true
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
    # mechanical = true
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
  # ===== NEW: 裂隙块（fractures）的等效性质：phi_eff = a*phi_f, k_eff = a*k_f =====

  [pf_temp_fracs]
    type = PorousFlowTemperature
    block = fractures
  []
  [pf_state_fracs]
    type = PorousFlow1PhaseFullySaturated
    porepressure = porepressure
    block = fractures
  []
  [pf_massfrac_fracs]
    type = PorousFlowMassFraction
    mass_fraction_vars = tracer_concentration
    block = fractures
  []
  [pf_fluid_fracs]
    type = PorousFlowSingleComponentFluid
    fp = water
    phase = 0
    block = fractures
  []
  [nearestqp_fracs]
    type = PorousFlowNearestQp
    block = fractures
  []


  [phi_fractures]
    type = PorousFlowPorosityConst
    block = fractures
    porosity = 6e-4              # <- 与 frac_phi_eff 对应（默认 a=6e-4, phi_f=1）
  []
  [k_fractures]
    type = PorousFlowPermeabilityConst
    block = fractures
    permeability = '1.8e-11 0 0   0 1.8e-11 0   0 0 1.8e-11'   # <- 与 frac_k_eff 对应（a*k_f）
  []
  [kth_fracs]
    type = PorousFlowThermalConductivityIdeal
    PorousFlowDictator = dictator
    dry_thermal_conductivity = '0 0 0   0 0 0   0 0 0'
    block = fractures
  []
  [diff_fractures]
    type = PorousFlowDiffusivityConst
    block = fractures
    diffusion_coeff = '1e-9 1e-9'
    tortuosity = 1.0
  []
  [diff_matrix]
    type = PorousFlowDiffusivityConst
    block = 'cap1 cap2 block'
    diffusion_coeff = '1e-9 1e-9'
    tortuosity = 0.2
  []
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
  [inj_flux]
    type = FunctionValuePostprocessor
    function = inj_flux
    execute_on = 'initial timestep_begin'
  []
  [power_prod]
    type      = SideIntegralVariablePostprocessor
    boundary  = production
    variable  = heat_outflow          # W/m²
    execute_on = 'TIMESTEP_END'
  []
[]


[BCs]
  [m_in]
    type = PorousFlowSink
    boundary = injection
    variable = porepressure
    flux_function = inj_flux
    fluid_phase = 0
  []
  [in_T]
    type = PorousFlowEnthalpySink
    variable = temperature
    boundary = injection
    T_in = 300
    fp = water
    flux_function = inj_flux
    porepressure_var = porepressure
  []
  [m_out]
    type = PorousFlowSink
    boundary = production
    variable = porepressure
    flux_function = prod_flux
    fluid_phase = 0
    use_mobility = False
    use_relperm = False
  []
  [m_out_heat]
    type = PorousFlowSink
    boundary = production
    variable = temperature
    flux_function = prod_flux
    fluid_phase = 0
    use_enthalpy = true
    save_in = heat_outflow
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
  # [outflow]
  #   type = PorousFlowOutflowBC
  #   PorousFlowDictator = dictator
  #   gravity = '0 -9.8 0'
  #   boundary = 'l r b t'
  #   flux_type = fluid
  #   variable = porepressure
  #   include_relperm = False
  # []
  [T_BC]
    type = DirichletBC
    variable = temperature
    boundary = 'l r t b'
    value = ${T_ini_bc}
  []
  [pp_BC]
    type = DirichletBC
    variable = porepressure
    boundary = 'l r t b'
    value = ${pp_ini_bc}
  []
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
  dt = 20
  [TimeStepper]
    type = IterationAdaptiveDT
    growth_factor = 1.5
    optimal_iterations = 12
    cutback_factor = 0.7
    dt = 20
  []
  end_time = 1e9 # near 32 years
[]

[VectorPostprocessors]
  [point_sample1]
    type = PositionsFunctorValueSampler
    discontinuous = False
    functors = 'disp_x disp_y porepressure temperature'
    # functors = 'porepressure temperature'
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
