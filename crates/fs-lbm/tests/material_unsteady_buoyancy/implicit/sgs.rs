//! Native implicit material/SGS composition: physical heat and work oracles.
use super::*;
use fs_lbm::conjugate::{
    AffineSpecificEnthalpy, ConvectionScheme, FvSgsConfig, TransientConfig, march_energy,
};
use fs_lbm::material::conjugate::buoyant::enthalpy::{
    ConvectionCaloricLaw, MaterialEnthalpyConfig, MaterialEnthalpyTransientConfig,
};
use fs_lbm::sgs3::{SgsModel, SgsUse};

fn sgs_controls(cs: f64) -> FvSgsConfig {
    FvSgsConfig {
        model: SgsModel::Smagorinsky { cs },
        filter_width_factor: 1.0,
        maximum_eddy_to_molecular: 100.0,
        usage: SgsUse::Exploration,
    }
}

fn shear_flow() -> SimpleConfig {
    let mut faces = [FvBoundary::wall(); 6];
    faces[3] = FvBoundary::Wall {
        velocity: [0.05, 0.0, 0.0],
    };
    faces[4] = FvBoundary::Symmetry;
    faces[5] = FvBoundary::Symmetry;
    SimpleConfig::new(faces)
}

#[test]
fn implicit_live_sgs_heat_uses_exactly_the_final_projected_flow_epoch() {
    let domain = VoxelDomain::new(4, 4, 1, 0.25).unwrap();
    let mut flow = shear_flow();
    flow.turbulent_prandtl = 0.63;
    let mut thermal = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    thermal.power_w = vec![0.01; domain.cell_count()];
    let saved = thermal.clone();
    let initial: Vec<_> = (0..domain.cell_count())
        .map(|c| 300.0 + 4.0 * domain.coords(c)[0] as f64 / 3.0)
        .collect();
    let reference = FluidTransportState::try_new(1.0, 2.0, 0.002, 0.01).unwrap();
    let bound = MaterialChtFluid::try_new(&reference).unwrap();
    let mut unsteady = UnsteadyConfig::new(0.2, 1);
    unsteady.scheme = TimeScheme::BackwardEuler;
    unsteady.inner_iterations = 1000;
    unsteady.inner_tolerance = 1e-11;
    let mut outer = endpoint_controls();
    outer.temperature_tolerance_k = 1e-8;
    outer.transport_relative_tolerance = 1e-8;
    let energy = EnergyConfig::default();
    let gate = CancelGate::new();
    let resolve = |_: usize, t: f64, _: &FluidTransportState| {
        FluidTransportState::try_new(1.0, 2.0, 0.002, 0.01 / (1.0 + 0.15 * (t - 300.0)))
    };
    let run = |cs| {
        bound
            .march_natural_convection_with_sgs_implicit(
                &domain,
                &[],
                &flow,
                &thermal,
                &energy,
                &initial,
                limits([0.0; 3]),
                &unsteady,
                sgs_controls(cs),
                outer,
                |_| 1.0,
                |_| 1.0,
                None,
                &gate,
                resolve,
            )
            .unwrap()
    };
    let result = run(0.5);
    let zero = run(0.0);
    let laminar = bound
        .march_natural_convection_implicit(
            &domain,
            &[],
            &flow,
            &thermal,
            &energy,
            &initial,
            limits([0.0; 3]),
            &unsteady,
            outer,
            |_| 1.0,
            |_| 1.0,
            None,
            &gate,
            resolve,
        )
        .unwrap();
    assert_eq!(zero.solution.march, laminar.solution.march);
    let march = &result.solution.march;
    assert!(result.records[0].passes > 1);
    assert_eq!(march.flow.records.len(), 1);
    assert_eq!(march.flow.averaged_steps, 1);
    assert_eq!(march.energy_records.len(), 1);
    assert_eq!(result.energy_solves, result.momentum_solves);
    assert!(
        march
            .flow
            .flow
            .eddy_viscosity_m2_s
            .iter()
            .any(|&nu| nu > 1e-7)
    );
    let velocity_difference = march
        .flow
        .flow
        .velocity_m_s
        .iter()
        .zip(&zero.solution.march.flow.flow.velocity_m_s)
        .flat_map(|(a, b)| a.iter().zip(b))
        .map(|(a, b)| (a - b).abs())
        .fold(0.0_f64, f64::max);
    assert!(velocity_difference > 1e-6);
    for (applied, endpoint) in result
        .solution
        .applied
        .iter()
        .zip(&result.solution.endpoint)
    {
        assert!(
            (applied.temperature_k - endpoint.temperature_k).abs() <= outer.temperature_tolerance_k
        );
        close(
            endpoint
                .material
                .resolved_transport()
                .thermal_conductivity_w_m_k(),
            0.002,
            1e-14,
        );
    }

    // Constant molecular k makes an independent native one-step heat replay
    // exact on the returned projected flow and its CURRENT eddy conductivity.
    // Zero/twice controls detect omitted or compounded turbulent heat transport.
    let replay = |scale: f64| {
        let mut setup = thermal.clone();
        setup.eddy_conductivity_w_m_k = march
            .flow
            .flow
            .heat_eddy_viscosity_m2_s
            .iter()
            .map(|nu_t| scale * 2.0 * nu_t / 0.63)
            .collect();
        march_energy(
            &domain,
            bound.properties(),
            &[],
            &march.flow.flow.field,
            &setup,
            &initial,
            |_| 1.0,
            &TransientConfig {
                time_step_s: unsteady.time_step_s,
                steps: 1,
                energy,
            },
            None,
            &gate,
        )
        .unwrap()
    };
    let error = |other: &[f64]| {
        march
            .temperature
            .iter()
            .zip(other)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0_f64, f64::max)
    };
    assert!(error(&replay(1.0).temperature) < 1e-8);
    assert!(error(&replay(0.0).temperature) > 1e-5);
    assert!(error(&replay(2.0).temperature) > 1e-5);
    let stored: f64 = march
        .temperature
        .iter()
        .zip(&initial)
        .enumerate()
        .map(|(cell, (t, old))| 2.0 * domain.volume(cell) * (t - old))
        .sum();
    let authored = unsteady.time_step_s * thermal.power_w.iter().sum::<f64>();
    close(stored, authored, 1e-8);
    close(march.energy_records[0].source_j, authored, 1e-12);
    assert!(march.energy_records[0].closure_j.abs() < 1e-8);
    assert_eq!(thermal, saved);
}

// The same exact synthetic source law used by the explicit SGS caloric
// regression. This is a constitutive numerical oracle, not measured fluid data.
struct Caloric;
impl ConvectionCaloricLaw for Caloric {
    fn tangent_at_temperature(
        &self,
        reference: f64,
        t: f64,
    ) -> Result<AffineSpecificEnthalpy, MaterialError> {
        Ok(AffineSpecificEnthalpy {
            temperature_k: t,
            enthalpy_j_kg: (t - reference) * (2.0 + 0.025 * (t + reference - 600.0)),
            heat_capacity_j_kg_k: cp(t),
        })
    }
}

fn cp(t: f64) -> f64 {
    2.0 + 0.05 * (t - 300.0)
}
fn h(t: f64) -> f64 {
    let theta = t - 300.0;
    2.0 * theta + 0.025 * theta * theta
}
fn molecular_k(t: f64) -> f64 {
    0.002 * (1.0 + 0.04 * (t - 300.0))
}
fn transport(t: f64, mu_slope: f64) -> Result<FluidTransportState, MaterialError> {
    FluidTransportState::try_new(
        1.0,
        cp(t),
        molecular_k(t),
        0.01 / (1.0 + mu_slope * (t - 300.0)),
    )
}
fn caloric_controls() -> MaterialEnthalpyTransientConfig {
    MaterialEnthalpyTransientConfig {
        caloric: MaterialEnthalpyConfig {
            reference_temperature_k: 300.0,
            maximum_flux_defect_w: 1e-9,
        },
        maximum_passes_per_step: 100,
        maximum_energy_solves: 1000,
        temperature_tolerance_k: 1e-9,
        maximum_storage_defect_w: 1e-8,
        maximum_energy_closure_j: 1e-8,
        temperature_relaxation: 1.0,
    }
}
fn thermal_pair() -> (
    VoxelDomain,
    ThermalSetup,
    UnsteadyConfig,
    BoussinesqConfig<3>,
    EnergyConfig,
) {
    let domain = VoxelDomain::new(2, 1, 1, 0.25).unwrap();
    let mut thermal = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    thermal.fixed_temperature = vec![(0, 320.0)];
    thermal.power_w = vec![0.0, 0.02];
    let mut unsteady = UnsteadyConfig::new(0.5, 1);
    unsteady.scheme = TimeScheme::BackwardEuler;
    unsteady.inner_iterations = 1000;
    unsteady.inner_tolerance = 1e-11;
    let mut bounds = limits([0.0; 3]);
    bounds.max_relative_specific_heat_departure = 0.9;
    let energy = EnergyConfig {
        scheme: ConvectionScheme::Upwind,
        tolerance: 1e-13,
        ..EnergyConfig::default()
    };
    (domain, thermal, unsteady, bounds, energy)
}

#[test]
fn implicit_molecular_k_and_source_cp_sgs_heat_close_the_independent_reciprocal_law() {
    let (domain, thermal, unsteady, bounds, energy) = thermal_pair();
    let reference = transport(300.0, 0.02).unwrap();
    let bound = MaterialChtFluid::try_new(&reference).unwrap();
    let caloric = Caloric;
    let result = bound
        .march_natural_convection_with_enthalpy_and_sgs_implicit(
            &domain,
            &[],
            &shear_flow(),
            &thermal,
            &energy,
            &[300.0; 2],
            bounds,
            &caloric,
            caloric_controls(),
            &unsteady,
            sgs_controls(0.5),
            endpoint_controls(),
            |_| 1.0,
            |_| 1.0,
            None,
            &CancelGate::new(),
            |_, t, _| transport(t, 0.02),
        )
        .unwrap();
    let solution = &result.solution;
    let march = &solution.convection.march;
    // These closed MAC cells have zero resolved volume flow. The original
    // wall derivative still gives analytic nonzero Smagorinsky shear.
    let nu_t = (0.5_f64 * 0.25).powi(2) * (0.05 / 0.25);
    for &actual in &march.flow.flow.eddy_viscosity_m2_s {
        close(actual, nu_t, 1e-12);
    }
    let heat = |t: f64| {
        let hot = molecular_k(320.0) + cp(320.0) * nu_t / 0.9;
        let cold = molecular_k(t) + cp(t) * nu_t / 0.9;
        0.0625 * (320.0 - t) / (0.125 / hot + 0.125 / cold)
    };
    let mut bracket = [300.0, 320.0];
    for _ in 0..80 {
        let t = 0.5 * (bracket[0] + bracket[1]);
        if domain.volume(1) * h(t) / unsteady.time_step_s > heat(t) + 0.02 {
            bracket[1] = t;
        } else {
            bracket[0] = t;
        }
    }
    let t = march.temperature[1];
    assert!((t - 0.5 * (bracket[0] + bracket[1])).abs() < 1e-7);
    close(
        domain.volume(1) * h(t) / unsteady.time_step_s - 0.02,
        heat(t),
        1e-8,
    );
    let stored = domain.volume(0) * h(320.0) + domain.volume(1) * h(t);
    close(march.energy_records[0].stored_energy_change_j, stored, 1e-9);
    close(
        stored + march.energy_records[0].boundary_outflow_j,
        march.energy_records[0].source_j,
        1e-8,
    );
    assert!(result.records[0].passes > 1);
    assert_eq!(solution.records.len(), 1);
    assert_eq!(march.flow.records.len(), 1);
    assert_eq!(result.energy_solves, solution.records[0].energy_solves);
    assert!(
        solution.records[0].maximum_sgs_heat_flux_defect_w
            <= caloric_controls().caloric.maximum_flux_defect_w
    );
    for sample in &solution.convection.endpoint {
        close(
            sample
                .material
                .resolved_transport()
                .thermal_conductivity_w_m_k(),
            molecular_k(sample.temperature_k),
            1e-14,
        );
    }
}

#[test]
fn implicit_sgs_retains_local_heat_cap_budget_and_exploration_refusals() {
    let (domain, thermal, unsteady, bounds, energy) = thermal_pair();
    let reference = transport(300.0, 0.02).unwrap();
    let bound = MaterialChtFluid::try_new(&reference).unwrap();
    let caloric = Caloric;
    let saved = thermal.clone();
    let run = |heat_controls, outer, sgs, mu_slope| {
        bound.march_natural_convection_with_enthalpy_and_sgs_implicit(
            &domain,
            &[],
            &shear_flow(),
            &thermal,
            &energy,
            &[300.0; 2],
            bounds,
            &caloric,
            heat_controls,
            &unsteady,
            sgs,
            outer,
            |_| 1.0,
            |_| 1.0,
            None,
            &CancelGate::new(),
            |_, t, _| transport(t, mu_slope),
        )
    };
    let mut one_pass = caloric_controls();
    one_pass.maximum_passes_per_step = 1;
    one_pass.temperature_tolerance_k = 100.0;
    one_pass.maximum_storage_defect_w = 100.0;
    one_pass.maximum_energy_closure_j = 100.0;
    assert!(
        matches!(run(one_pass, endpoint_controls(), sgs_controls(0.5), 0.02),
        Err(MaterialBuoyancyError::SgsHeatNotConverged { maximum_flux_defect_w, .. })
            if maximum_flux_defect_w > 1e-4)
    );
    let mut capped = sgs_controls(0.5);
    capped.maximum_eddy_to_molecular = 1.0;
    // Initial nu_t/nu=.3125 passes. The genuine hot source endpoint has
    // mu=.001 and ratio 3.125; physical admission must refuse that endpoint.
    assert!(matches!(
        run(caloric_controls(), endpoint_controls(), capped, 0.45),
        Err(MaterialBuoyancyError::Solver(ChtError::InvalidInput {
            field: "simple.sgs",
            ..
        }))
    ));
    // Deliberately admit the entire temperature jump and its 0.9 relative
    // viscosity defect in ONE outer pass. A second momentum solve cannot
    // rescue a missing actual-endpoint applicability check.
    let mut loose = endpoint_controls();
    loose.maximum_passes_per_step = 1;
    loose.temperature_tolerance_k = 100.0;
    loose.transport_relative_tolerance = 0.95;
    assert!(matches!(
        run(caloric_controls(), loose, capped, 0.45),
        Err(MaterialBuoyancyError::Solver(ChtError::InvalidInput {
            field: "simple.sgs",
            ..
        }))
    ));
    let mut acceptance = sgs_controls(0.5);
    acceptance.usage = SgsUse::AcceptanceRender;
    assert!(matches!(
        run(caloric_controls(), endpoint_controls(), acceptance, 0.02),
        Err(MaterialBuoyancyError::Solver(ChtError::InvalidInput {
            field: "simple.sgs",
            ..
        }))
    ));
    let mut outer = endpoint_controls();
    outer.maximum_energy_solves = 1;
    assert!(matches!(
        run(caloric_controls(), outer, sgs_controls(0.5), 0.02),
        Err(MaterialBuoyancyError::Solver(ChtError::InvalidInput {
            field: "unsteady.material_energy_solves",
            ..
        }))
    ));
    assert_eq!(thermal, saved);
    assert_eq!(bound.transport(), &reference);
}
