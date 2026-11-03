//! Synthetic manufactured transport and native caloric closure regressions.
//! These exercise numerical behavior, not experimental material qualification.

use super::*;
use fs_lbm::conjugate::{
    AdvectiveEnthalpyField, AffineSpecificEnthalpy, ConvectionScheme, FlowFace,
    FlowField, FluidProperties, simple_flow, solve_energy,
};
use fs_lbm::material::conjugate::buoyant::enthalpy::{
    ConvectionCaloricLaw, MaterialEnthalpyConfig, MaterialEnthalpyTransientConfig,
};

fn energy(scheme: ConvectionScheme) -> EnergyConfig {
    EnergyConfig { scheme, tolerance: 1e-13, ..EnergyConfig::default() }
}

fn properties() -> FluidProperties {
    FluidProperties {
        density_kg_m3: 1.0,
        // The supplied caloric field, not this unused reference Cp, owns heat.
        specific_heat_j_kg_k: 999.0,
        conductivity_w_m_k: 0.7,
        kinematic_viscosity_m2_s: 1.0,
    }
}

#[test]
fn g1_g3_graded_linear_enthalpy_transport_preserves_direction_rest_and_datum() {
    let domain = VoxelDomain::graded([
        vec![0.08, 0.13, 0.21, 0.11, 0.19, 0.09, 0.15, 0.17],
        vec![1.0], vec![1.0],
    ]).unwrap();
    let count = domain.cell_count();
    let length = domain.face_coord(0, count);
    let slope = 7.0;
    let cp = 2.0;
    let exact: Vec<_> = (0..count)
        .map(|cell| 300.0 + slope * domain.center(cell, 0, 0)[0]).collect();
    let capacity: Vec<_> = (0..count).map(|cell| cp * domain.volume(cell)).collect();
    let conductivity = vec![[0.7; 3]; count];
    for velocity in [-0.2, 0.0, 0.2] {
        let mut boundaries = [FlowFace::Wall; 6];
        boundaries[0] = FlowFace::Fixed;
        boundaries[1] = FlowFace::Fixed;
        let flow = FlowField::from_face_velocity(&domain, boundaries, |axis, _| {
            if axis == 0 { velocity } else { 0.0 }
        }).unwrap();
        let mut faces = [ThermalFace::Adiabatic; 6];
        faces[0] = ThermalFace::Outflow { backflow_temperature: 300.0 };
        faces[1] = ThermalFace::Outflow { backflow_temperature: 300.0 + slope * length };
        let mut thermal = ThermalSetup::new(faces);
        // For T=300+s*x, the exact integrated steady equation is
        // Q_cell=rho*u*A*Cp*s*width; constant conductive flux cancels.
        thermal.power_w = (0..count)
            .map(|cell| velocity * cp * slope * domain.widths(cell)[0]).collect();
        // Keep two end layers exact so the first-order open-face convention
        // does not enter the manufactured interior truncation comparison.
        thermal.fixed_temperature = [0, 1, count - 2, count - 1]
            .into_iter().map(|cell| (cell, exact[cell])).collect();
        for datum in [0.0, 131_072.0] {
            let tangents: Vec<_> = exact.iter().map(|&t| Some(AffineSpecificEnthalpy {
                temperature_k: t,
                enthalpy_j_kg: datum + cp * (t - 300.0),
                heat_capacity_j_kg_k: cp,
            })).collect();
            let field = AdvectiveEnthalpyField {
                cells: &tangents,
                boundary_j_kg: [
                    Some(datum), Some(datum + cp * slope * length),
                    None, None, None, None,
                ],
            };
            let run = |scheme| energy(scheme).step_with_enthalpy(
                &domain, &properties(), &[], &flow, &thermal, &conductivity,
                field, &capacity, &exact, 0.1, &CancelGate::new(),
            ).unwrap();
            let (limited, receipt) = run(ConvectionScheme::EnthalpyMinmod);
            let (upwind, _) = run(ConvectionScheme::Upwind);
            let error = |temperature: &[f64]| temperature.iter().zip(&exact)
                .map(|(a, b)| (a - b).abs()).fold(0.0_f64, f64::max);
            assert!(error(&limited.temperature) < 3e-8);
            assert!(receipt.closure_j.abs() < 1e-8);
            if velocity != 0.0 {
                assert!(error(&upwind.temperature) > 1e-4,
                    "graded donor-cell transport must distinguish this manufactured profile");
            } else {
                assert_eq!(limited.temperature, upwind.temperature);
            }
        }
    }
}

struct QuadraticCaloric;

fn enthalpy(t: f64) -> f64 {
    let theta = t - 300.0;
    2.0 * theta + 0.025 * theta * theta
}

impl ConvectionCaloricLaw for QuadraticCaloric {
    fn tangent_at_temperature(&self, reference: f64, t: f64)
        -> Result<AffineSpecificEnthalpy, MaterialError>
    {
        Ok(AffineSpecificEnthalpy {
            temperature_k: t,
            enthalpy_j_kg: (t - reference)
                * (2.0 + 0.025 * (t + reference - 600.0)),
            heat_capacity_j_kg_k: 2.0 + 0.05 * (t - 300.0),
        })
    }
}

fn caloric_controls(reference_temperature_k: f64) -> MaterialEnthalpyTransientConfig {
    MaterialEnthalpyTransientConfig {
        caloric: MaterialEnthalpyConfig {
            reference_temperature_k, maximum_flux_defect_w: 1e-9,
        },
        maximum_passes_per_step: 100,
        maximum_energy_solves: 1000,
        temperature_tolerance_k: 1e-9,
        maximum_storage_defect_w: 1e-8,
        maximum_energy_closure_j: 1e-8,
        temperature_relaxation: 1.0,
    }
}

#[test]
fn g1_g3_nonlinear_material_march_closes_native_storage_and_limited_endpoint_heat() {
    let domain = VoxelDomain::new(6, 1, 1, 0.25).unwrap();
    let reference = FluidTransportState::try_new(1.0, 2.0, 0.02, 0.1).unwrap();
    let fluid = MaterialChtFluid::try_new(&reference).unwrap();
    let mut faces = [FvBoundary::Symmetry; 6];
    faces[0] = FvBoundary::Inlet { velocity: [0.2, 0.0, 0.0] };
    faces[1] = FvBoundary::Outlet;
    let mut flow = SimpleConfig::new(faces);
    flow.pressure_tolerance = 1e-13;
    let mut faces = [ThermalFace::Adiabatic; 6];
    faces[0] = ThermalFace::Inflow { temperature: 300.0 };
    faces[1] = ThermalFace::Outflow { backflow_temperature: 300.0 };
    let thermal = ThermalSetup::new(faces);
    let saved = thermal.clone();
    let initial = [300.0, 302.0, 307.0, 309.0, 308.0, 305.0];
    let mut unsteady = UnsteadyConfig::new(0.1, 3);
    unsteady.scheme = TimeScheme::BackwardEuler;
    unsteady.inner_iterations = 1000;
    unsteady.inner_tolerance = 1e-11;
    let mut bounds = limits([0.0; 3]);
    bounds.max_relative_specific_heat_departure = 0.5;
    let caloric = QuadraticCaloric;
    let run = |scheme, datum| fluid.march_natural_convection_with_enthalpy(
        &domain, &[], &flow, &thermal, &energy(scheme), &initial, bounds,
        &caloric, caloric_controls(datum), &unsteady, |_| 1.0, |_| 1.0,
        None, &CancelGate::new(), |_, t, _| {
            FluidTransportState::try_new(1.0, 2.0 + 0.05 * (t - 300.0), 0.02, 0.1)
        },
    ).unwrap();
    let limited = run(ConvectionScheme::EnthalpyMinmod, 300.0);
    let shifted = run(ConvectionScheme::EnthalpyMinmod, 280.0);
    let upwind = run(ConvectionScheme::Upwind, 300.0);
    let march = &limited.convection.march;
    assert!(march.temperature.iter().zip(&upwind.convection.march.temperature)
        .any(|(a, b)| (a - b).abs() > 1e-5));
    for (actual, other) in march.temperature.iter().zip(&shifted.convection.march.temperature) {
        assert!((actual - other).abs() < 1e-8);
    }
    let stored: f64 = march.temperature.iter().zip(initial).enumerate()
        .map(|(cell, (&t, old))| domain.volume(cell) * (enthalpy(t) - enthalpy(old)))
        .sum();
    let recorded: f64 = march.energy_records.iter().map(|r| r.stored_energy_change_j).sum();
    let outflow: f64 = march.energy_records.iter().map(|r| r.boundary_outflow_j).sum();
    assert!((stored - recorded).abs() < 1e-8);
    assert!((stored + outflow).abs() < 3e-8);
    assert_eq!(limited.records.len(), unsteady.steps);
    assert_eq!(march.flow.records.len(), unsteady.steps);
    for record in &limited.records {
        assert!(record.passes > 1);
        assert!(record.advection.maximum_face_flux_defect_w <= 1e-9);
        assert!(record.advection.maximum_cell_flux_defect_w <= 1e-9);
        assert!(record.maximum_storage_defect_w <= 1e-8);
    }
    for record in &march.energy_records {
        assert!(record.closure_j.abs() < 1e-8);
        assert_eq!(record.source_j, 0.0);
    }
    for (cell, endpoint) in limited.enthalpy.endpoint.iter().enumerate() {
        close(endpoint.unwrap().enthalpy_j_kg, enthalpy(march.temperature[cell]), 1e-10);
    }
    close(limited.convection.endpoint_time_s, 0.3, 1e-13);
    assert_eq!(thermal, saved);
    assert_eq!(initial, [300.0, 302.0, 307.0, 309.0, 308.0, 305.0]);
}

#[test]
fn g0_limited_enthalpy_cannot_silently_select_plain_heat_or_momentum_transport() {
    let domain = VoxelDomain::new(3, 1, 1, 1.0).unwrap();
    let flow = FlowField::from_face_velocity(
        &domain, [FlowFace::Wall; 6], |_, _| 0.0,
    ).unwrap();
    let thermal = ThermalSetup::new([ThermalFace::Temperature(300.0); 6]);
    assert!(matches!(solve_energy(
        &domain, &properties(), &[], &flow, &thermal,
        &energy(ConvectionScheme::EnthalpyMinmod), &CancelGate::new(),
    ), Err(ChtError::InvalidInput { field: "energy.scheme", .. })));
    let mut momentum = SimpleConfig::new([FvBoundary::wall(); 6]);
    momentum.scheme = ConvectionScheme::EnthalpyMinmod;
    assert!(matches!(simple_flow(
        &domain, &properties(), &momentum, &CancelGate::new(),
    ), Err(ChtError::InvalidInput { field: "simple.scheme", .. })));
}
