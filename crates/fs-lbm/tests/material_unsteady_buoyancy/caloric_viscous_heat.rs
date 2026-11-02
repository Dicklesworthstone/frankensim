//! Physical molecular work and original source calorimetry share one clock.
//! Synthetic laws below are independent numerical oracles, not material claims.

use super::*;
use fs_lbm::conjugate::{AffineSpecificEnthalpy, ViscousHeatControls};
use fs_lbm::material::conjugate::buoyant::enthalpy::ConvectionCaloricLaw;

struct Caloric;
impl ConvectionCaloricLaw for Caloric {
    fn tangent_at_temperature(&self, reference: f64, t: f64)
        -> Result<AffineSpecificEnthalpy, MaterialError> {
        Ok(AffineSpecificEnthalpy {
            temperature_k: t,
            enthalpy_j_kg: (t - reference) * (2.0 + 0.025 * (t + reference - 600.0)),
            heat_capacity_j_kg_k: 2.0 + 0.05 * (t - 300.0),
        })
    }
}
fn h(t: f64) -> f64 {
    let theta = t - 300.0;
    2.0 * theta + 0.025 * theta * theta
}
fn mu(t: f64) -> f64 { 1.0 / (1.0 + 0.2 * (t - 300.0)) }
fn source(t: f64) -> Result<FluidTransportState, MaterialError> {
    FluidTransportState::try_new(1.0, 2.0 + 0.05 * (t - 300.0), 0.01, mu(t))
}
fn heat_controls() -> ViscousHeatControls {
    ViscousHeatControls { maximum_reynolds: 0.1, power_tolerance_w: 1e-8 }
}
fn caloric_bounds() -> BoussinesqConfig<3> {
    let mut result = bounds();
    result.max_relative_specific_heat_departure = 0.5;
    result
}
fn wall_flow(speed: f64) -> SimpleConfig {
    let mut faces = [FvBoundary::wall(); 6];
    faces[3] = FvBoundary::Wall { velocity: [speed, 0.0, 0.0] };
    faces[4] = FvBoundary::Symmetry;
    faces[5] = FvBoundary::Symmetry;
    SimpleConfig::new(faces)
}
fn clock(dt: f64, steps: usize) -> UnsteadyConfig {
    let mut unsteady = UnsteadyConfig::new(dt, steps);
    unsteady.scheme = TimeScheme::BackwardEuler;
    unsteady.inner_iterations = 1000;
    unsteady.inner_tolerance = 1e-11;
    unsteady
}

#[test]
fn moving_wall_heat_follows_original_cp_and_mu_without_compounding() {
    // Both closed MAC cells have zero volume flow. The one free tangential
    // wall stencil has area dx*dz=.0625 and half-height=.125. Each cell's
    // parallel half therefore receives (0.0625/.125)/2 * mu_i * U_wall^2.
    let domain = VoxelDomain::new(2, 1, 1, 0.25).unwrap();
    let reference = source(300.0).unwrap();
    let bound = MaterialChtFluid::try_new(&reference).unwrap();
    let flow = wall_flow(0.1);
    let thermal = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    let initial = [300.0; 2];
    let unsteady = clock(5.0, 4);
    let gate = CancelGate::new();
    let run = |setup: &ThermalSetup| bound.march_natural_convection_with_enthalpy_and_viscous_heat(
        &domain, &[], &flow, setup, &energy(), &initial, caloric_bounds(),
        &Caloric, controls(), &unsteady, heat_controls(), |_| 1.0, |_| 0.0,
        None, &gate, |_, t, _| source(t),
    );
    let result = run(&thermal).unwrap();
    let march = &result.thermal.convection.march;
    let mut expected = 300.0;
    let mut specific_energy = 0.0;
    for (receipt, thermal_record) in result.mechanical_steps.iter().zip(&march.energy_records) {
        let per_cell_w = 0.25 * mu(expected) * 0.1_f64.powi(2);
        specific_energy += per_cell_w * unsteady.time_step_s / domain.volume(0);
        // Stable independent inverse of h=2*theta+.025*theta^2.
        expected = 300.0 + 2.0 * specific_energy / (2.0 + (4.0 + 0.1 * specific_energy).sqrt());
        for &heat in &receipt.power.cell_heat_w { close(heat, per_cell_w, 1e-10); }
        close(receipt.power.heat_w, 2.0 * per_cell_w, 1e-10);
        close(thermal_record.source_j, 2.0 * per_cell_w * unsteady.time_step_s, 1e-9);
        close(thermal_record.max_temperature_k, expected, 1e-8);
        close(thermal_record.boundary_outflow_j, 0.0, 1e-9);
    }
    for &t in &march.temperature { close(t, expected, 1e-8); }
    let stored = march.temperature.iter().enumerate().map(|(c, &t)| domain.volume(c) * h(t)).sum();
    close(stored, result.heating.generated_heat_j, 1e-8);
    close(result.heating.molecular_heat_j, result.heating.generated_heat_j, 1e-10);
    close(result.heating.supplied_wall_work_j, result.heating.generated_heat_j, 1e-9);
    assert_eq!(result.heating.substeps, unsteady.steps);
    close(result.heating.duration_s, 20.0, 1e-14);
    assert!(result.thermal.records.iter().all(|r| r.passes > 1));

    // Physical dissipation is independent of the authored source multiplier.
    let mut scheduled = thermal.clone();
    scheduled.power_w = vec![100.0; 2];
    let rescheduled = run(&scheduled).unwrap();
    assert_eq!(rescheduled.thermal.convection.march.temperature, march.temperature);
    assert_eq!(rescheduled.heating, result.heating);
    let unheated = bound.march_natural_convection_with_enthalpy(
        &domain, &[], &flow, &thermal, &energy(), &initial, caloric_bounds(),
        &Caloric, controls(), &unsteady, |_| 1.0, |_| 0.0, None, &gate,
        |_, t, _| source(t),
    ).unwrap();
    assert_eq!(unheated.convection.march.temperature, initial);

    // The physical owner refuses unsupported history and source endpoints;
    // no accepted march or work receipt escapes either failure.
    let mut bdf = unsteady;
    bdf.scheme = TimeScheme::Bdf2;
    assert!(bound.march_natural_convection_with_enthalpy_and_viscous_heat(
        &domain, &[], &flow, &thermal, &energy(), &initial, caloric_bounds(),
        &Caloric, controls(), &bdf, heat_controls(), |_| 1.0, |_| 1.0,
        None, &gate, |_, t, _| source(t),
    ).is_err());
    assert!(bound.march_natural_convection_with_enthalpy_and_viscous_heat(
        &domain, &[], &flow, &thermal, &energy(), &initial, caloric_bounds(),
        &Caloric, controls(), &unsteady, heat_controls(), |_| 1.0, |_| 1.0,
        None, &gate, |_, t, _| {
            if t > 300.0 + 1e-6 { Err(MaterialError::State { what: "endpoint refused".into() }) }
            else { source(t) }
        },
    ).is_err());
    assert_eq!(initial, [300.0; 2]);
    assert_eq!(bound.transport(), &reference);
    assert!(thermal.power_w.is_empty());
}

#[test]
fn molecular_heat_at_thermostatted_cells_is_rejected_once_with_local_mu() {
    let domain = VoxelDomain::new(2, 1, 1, 0.25).unwrap();
    let reference = source(300.0).unwrap();
    let bound = MaterialChtFluid::try_new(&reference).unwrap();
    let initial = [300.0, 304.0];
    let mut thermal = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    thermal.fixed_temperature = vec![(0, 300.0), (0, 300.0), (1, 304.0)];
    let unsteady = clock(0.5, 3);
    let result = bound.march_natural_convection_with_enthalpy_and_viscous_heat(
        &domain, &[], &wall_flow(0.1), &thermal, &energy(), &initial, caloric_bounds(),
        &Caloric, controls(), &unsteady, heat_controls(), |_| 1.0, |_| 1.0,
        None, &CancelGate::new(), |_, t, _| source(t),
    ).unwrap();
    let expected = initial.map(|t| 0.25 * mu(t) * 0.1_f64.powi(2));
    for (mechanical, heat) in result.mechanical_steps.iter()
        .zip(&result.thermal.convection.march.energy_records) {
        for (&actual, &expected) in mechanical.power.cell_heat_w.iter().zip(&expected) {
            close(actual, expected, 1e-10);
        }
        close(heat.stored_energy_change_j, 0.0, 1e-9);
        close(heat.source_j, expected.iter().sum::<f64>() * unsteady.time_step_s, 1e-9);
        close(heat.boundary_outflow_j, heat.source_j, 1e-9);
        close(heat.closure_j, 0.0, 1e-9);
    }
    assert_eq!(result.thermal.convection.march.temperature, initial);
}

#[test]
fn circulating_material_heat_retains_native_kinetic_history_and_excludes_time_loss() {
    let domain = VoxelDomain::new(4, 4, 1, 0.25).unwrap();
    let reference = source(300.0).unwrap();
    let bound = MaterialChtFluid::try_new(&reference).unwrap();
    let initial = vec![300.0; domain.cell_count()];
    let thermal = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    let unsteady = clock(0.1, 4);
    let result = bound.march_natural_convection_with_enthalpy_and_viscous_heat(
        &domain, &[], &wall_flow(0.025), &thermal, &energy(), &initial, caloric_bounds(),
        &Caloric, controls(), &unsteady, heat_controls(), |_| 1.0, |_| 1.0,
        None, &CancelGate::new(), |_, t, _| source(t),
    ).unwrap();
    assert!(result.heating.final_kinetic_energy_j > 1e-9);
    assert!(result.heating.backward_euler_loss_j > 1e-9);
    for pair in result.mechanical_steps.windows(2) {
        assert_eq!(pair[0].final_kinetic_energy_j.to_bits(), pair[1].initial_kinetic_energy_j.to_bits());
        assert_eq!(pair[0].end_s.to_bits(), pair[1].start_s.to_bits());
    }
    let march = &result.thermal.convection.march;
    let actual_storage: f64 = march.temperature.iter().enumerate()
        .map(|(cell, &t)| domain.volume(cell) * h(t)).sum();
    close(actual_storage, result.heating.generated_heat_j, 1e-8);
    close(march.energy_records.iter().map(|r| r.source_j).sum(), result.heating.generated_heat_j, 1e-10);
    let kinetic_change = result.heating.final_kinetic_energy_j - result.heating.initial_kinetic_energy_j;
    let supply = result.heating.supplied_wall_work_j + result.heating.supplied_body_work_j
        + result.heating.pressure_constraint_work_j;
    close(supply, actual_storage + kinetic_change + result.heating.backward_euler_loss_j, 1e-8);
    assert!((actual_storage - result.heating.generated_heat_j).abs()
        < 0.05 * result.heating.backward_euler_loss_j,
        "numerical BE loss must remain outside source enthalpy");
}

#[test]
fn molecular_heat_uses_the_existing_joint_fluid_solid_caloric_owner() {
    let fluid_card = card(0.004, 0.0);
    let fluid_source = law(&fluid_card, [290.0, 400.0]);
    let caloric = source_caloric(&fluid_card, &fluid_source, [290.0, 400.0]);
    let transport = |t| {
        let state = super::transport(&fluid_source, t)?;
        FluidTransportState::try_new(state.density_kg_m3(), state.specific_heat_j_kg_k(),
            state.thermal_conductivity_w_m_k(), mu(t))
    };
    let reference = transport(300.0).unwrap();
    let bound = MaterialChtFluid::try_new(&reference).unwrap();
    let solid_card = card(0.2, 0.01);
    let solid = law(&solid_card, [290.0, 400.0]);
    let (domain, _, _, _, _) = shear_cavity(true);
    let initial = vec![300.0; domain.cell_count()];
    let thermal = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    let result = bound.march_natural_convection_with_enthalpy_and_material_solids_and_viscous_heat(
        &domain, std::slice::from_ref(&solid), &wall_flow(0.01), &thermal, &energy(), &initial,
        bounds(), &caloric, controls(), SolidThermalConfig::default(), &clock(0.25, 3),
        heat_controls(), |_| 1.0, |_| 1.0, None, &CancelGate::new(), |_, t, _| transport(t),
    ).unwrap();
    let march = &result.thermal.fluid.convection.march;
    let actual_storage: f64 = march.temperature.iter().enumerate().map(|(cell, &t)| {
        let slope = if domain.is_fluid(cell) { 0.004 } else { 0.2 };
        3.0 * domain.volume(cell) * (specific_energy(t, slope) - specific_energy(300.0, slope))
    }).sum();
    close(actual_storage, result.heating.generated_heat_j, 1e-8);
    assert!(result.thermal.solid_records.iter().map(|r| r.solid_stored_energy_change_j).sum::<f64>() > 0.0);
    for step in &result.mechanical_steps {
        for (cell, &heat) in step.power.cell_heat_w.iter().enumerate() {
            if !domain.is_fluid(cell) { assert_eq!(heat, 0.0); }
        }
    }
}
