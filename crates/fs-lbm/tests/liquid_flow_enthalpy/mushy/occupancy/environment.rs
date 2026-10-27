//! Finite atmosphere energy reaches original phase topology and native momentum.
use super::*;
use fs_lbm::conjugate::enthalpy::advection::feedback::environment::{
    FiniteGasContact, FiniteGasFlowError, advance_finite_gas_with_material_flow,
};
use fs_material::thermochemical_gas::energy::{GasEnergySolve, IsochoricGasParcel};

const R: f64 = 8.314_462_618_153_24;

#[path = "../../../nsf3_thermochemical/source.rs"]
mod gas_source;

fn gas_controls() -> GasEnergySolve {
    GasEnergySolve {
        minimum_temperature_k: 200.0,
        maximum_temperature_k: 690.0,
        energy_tolerance_j_kg: 1e-7,
        maximum_iterations: 128,
    }
}

fn contact(id: u64) -> FiniteGasContact {
    FiniteGasContact {
        cell: CellId(id),
        area_m2: 1.0,
        distances_m: [0.1; 2],
    }
}

#[test]
fn finite_gas_cools_by_nasa_internal_energy_while_its_heat_melts_the_original_cells() {
    let card = source(100.0);
    let curve = chart(&card);
    let law = MushyPhaseLaw::try_new_with_immobile_solid(&curve, 1.0, 1.0, 1.0).unwrap();
    let domain = VoxelDomain::new(3, 3, 1, 1.0).unwrap();
    let ids = std::array::from_fn::<_, 9, _>(|c| Some(CellId(c as u64)));
    let gate = CancelGate::new();
    let config = closed();
    let old = geometry(&domain, &curve, [-20.0; 9]);
    let flow = MaterialFlowFeedback::try_new_phase(
        &old,
        &domain,
        input(&ids, &curve),
        &config,
        (),
        &gate,
        4,
        &law,
        |_, cell: &DepositedCell, _| resolve(&card, &law, cell.specific_h_j_kg),
    )
    .unwrap()
    .with_unsteady_momentum()
    .unwrap()
    .with_unsteady_viscous_heating(controls())
    .unwrap();
    // Nonzero formation reference must not become transferred sensible heat.
    let gas = IsochoricGasParcel::new(
        gas_source::model(12000.0, 250000.0, 100.0),
        gas_source::point(600.0, 100000.0),
        0.1,
    )
    .unwrap();
    let contacts: Vec<_> = (0..9).map(contact).collect();
    let dt = 0.25;
    let (next, accepted, cooled, report) = advance_finite_gas_with_material_flow(
        &old,
        dt,
        step(1.0),
        &flow,
        &gas,
        &contacts,
        gas_controls(),
    )
    .unwrap();
    assert_eq!(report.substeps, 1);
    assert_eq!(accepted.momentum_solves(), 0); // Initial exact solid has no invented mu.
    // Independent two-layer resistance and NASA polynomial: cp/R=3+0.001T,
    // so u/R=2T+0.0005T^2+a7. cv, not cp, owns this fixed-volume gas.
    let kg = 2.6 * (0.5 + 99000.0 / 249000.0);
    let g = 1.0 / (0.1 / 0.1 + 0.1 / kg);
    let q_cell = dt * g * (600.0 - 490.0);
    let q = 9.0 * q_cell;
    let target = 2.0 * 600.0 + 0.0005 * 600.0 * 600.0 - q / gas.mass_kg() * 0.028 / R;
    let expected_t = ((4.0 + 0.002 * target).sqrt() - 2.0) / 0.001;
    near(cooled.state().state().temperature, expected_t);
    near(
        cooled.state().state().pressure,
        100000.0 * expected_t / 600.0,
    );
    assert_eq!(cooled.mass_kg(), gas.mass_kg());
    assert_eq!(cooled.volume_m3(), gas.volume_m3());
    near(report.heat_into_condensed_j, q);
    near(report.gas_internal_energy_change_j, -q);
    near(report.thermal.external_heat_j, q);
    near(report.other_thermal_external_heat_j, 0.0);
    assert!(report.gas_balance_defect_j.abs() < 1e-7);
    assert!(report.combined_thermal_balance_defect_j.abs() < 1e-7);
    near(next.energy_j() - old.energy_j(), q);
    for (index, sample) in accepted.endpoint().iter().enumerate() {
        near(sample.specific_enthalpy_j_kg, -20.0 + q_cell);
        near(sample.temperature_k, 500.0); // Latent heat changes phi, not T.
        assert!(sample.transport.mobile_transport().is_some());
        let phi = (-20.0 + q_cell) / 100.0;
        near(
            report.endpoint[index].condensed_conductivity_w_m_k,
            0.1 * (1.0 + phi),
        );
        near(
            report.endpoint[index].film.state().pressure,
            cooled.state().state().pressure,
        );
        near(
            report.endpoint[index].film.state().temperature,
            0.5 * (500.0 + expected_t),
        );
    }
    assert_eq!(accepted.unsteady_time_s(), Some(dt));
    assert_eq!(
        accepted.unsteady_heating_audit().unwrap().generated_heat_j,
        0.0
    );
    assert_eq!(flow.unsteady_time_s(), Some(0.0));
    assert_eq!(old.cells()[&CellId(0)].specific_h_j_kg, -20.0);
    assert_eq!(gas.state().state().temperature, 600.0);
}

#[test]
fn changing_only_gas_transport_data_opens_phase_faces_and_changes_native_momentum() {
    let card = source(100.0);
    let curve = chart(&card);
    let law = MushyPhaseLaw::try_new_with_immobile_solid(&curve, 1.0, 1.0, 1.0).unwrap();
    let domain = VoxelDomain::new(3, 3, 1, 1.0).unwrap();
    let ids = std::array::from_fn::<_, 9, _>(|c| Some(CellId(c as u64)));
    let gate = CancelGate::new();
    let config = closed();
    let mut enthalpy = [50.0; 9];
    enthalpy[4] = -20.0;
    let old = geometry(&domain, &curve, enthalpy);
    let flow = MaterialFlowFeedback::try_new_phase(
        &old,
        &domain,
        input(&ids, &curve),
        &config,
        (),
        &gate,
        20,
        &law,
        |_, cell: &DepositedCell, _| resolve(&card, &law, cell.specific_h_j_kg),
    )
    .unwrap()
    .with_unsteady_body_force_heating(corner_force(), controls())
    .unwrap();
    let mut outputs = Vec::new();
    for scale in [1.0, 100.0] {
        let gas = IsochoricGasParcel::new(
            gas_source::model(0.0, 250000.0, scale),
            gas_source::point(600.0, 100000.0),
            1.0,
        )
        .unwrap();
        let (next, accepted, cooled, report) = advance_finite_gas_with_material_flow(
            &old,
            0.3,
            step(0.3),
            &flow,
            &gas,
            &[contact(4)],
            gas_controls(),
        )
        .unwrap();
        let mechanics = accepted.unsteady_heating_audit().unwrap();
        near(
            report.other_thermal_external_heat_j,
            mechanics.generated_heat_j,
        );
        assert!(report.combined_thermal_balance_defect_j.abs() < 1e-7);
        assert!(cooled.state().state().temperature < 600.0);
        assert_eq!(accepted.unsteady_time_s(), Some(next.time_s()));
        let opened = accepted.endpoint()[4]
            .transport
            .mobile_transport()
            .is_some();
        let native = &accepted.unsteady_flow().unwrap().field;
        let internal_flux = [
            native.flux_x(1, 1, 0),
            native.flux_x(2, 1, 0),
            native.flux_y(1, 1, 0),
            native.flux_y(1, 2, 0),
        ]
        .into_iter()
        .map(f64::abs)
        .fold(0.0, f64::max);
        outputs.push((
            opened,
            internal_flux,
            next.cells()[&CellId(4)].specific_h_j_kg,
        ));
    }
    assert!(!outputs[0].0);
    assert_eq!(outputs[0].1, 0.0);
    assert!(outputs[1].0);
    assert!(outputs[1].1 > 1e-10);
    assert!(outputs[1].2 > outputs[0].2);
    assert_eq!(flow.unsteady_time_s(), Some(0.0));
}

#[test]
fn late_gas_pressure_domain_failure_leaves_thermal_momentum_and_reservoir_unchanged() {
    let card = source(100.0);
    let curve = chart(&card);
    let law = MushyPhaseLaw::try_new_with_immobile_solid(&curve, 1.0, 1.0, 1.0).unwrap();
    let domain = VoxelDomain::new(3, 3, 1, 1.0).unwrap();
    let ids = std::array::from_fn::<_, 9, _>(|c| Some(CellId(c as u64)));
    let gate = CancelGate::new();
    let config = closed();
    let old = geometry(&domain, &curve, [300.0; 9]);
    let flow = MaterialFlowFeedback::try_new_phase(
        &old,
        &domain,
        input(&ids, &curve),
        &config,
        (),
        &gate,
        20,
        &law,
        |_, cell: &DepositedCell, _| resolve(&card, &law, cell.specific_h_j_kg),
    )
    .unwrap()
    .with_unsteady_momentum()
    .unwrap();
    let gas = IsochoricGasParcel::new(
        gas_source::model(0.0, 250000.0, 100.0),
        gas_source::point(300.0, 249000.0),
        1e-5,
    )
    .unwrap();
    let saved = gas.clone();
    let result = advance_finite_gas_with_material_flow(
        &old,
        0.01,
        step(0.01),
        &flow,
        &gas,
        &[contact(0)],
        gas_controls(),
    );
    assert!(matches!(result, Err(FiniteGasFlowError::Gas(_))));
    assert_eq!(gas, saved);
    assert_eq!(old.time_s(), 0.0);
    assert!(
        old.cells()
            .values()
            .all(|cell| cell.specific_h_j_kg == 300.0)
    );
    assert_eq!(flow.unsteady_time_s(), Some(0.0));
    assert!(flow.unsteady_step().is_none());
    gate.request();
    let cancelled = advance_finite_gas_with_material_flow(
        &old,
        0.01,
        step(0.01),
        &flow,
        &gas,
        &[contact(0)],
        gas_controls(),
    );
    assert!(matches!(
        cancelled,
        Err(FiniteGasFlowError::Flow(MaterialFlowError::Thermal(
            DepositionError::Cancelled
        )))
    ));
    assert_eq!(gas, saved);
}

#[test]
fn combined_finite_capacity_bound_prevents_thermal_overshoot_with_a_tiny_reservoir() {
    let card = source(100.0);
    let curve = chart(&card);
    let law = MushyPhaseLaw::try_new_with_immobile_solid(&curve, 1.0, 1.0, 1.0).unwrap();
    let domain = VoxelDomain::new(3, 3, 1, 1.0).unwrap();
    let ids = std::array::from_fn::<_, 9, _>(|c| Some(CellId(c as u64)));
    let gate = CancelGate::new();
    let config = closed();
    let old = geometry(&domain, &curve, [-20.0; 9]);
    let flow = MaterialFlowFeedback::try_new_phase(
        &old,
        &domain,
        input(&ids, &curve),
        &config,
        (),
        &gate,
        20,
        &law,
        |_, cell: &DepositedCell, _| resolve(&card, &law, cell.specific_h_j_kg),
    )
    .unwrap()
    .with_unsteady_momentum()
    .unwrap();
    let gas = IsochoricGasParcel::new(
        gas_source::model(0.0, 250000.0, 100.0),
        gas_source::point(600.0, 100000.0),
        1e-5,
    )
    .unwrap();
    let (next, _, cooled, report) = advance_finite_gas_with_material_flow(
        &old,
        0.01,
        step(0.01),
        &flow,
        &gas,
        &[contact(0)],
        gas_controls(),
    )
    .unwrap();
    assert!(report.substeps > 1);
    let t = next.temperature_k(CellId(0)).unwrap();
    assert!(t >= 490.0 && t < 500.0);
    assert!(cooled.state().state().temperature >= t);
    assert!(cooled.state().state().temperature < 600.0);
    assert!(report.combined_thermal_balance_defect_j.abs() < 1e-7);
}
