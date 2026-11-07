//! G1/G3: timed restrained thermal work with independent intact/rupture oracles.
use super::{ResolvedThermalBrittle, card, near, properties, query, scalar, settings};
use fs_matdb::MaterialCard;
use fs_material::brittle::source::expansion::rate::resolve_thermal_viscous_brittle;
use fs_material::brittle::source::expansion::resolve_thermal_brittle;
use fs_material::brittle::source::rate::{
    BRITTLE_NORMAL_VISCOSITY_PROPERTY, BRITTLE_TANGENTIAL_VISCOSITY_PROPERTY,
    BRITTLE_VISCOSITY_DIMS,
};
use fs_material::state_point::{
    INVERSE_TEMPERATURE_DIMS, LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
    MaterialPropertySelection,
};

const ETA: f64 = 1e8;
const WIDTH: f64 = 0.01;

fn thermal_card(alpha: f64) -> MaterialCard {
    let mut claims = properties();
    claims.push(scalar(
        LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
        INVERSE_TEMPERATURE_DIMS,
        alpha,
    ));
    for property in [
        BRITTLE_NORMAL_VISCOSITY_PROPERTY,
        BRITTLE_TANGENTIAL_VISCOSITY_PROPERTY,
    ] {
        claims.push(scalar(property, BRITTLE_VISCOSITY_DIMS, ETA));
    }
    card(claims)
}

fn source(card: &MaterialCard, reference: f64, current: f64) -> ResolvedThermalBrittle {
    resolve_thermal_viscous_brittle(
        card,
        &query(reference),
        &query(current),
        MaterialPropertySelection::SingleClaimOnly,
        settings(WIDTH),
        1.0,
    )
    .unwrap()
}

#[test]
fn restrained_intact_cooling_uses_total_mean_stress_and_inverse_duration_loss() {
    let card = thermal_card(1e-5);
    let start = source(&card, 400.0, 400.0);
    let cold = start.at_temperature(&card, 399.0).unwrap();
    let virgin = start.band().law().initial_state();
    let epsilon = 1e-5;
    // E=1e9, nu=.2: hydrostatic elastic stress E*epsilon/(1-2nu),
    // transverse normal mode [1,nu/(1-nu),nu/(1-nu)] = [1,.25,.25].
    let hydrostatic = 1e9 * epsilon / 0.6;
    let mode = 1.5 * epsilon;
    let storage = 1.5 * hydrostatic * epsilon;
    for dt in [0.02, 0.2, 2.0] {
        let step = start
            .transfer_temperature_over_time(&virgin, &cold, dt, 1e-6)
            .unwrap();
        let viscous = ETA * WIDTH * mode / dt;
        let loss = ETA * WIDTH * mode * mode / dt;
        near(step.path.endpoint.stress_pa[0], hydrostatic + viscous);
        near(
            step.path.endpoint.stress_pa[1],
            hydrostatic + 0.25 * viscous,
        );
        near(step.path.mean_stress_pa[0], 0.5 * hydrostatic + viscous);
        near(step.path.endpoint.viscous_dissipated_j_m3, loss);
        near(step.stored_change_j_m3, storage);
        near(step.eigenstrain_work_j_m3, storage + loss);
        near(step.coefficient.parameter_work_j_m3, 0.0);
        near(step.fracture_dissipated_j_m3, 0.0);
        near(step.balance_defect_j_m3, 0.0);
        near(step.eigenstrain_work_defect_j_m3, 0.0);
        let end = &step.path.endpoint.equilibrium.state;
        for value in cold.total_strain(end).unwrap() {
            near(value, 0.0);
        }
        let hold = cold
            .transfer_temperature_over_time(end, &cold, dt, 1e-6)
            .unwrap();
        assert_eq!(hold.path.endpoint.viscous_dissipated_j_m3, 0.0);
        near(hold.eigenstrain_work_j_m3, 0.0);
        near(hold.path.endpoint.stress_pa[0], hydrostatic);
    }
    assert_eq!(virgin.strain(), &[0.0; 6]);
    assert_eq!(virgin.damage(), 0.0);
}

#[test]
fn restrained_cooling_fractures_once_and_reheating_cannot_restore_dashpots() {
    let card = thermal_card(1e-4);
    let start = source(&card, 400.0, 400.0);
    let cold = start.at_temperature(&card, 200.0).unwrap();
    let virgin = start.band().law().initial_state();
    let dt = 1.0;
    let step = start
        .transfer_temperature_over_time(&virgin, &cold, dt, 1e-6)
        .unwrap();
    let state = &step.path.endpoint.equilibrium.state;
    assert!(state.failed());
    let normal_modulus = 1e9 * 0.8 / (1.2 * 0.6);
    let onset_jump = WIDTH * 1e6 / normal_modulus;
    let failure_jump = 2.0 * 100.0 / 1e6;
    let end_jump = WIDTH * 1.5 * 0.02;
    // Independent integral of bilinear surviving-bond integrity over opening.
    let integral =
        onset_jump * failure_jump / (failure_jump - onset_jump) * (failure_jump / onset_jump).ln();
    let loss = ETA / WIDTH * (end_jump / dt) * integral;
    near(step.path.endpoint.viscous_dissipated_j_m3, loss);
    near(step.fracture_dissipated_j_m3, 100.0 / WIDTH);
    near(
        step.eigenstrain_work_j_m3,
        step.stored_change_j_m3 + 100.0 / WIDTH + loss,
    );
    near(step.balance_defect_j_m3, 0.0);
    near(step.path.endpoint.stress_pa[0], 0.0);
    for value in cold.total_strain(state).unwrap() {
        near(value, 0.0);
    }

    // This rate exceeds the former adhesive bound by many orders of magnitude;
    // accepted complete failure has no remaining adhesive dashpot or rate cap.
    let reheated = cold
        .transfer_temperature_over_time(state, &start, 1e-12, 1e-6)
        .unwrap();
    assert!(reheated.path.endpoint.equilibrium.state.failed());
    assert_eq!(reheated.path.endpoint.viscous_dissipated_j_m3, 0.0);
    near(reheated.fracture_dissipated_j_m3, 0.0);
    near(
        reheated.path.endpoint.equilibrium.state.dissipated_j_m3(),
        100.0 / WIDTH,
    );
    near(reheated.eigenstrain_work_j_m3, -step.stored_change_j_m3);
    near(reheated.balance_defect_j_m3, 0.0);
    assert_eq!(virgin.damage(), 0.0);
}

#[test]
fn timed_restrained_transfer_refuses_invalid_paths_and_preserves_intrinsic_behavior() {
    let card = thermal_card(1e-5);
    let start = source(&card, 400.0, 400.0);
    let cold = start.at_temperature(&card, 300.0).unwrap();
    let virgin = start.band().law().initial_state();
    let saved = virgin;
    for dt in [0.0, -1.0, f64::NAN, f64::INFINITY, 1e-12] {
        assert!(
            start
                .transfer_temperature_over_time(&virgin, &cold, dt, 1e-6)
                .is_err()
        );
    }
    let foreign_reference = source(&card, 350.0, 300.0);
    assert!(
        start
            .transfer_temperature_over_time(&virgin, &foreign_reference, 1.0, 1e-6)
            .is_err()
    );
    assert!(
        start
            .transfer_temperature_over_time(&virgin, &cold, 1.0, 0.0)
            .is_err()
    );
    assert_eq!(virgin, saved);

    let intrinsic = resolve_thermal_brittle(
        &card,
        &query(400.0),
        &query(400.0),
        MaterialPropertySelection::SingleClaimOnly,
        settings(WIDTH),
    )
    .unwrap();
    let next = intrinsic.at_temperature(&card, 300.0).unwrap();
    let old = intrinsic.band().law().initial_state();
    let untimed = intrinsic.transfer_temperature(&old, &next, 1e-6).unwrap();
    let timed = intrinsic
        .transfer_temperature_over_time(&old, &next, 0.03, 1e-6)
        .unwrap();
    assert_eq!(timed.path.endpoint.equilibrium.state, untimed.state);
    assert_eq!(timed.path.endpoint.stress_pa, untimed.stress_pa);
    assert_eq!(timed.path.mean_stress_pa, untimed.mean_stress_pa);
    assert_eq!(timed.path.endpoint.viscous_dissipated_j_m3, 0.0);
    near(
        timed.fracture_dissipated_j_m3,
        untimed.damage_dissipated_j_m3,
    );
    near(timed.eigenstrain_work_j_m3, untimed.eigenstrain_work_j_m3);
}
