//! Material refresh retains the original magnetic effort and electrode history.
use super::*;
use fs_material::conductor::ElectricalPropertyBasis;

#[path = "../../winding_ohmic/common.rs"]
pub(super) mod source;

pub(super) fn temperature_sources(
    card: &MaterialCard,
    temperature_k: f64,
    basis: ElectricalPropertyBasis,
) -> (
    Vec<SpatialMagneticMaterial>,
    Vec<Option<ResolvedConductorState>>,
) {
    let material = source::bundle_at(card, &source::absolute_point(temperature_k), basis);
    let conductor =
        ResolvedConductorState::try_from_material_with_basis(material.clone(), basis).unwrap();
    (
        vec![SpatialMagneticMaterial::Linear(material); 6],
        vec![Some(conductor); 6],
    )
}

pub(super) fn temperature_problem(
    card: &MaterialCard,
    basis: ElectricalPropertyBasis,
) -> (
    SpatialVectorMagnetostatics,
    Vec<Option<ResolvedConductorState>>,
) {
    let (materials, sources) = temperature_sources(card, 300.0, basis);
    let (positions, tets) = grid(1);
    (
        SpatialVectorMagnetostatics::new(positions, tets, materials, assembly(), &mut || false)
            .unwrap(),
        sources,
    )
}

fn voltage_state(
    problem: &SpatialVectorMagnetostatics,
    conductors: Vec<Option<ResolvedConductorState>>,
) -> SpatialTerminalInduction<'_> {
    let initial = problem
        .solve(
            &[[0.0; 3]; 6],
            &trace(problem, [0.0; 3], false),
            static_controls(),
            &mut || false,
        )
        .unwrap();
    let field =
        SpatialInduction::new(&initial, conductors, preparation(), 2e-9, &mut || false).unwrap();
    SpatialTerminalInduction::new(field, &electrodes(problem), preparation(), &mut || false)
        .unwrap()
}

#[test]
fn voltage_temperature_refresh_preserves_old_magnetic_effort_and_parameter_work() {
    // Independent one-edge midpoint equation, with BOTH material endpoints:
    // (sigma1/5 + dt*K1/2)*a1 = (sigma1/5 - dt*K0/2)*a0 + dt*sigma1*U/6.
    // Resetting the before-state to the new material gives the wrong K0 term.
    for basis in [
        ElectricalPropertyBasis::Resistivity,
        ElectricalPropertyBasis::Conductivity,
    ] {
        let card = source::absolute_affine_magnetic_card(0.2, 0.002, 2.0, 0.003, basis);
        let (problem, sources) = temperature_problem(&card, basis);
        let mut state = voltage_state(&problem, sources);
        let (dt, u): (f64, f64) = (0.2, 1.1);
        let mut a = 0.0;
        let mut k0 = 2.0;
        for temperature in [300.0, 600.0, 450.0, 350.0] {
            let (materials, conductors) = temperature_sources(&card, temperature, basis);
            let sigma = 1.0 / conductors[0].as_ref().unwrap().resistivity_ohm_m();
            let k1 = 4.0 / (2.0 + 0.003 * (temperature - 300.0));
            let old = a;
            a = ((sigma / 5.0 - 0.5 * dt * k0) * old + dt * sigma * u / 6.0)
                / (sigma / 5.0 + 0.5 * dt * k1);
            let rate = (a - old) / dt;
            let candidate = state
                .step_voltage_with_materials_candidate(
                    &[u, 0.0],
                    &[[0.0; 3]; 6],
                    &[[0.0; 3]; 8],
                    dt,
                    materials,
                    conductors.clone(),
                    assembly(),
                    controls(),
                    &mut || false,
                )
                .unwrap();
            let step = candidate.step();
            near(
                candidate.state().potential_wb()[diagonal(&problem)],
                a,
                3e-9,
            );
            near(step.stored_before_j, 0.5 * k0 * old * old, 3e-9);
            near(step.stored_after_j, 0.5 * k1 * a * a, 3e-9);
            near(step.parameter_work_j, 0.5 * (k1 - k0) * old * a, 3e-9);
            near(
                step.terminals[0].incoming_current_a,
                sigma * (u - rate / 6.0),
                3e-9,
            );
            near(
                step.heat_j,
                dt * sigma * (u * u - u * rate / 3.0 + rate * rate / 5.0),
                3e-9,
            );
            near(step.energy_residual_j, 0.0, 3e-9);
            assert_eq!(candidate.state().conductors(), conductors);
            near(state.potential_wb()[diagonal(&problem)], old, 3e-9);
            state = candidate.into_parts().0;
            k0 = k1;
        }
    }
}

#[test]
fn current_temperature_refresh_solves_voltage_and_matches_independent_inductive_mode() {
    let basis = ElectricalPropertyBasis::Resistivity;
    let card = source::absolute_affine_magnetic_card(0.2, 0.002, 2.0, 0.003, basis);
    let (problem, sources) = temperature_problem(&card, basis);
    let mut voltage = voltage_state(&problem, sources);
    let mut current = voltage
        .clone()
        .into_current_driven(preparation(), &mut || false)
        .unwrap();
    let (dt, incoming): (f64, f64) = (0.125, 1.7);
    let mut a = 0.0;
    let mut k0 = 2.0;
    for temperature in [300.0, 600.0, 450.0] {
        let (materials, conductors) = temperature_sources(&card, temperature, basis);
        let sigma = 1.0 / conductors[0].as_ref().unwrap().resistivity_ohm_m();
        let k1 = 4.0 / (2.0 + 0.003 * (temperature - 300.0));
        let mass = sigma * (1.0 / 5.0 - 1.0 / 36.0);
        let old = a;
        a = ((mass - 0.5 * dt * k0) * old + dt * incoming / 6.0) / (mass + 0.5 * dt * k1);
        let u = incoming / sigma + (a - old) / (6.0 * dt);
        let prescribed = voltage
            .step_voltage_with_materials_candidate(
                &[u, 0.0],
                &[[0.0; 3]; 6],
                &[[0.0; 3]; 8],
                dt,
                materials.clone(),
                conductors.clone(),
                assembly(),
                controls(),
                &mut || false,
            )
            .unwrap();
        let solved = current
            .step_current_with_materials_candidate(
                &[incoming, -incoming],
                &[[0.0; 3]; 6],
                &[[0.0; 3]; 8],
                dt,
                materials,
                conductors,
                assembly(),
                controls(),
                &mut || false,
            )
            .unwrap();
        near(solved.state().potential_wb()[diagonal(&problem)], a, 3e-9);
        near(-solved.step().terminals[1].voltage_to_reference_v, u, 3e-9);
        near(
            solved.step().terminals[0].incoming_current_a,
            incoming,
            3e-9,
        );
        near(
            solved.step().parameter_work_j,
            0.5 * (k1 - k0) * old * a,
            3e-9,
        );
        near(solved.step().stored_before_j, 0.5 * k0 * old * old, 3e-9);
        near(solved.step().heat_j, prescribed.step().heat_j, 3e-9);
        near(solved.step().terminal_work_j, dt * u * incoming, 3e-9);
        for (&a, &b) in solved
            .state()
            .potential_wb()
            .iter()
            .zip(prescribed.state().potential_wb())
        {
            near(a, b, 3e-9);
        }
        voltage = prescribed.into_parts().0;
        current = solved.into_parts().0;
        k0 = k1;
    }
}

#[test]
fn electrode_temperature_refresh_refuses_source_or_support_changes_without_mutation() {
    let basis = ElectricalPropertyBasis::Resistivity;
    let card = source::absolute_affine_magnetic_card(0.2, 0.002, 2.0, 0.003, basis);
    let (problem, sources) = temperature_problem(&card, basis);
    let state = voltage_state(&problem, sources)
        .step_voltage_candidate(
            &[1.0, 0.0],
            &[[0.0; 3]; 6],
            &[[0.0; 3]; 8],
            0.1,
            controls(),
            &mut || false,
        )
        .unwrap()
        .into_parts()
        .0;
    let retained = state.potential_wb().to_vec();
    let current = state
        .clone()
        .into_current_driven(preparation(), &mut || false)
        .unwrap();
    let (materials, conductors) = temperature_sources(&card, 600.0, basis);
    let mut missing = conductors.clone();
    missing[0] = None;
    assert!(
        state
            .step_voltage_with_materials_candidate(
                &[1.0, 0.0],
                &[[0.0; 3]; 6],
                &[[0.0; 3]; 8],
                0.1,
                materials.clone(),
                missing,
                assembly(),
                controls(),
                &mut || false,
            )
            .is_err()
    );
    let changed = source::absolute_affine_magnetic_card(0.3, 0.002, 2.0, 0.003, basis);
    let (changed_materials, changed_conductors) = temperature_sources(&changed, 600.0, basis);
    assert!(
        current
            .step_current_with_materials_candidate(
                &[1.0, -1.0],
                &[[0.0; 3]; 6],
                &[[0.0; 3]; 8],
                0.1,
                changed_materials,
                changed_conductors,
                assembly(),
                controls(),
                &mut || false,
            )
            .is_err()
    );
    assert!(
        current
            .step_current_with_materials_candidate(
                &[1.0, -1.0],
                &[[0.0; 3]; 6],
                &[[0.0; 3]; 8],
                0.1,
                materials,
                conductors,
                assembly(),
                controls(),
                &mut || true,
            )
            .is_err()
    );
    assert_eq!(state.potential_wb(), retained);
    assert_eq!(current.potential_wb(), retained);
    assert_eq!(state.time_s(), 0.1);
    assert_eq!(current.time_s(), 0.1);
}
