//! Native source-conjugate reservoir evolution and independent bath accounting.
use super::*;
use fs_thermochem::{MassSpecificEnergyQuantityV1, SpeciesId};

// Rows are [Cp/R constant, Cp/R linear coefficient, M/R]. The entropy and
// bath chemical-potential formulas below integrate those polynomials directly.
const T0: f64 = 0.001;
const COUNTER: [[f64; 3]; 2] = [[3.5, 0.0, 1.0], [6.5, 0.0, 1.0]];
const VARIABLE: [[f64; 3]; 2] = [[3.5, 0.2, 1.0], [4.0, 0.4, 2.0]];

fn pressure(rows: [[f64; 3]; 2], rho: [f64; 2], t: f64) -> f64 {
    t * (rho[0] / rows[0][2] + rho[1] / rows[1][2])
}

fn fixture(
    rows: [[f64; 3]; 2],
    formation: f64,
    rho: [f64; 2],
    velocity: f64,
    temperature: f64,
) -> MixtureEuler1d {
    let g = MixtureGas::new(
        bank(vec![
            source("A", rows[0][2] * R, rows[0][0], rows[0][1], formation),
            source(
                "B",
                rows[1][2] * R,
                rows[1][0],
                rows[1][1],
                -2.0 * formation,
            ),
        ]),
        CaloricDomain::new(Temperature::new(T0), Temperature::new(50.0)),
        || false,
    )
    .unwrap();
    let states = vec![
        g.conserved(&rho, velocity, pressure(rows, rho, temperature))
            .unwrap();
        2
    ];
    // Two actual cells, each one metre wide: the bath distance is 0.5 m.
    MixtureEuler1d::new(
        g,
        Length::new(1.0),
        0.7,
        states,
        Boundaries::Ends {
            left: Boundary::Transmissive,
            right: Boundary::Transmissive,
        },
        || false,
    )
    .unwrap()
}

fn physical_entropy(flow: &MixtureEuler1d, rows: [[f64; 3]; 2]) -> f64 {
    (0..flow.cell_count())
        .map(|i| {
            let c = flow.cell(i).unwrap();
            let t = flow.gas().primitive(c, || false).unwrap().temperature;
            flow.mesh().volumes()[i]
                * rows
                    .iter()
                    .enumerate()
                    .map(|(k, &[a, b, mass])| {
                        let rho = c.partial_densities[k];
                        rho * ((a - 1.0) * t.ln() + b * t - (rho / mass).ln()) / mass
                    })
                    .sum::<f64>()
        })
        .sum()
}

fn bath_input(rows: [[f64; 3]; 2], rho: [f64; 2], t: f64, transfer: Totals) -> f64 {
    let chemical: f64 = rows
        .iter()
        .enumerate()
        .map(|(k, &[a, b, mass])| {
            let h = (a * t + 0.5 * b * t * t - (a - 1.0) * T0 - 0.5 * b * T0 * T0) / mass;
            let s = ((a - 1.0) * t.ln() + b * t - (rho[k] / mass).ln()) / mass;
            (h / t - s) * transfer.species_masses[k]
        })
        .sum();
    transfer.energy / t - chemical
}

#[test]
fn g1_hot_composition_bath_repairs_negative_production_in_both_orientations() {
    for end in 0..2 {
        let mut gauge_reference = None;
        for formation in [0.0, 1e18] {
            let mut flow = fixture(COUNTER, formation, [0.1, 0.9], 0.0, 1.0);
            let bath = flow.gas().conserved(&[0.9, 0.1], 0.0, 10.0).unwrap();
            let mut baths = [None; 2];
            baths[end] = Some(bath.as_cell());
            let original = flow.states.clone();
            let before = physical_entropy(&flow, COUNTER);
            let controls = DissipationControl {
                maximum_step_s: 1e-4,
                maximum_halvings: 0,
                ..Default::default()
            };
            let entropy_controls = MolecularEntropyControl {
                absolute_tolerance_j_k: 1e-14,
                ..Default::default()
            };
            if formation == 0.0 {
                // The engineering bath-enthalpy face has initial generation
                // -0.1758153103430687 W/K; halving time cannot fix its sign.
                let mut legacy = operator(&flow, 0.0, 0.0, 0.0, 0.01);
                let error = legacy
                    .advance_with_reservoirs_and_entropy_for(
                        &mut flow,
                        Time::new(1e-4),
                        controls,
                        entropy_controls,
                        baths,
                        || false,
                    )
                    .unwrap_err();
                match error {
                    DissipationError::StepFailure { last, .. } => {
                        assert!(matches!(*last, DissipationError::EntropyDecrease { .. }));
                    }
                    other => panic!("unexpected legacy boundary refusal: {other:?}"),
                }
                assert_eq!(flow.states, original);
            }
            let mut op = operator(&flow, 0.0, 0.0, 0.0, 0.01);
            let (report, entropy) = op
                .advance_entropy_consistent_with_reservoirs_for(
                    &mut flow,
                    Time::new(1e-4),
                    controls,
                    entropy_controls,
                    baths,
                    || false,
                )
                .unwrap();
            assert_eq!(report.accepted_steps, 1);
            assert_eq!(report.rejected_steps, 0);
            assert_eq!(flow.time().value(), 0.0);
            let transfer = if end == 0 {
                report.left_boundary_change
            } else {
                report.right_boundary_change
            };
            let input = bath_input(COUNTER, [0.9, 0.1], 10.0, transfer);
            let change = physical_entropy(&flow, COUNTER) - before;
            assert!(change < 0.0 && input < change);
            assert_eq!(entropy.boundary_entropy_input_j_k[1 - end], 0.0);
            assert!((entropy.boundary_entropy_input_j_k[end] - input).abs() < 2e-13);
            assert!((entropy.source_entropy_change_j_k - change).abs() < 2e-13);
            assert!((entropy.source_entropy_generation_j_k - (change - input)).abs() < 3e-13);
            assert!(entropy.entropy_generation_j_k > 0.0);
            // Independent binary Maxwell--Stefan and polynomial source SSP-RK2
            // calculation: j_A initially -0.01 ln(9), divided caloric energy
            // 0.16857730185383066 W/m² in the cold-to-hot orientation.
            assert!((entropy.entropy_generation_j_k - 6.758876630478756e-6).abs() < 3e-12);
            assert!((input + 1.0570280628339708e-5).abs() < 3e-12);
            assert!(flow.cell(end).unwrap().partial_densities[0] > 0.1);
            assert!(entropy.maximum_relative_species_residual < 1e-12);
            if let Some((states, account)) = &gauge_reference {
                assert_eq!(&flow.states, states);
                assert_eq!(&entropy, account);
            } else {
                assert!(report.absolute_energy_residual_j.abs() < 2e-12);
                gauge_reference = Some((flow.states.clone(), entropy));
            }
        }
    }
}

#[test]
fn g1_variable_calorics_reciprocal_heat_and_viscous_work_use_actual_bath_energy() {
    let mut gauge_reference = None;
    for (formation, heat_offset) in [(0.0, 0.0), (1e18, 1e12)] {
        let mut flow = fixture(VARIABLE, formation, [0.8, 0.2], 1.2, 2.0);
        let bath_rho = [0.3, 0.7];
        let bath = flow
            .gas()
            .conserved(&bath_rho, 0.0, pressure(VARIABLE, bath_rho, 7.0))
            .unwrap();
        let cross = SoretDufourTransport::new(
            law(&flow, 0.07, 0.02, 0.09, 0.03),
            "declared-reservoir-cross-heat",
            1,
            TransportDomain {
                minimum_temperature_k: T0,
                maximum_temperature_k: 50.0,
                minimum_pressure_pa: 1e-8,
                maximum_pressure_pa: 1e8,
            },
            vec![
                (
                    SpeciesId::new("A").unwrap(),
                    MassSpecificEnergyQuantityV1::new(125.0 + heat_offset),
                ),
                (
                    SpeciesId::new("B").unwrap(),
                    MassSpecificEnergyQuantityV1::new(-75.0 + heat_offset),
                ),
            ],
        )
        .unwrap();
        let mut op = DissipationOperator::new(&flow, cross, || false).unwrap();
        let before = physical_entropy(&flow, VARIABLE);
        let original_kinetic = kinetic(&flow);
        let (report, entropy) = op
            .advance_entropy_consistent_with_reservoirs_for(
                &mut flow,
                Time::new(1e-6),
                DissipationControl::default(),
                MolecularEntropyControl::default(),
                [None, Some(bath.as_cell())],
                || false,
            )
            .unwrap();
        assert_eq!(report.accepted_steps, 1);
        assert_eq!(report.rejected_steps, 0);
        let input = bath_input(VARIABLE, bath_rho, 7.0, report.right_boundary_change);
        let change = physical_entropy(&flow, VARIABLE) - before;
        assert!((entropy.source_entropy_change_j_k - change).abs() < 2e-12);
        assert!((entropy.boundary_entropy_input_j_k[1] - input).abs() < 2e-12);
        assert!((entropy.source_entropy_generation_j_k - (change - input)).abs() < 3e-12);
        assert!(report.right_boundary_change.energy > 0.0);
        assert!(report.right_boundary_change.momentum < 0.0);
        assert!(kinetic(&flow) < original_kinetic);
        let endpoint = flow
            .gas()
            .primitive(flow.cell(1).unwrap(), || false)
            .unwrap();
        assert!(endpoint.temperature > 2.0 && endpoint.velocity < 1.2);
        // Independent binary friction + Fourier + axial viscous initial rate:
        // 0.7*(105.00836742340749 + 0.32142857142857145
        //      + 0.10491428571428571) = 73.80429719638524 W/K.
        // Analytic polynomial calorics and actual two-stage bath transfers give:
        assert!((entropy.entropy_generation_j_k - 7.379850861894129e-5).abs() < 3e-11);
        assert!((input - 3.184227219979743e-5).abs() < 3e-11);
        // Interior Fourier+Dufour heat is NOT bath heat for divided enthalpies
        // and symmetric stress work. Substituting it gives about 32.64 W/K
        // instead of the actual 31.84 W/K bath input in this fixture.
        assert!((input / 1e-6 - 32.64062670039358).abs() > 0.7);
        assert!(entropy.maximum_relative_species_residual < 1e-12);
        assert_eq!(flow.time().value(), 0.0);
        if let Some((states, account)) = &gauge_reference {
            assert_eq!(&flow.states, states);
            assert_eq!(&entropy, account);
        } else {
            assert!(report.absolute_energy_residual_j.abs() < 2e-12);
            gauge_reference = Some((flow.states.clone(), entropy));
        }
    }
}

#[test]
fn g4_one_sided_support_and_late_cancel_preserve_state_and_allow_reuse() {
    let make = || fixture(COUNTER, 0.0, [0.5, 0.5], 0.0, 2.0);
    let mut flow = make();
    let original = flow.states.clone();
    let pure = flow.gas().conserved(&[1.0, 0.0], 0.0, 3.0).unwrap();
    let bath = flow.gas().conserved(&[0.8, 0.2], 0.0, 3.0).unwrap();
    let controls = DissipationControl {
        maximum_step_s: 5e-4,
        ..Default::default()
    };
    let mut op = operator(&flow, 0.0, 0.0, 0.02, 0.01);
    let error = op
        .advance_entropy_consistent_with_reservoirs_for(
            &mut flow,
            Time::new(0.002),
            controls,
            MolecularEntropyControl::default(),
            [Some(pure.as_cell()), None],
            || false,
        )
        .unwrap_err();
    assert!(matches!(
        error,
        DissipationError::Invalid(
            "entropy-conjugate diffusion requires common positive species support"
        )
    ));
    assert_eq!(flow.states, original);
    assert_eq!(flow.time().value(), 0.0);

    // Match the reference's initial caloric hints before counting checkpoints.
    let mut op = operator(&flow, 0.0, 0.0, 0.02, 0.01);
    let mut reference = make();
    let mut reference_op = operator(&reference, 0.0, 0.0, 0.02, 0.01);
    let mut polls = 0;
    reference_op
        .advance_entropy_consistent_with_reservoirs_for(
            &mut reference,
            Time::new(0.002),
            controls,
            MolecularEntropyControl::default(),
            [Some(bath.as_cell()), None],
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    assert!(polls > 2);
    let mut visits = 0;
    let result = op.advance_entropy_consistent_with_reservoirs_for(
        &mut flow,
        Time::new(0.002),
        controls,
        MolecularEntropyControl::default(),
        [Some(bath.as_cell()), None],
        || {
            visits += 1;
            visits >= polls
        },
    );
    assert_eq!(result, Err(DissipationError::Cancelled));
    assert_eq!(flow.states, original);
    assert_eq!(flow.time().value(), 0.0);
    let (_, entropy) = op
        .advance_entropy_consistent_with_reservoirs_for(
            &mut flow,
            Time::new(0.002),
            controls,
            MolecularEntropyControl::default(),
            [Some(bath.as_cell()), None],
            || false,
        )
        .unwrap();
    assert!(entropy.entropy_generation_j_k > 0.0);
    for (actual, expected) in flow.states.iter().zip(&reference.states) {
        assert!((actual - expected).abs() < 1e-11);
    }
}
