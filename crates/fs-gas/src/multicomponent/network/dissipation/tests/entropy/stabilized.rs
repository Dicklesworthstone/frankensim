//! Physical RKL endpoints, finite thermal storage and complete entropy admission.
use super::*;
use crate::multicomponent::dissipation::StabilizedDissipationControl;
use fs_qty::{Pressure, Temperature};
use fs_thermochem::equilibrium::composition::IdealGasEquilibriumV1;

fn controls(step: f64) -> StabilizedDissipationControl {
    StabilizedDissipationControl {
        transport: DissipationControl {
            maximum_step_s: step,
            ..Default::default()
        },
        ..Default::default()
    }
}

// Legendre polynomial, independent of the native conservative stage recurrence.
fn amplification(stages: usize, z: f64) -> f64 {
    let s = stages as f64;
    let x = 1.0 + 4.0 * z / (s * s + s - 2.0);
    let (mut old, mut previous) = (1.0, x);
    for j in 2..=stages {
        let j = j as f64;
        let next = ((2.0 * j - 1.0) * x * previous - (j - 1.0) * old) / j;
        old = previous;
        previous = next;
    }
    let b = (s * s + s - 2.0) / (2.0 * s * (s + 1.0));
    1.0 - b + b * previous
}

#[test]
fn g1_stiff_finite_heat_capacity_matches_legendre_mode_and_source_entropy() {
    let g = gas(1, false, 0.0);
    let volume = 0.0002;
    let duration = 0.2;
    let make = || {
        let network = two_pools(
            &g,
            g.conserved(&[1.0], 0.0, 2.1).unwrap(),
            g.conserved(&[1.0], 0.0, 1.9).unwrap(),
            volume,
        );
        let op = NetworkDissipationOperator::new(
            &network,
            law(&g, 0.0, 0.0, 0.03, 0.0),
            vec![Length::new(PLENUM_DISTANCE); 2],
            || false,
        )
        .unwrap();
        (network, op)
    };
    let (mut network, mut op) = make();
    let before = source_entropy(&network, &[(1.0, 2.5, 0.0)]);
    let (report, entropy) = op
        .advance_stabilized_entropy_consistent_for(
            &mut network,
            Time::new(duration),
            controls(duration),
            MolecularEntropyControl::default(),
            || false,
        )
        .unwrap();
    assert_eq!(report.transport.accepted_steps, 1);
    assert_eq!(report.transport.rejected_steps, 0);
    assert!(report.maximum_stages_used > 4);
    assert_eq!(report.operator_evaluations, report.maximum_stages_used + 1);
    assert!(report.maximum_scaled_parabolic_number <= 0.40000000000001);
    assert!(report.transport.maximum_parabolic_number > 10.0);
    let cb = (2.0 * AREA * WIDTH) * 2.5;
    let cp = volume * 2.5;
    let conductance = 2.0 * AREA * 0.03 / (0.5 * WIDTH + PLENUM_DISTANCE);
    let rate = conductance * (1.0 / cb + 1.0 / cp);
    let mean = (cb * 2.1 + cp * 1.9) / (cb + cp);
    let decay = amplification(report.maximum_stages_used, -rate * duration);
    let expected = [mean + (2.1 - mean) * decay, mean + (1.9 - mean) * decay];
    for cell in 0..2 {
        let q = network.branches()[0].cell(cell).unwrap();
        close(
            g.primitive(q, || false).unwrap().temperature,
            expected[0],
            3e-11,
        );
        assert_eq!(q.partial_densities, [1.0]);
    }
    close(
        g.primitive(network.plenums()[0].state(), || false)
            .unwrap()
            .temperature,
        expected[1],
        3e-11,
    );
    let expected_entropy = cb * (expected[0] / 2.1).ln() + cp * (expected[1] / 1.9).ln();
    assert!(expected_entropy > 0.0);
    close(entropy.entropy_change_j_k, expected_entropy, 3e-12);
    close(
        entropy.source_entropy_change_j_k,
        source_entropy(&network, &[(1.0, 2.5, 0.0)]) - before,
        3e-12,
    );
    close(
        entropy.minimum_step_entropy_change_j_k,
        expected_entropy,
        3e-12,
    );
    close(report.transport.absolute_energy_residual_j, 0.0, 3e-12);
    for port in op.accepted_port_exchanges().unwrap() {
        close(port.energy_to_plenum, 0.5 * cp * (expected[1] - 1.9), 3e-13);
        assert_eq!(port.mass_to_plenum(), 0.0);
    }
    assert_eq!(network.time().value(), 0.0);
    assert_eq!(op.last_report(), Some(report.transport));

    let (mut explicit, mut explicit_op) = make();
    let (explicit_report, _) = explicit_op
        .advance_entropy_consistent_for(
            &mut explicit,
            Time::new(duration),
            controls(duration).transport,
            MolecularEntropyControl::default(),
            || false,
        )
        .unwrap();
    // Actual spatial work at the same interval; this is a stability-work
    // comparison, not a claim of equal temporal error or wall-clock speedup.
    assert!(2 * report.operator_evaluations < 2 * explicit_report.accepted_steps);
}

#[test]
fn g1_conjugate_rkl_ports_keep_real_calorics_common_absence_and_reference_invariance() {
    let calorics = [(1.0, 2.5, 0.0), (1.0, 5.5, 0.0), (1.0, 3.0, 0.0)];
    let mut reference = None;
    for formation in [0.0, 1e18] {
        let mut network = contrasting_pools(formation, true, false);
        let initial = source_entropy(&network, &calorics);
        let mut op = diffusion(&network);
        let (report, entropy) = op
            .advance_stabilized_entropy_consistent_for(
                &mut network,
                Time::new(0.01),
                controls(0.01),
                MolecularEntropyControl::default(),
                || false,
            )
            .unwrap();
        let actual = source_entropy(&network, &calorics) - initial;
        assert!(actual > 0.0);
        close(entropy.source_entropy_change_j_k, actual, 3e-12);
        close(entropy.entropy_change_j_k, actual, 3e-12);
        assert_eq!(report.transport.after.species_masses[2], 0.0);
        assert_eq!(report.transport.boundary_change.energy, 0.0);
        assert!(entropy.maximum_relative_species_residual < 1e-13);
        if let Some((state, measured)) = &reference {
            assert_eq!(&fields(&network), state);
            assert_eq!(&entropy, measured);
        } else {
            close(report.transport.absolute_energy_residual_j, 0.0, 3e-12);
            reference = Some((fields(&network), entropy));
        }
    }
    let mut unsupported = contrasting_pools(0.0, false, true);
    let before = fields(&unsupported);
    let error = diffusion(&unsupported)
        .advance_stabilized_entropy_consistent_for(
            &mut unsupported,
            Time::new(0.01),
            controls(0.01),
            MolecularEntropyControl::default(),
            || false,
        )
        .unwrap_err();
    assert!(matches!(error, DissipationError::Invalid(_)));
    assert_eq!(fields(&unsupported), before);
}

struct SwitchingConductivity(ConstantTransport);
impl TransportLaw for SwitchingConductivity {
    fn source(&self) -> &IdealGasEquilibriumV1 {
        self.0.source()
    }
    fn evaluate(
        &mut self,
        temperature: Temperature,
        pressure: Pressure,
        fractions: &[f64],
        binary: &mut [f64],
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<TransportCoefficients, DissipationError> {
        let mut value = self
            .0
            .evaluate(temperature, pressure, fractions, binary, cancelled)?;
        // Positive smooth k(T); unequal capacities move the face temperature
        // and expose the negative RKL anchor weight without a constitutive fault.
        value.thermal_conductivity = ThermalConductivityQuantityV1::new(
            0.001 + 0.999 / (1.0 + (-80.0 * (temperature.value() - 5.4)).exp()),
        );
        Ok(value)
    }
}

fn switching_heat() -> (
    MixtureNetwork1d,
    NetworkDissipationOperator<SwitchingConductivity>,
) {
    let g = gas(1, false, 0.0);
    let branch = MixtureEuler1d::new(
        g.clone(),
        Length::new(0.5),
        0.7,
        vec![
            g.conserved(&[1.0], 0.0, 10.0).unwrap(),
            g.conserved(&[2.0], 0.0, 2.0).unwrap(),
        ],
        Boundaries::Ends {
            left: Boundary::Reflecting,
            right: Boundary::Reflecting,
        },
        || false,
    )
    .unwrap();
    let network = MixtureNetwork1d::new(vec![branch], vec![], vec![], || false).unwrap();
    let op = NetworkDissipationOperator::new(
        &network,
        SwitchingConductivity(law(&g, 0.0, 0.0, 1.0, 0.0)),
        vec![],
        || false,
    )
    .unwrap();
    (network, op)
}

#[test]
fn g4_admissible_legendre_stages_still_refuse_negative_physical_entropy_endpoints() {
    let (mut legacy, mut legacy_op) = switching_heat();
    let original_entropy = source_entropy(&legacy, &[(1.0, 2.5, 0.0)]);
    let ordinary = legacy_op
        .advance_stabilized_for(&mut legacy, Time::new(0.2), controls(0.2), || false)
        .unwrap();
    assert_eq!(ordinary.transport.rejected_steps, 0);
    assert!(source_entropy(&legacy, &[(1.0, 2.5, 0.0)]) < original_entropy - 1.0);

    let (mut network, mut op) = switching_heat();
    let before = fields(&network);
    let mut strict = controls(0.2);
    strict.transport.maximum_halvings = 0;
    let error = op
        .advance_stabilized_entropy_consistent_for(
            &mut network,
            Time::new(0.2),
            strict,
            MolecularEntropyControl::default(),
            || false,
        )
        .unwrap_err();
    assert!(matches!(error, DissipationError::StepFailure { last, .. }
        if matches!(*last, DissipationError::EntropyDecrease { .. })));
    assert_eq!(fields(&network), before);
    assert!(op.last_report().is_none());
    let (report, entropy) = op
        .advance_stabilized_entropy_consistent_for(
            &mut network,
            Time::new(0.2),
            controls(0.2),
            MolecularEntropyControl::default(),
            || false,
        )
        .unwrap();
    assert!(report.transport.rejected_steps > 0);
    assert_eq!(
        entropy.rejected_entropy_steps,
        report.transport.rejected_steps
    );
    assert!(entropy.minimum_step_entropy_change_j_k > 0.0);
    close(
        entropy.source_entropy_change_j_k,
        source_entropy(&network, &[(1.0, 2.5, 0.0)]) - original_entropy,
        3e-11,
    );
    assert!(entropy.entropy_change_j_k > 0.0);
    close(report.transport.absolute_energy_residual_j, 0.0, 3e-12);
    assert_eq!(network.time().value(), 0.0);
}
