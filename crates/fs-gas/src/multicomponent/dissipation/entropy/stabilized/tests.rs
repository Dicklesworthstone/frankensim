//! Physical endpoint entropy for the original RKL2 molecular recurrence.
//! Internal Legendre stages are not independent physical entropy substeps.

use crate::multicomponent::dissipation::{
    ConstantTransport, DissipationControl, DissipationError, DissipationOperator,
    MolecularEntropyControl, StabilizedDissipationControl, TransportCoefficients, TransportDomain,
    TransportLaw,
};
use crate::multicomponent::tests::{gas, source};
use crate::multicomponent::{Boundaries, Boundary, MixtureEuler1d, MixtureGas};
use fs_qty::{DynViscosity, Length, Pressure, Temperature, Time};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;
use fs_thermochem::equilibrium::composition::IdealGasEquilibriumV1;
use fs_thermochem::transport::kinetic::ThermalConductivityQuantityV1;

fn controls(step: f64) -> StabilizedDissipationControl {
    StabilizedDissipationControl {
        transport: DissipationControl {
            maximum_step_s: step,
            ..DissipationControl::default()
        },
        ..StabilizedDissipationControl::default()
    }
}

// All manufactured sources below have M=R, so p=T*sum(rho_i).
fn field(
    gas: MixtureGas,
    count: usize,
    boundaries: Boundaries,
    initial: impl Fn(usize) -> (Vec<f64>, f64),
) -> MixtureEuler1d {
    let states = (0..count)
        .map(|i| {
            let (rho, t) = initial(i);
            let pressure = t * rho.iter().sum::<f64>();
            gas.conserved(&rho, 0.0, pressure).unwrap()
        })
        .collect();
    MixtureEuler1d::new(
        gas,
        Length::new(1.0 / count as f64),
        0.7,
        states,
        boundaries,
        || false,
    )
    .unwrap()
}

fn law(flow: &MixtureEuler1d, conductivity: f64, diffusion: f64) -> ConstantTransport {
    let n = flow.gas().species_count();
    let mut binary = vec![diffusion; n * n];
    for i in 0..n {
        binary[i * n + i] = 0.0;
    }
    ConstantTransport::new(
        flow.gas().source().clone(),
        "rkl2-entropy-fixture",
        1,
        TransportDomain {
            minimum_temperature_k: 0.001,
            maximum_temperature_k: 50.0,
            minimum_pressure_pa: 1e-8,
            maximum_pressure_pa: 1e8,
        },
        TransportCoefficients {
            dynamic_viscosity: DynViscosity::new(0.0),
            bulk_viscosity: DynViscosity::new(0.0),
            thermal_conductivity: ThermalConductivityQuantityV1::new(conductivity),
        },
        binary,
    )
    .unwrap()
}

fn temperature(flow: &MixtureEuler1d, i: usize) -> f64 {
    flow.gas()
        .primitive(flow.cell(i).unwrap(), || false)
        .unwrap()
        .temperature
}

// Original ideal-mixture state entropy up to constant species reference terms.
// Those constants cancel because these fixtures conserve every species.
fn entropy(flow: &MixtureEuler1d, cv: &[f64]) -> f64 {
    (0..flow.cell_count())
        .map(|i| {
            let q = flow.cell(i).unwrap();
            let t = temperature(flow, i);
            flow.mesh().volumes()[i]
                * q.partial_densities
                    .iter()
                    .zip(cv)
                    .filter(|(rho, _)| **rho > 0.0)
                    .map(|(&rho, &cv)| rho * (cv * t.ln() - rho.ln()))
                    .sum::<f64>()
        })
        .sum()
}

fn heat(count: usize) -> MixtureEuler1d {
    field(
        gas(vec![source("A", R, 3.5, 0.0, 0.0)]),
        count,
        Boundaries::Periodic,
        |i| {
            let x = (i as f64 + 0.5) / count as f64;
            let checker = if i % 2 == 0 { 1.0 } else { -1.0 };
            (
                vec![1.0],
                2.0 + 0.2 * (core::f64::consts::TAU * x).sin() + 0.01 * checker,
            )
        },
    )
}

// Evaluate the Legendre polynomial independently of the state/flux recurrence.
fn amplification(stages: usize, z: f64) -> f64 {
    let s = stages as f64;
    let x = 1.0 + 4.0 * z / (s * s + s - 2.0);
    let (mut previous, mut value) = (1.0, x);
    for j in 2..=stages {
        let j = j as f64;
        let next = ((2.0 * j - 1.0) * x * value - (j - 1.0) * previous) / j;
        previous = value;
        value = next;
    }
    let b = (s * s + s - 2.0) / (2.0 * s * (s + 1.0));
    1.0 - b + b * value
}

#[test]
fn g1_multistage_heat_matches_the_spectrum_and_original_state_entropy_with_less_work() {
    let n = 128;
    let dt = 0.05;
    let mut flow = heat(n);
    let initial: Vec<_> = (0..n).map(|i| temperature(&flow, i)).collect();
    let mut op = DissipationOperator::new(&flow, law(&flow, 0.05, 0.0), || false).unwrap();
    let (report, account) = op
        .advance_stabilized_entropy_consistent_for(
            &mut flow,
            Time::new(dt),
            controls(dt),
            MolecularEntropyControl::default(),
            || false,
        )
        .unwrap();
    assert_eq!(report.transport.accepted_steps, 1);
    assert_eq!(report.transport.rejected_steps, 0);
    assert!(report.maximum_stages_used > 2);
    assert_eq!(report.operator_evaluations, report.maximum_stages_used + 1);
    assert!(report.maximum_scaled_parabolic_number <= 0.40000000000001);

    let scale = 4.0 * 0.02 * (n * n) as f64;
    let low = amplification(
        report.maximum_stages_used,
        -dt * scale * (core::f64::consts::PI / n as f64).sin().powi(2),
    );
    let high = amplification(report.maximum_stages_used, -dt * scale);
    let mut expected_entropy = 0.0;
    for (i, &old) in initial.iter().enumerate() {
        let x = (i as f64 + 0.5) / n as f64;
        let checker = if i % 2 == 0 { 1.0 } else { -1.0 };
        let expected = 2.0 + 0.2 * low * (core::f64::consts::TAU * x).sin() + 0.01 * high * checker;
        assert!((temperature(&flow, i) - expected).abs() < 3e-11);
        expected_entropy += (0.7 * 2.5 / n as f64) * (expected / old).ln();
    }
    assert!(expected_entropy > 0.0);
    assert!((account.entropy_change_j_k - expected_entropy).abs() < 5e-12);
    assert!((account.source_entropy_change_j_k - expected_entropy).abs() < 5e-12);
    assert!((account.minimum_step_entropy_change_j_k - expected_entropy).abs() < 5e-12);
    assert!(report.transport.absolute_energy_residual_j.abs() < 3e-12);
    assert_eq!(flow.time().value(), 0.0);

    let mut explicit = heat(n);
    let mut op = DissipationOperator::new(&explicit, law(&explicit, 0.05, 0.0), || false).unwrap();
    let (ordinary, _) = op
        .advance_entropy_consistent_for(
            &mut explicit,
            Time::new(dt),
            DissipationControl::default(),
            MolecularEntropyControl::default(),
            || false,
        )
        .unwrap();
    // Compare actual constitutive RHS work, not a wall-clock performance claim.
    assert!(2 * report.operator_evaluations < 2 * ordinary.accepted_steps);
}

fn mixture(formation: f64) -> MixtureEuler1d {
    field(
        gas(vec![
            source("A", R, 3.5, 0.0, formation),
            source("B", R, 6.5, 0.0, -2.0 * formation),
            source("C", R, 4.0, 0.0, 3.0 * formation),
        ]),
        2,
        Boundaries::Ends {
            left: Boundary::Transmissive,
            right: Boundary::Transmissive,
        },
        |i| {
            if i == 0 {
                (vec![0.1, 0.9, 0.0], 1.0)
            } else {
                (vec![0.9, 0.1, 0.0], 10.0)
            }
        },
    )
}

#[test]
fn g1_unequal_calorics_repair_the_negative_face_without_inventing_energy_or_species() {
    let mut legacy = mixture(0.0);
    let before = entropy(&legacy, &[2.5, 5.5, 3.0]);
    let mut op = DissipationOperator::new(&legacy, law(&legacy, 0.0, 0.01), || false).unwrap();
    op.advance_stabilized_for(&mut legacy, Time::new(0.1), controls(0.1), || false)
        .unwrap();
    assert!(entropy(&legacy, &[2.5, 5.5, 3.0]) - before < -0.003);

    let mut reference = None;
    for formation in [0.0, 1e18] {
        let mut flow = mixture(formation);
        let before = entropy(&flow, &[2.5, 5.5, 3.0]);
        let mut op = DissipationOperator::new(&flow, law(&flow, 0.0, 0.01), || false).unwrap();
        let (report, account) = op
            .advance_stabilized_entropy_consistent_for(
                &mut flow,
                Time::new(0.1),
                controls(0.1),
                MolecularEntropyControl {
                    absolute_tolerance_j_k: 1e-14,
                    ..Default::default()
                },
                || false,
            )
            .unwrap();
        assert_eq!(report.maximum_stages_used, 2);
        assert_eq!(report.transport.accepted_steps, 1);
        assert_eq!(report.transport.rejected_steps, 0);
        let measured = entropy(&flow, &[2.5, 5.5, 3.0]) - before;
        // Independent equal-M binary equations: j_A=-D*x*(1-x)*Delta(logit(x))/dx;
        // h_i*=Cp_i*TL*TR*ln(TR/TL)/(TR-TL)-Cv_i*Tmin. Apply the two-stage
        // Butcher update Y1=Y0+dt*f0/3; Y2=Y0+dt*(-f0/2+3*f(Y1)/2).
        assert!((measured - 0.0066120946863565955).abs() < 3e-11);
        assert!((temperature(&flow, 0) - 0.996095063603396).abs() < 3e-11);
        assert!((temperature(&flow, 1) - 9.965487298703065).abs() < 3e-11);
        assert!((account.entropy_change_j_k - measured).abs() < 3e-12);
        assert_eq!(report.transport.boundary_change.energy, 0.0);
        assert_eq!(report.transport.after.species_masses[2], 0.0);
        for i in 0..2 {
            assert_eq!(flow.cell(i).unwrap().partial_densities[2], 0.0);
        }
        if let Some((states, entropy)) = &reference {
            assert_eq!(&flow.states, states);
            assert_eq!(&account, entropy);
        } else {
            assert!(report.transport.absolute_energy_residual_j.abs() < 3e-12);
            reference = Some((flow.states.clone(), account));
        }
    }
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
        mole_fractions: &[f64],
        binary_diffusivities: &mut [f64],
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<TransportCoefficients, DissipationError> {
        let mut result = self.0.evaluate(
            temperature,
            pressure,
            mole_fractions,
            binary_diffusivities,
            cancelled,
        )?;
        // Positive smooth manufactured material law. Unequal capacities change
        // the face temperature during transfer, exposing RKL2's signed weights.
        result.thermal_conductivity = ThermalConductivityQuantityV1::new(
            0.001 + 0.999 / (1.0 + (-80.0 * (temperature.value() - 5.4)).exp()),
        );
        Ok(result)
    }
}

fn switching_heat() -> MixtureEuler1d {
    field(
        gas(vec![source("A", R, 3.5, 0.0, 0.0)]),
        2,
        Boundaries::Ends {
            left: Boundary::Transmissive,
            right: Boundary::Transmissive,
        },
        |i| {
            if i == 0 {
                (vec![1.0], 10.0)
            } else {
                (vec![2.0], 1.0)
            }
        },
    )
}

#[test]
fn g4_positive_caloric_and_parabolic_stages_still_need_endpoint_entropy_retries() {
    let mut legacy = switching_heat();
    let before = entropy(&legacy, &[2.5]);
    let mut op = DissipationOperator::new(
        &legacy,
        SwitchingConductivity(law(&legacy, 1.0, 0.0)),
        || false,
    )
    .unwrap();
    let ordinary = op
        .advance_stabilized_for(&mut legacy, Time::new(0.2), controls(0.2), || false)
        .unwrap();
    assert_eq!(ordinary.transport.rejected_steps, 0);
    assert!(ordinary.maximum_scaled_parabolic_number < 0.321);
    assert!(temperature(&legacy, 0) > 0.001 && temperature(&legacy, 1) > 0.001);
    assert!((entropy(&legacy, &[2.5]) - before + 2.097337615085594).abs() < 2e-8);

    let mut flow = switching_heat();
    let original = flow.states.clone();
    let mut op =
        DissipationOperator::new(&flow, SwitchingConductivity(law(&flow, 1.0, 0.0)), || false)
            .unwrap();
    let mut strict = controls(0.2);
    strict.transport.maximum_halvings = 0;
    let error = op
        .advance_stabilized_entropy_consistent_for(
            &mut flow,
            Time::new(0.2),
            strict,
            MolecularEntropyControl::default(),
            || false,
        )
        .unwrap_err();
    match error {
        DissipationError::StepFailure { last, .. } => {
            assert!(matches!(*last, DissipationError::EntropyDecrease { .. }));
        }
        other => panic!("expected actual endpoint entropy refusal, got {other:?}"),
    }
    assert_eq!(flow.states, original);
    assert_eq!(flow.time().value(), 0.0);

    let (report, account) = op
        .advance_stabilized_entropy_consistent_for(
            &mut flow,
            Time::new(0.2),
            controls(0.2),
            MolecularEntropyControl::default(),
            || false,
        )
        .unwrap();
    assert_eq!(report.transport.accepted_steps, 2);
    assert_eq!(report.transport.rejected_steps, 2);
    assert_eq!(account.rejected_entropy_steps, 2);
    assert!((temperature(&flow, 0) - 9.34773553841244).abs() < 2e-8);
    assert!((temperature(&flow, 1) - 1.3261322307937802).abs() < 2e-8);
    let measured = entropy(&flow, &[2.5]) - before;
    assert!((measured - 0.43494696812052647).abs() < 2e-8);
    assert!((account.entropy_change_j_k - measured).abs() < 3e-11);
    assert!(account.minimum_step_entropy_change_j_k > 0.0);
    assert!(report.transport.absolute_energy_residual_j.abs() < 3e-12);
}

#[test]
fn g4_late_cancel_and_rhs_exhaustion_leave_no_partial_entropy_or_field_commit() {
    let control = controls(0.02);
    let mut reference = heat(64);
    let mut op =
        DissipationOperator::new(&reference, law(&reference, 0.05, 0.0), || false).unwrap();
    let mut polls = 0;
    let (_, expected) = op
        .advance_stabilized_entropy_consistent_for(
            &mut reference,
            Time::new(0.08),
            control,
            MolecularEntropyControl::default(),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    assert!(polls > 2);

    let mut flow = heat(64);
    let original = flow.states.clone();
    let mut op = DissipationOperator::new(&flow, law(&flow, 0.05, 0.0), || false).unwrap();
    let mut count = 0;
    let error = op
        .advance_stabilized_entropy_consistent_for(
            &mut flow,
            Time::new(0.08),
            control,
            MolecularEntropyControl::default(),
            || {
                count += 1;
                count >= polls - 2
            },
        )
        .unwrap_err();
    assert_eq!(error, DissipationError::Cancelled);
    assert_eq!(flow.states, original);
    assert_eq!(flow.time().value(), 0.0);

    let error = op
        .advance_stabilized_entropy_consistent_for(
            &mut flow,
            Time::new(0.08),
            StabilizedDissipationControl {
                maximum_operator_evaluations: 1,
                ..control
            },
            MolecularEntropyControl::default(),
            || false,
        )
        .unwrap_err();
    assert!(matches!(error, DissipationError::OperatorWorkLimit { .. }));
    assert_eq!(flow.states, original);

    let (_, account) = op
        .advance_stabilized_entropy_consistent_for(
            &mut flow,
            Time::new(0.08),
            control,
            MolecularEntropyControl::default(),
            || false,
        )
        .unwrap();
    assert!(
        flow.states
            .iter()
            .zip(&reference.states)
            .all(|(a, b)| (a - b).abs() < 2e-12)
    );
    assert!((account.entropy_change_j_k - expected.entropy_change_j_k).abs() < 2e-12);
    let accepted = flow.states.clone();
    let (noop, account) = op
        .advance_stabilized_entropy_consistent_for(
            &mut flow,
            Time::new(0.0),
            control,
            MolecularEntropyControl::default(),
            || false,
        )
        .unwrap();
    assert_eq!(noop.operator_evaluations, 0);
    assert_eq!(account.entropy_change_j_k, 0.0);
    assert_eq!(flow.states, accepted);
}
