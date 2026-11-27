//! G0/G4: held observations retain the original physical law and solver tiers.

use super::*;
use fs_exec::solver::snapshot_v2::StateIdentityCharterV2;
use fs_ports::component::{EndOfStep, ObservationError, StageResidual, StateSignals};

fn signals() -> Vec<PortSpec> {
    vec![PortSpec::signal_out(
        StableId::new("state").unwrap(),
        NonZeroUsize::new(2).unwrap(),
    )]
}

fn observe(
    _: &SpringMassDamper,
    state: &Oscillator,
    out: &mut [f64],
    _: &Cx<'_>,
) -> Result<(), Cancelled> {
    out.copy_from_slice(&[state.position_m, state.velocity_m_s]);
    Ok(())
}

#[test]
fn appended_signals_observe_old_state_and_preserve_original_physics_and_tangent() {
    let gate = CancelGate::new_clock_free();
    for component in [SpringMassDamper::new(Fault::None), undamped()] {
        let observed = StateSignals::new(component, signals(), observe).unwrap();
        assert_eq!(observed.class(), observed.inner().class());
        let layout = PortLayout::new(observed.ports()).unwrap();
        assert_eq!((layout.input_len(), layout.output_len()), (2, 4));
        with_cx(&gate, |cx| {
            check_component(&observed, &case(), cx).expect("wrapped physical law conforms");
            check_port_tangent(&observed, &case(), 1e-7, 1e-12, cx).expect("wrapped tangent");
            let (mut plain, mut added) = (STALE, STALE);
            let (mut plain_out, mut added_out) = ([0.0; 2], [42.0; 4]);
            let native = observed
                .inner()
                .propose(&START, case().dt_s, &INPUTS, &mut plain, &mut plain_out, cx)
                .unwrap();
            let wrapped = observed
                .propose(&START, case().dt_s, &INPUTS, &mut added, &mut added_out, cx)
                .unwrap();
            assert_eq!(plain, added);
            assert_eq!(native, wrapped);
            assert_eq!(&added_out[..2], &plain_out);
            assert_eq!(&added_out[2..], &[START.position_m, START.velocity_m_s]);
            assert_ne!(&added_out[2..], &[added.position_m, added.velocity_m_s]);
            let mut tangent = [99.0; 4];
            observed
                .jvp(&START, case().dt_s, &INPUTS, &[1.0, -2.0], &mut tangent, cx)
                .unwrap();
            assert_eq!(&tangent[2..], &[0.0; 2]);
            let mut pullback = [99.0; 2];
            observed
                .vjp(
                    &START,
                    case().dt_s,
                    &INPUTS,
                    &[0.0, 0.0, 1e100, -1e100],
                    &mut pullback,
                    cx,
                )
                .unwrap();
            assert_eq!(pullback, [0.0; 2]);
        });
    }
}

#[test]
fn observations_reject_duplicate_ports_wrong_roles_and_unwritten_or_nonfinite_values() {
    for port in [
        PortSpec::signal_out(StableId::new("drive").unwrap(), NonZeroUsize::MIN),
        PortSpec::signal_in(StableId::new("measurement").unwrap(), NonZeroUsize::MIN),
    ] {
        assert!(
            StateSignals::new(SpringMassDamper::new(Fault::None), vec![port], observe).is_err()
        );
    }
    let gate = CancelGate::new_clock_free();
    for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let observed = StateSignals::new(
            SpringMassDamper::new(Fault::None),
            signals(),
            move |_: &SpringMassDamper, _: &Oscillator, out: &mut [f64], _: &Cx<'_>| {
                out[0] = value;
                out[1] = 0.0;
                Ok(())
            },
        )
        .unwrap();
        with_cx(&gate, |cx| {
            assert!(matches!(
                observed.propose(
                    &START,
                    case().dt_s,
                    &INPUTS,
                    &mut STALE.clone(),
                    &mut [3.0; 4],
                    cx
                ),
                Err(ObservationError::InvalidSignal { index: 0, .. })
            ));
        });
    }
    let partial = StateSignals::new(
        SpringMassDamper::new(Fault::None),
        signals(),
        |_: &SpringMassDamper, _: &Oscillator, out: &mut [f64], _: &Cx<'_>| {
            out[0] = 1.0;
            Ok(())
        },
    )
    .unwrap();
    with_cx(&gate, |cx| {
        assert!(matches!(
            partial.propose(
                &START,
                case().dt_s,
                &INPUTS,
                &mut STALE.clone(),
                &mut [3.0; 4],
                cx
            ),
            Err(ObservationError::InvalidSignal { index: 1, .. })
        ));
        let complete =
            StateSignals::new(SpringMassDamper::new(Fault::None), signals(), observe).unwrap();
        for length in [0, 3, 5] {
            assert!(matches!(
                complete.propose(
                    &START,
                    case().dt_s,
                    &INPUTS,
                    &mut STALE.clone(),
                    &mut vec![0.0; length],
                    cx
                ),
                Err(ObservationError::Contract(
                    ContractViolation::VectorLength { .. }
                ))
            ));
        }
    });
}

#[test]
fn callback_cancellation_is_not_published_as_a_valid_observation() {
    let gate = CancelGate::new_clock_free();
    let observed = StateSignals::new(
        SpringMassDamper::new(Fault::None),
        signals(),
        |_: &SpringMassDamper, _: &Oscillator, out: &mut [f64], _: &Cx<'_>| {
            out.fill(0.0);
            gate.request();
            Ok(())
        },
    )
    .unwrap();
    with_cx(&gate, |cx| {
        assert!(matches!(
            observed.propose(
                &START,
                case().dt_s,
                &INPUTS,
                &mut STALE.clone(),
                &mut [0.0; 4],
                cx
            ),
            Err(ObservationError::Cancelled)
        ));
    });
}

/// Existing oscillator physics with metadata and a residual tier to check delegation.
struct Described(SpringMassDamper);
impl Component for Described {
    type State = Oscillator;
    type Error = &'static str;
    const STATE_CHARTER: Option<StateIdentityCharterV2> = Some(StateIdentityCharterV2 {
        owner: "fs-ports::tests::state-signals",
        state_family: "oscillator",
        schema_grammar: "position and velocity",
        codec_grammar: "two finite f64 bit patterns",
        codec_version: 1,
    });
    fn ports(&self) -> &[PortSpec] {
        self.0.ports()
    }
    fn class(&self) -> Class {
        self.0.class()
    }
    fn stored_energy_j(&self, s: &Oscillator) -> Result<f64, Self::Error> {
        self.0.stored_energy_j(s).map_err(|_| "energy")
    }
    fn stored_mass_kg(&self, _: &Oscillator) -> Result<Option<f64>, Self::Error> {
        Ok(Some(self.0.mass_kg))
    }
    fn stored_entropy_j_per_k(&self, _: &Oscillator) -> Result<Option<f64>, Self::Error> {
        Ok(Some(17.0))
    }
    fn stored_carried_quantity(
        &self,
        _: &Oscillator,
        carried: usize,
    ) -> Result<Option<f64>, Self::Error> {
        Ok((carried == 0).then_some(3.0))
    }
    fn limit_dt_s(&self, _: &Oscillator, dt_s: f64) -> Result<f64, Self::Error> {
        Ok(dt_s.min(0.2))
    }
    fn next_event_s(&self, _: &Oscillator) -> Option<f64> {
        Some(0.15)
    }
    fn end_of_step(&self, _: usize) -> EndOfStep {
        EndOfStep::Untold
    }
    fn propose(
        &self,
        old: &Oscillator,
        dt_s: f64,
        inputs: &[f64],
        new: &mut Oscillator,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, Self::Error> {
        if inputs[0] < 0.0 {
            return Err("negative force refused by native owner");
        }
        if dt_s > 0.2 {
            return Ok(Proposal::Shorten { dt_s: 0.2 });
        }
        self.0
            .propose(old, dt_s, inputs, new, outputs, cx)
            .map_err(|_| "cancelled")
    }
    fn admit(&self, _: &Oscillator, _: &Oscillator, _: &Receipt) -> Result<(), Self::Error> {
        Err("native admit")
    }
    fn encode_state(&self, state: &Oscillator, words: &mut Vec<u64>) -> bool {
        words.extend([state.position_m.to_bits(), state.velocity_m_s.to_bits()]);
        true
    }
    fn decode_state(&self, words: &[u64]) -> Option<Oscillator> {
        if words.len() != 2 {
            return None;
        }
        let (position_m, velocity_m_s) = (f64::from_bits(words[0]), f64::from_bits(words[1]));
        (position_m.is_finite() && velocity_m_s.is_finite()).then_some(Oscillator {
            position_m,
            velocity_m_s,
        })
    }
    fn identity_words(&self, words: &mut Vec<u64>) {
        words.push(self.0.mass_kg.to_bits());
    }
}

fn described_observation(
    _: &Described,
    state: &Oscillator,
    out: &mut [f64],
    _: &Cx<'_>,
) -> Result<(), &'static str> {
    out.copy_from_slice(&[state.position_m, state.velocity_m_s]);
    Ok(())
}

#[test]
fn identities_codecs_inventories_events_and_native_refusals_are_retained() {
    let plain = StateSignals::new(
        Described(SpringMassDamper::new(Fault::None)),
        signals(),
        described_observation,
    )
    .unwrap();
    let mut words = vec![123];
    assert!(!plain.encode_state(&START, &mut words));
    assert_eq!(words, [123]);
    assert!(plain.decode_state(&[0, 0]).is_none());
    let wrapped = plain.with_observation_identity([7; 32]);
    assert!(wrapped.encode_state(&START, &mut words));
    assert_eq!(wrapped.decode_state(&words[1..]), Some(START));
    assert!(wrapped.decode_state(&words).is_none());
    let mut original = Vec::new();
    wrapped.inner().encode_state(&START, &mut original);
    assert_eq!(&words[1..], &original);
    assert_eq!(
        wrapped.stored_energy_j(&START).unwrap(),
        wrapped.inner().stored_energy_j(&START).unwrap()
    );
    assert_eq!(wrapped.stored_mass_kg(&START).unwrap(), Some(0.25));
    assert_eq!(wrapped.stored_entropy_j_per_k(&START).unwrap(), Some(17.0));
    assert_eq!(
        wrapped.stored_carried_quantity(&START, 0).unwrap(),
        Some(3.0)
    );
    assert_eq!(wrapped.stored_carried_quantity(&START, 1).unwrap(), None);
    assert_eq!(wrapped.limit_dt_s(&START, 1.0).unwrap(), 0.2);
    assert_eq!(wrapped.next_event_s(&START), Some(0.15));
    assert_eq!(wrapped.end_of_step(0), EndOfStep::Untold);
    let mut identity_a = Vec::new();
    wrapped.identity_words(&mut identity_a);
    let other = StateSignals::new(
        Described(SpringMassDamper::new(Fault::None)),
        signals(),
        described_observation,
    )
    .unwrap()
    .with_observation_identity([8; 32]);
    let mut identity_b = Vec::new();
    other.identity_words(&mut identity_b);
    assert_ne!(identity_a, identity_b);
    assert!(matches!(
        wrapped.admit(&START, &START, &Receipt::lossless(0.0, 0.0)),
        Err(ObservationError::Owner("native admit"))
    ));
    let gate = CancelGate::new_clock_free();
    let untouched = StateSignals::new(
        Described(SpringMassDamper::new(Fault::None)),
        signals(),
        |_: &Described, _: &Oscillator, _: &mut [f64], _: &Cx<'_>| {
            panic!("observer cannot run on refusal/shortening")
        },
    )
    .unwrap();
    with_cx(&gate, |cx| {
        assert!(matches!(
            untouched.propose(&START, 0.5, &INPUTS, &mut STALE.clone(), &mut [9.0; 4], cx),
            Ok(Proposal::Shorten { dt_s: 0.2 })
        ));
        assert!(matches!(
            untouched.propose(
                &START,
                0.01,
                &[-1.0, 300.0],
                &mut STALE.clone(),
                &mut [9.0; 4],
                cx
            ),
            Err(ObservationError::Owner(
                "negative force refused by native owner"
            ))
        ));
    });
}

impl StageResidual for Described {
    fn unknowns(&self) -> usize {
        2
    }
    fn initial_guess(
        &self,
        old: &Oscillator,
        _: f64,
        _: &[f64],
        z: &mut [f64],
    ) -> Result<(), Self::Error> {
        z.copy_from_slice(&[old.position_m, old.velocity_m_s]);
        Ok(())
    }
    fn residual(
        &self,
        old: &Oscillator,
        dt: f64,
        inputs: &[f64],
        z: &[f64],
        r: &mut [f64],
        out: &mut [f64],
        _: &Cx<'_>,
    ) -> Result<(), Self::Error> {
        let v = 0.5 * (old.velocity_m_s + z[1]);
        r[0] = z[0] - old.position_m - dt * v;
        r[1] = self.0.mass_kg * (z[1] - old.velocity_m_s)
            - dt * (inputs[0]
                - self.0.stiffness_n_m * 0.5 * (old.position_m + z[0])
                - self.0.damping_n_s_m * v);
        out.copy_from_slice(&[v, -self.0.damping_n_s_m * v * v]);
        Ok(())
    }
    fn jvp_z(
        &self,
        old: &Oscillator,
        dt: f64,
        _: &[f64],
        z: &[f64],
        dz: &[f64],
        dr: &mut [f64],
        out: &mut [f64],
        _: &Cx<'_>,
    ) -> Result<(), Self::Error> {
        dr[0] = dz[0] - 0.5 * dt * dz[1];
        dr[1] = 0.5 * dt * self.0.stiffness_n_m * dz[0]
            + (self.0.mass_kg + 0.5 * dt * self.0.damping_n_s_m) * dz[1];
        out.copy_from_slice(&[
            0.5 * dz[1],
            -self.0.damping_n_s_m * 0.5 * (old.velocity_m_s + z[1]) * dz[1],
        ]);
        Ok(())
    }
    fn vjp_z(
        &self,
        old: &Oscillator,
        dt: f64,
        _: &[f64],
        z: &[f64],
        bar_r: &[f64],
        bar_out: &[f64],
        bar_z: &mut [f64],
        _: &Cx<'_>,
    ) -> Result<(), Self::Error> {
        bar_z[0] = bar_r[0] + 0.5 * dt * self.0.stiffness_n_m * bar_r[1];
        bar_z[1] = -0.5 * dt * bar_r[0]
            + (self.0.mass_kg + 0.5 * dt * self.0.damping_n_s_m) * bar_r[1]
            + 0.5 * bar_out[0]
            - self.0.damping_n_s_m * 0.5 * (old.velocity_m_s + z[1]) * bar_out[1];
        Ok(())
    }
    fn accept(
        &self,
        old: &Oscillator,
        _: f64,
        _: &[f64],
        z: &[f64],
        new: &mut Oscillator,
        _: &Cx<'_>,
    ) -> Result<Receipt, Self::Error> {
        *new = Oscillator {
            position_m: z[0],
            velocity_m_s: z[1],
        };
        Ok(Receipt::lossless(
            self.stored_energy_j(new)? - self.stored_energy_j(old)?,
            1e-12,
        ))
    }
}

#[test]
fn residual_tier_retains_original_unknowns_and_both_derivative_actions() {
    let wrapped = StateSignals::new(
        Described(SpringMassDamper::new(Fault::None)),
        signals(),
        described_observation,
    )
    .unwrap();
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        assert_eq!(wrapped.unknowns(), 2);
        let mut z = [0.0; 2];
        wrapped
            .initial_guess(&START, 0.01, &INPUTS, &mut z)
            .unwrap();
        assert_eq!(z, [START.position_m, START.velocity_m_s]);
        z = [0.003, -0.2];
        let (mut r, mut expected_r) = ([0.0; 2], [0.0; 2]);
        let (mut out, mut expected_out) = ([9.0; 4], [0.0; 2]);
        wrapped
            .residual(&START, 0.01, &INPUTS, &z, &mut r, &mut out, cx)
            .unwrap();
        wrapped
            .inner()
            .residual(
                &START,
                0.01,
                &INPUTS,
                &z,
                &mut expected_r,
                &mut expected_out,
                cx,
            )
            .unwrap();
        assert_eq!(r, expected_r);
        assert_eq!(&out[..2], &expected_out);
        assert_eq!(&out[2..], &[START.position_m, START.velocity_m_s]);
        wrapped
            .jvp_z(
                &START,
                0.01,
                &INPUTS,
                &z,
                &[0.1, -0.3],
                &mut r,
                &mut out,
                cx,
            )
            .unwrap();
        wrapped
            .inner()
            .jvp_z(
                &START,
                0.01,
                &INPUTS,
                &z,
                &[0.1, -0.3],
                &mut expected_r,
                &mut expected_out,
                cx,
            )
            .unwrap();
        assert_eq!(r, expected_r);
        assert_eq!(&out[..2], &expected_out);
        assert_eq!(&out[2..], &[0.0; 2]);
        wrapped
            .vjp_z(
                &START,
                0.01,
                &INPUTS,
                &z,
                &[0.7, -0.9],
                &[1.0, -2.0, 1e30, -1e30],
                &mut r,
                cx,
            )
            .unwrap();
        wrapped
            .inner()
            .vjp_z(
                &START,
                0.01,
                &INPUTS,
                &z,
                &[0.7, -0.9],
                &[1.0, -2.0],
                &mut expected_r,
                cx,
            )
            .unwrap();
        assert_eq!(r, expected_r);
        let (mut state, mut expected) = (STALE, STALE);
        assert_eq!(
            wrapped
                .accept(&START, 0.01, &INPUTS, &z, &mut state, cx)
                .unwrap(),
            wrapped
                .inner()
                .accept(&START, 0.01, &INPUTS, &z, &mut expected, cx)
                .unwrap()
        );
        assert_eq!(state, expected);
    });
}
