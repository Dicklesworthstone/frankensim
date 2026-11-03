//! Opt-in jet drive of Auvray, Fabre and Lagrée (2012, JASA 131(2) 1574) as used by Auvray, Ernoult,
//! Fabre and Lagrée (2014, JASA 136(1) 389), in place of the island's band chain and authored gain.
//!
//! Card row `jet_transfer auvray-2014 BETA GAMMA SOURCE_TAG`. Taken as a whole, it is the model
//! those papers validated against measured oscillation thresholds:
//! - receptivity at the flue exit, IN PHASE with the transverse acoustic velocity
//!   (de la Cuadra; 2012 Eq. 9, 2014 Eq. 4): `eta_0(t) = h v_ac(t) / U`;
//! - flat spatial growth and pure convection (2014 Eq. 3, `alpha_i = beta / h`, `c_p = gamma U`):
//!   `eta(W, t) = exp(beta W / h) eta_0(t - W / (gamma U))`;
//! - a Bickley profile at the labium with half-width `b = 2h/5`, from momentum conservation
//!   between a Poiseuille exit and a Bickley jet of equal centreline speed (2014, after Eq. 5):
//!   the flow into the pipe is `Q_in = H b U (1 + tanh((eta - y_off) / b))`, a swing of `2 b H U`;
//! - the jet-drive source `Delta p = (rho delta_d / (W H)) dQ_in/dt` with Verge's effective
//!   distance between the two injection points `delta_d = (4/pi) sqrt(2 h W)` (2014 Eq. 5);
//! - `U = sqrt(2 P / rho)`, Bernoulli on the reservoir (2014 Eq. 21). The row therefore refuses
//!   every windway law: it states its own exit speed.
//!
//! `beta` and `gamma` come from the card, never from a default here. The mouth's separation loss
//! that is always on in Auvray's model (2014 Eq. 2, `alpha_vc`) is the card's own `mouth_loss`
//! row; this module does not add it.
//!
//! Citation for the card values: Auvray, Ernoult, Fabre and Lagrée, "Jet-drive/discrete-vortex
//! models' simulation of a recorder", J. Acoust. Soc. Am. 136(1), 389-400 (2014). Its Table I
//! (row JD) gives `b = 0.3` and `c_jd = 0.4`, i.e. `beta 0.3` and `gamma 0.4` on this row; Sec. II B
//! states them as the experimental values of de la Cuadra and of Blanc (`beta ~ 0.3`,
//! `c_jd ~ 0.4`), which depend on the channel's velocity profile and hold over a few jet
//! thicknesses `h`. Both values are from the paper, not estimates made here. The paper validated
//! the model on one Zen-On Bressan alto recorder with `W` 4 mm (`W/h` about 4); a card at a much
//! larger `W/h` (a flute lip jet, `W/h` 6 to 9) uses them outside that validation and should say so.

/// Bickley half-width over flue height, `b / h = 2/5` (Auvray et al. 2014, after Eq. 5).
pub const BICKLEY_HALF_WIDTH_OVER_HEIGHT: f64 = 0.4;

/// Why the row and a windway law cannot share one voice, in either attachment order.
pub const WINDWAY_REFUSAL: &str =
    "the Auvray jet takes U = sqrt(2P/rho) at the exit (Auvray 2014 Eq. 21) and refuses every windway law";

/// Card values of the row: the growth coefficient `beta` (`alpha_i = beta / h`), the convection
/// ratio `gamma` (`c_p = gamma U`) and the explicit source tag the card cites them by.
#[derive(Debug, Clone, PartialEq)]
pub struct AuvraySpec {
    /// Spatial growth coefficient: `alpha_i h`.
    pub beta: f64,
    /// Convection velocity over jet centreline velocity.
    pub gamma: f64,
    /// The card's citation token for both values.
    pub source: String,
}

impl AuvraySpec {
    /// Validate the row's fields. Nothing is defaulted.
    ///
    /// # Errors
    /// A non-finite or non-positive `beta`, a `gamma` outside `(0, 1]`, or an empty source tag.
    pub fn new(beta: f64, gamma: f64, source: &str) -> Result<Self, &'static str> {
        if !(beta.is_finite() && beta > 0.0) {
            return Err("Auvray jet growth coefficient beta must be finite and positive");
        }
        if !(gamma.is_finite() && gamma > 0.0 && gamma <= 1.0) {
            return Err("Auvray jet convection ratio gamma must lie in (0, 1]");
        }
        if source.trim().is_empty() {
            return Err("Auvray jet needs an explicit source tag");
        }
        Ok(Self { beta, gamma, source: source.to_owned() })
    }
}

/// The row bound to one flue geometry: every constant of the closed-form transfer.
#[derive(Debug, Clone, PartialEq)]
pub struct AuvrayJet {
    spec: AuvraySpec,
    flue_height_m: f64,
    cut_up_m: f64,
    growth: f64,
    half_width_m: f64,
    source_length_m: f64,
}

impl AuvrayJet {
    /// Bind the row to the flue height `h` and cut-up `W`.
    ///
    /// # Errors
    /// A non-positive or non-finite height or cut-up, or a growth `exp(beta W / h)` that is not finite.
    pub fn new(spec: AuvraySpec, flue_height_m: f64, cut_up_m: f64) -> Result<Self, &'static str> {
        if !(flue_height_m.is_finite() && flue_height_m > 0.0 && cut_up_m.is_finite() && cut_up_m > 0.0) {
            return Err("Auvray jet needs a finite positive flue height and cut-up");
        }
        let growth = (spec.beta * cut_up_m / flue_height_m).exp();
        let half_width_m = BICKLEY_HALF_WIDTH_OVER_HEIGHT * flue_height_m;
        let source_length_m = 4.0 / core::f64::consts::PI * (2.0 * flue_height_m * cut_up_m).sqrt();
        if !(growth.is_finite() && half_width_m > 0.0 && source_length_m.is_finite() && source_length_m > 0.0) {
            return Err("Auvray jet growth or source length is not finite");
        }
        Ok(Self { spec, flue_height_m, cut_up_m, growth, half_width_m, source_length_m })
    }

    /// The card row.
    #[must_use]
    pub const fn spec(&self) -> &AuvraySpec { &self.spec }

    /// Cut-up over flue height `W / h`, the station the fitted `beta` and `gamma` apply at.
    #[must_use]
    pub fn cut_up_over_height(&self) -> f64 { self.cut_up_m / self.flue_height_m }

    /// Spatial growth to the labium, `exp(beta W / h)`: a magnitude, no phase.
    #[must_use]
    pub const fn growth(&self) -> f64 { self.growth }

    /// Bickley half-width at the labium, `b = 2h/5` [m].
    #[must_use]
    pub const fn half_width_m(&self) -> f64 { self.half_width_m }

    /// Verge's source length `delta_d = (4/pi) sqrt(2 h W)` [m].
    #[must_use]
    pub const fn source_length_m(&self) -> f64 { self.source_length_m }

    /// Labium displacement per unit retarded transverse velocity, `h exp(beta W / h) / U` [s];
    /// zero for a jet at rest.
    #[must_use]
    pub fn receptivity_s(&self, jet_speed_m_s: f64) -> f64 {
        if jet_speed_m_s > 0.0 { self.flue_height_m * self.growth / jet_speed_m_s } else { 0.0 }
    }

    /// Convection delay to the labium, `W / (gamma U)` [s]; infinite for a jet at rest.
    #[must_use]
    pub fn transit_s(&self, jet_speed_m_s: f64) -> f64 {
        self.cut_up_m / (self.spec.gamma * jet_speed_m_s)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::{FlueControl, FlueGeometry, FlueVoice, forced_jet, network, nolle_jet};
    use crate::bernoulli_aperture::network::{NetworkNode, TubeNetworkSpec, TubeSection};
    use fs_exec::CancelGate;
    use fs_material::gas::{GasSpec, GasState};

    const TAG: &str = "numerical-fixture";

    fn spec(beta: f64, gamma: f64) -> AuvraySpec { AuvraySpec::new(beta, gamma, TAG).unwrap() }

    // Recorder-like station (W/h 4.34), the committed recorder cards' flue row. A numerical fixture.
    fn geometry() -> FlueGeometry {
        FlueGeometry { flue_width_m: 0.012, flue_height_m: 1.059e-3, cut_up_m: 4.6e-3,
            labium_offset_m: 1.455e-4, jet_angle_rad: 0.0,
            bore_area_m2: core::f64::consts::PI * 0.006 * 0.006, mouth_loss_rel: 0.02 }
    }
    fn gas() -> GasState {
        GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap()
    }
    fn voice(rate: u32) -> FlueVoice {
        voice_with(rate, geometry())
    }
    fn voice_with(rate: u32, geometry: FlueGeometry) -> FlueVoice {
        let gas = gas();
        let source = network::FlueNetworkSpec {
            graph: TubeNetworkSpec {
                nodes: vec![NetworkNode::Inlet, NetworkNode::Termination { reflection: -0.7 }],
                sections: vec![TubeSection { nodes: [0, 1], length_m: 0.3, radius_m: 0.006,
                    max_length_error_m: 0.005 }],
                sound_speed_m_s: gas.sound_speed, max_wave_memory_bytes: 1_000_000,
            },
            gases: vec![gas], viscothermal: vec![], minimum_jet_speed_m_s: 1.0,
            maximum_jet_mach: 0.3, max_jet_history_samples: 4096,
        };
        FlueVoice::from_network(Some(&fs_aeroac::jetcard::mint_tonal_interim_card().unwrap()),
            90.0, 1e-6, geometry, source, rate, 3, 0, &CancelGate::new_clock_free()).unwrap()
    }

    #[test]
    fn auvray_spec_refuses_unusable_fields_with_text_and_defaults_nothing() {
        const BETA: &str = "Auvray jet growth coefficient beta must be finite and positive";
        const GAMMA: &str = "Auvray jet convection ratio gamma must lie in (0, 1]";
        const SOURCE: &str = "Auvray jet needs an explicit source tag";
        for (beta, gamma, tag, text) in [(0.0, 0.4, TAG, BETA), (-0.0, 0.4, TAG, BETA), (-0.1, 0.4, TAG, BETA),
            (f64::NAN, 0.4, TAG, BETA), (f64::INFINITY, 0.4, TAG, BETA), (0.3, 0.0, TAG, GAMMA),
            (0.3, -0.4, TAG, GAMMA), (0.3, 1.5, TAG, GAMMA), (0.3, f64::NAN, TAG, GAMMA), (0.3, 0.4, "", SOURCE),
            (0.3, 0.4, "  ", SOURCE)] {
            assert_eq!(AuvraySpec::new(beta, gamma, tag), Err(text), "{beta} {gamma} {tag:?}");
        }
        // The paper's values (Auvray et al. 2014 Table I) and the edge gamma = 1 are accepted as given.
        let paper = AuvraySpec::new(0.3, 0.4, TAG).unwrap();
        assert_eq!((paper.beta, paper.gamma, paper.source.as_str()), (0.3, 0.4, TAG));
        assert_eq!(AuvraySpec::new(1e-9, 1.0, TAG).unwrap().gamma, 1.0);
        assert!(AuvrayJet::new(spec(0.3, 0.4), 0.0, 4.6e-3).is_err());
        assert!(AuvrayJet::new(spec(0.3, 0.4), 1e-3, f64::NAN).is_err());
        assert!(AuvrayJet::new(spec(1e4, 0.4), 1e-3, 1.0).is_err(), "an overflowing growth is refused");
    }

    #[test]
    fn auvray_constants_are_the_published_closed_forms_at_the_flute_station() {
        // Flute D4 card: h 0.9 mm, W 8 mm; beta 0.3. Literals from exp(8/3), 0.4 h and
        // (4/pi) sqrt(2 x 0.9e-3 x 8e-3), worked by hand, not by this module.
        let jet = AuvrayJet::new(spec(0.3, 0.4), 0.9e-3, 8e-3).unwrap();
        assert!((jet.growth() / 14.391_916_095_149_892 - 1.0).abs() < 1e-12);
        assert!((jet.half_width_m() / 3.6e-4 - 1.0).abs() < 1e-12);
        assert!((jet.source_length_m() / 4.831_604_362_030_756e-3 - 1.0).abs() < 1e-12);
        assert!((jet.cut_up_over_height() / (80.0 / 9.0) - 1.0).abs() < 1e-12);
        // In phase and frequency-flat: the receptivity is a real number, h e^(beta W/h) / U.
        assert!((jet.receptivity_s(20.0) / (0.9e-3 * 14.391_916_095_149_892 / 20.0) - 1.0).abs() < 1e-12);
        assert_eq!(jet.receptivity_s(0.0), 0.0);
        assert!((jet.transit_s(20.0) / (8e-3 / 8.0) - 1.0).abs() < 1e-12);
        assert!(jet.transit_s(0.0).is_infinite());
    }

    // The closed-form check of the row's transfer at two frequencies: a voice whose mouth-velocity
    // history is a pure tone v = A sin(w t) must cut the jet at eta_n = R A sin(w (t_n - tau)) with
    // R = h exp(beta W/h) / U and tau = W / (gamma U), U = sqrt(2 P / rho), both taken from the
    // published forms here, and the flow into the pipe must be H b U (1 + tanh((eta - y0)/b)),
    // b = 2h/5. At 1200 Hz a one-sample delay error moves eta by w dt = 0.16 of R A, fifty times
    // the linear-interpolation bound the analytic check allows; the exact check compares with the
    // linear interpolation of the same samples to rounding. The island chain's own band lag, its
    // frequency-dependent magnitude or its half-width h/2 each fail one of the three.
    #[test]
    fn auvray_transfer_is_flat_in_phase_and_delayed_by_w_over_gamma_u_at_two_frequencies() {
        let rate = 48_000;
        let (beta, gamma, pressure, amplitude) = (0.3, 0.4, 415.0, 0.5);
        let g = geometry();
        for frequency in [300.0, 1200.0] {
            let mut v = voice(rate).with_auvray_jet(spec(beta, gamma)).unwrap();
            v.island.seed_amplitude_rel = 0.0;
            v.apply(FlueControl::SetBlowingPressure(pressure)).unwrap();
            let dt = 1.0 / f64::from(rate);
            let u = (2.0 * pressure / gas().density).sqrt();
            assert_eq!(v.jet_speed_m_s(), u);
            let receptivity = g.flue_height_m * (beta * g.cut_up_m / g.flue_height_m).exp() / u;
            let tau = g.cut_up_m / (gamma * u);
            assert!((v.transit_s() / tau - 1.0).abs() < 1e-14);
            assert!((v.source_coeff / (gas().density * 4.0 / core::f64::consts::PI
                * (2.0 * g.flue_height_m * g.cut_up_m).sqrt() / g.mouth_area_m2()) - 1.0).abs() < 1e-14);
            let b = 0.4 * g.flue_height_m;
            assert!((v.cut_half_width_m() / b - 1.0).abs() < 1e-15);
            let omega = 2.0 * core::f64::consts::PI * frequency;
            let tone = |samples: f64| amplitude * (omega * samples * dt).sin();
            let len = v.history.len();
            assert_eq!(v.head, 0);
            for (j, slot) in v.history.iter_mut().enumerate() {
                *slot = tone(-((len - j) as f64));
            }
            v.history_moving = len as u64;
            let delay = tau / dt;
            assert!(delay > 18.0, "the fixture reads only prefilled history");
            let back = delay - 1.0;
            let whole = back.floor();
            let fraction = back - whole;
            for n in 0..16 {
                let mut seen = None;
                v.step_block_labium_observed(&mut [0.0], |_, labium| { seen = Some(labium); Ok(()) }).unwrap();
                let eta = seen.expect("one observed sample").eta_m;
                let n = f64::from(n);
                let interpolated = (1.0 - fraction) * tone(n - 1.0 - whole) + fraction * tone(n - 2.0 - whole);
                let exact = receptivity * interpolated;
                let analytic = receptivity * amplitude * (omega * (n * dt - tau)).sin();
                assert!((eta - exact).abs() <= 1e-10 * receptivity * amplitude,
                    "{frequency} Hz sample {n}: eta {eta} against {exact}");
                let bound = receptivity * amplitude * (omega * dt).powi(2) / 8.0;
                assert!((eta - analytic).abs() <= 1.01 * bound + 1e-12 * receptivity * amplitude,
                    "{frequency} Hz sample {n}: eta {eta} against the delayed tone {analytic}");
                let flow = b * g.flue_width_m * u * (1.0 + ((eta - g.effective_offset_m()) / b).tanh());
                assert!((v.flow_prev - flow).abs() <= 1e-12 * flow, "{frequency} Hz sample {n}");
            }
        }
    }

    // MB 15924 (owner order): the row's cut keeps the card's labium offset, (eta - y0) / b, and the offset is
    // what makes even harmonics. A tone v = A sin(w t) with w = 2 pi rate / 16 (one period is 16 samples, all
    // read from prefilled history as in the test above) drives the cut. Over that period the second harmonic of
    // the flow into the pipe, over its first, must equal the closed form's on the same eta,
    // H b U (1 + tanh((eta - y0) / b)), and be far from zero: about 0.13 here, at y0/b 0.34 and eta/b 1.4.
    // Planted negative: the same voice with y0 = 0. The cut is then odd about its mean, so every even harmonic
    // vanishes to rounding.
    #[test]
    fn auvray_cut_keeps_the_labium_offset_and_only_an_offset_makes_even_harmonics() {
        let rate = 48_000;
        let (beta, gamma, pressure, amplitude) = (0.3, 0.4, 415.0, 4.0);
        let dt = 1.0 / f64::from(rate);
        let omega = 2.0 * core::f64::consts::PI * f64::from(rate) / 16.0;
        let harmonic = |flows: &[f64], k: f64| {
            let (mut re, mut im) = (0.0_f64, 0.0_f64);
            for (n, q) in flows.iter().enumerate() {
                let phase = k * omega * n as f64 * dt;
                re += q * phase.cos();
                im -= q * phase.sin();
            }
            re.hypot(im)
        };
        let mut ratios = Vec::new();
        for offset in [geometry().labium_offset_m, 0.0] {
            let g = FlueGeometry { labium_offset_m: offset, ..geometry() };
            let mut v = voice_with(rate, g).with_auvray_jet(spec(beta, gamma)).unwrap();
            v.island.seed_amplitude_rel = 0.0;
            v.apply(FlueControl::SetBlowingPressure(pressure)).unwrap();
            let u = v.jet_speed_m_s();
            let receptivity = g.flue_height_m * (beta * g.cut_up_m / g.flue_height_m).exp() / u;
            let (tau, b) = (g.cut_up_m / (gamma * u), 0.4 * g.flue_height_m);
            let tone = |samples: f64| amplitude * (omega * samples * dt).sin();
            let len = v.history.len();
            assert_eq!(v.head, 0);
            for (j, slot) in v.history.iter_mut().enumerate() {
                *slot = tone(-((len - j) as f64));
            }
            v.history_moving = len as u64;
            let back = tau / dt - 1.0;
            assert!(back > 17.0, "the fixture reads only prefilled history");
            let (whole, fraction) = (back.floor(), back - back.floor());
            let (mut voice_flows, mut closed_flows) = (Vec::new(), Vec::new());
            for n in 0..16 {
                v.step_block(&mut [0.0]).unwrap();
                voice_flows.push(v.flow_prev);
                let n = f64::from(n);
                let eta = receptivity
                    * ((1.0 - fraction) * tone(n - 1.0 - whole) + fraction * tone(n - 2.0 - whole));
                closed_flows.push(b * g.flue_width_m * u * (1.0 + ((eta - offset) / b).tanh()));
            }
            let ratio = harmonic(&voice_flows, 2.0) / harmonic(&voice_flows, 1.0);
            let closed = harmonic(&closed_flows, 2.0) / harmonic(&closed_flows, 1.0);
            println!("y0 {offset:e} m (y0/b {:.3}): flow H2/H1 {ratio:e}, closed form {closed:e}", offset / b);
            ratios.push((ratio, closed));
        }
        let ((with_offset, closed), (without, _)) = (ratios[0], ratios[1]);
        assert!(with_offset > 0.05, "the card's offset must make a second harmonic: H2/H1 {with_offset:e}");
        assert!((with_offset / closed - 1.0).abs() < 1e-9,
            "the cut must be (eta - y0)/b: H2/H1 {with_offset:e} against the closed form's {closed:e}");
        assert!(without < 1e-12, "planted negative y0 = 0: H2/H1 {without:e} must vanish");
    }

    #[test]
    fn auvray_jet_refuses_late_or_repeated_attachment_other_jets_windways_and_short_history() {
        let row = || spec(0.3, 0.4);
        assert!(voice(48_000).with_auvray_jet(row()).unwrap().with_auvray_jet(row()).is_err());
        let mut warm = voice(48_000);
        warm.apply(FlueControl::SetBlowingPressure(415.0)).unwrap();
        warm.step_block(&mut [0.0]).unwrap();
        assert!(warm.with_auvray_jet(row()).is_err());
        // Beside the forced Bickley jet and the Nolle/Price transfer, in both orders.
        let forced = voice(48_000).with_forced_jet(forced_jet::COARSE_SPEC, 300.0).unwrap();
        assert!(forced.with_auvray_jet(row()).is_err());
        assert!(voice(48_000).with_auvray_jet(row()).unwrap()
            .with_forced_jet(forced_jet::COARSE_SPEC, 300.0).is_err());
        let nolle = nolle_jet::NolleSpec::price_2015(nolle_jet::ChannelClass::Short, 0.0007).unwrap();
        assert!(voice(48_000).with_nolle_jet(nolle).unwrap().with_auvray_jet(row()).is_err());
        assert!(voice(48_000).with_auvray_jet(row()).unwrap().with_nolle_jet(nolle).is_err());
        // U = sqrt(2 P / rho): every windway law and the card centreline speed are refused, in both orders.
        assert!(voice(48_000).with_viscous_windway(0.015).unwrap().with_auvray_jet(row()).is_err());
        assert!(voice(48_000).with_auvray_jet(row()).unwrap().with_viscous_windway(0.015).is_err());
        assert!(voice(48_000).with_developing_windway(0.015, 100.0).unwrap().with_auvray_jet(row()).is_err());
        assert!(voice(48_000).with_auvray_jet(row()).unwrap().with_developing_windway(0.015, 100.0).is_err());
        assert!(voice(48_000).with_unsteady_windway(0.015, 8, 8).unwrap().with_auvray_jet(row()).is_err());
        assert!(voice(48_000).with_auvray_jet(row()).unwrap().with_unsteady_windway(0.015, 8, 8).is_err());
        assert!(voice(48_000).with_card_centreline_speed(415.0).unwrap().with_auvray_jet(row()).is_err());
        assert!(voice(48_000).with_auvray_jet(row()).unwrap().with_card_centreline_speed(415.0).is_err());
        // The retained history must cover W / (gamma U_min): a slower convection grows it while the
        // voice is cold, and one past the declared budget is refused.
        let before = voice(48_000).history.len();
        let slow = voice(48_000).with_auvray_jet(spec(0.3, 0.2)).unwrap();
        let needed = (geometry().cut_up_m / (0.2 * 1.0) / (1.0 / 48_000.0)).ceil() as usize + 2;
        assert!(needed > before);
        assert_eq!(slow.history.len(), needed);
        assert!(slow.history.iter().all(|x| *x == 0.0));
        assert_eq!(voice(48_000).with_auvray_jet(row()).unwrap().history.len(), before);
        assert!(voice(48_000).with_auvray_jet(spec(0.3, 0.01)).is_err());
        // Without the row nothing changes: the half-width and the source coefficient are the island's.
        let plain = voice(48_000);
        assert!(plain.auvray_jet().is_none());
        assert_eq!(plain.cut_half_width_m(), 0.5 * geometry().flue_height_m);
        assert_eq!(plain.source_coeff, plain.rho * geometry().end_correction_m() / geometry().mouth_area_m2());
    }
}
