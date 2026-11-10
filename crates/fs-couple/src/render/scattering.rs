//! Couple the original moving reed to an affine same-sample acoustic boundary.
//!
//! For incoming b = r*a + h and original impedance Z, pressure and flow obey
//! p = (1+r)*a+h and Z*q = (1-r)*a-h. For -1<r<1 this is the same physical
//! load as p = Z_eff*q + 2*h_eff, with Z_eff = Z*(1+r)/(1-r) and
//! h_eff = h/(1-r). The existing reed/contact owner therefore needs no new
//! integration scheme, force law or fictitious delay. Only wave coordinates
//! change; the source gas, pressure area, body flow and material state do not.

use super::*;
use fs_vfit::scattering::ScatteringResponse;

impl ReedBoreVoice {
    // The massless owner already accepts a direct reflection. Its massive
    // counterpart accepts a delayed wave; use an equivalent load coordinate
    // instead of freezing the incoming wave at the previous accepted sample.
    pub(super) fn step_massive_scattering(
        &mut self,
        boundary: ScatteringResponse,
        mouth_pressure_pa: f64,
        body_flow_m3_s: f64,
    ) -> Result<(f64, f64, f64), AcousticRealizeError> {
        let r = boundary.direct;
        if !r.is_finite() || r <= -1.0 || r >= 1.0 || !boundary.history.is_finite() {
            return Err(AcousticRealizeError::Reed {
                what: "massive reed requires a finite strictly passive instantaneous wave boundary",
            });
        }
        // Preserve the original arithmetic when the boundary has no direct term.
        let (impedance, incoming) = if r == 0.0 {
            (self.zc, boundary.history)
        } else {
            (self.zc * ((1.0 + r) / (1.0 - r)), boundary.history / (1.0 - r))
        };
        if !impedance.is_finite() || impedance <= 0.0 || !incoming.is_finite() {
            return Err(AcousticRealizeError::Reed {
                what: "same-sample reed load is outside the finite positive-impedance domain",
            });
        }
        let closing_force_n = self.closing_force_n;
        let closing_force = |_next_y: f64| Ok(closing_force_n);
        let material_force: Option<&dyn Fn(f64) -> Result<f64, AcousticRealizeError>> =
            if closing_force_n != 0.0 {
                Some(&closing_force)
            } else {
                None
            };
        // The midpoint step with its rigid closure stop at zero opening (the default) or without it (`closure_stop off`).
        let midpoint = if self.closure_stop {
            crate::reed_bore::OpeningPropagator::Midpoint
        } else {
            crate::reed_bore::OpeningPropagator::MidpointFree
        };
        let (outgoing, opening, velocity) = if let Some(slit) = self.viscous_slit {
            let lay = self.lay.as_ref().expect("massive voices retain their physical lay");
            let result = if let Some(recovery) = self.staple_recovery {
                // Almeida, Vergez and Caussé 2007: the flow follows the slit law on the constriction drop (dp)_c, the
                // blade is driven by (dp)_c (eq. 16), and the bore input is raised by the recovered head (eq. 17).
                midpoint.step_with_junction(
                    self.reed, self.rho, impedance, incoming, mouth_pressure_pa,
                    self.reed_y, self.reed_v, self.dt, body_flow_m3_s, lay,
                    |opening| Ok(opening.max(0.0)), material_force,
                    |opening, moving| {
                        let junction = slit.characteristic_junction_with_recovery(
                            recovery, self.reed.width_m, opening, self.rho, impedance,
                            incoming, mouth_pressure_pa, moving,
                        )?;
                        Ok(crate::reed_bore::JunctionPressure {
                            outgoing: junction.outgoing_pa,
                            blade_drop_pa: Some(junction.constriction_drop_pa),
                        })
                    },
                )?
            } else {
                midpoint.step_with_junction(
                    self.reed, self.rho, impedance, incoming, mouth_pressure_pa,
                    self.reed_y, self.reed_v, self.dt, body_flow_m3_s, lay,
                    |opening| Ok(opening.max(0.0)), material_force,
                    |opening, moving| slit.characteristic_pressure(
                        self.reed.width_m, opening, self.rho, impedance, incoming,
                        mouth_pressure_pa, moving,
                    ),
                )?
            };
            if self.solver_mode == ReedSolverMode::FastNewton {
                self.fast_stats.fallback_samples += 1;
            }
            result
        } else if self.staple_recovery.is_some()
            || self.jet_contraction.is_some()
            || self.residual_flow_area_m2.is_some()
        {
            // The same staple law on the inviscid Bernoulli junction, for the cards without `slit_length`: Almeida's own
            // model neglects the viscous loss along the reed. Same midpoint step, lay and tongue force as the branches below.
            // A jet contraction `C_c` (alone, or with the staple) narrows the jet to `C_c w h`; without a staple the recovery
            // is the `C_P = 0` one, whose junction is the plain Bernoulli one and whose blade drop is the whole mouth-to-bore drop.
            let recovery = match self.staple_recovery {
                Some(recovery) => recovery,
                None => crate::bernoulli_aperture::viscous::StapleRecovery::new(0.0, 1.0)?,
            };
            let contraction = self.jet_contraction.unwrap_or(1.0);
            // The residual flow (`residual_flow`, Almeida, Vergez and Caussé 2007 Fig. 9): the junction sees the opening
            // `max(h, 0) + S / (C_c w)`, so its jet area is `C_c w max(h, 0) + S` and the flow never falls to zero. The blade's
            // own mechanics, lay contact and stop keep the blade opening. With no row there is no floor: the original closure.
            let residual_opening = self.residual_flow_area_m2.map(|area| area / (contraction * self.reed.width_m));
            let lay = self.lay.as_ref().expect("massive voices retain their physical lay");
            let result = midpoint.step_with_junction(
                self.reed, self.rho, impedance, incoming, mouth_pressure_pa,
                self.reed_y, self.reed_v, self.dt, body_flow_m3_s, lay,
                |opening| Ok(match residual_opening {
                    Some(floor) => opening.max(0.0) + floor,
                    None => opening.max(0.0),
                }),
                material_force,
                |opening, moving| {
                    let junction = crate::bernoulli_aperture::moving::characteristic_junction_with_recovery(
                        recovery, contraction, self.reed.width_m, opening, self.rho, impedance,
                        incoming, mouth_pressure_pa, moving,
                    )?;
                    Ok(crate::reed_bore::JunctionPressure {
                        outgoing: junction.outgoing_pa,
                        blade_drop_pa: Some(junction.constriction_drop_pa),
                    })
                },
            )?;
            if self.solver_mode == ReedSolverMode::FastNewton {
                self.fast_stats.fallback_samples += 1;
            }
            result
        } else if closing_force_n != 0.0 {
            let lay = self.lay.as_ref().expect("massive voices retain their physical lay");
            let result = crate::reed_bore::step_profiled_aperture(
                self.reed, self.rho, impedance, incoming, mouth_pressure_pa,
                self.reed_y, self.reed_v, self.dt, body_flow_m3_s, lay,
                |opening| Ok(opening.max(0.0)), material_force,
                midpoint,
            )?;
            if self.solver_mode == ReedSolverMode::FastNewton {
                self.fast_stats.fallback_samples += 1;
            }
            result
        } else if !self.closure_stop {
            // `step_massive_reed` with the lay is `contact::step`, the midpoint step with the stop; with `closure_stop off` the
            // same arithmetic runs through the profiled-aperture entry with the stop-free propagator.
            let lay = self.lay.as_ref().expect("massive voices retain their physical lay");
            let result = crate::reed_bore::step_profiled_aperture(
                self.reed, self.rho, impedance, incoming, mouth_pressure_pa,
                self.reed_y, self.reed_v, self.dt, body_flow_m3_s, lay,
                |opening| Ok(opening.max(0.0)), None,
                midpoint,
            )?;
            if self.solver_mode == ReedSolverMode::FastNewton {
                self.fast_stats.fallback_samples += 1;
            }
            result
        } else {
            crate::reed_bore::step_massive_reed(
                self.reed, self.rho, impedance, incoming, mouth_pressure_pa,
                self.reed_y, self.reed_v, self.dt, body_flow_m3_s,
                self.lay.as_ref(), self.solver_mode, &mut self.fast_stats,
            )?
        };
        // p = a_eff + h_eff = (1+r)*a + h, and h=(1-r)*h_eff,
        // hence a=(a_eff+r*h_eff)/(1+r). The history term adds, not subtracts.
        let outgoing = if r == 0.0 { outgoing }
            else { (outgoing + r * incoming) / (1.0 + r) };
        if !outgoing.is_finite() {
            return Err(AcousticRealizeError::Reed {
                what: "same-sample reed outgoing wave is not finite",
            });
        }
        Ok((outgoing, opening, velocity))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_material::gas::GasSpec;
    use fs_vfit::discretize::DelayedFilter;

    fn voice(direct: f64, massive: bool) -> ReedBoreVoice {
        voice_with(direct, BeatingReed {
            rest_opening_m: 0.0004, width_m: 0.013,
            closing_pressure_pa: 6000.0, blowing_pressure_pa: 100.0,
            attack_s: 0.0, mass_kg: if massive { 0.00001 } else { 0.0 },
            stiffness_n_m: if massive { 100.0 } else { 0.0 }, damping_ratio: 0.35,
        })
    }

    fn voice_with(direct: f64, reed: BeatingReed) -> ReedBoreVoice {
        let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap();
        ReedBoreVoice::with_inlet_line(
            || Ok(0.007), &gas, reed,
            PlateBank::default(), 1.0, 48000,
            |_| Ok((
                DelayedFilter::from_impulse_response(1.0/48000.0,
                    vec![direct, 0.0, 0.0, -0.1, 0.0, 0.0, 0.0, 0.02]).unwrap(),
                crate::driving_point::LineRealization {
                    n_fft: 8, tail_energy: 0.0, peak_shift_cents: 0.0,
                    pole_states: Vec::new(), rungs_tried: Vec::new(),
                },
            )),
        ).unwrap()
    }

    // G0: derive reflection independently from accepted outgoing waves. The
    // old implementation's stale reflection fails the source flow equation,
    // even with r[0]=0, when the first delayed echo returns.
    #[test]
    fn source_motion_flow_and_committed_reflection_share_one_sample() {
        for direct in [-0.7, 0.0, 0.65] {
            for massive in [false, true] {
                let mut v = voice(direct, massive);
                let mut history = vec![5.0]; // the original explicit priming push
                for _ in 0..160 {
                    let (old_y, old_v) = (v.reed_y, v.reed_v);
                    v.step_block(&mut [0.0]).unwrap();
                    history.push(v.p_plus_prev);
                    let back = |k: usize| history.len().checked_sub(k+1)
                        .map_or(0.0, |i| history[i]);
                    let incoming = direct*back(0)-0.1*back(3)+0.02*back(7);
                    let actual = v.line.incoming();
                    assert!((incoming-actual).abs() < 1e-12*(1.0+incoming.abs()));
                    let pressure = v.p_plus_prev+actual;
                    let flow = (v.p_plus_prev-actual)/v.zc;
                    let dp = v.reed.blowing_pressure_pa-pressure;
                    let (opening, swept) = if massive {
                        (f64::midpoint(old_y,v.reed_y).max(0.0),
                            crate::reed_bore::reed_swept_flow(
                                crate::reed_bore::reed_pressure_face(v.reed),
                                f64::midpoint(old_v,v.reed_v)))
                    } else {
                        (crate::reed_bore::aperture_of(v.reed).opening_m(dp).max(0.0),0.0)
                    };
                    let jet = crate::bernoulli_aperture::moving::volume_flow(
                        v.reed.width_m,opening,dp,v.rho);
                    let scale = jet.abs()+swept.abs()+flow.abs();
                    assert!((jet+swept-flow).abs() < 1e-10*scale+1e-18,
                        "r={direct} massive={massive} flow residual={:e}",jet+swept-flow);
                }
            }
        }
    }

    // The opt-in constriction noise (the brass lip row's law on the reed channel) is admitted only before the first sample,
    // is seeded, and changes the voice only once the mouthpiece flow has been forward; with no source the voice is unchanged.
    #[test]
    fn constriction_noise_is_admitted_before_the_first_sample_only_and_is_seeded() {
        use crate::turbulence_noise::ConstrictionNoiseCard;
        let card = |seed| ConstrictionNoiseCard { seed, duct_radius_m: 0.007 };
        let run = |source: Option<ConstrictionNoiseCard>| {
            let mut v = voice(-0.2, true);
            if let Some(c) = source { v.set_constriction_noise(c).unwrap(); }
            let mut out = vec![0.0; 800];
            v.step_block(&mut out).unwrap();
            assert!(v.set_constriction_noise(card(1)).is_err(), "an advanced voice refuses the source");
            out
        };
        let (plain, a, a2, b) = (run(None), run(Some(card(5))), run(Some(card(5))), run(Some(card(6))));
        assert_eq!(a, a2, "one seed, one stream");
        assert_ne!(a, b, "another seed, another stream");
        assert_ne!(a, plain, "the source reaches the sample");
        assert_eq!(plain[0].to_bits(), a[0].to_bits(), "silent until the last sample's flow was forward");
        assert!(ConstrictionNoiseCard { seed: 1, duct_radius_m: 0.0 }.validate().is_err());
        assert!(voice(-0.2, true).set_constriction_noise(ConstrictionNoiseCard { seed: 1, duct_radius_m: -1.0 }).is_err());
    }

    // ---- staple pressure recovery (Almeida, Vergez and Caussé 2007, eqs. 14-17) -----------------------------------
    // Opt-in; absent, or with C_P = 0, the voice is bit for bit the original; refused outside its domain, by name.

    /// The oboe description's reed: 357 um rest gap on 5.6 mm (2.0 mm^2 at rest), mass and stiffness as its card.
    fn oboe_reed() -> BeatingReed {
        BeatingReed {
            rest_opening_m: 3.57e-4, width_m: 0.0056, closing_pressure_pa: 17720.0,
            blowing_pressure_pa: 7400.0, attack_s: 0.008, mass_kg: 2.2e-5,
            stiffness_n_m: 3464.0, damping_ratio: 0.114,
        }
    }

    fn staple_voice(recovery: Option<(f64, f64)>) -> ReedBoreVoice {
        use crate::bernoulli_aperture::viscous::StapleRecovery;
        let mut v = voice_with(-0.2, oboe_reed());
        v.set_viscous_slit_length(2.5e-2).unwrap(); // the cards' slit_length
        if let Some((c_p, s_c)) = recovery {
            v.set_staple_recovery(StapleRecovery::new(c_p, s_c).unwrap()).unwrap();
        }
        v
    }

    fn staple_run(recovery: Option<(f64, f64)>, samples: usize) -> Vec<f64> {
        let mut v = staple_voice(recovery);
        let mut out = vec![0.0; samples];
        v.step_block(&mut out).unwrap();
        out
    }

    #[test]
    fn g0_zero_or_absent_staple_recovery_is_the_original_voice_bit_for_bit() {
        let absent = staple_run(None, 1200);
        let zero = staple_run(Some((0.0, 4e-6)), 1200);
        assert!(absent.iter().all(|x| x.is_finite()));
        assert_eq!(
            absent.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
            zero.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
        );
        assert!(staple_voice(Some((0.0, 4e-6))).staple_recovery().is_none(), "C_P = 0 installs nothing");
        assert_eq!(staple_voice(Some((0.8, 4e-6))).staple_recovery().map(|r| r.c_p()), Some(0.8));
    }

    #[test]
    fn g1_staple_recovery_reaches_the_sample_and_keeps_the_callback_partition() {
        let plain = staple_run(None, 2400);
        let staple = staple_run(Some((0.8, 4e-6)), 2400);
        assert!(staple.iter().all(|x| x.is_finite()));
        let scale = plain.iter().fold(0.0_f64, |m, x| m.max(x.abs()));
        let difference = plain.iter().zip(&staple).fold(0.0_f64, |m, (a, b)| m.max((a - b).abs()));
        assert!(difference > 1e-3 * scale, "the recovery must change the voice: {difference:e} of {scale:e}");
        let (mut a, mut b) = (staple_voice(Some((0.8, 4e-6))), staple_voice(Some((0.8, 4e-6))));
        let (mut x, mut y) = ([0.0; 600], [0.0; 600]);
        a.step_block(&mut x).unwrap();
        for chunk in y.chunks_mut(37) { b.step_block(chunk).unwrap(); }
        assert_eq!(x, y);
        assert_eq!((a.reed_y, a.reed_v), (b.reed_y, b.reed_v));
    }

    #[test]
    fn g0_staple_recovery_is_refused_for_a_massless_reed_a_fresh_voice_or_an_undefined_rest_opening() {
        use crate::bernoulli_aperture::viscous::{StapleRecovery, STAPLE_RECOVERY_DOMAIN};
        let recovery = StapleRecovery::new(0.8, 4e-6).unwrap();
        // With or without the viscous slit: the row works on both laws.
        for slit in [None, Some(2.5e-2)] {
            let prepared = |mut v: ReedBoreVoice| {
                if let Some(length) = slit { v.set_viscous_slit_length(length).unwrap(); }
                v
            };
            assert!(prepared(voice(-0.2, false)).set_staple_recovery(recovery).is_err(), "massless reed {slit:?}");
            let mut advanced = prepared(voice_with(-0.2, oboe_reed()));
            advanced.step_block(&mut [0.0]).unwrap();
            assert!(advanced.set_staple_recovery(recovery).is_err(), "advanced voice {slit:?}");
            // Rest area 0.4 mm x 13 mm = 5.2 mm^2 is above S_c / sqrt(C_P) = 4.47 mm^2: refused, by name.
            let mut wide = prepared(voice(-0.2, true));
            match wide.set_staple_recovery(recovery) {
                Err(AcousticRealizeError::Reed { what }) => assert_eq!(what, STAPLE_RECOVERY_DOMAIN),
                other => panic!("expected the staple domain refusal, got {other:?}"),
            }
            assert!(wide.staple_recovery().is_none());
            assert!(prepared(voice_with(-0.2, oboe_reed())).set_staple_recovery(recovery).is_ok(), "{slit:?}");
        }
    }

    // The same law on the inviscid Bernoulli junction, the one the oboe cards use (no `slit_length`).
    fn inviscid_voice(recovery: Option<(f64, f64)>) -> ReedBoreVoice {
        use crate::bernoulli_aperture::viscous::StapleRecovery;
        let mut v = voice_with(-0.2, oboe_reed());
        if let Some((c_p, s_c)) = recovery {
            v.set_staple_recovery(StapleRecovery::new(c_p, s_c).unwrap()).unwrap();
        }
        v
    }

    fn inviscid_run(recovery: Option<(f64, f64)>, samples: usize) -> Vec<f64> {
        let mut v = inviscid_voice(recovery);
        let mut out = vec![0.0; samples];
        v.step_block(&mut out).unwrap();
        out
    }

    #[test]
    fn g0_zero_or_absent_staple_recovery_is_the_original_inviscid_voice_bit_for_bit() {
        let absent = inviscid_run(None, 1200);
        let zero = inviscid_run(Some((0.0, 4e-6)), 1200);
        assert!(absent.iter().all(|x| x.is_finite()));
        assert_eq!(
            absent.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
            zero.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
        );
        assert!(inviscid_voice(Some((0.0, 4e-6))).staple_recovery().is_none());
    }

    #[test]
    fn g1_inviscid_staple_recovery_reaches_the_sample_and_tracks_the_viscous_law_with_a_negligible_slit() {
        let plain = inviscid_run(None, 2400);
        let staple = inviscid_run(Some((0.8, 4e-6)), 2400);
        assert!(staple.iter().all(|x| x.is_finite()));
        let peak = plain.iter().fold(0.0_f64, |m, x| m.max(x.abs()));
        let difference = plain.iter().zip(&staple).fold(0.0_f64, |m, (a, b)| m.max((a - b).abs()));
        assert!(difference > 1e-3 * peak, "the recovery must change the inviscid voice: {difference:e} of {peak:e}");
        // The callback partition cannot matter.
        let (mut a, mut b) = (inviscid_voice(Some((0.8, 4e-6))), inviscid_voice(Some((0.8, 4e-6))));
        let (mut x, mut y) = ([0.0; 600], [0.0; 600]);
        a.step_block(&mut x).unwrap();
        for chunk in y.chunks_mut(37) { b.step_block(chunk).unwrap(); }
        assert_eq!(x, y);
        // With a slit so short that its viscous term is negligible, the viscous law gives the same voice: the two junctions
        // solve one equation (the junction-level test in `bernoulli_aperture::moving` pins the agreement to 1e-9).
        let mut slit = inviscid_voice(None);
        slit.set_viscous_slit_length_fully_developed(1e-9).unwrap(); // the jet law below L / h = 3 is alpha = 0.45, not the inviscid law
        slit.set_staple_recovery(crate::bernoulli_aperture::viscous::StapleRecovery::new(0.8, 4e-6).unwrap()).unwrap();
        let mut viscous = vec![0.0; 600];
        slit.step_block(&mut viscous).unwrap();
        let drift = viscous.iter().zip(&staple[..600]).fold(0.0_f64, |m, (a, b)| m.max((a - b).abs()));
        assert!(drift <= 1e-5 * peak, "negligible slit {drift:e} of {peak:e}");
    }

    // ---- jet contraction (the vena contracta coefficient C_c on the inviscid reed law) -------------------------------
    // Opt-in; absent, or with C_c = 1, the voice is bit for bit the original; the jet leaves with the area C_c w h.

    fn contracted_voice(contraction: f64, staple: Option<(f64, f64)>) -> ReedBoreVoice {
        use crate::bernoulli_aperture::viscous::StapleRecovery;
        let mut reed = oboe_reed();
        reed.attack_s = 0.0;
        let mut v = voice_with(-0.2, reed);
        if let Some((c_p, s_c)) = staple {
            v.set_staple_recovery(StapleRecovery::new(c_p, s_c).unwrap()).unwrap();
        }
        v.set_jet_contraction(contraction).unwrap();
        v
    }

    #[test]
    fn g0_unit_or_absent_jet_contraction_is_the_original_voice_bit_for_bit() {
        let run = |contraction: Option<f64>, staple: Option<(f64, f64)>| {
            let mut v = inviscid_voice(staple);
            if let Some(c) = contraction { v.set_jet_contraction(c).unwrap(); }
            let mut out = vec![0.0; 1200];
            v.step_block(&mut out).unwrap();
            assert!(out.iter().all(|x| x.is_finite()));
            out.iter().map(|x| x.to_bits()).collect::<Vec<_>>()
        };
        for staple in [None, Some((0.8, 4e-6))] {
            assert_eq!(run(None, staple), run(Some(1.0), staple), "staple {staple:?}");
        }
        let mut v = inviscid_voice(None);
        v.set_jet_contraction(1.0).unwrap();
        assert_eq!(v.jet_contraction(), None, "C_c = 1 installs nothing");
        v.set_jet_contraction(0.611).unwrap();
        assert_eq!(v.jet_contraction(), Some(0.611));
    }

    #[test]
    fn g1_jet_contraction_closes_the_same_sample_flow_on_the_contracted_area_and_reaches_the_sample() {
        for contraction in [0.611, 0.85] {
            let mut v = contracted_voice(contraction, None);
            for _ in 0..400 {
                let (old_y, old_v) = (v.reed_y, v.reed_v);
                v.step_block(&mut [0.0]).unwrap();
                let incoming = v.line.incoming();
                let pressure = v.p_plus_prev + incoming;
                let flow = (v.p_plus_prev - incoming) / v.zc;
                let dp = v.reed.blowing_pressure_pa - pressure;
                let opening = f64::midpoint(old_y, v.reed_y).max(0.0);
                let swept = crate::reed_bore::reed_swept_flow(
                    crate::reed_bore::reed_pressure_face(v.reed), f64::midpoint(old_v, v.reed_v));
                // The jet passes the contracted area; the blade's swept flow is geometric.
                let jet = contraction * crate::bernoulli_aperture::moving::volume_flow(v.reed.width_m, opening, dp, v.rho);
                let scale = jet.abs() + swept.abs() + flow.abs();
                assert!((jet + swept - flow).abs() < 1e-10 * scale + 1e-18,
                    "C_c={contraction}: flow residual {:e}", jet + swept - flow);
            }
        }
        // It reaches the sample, with or without the staple, and the callback partition cannot matter.
        for staple in [None, Some((0.8, 4e-6))] {
            let mut plain = contracted_voice(1.0, staple);
            let mut narrow = contracted_voice(0.611, staple);
            let (mut a, mut b) = (vec![0.0; 2400], vec![0.0; 2400]);
            plain.step_block(&mut a).unwrap();
            narrow.step_block(&mut b).unwrap();
            assert!(b.iter().all(|x| x.is_finite()));
            let peak = a.iter().fold(0.0_f64, |m, x| m.max(x.abs()));
            let difference = a.iter().zip(&b).fold(0.0_f64, |m, (x, y)| m.max((x - y).abs()));
            assert!(difference > 1e-3 * peak, "staple {staple:?}: {difference:e} of {peak:e}");
            let (mut c, mut d) = (contracted_voice(0.611, staple), contracted_voice(0.611, staple));
            let (mut x, mut y) = ([0.0; 600], [0.0; 600]);
            c.step_block(&mut x).unwrap();
            for chunk in y.chunks_mut(37) { d.step_block(chunk).unwrap(); }
            assert_eq!(x, y);
            assert_eq!((c.reed_y, c.reed_v), (d.reed_y, d.reed_v));
        }
    }

    #[test]
    fn g0_jet_contraction_is_refused_outside_its_domain_and_moves_the_staple_domain_to_the_jet_area() {
        use crate::bernoulli_aperture::viscous::{StapleRecovery, STAPLE_RECOVERY_DOMAIN};
        for bad in [0.0, -0.5, 1.0 + 1e-9, 2.0, f64::NAN, f64::INFINITY] {
            assert!(voice_with(-0.2, oboe_reed()).set_jet_contraction(bad).is_err(), "C_c={bad}");
        }
        assert!(voice(-0.2, false).set_jet_contraction(0.6).is_err(), "massless reed");
        let mut advanced = voice_with(-0.2, oboe_reed());
        advanced.step_block(&mut [0.0]).unwrap();
        assert!(advanced.set_jet_contraction(0.6).is_err(), "advanced voice");
        // Defined on the inviscid law only, in either order with the viscous channel.
        let mut slit = voice_with(-0.2, oboe_reed());
        slit.set_viscous_slit_length(2.5e-2).unwrap();
        assert!(slit.set_jet_contraction(0.6).is_err());
        let mut late = voice_with(-0.2, oboe_reed());
        late.set_jet_contraction(0.6).unwrap();
        assert!(late.set_viscous_slit_length(2.5e-2).is_err());
        // Rest area 0.4 mm x 13 mm = 5.2 mm^2: the staple law refuses it at C_c = 1; the jet area 0.8 x 5.2 = 4.16 mm^2 is
        // inside S_c / sqrt(C_P) = 4.47 mm^2, and 0.9 x 5.2 = 4.68 mm^2 is not.
        let recovery = StapleRecovery::new(0.8, 4e-6).unwrap();
        let mut ok = voice(-0.2, true);
        ok.set_jet_contraction(0.8).unwrap();
        assert!(ok.set_staple_recovery(recovery).is_ok());
        let mut wide = voice(-0.2, true);
        wide.set_jet_contraction(0.9).unwrap();
        match wide.set_staple_recovery(recovery) {
            Err(AcousticRealizeError::Reed { what }) => assert_eq!(what, STAPLE_RECOVERY_DOMAIN),
            other => panic!("expected the staple domain refusal, got {other:?}"),
        }
    }

    // G5: the physical source, its history and receiver are independent of the
    // host callback partition for both existing source mechanisms.
    #[test]
    fn same_sample_reed_preserves_callback_partition_and_control_history() {
        for massive in [false,true] {
            let (mut a,mut b) = (voice(-0.4,massive),voice(-0.4,massive));
            for pressure in [100.0,250.0,0.0] {
                a.set_blowing_pressure(pressure); b.set_blowing_pressure(pressure);
                let (mut x,mut y) = ([0.0;127],[0.0;127]);
                a.step_block(&mut x).unwrap();
                for chunk in y.chunks_mut(19) { b.step_block(chunk).unwrap(); }
                assert_eq!(x,y);
                assert_eq!(a.line.history(),b.line.history());
                assert_eq!((a.reed_y,a.reed_v),(b.reed_y,b.reed_v));
            }
        }
    }

    // G1: independently reconstruct the original affine wave boundary and the
    // parallel-plate plus separation pressure law. A post-filtered inviscid
    // result, an old-sample gap or an omitted swept flow cannot satisfy it.
    /// The card's slit selector (the separated-jet law), or the fully developed law it selected before.
    fn select_slit(v: &mut ReedBoreVoice, length: f64, jet_law: bool) {
        if jet_law { v.set_viscous_slit_length(length).unwrap(); } else { v.set_viscous_slit_length_fully_developed(length).unwrap(); }
        assert_eq!(v.viscous_slit().unwrap().is_separated_jet(), jet_law);
    }

    #[test]
    fn viscous_source_closes_the_same_sample_flow_for_both_reed_mechanisms() { source_closes_the_same_sample_flow(false); }

    #[test]
    fn jet_source_closes_the_same_sample_flow_for_both_reed_mechanisms() { source_closes_the_same_sample_flow(true); }

    fn source_closes_the_same_sample_flow(jet_law: bool) {
        for direct in [-0.7, 0.0, 0.65] {
            for massive in [false, true] {
                let mut v = voice(direct, massive);
                select_slit(&mut v, 0.001, jet_law);
                let mu = v.viscous_slit().unwrap().dynamic_viscosity_pa_s();
                for _ in 0..160 {
                    let (old_y, old_v) = (v.reed_y, v.reed_v);
                    v.step_block(&mut [0.0]).unwrap();
                    let incoming = v.line.incoming();
                    let pressure = v.p_plus_prev + incoming;
                    let total = (v.p_plus_prev - incoming) / v.zc;
                    let dp = v.reed.blowing_pressure_pa - pressure;
                    let (opening, swept) = if massive {
                        (f64::midpoint(old_y, v.reed_y).max(0.0),
                            crate::reed_bore::reed_swept_flow(
                                crate::reed_bore::reed_pressure_face(v.reed),
                                f64::midpoint(old_v, v.reed_v)))
                    } else {
                        (crate::reed_bore::aperture_of(v.reed).opening_m(dp), 0.0)
                    };
                    assert!(opening > 0.0);
                    let velocity = (total - swept) / (v.reed.width_m * opening);
                    if jet_law {
                        // The same-sample flow is the slit's own flow at the same-sample drop (`evaluate` is the independent entry point).
                        let flow = v.viscous_slit().unwrap().evaluate(v.reed.width_m, opening, v.rho, dp).unwrap().flow_m3_s;
                        let expected = total - swept;
                        assert!((flow - expected).abs() <= 1e-9 * (flow.abs() + expected.abs()) + 1e-18,
                            "r={direct} massive={massive}: {flow} != {expected}");
                        continue;
                    }
                    let viscous = 12.0 * mu * 0.001 * velocity / opening.powi(2);
                    let separation = 0.5 * v.rho * velocity * velocity.abs();
                    assert!((viscous + separation - dp).abs()
                        < 1e-9 * (dp.abs() + viscous.abs() + separation.abs()) + 1e-12,
                        "r={direct} massive={massive}: {viscous} + {separation} != {dp}");
                    assert!(viscous * velocity >= 0.0);
                }
                if !massive {
                    assert!(v.junction_residual_max().unwrap() < 1e-16);
                }
            }
        }
    }

    // G1/G3: actual beating and reopening must keep the contact potential,
    // mechanical damping and jet/wall heat in one pressure-work balance.
    #[test]
    fn viscous_massive_reed_closes_and_reopens_without_losing_pressure_work() { massive_reed_closes_and_reopens(false); }

    #[test]
    fn jet_massive_reed_closes_and_reopens_without_losing_pressure_work() { massive_reed_closes_and_reopens(true); }

    fn massive_reed_closes_and_reopens(jet_law: bool) {
        use crate::unilateral_contact::distributed::ApertureContactStep;
        let mut v = voice(-0.2, true);
        select_slit(&mut v, 0.001, jet_law);
        let slit = *v.viscous_slit().unwrap();
        let (stiffness, damping) = crate::reed_bore::reed_structural(v.reed);
        let face = crate::reed_bore::reed_pressure_face(v.reed);
        for _ in 0..10 {
            let mut closed = false;
            let mut reopened = false;
            for pressure in [15000.0, 100.0] {
                v.set_blowing_pressure(pressure);
                for _ in 0..192 {
                    let (old_y, old_v) = (v.reed_y, v.reed_v);
                    v.step_block(&mut [0.0]).unwrap();
                    // A step the closure stop engaged ends shut (opening exactly 0, velocity 0) at the mean velocity that
                    // brings the reed to zero in the step; every other step is the free midpoint scheme.
                    let stopped = v.reed_y == 0.0;
                    let velocity = if stopped { (v.reed_y - old_y) / v.dt } else { f64::midpoint(old_v, v.reed_v) };
                    let gap = f64::midpoint(old_y, v.reed_y).max(0.0);
                    let dp = pressure - (v.p_plus_prev + v.line.incoming());
                    let flow = (v.p_plus_prev - v.line.incoming()) / v.zc;
                    let jet = slit.evaluate(v.reed.width_m, gap, v.rho, dp).unwrap();
                    let lay = v.lay.as_ref().unwrap();
                    let contact = ApertureContactStep::new(lay, old_y).unwrap()
                        .response(v.reed_y, velocity).unwrap();
                    let energy = |y: f64, speed: f64| {
                        0.5 * v.reed.mass_kg * speed * speed
                            + 0.5 * stiffness * (y - v.reed.rest_opening_m).powi(2)
                            + lay.stiffness() / (lay.alpha() + 1.0)
                                * (-lay.gaps()[0] - y).max(0.0).powf(lay.alpha() + 1.0)
                    };
                    let delta = energy(v.reed_y, v.reed_v) - energy(old_y, old_v);
                    let loss = v.dt * (damping * velocity * velocity
                        + contact.dissipated_power + jet.dissipated_power_w());
                    let work = v.dt * dp * flow;
                    let scale = delta.abs() + loss.abs() + work.abs();
                    let residual = delta + loss - work;
                    if stopped {
                        // The stop is rigid and inelastic: it removes energy at the impact, never adds any.
                        assert!(residual <= 2e-9 * scale + 2e-18,
                            "closure stop created energy {:e} at gap {:e}", residual, gap);
                    } else {
                        assert!(residual.abs() <= 2e-9 * scale + 2e-18,
                            "work residual {:e} at gap {:e}", residual, gap);
                    }
                    let swept = -face * velocity;
                    assert!((jet.flow_m3_s + swept - flow).abs()
                        < 1e-10 * (jet.flow_m3_s.abs() + swept.abs() + flow.abs()) + 1e-18);
                    if pressure > 100.0 { closed |= v.reed_y <= 0.0; }
                    else { reopened |= v.reed_y > 0.0; }
                }
            }
            assert!(closed && reopened, "fixture must exercise both physical branches");
        }
    }

    #[test]
    fn viscous_selection_preserves_blocks_and_declares_strict_fallback() {
        for massive in [false, true] {
            let (mut a, mut b) = (voice(-0.4, massive), voice(-0.4, massive));
            a.set_viscous_slit_length(0.001).unwrap();
            b.set_viscous_slit_length(0.001).unwrap();
            b.set_solver_mode(ReedSolverMode::FastNewton);
            assert_eq!(b.realized_resolution().step_tolerance, None);
            assert!(b.set_junction_step_tolerance(1e-8).is_err());
            for pressure in [100.0, 250.0, 0.0] {
                a.set_blowing_pressure(pressure);
                b.set_blowing_pressure(pressure);
                let (mut x, mut y) = ([0.0; 127], [0.0; 127]);
                a.step_block(&mut x).unwrap();
                for chunk in y.chunks_mut(19) { b.step_block(chunk).unwrap(); }
                assert_eq!(x, y);
                assert_eq!(a.line.history(), b.line.history());
                assert_eq!((a.reed_y, a.reed_v), (b.reed_y, b.reed_v));
            }
            assert!(b.fast_stats.fallback_samples > 0);
            let original = *a.viscous_slit().unwrap();
            assert!(original.is_separated_jet(), "the card's slit_length selects the separated-jet law");
            assert!(a.set_viscous_slit_length(0.002).is_err());
            assert_eq!(a.viscous_slit(), Some(&original));
        }
        for length in [0.0, -0.001, f64::NAN, f64::INFINITY] {
            let mut v = voice(0.0, false);
            assert!(v.set_viscous_slit_length(length).is_err());
            assert!(v.viscous_slit().is_none());
        }
    }

    #[test]
    fn singular_massive_load_is_not_replaced_with_a_delayed_or_clipped_boundary() {
        for direct in [-1.0,1.0,f64::NAN] {
            let mut v = voice(0.0,true);
            let history = v.line.history();
            let state = (v.reed_y,v.reed_v,v.sample_index);
            assert!(v.step_massive_scattering(ScatteringResponse {
                direct,history:0.0,
            },100.0,0.0).is_err());
            assert_eq!(history,v.line.history());
            assert_eq!(state,(v.reed_y,v.reed_v,v.sample_index));
        }
    }
}
