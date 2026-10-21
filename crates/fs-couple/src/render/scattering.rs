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
        let (outgoing, opening, velocity) = crate::reed_bore::step_massive_reed(
            self.reed, self.rho, impedance, incoming, mouth_pressure_pa,
            self.reed_y, self.reed_v, self.dt, body_flow_m3_s,
            self.lay.as_ref(), self.solver_mode, &mut self.fast_stats,
        )?;
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
        let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap();
        ReedBoreVoice::with_inlet_line(
            || Ok(0.007), &gas,
            BeatingReed {
                rest_opening_m: 0.0004, width_m: 0.013,
                closing_pressure_pa: 6000.0, blowing_pressure_pa: 100.0,
                attack_s: 0.0, mass_kg: if massive { 0.00001 } else { 0.0 },
                stiffness_n_m: if massive { 100.0 } else { 0.0 }, damping_ratio: 0.35,
            },
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
