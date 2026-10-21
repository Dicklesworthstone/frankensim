//! Read the affine output of the next DelayedFilter push without advancing it.
//! The existing filter remains the only state owner and sample commit path.

use super::*;
use crate::scattering::ScatteringResponse;

impl DelayedFilter {
    /// Current-sample affine output `incoming = direct*outgoing + history`.
    ///
    /// Unlike [`Self::incoming`], this observes the sample about to be accepted,
    /// not the previous output. An implicit physical junction solves against
    /// this relation before calling [`Self::push`] exactly once. The original
    /// FIR/IIR, explicit poles, compliance and dispersion memories are retained.
    /// Reading allocates nothing and changes no state. FIR reading is O(taps);
    /// the ordinary push still performs its original convolution/arithmetic.
    ///
    /// The two forms can differ by summation roundoff because the affine form
    /// collects the direct coefficient separately. This is not a passivity
    /// certificate, a response fit, or a change to the legacy push convention.
    ///
    /// # Errors
    /// Nonfinite pending history, direct coefficient, filter or dispersion output.
    pub fn boundary_response(&self) -> Result<ScatteringResponse, DiscretizeError> {
        let (mut direct, history) = if self.fir.is_empty() {
            // The bulk delay is at least two samples: its next read is already
            // known, even when the downstream IIR itself has feedthrough.
            let n = self.buf.len();
            let i1 = (self.write + n - self.delay_int) % n;
            let i0 = (i1 + n - 1) % n;
            let delayed = (1.0-self.frac)*self.buf[i1] + self.frac*self.buf[i0];
            (0.0, self.filter.validate_step(&self.state, delayed)?)
        } else {
            let mut history = 0.0;
            // write names the current unknown input. Tap k>=1 reads its
            // existing predecessor, in the same ascending tap order as push.
            // Reuse push's two contiguous runs instead of per-tap modulo.
            let split = (self.write+1).min(self.fir.len());
            let (head, tail) = self.fir.split_at(split);
            for (&tap, &value) in head[1..].iter().zip(self.buf[..self.write].iter().rev()) {
                history += tap*value;
            }
            for (&tap, &value) in tail.iter().zip(self.buf[self.write+1..].iter().rev()) {
                history += tap*value;
            }
            let dt = self.filter.t_s;
            let mut direct = self.fir[0];
            if self.compliance_gain > 0.0 {
                history += self.compliance_state;
                direct += 0.5*self.compliance_gain*dt;
            }
            for pole in &self.exact_poles {
                history += dt*pole.state.re;
                direct += dt*(0.5*pole.residue.re);
            }
            (direct, history)
        };
        require_finite_scalar("boundary_history", history)?;
        // The allpass recurrence is affine too. Reuse its original immutable
        // history evaluation; only its current-input derivative is new.
        let history = self.validate_dispersion(history)?;
        for _ in &self.disp_x1 { direct *= self.disp_a; }
        require_finite_scalar("boundary_direct", direct)?;
        Ok(ScatteringResponse { direct, history })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64) {
        assert!((a-b).abs() <= 4e-12*(1.0+a.abs()+b.abs()), "{a:e} != {b:e}");
    }

    fn exercise(mut line: DelayedFilter, samples: usize) {
        for i in 0..samples {
            let old = format!("{line:?}");
            let response = line.boundary_response().unwrap();
            assert_eq!(response, line.boundary_response().unwrap());
            assert_eq!(old, format!("{line:?}"));
            for outgoing in [-2.0,0.0,3.0] {
                // Probe the ORIGINAL committing owner on a test-only clone.
                let mut probe = line.clone();
                close(response.incoming(outgoing).unwrap(), probe.push(outgoing).unwrap());
            }
            let outgoing = ((17*i+3)%41) as f64/19.0-1.0;
            close(response.incoming(outgoing).unwrap(), line.push(outgoing).unwrap());
        }
    }

    #[test]
    fn g0_fir_preview_matches_original_convolution_through_every_ring_position() {
        for length in [4,7,19] {
            for direct in [-0.7,0.0,0.8] {
                let mut taps: Vec<_> = (0..length).map(|k| (k as f64+1.0)/100.0).collect();
                taps[0] = direct;
                exercise(DelayedFilter::from_impulse_response(0.001,taps).unwrap(),3*length);
            }
        }
    }

    #[test]
    fn g0_echo_is_visible_on_its_physical_sample_not_one_step_late() {
        let mut line = DelayedFilter::from_impulse_response(0.001,vec![0.0,0.0,-0.5,0.0]).unwrap();
        line.push(1.0).unwrap(); line.push(0.0).unwrap();
        assert_eq!(line.incoming(),0.0);
        assert_eq!(line.boundary_response().unwrap(),ScatteringResponse {direct:0.0,history:-0.5});
        assert_eq!(line.push(7.0).unwrap(),-0.5);
    }

    #[test]
    fn g0_preview_includes_compliance_poles_and_full_dispersion_history() {
        for sections in [0,1,3] {
            let line = DelayedFilter::from_impulse_response(1.0/48000.0,
                vec![0.4,0.0,-0.2,0.1,0.0]).unwrap()
                .with_exact_poles(12.0,&[
                    (1400.0,0.0,C64::new(0.7,-0.2)),
                    (5000.0,25.0,C64::new(-0.3,0.1)),
                ]).unwrap();
            let line = if sections == 0 {line} else {line.with_dispersion(-0.4,sections).unwrap()};
            exercise(line,160);
        }
    }

    #[test]
    fn g0_bulk_delay_iir_and_dispersion_have_no_current_input_term() {
        let filter = DigitalFilter {sections:vec![Biquad {
            b:[0.1,0.2,0.0],a:[-0.4,0.1],
        }],direct:0.2,t_s:0.001,prewarp:0.0};
        for delay in [2.0,3.25,7.5] {
            let line = DelayedFilter::new(delay,filter.clone()).unwrap()
                .with_dispersion(0.3,2).unwrap();
            assert_eq!(line.boundary_response().unwrap().direct,0.0);
            exercise(line,160);
        }
    }

    #[test]
    fn g0_nonfinite_pending_response_refuses_without_advancing() {
        let mut line = DelayedFilter::from_impulse_response(0.001,
            vec![0.0,f64::MAX,0.0,0.0]).unwrap();
        line.push(2.0).unwrap();
        let before = format!("{line:?}");
        assert!(line.boundary_response().is_err());
        assert_eq!(format!("{line:?}"),before);
    }
}
