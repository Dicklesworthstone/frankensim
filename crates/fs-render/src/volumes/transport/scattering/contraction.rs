//! Conditional contraction of the retained nonemitting scattering model.
//!
//! This is an ideal-transport envelope, not a certificate for floating-point
//! free flights, phase sampling, geometry, input calibration or spectral error.

use super::{
    MassAbsorptionSpectrum, PhaseFunction, ScatteringGrid, VolumeTransportError, checkpoint,
};
use fs_exec::Cx;
use fs_ivl::Interval;

impl ScatteringGrid<'_> {
    /// Maximum ideal scattering gain over the complete retained spectral laws.
    ///
    /// Density multiplies absorption and scattering equally, so it cancels in
    /// the single-scattering albedo. The returned gain includes an outward
    /// enclosure of the phase integral; normalized table knots are not assumed
    /// to integrate to exactly one after floating-point admission. A gain of
    /// one or more, a thermal source, or a nonrepresentable envelope supplies
    /// no contraction metadata. A physically empty scattering field has gain
    /// zero. The caller must separately admit all surface materials and sources.
    pub(crate) fn photometric_scattering_gain_upper(
        &self,
        cx: &Cx<'_>,
    ) -> Result<Option<f64>, VolumeTransportError> {
        checkpoint(cx)?;
        if self.has_thermal_emission() {
            return Ok(None);
        }
        if !self.has_scattering() {
            return Ok(Some(0.0));
        }
        let albedo = max_albedo_upper(&self.absorption.absorption, &self.scattering.0);
        if albedo >= 1.0 {
            return Ok(None);
        }
        let Some(phase_mass) = phase_integral_upper(self.phase_ref(), cx)? else {
            return Ok(None);
        };
        let gain = (Interval::point(albedo) * Interval::point(phase_mass)).hi();
        checkpoint(cx)?;
        Ok((gain.is_finite() && gain > 0.0 && gain < 1.0).then_some(gain))
    }
}

fn max_albedo_upper(
    absorption: &MassAbsorptionSpectrum,
    scattering: &MassAbsorptionSpectrum,
) -> f64 {
    let mut maximum = 0.0_f64;
    // Between successive knots in the combined tables, s/(a+s) is a ratio
    // of affine functions. Its derivative has constant sign wherever the
    // denominator is positive, so extrema occur at the combined knots.
    for &(wavelength, _) in absorption.knots[..absorption.count]
        .iter()
        .chain(scattering.knots[..scattering.count].iter())
    {
        let a = coefficient_enclosure(absorption, wavelength);
        let s = coefficient_enclosure(scattering, wavelength);
        if s.hi() == 0.0 {
            // Includes an exact zero-extinction knot. Adjacent positive
            // endpoints also cover the limiting albedo on either side.
            continue;
        }
        // The ratio increases with s and decreases with a. Using this
        // monotonicity avoids losing the shared numerator/denominator s.
        let numerator = Interval::point(s.hi());
        let denominator = Interval::point(a.lo()) + numerator;
        let ratio = numerator / denominator;
        // Passivity independently bounds every albedo by one, including
        // extreme coefficient ranges where interval division is uninformative.
        maximum = maximum.max(ratio.hi().min(1.0));
    }
    maximum
}

fn coefficient_enclosure(spectrum: &MassAbsorptionSpectrum, wavelength: f64) -> Interval {
    let knots = &spectrum.knots[..spectrum.count];
    let upper = knots.partition_point(|&(lambda, _)| lambda < wavelength);
    if upper < knots.len() && knots[upper].0 == wavelength {
        return Interval::point(knots[upper].1);
    }
    // Callers use wavelengths from admitted tables with the same endpoints.
    let (left, a) = knots[upper - 1];
    let (right, b) = knots[upper];
    if a == b {
        return Interval::point(a);
    }
    let fraction = (Interval::point(wavelength) - Interval::point(left))
        / (Interval::point(right) - Interval::point(left));
    let fraction = unit_interval(fraction);
    let complement = unit_interval(Interval::point(1.0) - fraction);
    let value = Interval::point(a) * complement + Interval::point(b) * fraction;
    // The represented linear law is a convex combination of its knots. This
    // exact range also removes spurious negative endpoints near underflow.
    Interval::new(value.lo().max(a.min(b)), value.hi().min(a.max(b)))
}

fn unit_interval(value: Interval) -> Interval {
    Interval::new(value.lo().max(0.0), value.hi().min(1.0))
}

fn pi_enclosure() -> Interval {
    Interval::new(
        fs_math::next_down(core::f64::consts::PI),
        fs_math::next_up(core::f64::consts::PI),
    )
}

fn phase_integral_upper(
    phase: &PhaseFunction,
    cx: &Cx<'_>,
) -> Result<Option<f64>, VolumeTransportError> {
    checkpoint(cx)?;
    phase.validate()?;
    let upper = match phase {
        PhaseFunction::Isotropic | PhaseFunction::HenyeyGreenstein { .. } => {
            // The ideal HG kernel has integral 4*pi for every admitted g.
            // Isotropic evaluation returns the native retained INV_FOUR_PI
            // coefficient exactly, without requiring access to a second copy.
            let amplitude = PhaseFunction::Isotropic.eval_cos(0.0)?;
            (Interval::point(4.0) * pi_enclosure() * Interval::point(amplitude)).hi()
        }
        PhaseFunction::Rayleigh => {
            // Integral of (1+mu^2) over solid angle is 16*pi/3. Evaluation at
            // zero retains the native rounded 0.75*INV_FOUR_PI amplitude.
            let amplitude = phase.eval_cos(0.0)?;
            (Interval::point(16.0) * pi_enclosure() * Interval::point(amplitude)
                / Interval::point(3.0))
            .hi()
        }
        PhaseFunction::Tabulated(table) => {
            return tabulated_phase_integral_upper(table.knots_per_sr(), cx)
                .map(|upper| upper.map(|value| value.max(1.0)));
        }
    };
    // Do not turn downward normalization roundoff into physical absorption.
    Ok((upper.is_finite() && upper > 0.0).then_some(upper.max(1.0)))
}

fn tabulated_phase_integral_upper(
    knots: &[(f64, f64)],
    cx: &Cx<'_>,
) -> Result<Option<f64>, VolumeTransportError> {
    let mut mass = Interval::point(0.0);
    for (index, pair) in knots.windows(2).enumerate() {
        if index % 256 == 0 {
            checkpoint(cx)?;
        }
        if pair[0].1 == 0.0 && pair[1].1 == 0.0 {
            continue;
        }
        let width = Interval::point(pair[1].0) - Interval::point(pair[0].0);
        // Half each endpoint before adding, avoiding an unnecessary overflow
        // when a narrow phase lobe has large, finite density coefficients.
        let mean = Interval::point(0.5) * Interval::point(pair[0].1)
            + Interval::point(0.5) * Interval::point(pair[1].1);
        mass = mass + width * mean;
    }
    let upper = (Interval::point(2.0) * pi_enclosure() * mass).hi();
    checkpoint(cx)?;
    Ok((upper.is_finite() && upper > 0.0).then_some(upper))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::volumes::VolumeGrid;
    use crate::volumes::transport::scattering::MassScatteringSpectrum;
    use crate::volumes::transport::{AbsorbingGrid, AbsorbingGridBudget};
    use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};

    fn with_cx(f: impl FnOnce(&Cx<'_>)) {
        let gate = CancelGate::new_clock_free();
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|arena| {
            f(&Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 0xC071,
                    kernel_id: 1,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            ));
        });
    }

    fn encloses_tightly(upper: f64, exact: f64) {
        assert!(upper >= exact, "{upper} underestimates {exact}");
        assert!(upper <= exact * (1.0 + 1.0e-13));
    }

    #[test]
    fn combined_knots_include_absorption_only_albedo_extrema() {
        let absorption =
            MassAbsorptionSpectrum::table(&[(380.0, 4.0), (480.0, 1.0), (780.0, 4.0)]).unwrap();
        let scattering =
            MassAbsorptionSpectrum::table(&[(380.0, 1.0), (580.0, 5.0), (780.0, 1.0)]).unwrap();
        // At 480 nm, s=3 and a=1, so omega=3/4. At the scattering-only
        // knot 580 nm, a=2 and omega=5/7; endpoints have omega=1/5.
        encloses_tightly(max_albedo_upper(&absorption, &scattering), 0.75);
    }

    #[test]
    fn a_shared_zero_knot_does_not_discard_neighboring_albedo_limits() {
        let absorption =
            MassAbsorptionSpectrum::table(&[(380.0, 1.0), (580.0, 0.0), (780.0, 2.0)]).unwrap();
        let scattering =
            MassAbsorptionSpectrum::table(&[(380.0, 3.0), (580.0, 0.0), (780.0, 6.0)]).unwrap();
        encloses_tightly(max_albedo_upper(&absorption, &scattering), 0.75);
        let tiny = MassAbsorptionSpectrum::constant(f64::from_bits(1)).unwrap();
        assert!(max_albedo_upper(&tiny, &tiny) >= 0.5);
    }

    #[test]
    fn analytic_and_retained_tabulated_phase_mass_is_enclosed() {
        with_cx(|cx| {
            for phase in [
                PhaseFunction::Isotropic,
                PhaseFunction::Rayleigh,
                PhaseFunction::henyey_greenstein(0.8).unwrap(),
                PhaseFunction::henyey_greenstein(-fs_math::next_down(1.0)).unwrap(),
                PhaseFunction::tabulated(&[
                    (-1.0, 0.0),
                    (-0.6, 0.1),
                    (0.2, 0.3),
                    (0.9, 8.0),
                    (1.0, 12.0),
                ])
                .unwrap(),
            ] {
                encloses_tightly(phase_integral_upper(&phase, cx).unwrap().unwrap(), 1.0);
            }
        });
    }

    #[test]
    fn phase_enclosure_integrates_retained_density_instead_of_assuming_unity() {
        with_cx(|cx| {
            // A deliberately unnormalized retained density p=1/4 sr^-1 has
            // integral pi. The integration helper must inspect coefficients.
            let upper = tabulated_phase_integral_upper(&[(-1.0, 0.25), (1.0, 0.25)], cx)
                .unwrap()
                .unwrap();
            encloses_tightly(upper, core::f64::consts::PI);
        });
    }

    #[test]
    fn native_nonemitting_gain_is_density_independent_and_refuses_conservative_or_lte_media() {
        with_cx(|cx| {
            let temperatures = [3000.0; 2];
            let mut positive_density_gain = None;
            for (density, absorption, scattering, thermal, expected) in [
                ([1.0, 1.0], 1.0, 3.0, false, Some(0.75)),
                ([0.0, 16.0], 1.0, 3.0, false, Some(0.75)),
                ([0.0, 0.0], 0.0, 3.0, false, Some(0.0)),
                ([1.0, 1.0], 1.0, 0.0, false, Some(0.0)),
                ([1.0, 1.0], 0.0, 3.0, false, None),
                ([1.0, 1.0], 1.0, 3.0, true, None),
                ([1.0, 1.0], 1.0, 0.0, true, None),
            ] {
                let medium = ScatteringGrid::try_new(
                    cx,
                    AbsorbingGrid::try_new(
                        cx,
                        VolumeGrid::new([1, 1, 2], &density, [0.0; 3], [1.0; 3]),
                        MassAbsorptionSpectrum::constant(absorption).unwrap(),
                        thermal.then_some(temperatures.as_slice()),
                        AbsorbingGridBudget {
                            max_grid_cells: 2,
                            max_steps_per_segment: 8,
                        },
                    )
                    .unwrap(),
                    MassScatteringSpectrum::constant(scattering).unwrap(),
                    PhaseFunction::Isotropic,
                )
                .unwrap();
                let gain = medium.photometric_scattering_gain_upper(cx).unwrap();
                match expected {
                    Some(0.0) => assert_eq!(gain.unwrap().to_bits(), 0),
                    Some(value) => {
                        encloses_tightly(gain.unwrap(), value);
                        if let Some(previous) = positive_density_gain {
                            assert_eq!(gain, Some(previous));
                        }
                        positive_density_gain = gain;
                    }
                    None => assert_eq!(gain, None),
                }
            }
        });
    }
}
