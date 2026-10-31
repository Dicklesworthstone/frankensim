//! Conditional contraction of the retained nonemitting scattering model.
//!
//! This is an ideal-transport envelope, not a certificate for floating-point
//! free flights, phase sampling, geometry, input calibration or spectral error.

use super::{
    MassAbsorptionSpectrum, MassScatteringSpectrum, PhaseFunction, ScatteringGrid,
    VolumeTransportError, checkpoint,
};
use fs_exec::Cx;
use fs_ivl::Interval;
use std::collections::{HashMap, HashSet};

impl MassScatteringSpectrum {
    /// Enclose the retained tabulated, Rayleigh, or combined mass spectrum.
    /// This shares the original coefficient owner with the contraction bound;
    /// no fitting of a wavelength-dependent law into a gray coefficient occurs.
    pub(crate) fn coefficient_enclosure(
        &self,
        wavelength: f64,
    ) -> Result<Interval, VolumeTransportError> {
        let mut value = self.spectrum.coefficient_enclosure(wavelength)?;
        if let Some(reference) = self.rayleigh_reference_nm {
            if value.hi() != 0.0 {
                let ratio = Interval::point(reference) / Interval::point(wavelength);
                let square = ratio * ratio;
                value = value * square * square;
            }
        }
        if let Some((reference, coefficient)) = self.rayleigh_addend {
            if coefficient != 0.0 {
                let ratio = Interval::point(reference) / Interval::point(wavelength);
                let square = ratio * ratio;
                let extra = Interval::point(coefficient) * square * square;
                value = if value.hi() == 0.0 { extra } else { value + extra };
            }
        }
        if !value.lo().is_finite() || !value.hi().is_finite() {
            return Err(VolumeTransportError::Unrepresentable);
        }
        Ok(Interval::new(value.lo().max(0.0), value.hi()))
    }
}

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
        if self.absorption.absorption_palette().is_none() && self.scattering_palette().is_none() {
            let albedo = max_albedo_upper(&self.absorption.absorption, &self.scattering);
            if albedo >= 1.0 {
                return Ok(None);
            }
            let Some(phase_mass) = phase_integral_upper(&self.phase, cx)? else {
                return Ok(None);
            };
            let gain = (Interval::point(albedo) * Interval::point(phase_mass)).hi();
            checkpoint(cx)?;
            return Ok((gain.is_finite() && gain > 0.0 && gain < 1.0).then_some(gain));
        }

        // Independent absorption/scattering maps define physical PAIRS, not
        // every Cartesian combination of their palette entries. Empty cells
        // and unreferenced materials cannot contribute ideal transport gain.
        // Cache each occupied pair and each retained angular integral once;
        // insertion/iteration order never participates in the maximum.
        let mut pairs = HashSet::new();
        let mut phase_masses = HashMap::new();
        let mut maximum = 0.0_f64;
        for (cell, &density) in self.absorption.grid.data.iter().enumerate() {
            if cell % 1024 == 0 {
                checkpoint(cx)?;
            }
            if density == 0.0 {
                continue;
            }
            let key = self.coefficient_key(cell)?;
            if pairs.contains(&key) {
                continue;
            }
            pairs
                .try_reserve(1)
                .map_err(|_| VolumeTransportError::Budget)?;
            pairs.insert(key);
            let albedo = max_albedo_upper(
                self.absorption.absorption_at_cell(cell)?,
                self.scattering_at_cell(cell)?,
            );
            if albedo == 0.0 {
                continue;
            }
            if albedo >= 1.0 {
                return Ok(None);
            }
            let phase_mass = if let Some(&mass) = phase_masses.get(&key.1) {
                mass
            } else {
                let Some(mass) = phase_integral_upper(self.phase_at_cell(cell)?, cx)? else {
                    return Ok(None);
                };
                phase_masses
                    .try_reserve(1)
                    .map_err(|_| VolumeTransportError::Budget)?;
                phase_masses.insert(key.1, mass);
                mass
            };
            let gain = (Interval::point(albedo) * Interval::point(phase_mass)).hi();
            if !gain.is_finite() || gain >= 1.0 {
                return Ok(None);
            }
            maximum = maximum.max(gain);
        }
        checkpoint(cx)?;
        Ok(Some(maximum))
    }
}

fn max_albedo_upper(
    absorption: &MassAbsorptionSpectrum,
    scattering: &MassScatteringSpectrum,
) -> f64 {
    let mut maximum = 0.0_f64;
    // On each combined-knot interval, a=A+B*lambda and
    // s=C+D*lambda+K/lambda^4, with a,s,K nonnegative. The sign of the
    // albedo derivative is h=s'*a-s*a'; h'=s''*a=20*K*a/lambda^6 >= 0.
    // Every interior stationary point is therefore a minimum (or a flat
    // region), so the maximum is still attained at a combined-table knot.
    // This includes pure linear, pure Rayleigh, and their physical mixtures.
    for &(wavelength, _) in absorption.knots[..absorption.count]
        .iter()
        .chain(scattering.knots().iter())
    {
        let a = coefficient_enclosure(absorption, wavelength);
        let s = scattering.coefficient_enclosure(wavelength)
            .unwrap_or(Interval::new(0.0, f64::INFINITY));
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
    spectrum.coefficient_enclosure(wavelength)
        .unwrap_or(Interval::new(0.0, f64::INFINITY))
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
        PhaseFunction::RayleighDepolarized { depolarization } => {
            let d = Interval::point(*depolarization);
            let amplitude = Interval::point(1.5 * (1.0 / (4.0 * core::f64::consts::PI)));
            (amplitude
                * Interval::point(4.0)
                * pi_enclosure()
                * ((Interval::point(1.0) + d) + (Interval::point(1.0) - d) / Interval::point(3.0))
                / (Interval::point(2.0) + d))
                .hi()
        }
        PhaseFunction::Tabulated(table) => {
            return tabulated_phase_integral_upper(table.knots_per_sr(), cx)
                .map(|upper| upper.map(|value| value.max(1.0)));
        }
        PhaseFunction::SpectralTable(table) => {
            // Convex spectral interpolation cannot exceed the largest
            // endpoint phase integral. No wavelength sampling is used here.
            let mut maximum = 1.0_f64;
            for (_, phase) in table.knots() {
                let Some(bound) = phase_integral_upper(phase, cx)? else {
                    return Ok(None);
                };
                maximum = maximum.max(bound);
            }
            maximum
        }
        PhaseFunction::SpectralMixture(mixture) => {
            // At every wavelength the angular law is a convex combination
            // weighted by component scattering, not by mass alone. Its
            // integral cannot exceed the largest retained component bound.
            // The unit baseline also covers the isotropic zero-scattering
            // fallback where no component has a defined conditional weight.
            let mut maximum = 1.0_f64;
            for component in mixture.components() {
                checkpoint(cx)?;
                if component.mass_weight == 0.0 {
                    continue;
                }
                let Some(bound) = phase_integral_upper(&component.phase, cx)? else {
                    return Ok(None);
                };
                maximum = maximum.max(bound);
            }
            maximum
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
    use crate::volumes::transport::scattering::{
        MassScatteringSpectrum, ScatteringMixtureComponent, ScatteringOptics,
    };
    use crate::volumes::transport::{
        AbsorbingGrid, AbsorbingGridBudget, CellPalette, CellPaletteBudget, OpticalMixtureBudget,
    };
    use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
    use std::sync::Arc;

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
            MassScatteringSpectrum::table(&[(380.0, 1.0), (580.0, 5.0), (780.0, 1.0)]).unwrap();
        // At 480 nm, s=3 and a=1, so omega=3/4. At the scattering-only
        // knot 580 nm, a=2 and omega=5/7; endpoints have omega=1/5.
        encloses_tightly(max_albedo_upper(&absorption, &scattering), 0.75);
    }

    #[test]
    fn a_shared_zero_knot_does_not_discard_neighboring_albedo_limits() {
        let absorption =
            MassAbsorptionSpectrum::table(&[(380.0, 1.0), (580.0, 0.0), (780.0, 2.0)]).unwrap();
        let scattering =
            MassScatteringSpectrum::table(&[(380.0, 3.0), (580.0, 0.0), (780.0, 6.0)]).unwrap();
        encloses_tightly(max_albedo_upper(&absorption, &scattering), 0.75);
        let tiny = MassAbsorptionSpectrum::constant(f64::from_bits(1)).unwrap();
        let tiny_scattering = MassScatteringSpectrum::constant(f64::from_bits(1)).unwrap();
        assert!(max_albedo_upper(&tiny, &tiny_scattering) >= 0.5);
    }

    #[test]
    fn analytic_and_retained_tabulated_phase_mass_is_enclosed() {
        with_cx(|cx| {
            for phase in [
                PhaseFunction::Isotropic,
                PhaseFunction::Rayleigh,
                PhaseFunction::rayleigh_depolarized(0.03).unwrap(),
                PhaseFunction::spectral_table(&[
                    (380.0, PhaseFunction::Rayleigh),
                    (780.0, PhaseFunction::henyey_greenstein(0.8).unwrap()),
                ])
                .unwrap(),
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
    fn rayleigh_albedo_envelope_covers_nonmonotone_absorption_ratios() {
        for coefficients in [[5.0, 0.1], [0.1, 5.0], [1.0, 1.0]] {
            let absorption = MassAbsorptionSpectrum::table(&[
                (380.0, coefficients[0]),
                (780.0, coefficients[1]),
            ])
            .unwrap();
            let scattering = MassScatteringSpectrum::rayleigh(550.0, 1.0).unwrap();
            let bound = max_albedo_upper(&absorption, &scattering);
            for i in 0..=1000 {
                let wavelength = 380.0 + 0.4 * f64::from(i);
                let s = scattering.at(wavelength).unwrap();
                let a = absorption.at(wavelength).unwrap();
                assert!(bound >= s / (a + s));
            }
        }
    }

    #[test]
    fn linear_plus_rayleigh_albedo_retains_tight_endpoints_across_an_interior_minimum() {
        let absorption = MassAbsorptionSpectrum::constant(1.0).unwrap();
        let mut scattering =
            MassScatteringSpectrum::table(&[(380.0, 380.0 / 150.0), (780.0, 780.0 / 150.0)])
                .unwrap();
        scattering.rayleigh_addend = Some((600.0, 1.0));
        // Independent law s=lambda/150+(600/lambda)^4 has its minimum at
        // 600 nm. The mixed law must retain both terms: bounding the linear
        // table alone would miss the large short-wavelength Rayleigh tail.
        let albedo = |wavelength: f64| {
            let s = wavelength / 150.0 + (600.0 / wavelength).powi(4);
            s / (1.0 + s)
        };
        assert!(albedo(600.0) < albedo(380.0));
        assert!(albedo(600.0) < albedo(780.0));
        let upper = max_albedo_upper(&absorption, &scattering);
        encloses_tightly(upper, albedo(380.0).max(albedo(780.0)));
        for sample in 0..=4000 {
            let wavelength = 380.0 + f64::from(sample) / 10.0;
            assert!(upper >= albedo(wavelength));
        }

        // A separate absorption knot creates additional extrema that the
        // scattering table does not contain. The same union-knot bound must
        // still enclose every independently evaluated mixed-law sample.
        let absorption =
            MassAbsorptionSpectrum::table(&[(380.0, 4.0), (480.0, 0.5), (780.0, 2.0)]).unwrap();
        let upper = max_albedo_upper(&absorption, &scattering);
        for sample in 0..=4000 {
            let wavelength = 380.0 + f64::from(sample) / 10.0;
            let a = if wavelength <= 480.0 {
                4.0 - (wavelength - 380.0) * 0.035
            } else {
                0.5 + (wavelength - 480.0) * 0.005
            };
            let s = wavelength / 150.0 + (600.0 / wavelength).powi(4);
            assert!(upper >= s / (a + s));
        }
    }

    #[test]
    fn spectral_mixture_envelope_covers_wavelength_dependent_normalized_phase() {
        with_cx(|cx| {
            let components = [
                ScatteringMixtureComponent {
                    mass_weight: 1.0,
                    scattering: MassScatteringSpectrum::rayleigh(550.0, 1.0).unwrap(),
                    phase: PhaseFunction::Rayleigh,
                },
                ScatteringMixtureComponent {
                    mass_weight: 3.0,
                    scattering: MassScatteringSpectrum::table(&[(380.0, 0.1), (780.0, 2.0)])
                        .unwrap(),
                    phase: PhaseFunction::Isotropic,
                },
                ScatteringMixtureComponent {
                    mass_weight: 0.0,
                    scattering: MassScatteringSpectrum::constant(100.0).unwrap(),
                    phase: PhaseFunction::henyey_greenstein(0.9).unwrap(),
                },
            ];
            let phase = PhaseFunction::scattering_mixture(
                cx,
                &components,
                OpticalMixtureBudget { max_components: 3 },
            )
            .unwrap();
            let upper = phase_integral_upper(&phase, cx).unwrap().unwrap();
            encloses_tightly(upper, 1.0);
            let expected = components[..2]
                .iter()
                .map(|component| phase_integral_upper(&component.phase, cx).unwrap().unwrap())
                .fold(1.0_f64, f64::max);
            assert_eq!(upper, expected);
            for wavelength in [380.0, 550.0, 780.0] {
                // These components are quadratic in mu, so Simpson's angular
                // integral is exact in real arithmetic for every spectral mix.
                let integral = 2.0 * core::f64::consts::PI / 3.0
                    * (phase.eval_at(wavelength, -1.0).unwrap()
                        + 4.0 * phase.eval_at(wavelength, 0.0).unwrap()
                        + phase.eval_at(wavelength, 1.0).unwrap());
                assert!((integral - 1.0).abs() < 2e-14);
                assert!(integral <= upper + 4.0 * f64::EPSILON);
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

    #[test]
    fn indexed_gain_uses_occupied_pairs_and_ignores_conservative_vacuum_cells() {
        with_cx(|cx| {
            for density in [[1.0, 1.0, 0.0, 1.0, 1.0], [1.0; 5]] {
                let budget = CellPaletteBudget {
                    max_models: 3,
                    max_cells: 5,
                };
                let absorption = CellPalette::try_new(
                    cx,
                    [1.0, 100.0, 0.0]
                        .map(|a| MassAbsorptionSpectrum::constant(a).unwrap())
                        .to_vec(),
                    Arc::from([0, 1, 2, 0, 1]),
                    budget,
                )
                .unwrap();
                let scattering = CellPalette::try_new(
                    cx,
                    [100.0, 1.0, 10.0]
                        .map(|s| ScatteringOptics {
                            scattering: MassScatteringSpectrum::constant(s).unwrap(),
                            phase: PhaseFunction::Rayleigh,
                        })
                        .to_vec(),
                    Arc::from([1, 0, 2, 1, 0]),
                    budget,
                )
                .unwrap();
                let m = ScatteringGrid::try_new_indexed(
                    cx,
                    AbsorbingGrid::try_new_indexed(
                        cx,
                        VolumeGrid::new([5, 1, 1], &density, [0.0; 3], [1.0; 3]),
                        absorption,
                        None,
                        AbsorbingGridBudget {
                            max_grid_cells: 5,
                            max_steps_per_segment: 16,
                        },
                    )
                    .unwrap(),
                    scattering,
                )
                .unwrap();
                let gain = m.photometric_scattering_gain_upper(cx).unwrap();
                if density[2] == 0.0 {
                    // Occupied (a,s) pairs are (1,1) and (100,100). Mixing
                    // unrelated extrema would incorrectly produce 100/101;
                    // counting the conservative vacuum cell would refuse.
                    encloses_tightly(gain.unwrap(), 0.5);
                } else {
                    assert_eq!(gain, None);
                }
            }
        });
    }

    #[test]
    fn direct_extinction_encloses_the_retained_linear_plus_molecular_sum() {
        let mut spectrum = MassScatteringSpectrum::constant(0.125).unwrap();
        spectrum.rayleigh_addend = Some((400.0, 0.5));
        for (lambda, expected) in [(400.0, 0.625), (600.0, 145.0 / 648.0)] {
            let interval = spectrum.coefficient_enclosure(lambda).unwrap();
            assert!(interval.lo() <= expected && interval.hi() >= expected);
            assert!(interval.width() < 1e-13);
        }
        assert_eq!(spectrum.coefficient_enclosure(800.0), Err(VolumeTransportError::Wavelength));
    }

}
