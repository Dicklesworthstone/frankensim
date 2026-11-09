//! Exact rectangular observation windows on the original infinite planes.
//!
//! At fixed depth every original component has real tangential wavevector q.
//! The window average of exp(-i (q_i-q_j).x) is the centre phase times
//! sinc((q_ix-q_jx) width/2) sinc((q_iy-q_jy) height/2). Integrating these
//! cross terms retains interference, including standing and evanescent waves.
//! This is a passive observation footprint, not an aperture boundary that
//! clips, diffracts, or changes the incident field.

use super::{
    AcousticSpectrumError, FieldSum, MAX_ANGULAR_SPECTRUM_COMPONENTS, PlaneField,
    PreparedAcousticSpectrum, PreparedLayerAcousticSpectrum, finite, poll, propagation,
};
use fs_math::{c64::C64, dd::Dd, det};

/// Maximum distinct component pairs, including each diagonal once.
pub const MAX_ACOUSTIC_WINDOW_PAIRS: usize =
    MAX_ANGULAR_SPECTRUM_COMPONENTS * (MAX_ANGULAR_SPECTRUM_COMPONENTS + 1) / 2;

/// Explicit allowance for the exact O(N²) coherent window integration.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct AcousticWindowBudget {
    /// At least N(N+1)/2 and at most [`MAX_ACOUSTIC_WINDOW_PAIRS`].
    /// Zero-pressure and repeated-direction components still count.
    pub max_pairs: usize,
}

/// Uniform rectangular sensor footprint parallel to the original interfaces.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RectangularAcousticWindow {
    /// Centre [m] in the corresponding point observation's coordinate frame.
    /// The depth is fixed across the window, including at layer interfaces.
    pub center_m: [f64; 3],
    /// Positive finite x/y side lengths [m]; no discretization is performed.
    pub extent_m: [f64; 2],
}

/// Spatially averaged phasors and separately integrated quadratic observables.
/// All amplitudes use the original complex PEAK-pressure convention.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AcousticWindowSample {
    /// Coherent spatial mean, as for an ideal uniform pressure-sensing face.
    pub mean_pressure_pa: C64,
    pub mean_particle_velocity_m_s: [C64; 3],
    /// Spatial mean of the time-averaged squared real pressure: <|p|²>/2.
    /// This generally differs from |<p>|²/2.
    pub mean_squared_pressure_pa2: f64,
    /// Spatial/time means <|v_j|²>/2 for the three velocity components.
    pub mean_squared_particle_velocity_m2_s2: [f64; 3],
    /// Spatial mean Re(<p conj(v)>)/2 [W/m²], not a product of mean fields.
    pub active_intensity_w_m2: [f64; 3],
    /// Spatial mean Im(<p conj(v)>)/2 [W/m²]; this is not absorbed heat.
    pub reactive_intensity_w_m2: [f64; 3],
}

/// The same window with separate local constitutive heat mechanisms averaged
/// over its area [W/m³]. It is not a volume-integrated heat or temperature.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AcousticLayerWindowSample {
    pub field: AcousticWindowSample,
    pub viscous_heating_w_m3: f64,
    /// Bulk-compliance loss, including pore diffusion and any declared
    /// molecular relaxation; not a gas/frame temperature partition.
    pub thermal_heating_w_m3: f64,
}

impl PreparedAcousticSpectrum<'_, '_> {
    /// Average the TOTAL incident plus reflected field over a finite sensor
    /// at z<=0, with the entrance-face origin used by [`Self::incident_field`].
    /// The exact integration uses O(N) scratch storage and budgeted O(N²)
    /// pair work. No material query, resampling or new boundary solve occurs.
    pub fn incident_window(
        &self,
        window: RectangularAcousticWindow,
        budget: AcousticWindowBudget,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<AcousticWindowSample, AcousticSpectrumError> {
        self.exterior_window(window, budget, false, cancelled)
    }

    /// Average the refracted field at z>=0 measured from the EXIT face.
    /// Critical and evanescent waves retain their original pressure, velocity,
    /// decay and reactive flow. This does not introduce finite-edge scattering.
    pub fn transmitted_window(
        &self,
        window: RectangularAcousticWindow,
        budget: AcousticWindowBudget,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<AcousticWindowSample, AcousticSpectrumError> {
        self.exterior_window(window, budget, true, cancelled)
    }

    fn exterior_window(
        &self,
        window: RectangularAcousticWindow,
        budget: AcousticWindowBudget,
        transmitted: bool,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<AcousticWindowSample, AcousticSpectrumError> {
        poll(cancelled)?;
        if (transmitted && window.center_m[2] < 0.0) || (!transmitted && window.center_m[2] > 0.0) {
            return Err(AcousticSpectrumError::Invalid(
                "acoustic window must lie on the selected exterior side",
            ));
        }
        observe(
            self,
            window,
            budget,
            |index| {
                let wave = &self.components[index];
                let mut sum = FieldSum([Dd::ZERO; 8]);
                if transmitted {
                    sum.add(
                        wave.transmitted,
                        propagation(wave, window.center_m, wave.receiving_normal_wavenumber)?,
                    )?;
                } else {
                    sum.add(
                        wave.incident,
                        propagation(
                            wave,
                            window.center_m,
                            C64::from_re(wave.incident_normal_wavenumber),
                        )?,
                    )?;
                    sum.add(
                        wave.reflected,
                        propagation(
                            wave,
                            window.center_m,
                            C64::from_re(-wave.incident_normal_wavenumber),
                        )?,
                    )?;
                }
                sum.finish_fields()
            },
            None,
            cancelled,
        )
        .map(|sample| sample.field)
    }
}

impl PreparedLayerAcousticSpectrum<'_, '_, '_> {
    /// Exact spatial-window observation inside this original selected layer.
    /// Depth is measured from its front face, as in [`Self::field`]. Uses the
    /// same bounded internal propagation, O(N) scratch and budgeted O(N²)
    /// cross terms; no material is queried or boundary problem solved again.
    pub fn window(
        &self,
        window: RectangularAcousticWindow,
        budget: AcousticWindowBudget,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<AcousticWindowSample, AcousticSpectrumError> {
        self.interior_window(window, budget, false, cancelled)
            .map(|sample| sample.field)
    }

    /// Average the actual coherent constitutive heat, including all velocity
    /// directions, with the same original frozen density and bulk modulus.
    /// Heat is integrated before spatial averaging, not reconstructed from
    /// mean pressure or velocity. A zero-thickness layer refuses heat.
    pub fn window_with_heat(
        &self,
        window: RectangularAcousticWindow,
        budget: AcousticWindowBudget,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<AcousticLayerWindowSample, AcousticSpectrumError> {
        self.interior_window(window, budget, true, cancelled)
    }

    fn interior_window(
        &self,
        window: RectangularAcousticWindow,
        budget: AcousticWindowBudget,
        heat: bool,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<AcousticLayerWindowSample, AcousticSpectrumError> {
        poll(cancelled)?;
        let spectrum = self.spectrum();
        let thickness = spectrum.stack.prepared_layers().layers[self.layer_index()]
            .source
            .thickness_m();
        if window.center_m[2] < 0.0 || window.center_m[2] > thickness {
            return Err(AcousticSpectrumError::Invalid(
                "acoustic window must lie within its selected layer",
            ));
        }
        if heat && thickness == 0.0 {
            return Err(AcousticSpectrumError::Invalid(
                "a zero-thickness layer has no interior heat density",
            ));
        }
        observe(
            spectrum,
            window,
            budget,
            |index| {
                let mut sum = FieldSum([Dd::ZERO; 8]);
                sum.add(
                    self.component_field(index, window.center_m[2])?,
                    propagation(
                        &spectrum.components[index],
                        [window.center_m[0], window.center_m[1], 0.0],
                        C64::ZERO,
                    )?,
                )?;
                sum.finish_fields()
            },
            heat.then_some(self.heat_scales()),
            cancelled,
        )
    }
}

struct Term {
    field: [C64; 4],
    // Weighted before squaring so weak fields in strongly dissipative material
    // retain representable heat even if their unweighted norm underflows.
    heat: [C64; 4],
}

fn observe(
    spectrum: &PreparedAcousticSpectrum<'_, '_>,
    window: RectangularAcousticWindow,
    budget: AcousticWindowBudget,
    mut component_field: impl FnMut(usize) -> Result<PlaneField, AcousticSpectrumError>,
    heat_scales: Option<[f64; 2]>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<AcousticLayerWindowSample, AcousticSpectrumError> {
    poll(cancelled)?;
    let count = spectrum.components.len();
    if budget.max_pairs == 0
        || budget.max_pairs > MAX_ACOUSTIC_WINDOW_PAIRS
        || count * (count + 1) / 2 > budget.max_pairs
    {
        return Err(AcousticSpectrumError::Invalid(
            "coherent acoustic window exceeds its explicit pair budget",
        ));
    }
    if window.center_m.iter().any(|x| !x.is_finite())
        || window.extent_m.iter().any(|x| !x.is_finite() || *x <= 0.0)
        || heat_scales.is_some_and(|s| s.iter().any(|x| !x.is_finite() || *x < 0.0))
    {
        return Err(AcousticSpectrumError::Invalid(
            "acoustic window needs finite coordinates, positive extents and passive heat",
        ));
    }
    let mut terms = Vec::with_capacity(count);
    let mut mean = FieldSum([Dd::ZERO; 8]);
    for (index, wave) in spectrum.components.iter().enumerate() {
        poll(cancelled)?;
        let field = component_field(index)?;
        mean.add(
            field,
            C64::from_re(kernel(wave.tangential_wavenumber, [0.0; 2], window.extent_m)?.to_f64()),
        )?;
        let values = [
            field.pressure,
            field.velocity[0],
            field.velocity[1],
            field.velocity[2],
        ];
        let heat = heat_scales.map_or([C64::ZERO; 4], |s| {
            [
                values[0].scale(s[1]),
                values[1].scale(s[0]),
                values[2].scale(s[0]),
                values[3].scale(s[0]),
            ]
        });
        if values.iter().chain(&heat).any(|&value| !finite(value)) {
            return Err(AcousticSpectrumError::Invalid(
                "unrepresentable acoustic window field",
            ));
        }
        terms.push(Term {
            field: values,
            heat,
        });
    }
    // [four mean squares, three complex intensities, thermal/viscous heat].
    // Products and summation stay compensated through all cancellation.
    let mut quadratic = [Dd::ZERO; 12];
    for (i, left) in terms.iter().enumerate() {
        for (j, right) in terms[..=i].iter().enumerate() {
            poll(cancelled)?;
            let weight = kernel(
                spectrum.components[i].tangential_wavenumber,
                spectrum.components[j].tangential_wavenumber,
                window.extent_m,
            )?;
            let diagonal = i == j;
            let square_weight = if diagonal { Dd::from_f64(0.5) } else { weight };
            let intensity_weight = Dd::from_f64(0.5) * weight;
            for n in 0..4 {
                quadratic[n] =
                    quadratic[n] + product(left.field[n], right.field[n])[0] * square_weight;
            }
            for axis in 0..3 {
                let a = product(left.field[0], right.field[axis + 1]);
                let b = if diagonal {
                    [Dd::ZERO; 2]
                } else {
                    product(right.field[0], left.field[axis + 1])
                };
                for part in 0..2 {
                    let n = 4 + 2 * axis + part;
                    quadratic[n] = quadratic[n] + (a[part] + b[part]) * intensity_weight;
                }
            }
            if heat_scales.is_some() {
                let heat_weight = if diagonal {
                    Dd::ONE
                } else {
                    Dd::from_f64(2.0) * weight
                };
                quadratic[10] =
                    quadratic[10] + product(left.heat[0], right.heat[0])[0] * heat_weight;
                for axis in 1..4 {
                    quadratic[11] =
                        quadratic[11] + product(left.heat[axis], right.heat[axis])[0] * heat_weight;
                }
            }
        }
    }
    let q = quadratic.map(Dd::to_f64);
    if q.iter().any(|x| !x.is_finite()) || [0, 1, 2, 3, 10, 11].iter().any(|&i| q[i] < 0.0) {
        return Err(AcousticSpectrumError::Invalid(
            "unrepresentable or negative acoustic window mean square or heat",
        ));
    }
    let mean = mean.finish_fields()?;
    poll(cancelled)?;
    Ok(AcousticLayerWindowSample {
        field: AcousticWindowSample {
            mean_pressure_pa: mean.pressure,
            mean_particle_velocity_m_s: mean.velocity,
            mean_squared_pressure_pa2: q[0],
            mean_squared_particle_velocity_m2_s2: [q[1], q[2], q[3]],
            active_intensity_w_m2: [q[4], q[6], q[8]],
            reactive_intensity_w_m2: [q[5], q[7], q[9]],
        },
        thermal_heating_w_m3: q[10],
        viscous_heating_w_m3: q[11],
    })
}

// Real and imaginary parts of a conj(b), with full-precision products before
// the pair sum. In particular, identical opposite sources cancel exactly.
fn product(a: C64, b: C64) -> [Dd; 2] {
    let [ar, ai, br, bi] = [a.re, a.im, b.re, b.im].map(Dd::from_f64);
    [ar * br + ai * bi, ai * br - ar * bi]
}

fn kernel(a: [f64; 2], b: [f64; 2], extent: [f64; 2]) -> Result<Dd, AcousticSpectrumError> {
    let arguments = [0, 1].map(|i| {
        (Dd::from_f64(a[i]) - Dd::from_f64(b[i])) * Dd::from_f64(extent[i]) * Dd::from_f64(0.5)
    });
    if arguments
        .iter()
        .any(|x| !x.hi.is_finite() || !x.lo.is_finite())
    {
        return Err(AcousticSpectrumError::Invalid(
            "unrepresentable acoustic window phase",
        ));
    }
    Ok(sinc(arguments[0]) * sinc(arguments[1]))
}

fn sinc(x: Dd) -> Dd {
    if x == Dd::ZERO {
        Dd::ONE
    } else if x.hi.abs() <= 0.5 {
        // Retain sinc(x)-1 before adding one. A narrow dark-fringe window
        // has positive O(x²) energy even when sinc(x) rounds to exactly 1.
        // The original DD trigonometry uses the same bounded series policy;
        // at |x|<=1/2 the omitted term here is below 2e-36.
        let minus_square = -(x * x);
        let mut term = minus_square / Dd::from_f64(6.0);
        let mut decrement = term;
        for n in 2..=12 {
            term = term * minus_square / Dd::from_f64(f64::from((2 * n) * (2 * n + 1)));
            decrement = decrement + term;
        }
        Dd::ONE + decrement
    } else {
        det::sin_cos_dd(x).0 / x
    }
}
