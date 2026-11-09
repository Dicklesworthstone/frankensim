//! One original-frequency admission and basis builder for both storage forms.

use super::*;

#[derive(Clone, Copy)]
pub(super) enum Units {
    Hertz,
    Radians,
}

impl Units {
    fn width(self, window: HarmonicWindow, frequency: Dd) -> Result<Dd, WindowError> {
        match self {
            Self::Hertz => window.width_hz(frequency, 2.0),
            Self::Radians => window.half_width_extended(frequency),
        }
    }

    fn rotation(self, window: HarmonicWindow, frequency: Dd) -> Result<[Dd; 2], WindowError> {
        match self {
            Self::Hertz => window.rotation_hz(frequency),
            Self::Radians => window.rotation_parts(frequency),
        }
    }
}

pub(super) struct Spec {
    minimum: f64,
    center: Dd,
    pub(super) modes: usize,
    pub(super) real_field: bool,
}

pub(super) fn admit(
    window: HarmonicWindow,
    frequencies: &[f64],
    units: Units,
    real_field: bool,
) -> Result<Option<Spec>, WindowError> {
    let maximum_tones = MAX_UNIFORM_HARMONIC_BAND_TONES / if real_field { 2 } else { 1 };
    if frequencies.len() > maximum_tones {
        return Ok(None);
    }
    if frequencies.iter().any(|f| !f.is_finite()) {
        return Err(WindowError::InvalidFrequency);
    }
    if window.profile != HarmonicWindowProfile::Uniform || frequencies.len() < 3 {
        return Ok(None);
    }
    let (minimum, maximum) = frequencies.iter().fold(
        (f64::INFINITY, f64::NEG_INFINITY),
        |(lo, hi), &frequency| (lo.min(frequency), hi.max(frequency)),
    );
    let center = if real_field {
        Dd::ZERO
    } else {
        checked(Dd::from_f64(maximum) - Dd::from_f64(minimum))? * HALF
    };
    let radius_frequency = if real_field {
        Dd::from_f64(minimum.abs().max(maximum.abs()))
    } else {
        center
    };
    let width = units.width(window, radius_frequency)?;
    if HALF.lt(width) {
        return Ok(None);
    }
    // A real field has the degrees of its mirrored +/-frequency expansion,
    // even though parity lets us retain only the original source columns.
    let minimum_modes =
        (frequencies.len() * if real_field { 2 } else { 1 }).min(MAX_UNIFORM_HARMONIC_BAND_MODES);
    let radius = width.to_f64().next_up();
    let mut leading = 1.0;
    let mut retained = MAX_UNIFORM_HARMONIC_BAND_MODES;
    for modes in 1..=MAX_UNIFORM_HARMONIC_BAND_MODES {
        let odd = (2 * modes + 1) as f64;
        leading *= radius / odd;
        let tail = odd * leading * leading / (1.0 - radius * radius / (odd * (odd + 2.0)));
        // Margin for estimating the tail in f64. At radius 0.5, 17 retained
        // modes have an exact-arithmetic omitted modal tail below 4.15e-50.
        if modes >= minimum_modes && tail <= 4e-50 {
            retained = modes;
            break;
        }
    }
    Ok(Some(Spec {
        minimum: if real_field { 0.0 } else { minimum },
        center,
        modes: retained,
        real_field,
    }))
}

pub(super) fn mode_count(
    window: HarmonicWindow,
    frequencies: &[f64],
    units: Units,
    real_field: bool,
) -> Result<Option<usize>, WindowError> {
    Ok(admit(window, frequencies, units, real_field)?.map(|spec| spec.modes))
}

pub(super) fn fill(
    window: HarmonicWindow,
    frequencies: &[f64],
    units: Units,
    spec: &Spec,
    rotations: &mut [[Dd; 2]],
    coefficients: &mut [Dd],
) -> Result<(), WindowError> {
    for (i, &frequency) in frequencies.iter().enumerate() {
        let relative = (Dd::from_f64(frequency) - Dd::from_f64(spec.minimum)) - spec.center;
        rotations[i] = units.rotation(window, relative)?;
        let a = units.width(window, relative)?;
        let negative_square = -(a * a);
        let mut leading = Dd::ONE;
        for l in 0..spec.modes {
            let mut series = Dd::ONE;
            let mut term = Dd::ONE;
            for k in 1..=16 {
                term = term * negative_square / Dd::from_f64((2 * k * (2 * l + 2 * k + 1)) as f64);
                series = series + term;
            }
            coefficients[l * frequencies.len() + i] = leading * series;
            leading = leading * a / Dd::from_f64((2 * l + 3) as f64);
        }
    }
    Ok(())
}

pub(super) fn prepare_owned(
    window: HarmonicWindow,
    frequencies: &[f64],
    units: Units,
    real_field: bool,
) -> Result<Option<PreparedUniformHarmonicBand>, WindowError> {
    let Some(spec) = admit(window, frequencies, units, real_field)? else {
        return Ok(None);
    };
    let mut rotations = vec![[Dd::ZERO; 2]; frequencies.len()];
    let mut coefficients = vec![Dd::ZERO; spec.modes * frequencies.len()];
    fill(
        window,
        frequencies,
        units,
        &spec,
        &mut rotations,
        &mut coefficients,
    )?;
    Ok(Some(PreparedUniformHarmonicBand {
        rotations,
        coefficients,
        modes: spec.modes,
        real_field: spec.real_field,
    }))
}

impl HarmonicWindow {
    /// Analytic-band admission retaining original radians/second, without Hz conversion.
    pub fn uniform_band_mode_count_rad_s(
        self,
        frequencies: &[f64],
    ) -> Result<Option<usize>, WindowError> {
        mode_count(self, frequencies, Units::Radians, false)
    }

    /// Cache an analytic band on its original radians/second source clock.
    pub fn prepare_uniform_band_rad_s(
        self,
        frequencies: &[f64],
    ) -> Result<Option<PreparedUniformHarmonicBand>, WindowError> {
        prepare_owned(self, frequencies, Units::Radians, false)
    }

    /// Real-peak admission: 3..32 tones and max(abs(2*pi*f))*duration/2 <= 0.5.
    /// The full absolute frequencies are retained, including sum-frequency products.
    pub fn uniform_real_band_mode_count_hz(
        self,
        frequencies: &[f64],
    ) -> Result<Option<usize>, WindowError> {
        mode_count(self, frequencies, Units::Hertz, true)
    }

    /// Real-peak admission on original rad/s: max(abs(omega))*duration/2 <= 0.5.
    pub fn uniform_real_band_mode_count_rad_s(
        self,
        frequencies: &[f64],
    ) -> Result<Option<usize>, WindowError> {
        mode_count(self, frequencies, Units::Radians, true)
    }

    /// Cache actual real-peak products on the original Hz clock; no half factor.
    pub fn prepare_uniform_real_band_hz(
        self,
        frequencies: &[f64],
    ) -> Result<Option<PreparedUniformHarmonicBand>, WindowError> {
        prepare_owned(self, frequencies, Units::Hertz, true)
    }

    /// Cache actual real-peak products on the original rad/s clock; no half factor.
    pub fn prepare_uniform_real_band_rad_s(
        self,
        frequencies: &[f64],
    ) -> Result<Option<PreparedUniformHarmonicBand>, WindowError> {
        prepare_owned(self, frequencies, Units::Radians, true)
    }
}
