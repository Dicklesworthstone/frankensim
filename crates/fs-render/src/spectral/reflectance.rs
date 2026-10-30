//! Declared passive visible-band reflectance and the shared fractional table.
//!
//! Bounds apply to the represented piecewise-linear law. They do not establish
//! measurement provenance, interpolation accuracy or a scattering model.

/// Integration domain (nm): the standard visible band.
pub const LAMBDA_MIN: f64 = 380.0;
/// Upper edge of the visible band (nm).
pub const LAMBDA_MAX: f64 = 780.0;

const MAX_KNOTS: usize = 16;

/// A declared passive reflectance spectrum with exact black and white support.
///
/// The law covers 380--780 nm with 2--16 strictly increasing wavelength knots.
/// Values are finite fractions in [0,1], interpolated linearly in wavelength.
/// No RGB fit, extrapolation, temperature or thermal emission is implied.
/// Admission does not make a supplied table a calibrated measurement or bound
/// its interpolation, wavelength-resolution or material-model error.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SpectralReflectance {
    table: VisibleFractionTable,
}

impl SpectralReflectance {
    /// Declare wavelength-independent reflectance within the visible band.
    pub fn constant(reflectance: f64) -> Result<Self, SpectralReflectanceError> {
        Self::table(&[(LAMBDA_MIN, reflectance), (LAMBDA_MAX, reflectance)])
    }

    /// Admit (wavelength_nm, reflectance) knots without sorting or clipping.
    ///
    /// Endpoints must be exactly 380 and 780 nm. At most 16 knots are copied
    /// into inline storage, so the immutable spectrum remains allocation-free.
    pub fn table(knots: &[(f64, f64)]) -> Result<Self, SpectralReflectanceError> {
        Ok(Self {
            table: VisibleFractionTable::table(knots).map_err(SpectralReflectanceError::from)?,
        })
    }

    /// Evaluate the declared reflectance at a wavelength in the closed band.
    pub fn at(&self, wavelength_nm: f64) -> Result<f64, SpectralReflectanceError> {
        self.table
            .at(wavelength_nm)
            .map_err(SpectralReflectanceError::from)
    }

    /// Declared knots, excluding unused inline storage.
    #[must_use]
    pub fn knots(&self) -> &[(f64, f64)] {
        self.table.knots()
    }
}

/// Why a declared passive visible-band reflectance cannot be used.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SpectralReflectanceError {
    /// Require 2--16 strictly increasing finite knots from 380 to 780 nm.
    Knots,
    /// Every reflectance must be finite and in [0,1].
    Reflectance,
    /// Evaluation requires a finite wavelength within 380--780 nm.
    Wavelength,
}

impl core::fmt::Display for SpectralReflectanceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Knots => f.write_str("invalid visible-band reflectance knots"),
            Self::Reflectance => f.write_str("reflectance must be finite and in [0,1]"),
            Self::Wavelength => f.write_str("wavelength outside the visible reflectance band"),
        }
    }
}

impl core::error::Error for SpectralReflectanceError {}

impl From<VisibleFractionTableError> for SpectralReflectanceError {
    fn from(error: VisibleFractionTableError) -> Self {
        match error {
            VisibleFractionTableError::Knots => Self::Knots,
            VisibleFractionTableError::Fraction => Self::Reflectance,
            VisibleFractionTableError::Wavelength => Self::Wavelength,
        }
    }
}

/// Numerical owner shared by semantically distinct reflectance/emissivity laws.
/// Values are stored directly, never as a rounded complement of another law.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct VisibleFractionTable {
    knots: [(f64, f64); MAX_KNOTS],
    len: u8,
    positive: bool,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum VisibleFractionTableError {
    Knots,
    Fraction,
    Wavelength,
}

impl VisibleFractionTable {
    pub(crate) fn table(knots: &[(f64, f64)]) -> Result<Self, VisibleFractionTableError> {
        if !(2..=MAX_KNOTS).contains(&knots.len())
            || knots[0].0 != LAMBDA_MIN
            || knots[knots.len() - 1].0 != LAMBDA_MAX
        {
            return Err(VisibleFractionTableError::Knots);
        }
        let mut accepted = Self {
            knots: [(0.0, 0.0); MAX_KNOTS],
            len: knots.len() as u8,
            positive: false,
        };
        for (index, &(wavelength_nm, value)) in knots.iter().enumerate() {
            if !wavelength_nm.is_finite()
                || !(LAMBDA_MIN..=LAMBDA_MAX).contains(&wavelength_nm)
                || (index > 0 && wavelength_nm <= knots[index - 1].0)
            {
                return Err(VisibleFractionTableError::Knots);
            }
            if !value.is_finite() || !(0.0..=1.0).contains(&value) {
                return Err(VisibleFractionTableError::Fraction);
            }
            // One physical zero has one declaration identity.
            let value = if value == 0.0 { 0.0 } else { value };
            accepted.knots[index] = (wavelength_nm, value);
            accepted.positive |= value > 0.0;
        }
        Ok(accepted)
    }

    pub(crate) fn at(&self, wavelength_nm: f64) -> Result<f64, VisibleFractionTableError> {
        if !wavelength_nm.is_finite() || !(LAMBDA_MIN..=LAMBDA_MAX).contains(&wavelength_nm) {
            return Err(VisibleFractionTableError::Wavelength);
        }
        let knots = self.knots();
        let index = knots.partition_point(|&(wavelength, _)| wavelength <= wavelength_nm) - 1;
        let (left_wavelength, left_value) = knots[index];
        if wavelength_nm == left_wavelength {
            return Ok(left_value);
        }
        let (right_wavelength, right_value) = knots[index + 1];
        if left_value == right_value {
            return Ok(left_value);
        }
        let fraction = (wavelength_nm - left_wavelength) / (right_wavelength - left_wavelength);
        // Preserve nonnegative support near zero without subtracting nearly
        // equal endpoint values or imposing an artificial density floor.
        Ok((1.0 - fraction) * left_value + fraction * right_value)
    }

    pub(crate) fn knots(&self) -> &[(f64, f64)] {
        &self.knots[..usize::from(self.len)]
    }

    pub(crate) const fn has_positive_value(&self) -> bool {
        self.positive
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exact_black_white_and_constant_reflectance() {
        for value in [0.0, 0.25, 1.0] {
            let curve = SpectralReflectance::constant(value).unwrap();
            for wavelength in [380.0, 382.5, 447.5, 555.0, 777.5, 780.0] {
                assert_eq!(curve.at(wavelength).unwrap().to_bits(), value.to_bits());
            }
            assert_eq!(curve.table.has_positive_value(), value > 0.0);
        }
        let zero = SpectralReflectance::constant(-0.0).unwrap();
        assert_eq!(zero, SpectralReflectance::constant(0.0).unwrap());
        assert!(zero.knots().iter().all(|&(_, value)| value.to_bits() == 0));
    }

    #[test]
    fn wavelength_linear_interpolation_uses_declared_knots() {
        let knots = [
            (380.0, 0.125),
            (450.0, 0.375),
            (550.0, 0.875),
            (780.0, 0.25),
        ];
        let curve = SpectralReflectance::table(&knots).unwrap();
        assert_eq!(curve.knots(), knots);
        for &(wavelength, value) in &knots {
            assert_eq!(curve.at(wavelength).unwrap(), value);
        }
        for (wavelength, expected) in [(415.0, 0.25), (500.0, 0.625), (665.0, 0.5625)] {
            assert_eq!(curve.at(wavelength).unwrap(), expected);
        }
    }

    #[test]
    fn tiny_reflectance_is_not_lost_by_complement_storage() {
        let tiny = 1.0e-300;
        assert_eq!(1.0 - (1.0 - tiny), 0.0);
        let constant = SpectralReflectance::constant(tiny).unwrap();
        assert_eq!(constant.at(555.0).unwrap(), tiny);
        let ramp = SpectralReflectance::table(&[(380.0, 0.0), (780.0, tiny)]).unwrap();
        assert_eq!(ramp.at(580.0).unwrap(), tiny / 2.0);
        assert!(ramp.table.has_positive_value());
    }

    #[test]
    fn narrow_zero_supported_band_stays_passive() {
        let curve = SpectralReflectance::table(&[
            (380.0, 0.0),
            (500.0, 0.0),
            (501.0, 1.0),
            (502.0, 0.0),
            (780.0, 0.0),
        ])
        .unwrap();
        for wavelength in [380.0, 499.99, 500.0, 502.0, 550.0, 780.0] {
            assert_eq!(curve.at(wavelength).unwrap(), 0.0);
        }
        assert_eq!(curve.at(500.25).unwrap(), 0.25);
        assert_eq!(curve.at(501.25).unwrap(), 0.75);
        for sample in 0..=1600 {
            let value = curve.at(380.0 + f64::from(sample) * 0.25).unwrap();
            assert!((0.0..=1.0).contains(&value));
        }
    }

    #[test]
    fn table_capacity_is_explicit_and_bounded() {
        let maximum = core::array::from_fn::<_, MAX_KNOTS, _>(|index| {
            let fraction = index as f64 / (MAX_KNOTS - 1) as f64;
            (380.0 + 400.0 * fraction, fraction)
        });
        assert_eq!(
            SpectralReflectance::table(&maximum).unwrap().knots().len(),
            MAX_KNOTS
        );
        let excessive =
            core::array::from_fn::<_, 17, _>(|index| (380.0 + 25.0 * index as f64, 0.5));
        for knots in [&[][..], &[(380.0, 0.0)][..], &excessive[..]] {
            assert_eq!(
                SpectralReflectance::table(knots),
                Err(SpectralReflectanceError::Knots)
            );
        }
    }

    #[test]
    fn malformed_knots_and_active_reflectance_refuse() {
        for knots in [
            vec![(381.0, 0.0), (780.0, 1.0)],
            vec![(380.0, 0.0), (779.0, 1.0)],
            vec![(380.0, 0.0), (500.0, 0.5), (500.0, 0.6), (780.0, 1.0)],
            vec![(380.0, 0.0), (600.0, 0.5), (500.0, 0.6), (780.0, 1.0)],
            vec![(380.0, 0.0), (f64::NAN, 0.5), (780.0, 1.0)],
        ] {
            assert_eq!(
                SpectralReflectance::table(&knots),
                Err(SpectralReflectanceError::Knots)
            );
        }
        for value in [-0.01, 1.01, f64::INFINITY, f64::NEG_INFINITY, f64::NAN] {
            assert_eq!(
                SpectralReflectance::constant(value),
                Err(SpectralReflectanceError::Reflectance)
            );
        }
    }

    #[test]
    fn evaluation_outside_declared_band_refuses() {
        let curve = SpectralReflectance::constant(0.5).unwrap();
        for wavelength in [379.0, 781.0, f64::NAN, f64::NEG_INFINITY, f64::INFINITY] {
            assert_eq!(
                curve.at(wavelength),
                Err(SpectralReflectanceError::Wavelength)
            );
        }
    }
}
