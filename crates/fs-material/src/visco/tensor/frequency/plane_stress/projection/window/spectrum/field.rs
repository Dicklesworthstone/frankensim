//! Shared positive real-field modes with original binary constitutive scales.

use super::*;

/// A linear copy of the shared temporal owner's centered basis. A material
/// channel is a binary-factored Sum, not necessarily a representable DD value;
/// retaining that scale through projection is essential before physical work,
/// duration and original relaxation time supply their compensating factors.
pub(super) struct FieldModes {
    coefficients: Vec<Dd>,
    sources: usize,
    modes: usize,
}

impl FieldModes {
    pub(super) fn new(
        frequencies: &[f64],
        window: HarmonicWindow,
    ) -> Result<Option<Self>, WindowError> {
        // Physical storage/heat/work are unweighted interval integrals. Only
        // the shared Uniform real factor applies; its absolute frequencies
        // retain both sum and beat terms without a real-peak half factor.
        let Some(band) = window.prepare_uniform_real_band_rad_s(frequencies)? else {
            return Ok(None);
        };
        let mut coefficients = Vec::with_capacity(band.mode_count() * frequencies.len());
        for mode in 0..band.mode_count() {
            for source in 0..frequencies.len() {
                coefficients.push(band.centered_coefficient_dd(mode, source)?);
            }
        }
        Ok(Some(Self {
            coefficients,
            sources: frequencies.len(),
            modes: band.mode_count(),
        }))
    }

    pub(super) fn heap_payload_bytes(&self) -> usize {
        self.coefficients.len() * core::mem::size_of::<Dd>()
    }

    fn project(&self, mode: usize, fields: &[Channel]) -> Sum {
        let mut total = Sum::default();
        let coefficients = &self.coefficients[mode * self.sources..(mode + 1) * self.sources];
        for (&coefficient, field) in coefficients.iter().zip(fields) {
            add(&mut total, scaled(field[mode % 2], coefficient));
        }
        total
    }

    pub(super) fn mean(&self, fields: &[Channel]) -> Sum {
        self.project(0, fields)
    }

    pub(super) fn add_product(
        &self,
        left: &[Channel],
        right: &[Channel],
        duration: f64,
        denominator: &[f64],
        total: &mut Sum,
    ) {
        let same_field = core::ptr::eq(left, right);
        for mode in 0..self.modes {
            let a = self.project(mode, left);
            let b = if same_field {
                a
            } else {
                self.project(mode, right)
            };
            weighted(
                total,
                a,
                b,
                Dd::from_f64((2 * mode + 1) as f64),
                &[duration],
                denominator,
            );
        }
    }
}
