//! Reusable interval packets over the ORIGINAL independently indexed optics.
//!
//! Uniform coefficients retain the single-column fast path. Indexed absorption
//! and scattering retain separate tables, not a Cartesian product or a copied
//! cell field. Every ray reads the original cell maps, sums optical depth by
//! wavelength, and only then applies exponential attenuation.

use super::{
    AbsorbingGrid, GridEnclosureBudget, Interval, MassAbsorptionSpectrum, Point3,
    VolumeTransportError, attenuate, checkpoint, exp_neg, finite, nonnegative,
};
use crate::volumes::transport::scattering::{MassScatteringSpectrum, ScatteringGrid};
use fs_exec::Cx;

/// Default cap on retained coefficient intervals, across both optical tables.
/// At 16 bytes per interval this permits 5 MiB of initialized table entries.
/// Borrowed fields, object/stack storage, spare capacity and allocator overhead
/// are not part of this coefficient-count budget.
pub const MAX_PREPARED_EXTINCTION_VALUES: usize = 327_680;

enum Coefficients<const N: usize> {
    Uniform([Interval; N]),
    Indexed {
        absorption: Vec<[Interval; N]>,
        scattering: Vec<[Interval; N]>,
    },
}

/// Prepared spectral enclosures immutably bound to the original native medium.
///
/// Density, cell maps and spectral/phase declarations are never copied or
/// replaced. Preparation alone allocates bounded spectral tables; segment
/// queries allocate nothing. This bounds the retained numerical model, not
/// calibration, spatial/spectral discretization, LTE emission or in-scattering.
pub struct PreparedGridExtinction<'g, 'field, const N: usize> {
    grid: &'g AbsorbingGrid<'field>,
    scattering: Option<&'g ScatteringGrid<'field>>,
    coefficients: Coefficients<N>,
}

impl<'field> AbsorbingGrid<'field> {
    /// Prepare wavelength packets once for any number of finite ray segments.
    /// `max_coefficient_values` bounds all retained coefficient intervals before
    /// allocation. Uniform optics need N values; indexed optics need N values
    /// per declared absorption model, including unused entries.
    pub fn prepare_extinction_enclosure<'g, const N: usize>(
        &'g self,
        cx: &Cx<'_>,
        wavelengths_nm: &[f64; N],
        max_coefficient_values: usize,
    ) -> Result<PreparedGridExtinction<'g, 'field, N>, VolumeTransportError> {
        PreparedGridExtinction::prepare(cx, self, None, wavelengths_nm, max_coefficient_values)
    }
}

impl<'field> ScatteringGrid<'field> {
    /// Enclose a global extinction law. Indexed optics deliberately refuse
    /// `Palette`: they have no representative global coefficient packet.
    /// Use `prepare_extinction_enclosure` for their original cell-selected laws.
    pub fn extinction_coefficients_enclosure<const N: usize>(
        &self,
        wavelengths_nm: &[f64; N],
    ) -> Result<[Interval; N], VolumeTransportError> {
        if self.scattering_palette().is_some() { return Err(VolumeTransportError::Palette); }
        let mut values = self.absorption_grid().extinction_coefficients_enclosure(wavelengths_nm)?;
        let scattering = self.scattering_at_cell(0)?;
        for (value, &wavelength) in values.iter_mut().zip(wavelengths_nm) {
            *value = positive_sum(*value, scattering.coefficient_enclosure(wavelength)?)?;
        }
        Ok(values)
    }

    /// Prepare the original absorption and scattering maps independently.
    /// Uniform optics need N intervals. Indexed optics need N per absorption
    /// and scattering model; no Cartesian palette-pair table is allocated.
    pub fn prepare_extinction_enclosure<'g, const N: usize>(
        &'g self,
        cx: &Cx<'_>,
        wavelengths_nm: &[f64; N],
        max_coefficient_values: usize,
    ) -> Result<PreparedGridExtinction<'g, 'field, N>, VolumeTransportError> {
        PreparedGridExtinction::prepare(
            cx, self.absorption_grid(), Some(self), wavelengths_nm, max_coefficient_values,
        )
    }

    /// Enclose unscattered beam transmission with the default coefficient cap.
    /// Repeated rays should retain a prepared packet to amortize preparation.
    /// Phase response and LTE emission never enter direct transmittance.
    pub fn transmittance_enclosure<const N: usize>(
        &self,
        cx: &Cx<'_>,
        start: Point3,
        end: Point3,
        wavelengths_nm: &[f64; N],
        budget: GridEnclosureBudget,
    ) -> Result<[Interval; N], VolumeTransportError> {
        self.prepare_extinction_enclosure(cx, wavelengths_nm, MAX_PREPARED_EXTINCTION_VALUES)?
            .transmittance(cx, start, end, budget)
    }
}

impl<'g, 'field, const N: usize> PreparedGridExtinction<'g, 'field, N> {
    fn prepare(
        cx: &Cx<'_>,
        grid: &'g AbsorbingGrid<'field>,
        scattering: Option<&'g ScatteringGrid<'field>>,
        wavelengths: &[f64; N],
        max_values: usize,
    ) -> Result<Self, VolumeTransportError> {
        checkpoint(cx)?;
        if N == 0 { return Err(VolumeTransportError::Sample); }
        let indexed = grid.absorption_palette().is_some()
            || scattering.is_some_and(|s| s.scattering_palette().is_some());
        if !indexed {
            if N > max_values { return Err(VolumeTransportError::Budget); }
            let coefficients = match scattering {
                Some(medium) => medium.extinction_coefficients_enclosure(wavelengths)?,
                None => grid.extinction_coefficients_enclosure(wavelengths)?,
            };
            checkpoint(cx)?;
            return Ok(Self { grid, scattering, coefficients: Coefficients::Uniform(coefficients) });
        }
        let a_count = grid.absorption_palette().map_or(1, |p| p.values().len());
        let s_count = scattering.map_or(0, |s| s.scattering_palette().map_or(1, |p| p.values().len()));
        let count = a_count.checked_add(s_count).and_then(|n| n.checked_mul(N))
            .ok_or(VolumeTransportError::Budget)?;
        if count > max_values { return Err(VolumeTransportError::Budget); }
        let mut a = Vec::new();
        let mut s = Vec::new();
        a.try_reserve_exact(a_count).map_err(|_| VolumeTransportError::Budget)?;
        s.try_reserve_exact(s_count).map_err(|_| VolumeTransportError::Budget)?;
        for index in 0..a_count {
            checkpoint(cx)?;
            let spectrum = match grid.absorption_palette() {
                Some(palette) => &palette.values()[index],
                None => grid.absorption_at_cell(0)?,
            };
            a.push(absorption_packet(spectrum, wavelengths)?);
        }
        if let Some(medium) = scattering {
            for index in 0..s_count {
                checkpoint(cx)?;
                let spectrum = match medium.scattering_palette() {
                    Some(palette) => &palette.values()[index].scattering,
                    None => medium.scattering_at_cell(0)?,
                };
                s.push(scattering_packet(spectrum, wavelengths)?);
            }
        }
        checkpoint(cx)?;
        Ok(Self { grid, scattering, coefficients: Coefficients::Indexed { absorption: a, scattering: s } })
    }

    /// Bound attenuation over the finite endpoint segment without allocation.
    /// Original cell and caller node-work budgets remain binding. Cancellation
    /// and exhausted/unsupported arithmetic never return a partial packet.
    pub fn transmittance(
        &self,
        cx: &Cx<'_>,
        start: Point3,
        end: Point3,
        budget: GridEnclosureBudget,
    ) -> Result<[Interval; N], VolumeTransportError> {
        checkpoint(cx)?;
        let (absorption, scattering) = match &self.coefficients {
            Coefficients::Uniform(coefficients) => {
                let column = self.grid.column_density_enclosure(cx, start, end, budget)?;
                let values = attenuate(column.kg_per_m2, coefficients)?;
                checkpoint(cx)?;
                return Ok(values);
            }
            Coefficients::Indexed { absorption, scattering } => (absorption, scattering),
        };
        let mut depth = [Interval::point(0.0); N];
        self.grid.visit_column(cx, start, end, budget, |cell, column| {
            let a_index = self.grid.absorption_palette().map_or(0, |p| p.cell_indices()[cell] as usize);
            let a = &absorption[a_index];
            let s = self.scattering.map(|medium| {
                let index = medium.scattering_palette().map_or(0, |p| p.cell_indices()[cell] as usize);
                &scattering[index]
            });
            for lane in 0..N {
                if lane % 64 == 0 { checkpoint(cx)?; }
                let extinction = match s {
                    Some(s) => positive_sum(a[lane], s[lane])?,
                    None => a[lane],
                };
                if extinction.hi() == 0.0 { continue; }
                depth[lane] = nonnegative(depth[lane] + column * extinction);
                if !finite(depth[lane]) { return Err(VolumeTransportError::Unrepresentable); }
            }
            Ok(())
        })?;
        let mut result = [Interval::point(1.0); N];
        for (lane, tau) in depth.into_iter().enumerate() {
            if lane % 64 == 0 { checkpoint(cx)?; }
            result[lane] = exp_neg(tau).ok_or(VolumeTransportError::Unrepresentable)?;
        }
        checkpoint(cx)?;
        Ok(result)
    }
}

fn positive_sum(a: Interval, b: Interval) -> Result<Interval, VolumeTransportError> {
    let value = if a.hi() == 0.0 { b } else if b.hi() == 0.0 { a } else { nonnegative(a + b) };
    if finite(value) { Ok(value) } else { Err(VolumeTransportError::Unrepresentable) }
}

fn absorption_packet<const N: usize>(
    spectrum: &MassAbsorptionSpectrum, wavelengths: &[f64; N],
) -> Result<[Interval; N], VolumeTransportError> {
    let mut values = [Interval::point(0.0); N];
    for (value, &wavelength) in values.iter_mut().zip(wavelengths) {
        *value = spectrum.coefficient_enclosure(wavelength)?;
    }
    Ok(values)
}

fn scattering_packet<const N: usize>(
    spectrum: &MassScatteringSpectrum, wavelengths: &[f64; N],
) -> Result<[Interval; N], VolumeTransportError> {
    let mut values = [Interval::point(0.0); N];
    for (value, &wavelength) in values.iter_mut().zip(wavelengths) {
        *value = spectrum.coefficient_enclosure(wavelength)?;
    }
    Ok(values)
}
