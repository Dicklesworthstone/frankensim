//! Source-amplitude derivatives of the actual finite-area, finite-time energy.
//!
//! The original modal reconstruction supplies one complex field basis per
//! authored control bank. Controls act BEFORE real-field products. In a short
//! shutter the carrier-sum terms matter, so the Hessian is real-linear in
//! complex Cartesian coordinates, not generally a complex Hermitian matrix.

use super::*;

#[path = "control/action.rs"]
mod action;
use action::{EnergyAction, SourceGroup};

/// Whole preparation and repeated-query bounds for coherent source control.
#[derive(Clone, Copy, Debug)]
pub struct AngularSpectrumApertureControlLimits {
    /// At most 64 independently weighted original coherent banks.
    pub maximum_controls: usize,
    /// N times the total original line count, including the fixed bank.
    /// Bounds the one-time original modal reconstruction; no FFT is added.
    pub maximum_point_line_pairs: usize,
    /// N times the sum of distinct-frequency counts of all original banks.
    /// Bounds retained basis fields and each forward/source-adjoint sweep.
    pub maximum_source_field_entries: usize,
    /// N*F*(width+height+F), where F is the union of original frequencies.
    /// Bounds the three separable kernel passes on six real field components.
    pub maximum_kernel_pairs: usize,
    /// Peak logical storage, including original reconstruction scratch,
    /// sparse source fields, union-frequency fields, kernels and gradient.
    /// Borrowed owners and allocator overhead are excluded.
    pub maximum_workspace_bytes: usize,
}

/// The same physical energy evaluation and its Cartesian source derivative.
/// For a source multiplier a=x+i*y, re=dJ/dx and im=dJ/dy. J is mean field
/// energy density in J/m³; each multiplier is dimensionless and multiplies
/// the original bank's SI peak fields on its original common phase clock.
pub struct AngularSpectrumApertureEnergyGradient<'a> {
    pub mean_field_energy_density_j_m3: f64,
    pub source_gradient: &'a [C64],
}

/// Reusable finite-footprint optical energy objective for coherent controls.
///
/// The observed real field is the fixed original bank plus sum_s a_s F_s,
/// where each F_s is an original coherent exposure bank. The objective is
/// <(epsilon E.E + mu H.H)/2>, averaged over the original rectangle and
/// normalized shutter. It is an energy DENSITY, not absorbed heat, a body
/// energy or detector power. The original observed half-space must admit its
/// lossless nondispersive storage law. A lossy exit does not forbid a valid
/// left-side observation.
///
/// Preparation retains sparse original modal fields and the existing exact
/// aperture/shutter kernels. Queries combine the fields first, apply that
/// real quadratic operator, then pull back to the original control order.
/// This preserves cancellation before products, needs no dense source
/// Hessian, and performs no repeated Maxwell solve, FFT or allocation.
/// Separable x/y/time passes cost O(N*F*(width+height+F) + N*sum_s F_s)
/// per value/gradient or Hessian action. Even/odd channels keep small positive
/// window variances without subtracting nearly equal sum/difference kernels.
///
/// All banks retain their original material, geometry and source bindings.
/// They are deterministic coherent controls; independently sampled covariance
/// factors do not acquire mutual phase coherence by becoming control inputs.
/// The multipliers describe authored spectral amplitudes; this is neither a
/// causal time-domain controller nor a material/geometry derivative. The
/// rectangle observes the original periodic continuation, with no added
/// isolated-aperture diffraction or cubature/convergence certificate.
/// DD arithmetic does not certify relative accuracy at every near-null field.
pub struct PreparedAngularSpectrumApertureEnergyControl<'e, 'b, 's, 'm> {
    fixed: &'e AngularSpectrumExposure<'b, 's, 'm>,
    controls: &'e [AngularSpectrumExposure<'b, 's, 'm>],
    aperture: AngularSpectrumAperture,
    window: HarmonicWindow,
    distance: f64,
    action: EnergyAction,
    workspace_bytes: usize,
    source_field_entries: usize,
    kernel_pairs: usize,
}

impl<'e, 'b, 's, 'm> PreparedAngularSpectrumApertureEnergyControl<'e, 'b, 's, 'm> {
    /// Bind one-sided observation semantics of the original exposure owners.
    /// An empty fixed bank permits pure source design; zero controls preserve
    /// the fixed observation. Equal-frequency banks interfere coherently.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        fixed: &'e AngularSpectrumExposure<'b, 's, 'm>,
        controls: &'e [AngularSpectrumExposure<'b, 's, 'm>],
        side: AngularSpectrumSide,
        distance_from_face_m: f64,
        aperture: AngularSpectrumAperture,
        window: HarmonicWindow,
        limits: AngularSpectrumApertureControlLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, AngularSpectrumError> {
        Self::prepare(
            fixed,
            controls,
            ExposureSide::OneSided(side),
            distance_from_face_m,
            aperture,
            window,
            limits,
            &mut cancelled,
        )
    }

    /// Two-face counterpart. A control may correlate both physical faces;
    /// its complete incoming/outgoing field is differentiated together.
    #[allow(clippy::too_many_arguments)]
    pub fn new_illumination(
        fixed: &'e AngularSpectrumExposure<'b, 's, 'm>,
        controls: &'e [AngularSpectrumExposure<'b, 's, 'm>],
        side: AngularSpectrumIlluminationSide,
        distance_from_face_m: f64,
        aperture: AngularSpectrumAperture,
        window: HarmonicWindow,
        limits: AngularSpectrumApertureControlLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, AngularSpectrumError> {
        Self::prepare(
            fixed,
            controls,
            ExposureSide::TwoSided(side),
            distance_from_face_m,
            aperture,
            window,
            limits,
            &mut cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn prepare(
        fixed: &'e AngularSpectrumExposure<'b, 's, 'm>,
        controls: &'e [AngularSpectrumExposure<'b, 's, 'm>],
        side: ExposureSide,
        distance: f64,
        aperture: AngularSpectrumAperture,
        window: HarmonicWindow,
        limits: AngularSpectrumApertureControlLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, AngularSpectrumError> {
        poll(cancelled)?;
        if controls.len() > 64 || controls.len() > limits.maximum_controls {
            return Err(bad("optical control count exceeds its budget"));
        }
        if !distance.is_finite() || distance < 0.0 {
            return Err(bad(
                "optical control distance must be finite and nonnegative",
            ));
        }
        let anchor = fixed.anchor;
        let n = anchor.modes.len();
        let mut count = 0;
        let mut lines = 0;
        let mut original_workspace = 0;
        for bank in core::iter::once(fixed).chain(controls.iter()) {
            poll(cancelled)?;
            if bank.anchor.grid() != anchor.grid()
                || !bank
                    .anchor
                    .source()
                    .same_sources_and_geometry(anchor.source())
            {
                return Err(bad(
                    "optical controls require the same original materials, layers and grid",
                ));
            }
            count = add(count, bank.frequency_count())?;
            lines = add(lines, bank.line_count())?;
            original_workspace = original_workspace.max(bank.workspace_bytes());
        }
        let source_field_entries = mul(n, count)?;
        let preliminary = add(
            core::mem::size_of::<Self>(),
            mul(count, core::mem::size_of::<f64>())?,
        )?;
        if mul(n, lines)? > limits.maximum_point_line_pairs
            || source_field_entries > limits.maximum_source_field_entries
            || preliminary > limits.maximum_workspace_bytes
        {
            return Err(bad(
                "optical control reconstruction exceeds its work or storage budget",
            ));
        }
        let medium = side.momentum_material(anchor.source())?;
        let epsilon = linear_electric_permittivity(medium).map_err(ReflectionError::from)?;
        let mu = linear_magnetic_permeability(medium).map_err(ReflectionError::from)?;
        let mut frequencies = reserve(count)?;
        for bank in core::iter::once(fixed).chain(controls.iter()) {
            for (i, &omega) in bank.frequencies.iter().enumerate() {
                if i % 64 == 0 {
                    poll(cancelled)?;
                }
                frequencies.push(omega);
            }
        }
        frequencies.sort_unstable_by(f64::total_cmp);
        frequencies.dedup();
        poll(cancelled)?;
        let f = frequencies.len();
        let nf = mul(n, f)?;
        let grid = anchor.grid();
        let kernel_pairs = mul(nf, add(add(grid.width, grid.height)?, f)?)?;
        let axis_pairs = if f == 0 {
            0
        } else {
            add(mul(grid.width, grid.width)?, mul(grid.height, grid.height)?)?
        };
        let mut workspace_bytes = add(preliminary, original_workspace)?;
        for bytes in [
            mul(count, core::mem::size_of::<SourceGroup>())?,
            mul(source_field_entries, core::mem::size_of::<FieldParts>())?,
            mul(mul(nf, 4)?, core::mem::size_of::<FieldParts>())?,
            mul(
                add(axis_pairs, mul(f, f)?)?,
                core::mem::size_of::<RealKernel>(),
            )?,
            mul(controls.len(), core::mem::size_of::<C64>())?,
        ] {
            workspace_bytes = add(workspace_bytes, bytes)?;
        }
        if kernel_pairs > limits.maximum_kernel_pairs
            || workspace_bytes > limits.maximum_workspace_bytes
        {
            return Err(bad(
                "optical control kernels exceed their work or storage budget",
            ));
        }
        let spatial = if f == 0 {
            None
        } else {
            Some(Spatial::<RealKernel>::prepare(anchor, aperture, cancelled)?)
        };
        let temporal = pair_kernels(f, window, |i| -frequencies[i], cancelled)?;
        let mut groups = reserve(count)?;
        for (bank_index, bank) in core::iter::once(fixed).chain(controls.iter()).enumerate() {
            for group in 0..bank.frequency_count() {
                poll(cancelled)?;
                let omega = bank.frequencies[group];
                let frequency = frequencies
                    .binary_search_by(|v| v.total_cmp(&omega))
                    .map_err(|_| bad("original optical control frequency is unavailable"))?;
                let mut fields = bank.frequency_modal_parts(group, side, distance, cancelled)?;
                center_fields(
                    &mut fields,
                    omega,
                    |i| anchor.modes[i].kx,
                    |i| anchor.modes[i].ky,
                    aperture,
                    window,
                    cancelled,
                )?;
                groups.push(SourceGroup {
                    control: bank_index.checked_sub(1),
                    frequency,
                    fields,
                });
            }
        }
        let action = EnergyAction::new(
            controls.len(),
            n,
            frequencies,
            spatial,
            temporal,
            groups,
            epsilon,
            mu,
        )?;
        poll(cancelled)?;
        Ok(Self {
            fixed,
            controls,
            aperture,
            window,
            distance,
            action,
            workspace_bytes,
            source_field_entries,
            kernel_pairs,
        })
    }

    pub fn fixed_bank(&self) -> &AngularSpectrumExposure<'b, 's, 'm> {
        self.fixed
    }
    pub fn control_banks(&self) -> &[AngularSpectrumExposure<'b, 's, 'm>] {
        self.controls
    }
    pub fn aperture(&self) -> AngularSpectrumAperture {
        self.aperture
    }
    pub fn window(&self) -> HarmonicWindow {
        self.window
    }
    pub fn distance_from_face_m(&self) -> f64 {
        self.distance
    }
    pub fn frequencies_rad_s(&self) -> &[f64] {
        &self.action.frequencies
    }
    pub fn workspace_bytes(&self) -> usize {
        self.workspace_bytes
    }
    pub fn source_field_entries(&self) -> usize {
        self.source_field_entries
    }
    pub fn kernel_pairs(&self) -> usize {
        self.kernel_pairs
    }

    /// Allocate nothing; return energy and its complete original-order source
    /// derivative. On refusal no result is published and the owner is reusable.
    pub fn value_gradient(
        &mut self,
        coefficients: &[C64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AngularSpectrumApertureEnergyGradient<'_>, AngularSpectrumError> {
        let energy = self.action.evaluate(coefficients, true, &mut cancelled)?;
        Ok(AngularSpectrumApertureEnergyGradient {
            mean_field_energy_density_j_m3: energy,
            source_gradient: &self.action.gradient,
        })
    }

    /// Matrix-free action of the real Cartesian Hessian on a complex source
    /// direction. Fixed illumination does not enter this action. In general
    /// H(i*v) != i*H(v), because an actual finite shutter retains carrier sums.
    pub fn hessian_vector(
        &mut self,
        direction: &[C64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<&[C64], AngularSpectrumError> {
        self.action.evaluate(direction, false, &mut cancelled)?;
        Ok(&self.action.gradient)
    }
}

fn add(a: usize, b: usize) -> Result<usize, AngularSpectrumError> {
    a.checked_add(b)
        .ok_or_else(|| bad("optical control size overflow"))
}
fn mul(a: usize, b: usize) -> Result<usize, AngularSpectrumError> {
    a.checked_mul(b)
        .ok_or_else(|| bad("optical control size overflow"))
}
fn reserve<T>(count: usize) -> Result<Vec<T>, AngularSpectrumError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(count)
        .map_err(|_| bad("optical control allocation failed"))?;
    Ok(values)
}
fn filled<T: Clone>(count: usize, value: T) -> Result<Vec<T>, AngularSpectrumError> {
    let mut values = reserve(count)?;
    values.resize(count, value);
    Ok(values)
}

#[cfg(test)]
#[path = "control/tests.rs"]
mod tests;
