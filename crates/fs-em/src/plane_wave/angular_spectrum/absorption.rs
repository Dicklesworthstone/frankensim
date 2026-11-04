//! Source-bound dielectric/Joule heat in finite films and the exit half-space.
//!
//! The exp(+i omega t) convention gives q = -omega Im(epsilon) |E|²/2.
//! All electric components and coherent spatial interference contribute. The
//! material's permittivity already includes its admitted Ohmic, relaxation or
//! optical loss; conductivity must not be added a second time.

use super::*;

mod kernels;
use kernels::{electric_heat_density, mode_absorbed_power};

/// A coherent electric/magnetic field plane and its physical absorbed heat.
/// Heat density is cycle-mean W/m³. Its transverse integral is W/m of depth,
/// not the total watts deposited throughout a finite or semi-infinite volume.
#[derive(Debug)]
pub struct AngularSpectrumAbsorptionField {
    field: AngularSpectrumField,
    absorbed_power_density_w_m3: Vec<f64>,
    absorbed_power_per_depth_w_m: f64,
}

impl AngularSpectrumAbsorptionField {
    pub fn field(&self) -> &AngularSpectrumField {
        &self.field
    }
    pub fn absorbed_power_density_w_m3(&self) -> &[f64] {
        &self.absorbed_power_density_w_m3
    }
    pub fn absorbed_power_per_depth_w_m(&self) -> f64 {
        self.absorbed_power_per_depth_w_m
    }

    // Both one- and two-sided film illumination, and the transmitted substrate,
    // synthesize the complete coherent E before applying constitutive loss.
    pub(super) fn from_field(
        field: AngularSpectrumField,
        heat_coefficient: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, AngularSpectrumError> {
        let n = field.electric_v_per_m().len();
        let mut absorbed_power_density_w_m3 = Vec::with_capacity(n);
        let mut sum = 0.0;
        let heat_scale = det::sqrt(heat_coefficient);
        for (index, e) in field.electric_v_per_m().iter().enumerate() {
            if index % 256 == 0 {
                poll(&mut cancelled)?;
            }
            let heat = electric_heat_density(e, heat_scale);
            if !heat.is_finite() || heat < 0.0 {
                return Err(bad("local dielectric heat is unrepresentable"));
            }
            absorbed_power_density_w_m3.push(heat);
            sum += heat;
        }
        let absorbed_power_per_depth_w_m = sum * field.grid.spacing_m[0] * field.grid.spacing_m[1];
        if !absorbed_power_per_depth_w_m.is_finite() {
            return Err(bad("integrated dielectric heat is unrepresentable"));
        }
        poll(&mut cancelled)?;
        Ok(Self {
            field,
            absorbed_power_density_w_m3,
            absorbed_power_per_depth_w_m,
        })
    }
}

impl PlanarAngularSpectrum<'_, '_> {
    /// Logical workspace for one transmitted E/H/heat plane, including this
    /// beam's retained modes and FFT workspace. Other retained query results,
    /// source storage and allocator overhead are excluded as in `prepare`.
    pub fn transmitted_absorption_workspace_bytes(&self) -> Result<usize, AngularSpectrumError> {
        self.modes
            .len()
            .checked_mul(core::mem::size_of::<f64>())
            .and_then(|bytes| bytes.checked_add(self.workspace_bytes))
            .ok_or_else(|| bad("transmitted heat workspace overflow"))
    }

    /// Reconstruct deposited heat at a depth measured INTO the exit half-space
    /// from the last film face. Uses the original complex permittivity and the
    /// same cached coherent E/H reconstruction as `field(Transmitted, ...)`.
    /// This also handles Ohmic skin loss, dispersive dielectric/optical loss,
    /// longitudinal TM fields and transmitted evanescent modes. A lossless
    /// evanescent field deposits exactly zero heat despite its spatial decay.
    ///
    /// Admission precedes allocation. Each query requires six existing FFTs
    /// and O(N) heat work, with no additional material or multilayer solve.
    pub fn transmitted_absorption_field(
        &self,
        depth_m: f64,
        limits: AngularSpectrumLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AngularSpectrumAbsorptionField, AngularSpectrumError> {
        poll(&mut cancelled)?;
        if self.modes.len() > limits.maximum_grid_points
            || self.transmitted_absorption_workspace_bytes()? > limits.maximum_workspace_bytes
            || self.field_fft_work > limits.maximum_fft_work
        {
            return Err(bad("transmitted heat exceeds its work or workspace budget"));
        }
        let field = self.field(AngularSpectrumSide::Transmitted, depth_m, &mut cancelled)?;
        AngularSpectrumAbsorptionField::from_field(
            field,
            -0.5 * self.exit_omega_epsilon.im,
            &mut cancelled,
        )
    }

    /// Watts absorbed between two depths in the original transmitted medium,
    /// integrated over the entire periodic transverse cell. The start must be
    /// finite and nonnegative; the end may be +infinity for a substrate tail.
    /// Finite-film absorption remains separately reported by `power()`.
    ///
    /// Cell orthogonality removes cross-mode terms only after the transverse
    /// integral. Each mode deposits P(0) exp(-2 alpha z0) [-expm1(-2 alpha dz)],
    /// alpha = -Im(kz), so even weak absorption survives without subtracting
    /// nearly equal powers. The analytic O(N) integral allocates no workspace
    /// and needs neither depth quadrature nor FFTs. A finite detector cannot
    /// use this whole-cell result to discard coherent spatial interference.
    pub fn transmitted_absorbed_power_w(
        &self,
        depth_start_m: f64,
        depth_end_m: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<f64, AngularSpectrumError> {
        poll(&mut cancelled)?;
        if !depth_start_m.is_finite()
            || depth_start_m < 0.0
            || depth_end_m.is_nan()
            || depth_end_m < depth_start_m
        {
            return Err(bad(
                "absorption depths must be ordered, nonnegative and start finite",
            ));
        }
        if depth_start_m == depth_end_m || self.exit_omega_epsilon.im == 0.0 {
            return Ok(0.0);
        }
        let mut absorbed = 0.0;
        for (index, mode) in self.modes.iter().enumerate() {
            if index % 256 == 0 {
                poll(&mut cancelled)?;
            }
            if mode.transmitted_power_w == 0.0 {
                continue;
            }
            if mode.exit_kz.im >= 0.0 {
                return Err(bad("absorbing substrate attenuation is unrepresentable"));
            }
            absorbed += mode_absorbed_power(
                mode.transmitted_power_w,
                mode.exit_kz.im,
                depth_start_m,
                depth_end_m,
            );
        }
        if !absorbed.is_finite() || absorbed < 0.0 {
            return Err(bad("transmitted absorbed power is unrepresentable"));
        }
        poll(&mut cancelled)?;
        Ok(absorbed)
    }
}
