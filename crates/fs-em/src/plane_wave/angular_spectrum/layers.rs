//! Coherent diffracting fields and physical dielectric heat INSIDE a film.
//!
//! A selected layer keeps two original E/H responses per spatial mode. Depth
//! planes reuse the prepared beam's FFTs; the coherent E field is synthesized
//! BEFORE local heat is formed. A sum of mode intensities would lose spatial
//! interference and cannot represent the requested absorption distribution.

use super::modes::fields_from_tangential;
use super::*;
use crate::plane_wave::layers::{FIELD_PREPARATION_WORKSPACE_BYTES, PreparedLayerField};

/// Preparation visits two polarizations of every mode and every source layer.
/// The byte cap includes the borrowed beam's numerical workspace, the selected
/// layer cache, fixed preparation scratch and one returned E/H/heat plane.
/// Other simultaneously retained planes and allocator overhead are excluded.
#[derive(Clone, Copy, Debug)]
pub struct AngularSpectrumLayerLimits {
    pub maximum_layer_mode_pairs: usize,
    pub maximum_workspace_bytes: usize,
    pub maximum_fft_work: u64,
}

/// A plane inside the selected physical film. Heat is cycle-mean W/m³; its
/// transverse integral is W/m of depth, not the film's total absorbed watts.
#[derive(Debug)]
pub struct AngularSpectrumLayerField {
    field: AngularSpectrumField,
    absorbed_power_density_w_m3: Vec<f64>,
    absorbed_power_per_depth_w_m: f64,
}
impl AngularSpectrumLayerField {
    pub fn field(&self) -> &AngularSpectrumField {
        &self.field
    }
    pub fn absorbed_power_density_w_m3(&self) -> &[f64] {
        &self.absorbed_power_density_w_m3
    }
    pub fn absorbed_power_per_depth_w_m(&self) -> f64 {
        self.absorbed_power_per_depth_w_m
    }
}

/// Retains only one selected layer, O(grid points) storage rather than the
/// product of all layers and grid points. Each additional depth uses six
/// existing FFTs and linear modal work, without another film/material solve.
pub struct PreparedAngularSpectrumLayer<'b, 's, 'm> {
    beam: &'b PlanarAngularSpectrum<'s, 'm>,
    layer_index: usize,
    modes: Vec<[PreparedLayerField<'s, 'm>; 2]>,
    omega_epsilon: C64,
    omega_mu: f64,
    heat_coefficient: f64,
    absorbed_power_w: f64,
    workspace_bytes: usize,
    layer_mode_pairs: usize,
}

impl<'s, 'm> PlanarAngularSpectrum<'s, 'm> {
    /// Cache a film's internal TE/TM boundary response for this actual coherent
    /// beam. The original stack's source and layer ordering remain unchanged.
    pub fn prepare_layer(
        &self,
        layer_index: usize,
        limits: AngularSpectrumLayerLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedAngularSpectrumLayer<'_, 's, 'm>, AngularSpectrumError> {
        poll(&mut cancelled)?;
        let layer = self
            .source
            .layer(layer_index)
            .ok_or_else(|| bad("selected field layer is outside the original stack"))?;
        let n = self.modes.len();
        let layer_mode_pairs = n
            .checked_mul(2)
            .and_then(|v| v.checked_mul(self.source.layer_count()))
            .ok_or_else(|| bad("internal field preparation work overflow"))?;
        let per_point =
            core::mem::size_of::<[PreparedLayerField<'_, '_>; 2]>() + core::mem::size_of::<f64>();
        let workspace_bytes = n
            .checked_mul(per_point)
            .and_then(|v| v.checked_add(self.workspace_bytes))
            .and_then(|v| v.checked_add(FIELD_PREPARATION_WORKSPACE_BYTES))
            .ok_or_else(|| bad("internal field workspace overflow"))?;
        if layer_mode_pairs > limits.maximum_layer_mode_pairs
            || workspace_bytes > limits.maximum_workspace_bytes
            || self.field_fft_work > limits.maximum_fft_work
        {
            return Err(bad("internal field exceeds its work or workspace budget"));
        }
        let omega = self.source.angular_frequency_rad_s();
        let epsilon = self
            .source
            .layer_permittivity_f_m(layer_index)
            .ok_or_else(|| bad("selected field layer is outside the original stack"))?;
        let omega_epsilon = epsilon.scale(omega);
        let omega_mu = omega
            * linear_magnetic_permeability(layer.response.material())
                .map_err(ReflectionError::from)?;
        let heat_coefficient = -0.5 * omega_epsilon.im;
        if !finite(omega_epsilon)
            || omega_epsilon == C64::ZERO
            || !omega_mu.is_finite()
            || omega_mu <= 0.0
            || !heat_coefficient.is_finite()
            || heat_coefficient < 0.0
        {
            return Err(bad(
                "internal constitutive field coefficients are unrepresentable",
            ));
        }
        let mut modes = Vec::with_capacity(n);
        let area = (self.grid.width as f64 * self.grid.spacing_m[0])
            * (self.grid.height as f64 * self.grid.spacing_m[1]);
        let mut absorbed_power_w = 0.0;
        for mode in &self.modes {
            // One source response is capped at 64 original layers.
            poll(&mut cancelled)?;
            let cosine = mode.incidence_cosine;
            let te = self
                .source
                .prepare_layer_field(layer_index, cosine, Polarization::Te)?;
            poll(&mut cancelled)?;
            let tm = self
                .source
                .prepare_layer_field(layer_index, cosine, Polarization::Tm)?;
            let amplitudes = polarized(mode.incident_e, basis(mode.kx, mode.ky));
            absorbed_power_w += area
                * (amplitudes[0].norm_sq() * te.absorbed_power_w_m2()
                    + amplitudes[1].norm_sq() * tm.absorbed_power_w_m2());
            modes.push([te, tm]);
        }
        if !absorbed_power_w.is_finite() || absorbed_power_w < 0.0 {
            return Err(bad("selected layer beam absorption is unrepresentable"));
        }
        poll(&mut cancelled)?;
        Ok(PreparedAngularSpectrumLayer {
            beam: self,
            layer_index,
            modes,
            omega_epsilon,
            omega_mu,
            heat_coefficient,
            absorbed_power_w,
            workspace_bytes,
            layer_mode_pairs,
        })
    }
}

impl<'b, 's, 'm> PreparedAngularSpectrumLayer<'b, 's, 'm> {
    pub fn beam(&self) -> &'b PlanarAngularSpectrum<'s, 'm> {
        self.beam
    }
    pub fn layer_index(&self) -> usize {
        self.layer_index
    }
    /// Physical watts in this film from the original positive loss modes.
    /// No depth quadrature or subtraction of nearly equal fluxes is involved.
    pub fn absorbed_power_w(&self) -> f64 {
        self.absorbed_power_w
    }
    pub fn workspace_bytes(&self) -> usize {
        self.workspace_bytes
    }
    pub fn preparation_layer_mode_pairs(&self) -> usize {
        self.layer_mode_pairs
    }
    /// Full coherent E/H and dielectric heating at depth from the selected
    /// film's first face. At either interface, longitudinal E belongs to THIS
    /// material; tangential E/H and normal D/B obey their own continuity laws.
    pub fn field(
        &self,
        depth_m: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AngularSpectrumLayerField, AngularSpectrumError> {
        poll(&mut cancelled)?;
        let thickness = self
            .beam
            .source
            .layer(self.layer_index)
            .unwrap()
            .thickness_m;
        if !depth_m.is_finite() || depth_m < 0.0 || depth_m > thickness {
            return Err(bad("field depth is outside the selected layer"));
        }
        let n = self.modes.len();
        let mut electric = Vec::with_capacity(n);
        let mut magnetic = Vec::with_capacity(n);
        for (index, (mode, prepared)) in self.beam.modes.iter().zip(&self.modes).enumerate() {
            if index % 64 == 0 {
                poll(&mut cancelled)?;
            }
            let p = basis(mode.kx, mode.ky);
            let amplitude = polarized(mode.incident_e, p);
            let te = prepared[0].field(depth_m)?;
            let tm = prepared[1].field(depth_m)?;
            let et = [
                te.electric_v_per_m * amplitude[0],
                tm.electric_v_per_m * amplitude[1],
            ];
            let ht = [
                te.conjugate_magnetic_a_per_m * amplitude[0],
                tm.conjugate_magnetic_a_per_m * amplitude[1],
            ];
            // Reconstruct from tangential E AND H, so critical kz=0 needs no
            // singular division. Source epsilon and mu remain independent.
            let (e, h) =
                fields_from_tangential(mode.kx, mode.ky, self.omega_epsilon, self.omega_mu, et, ht)
                    .ok_or_else(|| bad("internal beam field is unrepresentable"))?;
            electric.push(e);
            magnetic.push(h);
        }
        let field = self.beam.synthesize(electric, magnetic, &mut cancelled)?;
        let mut absorbed_power_density_w_m3 = Vec::with_capacity(n);
        let mut sum = 0.0;
        for (index, e) in field.electric_v_per_m().iter().enumerate() {
            if index % 256 == 0 {
                poll(&mut cancelled)?;
            }
            let heat = if self.heat_coefficient == 0.0 {
                0.0
            } else {
                self.heat_coefficient * e.iter().map(|v| v.norm_sq()).sum::<f64>()
            };
            if !heat.is_finite() || heat < 0.0 {
                return Err(bad("local dielectric heat is unrepresentable"));
            }
            absorbed_power_density_w_m3.push(heat);
            sum += heat;
        }
        let absorbed_power_per_depth_w_m =
            sum * self.beam.grid.spacing_m[0] * self.beam.grid.spacing_m[1];
        if !absorbed_power_per_depth_w_m.is_finite() {
            return Err(bad("integrated dielectric heat is unrepresentable"));
        }
        poll(&mut cancelled)?;
        Ok(AngularSpectrumLayerField {
            field,
            absorbed_power_density_w_m3,
            absorbed_power_per_depth_w_m,
        })
    }
}
