//! Coherent diffracting fields and physical dielectric heat INSIDE a film.
//!
//! A selected layer keeps two original E/H responses per spatial mode. Depth
//! planes reuse the prepared beam's FFTs; the coherent E field is synthesized
//! BEFORE local heat is formed. A sum of mode intensities would lose spatial
//! interference and cannot represent the requested absorption distribution.

use super::modes::fields_from_tangential;
use super::*;
use crate::plane_wave::layers::{FIELD_PREPARATION_WORKSPACE_BYTES, PreparedLayerField};

// Unsupported zero source bins retain their FFT positions without inventing
// an incidence angle or evaluating an unused constitutive response.
pub(super) type PreparedMode<'s, 'm> = Option<[PreparedLayerField<'s, 'm>; 2]>;

/// Preparation work is bounded by two polarizations of every mode/source layer.
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
pub type AngularSpectrumLayerField = AngularSpectrumAbsorptionField;

/// Retains only one selected layer, O(grid points) storage rather than the
/// product of all layers and grid points. Each additional depth uses six
/// existing FFTs and linear modal work, without another film/material solve.
pub struct PreparedAngularSpectrumLayer<'b, 's, 'm> {
    beam: &'b PlanarAngularSpectrum<'s, 'm>,
    layer_index: usize,
    modes: Vec<PreparedMode<'s, 'm>>,
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
        let per_point = core::mem::size_of::<PreparedMode<'_, '_>>() + core::mem::size_of::<f64>();
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
            if mode.is_dormant() {
                modes.push(None);
                continue;
            }
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
            modes.push(Some([te, tm]));
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

    /// Cycle-mean watts absorbed between two depths in this original film,
    /// integrated continuously over the beam's complete periodic cell.
    ///
    /// The native positive film-loss operator integrates each depth interval
    /// exactly, including attenuation and counterpropagating interference.
    /// Distinct transverse Fourier modes are orthogonal on the full cell;
    /// TE and TM electric fields are orthogonal within each mode. Their powers
    /// therefore add only AFTER each mode's coherent interior loss is formed.
    /// This is not a local sample or a finite transverse-aperture average.
    ///
    /// Depths are measured from the selected layer's first face and must obey
    /// `0 <= start <= end <= thickness`. Equal depths absorb exactly zero.
    /// The query reuses the prepared material/field responses, takes linear
    /// work in the admitted grid size, and allocates no FFT or depth buffer.
    pub fn absorbed_power_in_depth_interval_w(
        &self,
        start_depth_m: f64,
        end_depth_m: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<f64, AngularSpectrumError> {
        poll(&mut cancelled)?;
        let thickness_m = self
            .beam
            .source
            .layer(self.layer_index)
            .unwrap()
            .thickness_m;
        if !start_depth_m.is_finite()
            || !end_depth_m.is_finite()
            || start_depth_m < 0.0
            || end_depth_m < start_depth_m
            || end_depth_m > thickness_m
        {
            return Err(bad(
                "depth interval must be ordered and inside the original film",
            ));
        }
        let grid = self.beam.grid;
        let area =
            (grid.width as f64 * grid.spacing_m[0]) * (grid.height as f64 * grid.spacing_m[1]);
        let mut sum = fs_math::dd::Dd::ZERO;
        for (index, (mode, prepared)) in self.beam.modes.iter().zip(&self.modes).enumerate() {
            if index % 64 == 0 {
                poll(&mut cancelled)?;
            }
            let Some(prepared) = prepared else {
                continue;
            };
            let amplitudes = polarized(mode.incident_e, basis(mode.kx, mode.ky));
            for p in 0..2 {
                let absorbed = prepared[p]
                    .absorbed_power_in_depth_interval_w_m2(start_depth_m, end_depth_m)?;
                sum = sum + fs_math::dd::Dd::from_f64(area * (amplitudes[p].norm_sq() * absorbed));
            }
        }
        let power = sum.to_f64();
        if !power.is_finite() || power < 0.0 {
            return Err(bad(
                "selected layer depth-interval absorption is unrepresentable",
            ));
        }
        poll(&mut cancelled)?;
        Ok(power)
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
        self.field_for_beam(self.beam, depth_m, &mut cancelled)
    }

    /// Resolve complete fields and local constitutive heat at arbitrary
    /// transverse positions and depths inside this prepared original film.
    /// Results retain input order. The cached native TE/TM film fields supply
    /// each mode; no material law, film solve or FFT interpolation is repeated.
    pub fn probe_fields(
        &self,
        observations: &[AngularSpectrumLayerObservation],
        limits: AngularSpectrumProbeLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<AngularSpectrumLayerProbe>, AngularSpectrumError> {
        poll(&mut cancelled)?;
        if observations.len() > limits.maximum_points {
            return Err(bad("probe count exceeds its budget"));
        }
        let thickness = self
            .beam
            .source
            .layer(self.layer_index)
            .unwrap()
            .thickness_m;
        for (index, observation) in observations.iter().enumerate() {
            if index % 256 == 0 {
                poll(&mut cancelled)?;
            }
            probes::admit_position(observation.position_m)?;
            if !observation.depth_m.is_finite()
                || observation.depth_m < 0.0
                || observation.depth_m > thickness
            {
                return Err(bad("probe depth is outside the selected layer"));
            }
        }
        let active = probes::admit::<AngularSpectrumLayerProbe>(
            self.beam,
            observations.len(),
            self.workspace_bytes,
            limits,
            &mut cancelled,
        )?;
        let heat_scale = det::sqrt(self.heat_coefficient);
        let mut output = Vec::with_capacity(observations.len());
        for observation in observations {
            let field = probes::synthesize(
                self.beam,
                &active,
                observation.position_m,
                |index| self.mode_field_for_beam(self.beam, index, observation.depth_m),
                &mut cancelled,
            )?;
            let absorbed_power_density_w_m3 =
                absorption::electric_heat_density(&field.electric_v_per_m, heat_scale);
            if !absorbed_power_density_w_m3.is_finite() || absorbed_power_density_w_m3 < 0.0 {
                return Err(bad("local probe heat is unrepresentable"));
            }
            output.push(AngularSpectrumLayerProbe {
                field,
                absorbed_power_density_w_m3,
            });
        }
        poll(&mut cancelled)?;
        Ok(output)
    }

    // A coherence factor changes only the incident amplitudes. The source,
    // spatial bins and their TE/TM boundary responses must remain identical.
    fn admit_beam(&self, beam: &PlanarAngularSpectrum<'s, 'm>) -> Result<(), AngularSpectrumError> {
        let same_source = core::ptr::eq(beam.source(), self.beam.source())
            || (beam.source().angular_frequency_rad_s()
                == self.beam.source().angular_frequency_rad_s()
                && beam.source().same_sources_and_geometry(self.beam.source()));
        if !same_source
            || beam.grid() != self.beam.grid()
            || beam.carrier_rad_m != self.beam.carrier_rad_m
        {
            return Err(bad(
                "layer factors require the same original source, grid and carrier",
            ));
        }
        Ok(())
    }

    pub(super) fn absorbed_power_for_beam(
        &self,
        beam: &PlanarAngularSpectrum<'s, 'm>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<f64, AngularSpectrumError> {
        self.admit_beam(beam)?;
        let area = (beam.grid.width as f64 * beam.grid.spacing_m[0])
            * (beam.grid.height as f64 * beam.grid.spacing_m[1]);
        let mut sum = fs_math::dd::Dd::ZERO;
        for (index, (mode, prepared)) in beam.modes.iter().zip(&self.modes).enumerate() {
            if index % 64 == 0 {
                poll(&mut cancelled)?;
            }
            let Some(prepared) = prepared else {
                continue;
            };
            let amplitudes = polarized(mode.incident_e, basis(mode.kx, mode.ky));
            for p in 0..2 {
                sum = sum
                    + fs_math::dd::Dd::from_f64(
                        area * amplitudes[p].norm_sq() * prepared[p].absorbed_power_w_m2(),
                    );
            }
        }
        let power = sum.to_f64();
        if !power.is_finite() || power < 0.0 {
            return Err(bad("coherence factor film absorption is unrepresentable"));
        }
        Ok(power)
    }

    pub(super) fn modal_field_for_beam(
        &self,
        beam: &PlanarAngularSpectrum<'s, 'm>,
        depth_m: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(Vec<[C64; 3]>, Vec<[C64; 3]>), AngularSpectrumError> {
        poll(&mut cancelled)?;
        self.admit_beam(beam)?;
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
        for index in 0..n {
            if index % 64 == 0 {
                poll(&mut cancelled)?;
            }
            let (e, h) = self.mode_field_for_beam(beam, index, depth_m)?;
            electric.push(e);
            magnetic.push(h);
        }
        Ok((electric, magnetic))
    }

    // FFT planes and nonuniform probes use precisely the same native interior
    // field reconstruction. Their callers admit depth and source identity.
    fn mode_field_for_beam(
        &self,
        beam: &PlanarAngularSpectrum<'s, 'm>,
        index: usize,
        depth_m: f64,
    ) -> Result<([C64; 3], [C64; 3]), AngularSpectrumError> {
        let mode = &beam.modes[index];
        let Some(prepared) = &self.modes[index] else {
            return Ok(([C64::ZERO; 3], [C64::ZERO; 3]));
        };
        let amplitude = polarized(mode.incident_e, basis(mode.kx, mode.ky));
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
        // Tangential E AND H retain exact-critical kz=0 without division.
        fields_from_tangential(mode.kx, mode.ky, self.omega_epsilon, self.omega_mu, et, ht)
            .ok_or_else(|| bad("internal beam field is unrepresentable"))
    }

    pub(super) fn field_for_beam(
        &self,
        beam: &PlanarAngularSpectrum<'s, 'm>,
        depth_m: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AngularSpectrumLayerField, AngularSpectrumError> {
        let (electric, magnetic) = self.modal_field_for_beam(beam, depth_m, &mut cancelled)?;
        let field = beam.synthesize(electric, magnetic, &mut cancelled)?;
        AngularSpectrumLayerField::from_field(field, self.heat_coefficient, &mut cancelled)
    }
}
