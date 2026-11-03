//! Simultaneous coherent vector beams at both faces of the original films.
//!
//! The left beam owns the existing periodic grid, Fourier plans and source.
//! Right-face samples share its frequency, transverse coordinates and phase
//! clock. Each Fourier bin uses the original two-port illumination solve;
//! amplitudes interfere before outgoing power, interior fields or heat form.
//! Both exterior light cones must contain the entire grid. No bin is removed.
//!
//! The cache retains both native source columns for every bin and film. This
//! costs O(points * layers) storage once; subsequent exterior/interior planes
//! require only O(points) field work and the same six existing Fourier passes.

use super::*;
use crate::plane_wave::layers::{
    FIELD_PREPARATION_WORKSPACE_BYTES, MAX_LAYERS, PreparedLayerIllumination, PreparedLayerPoint,
};
use fs_math::dd::Dd;

/// Explicit additional preparation work and total logical numerical storage.
/// Storage includes the borrowed beam, both native film-response columns per
/// polarization/bin, conservative preparation scratch, and one returned field
/// with local heat. Input buffers, other retained results and allocator
/// overhead are excluded. FFT work bounds either preparation or one query.
#[derive(Clone, Copy, Debug)]
pub struct AngularSpectrumIlluminationLimits {
    pub maximum_layer_mode_pairs: usize,
    pub maximum_workspace_bytes: usize,
    pub maximum_fft_work: u64,
}

/// Cycle-mean watts over the declared periodic cell. Ports are left then right;
/// both incoming and outgoing entries are positive magnitudes. An outgoing
/// wave contains the coherent sum of reflection and opposite-face transmission.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct AngularSpectrumIlluminationPower {
    pub incoming_w: [f64; 2],
    pub outgoing_w: [f64; 2],
    pub film_absorbed_w: f64,
}
impl AngularSpectrumIlluminationPower {
    pub fn balance_defect_w(self) -> f64 {
        self.incoming_w[0] + self.incoming_w[1]
            - self.outgoing_w[0]
            - self.outgoing_w[1]
            - self.film_absorbed_w
    }
}

/// Distances are outward from the selected physical face, as in the original
/// one-sided beam. Total fields retain incoming/outgoing standing-wave terms.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AngularSpectrumIlluminationSide {
    IncomingLeft,
    OutgoingLeft,
    TotalLeft,
    IncomingRight,
    OutgoingRight,
    TotalRight,
}

struct IlluminatedMode<'s, 'm> {
    prepared: [PreparedLayerIllumination<'s, 'm>; 2],
    // TE/TM major, left/right incoming face minor, on fixed tangential axes.
    input: [[C64; 2]; 2],
    // Left/right outgoing face major, TE/TM minor.
    outgoing: [[C64; 2]; 2],
}

/// An immutable two-sided illumination of the exact original stack and grid.
/// It borrows the left beam and reuses its FFT plans. Neither incoming source
/// is changed. A failure or cancellation publishes no partially prepared owner.
pub struct PreparedAngularSpectrumIllumination<'b, 's, 'm> {
    beam: &'b PlanarAngularSpectrum<'s, 'm>,
    modes: Vec<IlluminatedMode<'s, 'm>>,
    layer_power_w: Vec<f64>,
    power: AngularSpectrumIlluminationPower,
    workspace_bytes: usize,
    layer_mode_pairs: usize,
}

impl<'s, 'm> PlanarAngularSpectrum<'s, 'm> {
    /// Add a coherent beam incident from the right face. Its Ex/Ey samples use
    /// the SAME x/y axes, grid origin, frequency and exp(+i omega t) clock as
    /// this left beam. Samples are peak V/m at the right physical face.
    ///
    /// A lossy, critical or evanescent exterior refuses even for a dark input
    /// bin: the complete declared two-port source space must be propagating.
    pub fn prepare_illumination(
        &self,
        right_tangential_e_v_per_m: &[[C64; 2]],
        limits: AngularSpectrumIlluminationLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedAngularSpectrumIllumination<'_, 's, 'm>, AngularSpectrumError> {
        poll(&mut cancelled)?;
        let n = self.modes.len();
        let layers = self.source.layer_count();
        if right_tangential_e_v_per_m.len() != n {
            return Err(bad("right-face Ex/Ey must cover the original beam grid"));
        }
        let layer_mode_pairs = n
            .checked_mul(2)
            .and_then(|v| v.checked_mul(layers))
            .ok_or_else(|| bad("two-sided beam preparation work overflow"))?;
        // Each native polarization retains loss_modes and interior_inputs,
        // each with four complex values per film. Vec headers are in the mode.
        let per_point = layers
            .checked_mul(4 * core::mem::size_of::<[[C64; 2]; 2]>())
            .and_then(|v| v.checked_add(core::mem::size_of::<IlluminatedMode<'_, '_>>()))
            .and_then(|v| {
                v.checked_add(2 * core::mem::size_of::<C64>() + core::mem::size_of::<f64>())
            })
            .ok_or_else(|| bad("two-sided beam workspace overflow"))?;
        let fixed = 2 * FIELD_PREPARATION_WORKSPACE_BYTES
            + core::mem::size_of::<[Dd; MAX_LAYERS]>()
            + core::mem::size_of::<[f64; MAX_LAYERS]>()
            + 2 * core::mem::size_of::<PreparedLayerPoint<'_, '_, '_>>();
        let workspace_bytes = n
            .checked_mul(per_point)
            .and_then(|v| v.checked_add(self.workspace_bytes))
            .and_then(|v| v.checked_add(fixed))
            .and_then(|v| v.checked_add(layers * core::mem::size_of::<f64>()))
            .ok_or_else(|| bad("two-sided beam workspace overflow"))?;
        if layer_mode_pairs > limits.maximum_layer_mode_pairs
            || workspace_bytes > limits.maximum_workspace_bytes
            || self.field_fft_work > limits.maximum_fft_work
        {
            return Err(bad("two-sided beam exceeds its work or workspace budget"));
        }
        for (index, (sample, mode)) in right_tangential_e_v_per_m
            .iter()
            .zip(&self.modes)
            .enumerate()
        {
            if index % 256 == 0 {
                poll(&mut cancelled)?;
            }
            if !sample.iter().all(|&z| finite(z)) {
                return Err(bad("right-face Ex/Ey must be finite"));
            }
            if mode.exit_kz.im != 0.0
                || mode.exit_kz.re <= 0.0
                || self.source.exit_permittivity_f_m().im != 0.0
            {
                return Err(bad(
                    "two-sided beam requires both complete exterior light cones",
                ));
            }
        }
        let mut ex = Vec::with_capacity(n);
        let mut ey = Vec::with_capacity(n);
        for (index, sample) in right_tangential_e_v_per_m.iter().enumerate() {
            if index % 256 == 0 {
                poll(&mut cancelled)?;
            }
            ex.push(fft_complex(sample[0]));
            ey.push(fft_complex(sample[1]));
        }
        let mut line = vec![FftComplex::default(); self.grid.width.max(self.grid.height)];
        let mut scratch = line.clone();
        for field in [&mut ex, &mut ey] {
            transform(
                self.grid,
                &self.x_fft,
                &self.y_fft,
                field,
                &mut line,
                &mut scratch,
                true,
                &mut cancelled,
            )?;
        }
        let area = (self.grid.width as f64 * self.grid.spacing_m[0])
            * (self.grid.height as f64 * self.grid.spacing_m[1]);
        let mut modes = Vec::with_capacity(n);
        let mut sums = [Dd::ZERO; 5];
        let mut layer_sums = [Dd::ZERO; MAX_LAYERS];
        let mut heat = [0.0; MAX_LAYERS];
        for (index, mode) in self.modes.iter().enumerate() {
            poll(&mut cancelled)?;
            let p = basis(mode.kx, mode.ky);
            let left = polarized(mode.incident_e, p);
            let right = polarized([complex(ex[index]), complex(ey[index])], p);
            let input = [[left[0], right[0]], [left[1], right[1]]];
            let te = self
                .source
                .prepare_illumination(mode.incidence_cosine, Polarization::Te)?;
            poll(&mut cancelled)?;
            let tm = self
                .source
                .prepare_illumination(mode.incidence_cosine, Polarization::Tm)?;
            let prepared = [te, tm];
            let mut outgoing = [[C64::ZERO; 2]; 2];
            for polarization in 0..2 {
                let response =
                    prepared[polarization].evaluate(input[polarization], &mut heat[..layers])?;
                for side in 0..2 {
                    outgoing[side][polarization] = response.outgoing_tangential_e_v_m[side];
                    sums[side] =
                        sums[side] + Dd::from_f64(area * response.power.incoming_w_m2[side]);
                    sums[side + 2] =
                        sums[side + 2] + Dd::from_f64(area * response.power.outgoing_w_m2[side]);
                }
                sums[4] = sums[4] + Dd::from_f64(area * response.power.absorbed_w_m2);
                for (sum, &value) in layer_sums.iter_mut().zip(&heat[..layers]) {
                    *sum = *sum + Dd::from_f64(area * value);
                }
            }
            modes.push(IlluminatedMode {
                prepared,
                input,
                outgoing,
            });
        }
        let values = sums.map(Dd::to_f64);
        let layer_power_w: Vec<_> = layer_sums[..layers].iter().map(|s| s.to_f64()).collect();
        if values
            .iter()
            .chain(&layer_power_w)
            .any(|v| !v.is_finite() || *v < 0.0)
        {
            return Err(bad("two-sided beam power is unrepresentable"));
        }
        poll(&mut cancelled)?;
        Ok(PreparedAngularSpectrumIllumination {
            beam: self,
            modes,
            layer_power_w,
            power: AngularSpectrumIlluminationPower {
                incoming_w: [values[0], values[1]],
                outgoing_w: [values[2], values[3]],
                film_absorbed_w: values[4],
            },
            workspace_bytes,
            layer_mode_pairs,
        })
    }
}

impl<'b, 's, 'm> PreparedAngularSpectrumIllumination<'b, 's, 'm> {
    pub fn left_beam(&self) -> &'b PlanarAngularSpectrum<'s, 'm> {
        self.beam
    }
    pub fn source(&self) -> &'s PlanarDielectricLayers<'m> {
        self.beam.source()
    }
    pub fn grid(&self) -> AngularSpectrumGrid {
        self.beam.grid()
    }
    pub fn power(&self) -> AngularSpectrumIlluminationPower {
        self.power
    }
    pub fn layer_absorbed_power_w(&self, layer: usize) -> Option<f64> {
        self.layer_power_w.get(layer).copied()
    }
    pub fn workspace_bytes(&self) -> usize {
        self.workspace_bytes
    }
    pub fn preparation_layer_mode_pairs(&self) -> usize {
        self.layer_mode_pairs
    }

    /// Complete Maxwell fields outside either face, with physical signed flux.
    /// The same native six-transform synthesis handles every source/side.
    pub fn field(
        &self,
        side: AngularSpectrumIlluminationSide,
        distance_from_face_m: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AngularSpectrumField, AngularSpectrumError> {
        let (electric, magnetic) = self.modal_field(side, distance_from_face_m, &mut cancelled)?;
        self.beam.synthesize(electric, magnetic, &mut cancelled)
    }

    // Finite-window banks combine equal-frequency sources before one spatial
    // synthesis, retaining the original native two-port field operator.
    pub(super) fn modal_field(
        &self,
        side: AngularSpectrumIlluminationSide,
        distance_from_face_m: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(Vec<[C64; 3]>, Vec<[C64; 3]>), AngularSpectrumError> {
        poll(&mut cancelled)?;
        if !distance_from_face_m.is_finite() || distance_from_face_m < 0.0 {
            return Err(bad("field distance must be finite and nonnegative"));
        }
        use AngularSpectrumIlluminationSide as S;
        let mut electric = Vec::with_capacity(self.modes.len());
        let mut magnetic = Vec::with_capacity(self.modes.len());
        for (index, (mode, illumination)) in self.beam.modes.iter().zip(&self.modes).enumerate() {
            if index % 64 == 0 {
                poll(&mut cancelled)?;
            }
            let p = basis(mode.kx, mode.ky);
            let mut e = [C64::ZERO; 3];
            let mut h = [C64::ZERO; 3];
            let mut add =
                |kz: f64, omega_mu: f64, input, distance| -> Result<(), AngularSpectrumError> {
                    let kz = C64::from_re(kz);
                    let phase = propagation(kz, distance)
                        .ok_or_else(|| bad("two-sided propagation phase is unrepresentable"))?;
                    let (em, hm) = fields(mode.kx, mode.ky, kz, omega_mu, input)
                        .ok_or_else(|| bad("two-sided Maxwell field is unrepresentable"))?;
                    for c in 0..3 {
                        e[c] = e[c] + em[c] * phase;
                        h[c] = h[c] + hm[c] * phase;
                    }
                    Ok(())
                };
            if matches!(side, S::IncomingLeft | S::TotalLeft) {
                add(
                    mode.incident_kz,
                    self.beam.incident_omega_mu,
                    mode.incident_e,
                    -distance_from_face_m,
                )?;
            }
            if matches!(side, S::OutgoingLeft | S::TotalLeft) {
                add(
                    -mode.incident_kz,
                    self.beam.incident_omega_mu,
                    cartesian(illumination.outgoing[0], p),
                    -distance_from_face_m,
                )?;
            }
            if matches!(side, S::IncomingRight | S::TotalRight) {
                add(
                    -mode.exit_kz.re,
                    self.beam.exit_omega_mu,
                    cartesian([illumination.input[0][1], illumination.input[1][1]], p),
                    distance_from_face_m,
                )?;
            }
            if matches!(side, S::OutgoingRight | S::TotalRight) {
                add(
                    mode.exit_kz.re,
                    self.beam.exit_omega_mu,
                    cartesian(illumination.outgoing[1], p),
                    distance_from_face_m,
                )?;
            }
            electric.push(e);
            magnetic.push(h);
        }
        Ok((electric, magnetic))
    }

    /// Coherent fields and local dielectric heating inside an original film.
    /// Native two-source columns retain counterpropagating interference; all
    /// spatial/polarization modes combine BEFORE physical local heat is formed.
    pub fn layer_field(
        &self,
        layer_index: usize,
        depth_m: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AngularSpectrumLayerField, AngularSpectrumError> {
        let (electric, magnetic) = self.modal_layer_field(layer_index, depth_m, &mut cancelled)?;
        let omega = self.source().angular_frequency_rad_s();
        let epsilon = self.source().layer_permittivity_f_m(layer_index).unwrap();
        let heat_coefficient = -0.5 * epsilon.scale(omega).im;
        let field = self.beam.synthesize(electric, magnetic, &mut cancelled)?;
        AngularSpectrumLayerField::from_field(field, heat_coefficient, &mut cancelled)
    }

    pub(super) fn modal_layer_field(
        &self,
        layer_index: usize,
        depth_m: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(Vec<[C64; 3]>, Vec<[C64; 3]>), AngularSpectrumError> {
        poll(&mut cancelled)?;
        let layer = self
            .source()
            .layer(layer_index)
            .ok_or_else(|| bad("selected field layer is outside the original stack"))?;
        if layer.thickness_m == 0.0
            || !depth_m.is_finite()
            || depth_m < 0.0
            || depth_m > layer.thickness_m
        {
            return Err(bad(
                "interior depth needs a positive film and lies outside it",
            ));
        }
        let omega = self.source().angular_frequency_rad_s();
        let epsilon = self.source().layer_permittivity_f_m(layer_index).unwrap();
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
        let mut electric = Vec::with_capacity(self.modes.len());
        let mut magnetic = Vec::with_capacity(self.modes.len());
        for (index, (mode, illumination)) in self.beam.modes.iter().zip(&self.modes).enumerate() {
            if index % 64 == 0 {
                poll(&mut cancelled)?;
            }
            let te = illumination.prepared[0]
                .prepare_point(layer_index, depth_m)?
                .evaluate(illumination.input[0])?;
            let tm = illumination.prepared[1]
                .prepare_point(layer_index, depth_m)?
                .evaluate(illumination.input[1])?;
            let (e, h) = fields_from_tangential(
                mode.kx,
                mode.ky,
                omega_epsilon,
                omega_mu,
                [te.electric_v_m[1], tm.electric_v_m[0]],
                [-te.magnetic_a_m[0], tm.magnetic_a_m[1]],
            )
            .ok_or_else(|| bad("two-sided interior Maxwell field is unrepresentable"))?;
            electric.push(e);
            magnetic.push(h);
        }
        Ok((electric, magnetic))
    }

    /// Net cycle-mean force [N] on the stack over the whole periodic cell.
    /// Uses each original lossless exterior's Maxwell stress and positive
    /// film-loss modes, including transmitted refraction between unlike media.
    /// No FFT, allocation, striction or mechanical-work credit is added.
    /// The retained native loss contractions cost O(points * layers).
    pub fn radiation_force(
        &self,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<[f64; 3], AngularSpectrumError> {
        poll(&mut cancelled)?;
        self.source().admit_refracting_momentum_exterior()?;
        let grid = self.grid();
        let area =
            (grid.width as f64 * grid.spacing_m[0]) * (grid.height as f64 * grid.spacing_m[1]);
        let mut sums = [Dd::ZERO; 3];
        for (mode, illumination) in self.beam.modes.iter().zip(&self.modes) {
            poll(&mut cancelled)?;
            let p = basis(mode.kx, mode.ky);
            for polarization in 0..2 {
                let pressure = illumination.prepared[polarization]
                    .radiation_pressure(illumination.input[polarization])?;
                for (sum, force) in
                    sums.iter_mut()
                        .zip([pressure[0] * p[0], pressure[0] * p[1], pressure[2]])
                {
                    *sum = *sum + Dd::from_f64(area * force);
                }
            }
        }
        let force = sums.map(Dd::to_f64);
        if force.iter().any(|v| !v.is_finite()) {
            return Err(bad("two-sided beam force is unrepresentable"));
        }
        poll(&mut cancelled)?;
        Ok(force)
    }
}
