//! Fields inside a selected original film, including critical and opaque limits.
//!
//! The same bounded slab/state recursion supplies both boundary traces. The
//! interior uses waves decaying from opposite faces, or an entire E/H
//! propagator near kz=0. No division by a vanishing normal wave number and no
//! growing evanescent transfer matrix is needed at a critical film.
//! Depth-cell heat reuses the same positive slab loss modes between its two
//! field traces. It needs neither depth quadrature nor a flux subtraction.

use super::*;

mod kernel;
use kernel::LayerFieldKernel;

pub(in crate::plane_wave) const PREPARATION_WORKSPACE_BYTES: usize =
    core::mem::size_of::<[Slab; MAX_LAYERS]>()
        + core::mem::size_of::<[[C64; 2]; MAX_LAYERS + 1]>()
        + core::mem::size_of::<[C64; MAX_LAYERS]>();

/// Total peak tangential fields for ONE V/m of incident tangential electric
/// field. Multiply both by the actual complex incident amplitude. The H axis
/// is power-conjugate: -H_x for TE (E_y), H_y for TM (E_x), with +z into the
/// stack. Normal flux is Re(E conj(H))/2, including standing-wave terms.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LayerTangentialField {
    pub electric_v_per_m: C64,
    pub conjugate_magnetic_a_per_m: C64,
}

/// One source-bound selected film and polarization. Preparation is O(layers)
/// with fixed stack storage; repeated depth queries allocate nothing and do
/// not reevaluate a material law or solve the stack again.
pub struct PreparedLayerField<'s, 'm> {
    source: &'s PlanarDielectricLayers<'m>,
    layer_index: usize,
    polarization: Polarization,
    incidence_cosine: f64,
    admittance_s: f64,
    kernel: LayerFieldKernel,
    absorbed_power_w_m2: f64,
}

impl<'m> PlanarDielectricLayers<'m> {
    /// Prepare the field of a unit tangential-E plane wave incident from the
    /// first face, using the original exp(+i omega t) convention. Lossy and
    /// evanescent exits remain admitted; zero-thickness films are identities.
    pub fn prepare_layer_field(
        &self,
        layer_index: usize,
        incidence_cosine: f64,
        polarization: Polarization,
    ) -> Result<PreparedLayerField<'_, 'm>, ReflectionError> {
        let selected = self
            .layers
            .get(layer_index)
            .ok_or_else(|| bad("selected field layer is outside the original stack"))?;
        let c = positive(incidence_cosine)?;
        if c > 1.0 {
            return Err(bad("incidence cosine must lie in (0,1]"));
        }
        let epsilon_i = self.incident_epsilon_f_m;
        let mu_i = self.incident_mu_h_m;
        let normal_admittance = positive(det::sqrt(epsilon_i) / det::sqrt(mu_i))?;
        let admittance = positive(match polarization {
            Polarization::Te => normal_admittance * c,
            Polarization::Tm => normal_admittance / c,
        })?;
        let mut slabs = [Slab::IDENTITY; MAX_LAYERS];
        let mut states = [[C64::ZERO; 2]; MAX_LAYERS + 1];
        let mut gains = [C64::ONE; MAX_LAYERS];
        let n = self.layers.len();
        states[n] = admittance_state(self.exit_wave, c, polarization)?;
        for i in (0..n).rev() {
            let layer = self.layers[i];
            slabs[i] = Slab::new(
                layer.wave,
                self.wave_number,
                c,
                polarization,
                layer.source.thickness_m,
            )?;
            (states[i], gains[i]) = if layer.source.thickness_m == 0.0 {
                (states[i + 1], C64::ONE)
            } else {
                propagate_state(slabs[i], states[i + 1])?
            };
        }
        let mut incoming = C64::ONE;
        for &gain in &gains[..layer_index] {
            incoming = incoming * gain;
        }
        let outgoing = incoming * gains[layer_index];
        let trace = |state: [C64; 2], amplitude: C64| {
            let sum = state[0] + state[1];
            [
                amplitude * (state[0] / sum).scale(2.0),
                (amplitude * (state[1] / sum)).scale(2.0 * admittance),
            ]
        };
        let left = trace(states[layer_index], incoming);
        let right = trace(states[layer_index + 1], outgoing);
        let tail = states[layer_index + 1];
        let opposite = ((tail[0] - tail[1]) / (tail[0] + tail[1])) * outgoing;
        let slab = slabs[layer_index];
        // Keep the original positive loss modes. Tiny layer heat is never a
        // difference of nearly equal boundary powers.
        let even = (incoming + opposite).scale(det::sqrt(0.5 * slab.even_loss));
        let odd = (incoming - opposite).scale(det::sqrt(0.5 * slab.odd_loss));
        let absorbed_power_w_m2 = 0.5 * admittance * (even.norm_sq() + odd.norm_sq());
        let omega_mu = positive(
            self.omega * linear_magnetic_permeability(selected.source.response.material())?,
        )?;
        let omega_epsilon = selected.wave.epsilon_f_m.scale(self.omega);
        let normal = selected.wave.normal(c)?;
        let kz = normal.scale(self.wave_number);
        // Match the slab's full q*k*d phase instead of rounding q*k before
        // reconstructing an interior field from its corrected face traces.
        let phase_wavenumber = Dd::from_f64(normal.re) * Dd::from_f64(self.wave_number);
        let (a, b) = match polarization {
            Polarization::Te => (C64::from_re(omega_mu), (kz * kz).scale(1.0 / omega_mu)),
            Polarization::Tm => ((kz * kz) / omega_epsilon, omega_epsilon),
        };
        let kernel = LayerFieldKernel::new(
            selected.source.thickness_m,
            kz,
            phase_wavenumber,
            a,
            b,
            left,
            right,
        )
        .ok_or_else(|| bad("selected layer field is unrepresentable"))?;
        if !absorbed_power_w_m2.is_finite() || absorbed_power_w_m2 < 0.0 {
            return Err(bad("selected layer absorbed power is unrepresentable"));
        }
        Ok(PreparedLayerField {
            source: self,
            layer_index,
            polarization,
            incidence_cosine: c,
            admittance_s: admittance,
            kernel,
            absorbed_power_w_m2,
        })
    }
}

impl<'s, 'm> PreparedLayerField<'s, 'm> {
    pub fn source(&self) -> &'s PlanarDielectricLayers<'m> {
        self.source
    }
    pub fn layer_index(&self) -> usize {
        self.layer_index
    }
    pub fn polarization(&self) -> Polarization {
        self.polarization
    }
    /// Physical heat integrated through this film for unit incident peak E.
    pub fn absorbed_power_w_m2(&self) -> f64 {
        self.absorbed_power_w_m2
    }
    /// Physical cycle-mean heat in a depth interval [W/m²] for unit incident
    /// peak tangential E. Both depths are measured from this film's left face;
    /// 0 <= start <= end <= thickness. An empty interval has exactly zero heat.
    ///
    /// The original slab's positive dissipative modes integrate the local
    /// dielectric loss exactly, including standing waves and longitudinal TM
    /// fields. No depth mesh, flux subtraction, new material evaluation or
    /// stack solve is needed. The query is allocation-free and O(1). Multiply
    /// by the actual incident peak amplitude squared. A thermal consumer owns
    /// the area/time conversion and deposits this energy only once.
    pub fn absorbed_power_in_depth_interval_w_m2(
        &self,
        start_depth_m: f64,
        end_depth_m: f64,
    ) -> Result<f64, ReflectionError> {
        let layer = self.source.layers[self.layer_index];
        let width = admit_depth_interval(layer.source.thickness_m, start_depth_m, end_depth_m)?;
        if width == 0.0 || layer.wave.epsilon_f_m.im == 0.0 {
            return Ok(0.0);
        }
        if start_depth_m == 0.0 && end_depth_m == layer.source.thickness_m {
            return Ok(self.absorbed_power_w_m2);
        }
        let slab = Slab::new_with_phase_thickness(
            layer.wave,
            self.source.wave_number,
            self.incidence_cosine,
            self.polarization,
            width,
            Dd::from_f64(end_depth_m) - Dd::from_f64(start_depth_m),
        )?;
        let trace = |depth| {
            self.kernel
                .field(depth)
                .ok_or_else(|| bad("depth-interval field trace is unrepresentable"))
        };
        let modes = interval_loss_transfer(
            slab,
            self.admittance_s,
            trace(start_depth_m)?,
            trace(end_depth_m)?,
        )?;
        let heat = modes[0].norm_sq() + modes[1].norm_sq();
        if !heat.is_finite() {
            return Err(bad("depth-interval absorbed power is unrepresentable"));
        }
        Ok(heat)
    }
    /// Depth in metres from this film's incident face, including both faces.
    pub fn field(&self, depth_m: f64) -> Result<LayerTangentialField, ReflectionError> {
        let [electric_v_per_m, conjugate_magnetic_a_per_m] = self
            .kernel
            .field(depth_m)
            .ok_or_else(|| bad("field depth is outside the layer or is unrepresentable"))?;
        Ok(LayerTangentialField {
            electric_v_per_m,
            conjugate_magnetic_a_per_m,
        })
    }
}

pub(super) fn admit_depth_interval(
    thickness_m: f64,
    start_depth_m: f64,
    end_depth_m: f64,
) -> Result<f64, ReflectionError> {
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
    Ok(end_depth_m - start_depth_m)
}

// Restrict the already solved film to these two traces. The artificial cut
// faces use the SAME real reference admittance as every original slab; they
// add no material interface. Contract incoming reference waves through the
// positive slab eigenmodes, preserving weak heat without subtracting fluxes.
pub(super) fn interval_loss_transfer(
    slab: Slab,
    admittance_s: f64,
    left: [C64; 2],
    right: [C64; 2],
) -> Result<[C64; 2], ReflectionError> {
    let incoming = (left[0] + divided(left[1], admittance_s)).scale(0.5);
    let opposite = (right[0] - divided(right[1], admittance_s)).scale(0.5);
    let root_reference = 0.5 * det::sqrt(admittance_s);
    let modes = [
        (incoming + opposite).scale(root_reference * det::sqrt(slab.even_loss)),
        (incoming - opposite).scale(root_reference * det::sqrt(slab.odd_loss)),
    ];
    if modes.iter().any(|&z| !finite(z)) {
        return Err(bad("depth-interval loss transfer is unrepresentable"));
    }
    Ok(modes)
}
