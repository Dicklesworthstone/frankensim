//! Fields inside a selected original film, including critical and opaque limits.
//!
//! The same bounded slab/state recursion supplies both boundary traces. The
//! interior uses waves decaying from opposite faces, or an entire E/H
//! propagator near kz=0. No division by a vanishing normal wave number and no
//! growing evanescent transfer matrix is needed at a critical film.

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
