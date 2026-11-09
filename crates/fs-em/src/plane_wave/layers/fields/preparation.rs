//! One native backward solve, shared by selected and whole-stack field owners.
use super::*;

pub(super) struct LayerFieldPreparation<'s, 'm> {
    source: &'s PlanarDielectricLayers<'m>,
    polarization: Polarization,
    incidence_cosine: f64,
    admittance_s: f64,
    slabs: [Slab; MAX_LAYERS],
    states: [[C64; 2]; MAX_LAYERS + 1],
    // Native forward gains are overwritten by their original cumulative input.
    incoming: [C64; MAX_LAYERS + 1],
}

impl<'s, 'm> LayerFieldPreparation<'s, 'm> {
    pub(super) fn new(
        source: &'s PlanarDielectricLayers<'m>,
        incidence_cosine: f64,
        polarization: Polarization,
    ) -> Result<Self, ReflectionError> {
        let c = positive(incidence_cosine)?;
        if c > 1.0 {
            return Err(bad("incidence cosine must lie in (0,1]"));
        }
        let normal_admittance =
            positive(det::sqrt(source.incident_epsilon_f_m) / det::sqrt(source.incident_mu_h_m))?;
        let admittance = positive(match polarization {
            Polarization::Te => normal_admittance * c,
            Polarization::Tm => normal_admittance / c,
        })?;
        let mut slabs = [Slab::IDENTITY; MAX_LAYERS];
        let mut states = [[C64::ZERO; 2]; MAX_LAYERS + 1];
        let mut incoming = [C64::ONE; MAX_LAYERS + 1];
        let n = source.layers.len();
        states[n] = admittance_state(source.exit_wave, c, polarization)?;
        for i in (0..n).rev() {
            let layer = source.layers[i];
            slabs[i] = Slab::new(
                layer.wave,
                source.wave_number,
                c,
                polarization,
                layer.source.thickness_m,
            )?;
            (states[i], incoming[i + 1]) = if layer.source.thickness_m == 0.0 {
                (states[i + 1], C64::ONE)
            } else {
                propagate_state(slabs[i], states[i + 1])?
            };
        }
        for i in 0..n {
            incoming[i + 1] = incoming[i] * incoming[i + 1];
            if !finite(incoming[i + 1]) {
                return Err(bad("one-sided layer field gain is unrepresentable"));
            }
        }
        Ok(Self {
            source,
            polarization,
            incidence_cosine: c,
            admittance_s: admittance,
            slabs,
            states,
            incoming,
        })
    }

    pub(super) fn admittance_s(&self) -> f64 {
        self.admittance_s
    }

    fn trace(&self, face: usize) -> [C64; 2] {
        let state = self.states[face];
        let sum = state[0] + state[1];
        let amplitude = self.incoming[face];
        [
            amplitude * (state[0] / sum).scale(2.0),
            (amplitude * (state[1] / sum)).scale(2.0 * self.admittance_s),
        ]
    }

    pub(super) fn boundary_traces(&self) -> [[C64; 2]; 2] {
        [self.trace(0), self.trace(self.source.layers.len())]
    }

    /// O(1) selection after the original bounded stack sweeps.
    pub(super) fn field(
        &self,
        layer_index: usize,
    ) -> Result<PreparedLayerField<'s, 'm>, ReflectionError> {
        let selected = self
            .source
            .layers
            .get(layer_index)
            .ok_or_else(|| bad("selected field layer is outside the original stack"))?;
        let tail = self.states[layer_index + 1];
        let incoming = self.incoming[layer_index];
        let outgoing = self.incoming[layer_index + 1];
        let opposite = ((tail[0] - tail[1]) / (tail[0] + tail[1])) * outgoing;
        let slab = self.slabs[layer_index];
        let even = loss_amplitude(incoming + opposite, slab.even_loss);
        let odd = loss_amplitude(incoming - opposite, slab.odd_loss);
        let heat = 0.5 * self.admittance_s * (even.norm_sq() + odd.norm_sq());
        if !heat.is_finite() || heat < 0.0 {
            return Err(bad("selected layer absorbed power is unrepresentable"));
        }
        let omega_mu = positive(
            self.source.omega * linear_magnetic_permeability(selected.source.response.material())?,
        )?;
        let omega_epsilon = selected.wave.epsilon_f_m.scale(self.source.omega);
        let normal = selected.wave.normal(self.incidence_cosine)?;
        let kz = normal.scale(self.source.wave_number);
        let phase_wavenumber = Dd::from_f64(normal.re) * Dd::from_f64(self.source.wave_number);
        let (a, b) = match self.polarization {
            Polarization::Te => (C64::from_re(omega_mu), (kz * kz).scale(1.0 / omega_mu)),
            Polarization::Tm => ((kz * kz) / omega_epsilon, omega_epsilon),
        };
        let kernel = LayerFieldKernel::new(
            selected.source.thickness_m,
            kz,
            phase_wavenumber,
            a,
            b,
            self.trace(layer_index),
            self.trace(layer_index + 1),
        )
        .ok_or_else(|| bad("selected layer field is unrepresentable"))?;
        Ok(PreparedLayerField {
            source: self.source,
            layer_index,
            polarization: self.polarization,
            incidence_cosine: self.incidence_cosine,
            admittance_s: self.admittance_s,
            kernel,
            absorbed_power_w_m2: heat,
        })
    }
}
