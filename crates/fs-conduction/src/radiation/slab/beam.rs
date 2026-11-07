//! Collimated irradiation through the original incoherent Fresnel slab.
//!
//! Both exterior media have index one. Input powers are already projected onto
//! physical slab area; neither a second cosine nor a diffuse angular average
//! is applied. Original layer sources and temperature/band support are retained.

use super::{
    AbsorbingLayer, AbsorbingSlab, Arc, ConductionError, MAX_SLAB_LAYERS, add, checkpoint,
    decayed_fraction, finite, invalid, layered, nonnegative, refraction, remaining_fraction,
};

/// Prepared oblique beam response for the original ordered material layers.
///
/// The s fraction is incident POWER in the plane-parallel incidence frame.
/// Phase between s and p cannot change absorption in these isotropic layers.
/// Counter-propagating beams are mutually incoherent. Smooth real-index
/// interfaces, Beer absorption and every internal return use the same positive
/// elimination as the diffuse thermal owner. No quadrature is needed.
#[derive(Clone, Debug, PartialEq)]
pub struct SlabBeamResponse {
    layers: Arc<[AbsorbingLayer]>,
    direction_cosine: f64,
    s_power_fraction: f64,
    reflectance: [f64; 2],
    transmittance: f64,
    absorptance: Arc<[[f64; 2]]>,
    partition_residual: [f64; 2],
}

/// Outgoing beam power per slab area, with no thermal-emission contribution.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SlabBeamExchange {
    /// Reflected power leaving [front, back] [W/m²].
    pub reflected_w_m2: [f64; 2],
    /// Transmitted power leaving [front, back] [W/m²], from the opposite face.
    pub transmitted_w_m2: [f64; 2],
    /// Total positive absorbed power assigned to the original layers [W/m²].
    pub absorbed_w_m2: f64,
    /// Signed incoming minus outgoing minus absorbed arithmetic residual [W/m²].
    pub balance_residual_w_m2: f64,
}

impl SlabBeamResponse {
    /// Prepare a beam-only stack in O(N), without preparing diffuse coefficients.
    /// `direction_cosine` is in (0,1]; `s_power_fraction` is in [0,1].
    pub fn new(
        layers: Vec<AbsorbingLayer>,
        direction_cosine: f64,
        s_power_fraction: f64,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, ConductionError> {
        Self::prepare(layers.into(), direction_cosine, s_power_fraction, cancelled)
    }

    fn prepare(
        layers: Arc<[AbsorbingLayer]>,
        direction_cosine: f64,
        s_power_fraction: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, ConductionError> {
        checkpoint(&mut cancelled, 0)?;
        let n = layers.len();
        if n == 0
            || n > MAX_SLAB_LAYERS
            || !direction_cosine.is_finite()
            || direction_cosine <= 0.0
            || direction_cosine > 1.0
            || !s_power_fraction.is_finite()
            || !(0.0..=1.0).contains(&s_power_fraction)
        {
            return Err(invalid(
                "need 1..=64 sourced layers, beam cosine in (0,1] and s power fraction in [0,1]",
            ));
        }
        let mut indices = [1.0; MAX_SLAB_LAYERS + 2];
        let mut mu = [direction_cosine; MAX_SLAB_LAYERS + 2];
        let mut t = [0.0; MAX_SLAB_LAYERS];
        let mut a = [0.0; MAX_SLAB_LAYERS];
        for (i, layer) in layers.iter().enumerate() {
            let index = layer.refractive_index();
            indices[i + 1] = index;
            // Positive Snell discriminant; the matched-index branch also
            // retains direction cosines whose square would underflow.
            mu[i + 1] = if index == 1.0 {
                direction_cosine
            } else {
                (((index - 1.0) / index) * ((index + 1.0) / index)
                    + (direction_cosine / index).powi(2))
                .sqrt()
            };
            let depth = layer.optical_depth() / mu[i + 1];
            t[i] = remaining_fraction(depth);
            a[i] = decayed_fraction(depth);
        }
        let channels = 3 + 2 * n;
        let mut output = [0.0; 3 + 2 * MAX_SLAB_LAYERS];
        let mut correction = output;
        for (pol, weight) in [s_power_fraction, 1.0 - s_power_fraction]
            .into_iter()
            .enumerate()
        {
            checkpoint(&mut cancelled, pol + 1)?;
            if weight == 0.0 {
                continue;
            }
            let mut r = [0.0; MAX_SLAB_LAYERS + 1];
            let mut e = [1.0; MAX_SLAB_LAYERS + 1];
            for i in 0..=n {
                if indices[i] != indices[i + 1] {
                    (r[i], e[i]) = if pol == 0 {
                        refraction::power_fractions(indices[i] * mu[i], indices[i + 1] * mu[i + 1])
                    } else {
                        refraction::power_fractions(indices[i] * mu[i + 1], indices[i + 1] * mu[i])
                    };
                }
            }
            layered::accumulate(
                &t[..n],
                &a[..n],
                &r[..=n],
                &e[..=n],
                true,
                false,
                weight,
                &mut output[..channels],
                &mut correction[..channels],
            );
        }
        let reflectance = [output[channels - 2], output[channels - 1]];
        let transmittance = output[0];
        let absorptance: Arc<[[f64; 2]]> = (0..n)
            .map(|i| [output[1 + 2 * i], output[2 + 2 * i]])
            .collect();
        let mut partition_residual = [0.0; 2];
        let allowance = 256.0 * f64::EPSILON * n as f64;
        if output[..channels]
            .iter()
            .any(|x| !nonnegative(*x) || *x > 1.0 + allowance)
        {
            return Err(invalid("unrepresentable directional slab coefficient"));
        }
        for side in 0..2 {
            let mut residual = -1.0;
            let mut correction = 0.0;
            for coefficient in [reflectance[side], transmittance]
                .into_iter()
                .chain(absorptance.iter().map(|a| a[side]))
            {
                add(&mut residual, &mut correction, coefficient);
            }
            if residual.abs() > allowance {
                return Err(invalid("directional slab partition does not close"));
            }
            partition_residual[side] = residual;
        }
        checkpoint(&mut cancelled, 3)?;
        Ok(Self {
            layers,
            direction_cosine,
            s_power_fraction,
            reflectance,
            transmittance,
            absorptance,
            partition_residual,
        })
    }

    pub fn layers(&self) -> &[AbsorbingLayer] {
        &self.layers
    }
    pub const fn direction_cosine(&self) -> f64 {
        self.direction_cosine
    }
    pub const fn s_power_fraction(&self) -> f64 {
        self.s_power_fraction
    }
    pub const fn reflectances(&self) -> [f64; 2] {
        self.reflectance
    }
    pub const fn transmittance(&self) -> f64 {
        self.transmittance
    }
    /// Layer-major absorption fractions for [front, back] incident beams.
    pub fn layer_absorptance(&self) -> &[[f64; 2]] {
        &self.absorptance
    }
    pub const fn partition_residual(&self) -> [f64; 2] {
        self.partition_residual
    }

    /// Apply nonnegative incident fluxes [front, back] in W/m² of slab area.
    /// All internal returns are already included. No allocation, quadrature or
    /// emission solve occurs. `layer_heat_w_m2` is published only on success.
    pub fn exchange(
        &self,
        irradiation_w_m2: [f64; 2],
        layer_heat_w_m2: &mut [f64],
    ) -> Result<SlabBeamExchange, ConductionError> {
        let n = self.layers.len();
        if layer_heat_w_m2.len() != n || irradiation_w_m2.iter().any(|x| !nonnegative(*x)) {
            return Err(invalid(
                "need nonnegative beam fluxes and one output per original layer",
            ));
        }
        let reflected_w_m2 =
            std::array::from_fn(|side| self.reflectance[side] * irradiation_w_m2[side]);
        let transmitted_w_m2 =
            std::array::from_fn(|side| self.transmittance * irradiation_w_m2[1 - side]);
        let mut heat = [0.0; MAX_SLAB_LAYERS];
        let mut absorbed_w_m2 = 0.0;
        let mut correction = 0.0;
        for i in 0..n {
            heat[i] = finite(
                self.absorptance[i][0] * irradiation_w_m2[0]
                    + self.absorptance[i][1] * irradiation_w_m2[1],
            )?;
            add(&mut absorbed_w_m2, &mut correction, heat[i]);
        }
        let mut residual = 0.0;
        let mut correction = 0.0;
        // Cancel each side's outgoing terms before combining large inputs.
        for side in 0..2 {
            for value in [
                irradiation_w_m2[side],
                -reflected_w_m2[side],
                -transmitted_w_m2[1 - side],
            ] {
                add(&mut residual, &mut correction, finite(value)?);
            }
        }
        add(&mut residual, &mut correction, -finite(absorbed_w_m2)?);
        let balance_residual_w_m2 = finite(residual)?;
        layer_heat_w_m2.copy_from_slice(&heat[..n]);
        Ok(SlabBeamExchange {
            reflected_w_m2,
            transmitted_w_m2,
            absorbed_w_m2,
            balance_residual_w_m2,
        })
    }
}

impl AbsorbingSlab {
    /// Reuse the exact material stack for a prepared polarized incident beam.
    pub fn beam_response(
        &self,
        direction_cosine: f64,
        s_power_fraction: f64,
        cancelled: impl FnMut() -> bool,
    ) -> Result<SlabBeamResponse, ConductionError> {
        SlabBeamResponse::prepare(
            self.layers.clone(),
            direction_cosine,
            s_power_fraction,
            cancelled,
        )
    }

    /// Add actual collimated heating to the original diffuse/emitting layer heat.
    /// The beam must retain exactly this material stack and the same gray/band
    /// scope. Beam outgoing powers stay directional and are returned separately;
    /// they must not silently enter a diffuse view-factor closure.
    pub fn layer_net_fluxes_with_beam(
        &self,
        diffuse_irradiation: [f64; 2],
        blackbody_exitance: &[f64],
        beam: &SlabBeamResponse,
        beam_irradiation: [f64; 2],
        out: &mut [f64],
    ) -> Result<SlabBeamExchange, ConductionError> {
        let n = self.layer_count();
        if self.layers != beam.layers || out.len() != n {
            return Err(invalid(
                "beam and thermal exchange must share the original ordered layers",
            ));
        }
        let mut thermal = [0.0; MAX_SLAB_LAYERS];
        let mut directed = [0.0; MAX_SLAB_LAYERS];
        self.layer_net_fluxes(diffuse_irradiation, blackbody_exitance, &mut thermal[..n])?;
        let result = beam.exchange(beam_irradiation, &mut directed[..n])?;
        for i in 0..n {
            thermal[i] = finite(thermal[i] + directed[i])?;
        }
        out.copy_from_slice(&thermal[..n]);
        Ok(result)
    }
}
