//! Collimated irradiation through the original incoherent Fresnel slab.
//!
//! Both faces share one explicitly retained exterior dielectric; the original
//! constructor uses index one. Input powers are already projected onto physical
//! slab area; neither a second cosine nor a diffuse angular average is applied.
//! Original layer sources and temperature/band support are retained.

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
    exterior_refractive_index: f64,
    exterior_index_source: Option<Arc<str>>,
    direction_cosine: f64,
    s_power_fraction: f64,
    reflectance: [f64; 2],
    transmittance: f64,
    absorptance: Arc<[[f64; 2]]>,
    layer_fields: Arc<[BeamLayerFields]>,
    partition_residual: [f64; 2],
}

#[derive(Clone, Copy, Debug, PartialEq)]
struct BeamLayerFields {
    direction_cosine: f64,
    // [forward@front, backward@back] per unit [front, back] irradiation.
    endpoints: [[f64; 2]; 2],
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
        Self::prepare(
            layers.into(),
            1.0,
            None,
            direction_cosine,
            s_power_fraction,
            cancelled,
        )
    }

    /// Prepare illumination from the same lossless dielectric on both faces.
    ///
    /// The sourced real exterior index must be >= 1 and support the original
    /// layers' temperature and gray/band scope. `direction_cosine` is measured
    /// IN THAT MEDIUM; incident W/m² is physical flux through the slab face.
    /// Every layer retains its absolute index, Snell angle, polarized Fresnel
    /// reflection and all propagating internal returns. No index-squared
    /// radiance factor belongs in this already power-normalized input.
    ///
    /// A layer outside the propagating cone blocks transport by total internal
    /// reflection. Its ray field and heat are zero. Evanescent absorption and
    /// frustrated transmission through a thin gap require a coherent complex-
    /// index wave model and are not represented by this geometric-optics owner.
    /// Different exterior media and thermal emission in the surrounding medium
    /// are outside this beam-only entry point.
    pub fn new_in_medium(
        layers: Vec<AbsorbingLayer>,
        exterior_refractive_index: f64,
        exterior_index_source: impl Into<String>,
        direction_cosine: f64,
        s_power_fraction: f64,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, ConductionError> {
        let source = exterior_index_source.into();
        if !exterior_refractive_index.is_finite()
            || exterior_refractive_index < 1.0
            || !(exterior_refractive_index * exterior_refractive_index).is_finite()
            || source.trim().is_empty()
        {
            return Err(invalid(
                "need a sourced real exterior index >= 1 with representable n²",
            ));
        }
        Self::prepare(
            layers.into(),
            exterior_refractive_index,
            Some(source.into()),
            direction_cosine,
            s_power_fraction,
            cancelled,
        )
    }

    fn prepare(
        layers: Arc<[AbsorbingLayer]>,
        exterior_refractive_index: f64,
        exterior_index_source: Option<Arc<str>>,
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
        let mut indices = [exterior_refractive_index; MAX_SLAB_LAYERS + 2];
        let mut mu = [direction_cosine; MAX_SLAB_LAYERS + 2];
        let mut t = [0.0; MAX_SLAB_LAYERS];
        let mut a = [0.0; MAX_SLAB_LAYERS];
        for (i, layer) in layers.iter().enumerate() {
            let index = layer.refractive_index();
            indices[i + 1] = index;
            // n*sin(theta) is conserved across the whole stack. A matched
            // layer retains even a cosine whose square would underflow.
            mu[i + 1] = if index == exterior_refractive_index || direction_cosine == 1.0 {
                direction_cosine
            } else if index < exterior_refractive_index {
                let ratio = index / exterior_refractive_index;
                let sin_squared = (1.0 - direction_cosine) * (1.0 + direction_cosine);
                let discriminant = ratio.mul_add(ratio, -sin_squared);
                if discriminant > 0.0 {
                    discriminant.sqrt() / ratio
                } else {
                    0.0
                }
            } else {
                (((index - exterior_refractive_index) / index)
                    * ((index + exterior_refractive_index) / index)
                    + (direction_cosine * (exterior_refractive_index / index)).powi(2))
                .sqrt()
            };
            if mu[i + 1] == 0.0 {
                // Unavailable propagating channel, not an opaque absorbing
                // material. Neighboring Fresnel interfaces reflect all power.
                continue;
            }
            let depth = optical_path(layer.absorption_m_inv(), layer.thickness_m(), mu[i + 1]);
            t[i] = remaining_fraction(depth);
            a[i] = decayed_fraction(depth);
        }
        let channels = 3 + 2 * n;
        let mut output = [0.0; 3 + 2 * MAX_SLAB_LAYERS];
        let mut correction = output;
        let mut fields = [[[0.0; 2]; 2]; MAX_SLAB_LAYERS];
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
                if mu[i] == 0.0 || mu[i + 1] == 0.0 {
                    (r[i], e[i]) = (1.0, 0.0);
                } else if indices[i] != indices[i + 1] {
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
                Some(&mut fields[..n]),
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
            || fields[..n]
                .iter()
                .flatten()
                .flatten()
                .any(|x| !nonnegative(*x))
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
            exterior_refractive_index,
            exterior_index_source,
            direction_cosine,
            s_power_fraction,
            reflectance,
            transmittance,
            absorptance,
            layer_fields: (0..n)
                .map(|i| BeamLayerFields {
                    direction_cosine: mu[i + 1],
                    endpoints: fields[i],
                })
                .collect(),
            partition_residual,
        })
    }

    pub fn layers(&self) -> &[AbsorbingLayer] {
        &self.layers
    }
    /// Absolute real index of the common surrounding dielectric.
    pub const fn exterior_refractive_index(&self) -> f64 {
        self.exterior_refractive_index
    }
    /// Original exterior optical source; absent only for the index-one constructor.
    pub fn exterior_index_source(&self) -> Option<&str> {
        self.exterior_index_source.as_deref()
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

    /// Original Snell-adjusted direction cosine inside a physical layer.
    /// Zero denotes an unavailable propagating channel under total reflection.
    pub fn layer_direction_cosine(&self, layer: usize) -> Option<f64> {
        self.layer_fields
            .get(layer)
            .map(|field| field.direction_cosine)
    }

    /// Nonnegative [forward, backward] directional fluxes [W/m² of slab area]
    /// inside the original layer, at depth measured from its front face.
    /// Forward points toward the slab back; backward points toward its front.
    /// Their difference is signed normal flux. At an interface, the requested
    /// layer selects the material-side values; no interface is averaged.
    ///
    /// Both incident fluxes are already projected onto slab area and mutually
    /// incoherent. All s/p-weighted internal returns are cached. Evaluation is
    /// O(1), without allocation, quadrature, a ray solve or thermal emission.
    pub fn layer_fluxes_at_depth(
        &self,
        layer: usize,
        depth_m: f64,
        irradiation_w_m2: [f64; 2],
    ) -> Result<[f64; 2], ConductionError> {
        let (material, field) = self.admit_depth_interval(layer, [depth_m; 2], irradiation_w_m2)?;
        let distances = [depth_m, material.thickness_m() - depth_m];
        let mut result = [0.0; 2];
        for direction in 0..2 {
            let attenuation = remaining_fraction(optical_path(
                material.absorption_m_inv(),
                distances[direction],
                field.direction_cosine,
            ));
            result[direction] = finite(
                field.endpoints[direction][0] * attenuation * irradiation_w_m2[0]
                    + field.endpoints[direction][1] * attenuation * irradiation_w_m2[1],
            )?;
        }
        Ok(result)
    }

    /// Positive beam power absorbed in the depth interval [lo, hi] [W/m²],
    /// measured inside one original layer. Dividing by positive hi-lo gives
    /// cell-average volumetric heating [W/m³]. A zero-width interval is zero.
    ///
    /// Exact Beer interval integration uses expm1 and positive endpoint fields,
    /// never a difference of nearly equal incoming/outgoing fluxes. Narrow and
    /// weakly absorbing cells retain heat without subdividing the optical stack.
    /// Preparation and query complexity are unchanged from the point query.
    pub fn layer_absorbed_flux(
        &self,
        layer: usize,
        interval_m: [f64; 2],
        irradiation_w_m2: [f64; 2],
    ) -> Result<f64, ConductionError> {
        let (material, field) = self.admit_depth_interval(layer, interval_m, irradiation_w_m2)?;
        let fractions = interval_absorptance(material, field, interval_m);
        finite(fractions[0] * irradiation_w_m2[0] + fractions[1] * irradiation_w_m2[1])
    }

    // Share exact positive interval integration with the diffuse preparation.
    pub(super) fn layer_interval_absorptance(
        &self,
        layer: usize,
        interval_m: [f64; 2],
    ) -> Result<[f64; 2], ConductionError> {
        let (material, field) = self.admit_depth_interval(layer, interval_m, [0.0; 2])?;
        let fractions = interval_absorptance(material, field, interval_m);
        Ok([finite(fractions[0])?, finite(fractions[1])?])
    }

    fn admit_depth_interval(
        &self,
        layer: usize,
        interval: [f64; 2],
        irradiation: [f64; 2],
    ) -> Result<(&AbsorbingLayer, &BeamLayerFields), ConductionError> {
        let material = self
            .layers
            .get(layer)
            .ok_or_else(|| invalid("beam depth query needs an original layer index"))?;
        if interval
            .iter()
            .chain(&irradiation)
            .any(|x| !nonnegative(*x))
            || interval[1] < interval[0]
            || interval[1] > material.thickness_m()
        {
            return Err(invalid(
                "beam depth query needs an ordered interval within the layer and nonnegative finite irradiation",
            ));
        }
        Ok((material, &self.layer_fields[layer]))
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

fn interval_absorptance(
    material: &AbsorbingLayer,
    field: &BeamLayerFields,
    [lo, hi]: [f64; 2],
) -> [f64; 2] {
    let alpha = material.absorption_m_inv();
    let mu = field.direction_cosine;
    let absorbed = decayed_fraction(optical_path(alpha, hi - lo, mu));
    let forward = remaining_fraction(optical_path(alpha, lo, mu));
    let backward = remaining_fraction(optical_path(alpha, material.thickness_m() - hi, mu));
    core::array::from_fn(|side| {
        absorbed * (field.endpoints[0][side] * forward + field.endpoints[1][side] * backward)
    })
}

// Use original finite coefficients, not possibly infinite total optical depth
// times a depth fraction. Divide the smaller factor first to retain weak
// grazing attenuation. Overflow is an opaque optical path, not a material input.
fn optical_path(absorption: f64, distance: f64, mu: f64) -> f64 {
    if absorption == 0.0 || distance == 0.0 {
        0.0
    } else {
        (absorption.min(distance) / mu) * absorption.max(distance)
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
            1.0,
            None,
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
        if self.layers != beam.layers || out.len() != n || beam.exterior_refractive_index != 1.0 {
            return Err(invalid(
                "beam and thermal exchange must share the original ordered layers and index-one exterior",
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
