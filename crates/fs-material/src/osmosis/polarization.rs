//! Steady concentration polarization on both sides of a selective membrane.
//!
//! Each optional film integrates `Js = Jv*c - D*dc/dx` with constant `D`,
//! thickness, and the SAME solvent and solute flux as the membrane. Its
//! permeance is `k = D / thickness` [m/s]. The original membrane integrates
//! `Js = (1-sigma)*Jv*c - P*dc/ds`; composing these three exact solutions gives
//! nonnegative bulk donor coefficients and both membrane-face concentrations.
//! The membrane water law must use those FACE concentrations, so its water
//! flux generally requires a scalar implicit solve owned by the consumer.
//!
//! This is the stagnant-film, isothermal, dilute-solution approximation with
//! negligible film storage and no additional hydraulic pressure drop. The
//! caller supplies species-specific film coefficients and their applicability;
//! no flow geometry, turbulence correlation, salt dissociation, ionic field,
//! or solute boundary condition is inferred. An absent film is an exact
//! well-mixed boundary. A present film always has positive finite permeance.
//! Fully rejected solutes remain exactly sealed; polarization never creates
//! an artificial solute leak. Large Peclet numbers use positive scaled sums
//! and logarithmic reconstruction rather than capped exponents or floors.
//!
//! Primary film formulation: Biesheuvel et al. (2023), equations 1--3,
//! <https://arxiv.org/html/2308.14233v2>. Feed/draw polarization and its effect
//! on osmotic driving force: McCutcheon and Elimelech (2006),
//! <https://doi.org/10.1016/j.memsci.2006.07.049>. The membrane itself remains
//! the original Spiegler--Kedem owner in [`super::MembraneTransport`].

use super::{MembraneTransport, OsmosisError, Result, diffusive_exchange, finite, positive, shape};
use fs_math::det;

/// Species-specific diffusive film permeance `D / thickness` [m/s].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FilmPermeance(f64);

impl FilmPermeance {
    /// Admit a strictly positive, finite coefficient in m/s.
    /// Use `None` at preparation for an absent, perfectly mixed film.
    pub fn new(value_m_s: f64) -> Result<Self> {
        if positive(value_m_s) {
            Ok(Self(value_m_s))
        } else {
            Err(OsmosisError::Invalid(
                "film permeance must be finite and strictly positive m/s",
            ))
        }
    }

    /// Construct `k = D / thickness` from diffusivity [m²/s] and thickness [m].
    pub fn from_diffusion_layer(diffusivity_m2_s: f64, thickness_m: f64) -> Result<Self> {
        if !positive(diffusivity_m2_s) || !positive(thickness_m) {
            return Err(OsmosisError::Invalid(
                "film diffusivity and thickness must be finite and strictly positive",
            ));
        }
        let value = diffusivity_m2_s / thickness_m;
        if positive(value) {
            Ok(Self(value))
        } else {
            Err(OsmosisError::Arithmetic(
                "diffusion-layer permeance is not representable",
            ))
        }
    }

    /// Numerical coherent-SI value [m/s].
    pub const fn value(self) -> f64 {
        self.0
    }
}

/// Original membrane and one optional external film per species on each side.
#[derive(Clone, Debug, PartialEq)]
pub struct PolarizedMembraneTransport {
    membrane: MembraneTransport,
    left: Box<[Option<FilmPermeance>]>,
    right: Box<[Option<FilmPermeance>]>,
}

impl PolarizedMembraneTransport {
    /// Prepare the three-layer law on the original membrane's species order.
    /// Allocations occur only here; each later coefficient/reconstruction call
    /// uses fixed-size local storage and performs no allocation.
    pub fn new(
        membrane: MembraneTransport,
        left: Vec<Option<FilmPermeance>>,
        right: Vec<Option<FilmPermeance>>,
    ) -> Result<Self> {
        shape(membrane.solutes().len(), left.len())?;
        shape(membrane.solutes().len(), right.len())?;
        Ok(Self {
            membrane,
            left: left.into_boxed_slice(),
            right: right.into_boxed_slice(),
        })
    }

    /// The original membrane; evaluate its water law at reconstructed faces.
    pub const fn membrane(&self) -> &MembraneTransport {
        &self.membrane
    }

    /// Owned boxed preparation payload in bytes, excluding `size_of::<Self>()`.
    /// Includes the original membrane's solute array and both film arrays.
    pub fn heap_bytes(&self) -> usize {
        core::mem::size_of_val(self.membrane.solutes())
            + core::mem::size_of_val(self.left.as_ref())
            + core::mem::size_of_val(self.right.as_ref())
    }

    /// Prepare bulk donor coefficients and reusable positive face weights at
    /// one trial solvent flux [m/s], with positive flux directed left to right.
    ///
    /// The consumer can eliminate its finite solute inventories using the
    /// donor coefficients, then recover membrane-face concentrations from the
    /// SAME prepared trial without re-evaluating the layer exponential fits.
    pub fn coefficients(
        &self,
        index: usize,
        volume_flux_m_s: f64,
    ) -> Result<PolarizedSoluteCoefficients> {
        let &(reflection, permeance) =
            self.membrane
                .solutes()
                .get(index)
                .ok_or(OsmosisError::Invalid(
                    "solute index is outside the prepared polarization layout",
                ))?;
        if !volume_flux_m_s.is_finite() {
            return Err(OsmosisError::Invalid(
                "polarization solvent flux must be finite m/s",
            ));
        }
        let left = self.left[index];
        let right = self.right[index];
        if left.is_none() && right.is_none() {
            let donor = self
                .membrane
                .solute_stream_coefficients(index, volume_flux_m_s)?;
            return Ok(PolarizedSoluteCoefficients::identity_faces(donor));
        }

        let speed = volume_flux_m_s.abs();
        let membrane_speed = (1.0 - reflection.value()) * speed;
        if permeance.value() == 0.0 && membrane_speed == 0.0 {
            // Js=0 disconnects the two film concentration problems. Keeping
            // each exponential as a weight permits tiny bulk concentrations
            // to compensate a large enrichment exponent without overflow.
            return Ok(PolarizedSoluteCoefficients {
                donor: (0.0, 0.0),
                left_face: [
                    Weight::from_log(left.map_or(0.0, |k| volume_flux_m_s / k.0)),
                    Weight::ZERO,
                ],
                right_face: [
                    Weight::ZERO,
                    Weight::from_log(right.map_or(0.0, |k| -volume_flux_m_s / k.0)),
                ],
            });
        }

        // All layers' advective velocities have the same sign. Orient them
        // downstream first, which keeps every attenuation exponent <=0.
        let reverse = volume_flux_m_s < 0.0;
        let (upstream, downstream) = if reverse {
            (right, left)
        } else {
            (left, right)
        };
        let layers = [
            Layer::film(upstream, speed)?,
            Layer::new(permeance.value(), membrane_speed)?,
            Layer::film(downstream, speed)?,
        ];
        let mut result = match ordinary_coefficients(&layers) {
            Some(result) => result,
            None => logarithmic_coefficients(&layers)?,
        };
        // Preserve absent-film concentration identity exactly, including when
        // the other side needs logarithmic reconstruction.
        if upstream.is_none() {
            result.left_face = [Weight::ONE, Weight::ZERO];
        }
        if downstream.is_none() {
            result.right_face = [Weight::ZERO, Weight::ONE];
        }
        if reverse {
            result = result.reversed();
        }
        Ok(result)
    }
}

/// Bulk solute transport and both membrane faces at one fixed water-flux trial.
/// This value has no heap storage and is reusable for different bulk endpoints.
#[derive(Clone, Copy, Debug)]
pub struct PolarizedSoluteCoefficients {
    donor: (f64, f64),
    left_face: [Weight; 2],
    right_face: [Weight; 2],
}

impl PolarizedSoluteCoefficients {
    /// Nonnegative `(k_left,k_right)` [m/s] for
    /// `Js = k_left*c_bulk_left - k_right*c_bulk_right` [mol/(m² s)].
    pub const fn donor_coefficients(&self) -> (f64, f64) {
        self.donor
    }

    /// Recover `(c_membrane_left,c_membrane_right)` [mol/m³] from nonnegative
    /// finite bulk concentrations. These, not bulk concentrations, enter the
    /// original membrane's osmotic water-driving law.
    ///
    /// Positive reconstruction weights avoid cancellation and preserve exact
    /// empty endpoints. An unrepresentable enriched face refuses; neither
    /// concentration floors nor exponent caps modify the constitutive law.
    pub fn membrane_face_concentrations(
        &self,
        bulk_left_mol_m3: f64,
        bulk_right_mol_m3: f64,
    ) -> Result<(f64, f64)> {
        if !super::nonnegative(bulk_left_mol_m3) || !super::nonnegative(bulk_right_mol_m3) {
            return Err(OsmosisError::Invalid(
                "polarization bulk concentrations must be finite and nonnegative",
            ));
        }
        let left = reconstruct(self.left_face, bulk_left_mol_m3, bulk_right_mol_m3)?;
        let right = reconstruct(self.right_face, bulk_left_mol_m3, bulk_right_mol_m3)?;
        Ok((left, right))
    }

    fn identity_faces(donor: (f64, f64)) -> Self {
        Self {
            donor,
            left_face: [Weight::ONE, Weight::ZERO],
            right_face: [Weight::ZERO, Weight::ONE],
        }
    }

    fn reversed(self) -> Self {
        Self {
            donor: (self.donor.1, self.donor.0),
            left_face: [self.right_face[1], self.right_face[0]],
            right_face: [self.left_face[1], self.left_face[0]],
        }
    }
}

// Most weights are ordinary finite doubles. Only weights outside the normal
// range retain their logarithm, allowing weight*bulk to remain representable
// even when the weight alone would overflow or round to zero.
#[derive(Clone, Copy, Debug)]
enum Weight {
    Value(f64),
    Log(f64),
}

impl Weight {
    const ZERO: Self = Self::Value(0.0);
    const ONE: Self = Self::Value(1.0);

    fn from_log(logarithm: f64) -> Self {
        if logarithm == f64::NEG_INFINITY {
            return Self::ZERO;
        }
        let value = det::exp(logarithm);
        if value.is_normal() {
            Self::Value(value)
        } else {
            Self::Log(logarithm)
        }
    }

    fn times(self, concentration: f64) -> Result<f64> {
        if concentration == 0.0 {
            return Ok(0.0);
        }
        let value = match self {
            Self::Value(weight) => weight * concentration,
            Self::Log(logarithm) => det::exp(logarithm + det::ln(concentration)),
        };
        finite(
            value,
            "polarized membrane-face concentration is not representable",
        )
    }
}

fn reconstruct(weights: [Weight; 2], left: f64, right: f64) -> Result<f64> {
    finite(
        weights[0].times(left)? + weights[1].times(right)?,
        "polarized membrane-face concentration is not representable",
    )
}

// Backward relation for one oriented layer: c_in = r*c_out + s*Js,
// r=exp(-Pe), s=1/k_upwind. An absent film has r=1,s=0. Store the
// original Pe as well as rounded coefficients so an underflowed reverse donor
// cannot destroy an otherwise representable composite resistance or face.
#[derive(Clone, Copy)]
struct Layer {
    donor: f64,
    exchange: f64,
    speed: f64,
    peclet: f64,
    absent: bool,
}

impl Layer {
    fn film(film: Option<FilmPermeance>, speed: f64) -> Result<Self> {
        match film {
            Some(k) => Self::new(k.0, speed),
            None => Ok(Self {
                donor: 1.0,
                exchange: 1.0,
                speed: 0.0,
                peclet: 0.0,
                absent: true,
            }),
        }
    }

    fn new(permeance: f64, speed: f64) -> Result<Self> {
        let exchange = diffusive_exchange(permeance, speed)?;
        Ok(Self {
            donor: exchange + speed,
            exchange,
            speed,
            peclet: if speed == 0.0 { 0.0 } else { speed / permeance },
            absent: false,
        })
    }

    fn log_resistance(self) -> f64 {
        if self.absent {
            f64::NEG_INFINITY
        } else if self.donor.is_finite() {
            -det::ln(self.donor)
        } else {
            // The sum of two finite local donor terms may exceed f64::MAX
            // even when the composite membrane law remains representable.
            let large = self.exchange.max(self.speed);
            let small = self.exchange.min(self.speed);
            -(det::ln(large) + det::ln(1.0 + small / large))
        }
    }
}

// Ordinary conditions need only the original layer fits and positive scaled
// arithmetic. The thresholds select a computational path, not a physical
// cap: all remaining valid cases use the same equations in logarithms below.
fn ordinary_coefficients(layers: &[Layer; 3]) -> Option<PolarizedSoluteCoefficients> {
    if layers.iter().map(|layer| layer.peclet).sum::<f64>() > 400.0 {
        return None;
    }
    let mut resistance = [0.0; 3];
    let mut attenuation = [1.0; 3];
    for (index, layer) in layers.iter().enumerate() {
        if !layer.absent {
            resistance[index] = 1.0 / layer.donor;
            attenuation[index] = layer.exchange / layer.donor;
            if !resistance[index].is_normal() || !attenuation[index].is_normal() {
                return None;
            }
        }
    }
    let scale = resistance[0].max(resistance[1]).max(resistance[2]);
    for value in &mut resistance {
        if *value > 0.0 {
            *value /= scale;
            if *value < 1e-100 {
                return None;
            }
        }
    }
    let [left, membrane, right] = resistance;
    let [rl, rm, rr] = attenuation;
    let total = left + rl * membrane + rl * rm * right;
    let upwind = (1.0 / scale) / total;
    if !positive(upwind) {
        return None;
    }
    Some(PolarizedSoluteCoefficients {
        donor: (upwind, upwind * (rl * rm * rr)),
        left_face: [
            Weight::Value((membrane + rm * right) / total),
            Weight::Value((rm * rr * left) / total),
        ],
        right_face: [
            Weight::Value(right / total),
            Weight::Value(rr * (left + rl * membrane) / total),
        ],
    })
}

fn logarithmic_coefficients(layers: &[Layer; 3]) -> Result<PolarizedSoluteCoefficients> {
    let [left, membrane, right] = layers.map(Layer::log_resistance);
    let [pl, pm, pr] = layers.map(|layer| layer.peclet);
    let total = log_sum([left, membrane - pl, right - pl - pm]);
    let upwind = det::exp(-total);
    if !positive(upwind) {
        return Err(OsmosisError::Arithmetic(
            "composite solute donor coefficient is not representable",
        ));
    }
    let backward = finite(
        det::exp(-total - pl - pm - pr),
        "composite reverse donor coefficient is not representable",
    )?;
    Ok(PolarizedSoluteCoefficients {
        donor: (upwind, backward),
        left_face: [
            Weight::from_log(log_sum([membrane, right - pm, f64::NEG_INFINITY]) - total),
            Weight::from_log(left - pm - pr - total),
        ],
        right_face: [
            Weight::from_log(right - total),
            Weight::from_log(log_sum([left, membrane - pl, f64::NEG_INFINITY]) - pr - total),
        ],
    })
}

fn log_sum(values: [f64; 3]) -> f64 {
    let maximum = values[0].max(values[1]).max(values[2]);
    if maximum == f64::NEG_INFINITY {
        return maximum;
    }
    maximum
        + det::ln(
            det::exp(values[0] - maximum)
                + det::exp(values[1] - maximum)
                + det::exp(values[2] - maximum),
        )
}

#[cfg(test)]
mod tests;
