//! Spatial derivatives of the original nonlinear generated-conduction profile.
//!
//! Fourier's law gives T'=-q/k(T). Differentiating with q'=s gives
//! T''=-(s+k'(T)*T'^2)/k(T), including curvature without internal generation.
//! Original source knots are never smoothed: T and T' remain defined at a
//! continuous-knot crossing, while a nonzero-flux T'' requires equal slopes.

use super::{NonlinearGeneratedLayeredWallSolution, checkpoint, finite, invalid, positive};
use crate::ConductionError;
use crate::layered::PlanarThermalLayer;
use crate::layered::geometry::LayerShape;
use fs_exec::Cx;

/// Temperature and first spatial derivative in one original physical layer.
/// Unlike curvature, these remain defined when continuous k(T) changes slope.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NonlinearGeneratedLayerTemperatureSample {
    pub temperature_k: f64,
    /// Signed left-to-right Fourier flux, W/m².
    pub heat_flux_w_per_m2: f64,
    /// Original dT/dx, K/m.
    pub gradient_k_per_m: f64,
}

/// A temperature sample plus its unique analytic curvature, K/m².
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NonlinearGeneratedLayerTemperatureJet {
    pub sample: NonlinearGeneratedLayerTemperatureSample,
    pub curvature_k_per_m2: f64,
}

/// Immutable borrow of one solved nonlinear planar layer and original source.
/// Queries reuse the original exact segment inversion and its cancellation
/// checks. They allocate nothing, with work bounded by retained source knots;
/// they do not perform another wall solve or interpolate between face values.
#[derive(Clone, Copy, Debug)]
pub struct NonlinearGeneratedLayerTemperatureProfile<'a> {
    solution: &'a NonlinearGeneratedLayeredWallSolution,
    layer: &'a PlanarThermalLayer,
    layer_index: usize,
    source_w_m3: f64,
}

impl NonlinearGeneratedLayeredWallSolution {
    /// Borrow one original planar layer for temperature/gradient/curvature
    /// observations. The private retained Kirchhoff law still controls material
    /// temperature admission; contacts and source knots are not reconstructed.
    pub fn temperature_profile(
        &self,
        cx: &Cx<'_>,
        layer_index: usize,
    ) -> Result<NonlinearGeneratedLayerTemperatureProfile<'_>, ConductionError> {
        checkpoint(cx, layer_index)?;
        let layer = self
            .wall
            .layers
            .get(layer_index)
            .ok_or_else(|| invalid("nonlinear temperature layer index is out of range"))?;
        let row = self
            .layer_temperatures
            .get(layer_index)
            .ok_or_else(|| invalid("nonlinear layer temperature row is missing"))?;
        let source_w_m3 = *self
            .volumetric_generation_w_m3
            .get(layer_index)
            .ok_or_else(|| invalid("nonlinear layer generation row is missing"))?;
        let law = self
            .laws
            .get(layer_index)
            .ok_or_else(|| invalid("nonlinear layer source integral is missing"))?;
        let shape = self
            .shapes
            .get(layer_index)
            .ok_or_else(|| invalid("nonlinear layer geometry is missing"))?;
        if *shape != LayerShape::planar(layer.thickness_m()) {
            return Err(invalid(
                "planar temperature derivatives require original planar geometry",
            ));
        }
        finite("nonlinear volumetric generation", source_w_m3)?;
        let [low, high] = law.interval();
        for temperature in row.face_temperatures_k {
            positive("nonlinear layer face temperature", temperature)?;
            crate::material::TemperatureSpan::Sampled { low, high }.check(temperature)?;
            layer.conductivity().span().check(temperature)?;
        }
        for flux in row.face_fluxes_w_per_m2 {
            finite("nonlinear layer face flux", flux)?;
        }
        Ok(NonlinearGeneratedLayerTemperatureProfile {
            solution: self,
            layer,
            layer_index,
            source_w_m3,
        })
    }

    /// First-order observation; a conductivity slope change does not remove a
    /// well-defined temperature or Fourier gradient.
    pub fn temperature_gradient_at(
        &self,
        cx: &Cx<'_>,
        layer_index: usize,
        offset_m: f64,
    ) -> Result<NonlinearGeneratedLayerTemperatureSample, ConductionError> {
        self.temperature_profile(cx, layer_index)?
            .temperature_gradient_at(cx, offset_m)
    }

    /// Second-order observation, refusing a genuine interior curvature jump.
    pub fn temperature_jet_at(
        &self,
        cx: &Cx<'_>,
        layer_index: usize,
        offset_m: f64,
    ) -> Result<NonlinearGeneratedLayerTemperatureJet, NonlinearTemperatureProfileError> {
        self.temperature_profile(cx, layer_index)?
            .temperature_jet_at(cx, offset_m)
    }
}

impl NonlinearGeneratedLayerTemperatureProfile<'_> {
    pub const fn solution(&self) -> &NonlinearGeneratedLayeredWallSolution {
        self.solution
    }
    pub const fn layer_index(&self) -> usize {
        self.layer_index
    }
    pub const fn thickness_m(&self) -> f64 {
        self.layer.thickness_m()
    }

    /// Original solved T and T'=-q/k(T). Endpoint T and q retain the stored
    /// face values exactly. Source support and cancellation remain enforced by
    /// the original nonlinear temperature owner at every query.
    pub fn temperature_gradient_at(
        &self,
        cx: &Cx<'_>,
        offset_m: f64,
    ) -> Result<NonlinearGeneratedLayerTemperatureSample, ConductionError> {
        self.sample_with_conductivity(cx, offset_m)
            .map(|(sample, _)| sample)
    }

    // The second-order query reuses this same source evaluation instead of
    // repeating the conductivity table lookup at every optical stage.
    fn sample_with_conductivity(
        &self,
        cx: &Cx<'_>,
        offset_m: f64,
    ) -> Result<(NonlinearGeneratedLayerTemperatureSample, f64), ConductionError> {
        let temperature_k = self
            .solution
            .temperature_at(cx, self.layer_index, offset_m)?;
        positive("nonlinear local temperature", temperature_k)?;
        let row = &self.solution.layer_temperatures[self.layer_index];
        let heat_flux_w_per_m2 = if offset_m == 0.0 {
            row.face_fluxes_w_per_m2[0]
        } else if offset_m == self.layer.thickness_m() {
            row.face_fluxes_w_per_m2[1]
        } else {
            finite(
                "nonlinear local heat flux",
                self.source_w_m3
                    .mul_add(offset_m, row.face_fluxes_w_per_m2[0]),
            )?
        };
        let conductivity = self.layer.conductivity().eval(temperature_k)?;
        positive("nonlinear local conductivity", conductivity)?;
        let gradient_k_per_m = finite(
            "nonlinear temperature gradient",
            -heat_flux_w_per_m2 / conductivity,
        )?;
        if heat_flux_w_per_m2 != 0.0 && gradient_k_per_m == 0.0 {
            return Err(invalid(
                "nonzero nonlinear temperature gradient underflowed",
            ));
        }
        Ok((
            NonlinearGeneratedLayerTemperatureSample {
                temperature_k,
                heat_flux_w_per_m2,
                gradient_k_per_m,
            },
            conductivity,
        ))
    }

    /// Original T''=-(s+k'(T)*T'^2)/k(T). At an interior source knot, unequal
    /// slopes imply a curvature jump only when q is nonzero. Exact q=0 removes
    /// the k' term and admits the unique curvature -s/k. Physical layer faces
    /// use the derivative from inside this layer, not a table's arbitrary
    /// right-sided convention. No finite differences or smoothing are applied.
    pub fn temperature_jet_at(
        &self,
        cx: &Cx<'_>,
        offset_m: f64,
    ) -> Result<NonlinearGeneratedLayerTemperatureJet, NonlinearTemperatureProfileError> {
        let (sample, conductivity) = self.sample_with_conductivity(cx, offset_m)?;
        let curvature_k_per_m2 = if sample.heat_flux_w_per_m2 == 0.0 {
            -self.source_w_m3 / conductivity
        } else {
            let slope = self.curvature_slope(offset_m, sample)?;
            // Use T' rather than q²/k³ so a representable Fourier gradient
            // does not require squaring a potentially large heat flux.
            (-slope * sample.gradient_k_per_m).mul_add(sample.gradient_k_per_m, -self.source_w_m3)
                / conductivity
        };
        finite("nonlinear temperature curvature", curvature_k_per_m2)?;
        Ok(NonlinearGeneratedLayerTemperatureJet {
            sample,
            curvature_k_per_m2,
        })
    }

    fn curvature_slope(
        &self,
        offset_m: f64,
        sample: NonlinearGeneratedLayerTemperatureSample,
    ) -> Result<f64, NonlinearTemperatureProfileError> {
        let table = self.layer.conductivity();
        let knots = table.knots();
        let at = knots.partition_point(|&(temperature, _)| temperature < sample.temperature_k);
        let slope = if at > 0 && at + 1 < knots.len() && knots[at].0 == sample.temperature_k {
            let left = (knots[at].1 - knots[at - 1].1) / (knots[at].0 - knots[at - 1].0);
            let right = (knots[at + 1].1 - knots[at].1) / (knots[at + 1].0 - knots[at].0);
            if offset_m == 0.0 {
                if sample.heat_flux_w_per_m2 < 0.0 {
                    right
                } else {
                    left
                }
            } else if offset_m == self.layer.thickness_m() {
                if sample.heat_flux_w_per_m2 < 0.0 {
                    left
                } else {
                    right
                }
            } else if left == right {
                left
            } else {
                return Err(NonlinearTemperatureProfileError::UndefinedCurvature {
                    temperature_k: sample.temperature_k,
                });
            }
        } else {
            table.derivative(sample.temperature_k)?
        };
        Ok(finite("nonlinear conductivity slope", slope)?)
    }
}

/// Original source refusal, or absence of a unique second derivative.
#[derive(Clone, Debug, PartialEq)]
pub enum NonlinearTemperatureProfileError {
    Conduction(ConductionError),
    UndefinedCurvature { temperature_k: f64 },
}
impl From<ConductionError> for NonlinearTemperatureProfileError {
    fn from(error: ConductionError) -> Self {
        Self::Conduction(error)
    }
}
impl core::fmt::Display for NonlinearTemperatureProfileError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "nonlinear generated temperature profile: {self:?}")
    }
}
impl core::error::Error for NonlinearTemperatureProfileError {}
