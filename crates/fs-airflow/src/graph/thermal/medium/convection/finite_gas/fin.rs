//! Original material and gas films on quasi-steady constant-section fins.
//!
//! The base body retains its original native thermal storage. Fins here are
//! massless steady attachments: their axial temperature profile is eliminated,
//! not given a second heat or storage integrator. A uniform coefficient is
//! sampled from the original film at the actual base/gas trial, and constant k
//! is sampled at the base. An explicit bound admits the neglected k variation
//! over the actual tip-to-base interval, including all interior source knots.
//! This is not a distributed nonlinear k(T)/h(T) fin solve or a fin-array flow
//! correlation. Fin mass, radiation, contact resistance and mutual shielding
//! require separately resolved models.

/// Distributed source-dependent conductivity and film properties along a fin.
pub mod nonlinear;

use super::conjugate::{ConjugateFilmSelection, ConjugateGasFilm};
use super::film::{FilmGasSource, FilmGasSurface};
use super::*;
use fs_convection::fin::{FinError, FinTip, StraightFinGeometry, StraightFinResponse};
use fs_matdb::PropertyLineSlice;
use fs_material::phase::SolidLiquidPhase;
use fs_material::thermal::conduction::source::ResolvedConductionLayer;
use fs_qty::Temperature;
use fs_thermal::enthalpy::DepositionError;
use std::sync::Arc;

/// Tip convection uses the same held film h as the exposed sides and bare base.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FinTipExposure {
    Adiabatic,
    Exposed,
}

/// One repeated fin geometry, actual number of fins and the remaining bare base.
/// The covered root footprints are excluded from convection. They remain part
/// of the base body's surface used by its original native Biot/storage owner.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StraightFinArray {
    fin: StraightFinGeometry,
    count: u32,
    bare_base_area_m2: f64,
    tip: FinTipExposure,
    base_area_m2: f64,
    wetted_area_m2: f64,
}
impl StraightFinArray {
    pub fn try_new(
        fin: StraightFinGeometry,
        count: u32,
        bare_base_area_m2: f64,
        tip: FinTipExposure,
    ) -> Result<Self, SourceFinError> {
        if count == 0 || !bare_base_area_m2.is_finite() || bare_base_area_m2 < 0.0 {
            return Err(SourceFinError::Invalid(
                "fins need positive count and finite nonnegative bare base area",
            ));
        }
        let n = f64::from(count);
        let base = bare_base_area_m2 + n * fin.cross_section_m2();
        let wetted = bare_base_area_m2
            + n * (fin.lateral_area_m2()
                + if tip == FinTipExposure::Exposed {
                    fin.cross_section_m2()
                } else {
                    0.0
                });
        if !base.is_finite() || base <= 0.0 || !wetted.is_finite() || wetted <= 0.0 {
            return Err(SourceFinError::Fin(FinError::Arithmetic));
        }
        Ok(Self {
            fin,
            count,
            bare_base_area_m2,
            tip,
            base_area_m2: base,
            wetted_area_m2: wetted,
        })
    }
    pub const fn fin(self) -> StraightFinGeometry {
        self.fin
    }
    pub const fn count(self) -> u32 {
        self.count
    }
    pub const fn bare_base_area_m2(self) -> f64 {
        self.bare_base_area_m2
    }
    pub const fn tip(self) -> FinTipExposure {
        self.tip
    }
    pub const fn base_area_m2(self) -> f64 {
        self.base_area_m2
    }
    pub const fn wetted_area_m2(self) -> f64 {
        self.wetted_area_m2
    }
}

/// Explicit model-approximation limits; no default empirical allowance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SourceFinLimits {
    pub maximum_transverse_biot: f64,
    pub maximum_relative_conductivity_variation: f64,
}

#[derive(Debug)]
pub enum SourceFinError {
    Invalid(&'static str),
    Film(MaterialConvectionError),
    Fin(FinError),
    Thermal(DepositionError),
    TransverseBiot { actual: f64, maximum: f64 },
    ConductivityVariation { actual: f64, maximum: f64 },
}
impl core::fmt::Display for SourceFinError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "source fin convection: {self:?}")
    }
}
impl std::error::Error for SourceFinError {}
impl From<MaterialConvectionError> for SourceFinError {
    fn from(e: MaterialConvectionError) -> Self {
        Self::Film(e)
    }
}
impl From<FinError> for SourceFinError {
    fn from(e: FinError) -> Self {
        Self::Fin(e)
    }
}
impl From<DepositionError> for SourceFinError {
    fn from(e: DepositionError) -> Self {
        Self::Thermal(e)
    }
}

/// Source-bearing fin closure. The film's raw h*A is retained as evidence only;
/// total_conductance_w_k is the SINGLE heat connection applied by the owner.
#[derive(Debug)]
pub struct SourceFinSample<S> {
    pub film: S,
    pub conductivity_source: Arc<PropertyLineSlice>,
    pub geometry: StraightFinArray,
    pub base_temperature_k: f64,
    pub bulk_temperature_k: f64,
    pub tip_temperature_k: f64,
    pub fin: StraightFinResponse,
    pub relative_conductivity_variation: f64,
    /// Largest h A/(k(T) P) on the admitted original source interval.
    /// The fin response separately reports its held-base-k value.
    pub maximum_source_transverse_biot: f64,
    pub total_conductance_w_k: f64,
    /// Signed heat into the gas, including the bare base and all fins once.
    pub applied_heat_into_gas_w: f64,
}

/// Original immutable film and solid conductivity source with explicit geometry.
/// Generic over the SAME film policy as the existing conducting-surface owner,
/// including original specialized tube/sphere and prescribed-flow selections.
pub struct SourceFinGasSurface<'a, 'm, F: ConjugateGasFilm> {
    film: F,
    material: &'a ResolvedConductionLayer<'m>,
    geometry: StraightFinArray,
    limits: SourceFinLimits,
}
impl<'a, 'm, F: ConjugateGasFilm> SourceFinGasSurface<'a, 'm, F> {
    pub fn new(
        film: F,
        material: &'a ResolvedConductionLayer<'m>,
        geometry: StraightFinArray,
        limits: SourceFinLimits,
    ) -> Result<Self, SourceFinError> {
        if material.phase() != SolidLiquidPhase::Solid {
            return Err(SourceFinError::Invalid(
                "a fixed fin requires an original all-solid conductivity branch",
            ));
        }
        if [
            limits.maximum_transverse_biot,
            limits.maximum_relative_conductivity_variation,
        ]
        .into_iter()
        .any(|v| !v.is_finite() || v < 0.0)
        {
            return Err(SourceFinError::Invalid(
                "finite nonnegative fin approximation limits required",
            ));
        }
        Ok(Self {
            film,
            material,
            geometry,
            limits,
        })
    }
    pub const fn film(&self) -> &F {
        &self.film
    }
    pub const fn material(&self) -> &ResolvedConductionLayer<'m> {
        self.material
    }
    pub const fn geometry(&self) -> StraightFinArray {
        self.geometry
    }
    pub const fn limits(&self) -> SourceFinLimits {
        self.limits
    }

    /// Resolve one original bulk/base state, retaining all original gas and solid
    /// evidence. The exposed-film area must match actual wetted fins plus bare
    /// base; the covered root footprints are never added to that area.
    ///
    /// k is held at its base value for this axial BVP. The entire resulting
    /// monotone temperature interval is admitted against the original source.
    /// Work is bounded by its original 4096-knot cap, with checkpoints every 64
    /// visited interior knots. No source interpolation or material name is added.
    pub fn resolve(
        &self,
        gas: &F::State,
        base: Temperature,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<ConductanceSample<SourceFinSample<F::Sample>>, SourceFinError> {
        checkpoint()?;
        let g = self.geometry;
        let area = self.film.admitted_area_m2(gas)?;
        if area != g.wetted_area_m2 {
            return Err(SourceFinError::Invalid(
                "fin film must cover the actual bare-base plus wetted-fin area exactly",
            ));
        }
        let curve = self.material.conductivity_curve();
        let k = curve.conductivity_w_m_k(base.value())?;
        let film = self.film.resolve_film_at_surface(gas, base)?;
        let h = film.conductance_w_k / area;
        let response = StraightFinResponse::try_new(
            g.fin,
            k,
            h,
            match g.tip {
                FinTipExposure::Adiabatic => FinTip::Adiabatic,
                FinTipExposure::Exposed => FinTip::Convection {
                    coefficient_w_m2_k: h,
                },
            },
        )?;
        if response.transverse_biot() > self.limits.maximum_transverse_biot {
            return Err(SourceFinError::TransverseBiot {
                actual: response.transverse_biot(),
                maximum: self.limits.maximum_transverse_biot,
            });
        }
        let bulk = gas.temperature_k();
        let tip = response.temperature_k_at(g.fin.length_m(), base.value(), bulk)?;
        let tip_k = curve.conductivity_w_m_k(tip)?;
        let mut variation = (tip_k - k).abs() / k;
        let mut minimum_conductivity = k.min(tip_k);
        let lo = tip.min(base.value());
        let hi = tip.max(base.value());
        let knots = curve.knots();
        if hi > lo {
            let first = knots.partition_point(|v| v[0] <= lo);
            let last = knots.partition_point(|v| v[0] < hi);
            for (index, knot) in knots[first..last].iter().enumerate() {
                if index % 64 == 0 {
                    checkpoint()?;
                }
                variation = variation.max((knot[1] - k).abs() / k);
                minimum_conductivity = minimum_conductivity.min(knot[1]);
            }
        }
        if !variation.is_finite() || variation > self.limits.maximum_relative_conductivity_variation
        {
            return Err(SourceFinError::ConductivityVariation {
                actual: variation,
                maximum: self.limits.maximum_relative_conductivity_variation,
            });
        }
        let maximum_source_biot =
            (h / minimum_conductivity) * (g.fin.cross_section_m2() / g.fin.exposed_perimeter_m());
        if !maximum_source_biot.is_finite()
            || maximum_source_biot > self.limits.maximum_transverse_biot
        {
            return Err(SourceFinError::TransverseBiot {
                actual: maximum_source_biot,
                maximum: self.limits.maximum_transverse_biot,
            });
        }
        let conductance = h * g.bare_base_area_m2 + f64::from(g.count) * response.conductance_w_k();
        let heat = conductance * (base.value() - bulk);
        if !conductance.is_finite()
            || conductance < 0.0
            || !heat.is_finite()
            || (h > 0.0 && conductance == 0.0)
            || (conductance > 0.0 && base.value() != bulk && heat == 0.0)
        {
            return Err(SourceFinError::Fin(FinError::Arithmetic));
        }
        checkpoint()?;
        Ok(ConductanceSample {
            conductance_w_k: conductance,
            source: SourceFinSample {
                film: film.source,
                conductivity_source: Arc::clone(self.material.source_line()),
                geometry: g,
                base_temperature_k: base.value(),
                bulk_temperature_k: bulk,
                tip_temperature_k: tip,
                fin: response,
                relative_conductivity_variation: variation,
                maximum_source_transverse_biot: maximum_source_biot,
                total_conductance_w_k: conductance,
                applied_heat_into_gas_w: heat,
            },
        })
    }

    /// Reuse this source closure on the original single-, multiple-body or
    /// adaptive wall-feedback residual. Native area is the BASE surface,
    /// including covered root footprints, not the larger fin wetted area.
    pub fn resolve_isochoric_coefficient(
        &self,
        cx: &Cx<'_>,
        gas: &F::State,
        body_surface_area_m2: f64,
        base: Temperature,
    ) -> Result<ConvectionCoefficient<SourceFinSample<F::Sample>>, SourceFinError> {
        if body_surface_area_m2 != self.geometry.base_area_m2 {
            return Err(SourceFinError::Invalid(
                "native fin body must cover the original base surface exactly",
            ));
        }
        let sample = self.resolve(gas, base, || {
            cx.checkpoint().map_err(|_| DepositionError::Cancelled)
        })?;
        let coefficient = sample.conductance_w_k / body_surface_area_m2;
        if !coefficient.is_finite() || coefficient < 0.0 {
            return Err(SourceFinError::Fin(FinError::Arithmetic));
        }
        Ok(ConvectionCoefficient {
            coefficient_w_per_m2_k: coefficient,
            source: sample.source,
        })
    }

    /// Original finite gas/body transaction with nonlinear base/film feedback.
    /// The body owns only its declared inventory; no fin mass or enthalpy is
    /// inferred. Native body radiation must be disabled because its uniform
    /// base surface does not represent the eliminated fin temperature profile.
    /// Initial and every trial/actual endpoint keep their original source gates.
    pub fn isochoric_body_candidate<G: IsochoricGasStorage<State = F::State>>(
        &self,
        cx: &Cx<'_>,
        environment: &FiniteGasThermalEnvironment<'_, G>,
        body: UniformEnthalpyStepInput<'_>,
        wall: WallConvectionSolve,
    ) -> Result<
        UniformEnthalpyStep<
            FiniteGasThermalReport<WallConvectionSample<SourceFinSample<F::Sample>>, G>,
        >,
        ConvectionError<WallConvectionError<SourceFinError>>,
    > {
        let source_error = |e| match e {
            SourceFinError::Thermal(DepositionError::Cancelled) => {
                ConvectionError::Thermal(ThermalPistonError::Cancelled)
            }
            other => ConvectionError::Convection(WallConvectionError::Source(other)),
        };
        cx.checkpoint()
            .map_err(|_| ConvectionError::Thermal(ThermalPistonError::Cancelled))?;
        if !environment.body.transport.radiation_is_disabled() {
            return Err(source_error(SourceFinError::Invalid(
                "quasi-steady fins require disabled body-node radiation",
            )));
        }
        self.resolve_isochoric_coefficient(
            cx,
            environment.gas.state(),
            body.surface_area_m2,
            Temperature::new(body.initial.temperature_k()),
        )
        .map_err(source_error)?;
        environment
            .advance_with_wall_convection(cx, body, wall, |gas, t| {
                self.resolve_isochoric_coefficient(
                    cx,
                    gas,
                    body.surface_area_m2,
                    Temperature::new(t),
                )
            })
            .map_err(|e| match e {
                ConvectionError::Convection(WallConvectionError::Source(
                    SourceFinError::Thermal(DepositionError::Cancelled),
                )) => ConvectionError::Thermal(ThermalPistonError::Cancelled),
                other => other,
            })
    }
}

impl<'g, M: FilmGasSource> FilmGasSurface<'g, M> {
    /// Attach source-owned straight fins to this original film selection.
    pub fn with_straight_fins<'a, 'm>(
        self,
        material: &'a ResolvedConductionLayer<'m>,
        geometry: StraightFinArray,
        limits: SourceFinLimits,
    ) -> Result<SourceFinGasSurface<'a, 'm, ConjugateFilmSelection<'g, M>>, SourceFinError> {
        SourceFinGasSurface::new(
            ConjugateFilmSelection::new(self, None),
            material,
            geometry,
            limits,
        )
    }
}
