//! Gray, opaque liquid surfaces radiating to prescribed isothermal surroundings.
//!
//! The same extensive-enthalpy owner integrates radiation, gas-film exchange,
//! contact conduction and authored heat before native momentum/VOF transport.
//! Geometry is declared at the accepted configuration, not inferred from a fill
//! fraction. Surface temperature is the resolved cell temperature: the caller
//! must resolve the liquid thermal boundary layer or justify this lumped value.
//! `Material` emissivity comes from the original liquid source bundle; the
//! separately solved interface path evaluates that same claim at interface T.

use super::*;
use crate::conjugate::STEFAN_BOLTZMANN;
use fs_material::state_point::{
    MaterialPropertySelection, ResolvedScalarProperty, resolve_material_state_point,
};
use fs_material::thermal::source::{
    EMISSIVITY_DIMS, SURFACE_EMISSIVITY_PROPERTY, ThermalMaterialPoint,
};

mod body;
pub(super) use body::BodyRadiation;
pub use body::{
    FreeSurfaceBodyRadiationExchange, FreeSurfaceBodyRadiationFace, FreeSurfaceBodyRadiationSample,
};

pub(super) fn admit_emissivity(value: f64) -> Result<(), Error> {
    if !value.is_finite() || !(0.0..=1.0).contains(&value) {
        return Err(Error::Invalid(
            "total hemispherical emissivity must lie in [0,1]",
        ));
    }
    Ok(())
}

// One fourth-power kernel for source-backed solid and liquid surfaces. The
// exact positive secant preserves small fluxes near thermal equilibrium and
// participates in the same combined thermal CFL as all other heat paths.
fn radiant_conductance(
    emissivity: f64,
    area: f64,
    temperature: f64,
    bath: f64,
) -> Result<f64, Error> {
    if area == 0.0 || emissivity == 0.0 {
        return Ok(0.0);
    }
    let value = environment::finite(
        (emissivity * STEFAN_BOLTZMANN * area)
            * (temperature + bath)
            * (temperature * temperature + bath * bath),
    )?;
    if value <= 0.0 {
        return Err(Error::Invalid("radiant conductance is unrepresentable"));
    }
    Ok(value)
}

/// The explicitly selected diffuse-gray surface law. Geometry and radiation
/// reservoirs remain boundary data; a material law follows the original liquid.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum FreeSurfaceEmissivity {
    /// A prescribed constant in [0, 1], including zero to omit radiation.
    Prescribed(f64),
    /// The canonical total hemispherical emissivity already present in the
    /// original `LiquidTemperatureLaw` bundle. Its selected claim is retained
    /// at every bulk or solved-surface temperature; there is no fresh selection.
    /// A material value must lie in [0, 1]. A known zero is retained as a
    /// reflecting surface; missing/unsupported data still refuses.
    Material,
}

impl From<f64> for FreeSurfaceEmissivity {
    fn from(value: f64) -> Self {
        Self::Prescribed(value)
    }
}

impl FreeSurfaceEmissivity {
    pub(super) fn admit(self) -> Result<(), Error> {
        match self {
            Self::Prescribed(value) if !value.is_finite() || !(0.0..=1.0).contains(&value) => Err(
                Error::Invalid("prescribed total emissivity must lie in [0,1]"),
            ),
            _ => Ok(()),
        }
    }

    fn original<'a>(
        self,
        law: &'a LiquidTemperatureLaw<'_>,
    ) -> Result<Option<&'a ResolvedScalarProperty>, Error> {
        self.admit()?;
        if self == Self::Material {
            let original = law
                .reference_state()
                .material()
                .property(SURFACE_EMISSIVITY_PROPERTY)
                .ok_or(Error::Invalid(
                    "material radiation needs the emissivity claim in the original liquid bundle",
                ))?;
            if original.requirement().dims() != EMISSIVITY_DIMS {
                return Err(Error::Invalid(
                    "total hemispherical emissivity must be dimensionless",
                ));
            }
            Ok(Some(original))
        } else {
            Ok(None)
        }
    }

    /// Source validity bounds narrow the existing interface bracket. The exact
    /// pinned query below still owns curve/hull support and may refuse a trial.
    pub(super) fn temperature_bounds(
        self,
        law: &LiquidTemperatureLaw<'_>,
    ) -> Result<Option<(f64, f64)>, Error> {
        let Some(source) = self.original(law)? else {
            return Ok(None);
        };
        let (mut lo, mut hi) = source
            .answer()
            .evidence
            .model
            .validity
            .bound(law.temperature_axis())
            .unwrap_or((0.0, f64::MAX));
        let claim = law
            .card()
            .claims()
            .claim(source.answer().receipt.selected)
            .ok_or(Error::Invalid(
                "material radiation lost its retained source claim",
            ))?;
        if let fs_matdb::PropertyValue::Curve {
            abscissa, knots, ..
        } = &claim.value
        {
            if abscissa == law.temperature_axis() {
                lo = lo.max(knots[0].0);
                hi = hi.min(knots[knots.len() - 1].0);
            }
        }
        Ok(Some((lo, hi)))
    }

    pub(super) fn resolve(
        self,
        law: &LiquidTemperatureLaw<'_>,
        liquid: &ResolvedLiquidState,
        temperature_k: f64,
    ) -> Result<(f64, Option<ResolvedScalarProperty>), Error> {
        let Some(original) = self.original(law)? else {
            let Self::Prescribed(value) = self else {
                unreachable!("material has a source")
            };
            return Ok((value, None));
        };
        let point = liquid
            .material()
            .query()
            .with(law.temperature_axis(), temperature_k)
            .map_err(|_| Error::Invalid("invalid emissivity surface-temperature coordinate"))?;
        let state = resolve_material_state_point(
            law.card(),
            &point,
            std::slice::from_ref(original.requirement()),
            MaterialPropertySelection::PinnedByProperty(vec![(
                SURFACE_EMISSIVITY_PROPERTY.into(),
                original.answer().receipt.selected,
            )]),
        )?;
        let emissivity = state
            .property(SURFACE_EMISSIVITY_PROPERTY)
            .ok_or(Error::Invalid(
                "material radiation lost its original emissivity source",
            ))?;
        ThermalMaterialPoint {
            material: law.card(),
            context: &point,
            temperature_axis: law.temperature_axis(),
            liquid_fraction_axis: None,
        }
        .admit_static(emissivity)?;
        let value = emissivity.value_si();
        admit_emissivity(value)?;
        Ok((value, Some(emissivity.clone())))
    }
}

/// One exposed patch facing a black isothermal surrounding hemisphere.
///
/// This is a diffuse-gray, opaque-surface approximation with view factor one
/// to the declared surroundings. It does not infer visibility, emissivity,
/// surface-to-surface exchange, gas absorption or evaporation from a liquid's
/// name. A gas film may exchange heat through the SAME patch in parallel.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FreeSurfaceRadiationFace {
    /// Positive-mass liquid cell on the accepted free surface.
    pub liquid: [usize; 3],
    /// Adjacent in-domain dry gas cell; wet cells and permanent walls refuse.
    pub gas: [usize; 3],
    /// Radiating area [m2], bounded by one Cartesian face. Zero disables heat.
    pub area_m2: f64,
    /// Prescribed or original material total hemispherical gray emissivity.
    pub emissivity: FreeSurfaceEmissivity,
    /// Absolute temperature of the prescribed surroundings [K], positive.
    pub surroundings_k: f64,
}

/// The physical flux and secant conductance used by one thermal substep.
#[derive(Clone, Debug, PartialEq)]
pub struct FreeSurfaceRadiationSample {
    /// Declared exposed geometry and radiative bath for this sample.
    pub face: FreeSurfaceRadiationFace,
    /// Accepted caloric temperature at the sampling epoch [K].
    pub liquid_temperature_k: f64,
    /// Actual coefficient at this temperature, not a reference-temperature cache.
    pub emissivity: f64,
    /// Original source value and receipt when material radiation was selected.
    pub emissivity_source: Option<ResolvedScalarProperty>,
    /// eps sigma A (T + Ts) (T^2 + Ts^2) [W/K]. This is the exact
    /// secant of the fourth-power law, not a linearized radiation model.
    pub conductance_w_k: f64,
    /// eps sigma A (Ts^4 - T^4) [W], positive into liquid.
    pub heat_into_liquid_w: f64,
}

/// Signed radiation over the same accepted thermal and native flow interval.
#[derive(Clone, Debug)]
pub struct FreeSurfaceRadiationExchange {
    /// Already included in thermal/flow external heat; never add it twice.
    pub heat_to_liquid_j: f64,
    /// Fluxes held during the last actually accepted thermal substep.
    pub last_applied: Vec<FreeSurfaceRadiationSample>,
    /// Supported thermal endpoint before native enthalpy advection.
    pub thermal_endpoint: Vec<FreeSurfaceRadiationSample>,
}

impl MaterialFreeSurfaceEnthalpy3<'_, '_> {
    pub(super) fn admit_radiation(
        &self,
        radiation: &[FreeSurfaceRadiationFace],
        cell_id: &impl Fn([usize; 3]) -> Result<CellId, Error>,
        dx: f64,
        gate: &CancelGate,
    ) -> Result<Vec<CellId>, Error> {
        let mut faces = BTreeSet::new();
        let mut endpoints = Vec::with_capacity(radiation.len());
        for &face in radiation {
            poll(gate)?;
            face.emissivity.admit()?;
            let id = cell_id(face.liquid)?;
            let separation = (0..3).fold(0usize, |sum, axis| {
                sum.saturating_add(face.liquid[axis].abs_diff(face.gas[axis]))
            });
            let [Ok(x), Ok(y), Ok(z)] = face.gas.map(i64::try_from) else {
                return Err(Error::Invalid(
                    "radiating gas coordinate is unrepresentable",
                ));
            };
            if separation != 1
                || self.flow.lattice.transport_wall_at([x, y, z])
                || self
                    .flow
                    .lattice
                    .sample_at(face.gas[0], face.gas[1], face.gas[2])
                    .is_some_and(|cell| cell.fill_fraction != 0.0)
                || !faces.insert((id, face.gas))
                || !face.area_m2.is_finite()
                || face.area_m2 < 0.0
                || face.area_m2 > dx * dx
                || !face.surroundings_k.is_finite()
                || face.surroundings_k <= 0.0
            {
                return Err(Error::Invalid(
                    "radiation needs a unique exposed face, supported area/emissivity and positive kelvin surroundings",
                ));
            }
            endpoints.push(id);
        }
        Ok(endpoints)
    }

    /// Advance radiant and gas-film heat, finite wall contacts, conductivity,
    /// optional molecular heat and native VOF motion on one physical interval.
    ///
    /// Radiation and gas films may occupy the same exposed patch; their watts
    /// enter separate reports and the external heat account exactly once each.
    /// Finite wall contacts and thermal interface laws retain their existing
    /// owners. If electrical heating is attached, its original owner also
    /// participates. The radiative bath is external; it is not finite gas energy.
    ///
    /// Fourth-power radiation and an explicitly selected material emissivity
    /// are reevaluated at every thermal substep. The original emissivity claim
    /// must already belong to the complete liquid law. All nonthermal source
    /// coordinates and claim identity survive; missing support refuses.
    /// Its
    /// positive secant conductance is added to gas-film and internal-contact
    /// conductances in a SINGLE monotonicity bound. Strong radiation therefore
    /// shortens the step instead of overshooting the bath or clipping energy.
    /// This remains explicit first-order thermal/flow splitting and needs time
    /// refinement for accuracy. Exhausted budgets, unsupported material states,
    /// cancellation and native flow refusals publish no partial histories.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_radiative_environment(
        &mut self,
        faces: &[FreeSurfaceConductionFace],
        heat: &[FreeSurfaceCellHeat],
        controls: FreeSurfaceConductionControls,
        bodies: &mut [FreeSurfaceThermalBody<'_>],
        contacts: &[FreeSurfaceBodyContact],
        interfaces: &[FreeSurfaceBodyInterface<'_>],
        environment: &[FreeSurfaceGasContact<'_>],
        radiation: &[FreeSurfaceRadiationFace],
        with_viscous_heating: bool,
        gate: &CancelGate,
    ) -> Result<viscous::FreeSurfaceViscousConductionReport, Error> {
        let applied = if with_viscous_heating {
            self.viscous_sources(gate)?
        } else {
            Vec::new()
        };
        let conduction = self.advance_with_radiative_environment(
            faces,
            heat,
            controls,
            &applied,
            bodies,
            contacts,
            interfaces,
            environment,
            radiation,
            None,
            gate,
        )?;
        Ok(viscous::FreeSurfaceViscousConductionReport {
            conduction,
            applied,
        })
    }
}

impl Conduction<'_, '_> {
    pub(super) fn sample_radiation(
        &self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Vec<FreeSurfaceRadiationSample>, Error> {
        let mut samples = Vec::with_capacity(self.radiation.len());
        let mut covered = BTreeMap::new();
        for (&face, &id) in self.radiation.iter().zip(&self.radiation_endpoints) {
            check()?;
            let area = covered.entry(id).or_insert(0.0);
            *area = environment::finite(*area + face.area_m2)?;
            if *area > thermal.exposed_area_m2(id)? {
                return Err(Error::Invalid("radiation exceeds exposed liquid area"));
            }
            let temperature = thermal.temperature_k(id)?;
            let (emissivity, emissivity_source) = match face.emissivity {
                FreeSurfaceEmissivity::Prescribed(value) => (value, None),
                FreeSurfaceEmissivity::Material => {
                    let liquid =
                        self.sample_liquid(id, thermal.cells()[&id].specific_h_j_kg, check)?;
                    face.emissivity.resolve(self.law, &liquid, temperature)?
                }
            };
            let conductance =
                radiant_conductance(emissivity, face.area_m2, temperature, face.surroundings_k)?;
            samples.push(FreeSurfaceRadiationSample {
                face,
                liquid_temperature_k: temperature,
                emissivity,
                emissivity_source,
                conductance_w_k: conductance,
                heat_into_liquid_w: environment::finite(
                    conductance * (face.surroundings_k - temperature),
                )?,
            });
        }
        Ok(samples)
    }
}
