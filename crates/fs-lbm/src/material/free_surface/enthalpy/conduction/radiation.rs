//! Gray, opaque liquid surfaces radiating to prescribed isothermal surroundings.
//!
//! The same extensive-enthalpy owner integrates radiation, gas-film exchange,
//! contact conduction and authored heat before native momentum/VOF transport.
//! Geometry is declared at the accepted configuration, not inferred from a fill
//! fraction. Surface temperature is the resolved cell temperature: the caller
//! must resolve the liquid thermal boundary layer or justify this lumped value.

use super::*;
use crate::conjugate::radiation::STEFAN_BOLTZMANN;

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
    /// Prescribed total hemispherical gray emissivity, in [0, 1].
    pub emissivity: f64,
    /// Absolute temperature of the prescribed surroundings [K], positive.
    pub surroundings_k: f64,
}

/// The physical flux and secant conductance used by one thermal substep.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FreeSurfaceRadiationSample {
    /// Declared exposed geometry and radiative bath for this sample.
    pub face: FreeSurfaceRadiationFace,
    /// Accepted caloric temperature at the sampling epoch [K].
    pub liquid_temperature_k: f64,
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
                || !face.emissivity.is_finite()
                || !(0.0..=1.0).contains(&face.emissivity)
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
    /// Fourth-power radiation is reevaluated at every thermal substep. Its
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
            // Factoring preserves small fluxes near equilibrium. No difference
            // of almost equal fourth powers and no empirical clipping occurs.
            let conductance = if face.area_m2 == 0.0 || face.emissivity == 0.0 {
                0.0
            } else {
                let bath = face.surroundings_k;
                let value = environment::finite(
                    (face.emissivity * STEFAN_BOLTZMANN * face.area_m2)
                        * (temperature + bath)
                        * (temperature * temperature + bath * bath),
                )?;
                if value <= 0.0 {
                    return Err(Error::Invalid("radiant conductance is unrepresentable"));
                }
                value
            };
            samples.push(FreeSurfaceRadiationSample {
                face,
                liquid_temperature_k: temperature,
                conductance_w_k: conductance,
                heat_into_liquid_w: environment::finite(
                    conductance * (face.surroundings_k - temperature),
                )?,
            });
        }
        Ok(samples)
    }
}
