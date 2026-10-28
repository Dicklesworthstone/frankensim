//! Heat from a declared environmental gas film reaches native liquid motion.
//!
//! This is a continuum, quasi-steady two-layer Fourier closure with prescribed
//! gas temperature, pressure and composition. It does not infer a film thickness
//! or a convective coefficient from a gas name. The original material owners
//! evaluate liquid conductivity at the cell temperature and gas conductivity
//! at the mean cell/bath temperature. The original thermal timeline integrates
//! that heat before native VOF/enthalpy transport and viscosity refresh.
//! The optional finite atmosphere supplies evolving gas T and EOS pressure
//! from its original internal-energy owner while retaining the same film law.

use super::*;
use fs_matdb::QueryPoint;
use fs_material::thermochemical_gas::{ResolvedThermochemicalGasState, ThermochemicalGasModel};

/// One explicitly exposed liquid/gas patch at the accepted flow configuration.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FreeSurfaceGasFace {
    pub liquid: [usize; 3],
    /// Face-adjacent in-domain gas cell for Cartesian exchange. Resolved Local
    /// surface exchange instead names this partial cell's gas subvolume using
    /// the same coordinate as `liquid`; other wet-cell bindings refuse.
    pub gas: [usize; 3],
    /// Physical patch area [m2]. Cartesian exchange admits at most one cell
    /// face and zero insulates. Resolved graph/Local surface exchange requires
    /// its complete positive reconstructed patch area.
    pub area_m2: f64,
    /// Liquid-centroid-to-interface and gas-film thicknesses [m]. Both are
    /// positive. Cartesian exchange bounds the liquid distance by cell width;
    /// resolved surface exchange retains an explicit positive normal distance.
    pub distances_m: [f64; 2],
}

/// Source-owned gas bath on an actual exposed patch. All context axes except
/// the film evaluation temperature stay fixed during this flow interval.
#[derive(Clone, Copy, Debug)]
pub struct FreeSurfaceGasContact<'a> {
    pub face: FreeSurfaceGasFace,
    pub model: &'a ThermochemicalGasModel,
    /// Typed bath temperature and pressure plus any original source conditions.
    /// The bath is an external reservoir; its mass/energy are not simulated.
    pub context: &'a QueryPoint,
}

/// Physical coefficients and sources actually used for one exposed patch.
#[derive(Clone, Debug, PartialEq)]
pub struct FreeSurfaceGasSample {
    pub face: FreeSurfaceGasFace,
    pub bath_temperature_k: f64,
    pub liquid_temperature_k: f64,
    pub liquid: ResolvedLiquidState,
    /// Full species, EOS and transport snapshot at the declared film temperature.
    pub gas: ResolvedThermochemicalGasState,
    pub conductance_w_k: f64,
    /// Positive into the liquid, negative into the gas bath [W]. A finite
    /// atmosphere stores the exact opposite heat on the same interval.
    pub heat_into_liquid_w: f64,
}

#[derive(Clone, Debug)]
pub struct FreeSurfaceGasExchange {
    /// Signed bath supply [J], already included in the thermal and flow reports'
    /// external heat. With `finite_gas`, its opposite is stored in the gas and
    /// this transfer is internal to their union. Never add it twice.
    pub heat_to_liquid_j: f64,
    pub last_applied: Vec<FreeSurfaceGasSample>,
    /// Supported endpoint before the subsequent native advection/momentum step.
    pub thermal_endpoint: Vec<FreeSurfaceGasSample>,
}

pub(super) fn finite(value: f64) -> Result<f64, Error> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(Error::Invalid(
            "unrepresentable environmental heat exchange",
        ))
    }
}

pub(super) fn same_sources(
    a: &ResolvedThermochemicalGasState,
    b: &ResolvedThermochemicalGasState,
) -> bool {
    a.components().len() == b.components().len()
        && a.components().iter().zip(b.components()).all(|(a, b)| {
            a.card_identity() == b.card_identity()
                && a.properties().len() == b.properties().len()
                && a.properties().iter().zip(b.properties()).all(|(a, b)| {
                    a.requirement() == b.requirement()
                        && a.answer().receipt.selected == b.answer().receipt.selected
                })
        })
}

impl MaterialFreeSurfaceEnthalpy3<'_, '_> {
    pub(super) fn admit_environment(
        &self,
        environment: &[FreeSurfaceGasContact<'_>],
        cell_id: &impl Fn([usize; 3]) -> Result<CellId, Error>,
        dx: f64,
        gate: &CancelGate,
    ) -> Result<Vec<CellId>, Error> {
        let mut faces = BTreeSet::new();
        let mut endpoints = Vec::with_capacity(environment.len());
        for contact in environment {
            poll(gate)?;
            let face = contact.face;
            let id = cell_id(face.liquid)?;
            let separation = (0..3).fold(0usize, |sum, axis| {
                sum.saturating_add(face.liquid[axis].abs_diff(face.gas[axis]))
            });
            let gas = face.gas.map(i64::try_from);
            let [Ok(x), Ok(y), Ok(z)] = gas else {
                return Err(Error::Invalid(
                    "environment gas coordinate is unrepresentable",
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
                || face.distances_m.iter().any(|d| !d.is_finite() || *d <= 0.0)
                || face.distances_m[0] > dx
            {
                return Err(Error::Invalid(
                    "environment needs one unique exposed face and an adjacent gas cell",
                ));
            }
            // Bath support is required even if the interpolated film happens
            // to lie inside the source domain or the declared area is zero.
            contact
                .model
                .evaluate(contact.context)
                .map_err(MaterialLbmError::from)?;
            endpoints.push(id);
        }
        Ok(endpoints)
    }

    /// Advance liquid conduction, finite solid contacts and prescribed gas-film
    /// exchange on one thermal timeline, then the SAME native flow and enthalpy.
    /// Supply current exposed geometry each call; a liquid-covered gas neighbor
    /// refuses instead of exchanging heat through a stale free surface.
    ///
    /// `with_viscous_heating` selects the original molecular stress-power owner.
    /// The environment is a reservoir: bath heat is external, wall contact heat
    /// is internal, and molecular heating has its separate existing account.
    /// The gas and liquid source claims used at the initial film state must
    /// remain supported through every substep and endpoint. No source switches,
    /// evaporation, gas mass transport, gas pressure loading, surface tension,
    /// radiation, rarefied-gas transfer or resolved gas convection are supplied.
    ///
    /// This is the existing first-order thermal/flow split. Conductances are
    /// frozen over each explicit substep, with a combined monotonicity bound;
    /// nonlinear source accuracy still requires timestep refinement. Any error
    /// leaves the complete liquid and finite-wall histories unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_gas_environment(
        &mut self,
        faces: &[FreeSurfaceConductionFace],
        heat: &[FreeSurfaceCellHeat],
        controls: FreeSurfaceConductionControls,
        bodies: &mut [FreeSurfaceThermalBody<'_>],
        contacts: &[FreeSurfaceBodyContact],
        environment: &[FreeSurfaceGasContact<'_>],
        with_viscous_heating: bool,
        gate: &CancelGate,
    ) -> Result<viscous::FreeSurfaceViscousConductionReport, Error> {
        let applied = if with_viscous_heating {
            self.viscous_sources(gate)?
        } else {
            Vec::new()
        };
        let conduction = self.advance_with_environment(
            faces,
            heat,
            controls,
            &applied,
            bodies,
            contacts,
            environment,
            gate,
        )?;
        Ok(viscous::FreeSurfaceViscousConductionReport {
            conduction,
            applied,
        })
    }
}

impl Conduction<'_, '_> {
    pub(super) fn sample_environment(
        &self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Vec<FreeSurfaceGasSample>, Error> {
        let mut samples = Vec::with_capacity(self.environment.len());
        let mut covered = BTreeMap::new();
        for (index, (contact, &id)) in self
            .environment
            .iter()
            .zip(&self.environment_endpoints)
            .enumerate()
        {
            check()?;
            let area = covered.entry(id).or_insert(0.0);
            *area = finite(*area + contact.face.area_m2)?;
            if *area > thermal.exposed_area_m2(id)? {
                return Err(Error::Invalid("environment exceeds exposed liquid area"));
            }
            let h = thermal.cells()[&id].specific_h_j_kg;
            let liquid = self.sample_liquid(id, h, check)?;
            let temperature = thermal.temperature_k(id)?;
            let (model, context) = self.finite_atmosphere.as_ref().map_or(
                (contact.model, contact.context),
                |gas| (gas.parcel.model(), gas.parcel.point()),
            );
            let bath = context.axes()["temperature"];
            let film = 0.5 * temperature + 0.5 * bath;
            let point = context
                .clone()
                .with("temperature", film)
                .map_err(|_| Error::Invalid("invalid environmental film temperature"))?;
            let gas = model.evaluate(&point)
                .map_err(MaterialLbmError::from)?;
            if self
                .environment_reference
                .get(index)
                .is_some_and(|start| !same_sources(&start.gas, &gas))
            {
                return Err(Error::Invalid(
                    "environmental film left its original transport source claims",
                ));
            }
            let conductance = series_conductance_w_k(
                contact.face.area_m2,
                contact.face.distances_m,
                [
                    liquid.thermal_conductivity_w_m_k(),
                    gas.state().thermal_conductivity,
                ],
            )?;
            samples.push(FreeSurfaceGasSample {
                face: contact.face,
                bath_temperature_k: bath,
                liquid_temperature_k: temperature,
                liquid,
                gas,
                conductance_w_k: conductance,
                heat_into_liquid_w: finite(conductance * (bath - temperature))?,
            });
        }
        Ok(samples)
    }

    pub(super) fn environment_timestep(
        &self,
        thermal: &DepositionSimulation<'_>,
        proposed_s: f64,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<f64, Error> {
        let mut outgoing = BTreeMap::new();
        for (&id, sample) in self
            .environment_endpoints
            .iter()
            .zip(&self.environment_last)
        {
            check()?;
            if sample.conductance_w_k == 0.0 {
                continue;
            }
            let g = outgoing.entry(id).or_insert(0.0);
            *g = finite(*g + sample.conductance_w_k)?;
        }
        for (&id, sample) in self
            .radiation_endpoints
            .iter()
            .zip(&self.radiation_last)
        {
            check()?;
            if sample.conductance_w_k != 0.0 {
                let g = outgoing.entry(id).or_insert(0.0);
                *g = finite(*g + sample.conductance_w_k)?;
            }
        }
        for (&id, sample) in self.surface_endpoints.iter().zip(&self.surface_last) {
            check()?;
            let g = outgoing.entry(id).or_insert(0.0);
            *g = finite(*g + sample.liquid_decay_conductance_w_k)?;
        }
        if outgoing.is_empty() {
            return Ok(proposed_s);
        }
        // The thermal owner has already installed this substep's material
        // conductances. Include them WITH gas-film and radiative decay, not as
        // separate CFL minima, so their sum cannot overshoot a passive bath.
        for link in thermal.links().values() {
            check()?;
            for id in [link.a, link.b] {
                if let Some(g) = outgoing.get_mut(&id) {
                    *g = finite(*g + link.conductance_w_k)?;
                }
            }
        }
        let mut dt = proposed_s;
        for (id, g) in outgoing {
            check()?;
            let decay = finite(
                g * self.caloric.maximum_temperature_slope() / thermal.cells()[&id].mass_kg,
            )?;
            if decay > 0.0 {
                // A half-capacity numerical margin preserves a convex thermal
                // update for frozen coefficients. It is not a material law.
                dt = dt.min(0.5 / decay);
            }
        }
        if let Some(gas) = &self.finite_atmosphere {
            let mut conductance = 0.0;
            for sample in &self.environment_last {
                check()?;
                conductance = finite(conductance + sample.conductance_w_k)?;
            }
            // Eliminating the surface also couples gas to the radiative bath.
            // Its finite capacity must include BOTH gas-side decay channels.
            for sample in &self.surface_last {
                check()?;
                conductance = finite(conductance + sample.gas_decay_conductance_w_k)?;
            }
            let decay = finite(gas.temperature_slope * conductance / gas.parcel.mass_kg())?;
            if decay > 0.0 {
                // This is an additional finite-capacity row of the SAME heat
                // operator. Liquid rows already include film + contact loss.
                dt = dt.min(0.5 / decay);
            }
        }
        Ok(dt)
    }
}
