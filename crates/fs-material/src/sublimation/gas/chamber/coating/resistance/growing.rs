//! Mass-derived area and thermal resistance in the original finite-gas history.
//!
//! The original phase geometry owns the coating mass, thickness and exposed
//! area. The original frozen path queries conductivity on the deposit's
//! specific-enthalpy chart at each interval start. The integrated path instead
//! evaluates the original conductivity integral across the actual temperature
//! drop at every thermal trial, including its finite-contact face solve.
//! Geometry is refreshed between physical half steps. The same two-temperature
//! chamber owns molecular streams, heat storage and one gas reconstruction.
//!
//! This is a shape-preserving, constant-density two-node approximation, not a
//! resolved temperature profile, pore model or moving gas-volume solve. The
//! gas's supplied volume is the accessible volume after fixed supports; the
//! total deposit volume must stay below the explicit dilute-volume ceiling.
//! That ceiling bounds omitted displacement, not pressure-work or total error.
//! Both initial and final films remain positive; birth and extinction need a
//! separate interface transition.

use super::super::super::Total;
use super::advance::{
    GasResistiveAdvanceCandidate, GasResistiveAdvanceError, GasResistiveAdvanceSolve,
};
use super::conductive::ConductiveCoatingInterface;
use super::driven::GasResistivePower;
use super::{
    ClosedGasParcel, GasResistiveDrive, GasResistiveSurface, GasSurfaceError,
    ResistiveSurfaceState, ThermalAccommodation, ThermochemicalGasModel, poll,
};
use crate::thermal::conduction::{ConductivityCurve, LayerConductionConfig};
use fs_thermal::enthalpy::DepositionError;
use fs_thermal::enthalpy::conduction::ConductivityLaw;
use fs_thermal::sublimation::coating::resistance::ResistiveInterfaceLaw;
use fs_thermal::sublimation::coating::resistance::geometry::{
    ConductiveCoatingEndpoint, ConductiveCoatingGeometry,
};
use std::cell::RefCell;
use std::collections::BTreeMap;

/// The original material transport and explicit supported-layer geometry.
#[derive(Clone, Copy)]
pub struct GasGrowingCoating<'a> {
    pub geometry: ConductiveCoatingGeometry,
    /// Must use the SAME deposit-specific enthalpy chart as the surface.
    /// Source-backed conduction laws retain their own selected material data.
    pub conductivity: &'a dyn ConductivityLaw,
    /// Per contact area, not per exposed outer-shell area [m2 K/W].
    pub contact_resistance_m2_k_w: f64,
}

/// A growing layer with its original conductivity integral across the actual
/// temperature drop. No endpoint conductivity is frozen in its thermal solve.
#[derive(Clone, Copy)]
pub struct GasIntegratedGrowingCoating<'a> {
    pub geometry: ConductiveCoatingGeometry,
    pub conductivity: &'a ConductivityCurve,
    pub contact_resistance_m2_k_w: f64,
    /// Bounds each original material-face/contact solve inside a thermal trial.
    pub interface_config: LayerConductionConfig,
}

trait GrowingModel {
    const NONLINEAR_INTERFACE: bool = false;
    fn geometry(&self) -> ConductiveCoatingGeometry;
    fn sample(
        &self,
        surface: &GasResistiveSurface<'_>,
        state: ResistiveSurfaceState,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ConductiveCoatingEndpoint, GasResistiveAdvanceError>;
    fn interface(
        &self,
        _mass_kg: f64,
    ) -> Result<Option<ConductiveCoatingInterface<'_>>, GasResistiveAdvanceError> {
        Ok(None)
    }
}

impl GrowingModel for GasGrowingCoating<'_> {
    fn geometry(&self) -> ConductiveCoatingGeometry {
        self.geometry
    }
    fn sample(
        &self,
        surface: &GasResistiveSurface<'_>,
        state: ResistiveSurfaceState,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ConductiveCoatingEndpoint, GasResistiveAdvanceError> {
        poll(cancelled)?;
        sample_frozen(self, surface, state)
    }
}

impl GrowingModel for GasIntegratedGrowingCoating<'_> {
    const NONLINEAR_INTERFACE: bool = true;
    fn geometry(&self) -> ConductiveCoatingGeometry {
        self.geometry
    }
    fn interface(
        &self,
        mass_kg: f64,
    ) -> Result<Option<ConductiveCoatingInterface<'_>>, GasResistiveAdvanceError> {
        Ok(Some(ConductiveCoatingInterface::new(
            self.geometry,
            mass_kg,
            self.conductivity,
            self.contact_resistance_m2_k_w,
            self.interface_config,
        )?))
    }
    fn sample(
        &self,
        surface: &GasResistiveSurface<'_>,
        state: ResistiveSurfaceState,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ConductiveCoatingEndpoint, GasResistiveAdvanceError> {
        poll(cancelled)?;
        let deposit_t = (surface.deposit_temperature)(state.deposit.specific_enthalpy_j_kg)
            .map_err(GasSurfaceError::from)?;
        let support_t = (surface.substrate_temperature)(state.substrate_specific_enthalpy_j_kg)
            .map_err(GasSurfaceError::from)?;
        self.geometry
            .admit_density(
                surface
                    .law
                    .properties(deposit_t)
                    .map_err(GasSurfaceError::from)?
                    .solid_density_kg_m3,
            )
            .map_err(GasSurfaceError::from)?;
        let interface = ConductiveCoatingInterface::new(
            self.geometry,
            state.deposit.mass_kg,
            self.conductivity,
            self.contact_resistance_m2_k_w,
            self.interface_config,
        )?;
        let response = interface.evaluate([support_t, deposit_t], || {
            if cancelled() {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        })?;
        let face_t = response.material_face_temperature_k;
        // Average the ORIGINAL conductivity across only the material segment.
        // Avoid subtracting near-equal total/contact resistances when the
        // contact dominates. The equal-temperature limit uses its exact knot law.
        let average_k = if face_t == deposit_t {
            self.conductivity.conductivity_w_m_k(deposit_t)?
        } else {
            self.conductivity.integral_w_m(deposit_t, face_t)? / (face_t - deposit_t)
        };
        self.geometry
            .at_mass(
                state.deposit.mass_kg,
                average_k,
                self.contact_resistance_m2_k_w,
            )
            .map_err(GasSurfaceError::from)
            .map_err(Into::into)
    }
}

/// Temporal controls and the explicit fixed-gas-volume approximation domain.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasGrowingCoatingSolve {
    pub advance: GasResistiveAdvanceSolve,
    /// Sum of deposit volumes / fixed accessible gas volume, in (0, 0.01].
    /// Fixed substrates are already excluded from the supplied gas volume.
    pub maximum_deposit_volume_fraction: f64,
    pub absolute_extent_tolerance_m: f64,
    pub absolute_conductance_tolerance_w_k: f64,
}

/// Actual geometry and transport at the ends of the accepted history.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasGrowingCoatingHistory {
    pub id: u64,
    pub initial: ConductiveCoatingEndpoint,
    pub final_geometry: ConductiveCoatingEndpoint,
    /// Exact mass-derived shell displacement; positive retreat, negative growth.
    pub normal_recession_m: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct GasGrowingCoatingCandidate<'catalog> {
    history: GasResistiveAdvanceCandidate<'catalog>,
    geometry: Vec<GasGrowingCoatingHistory>,
}
impl<'catalog> GasGrowingCoatingCandidate<'catalog> {
    #[must_use]
    pub const fn history(&self) -> &GasResistiveAdvanceCandidate<'catalog> {
        &self.history
    }
    #[must_use]
    pub fn geometry(&self) -> &[GasGrowingCoatingHistory] {
        &self.geometry
    }
    #[must_use]
    pub const fn gas(&self) -> &ClosedGasParcel {
        self.history.gas()
    }
    #[must_use]
    pub fn into_parts(
        self,
    ) -> (
        GasResistiveAdvanceCandidate<'catalog>,
        Vec<GasGrowingCoatingHistory>,
    ) {
        (self.history, self.geometry)
    }
}

fn sample_frozen(
    model: &GasGrowingCoating<'_>,
    surface: &GasResistiveSurface<'_>,
    state: ResistiveSurfaceState,
) -> Result<ConductiveCoatingEndpoint, GasResistiveAdvanceError> {
    let temperature = (surface.deposit_temperature)(state.deposit.specific_enthalpy_j_kg)
        .map_err(GasSurfaceError::from)?;
    let phase = surface
        .law
        .properties(temperature)
        .map_err(GasSurfaceError::from)?;
    model
        .geometry
        .admit_density(phase.solid_density_kg_m3)
        .map_err(GasSurfaceError::from)?;
    let conductivity = model
        .conductivity
        .conductivity_w_m_k(state.deposit.specific_enthalpy_j_kg)
        .map_err(GasResistiveAdvanceError::Thermal)?;
    model
        .geometry
        .at_mass(
            state.deposit.mass_kg,
            conductivity,
            model.contact_resistance_m2_k_w,
        )
        .map_err(GasSurfaceError::from)
        .map_err(Into::into)
}

fn admit_volume(
    volume_m3: f64,
    gas: &ClosedGasParcel,
    maximum: f64,
) -> Result<(), GasResistiveAdvanceError> {
    let fraction = volume_m3 / gas.volume_m3();
    if !fraction.is_finite() || fraction < 0.0 || fraction > maximum {
        return Err(GasResistiveAdvanceError::Invalid(
            "deposit volume exceeds the declared fixed-gas-volume approximation",
        ));
    }
    Ok(())
}

fn ratio(
    coarse: f64,
    fine: f64,
    initial: f64,
    absolute: f64,
    relative: f64,
) -> Result<f64, GasResistiveAdvanceError> {
    let scale = absolute + relative * initial.abs().max(fine.abs());
    let result = (coarse - fine).abs() / scale;
    if !scale.is_finite() || scale <= 0.0 || !result.is_finite() {
        return Err(GasResistiveAdvanceError::Invalid(
            "growing coating error scale is unrepresentable",
        ));
    }
    Ok(result)
}

/// Advance one growing, thermally resistive coating system and its finite gas.
///
/// Layer models must cover the exact surface roster. Area and conductance in
/// the input surface records are replaced by mass-derived values before
/// admission. The masses and both original caloric coordinates are retained.
/// The selected models are sorted once; each trial uses a linear walk with
/// borrowed material laws, without repeated material-source resolution.
///
/// Only the two-half-step conservative history is accepted. Besides the
/// controller's two temperatures, independent storage increments and species
/// inventories, both coating thickness and thermal conductance participate in
/// temporal error control. Invalid source states and cancellation publish no
/// partial geometry, gas or heat history.
#[allow(clippy::too_many_arguments)]
pub fn advance_growing_resistive_chamber<'catalog>(
    surfaces: &[GasResistiveSurface<'_>],
    layers: &BTreeMap<u64, GasGrowingCoating<'_>>,
    accommodation: Option<&BTreeMap<u64, ThermalAccommodation>>,
    films: Option<&BTreeMap<u64, f64>>,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    duration_s: f64,
    controls: GasGrowingCoatingSolve,
    cancelled: impl FnMut() -> bool,
) -> Result<GasGrowingCoatingCandidate<'catalog>, GasResistiveAdvanceError> {
    let mut films = films;
    advance_growing(
        surfaces,
        layers,
        accommodation,
        &mut films,
        None,
        gas,
        catalog,
        duration_s,
        controls,
        cancelled,
    )
}

/// Evolve actual coating mass/area with the original nonlinear conductivity
/// integral and an optional external support/deposit power supply.
///
/// Each interval freezes only geometry. Every thermal endpoint evaluates the
/// original full conductivity curve and finite-contact face solve at BOTH
/// current body temperatures. Geometry and secant-circuit metadata are refreshed
/// after every physical half step and join the shared temporal error estimate.
///
/// Source powers are constant over this requested duration; only accepted
/// physical half steps count toward supplied heat. An authored power interface
/// is refused because this layer already owns the support/deposit heat law.
/// Supply piecewise power schedules as successive explicit time windows.
#[allow(clippy::too_many_arguments)]
pub fn advance_growing_integrated_chamber<'catalog>(
    surfaces: &[GasResistiveSurface<'_>],
    layers: &BTreeMap<u64, GasIntegratedGrowingCoating<'_>>,
    accommodation: Option<&BTreeMap<u64, ThermalAccommodation>>,
    films: Option<&BTreeMap<u64, f64>>,
    powers: Option<&BTreeMap<u64, GasResistivePower<'_>>>,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    duration_s: f64,
    controls: GasGrowingCoatingSolve,
    cancelled: impl FnMut() -> bool,
) -> Result<GasGrowingCoatingCandidate<'catalog>, GasResistiveAdvanceError> {
    let mut films = films;
    advance_growing(
        surfaces,
        layers,
        accommodation,
        &mut films,
        powers,
        gas,
        catalog,
        duration_s,
        controls,
        cancelled,
    )
}

/// Evolve an original enthalpy-chart conductivity law with live gas films.
///
/// This retains the frozen-within-interval conductivity model while refreshing
/// gas-film transport after each mass-owned area rebound. Optional prescribed
/// powers use the same accepted-history accounting as the integrated path.
#[allow(clippy::too_many_arguments)]
pub fn advance_growing_resistive_chamber_resolved<'catalog, F, E>(
    surfaces: &[GasResistiveSurface<'_>],
    layers: &BTreeMap<u64, GasGrowingCoating<'_>>,
    accommodation: Option<&BTreeMap<u64, ThermalAccommodation>>,
    mut resolver: F,
    powers: Option<&BTreeMap<u64, GasResistivePower<'_>>>,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    duration_s: f64,
    controls: GasGrowingCoatingSolve,
    cancelled: impl FnMut() -> bool,
) -> Result<GasGrowingCoatingCandidate<'catalog>, super::film::GasResistiveFilmError<E>>
where
    F: FnMut(
        &[GasResistiveSurface<'_>],
        &ClosedGasParcel,
        &mut BTreeMap<u64, f64>,
        &mut dyn FnMut() -> bool,
    ) -> Result<(), E>,
{
    super::film::with_resolved_films(&mut resolver, |films| {
        advance_growing(
            surfaces,
            layers,
            accommodation,
            films,
            powers,
            gas,
            catalog,
            duration_s,
            controls,
            cancelled,
        )
    })
}

/// Evolve a growing source-integrated coating with live gas-film transport.
///
/// The resolver runs AFTER the original geometry has rebound each trial's
/// exposed area. It receives the current finite gas and both original body
/// states, including on the second half interval and zero-time admission.
/// Source knots and regime checks can therefore follow evolving pressure,
/// composition and temperatures. The same controller compares complete
/// coupled full/half histories; no film state or rejected energy is committed.
///
/// Resolver output/error semantics are those of
/// [super::film::advance_resistive_chamber_adaptive_resolved].
#[allow(clippy::too_many_arguments)]
pub fn advance_growing_integrated_chamber_resolved<'catalog, F, E>(
    surfaces: &[GasResistiveSurface<'_>],
    layers: &BTreeMap<u64, GasIntegratedGrowingCoating<'_>>,
    accommodation: Option<&BTreeMap<u64, ThermalAccommodation>>,
    mut resolver: F,
    powers: Option<&BTreeMap<u64, GasResistivePower<'_>>>,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    duration_s: f64,
    controls: GasGrowingCoatingSolve,
    cancelled: impl FnMut() -> bool,
) -> Result<GasGrowingCoatingCandidate<'catalog>, super::film::GasResistiveFilmError<E>>
where
    F: FnMut(
        &[GasResistiveSurface<'_>],
        &ClosedGasParcel,
        &mut BTreeMap<u64, f64>,
        &mut dyn FnMut() -> bool,
    ) -> Result<(), E>,
{
    super::film::with_resolved_films(&mut resolver, |films| {
        advance_growing(
            surfaces,
            layers,
            accommodation,
            films,
            powers,
            gas,
            catalog,
            duration_s,
            controls,
            cancelled,
        )
    })
}

#[allow(clippy::too_many_arguments)]
fn advance_growing<'catalog, M: GrowingModel>(
    surfaces: &[GasResistiveSurface<'_>],
    layers: &BTreeMap<u64, M>,
    accommodation: Option<&BTreeMap<u64, ThermalAccommodation>>,
    films: &mut dyn super::film::TrialFilms,
    powers: Option<&BTreeMap<u64, GasResistivePower<'_>>>,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    duration_s: f64,
    controls: GasGrowingCoatingSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<GasGrowingCoatingCandidate<'catalog>, GasResistiveAdvanceError> {
    poll(&mut cancelled)?;
    if surfaces.len() > controls.advance.step.chamber.maximum_surfaces
        || layers.len() != surfaces.len()
        || !controls.maximum_deposit_volume_fraction.is_finite()
        || controls.maximum_deposit_volume_fraction <= 0.0
        || controls.maximum_deposit_volume_fraction > 0.01
        || !controls.absolute_extent_tolerance_m.is_finite()
        || controls.absolute_extent_tolerance_m <= 0.0
        || !controls.absolute_conductance_tolerance_w_k.is_finite()
        || controls.absolute_conductance_tolerance_w_k <= 0.0
    {
        return Err(GasResistiveAdvanceError::Invalid(
            "complete growing-coating models, dilute volume cap and positive geometry tolerances required",
        ));
    }
    if let Some(powers) = powers {
        if powers.len() > surfaces.len()
            || powers
                .iter()
                .any(|(id, power)| !layers.contains_key(id) || power.interface.is_some())
        {
            return Err(GasResistiveAdvanceError::Invalid(
                "growing layer power has an unknown surface or a competing interface owner",
            ));
        }
        for power in powers.values() {
            power.for_interval(duration_s)?;
        }
    }
    let mut initial = surfaces.to_vec();
    initial.sort_by_key(|surface| surface.id);
    if initial.windows(2).any(|pair| pair[0].id == pair[1].id) {
        return Err(GasResistiveAdvanceError::Invalid(
            "duplicate growing coating identity",
        ));
    }
    let mut models = Vec::with_capacity(initial.len());
    let mut initial_geometry = Vec::with_capacity(initial.len());
    let mut volume = Total::default();
    for surface in &mut initial {
        poll(&mut cancelled)?;
        let model = layers
            .get(&surface.id)
            .ok_or(GasResistiveAdvanceError::Invalid(
                "growing coating model roster differs from surfaces",
            ))?;
        let geometry = model.sample(surface, surface.state, &mut cancelled)?;
        volume.add(geometry.front.volume_m3)?;
        surface.state.deposit.area_m2 = geometry.front.area_m2;
        surface.conductance_w_k = geometry.conductance_w_k;
        models.push(model);
        initial_geometry.push(geometry);
    }
    admit_volume(
        volume.value()?,
        gas,
        controls.maximum_deposit_volume_fraction,
    )?;
    let history = super::advance::advance_with(
        &initial,
        gas,
        catalog,
        duration_s,
        controls.advance,
        |rows, gas, dt, cancelled| {
            let mut rebound = rows.to_vec();
            let mut volume = Total::default();
            for ((row, model), original) in rebound.iter_mut().zip(&models).zip(&initial) {
                poll(cancelled)?;
                if row.id != original.id {
                    return Err(GasResistiveAdvanceError::Invalid(
                        "growing coating trial changed canonical identity",
                    ));
                }
                let geometry = model.sample(row, row.state, cancelled)?;
                volume.add(geometry.front.volume_m3)?;
                row.state.deposit.area_m2 = geometry.front.area_m2;
                row.conductance_w_k = geometry.conductance_w_k;
            }
            admit_volume(
                volume.value()?,
                gas,
                controls.maximum_deposit_volume_fraction,
            )?;
            let films = films.resolve(&rebound, gas, cancelled)?;
            let mut candidate = if !M::NONLINEAR_INTERFACE && powers.is_none() {
                // Preserve the original frozen path without allocating any
                // interface objects or source maps inside its physical trials.
                super::advance(
                    &rebound,
                    accommodation,
                    films,
                    None,
                    gas,
                    catalog,
                    dt,
                    controls.advance.step,
                    false,
                    &mut *cancelled,
                )?
            } else {
                // Share the SAME caller checkpoint with the bounded original
                // inner contact solve. Borrows end at every poll; no source or
                // caller state is cloned or advanced on a rejected trial.
                let gate = RefCell::new(&mut *cancelled);
                let checkpoint = || {
                    if (**gate.borrow_mut())() {
                        Err(DepositionError::Cancelled)
                    } else {
                        Ok(())
                    }
                };
                let mut interfaces = Vec::with_capacity(models.len());
                for (model, row) in models.iter().zip(&rebound) {
                    checkpoint()?;
                    interfaces.push(
                        model
                            .interface(row.state.deposit.mass_kg)?
                            .map(|interface| interface.with_checkpoint(&checkpoint)),
                    );
                }
                let mut drives = BTreeMap::new();
                for (row, interface) in rebound.iter().zip(&interfaces) {
                    checkpoint()?;
                    let mut drive = match powers.and_then(|powers| powers.get(&row.id)) {
                        Some(power) => power.for_interval(dt)?,
                        None => GasResistiveDrive {
                            substrate_heat_j: 0.0,
                            deposit_heat_j: 0.0,
                            interface: None,
                        },
                    };
                    drive.interface = interface
                        .as_ref()
                        .map(|law| law as &dyn ResistiveInterfaceLaw);
                    if drive.interface.is_some() || powers.is_some_and(|p| p.contains_key(&row.id))
                    {
                        drives.insert(row.id, drive);
                    }
                }
                super::advance(
                    &rebound,
                    accommodation,
                    films,
                    if drives.is_empty() {
                        None
                    } else {
                        Some(&drives)
                    },
                    gas,
                    catalog,
                    dt,
                    controls.advance.step,
                    false,
                    || (**gate.borrow_mut())(),
                )?
            };
            let mut volume = Total::default();
            for ((endpoint, row), model) in candidate.surfaces.iter_mut().zip(&rebound).zip(&models)
            {
                poll(cancelled)?;
                if endpoint.id != row.id {
                    return Err(GasResistiveAdvanceError::Invalid(
                        "growing coating endpoint changed canonical identity",
                    ));
                }
                let geometry = model.sample(row, endpoint.state, cancelled)?;
                volume.add(geometry.front.volume_m3)?;
                endpoint.state.deposit.area_m2 = geometry.front.area_m2;
                endpoint.exchange.solid.area_m2 = geometry.front.area_m2;
                endpoint.exchange.recession_m = model
                    .geometry()
                    .geometry()
                    .recession_m(row.state.deposit.mass_kg, endpoint.state.deposit.mass_kg)
                    .map_err(GasSurfaceError::from)?;
            }
            admit_volume(
                volume.value()?,
                candidate.gas(),
                controls.maximum_deposit_volume_fraction,
            )?;
            Ok(candidate)
        },
        |before, coarse, fine, relative, cancelled| {
            let mut error: f64 = 0.0;
            for (((row, coarse), fine), model) in before
                .iter()
                .zip(coarse.surfaces())
                .zip(fine.surfaces())
                .zip(&models)
            {
                poll(cancelled)?;
                let initial = model.sample(row, row.state, cancelled)?;
                let coarse = model.sample(row, coarse.state, cancelled)?;
                let fine = model.sample(row, fine.state, cancelled)?;
                error = error
                    .max(ratio(
                        coarse.front.extent_m,
                        fine.front.extent_m,
                        initial.front.extent_m,
                        controls.absolute_extent_tolerance_m,
                        relative,
                    )?)
                    .max(ratio(
                        coarse.conductance_w_k,
                        fine.conductance_w_k,
                        initial.conductance_w_k,
                        controls.absolute_conductance_tolerance_w_k,
                        relative,
                    )?);
            }
            Ok(error)
        },
        &mut cancelled,
    )?;
    if let Some(powers) = powers {
        super::driven::admit_power_history(
            powers,
            &history,
            duration_s,
            controls.advance.step.chamber.surface.gas.energy_tolerance_j,
            &mut cancelled,
        )?;
    }
    let mut geometry = Vec::with_capacity(initial.len());
    for (((row, initial_geometry), endpoint), model) in initial
        .iter()
        .zip(initial_geometry)
        .zip(history.endpoint().surfaces())
        .zip(models)
    {
        poll(&mut cancelled)?;
        geometry.push(GasGrowingCoatingHistory {
            id: row.id,
            initial: initial_geometry,
            final_geometry: model.sample(row, endpoint.state, &mut cancelled)?,
            normal_recession_m: model
                .geometry()
                .geometry()
                .recession_m(row.state.deposit.mass_kg, endpoint.state.deposit.mass_kg)
                .map_err(GasSurfaceError::from)?,
        });
    }
    poll(&mut cancelled)?;
    Ok(GasGrowingCoatingCandidate { history, geometry })
}

#[cfg(test)]
#[path = "growing/tests.rs"]
mod tests;

#[cfg(test)]
#[path = "growing/integrated_tests.rs"]
mod integrated_tests;
