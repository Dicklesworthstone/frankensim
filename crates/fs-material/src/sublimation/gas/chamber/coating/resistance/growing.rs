//! Mass-derived area and thermal resistance in the original finite-gas history.
//!
//! The original phase geometry owns the coating mass, thickness and exposed
//! area. The borrowed conductivity law is queried on the deposit's original
//! specific-enthalpy chart. Each physical interval freezes that conductivity
//! and the series layer/contact conductance at its starting state; the next
//! half step refreshes them. The same two-temperature chamber owns all molecular
//! streams, heat storage and the one gas reconstruction.
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
use super::{
    ClosedGasParcel, GasResistiveSurface, GasSurfaceError, ResistiveSurfaceState,
    ThermalAccommodation, ThermochemicalGasModel, poll,
};
use fs_thermal::enthalpy::conduction::ConductivityLaw;
use fs_thermal::sublimation::coating::resistance::geometry::{
    ConductiveCoatingEndpoint, ConductiveCoatingGeometry,
};
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

fn sample(
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
        let geometry = sample(model, surface, surface.state)?;
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
                let geometry = sample(model, row, row.state)?;
                volume.add(geometry.front.volume_m3)?;
                row.state.deposit.area_m2 = geometry.front.area_m2;
                row.conductance_w_k = geometry.conductance_w_k;
            }
            admit_volume(
                volume.value()?,
                gas,
                controls.maximum_deposit_volume_fraction,
            )?;
            let mut candidate = super::advance(
                &rebound,
                accommodation,
                films,
                gas,
                catalog,
                dt,
                controls.advance.step,
                false,
                &mut *cancelled,
            )?;
            let mut volume = Total::default();
            for ((endpoint, row), model) in candidate.surfaces.iter_mut().zip(&rebound).zip(&models)
            {
                poll(cancelled)?;
                if endpoint.id != row.id {
                    return Err(GasResistiveAdvanceError::Invalid(
                        "growing coating endpoint changed canonical identity",
                    ));
                }
                let geometry = sample(model, row, endpoint.state)?;
                volume.add(geometry.front.volume_m3)?;
                endpoint.state.deposit.area_m2 = geometry.front.area_m2;
                endpoint.exchange.solid.area_m2 = geometry.front.area_m2;
                endpoint.exchange.recession_m = model
                    .geometry
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
                let initial = sample(model, row, row.state)?;
                let coarse = sample(model, row, coarse.state)?;
                let fine = sample(model, row, fine.state)?;
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
            final_geometry: sample(model, row, endpoint.state)?,
            normal_recession_m: model
                .geometry
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
