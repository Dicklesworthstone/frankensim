//! One physical residual shared by numerical trials and final admission.

use std::cell::RefCell;

use fs_qty::Temperature;
use fs_thermal::enthalpy::evaporation::implicit::{
    EvaporatingCellBinding, EvaporatingNetworkBackwardEulerStep,
    EvaporatingNetworkVaporEnthalpyTrial, VaporEnthalpyStepCoefficients,
};
use fs_thermal::evaporation::GasFilm;
use fs_thermochem::{SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K};

use super::super::{admit_liquid_source, chamber::Total, smaller_step};
use super::*;
use crate::thermochemical_gas::energy::parcel::transfer::species::boundary::{
    SpeciesBoundaryRelativeTransfer, admit_catalog,
};

pub(super) struct Prepared<'s, 'm, 'l, 'catalog> {
    initial: &'s DepositionSimulation<'m>,
    step: EvaporatingNetworkBackwardEulerStep<'s, 'm, 'l, dyn EvaporationLaw + 'l>,
    surfaces: BTreeMap<CellId, ImplicitGasSurface<'l>>,
    gas: &'s ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    controls: ImplicitGasControls,
    species: BTreeMap<SpeciesId, Volatile>,
    inventory: BTreeMap<SpeciesId, f64>,
    temperature_scale_j_k: f64,
}

struct Volatile {
    initial_mass_kg: f64,
    gas_constant_j_kg_k: f64,
}

#[derive(Default)]
struct SpeciesSums {
    mass: Total,
    condensation: Total,
    relative_energy: Total,
}

pub(super) struct Trial<'m, T> {
    thermal: EvaporatingNetworkVaporEnthalpyTrial<'m, T>,
    sensible: BTreeMap<CellId, f64>,
    sensible_total_j: f64,
    transfers: Vec<SpeciesBoundaryRelativeTransfer>,
    gas_mass_residual_kg: f64,
    gas_energy_residual_j: f64,
    coupled_mass_residual_kg: f64,
    energy_residual_j: f64,
    residual: Vec<f64>,
}
impl<T> Trial<'_, T> {
    pub(super) fn residual_j(&self) -> Vec<f64> {
        self.residual.clone()
    }
}

impl<'s, 'm, 'l, 'catalog> Prepared<'s, 'm, 'l, 'catalog> {
    #[allow(clippy::too_many_arguments)]
    pub(super) fn new(
        initial: &'s DepositionSimulation<'m>,
        surfaces: &[ImplicitGasSurface<'l>],
        gas: &'s ClosedGasParcel,
        catalog: &'catalog ThermochemicalGasModel,
        target_s: f64,
        controls: ImplicitGasControls,
        resolve_films: &mut impl FnMut(
            &DepositionSimulation<'m>,
            &ClosedGasParcel,
            &mut BTreeMap<CellId, GasFilm>,
        ) -> Result<(), ImplicitGasError>,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, ImplicitGasError> {
        checkpoint()?;
        controls.validate(target_s - initial.time_s())?;
        if initial
            .pending_events()
            .first()
            .is_some_and(|event| event.time_s < target_s)
        {
            return Err(DepositionError::InvalidInput(
                "process pending events or split implicit step at the next event",
            )
            .into());
        }
        if surfaces.len() > initial.cells().len() {
            return Err(ImplicitGasError::Invalid(
                "select one or more active liquid cells",
            ));
        }
        admit_catalog(gas, catalog)?;
        let inventory = gas.species_inventory_kg()?;
        let state = gas.state().state();
        let temperature_scale_j_k = finite(gas.mass_kg() * state.specific_heat_cp / state.gamma)?;
        if !positive(temperature_scale_j_k)
            || !(controls.phase.gas.gas.minimum_temperature_k
                ..=controls.phase.gas.gas.maximum_temperature_k)
                .contains(&state.temperature)
        {
            return Err(ImplicitGasError::Invalid(
                "initial gas must lie inside its admitted caloric bracket",
            ));
        }
        let mut ordered = BTreeMap::new();
        let mut species = BTreeMap::new();
        for surface in surfaces {
            checkpoint()?;
            if !nonnegative(surface.mass_transfer_coefficient_m_s)
                || !nonnegative(surface.heat_transfer_coefficient_w_m2_k)
                || ordered.insert(surface.cell_id, *surface).is_some()
            {
                return Err(ImplicitGasError::Invalid(
                    "unique liquid cells and finite nonnegative film coefficients required",
                ));
            }
            let cell = initial
                .cells()
                .get(&surface.cell_id)
                .ok_or(ImplicitGasError::Invalid(
                    "selected liquid cell must be active",
                ))?;
            let source = admit_liquid_source(
                surface.law,
                surface.binding,
                catalog,
                initial.temperature_k(surface.cell_id)?,
                cell.specific_h_j_kg,
                controls.phase,
            )?;
            species
                .entry(surface.binding.species.clone())
                .or_insert(Volatile {
                    initial_mass_kg: inventory
                        .get(&surface.binding.species)
                        .copied()
                        .unwrap_or(0.0),
                    gas_constant_j_kg_k: finite(
                        UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / source.molecular_mass,
                    )?,
                });
        }
        let initial_films = films(&ordered, initial, gas, &species, resolve_films, checkpoint)?;
        let bindings: Vec<_> = ordered
            .iter()
            .map(|(&id, surface)| EvaporatingCellBinding {
                id,
                law: surface.law,
                film: initial_films[&id],
                pure_liquid_enthalpy_range_j_kg: surface.pure_liquid_enthalpy_range_j_kg,
            })
            .collect();
        let step = EvaporatingNetworkBackwardEulerStep::new(
            initial,
            &bindings,
            target_s,
            &mut *checkpoint,
        )?;
        checkpoint()?;
        Ok(Self {
            initial,
            step,
            surfaces: ordered,
            gas,
            catalog,
            controls,
            species,
            inventory,
            temperature_scale_j_k,
        })
    }

    pub(super) fn dimension(&self) -> usize {
        self.thermal_dimension() + self.species.len() + usize::from(!self.species.is_empty())
    }

    pub(super) fn thermal_dimension(&self) -> usize {
        self.step.cell_ids().len()
    }

    /// Same representability projection as native heat, extended only to the
    /// species-mass and gas-temperature coordinates. Residuals retain every mass/J.
    pub(super) fn project(&self, coordinates: &[f64]) -> Vec<f64> {
        let n = self.thermal_dimension();
        if coordinates.len() != self.dimension() {
            return coordinates.to_vec();
        }
        let mut projected = self.step.resolvable_increments_j(&coordinates[..n]);
        if self.species.is_empty() {
            return projected;
        }
        for (index, species) in self.species.values().enumerate() {
            let coordinate = coordinates[n + index];
            projected.push(
                if species.initial_mass_kg + coordinate / self.controls.mass_coordinate_scale_j_kg
                    == species.initial_mass_kg
                {
                    0.0
                } else {
                    coordinate
                },
            );
        }
        let t0 = self.gas.state().state().temperature;
        let temperature_coordinate = coordinates[n + self.species.len()];
        projected.push(
            if t0 + temperature_coordinate / self.temperature_scale_j_k == t0 {
                0.0
            } else {
                temperature_coordinate
            },
        );
        projected
    }

    /// A source-owned mass predictor puts initially dry evaporation inside the
    /// feasible inventory domain. No trace gas, clipped donor or heat-only
    /// temperature endpoint is introduced. This trial spends the common budget.
    pub(super) fn initial_guess<T>(&self, predictor: &Trial<'_, T>) -> Vec<f64> {
        let mut guess = vec![0.0; self.dimension()];
        for (index, transfer) in predictor.transfers.iter().enumerate() {
            guess[self.thermal_dimension() + index] =
                transfer.mass_into_gas_kg * self.controls.mass_coordinate_scale_j_kg;
        }
        self.project(&guess)
    }

    pub(super) fn evaluate<T>(
        &self,
        coordinates: &[f64],
        resolve: &mut impl FnMut(
            &DepositionSimulation<'m>,
            &BTreeMap<CellId, LiquidSurfaceState>,
            &ClosedGasParcel,
        )
            -> Result<BoundaryExchangeStepCoefficients<T>, ImplicitGasError>,
        resolve_films: &mut impl FnMut(
            &DepositionSimulation<'m>,
            &ClosedGasParcel,
            &mut BTreeMap<CellId, GasFilm>,
        ) -> Result<(), ImplicitGasError>,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Trial<'m, T>, ImplicitGasError> {
        checkpoint()?;
        let n = self.thermal_dimension();
        if coordinates.len() != self.dimension() || !coordinates.iter().all(|x| x.is_finite()) {
            return Err(ImplicitGasError::Invalid(
                "finite coupled coordinate vector required",
            ));
        }
        if self.species.is_empty() {
            return self.evaluate_at_gas(coordinates, self.gas, resolve, resolve_films, checkpoint);
        }
        let t = finite(
            self.gas.state().state().temperature
                + coordinates[n + self.species.len()] / self.temperature_scale_j_k,
        )?;
        let mut target = None;
        for (index, (id, species)) in self.species.iter().enumerate() {
            let mass = finite(
                species.initial_mass_kg
                    + coordinates[n + index] / self.controls.mass_coordinate_scale_j_kg,
            )?;
            if mass < 0.0 {
                return Err(ImplicitGasError::Invalid(
                    "negative proposed vapor inventory",
                ));
            }
            if mass != species.initial_mass_kg {
                target
                    .get_or_insert_with(|| self.inventory.clone())
                    .insert(id.clone(), mass);
            }
        }
        let gas = if let Some(target) = target {
            self.gas
                .at_species_inventory_temperature(self.catalog, &target, t)?
        } else if t == self.gas.state().state().temperature {
            self.gas.clone()
        } else {
            self.gas
                .at_model_mass_temperature(self.gas.model().clone(), self.gas.mass_kg(), t)?
        };
        self.evaluate_at_gas(&coordinates[..n], &gas, resolve, resolve_films, checkpoint)
    }

    pub(super) fn evaluate_at_gas<T>(
        &self,
        thermal_coordinates: &[f64],
        gas: &ClosedGasParcel,
        resolve: &mut impl FnMut(
            &DepositionSimulation<'m>,
            &BTreeMap<CellId, LiquidSurfaceState>,
            &ClosedGasParcel,
        )
            -> Result<BoundaryExchangeStepCoefficients<T>, ImplicitGasError>,
        resolve_films: &mut impl FnMut(
            &DepositionSimulation<'m>,
            &ClosedGasParcel,
            &mut BTreeMap<CellId, GasFilm>,
        ) -> Result<(), ImplicitGasError>,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Trial<'m, T>, ImplicitGasError> {
        checkpoint()?;
        let t0 = self.gas.state().state().temperature;
        let tg = gas.state().state().temperature;
        let dt = self.step.duration_s();
        if !(self.controls.phase.gas.gas.minimum_temperature_k
            ..=self.controls.phase.gas.gas.maximum_temperature_k)
            .contains(&tg)
        {
            return Err(ImplicitGasError::Invalid(
                "gas trial left its admitted caloric bracket",
            ));
        }
        if (tg - t0).abs() > self.controls.phase.maximum_temperature_increment_k {
            return Err(
                smaller_step(dt, "implicit gas temperature departure exceeds budget").into(),
            );
        }
        let poll = RefCell::new(checkpoint);
        let mut sensible = BTreeMap::new();
        let mut species_sums: BTreeMap<_, _> = self
            .species
            .keys()
            .map(|id| (id.clone(), SpeciesSums::default()))
            .collect();
        let mut heat_total = Total::default();
        let thermal = self.step.evaluate_with_vapor_enthalpy(
            thermal_coordinates,
            |thermal| {
                films(
                    &self.surfaces,
                    thermal,
                    gas,
                    &self.species,
                    resolve_films,
                    &mut || (poll.borrow_mut())(),
                )
            },
            |thermal, liquids, _| {
                (poll.borrow_mut())()?;
                let mut heat = resolve(thermal, liquids, gas)?;
                let mut excess = BTreeMap::new();
                for (&id, surface) in &self.surfaces {
                    (poll.borrow_mut())()?;
                    let species = &surface.binding.species;
                    let liquid = liquids[&id];
                    let old = &self.initial.cells()[&id];
                    let tl = thermal.temperature_k(id)?;
                    if (tl - self.initial.temperature_k(id)?).abs()
                        > self.controls.phase.maximum_temperature_increment_k
                    {
                        return Err(smaller_step(
                            dt,
                            "implicit liquid temperature departure exceeds budget",
                        )
                        .into());
                    }
                    let source = admit_liquid_source(
                        surface.law,
                        surface.binding,
                        self.catalog,
                        tl,
                        liquid.specific_enthalpy_j_kg,
                        self.controls.phase,
                    )?;
                    let dm = finite(old.mass_kg - liquid.mass_kg)?;
                    if dm > self.controls.phase.maximum_inventory_fraction * old.mass_kg {
                        return Err(smaller_step(
                            dt,
                            "implicit evaporation exceeds its original liquid donor budget",
                        )
                        .into());
                    }
                    let sums = species_sums
                        .get_mut(species)
                        .expect("prepared volatile species");
                    if dm < 0.0 {
                        sums.condensation.add(-dm)?;
                    }
                    sums.mass.add(dm)?;
                    let upstream = if dm < 0.0 { tg } else { tl };
                    let incoming_sensible = if dm < 0.0 {
                        self.catalog
                            .species_specific_enthalpy_change_j_kg(species, tl, tg)?
                    } else {
                        0.0
                    };
                    excess.insert(id, finite(source.alignment + incoming_sensible)?);
                    sums.relative_energy.add(finite(
                        dm * self
                            .catalog
                            .species_specific_enthalpy_change_j_kg(species, t0, upstream)?,
                    )?)?;
                    let power = finite(
                        surface.heat_transfer_coefficient_w_m2_k * liquid.area_m2 * (tg - tl),
                    )?;
                    let q = finite(dt * power)?;
                    sensible.insert(id, q);
                    heat_total.add(q)?;
                    let slot = heat.transport.heat.bulk.sources_w.entry(id).or_default();
                    *slot = finite(*slot + power)?;
                }
                for (id, sums) in &species_sums {
                    if sums.condensation.value()?
                        > self.controls.phase.maximum_inventory_fraction
                            * self.species[id].initial_mass_kg
                    {
                        return Err(smaller_step(
                            dt,
                            "joint condensation exceeds the original gas vapor donor budget",
                        )
                        .into());
                    }
                }
                Ok::<_, ImplicitGasError>(VaporEnthalpyStepCoefficients {
                    heat,
                    vapor_enthalpy_excess_j_kg: excess,
                })
            },
            || (poll.borrow_mut())(),
        )?;
        let mut transfers = Vec::with_capacity(self.species.len());
        let mut mass = Total::default();
        let mut relative = Total::default();
        let mut pressure_flow = Total::default();
        // Retain zero-net rows: opposite streams can exchange energy without
        // changing this species' inventory. Species order also owns coordinates.
        for (id, sums) in species_sums {
            let dm = sums.mass.value()?;
            let relative_energy = sums.relative_energy.value()?;
            mass.add(dm)?;
            relative.add(relative_energy)?;
            pressure_flow.add(finite(dm * self.species[&id].gas_constant_j_kg_k * t0)?)?;
            transfers.push(SpeciesBoundaryRelativeTransfer {
                species: id,
                mass_into_gas_kg: dm,
                relative_energy_into_gas_j: relative_energy,
            });
        }
        let transferred_mass_kg = mass.value()?;
        let relative_energy_into_gas_j = relative.value()?;
        let sensible_total_j = heat_total.value()?;
        let actual = gas.species_inventory_kg()?;
        let total_mass_residual =
            finite((gas.mass_kg() - self.gas.mass_kg()) - transferred_mass_kg)?;
        let mut gas_mass_residual_kg = total_mass_residual.abs();
        let mut residual = thermal.residual_j().to_vec();
        for transfer in &transfers {
            let id = &transfer.species;
            let vapor_residual = finite(
                (actual.get(id).copied().unwrap_or(0.0) - self.species[id].initial_mass_kg)
                    - transfer.mass_into_gas_kg,
            )?;
            gas_mass_residual_kg = gas_mass_residual_kg.max(vapor_residual.abs());
            residual.push(finite(
                vapor_residual * self.controls.mass_coordinate_scale_j_kg,
            )?);
        }
        for source in self.catalog.thermodynamics().components() {
            (poll.borrow_mut())()?;
            let id = source.species();
            if !self.species.contains_key(id) {
                let defect = finite(
                    actual.get(id).copied().unwrap_or(0.0)
                        - self.inventory.get(id).copied().unwrap_or(0.0),
                )?;
                gas_mass_residual_kg = gas_mass_residual_kg.max(defect.abs());
            }
        }
        let increment = gas
            .model()
            .thermodynamics()
            .specific_internal_energy_change(Temperature::new(t0), Temperature::new(tg))
            .map_err(|error| {
                GasEnergyError::Material(crate::MaterialError::Parameters {
                    what: format!("implicit gas original caloric increment: {error}"),
                })
            })?
            .value();
        let gas_energy_residual_j = sum([
            finite(gas.mass_kg() * increment)?,
            -pressure_flow.value()?,
            -relative_energy_into_gas_j,
            sensible_total_j,
        ])?;
        let coupled_mass_residual_kg = sum([
            thermal.audit().thermal.mass_balance_residual_kg,
            total_mass_residual,
        ])?;
        let energy_residual_j = sum([
            thermal.audit().thermal.energy_residual_j,
            gas_energy_residual_j,
        ])?;
        if !self.species.is_empty() {
            residual.push(gas_energy_residual_j);
        }
        (poll.borrow_mut())()?;
        Ok(Trial {
            thermal,
            sensible,
            sensible_total_j,
            transfers,
            gas_mass_residual_kg,
            gas_energy_residual_j,
            coupled_mass_residual_kg,
            energy_residual_j,
            residual,
        })
    }

    pub(super) fn check_balance<T>(&self, trial: &Trial<'_, T>) -> Result<(), ImplicitGasError> {
        let mass_tolerance = self.controls.phase.gas.mass_tolerance_kg;
        if trial.gas_mass_residual_kg > mass_tolerance {
            return Err(ImplicitGasError::GasMassBalance {
                residual_kg: trial.gas_mass_residual_kg,
                tolerance_kg: mass_tolerance,
            });
        }
        let coupled_mass_tolerance = mass_tolerance.min(self.controls.thermal.mass_tolerance_kg);
        if trial.coupled_mass_residual_kg.abs() > coupled_mass_tolerance {
            return Err(ImplicitGasError::GasMassBalance {
                residual_kg: trial.coupled_mass_residual_kg,
                tolerance_kg: coupled_mass_tolerance,
            });
        }
        if trial.gas_energy_residual_j.abs() > self.controls.phase.gas.energy_tolerance_j {
            return Err(ImplicitGasError::GasEnergyBalance {
                residual_j: trial.gas_energy_residual_j,
                tolerance_j: self.controls.phase.gas.energy_tolerance_j,
            });
        }
        if trial.energy_residual_j.abs() > self.controls.thermal.total_energy_tolerance_j {
            return Err(ImplicitGasError::CoupledEnergyBalance {
                residual_j: trial.energy_residual_j,
                tolerance_j: self.controls.thermal.total_energy_tolerance_j,
            });
        }
        Ok(())
    }

    pub(super) fn reconstruct<T>(
        &self,
        trial: &Trial<'_, T>,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<GasSpeciesBoundaryCandidate<'catalog>, ImplicitGasError> {
        checkpoint()?;
        let mut interrupted = None;
        let candidate = self.gas.exchange_species_boundary_relative_with_energy(
            self.catalog,
            &trial.transfers,
            |_| Ok(-trial.sensible_total_j),
            self.controls.phase.gas,
            || match checkpoint() {
                Ok(()) => false,
                Err(error) => {
                    interrupted = Some(error);
                    true
                }
            },
        );
        if let Some(error) = interrupted {
            return Err(error.into());
        }
        let candidate = candidate?;
        checkpoint()?;
        Ok(candidate)
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn accept<T>(
        &self,
        trial: Trial<'m, T>,
        gas: GasSpeciesBoundaryCandidate<'catalog>,
        newton_iterations: usize,
        residual_evaluations: usize,
        maximum_accepted_residual_norm_j: f64,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<ImplicitGasCandidate<'m, 'catalog, T>, ImplicitGasError> {
        checkpoint()?;
        self.check_balance(&trial)?;
        let external_heat_j = sum([
            trial.thermal.audit().thermal.external_heat_j,
            -trial.sensible_total_j,
        ])?;
        let (thermal, payload, audit) = trial.thermal.accept(
            self.controls.thermal.cell_energy_tolerance_j,
            self.controls.thermal.total_energy_tolerance_j,
            self.controls.thermal.mass_tolerance_kg,
            &mut *checkpoint,
        )?;
        let report = ImplicitGasReport {
            thermal: audit,
            sensible_heat_into_liquid_j: trial.sensible,
            species_transferred_mass_kg: trial
                .transfers
                .into_iter()
                .map(|transfer| (transfer.species, transfer.mass_into_gas_kg))
                .collect(),
            gas_mass_residual_kg: trial.gas_mass_residual_kg,
            gas_energy_residual_j: trial.gas_energy_residual_j,
            mass_residual_kg: trial.coupled_mass_residual_kg,
            energy_residual_j: trial.energy_residual_j,
            external_heat_j,
            newton_iterations,
            residual_evaluations,
            gas_iterations: gas.step().iterations,
            maximum_accepted_residual_norm_j,
        };
        checkpoint()?;
        Ok(ImplicitGasCandidate {
            thermal,
            gas,
            payload,
            report,
        })
    }
}

fn films<'m>(
    surfaces: &BTreeMap<CellId, ImplicitGasSurface<'_>>,
    thermal: &DepositionSimulation<'m>,
    gas: &ClosedGasParcel,
    species: &BTreeMap<SpeciesId, Volatile>,
    resolve: &mut impl FnMut(
        &DepositionSimulation<'m>,
        &ClosedGasParcel,
        &mut BTreeMap<CellId, GasFilm>,
    ) -> Result<(), ImplicitGasError>,
    checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
) -> Result<BTreeMap<CellId, GasFilm>, ImplicitGasError> {
    checkpoint()?;
    let state = gas.state().state();
    let mut partial_pressures = BTreeMap::new();
    for id in species.keys() {
        checkpoint()?;
        partial_pressures.insert(id.clone(), gas.species_partial_pressure_pa(id)?);
    }
    let mut films: BTreeMap<_, _> = surfaces
        .iter()
        .map(|(&id, surface)| {
            (
                id,
                GasFilm {
                    total_pressure_pa: state.pressure,
                    vapor_partial_pressure_pa: partial_pressures[&surface.binding.species],
                    temperature_k: state.temperature,
                    mass_transfer_coefficient_m_s: surface.mass_transfer_coefficient_m_s,
                },
            )
        })
        .collect();
    resolve(thermal, gas, &mut films)?;
    checkpoint()?;
    if !films.keys().eq(surfaces.keys()) {
        return Err(ImplicitGasError::Invalid(
            "film resolver must preserve the exact active liquid cell roster",
        ));
    }
    for (&id, film) in &films {
        checkpoint()?;
        let species = &surfaces[&id].binding.species;
        if film.total_pressure_pa != state.pressure
            || film.vapor_partial_pressure_pa != partial_pressures[species]
        {
            return Err(ImplicitGasError::Invalid(
                "film resolver must preserve the original gas total and volatile partial pressures",
            ));
        }
        film.validate().map_err(GasEvaporationError::from)?;
    }
    Ok(films)
}

fn finite(value: f64) -> Result<f64, ImplicitGasError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(ImplicitGasError::Invalid(
            "implicit finite-gas arithmetic is unrepresentable",
        ))
    }
}
fn sum(values: impl IntoIterator<Item = f64>) -> Result<f64, ImplicitGasError> {
    let mut total = Total::default();
    for value in values {
        total.add(value)?;
    }
    Ok(total.value()?)
}
