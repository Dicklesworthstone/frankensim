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
    species: Option<SpeciesId>,
    inventory: BTreeMap<SpeciesId, f64>,
    old_vapor_mass_kg: f64,
    vapor_gas_constant_j_kg_k: f64,
    temperature_scale_j_k: f64,
}

pub(super) struct Trial<'m, T> {
    thermal: EvaporatingNetworkVaporEnthalpyTrial<'m, T>,
    sensible: BTreeMap<CellId, f64>,
    sensible_total_j: f64,
    transferred_mass_kg: f64,
    relative_energy_into_gas_j: f64,
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
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, ImplicitGasError> {
        checkpoint()?;
        controls.validate(target_s - initial.time_s())?;
        if surfaces.len() > initial.cells().len() {
            return Err(ImplicitGasError::Invalid(
                "select one or more active liquid cells",
            ));
        }
        admit_catalog(gas, catalog)?;
        let species = surfaces
            .first()
            .map(|surface| surface.binding.species.clone());
        let inventory = gas.species_inventory_kg()?;
        let old_vapor_mass_kg = species
            .as_ref()
            .and_then(|species| inventory.get(species))
            .copied()
            .unwrap_or(0.0);
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
        let mut molecular_mass = 0.0;
        for surface in surfaces {
            checkpoint()?;
            if Some(&surface.binding.species) != species.as_ref()
                || !nonnegative(surface.mass_transfer_coefficient_m_s)
                || !nonnegative(surface.heat_transfer_coefficient_w_m2_k)
                || ordered.insert(surface.cell_id, *surface).is_some()
            {
                return Err(ImplicitGasError::Invalid(
                    "unique liquid cells, one volatile identity and finite nonnegative film coefficients required",
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
            molecular_mass = source.molecular_mass;
        }
        let initial_films = films(&ordered, gas, species.as_ref())?;
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
        let vapor_gas_constant_j_kg_k = if species.is_some() {
            finite(UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / molecular_mass)?
        } else {
            0.0
        };
        Ok(Self {
            initial,
            step,
            surfaces: ordered,
            gas,
            catalog,
            controls,
            species,
            inventory,
            old_vapor_mass_kg,
            temperature_scale_j_k,
            vapor_gas_constant_j_kg_k,
        })
    }

    pub(super) fn dimension(&self) -> usize {
        self.thermal_dimension() + if self.species.is_some() { 2 } else { 0 }
    }

    pub(super) fn thermal_dimension(&self) -> usize {
        self.step.cell_ids().len()
    }

    /// Same representability projection as native heat, extended only to the
    /// two numerical gas coordinates. Residuals still retain every mass/J.
    pub(super) fn project(&self, coordinates: &[f64]) -> Vec<f64> {
        let n = self.thermal_dimension();
        if coordinates.len() != self.dimension() {
            return coordinates.to_vec();
        }
        let mut projected = self.step.resolvable_increments_j(&coordinates[..n]);
        if self.species.is_none() {
            return projected;
        }
        projected.push(
            if self.old_vapor_mass_kg + coordinates[n] / self.controls.mass_coordinate_scale_j_kg
                == self.old_vapor_mass_kg
            {
                0.0
            } else {
                coordinates[n]
            },
        );
        let t0 = self.gas.state().state().temperature;
        projected.push(
            if t0 + coordinates[n + 1] / self.temperature_scale_j_k == t0 {
                0.0
            } else {
                coordinates[n + 1]
            },
        );
        projected
    }

    /// A source-owned mass predictor puts initially dry evaporation inside the
    /// feasible inventory domain. No trace gas, clipped donor or heat-only
    /// temperature endpoint is introduced. This trial spends the common budget.
    pub(super) fn initial_guess<T>(&self, predictor: &Trial<'_, T>) -> Vec<f64> {
        let mut guess = vec![0.0; self.dimension()];
        if self.species.is_some() {
            guess[self.thermal_dimension()] =
                predictor.transferred_mass_kg * self.controls.mass_coordinate_scale_j_kg;
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
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Trial<'m, T>, ImplicitGasError> {
        checkpoint()?;
        let n = self.thermal_dimension();
        if coordinates.len() != self.dimension() || !coordinates.iter().all(|x| x.is_finite()) {
            return Err(ImplicitGasError::Invalid(
                "finite coupled coordinate vector required",
            ));
        }
        let Some(species) = &self.species else {
            return self.evaluate_at_gas(coordinates, self.gas, resolve, checkpoint);
        };
        let t = finite(
            self.gas.state().state().temperature + coordinates[n + 1] / self.temperature_scale_j_k,
        )?;
        let mass = finite(
            self.old_vapor_mass_kg + coordinates[n] / self.controls.mass_coordinate_scale_j_kg,
        )?;
        if mass < 0.0 {
            return Err(ImplicitGasError::Invalid(
                "negative proposed vapor inventory",
            ));
        }
        let gas = if mass == self.old_vapor_mass_kg {
            if t == self.gas.state().state().temperature {
                self.gas.clone()
            } else {
                self.gas.at_model_mass_temperature(
                    self.gas.model().clone(),
                    self.gas.mass_kg(),
                    t,
                )?
            }
        } else {
            let mut target = self.inventory.clone();
            target.insert(species.clone(), mass);
            self.gas
                .at_species_inventory_temperature(self.catalog, &target, t)?
        };
        self.evaluate_at_gas(&coordinates[..n], &gas, resolve, checkpoint)
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
        let mut mass = Total::default();
        let mut condensation = Total::default();
        let mut relative = Total::default();
        let mut heat_total = Total::default();
        let thermal = self.step.evaluate_with_vapor_enthalpy(
            thermal_coordinates,
            |_| films(&self.surfaces, gas, self.species.as_ref()),
            |thermal, liquids, _| {
                (poll.borrow_mut())()?;
                let mut heat = resolve(thermal, liquids, gas)?;
                let mut excess = BTreeMap::new();
                for (&id, surface) in &self.surfaces {
                    (poll.borrow_mut())()?;
                    let species = self.species.as_ref().ok_or(ImplicitGasError::Invalid(
                        "active liquid surface requires its volatile identity",
                    ))?;
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
                    if dm < 0.0 {
                        condensation.add(-dm)?;
                    }
                    mass.add(dm)?;
                    let upstream = if dm < 0.0 { tg } else { tl };
                    let incoming_sensible = if dm < 0.0 {
                        self.catalog
                            .species_specific_enthalpy_change_j_kg(species, tl, tg)?
                    } else {
                        0.0
                    };
                    excess.insert(id, finite(source.alignment + incoming_sensible)?);
                    relative.add(finite(
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
                if condensation.value()?
                    > self.controls.phase.maximum_inventory_fraction * self.old_vapor_mass_kg
                {
                    return Err(smaller_step(
                        dt,
                        "joint condensation exceeds the original gas vapor donor budget",
                    )
                    .into());
                }
                Ok::<_, ImplicitGasError>(VaporEnthalpyStepCoefficients {
                    heat,
                    vapor_enthalpy_excess_j_kg: excess,
                })
            },
            || (poll.borrow_mut())(),
        )?;
        let transferred_mass_kg = mass.value()?;
        let relative_energy_into_gas_j = relative.value()?;
        let sensible_total_j = heat_total.value()?;
        let actual = gas.species_inventory_kg()?;
        let vapor_residual = if let Some(species) = &self.species {
            finite(
                (actual.get(species).copied().unwrap_or(0.0) - self.old_vapor_mass_kg)
                    - transferred_mass_kg,
            )?
        } else {
            0.0
        };
        let total_mass_residual =
            finite((gas.mass_kg() - self.gas.mass_kg()) - transferred_mass_kg)?;
        let mut gas_mass_residual_kg = vapor_residual.abs().max(total_mass_residual.abs());
        for source in self.catalog.thermodynamics().components() {
            (poll.borrow_mut())()?;
            let id = source.species();
            if Some(id) != self.species.as_ref() {
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
            -finite(transferred_mass_kg * self.vapor_gas_constant_j_kg_k * t0)?,
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
        let mut residual = thermal.residual_j().to_vec();
        if self.species.is_some() {
            residual.push(finite(
                vapor_residual * self.controls.mass_coordinate_scale_j_kg,
            )?);
            residual.push(gas_energy_residual_j);
        }
        (poll.borrow_mut())()?;
        Ok(Trial {
            thermal,
            sensible,
            sensible_total_j,
            transferred_mass_kg,
            relative_energy_into_gas_j,
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
        let transfer = self
            .species
            .as_ref()
            .map(|species| SpeciesBoundaryRelativeTransfer {
                species: species.clone(),
                mass_into_gas_kg: trial.transferred_mass_kg,
                relative_energy_into_gas_j: trial.relative_energy_into_gas_j,
            });
        let candidate = self.gas.exchange_species_boundary_relative_with_energy(
            self.catalog,
            transfer.as_slice(),
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

fn films(
    surfaces: &BTreeMap<CellId, ImplicitGasSurface<'_>>,
    gas: &ClosedGasParcel,
    species: Option<&SpeciesId>,
) -> Result<BTreeMap<CellId, GasFilm>, ImplicitGasError> {
    let Some(species) = species else {
        return Ok(BTreeMap::new());
    };
    let state = gas.state().state();
    let vapor_partial_pressure_pa = gas.species_partial_pressure_pa(species)?;
    Ok(surfaces
        .iter()
        .map(|(&id, surface)| {
            (
                id,
                GasFilm {
                    total_pressure_pa: state.pressure,
                    vapor_partial_pressure_pa,
                    temperature_k: state.temperature,
                    mass_transfer_coefficient_m_s: surface.mass_transfer_coefficient_m_s,
                },
            )
        })
        .collect())
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
