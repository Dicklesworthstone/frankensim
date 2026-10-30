//! Thermal-owner admission for externally solved solid/gas exchange.
//!
//! The external owner supplies a pure-solid enthalpy interval and the complete
//! expected cell. This owner retains identities, caloric laws, reference
//! geometry, contacts, clock and schedule; only accepted mass and h can change.
//! Phase identity is caller-admitted: a temperature or an interval alone cannot
//! prove that an arbitrary material chart describes a pure solid.

use std::collections::BTreeSet;

use crate::enthalpy::{
    DepositedCell, DepositionBudget, DepositionError, DepositionSimulation, ThermalStorageBasis,
    TimedEvent,
};
use crate::sublimation::SolidSurfaceState;

/// One prepared endpoint bound to the unchanged original thermal cell.
#[derive(Clone, Debug, PartialEq)]
pub struct SolidCellExchange {
    /// Full post-heat cell, including its identities, geometry and boundary.
    /// Any stale mass, enthalpy or other field refuses the complete batch.
    pub expected: DepositedCell,
    /// Exposure owned by the preceding physical interval. It may differ from
    /// current exposure when a contact event occurred at that interval's end.
    pub interval_exposed_area_m2: f64,
    /// Accepted mass/h; its area must equal the supplied interval exposure.
    pub replacement: SolidSurfaceState,
    /// Inclusive caller-admitted pure-solid branch, in the ORIGINAL h datum.
    pub pure_solid_enthalpy_range_j_kg: [f64; 2],
}

impl DepositionSimulation<'_> {
    /// Unconsumed original events, in their authored equal-time order.
    #[must_use]
    pub fn pending_events(&self) -> &[TimedEvent] {
        &self.schedule[self.cursor..]
    }

    /// Original per-advance work ceilings; a coupled owner must not renew them
    /// each time it splits one requested advance into physical intervals.
    #[must_use]
    pub const fn work_budget(&self) -> DepositionBudget {
        self.budget
    }

    /// Evaluate this active cell's ORIGINAL caloric chart at a proposed h.
    /// No alternate heat capacity, extrapolation or state mutation is supplied.
    pub fn temperature_at_specific_enthalpy_k(
        &self,
        id: crate::enthalpy::CellId,
        specific_h_j_kg: f64,
    ) -> Result<f64, DepositionError> {
        if !specific_h_j_kg.is_finite() {
            return Err(DepositionError::MaterialDomain);
        }
        let cell = self.cells.get(&id).ok_or(DepositionError::InvalidInput(
            "external solid exchange requires an active cell",
        ))?;
        let t = self.materials[&cell.material_id].temperature_k(specific_h_j_kg)?;
        if !crate::enthalpy::positive(t) {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(t)
    }

    /// Atomically accept prepared endpoints on the existing thermal storage.
    ///
    /// The caller owns interval timing, the physical pure-solid admission and
    /// joint gas/solid conservation. Process heat and endpoint events before
    /// capturing the expected cell, but retain the START interval area.
    /// End events cannot retroactively change exposure during that interval.
    ///
    /// Every expected cell, both caloric endpoints, original-geometry area,
    /// positive remaining inventory and the final cancellation gate are checked
    /// before any write. No mesh, contact, material, boundary or clock is changed.
    /// A surrounding owner can mutate one private trial without cloning the
    /// entire graph again for each phase substep.
    pub fn apply_solid_exchanges(
        &mut self,
        exchanges: &[SolidCellExchange],
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), DepositionError> {
        checkpoint()?;
        if exchanges.len() > self.budget.max_cells {
            return Err(DepositionError::Budget("external solid exchange cells"));
        }
        let mut ids = BTreeSet::new();
        for (index, exchange) in exchanges.iter().enumerate() {
            if index % 256 == 0 {
                checkpoint()?;
            }
            let id = exchange.expected.id;
            let cell = self.cells.get(&id).ok_or(DepositionError::InvalidInput(
                "external solid exchange requires an active cell",
            ))?;
            if !ids.insert(id) || cell != &exchange.expected {
                return Err(DepositionError::InvalidInput(
                    "duplicate or stale external solid exchange",
                ));
            }
            if self.material_storage_basis(cell.material_id)? != ThermalStorageBasis::Enthalpy {
                return Err(DepositionError::InvalidInput(
                    "external solid exchange requires enthalpy storage and an explicit open-mass model",
                ));
            }
            let [lo, hi] = exchange.pure_solid_enthalpy_range_j_kg;
            let area = exchange.interval_exposed_area_m2;
            let next = exchange.replacement;
            if !lo.is_finite()
                || !hi.is_finite()
                || lo > hi
                || !(lo..=hi).contains(&cell.specific_h_j_kg)
                || !(lo..=hi).contains(&next.specific_enthalpy_j_kg)
                || !crate::enthalpy::positive(next.mass_kg)
                || !crate::enthalpy::nonnegative(area)
                || area > cell.geometry.surface_area_m2
                || next.area_m2 != area
            {
                return Err(DepositionError::InvalidInput(
                    "external solid endpoint leaves admitted branch, inventory or reference area",
                ));
            }
            self.temperature_at_specific_enthalpy_k(id, cell.specific_h_j_kg)?;
            self.temperature_at_specific_enthalpy_k(id, next.specific_enthalpy_j_kg)?;
        }
        checkpoint()?;
        // No fallible work or callbacks after this gate.
        for exchange in exchanges {
            let cell = self
                .cells
                .get_mut(&exchange.expected.id)
                .expect("all external solid endpoints validated");
            cell.mass_kg = exchange.replacement.mass_kg;
            cell.specific_h_j_kg = exchange.replacement.specific_enthalpy_j_kg;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::enthalpy::{
        Boundary, CellGeometry, CellId, Event, MaterialId, Reservoir, StepConfig, ThermalLaw,
    };

    struct Caloric;
    impl ThermalLaw for Caloric {
        fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
            let t = 250.0 + h / 2000.0;
            if (200.0..=300.0).contains(&t) {
                Ok(t)
            } else {
                Err(DepositionError::MaterialDomain)
            }
        }
        fn maximum_temperature_slope(&self) -> f64 {
            0.0005
        }
    }
    static CALORIC: Caloric = Caloric;

    fn state() -> DepositionSimulation<'static> {
        let events = [1, 2].map(|id| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(id),
                segment_id: id,
                body_id: id,
                material_id: MaterialId(1),
                mass_kg: 0.01,
                specific_h_j_kg: 0.0,
                geometry: CellGeometry {
                    start_m: [id as f64, 0.0, 0.0],
                    end_m: [id as f64, 0.01, 0.0],
                    volume_m3: 1e-5,
                    surface_area_m2: 0.002,
                },
                boundary: Boundary {
                    ambient_k: 250.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        });
        DepositionSimulation::new(
            events.to_vec(),
            &[(MaterialId(1), &CALORIC)],
            DepositionBudget {
                max_cells: 2,
                max_links: 0,
                max_events: 2,
                max_substeps: 4,
            },
        )
        .unwrap()
        .advance_to(
            0.0,
            StepConfig {
                max_dt_s: 1.0,
                energy_tolerance_j: 1e-10,
            },
            || Ok(()),
        )
        .unwrap()
        .0
    }
    fn exchange(state: &DepositionSimulation<'_>, id: u64) -> SolidCellExchange {
        SolidCellExchange {
            expected: state.cells()[&CellId(id)].clone(),
            interval_exposed_area_m2: 0.001,
            replacement: SolidSurfaceState {
                mass_kg: 0.009,
                specific_enthalpy_j_kg: -1000.0,
                area_m2: 0.001,
            },
            pure_solid_enthalpy_range_j_kg: [-50_000.0, 30_000.0],
        }
    }

    #[test]
    fn external_endpoint_changes_only_mass_and_original_caloric_enthalpy() {
        let mut state = state();
        let update = exchange(&state, 1);
        let mut expected = update.expected.clone();
        expected.mass_kg = update.replacement.mass_kg;
        expected.specific_h_j_kg = update.replacement.specific_enthalpy_j_kg;
        state.apply_solid_exchanges(&[update], || Ok(())).unwrap();
        assert_eq!(state.cells()[&CellId(1)], expected);
        assert_eq!(state.temperature_k(CellId(1)), Ok(249.5));
        assert_eq!(state.time_s(), 0.0);
    }

    #[test]
    fn internal_energy_endpoint_refuses_the_entire_external_exchange() {
        struct InternalEnergy;
        impl ThermalLaw for InternalEnergy {
            fn storage_basis(&self) -> ThermalStorageBasis {
                ThermalStorageBasis::InternalEnergy
            }
            fn temperature_k(&self, u: f64) -> Result<f64, DepositionError> {
                CALORIC.temperature_k(u)
            }
            fn maximum_temperature_slope(&self) -> f64 {
                CALORIC.maximum_temperature_slope()
            }
        }
        let internal = InternalEnergy;
        let events = state()
            .cells()
            .values()
            .cloned()
            .map(|mut cell| {
                if cell.id == CellId(2) {
                    cell.material_id = MaterialId(2);
                }
                TimedEvent {
                    time_s: 0.0,
                    event: Event::Birth(cell),
                }
            })
            .collect();
        let mut state = DepositionSimulation::new(
            events,
            &[(MaterialId(1), &CALORIC), (MaterialId(2), &internal)],
            DepositionBudget {
                max_cells: 2,
                max_links: 0,
                max_events: 2,
                max_substeps: 4,
            },
        )
        .unwrap()
        .advance_to(
            0.0,
            StepConfig {
                max_dt_s: 1.0,
                energy_tolerance_j: 1e-10,
            },
            || Ok(()),
        )
        .unwrap()
        .0;
        let before = state.cells().clone();
        let updates = [exchange(&state, 1), exchange(&state, 2)];
        assert_eq!(
            state.apply_solid_exchanges(&updates, || Ok(())),
            Err(DepositionError::InvalidInput(
                "external solid exchange requires enthalpy storage and an explicit open-mass model",
            )),
        );
        assert_eq!(state.cells(), &before);
        assert_eq!(state.time_s(), 0.0);
    }

    #[test]
    fn start_exposure_survives_endpoint_covering() {
        let mut state = state();
        let cell = state.cells.get_mut(&CellId(1)).unwrap();
        cell.boundary.reservoir = Some(Reservoir {
            area_m2: 0.002,
            conductance_w_k: 0.0,
            temperature_k: 250.0,
        });
        assert_eq!(state.exposed_area_m2(CellId(1)), Ok(0.0));
        let update = exchange(&state, 1);
        state.apply_solid_exchanges(&[update], || Ok(())).unwrap();
        assert_eq!(state.cells()[&CellId(1)].mass_kg, 0.009);
        assert_eq!(state.exposed_area_m2(CellId(1)), Ok(0.0));
    }

    #[test]
    fn a_bad_later_endpoint_or_stale_geometry_publishes_nothing() {
        let mut state = state();
        let before = state.cells().clone();
        let a = exchange(&state, 1);
        for bad in [0, 1, 2, 3] {
            let mut b = exchange(&state, 2);
            match bad {
                0 => b.replacement.mass_kg = 0.0,
                1 => b.replacement.specific_enthalpy_j_kg = 40_000.0,
                2 => b.expected.geometry.volume_m3 *= 2.0,
                _ => b.replacement.area_m2 *= 2.0,
            }
            assert!(
                state
                    .apply_solid_exchanges(&[a.clone(), b], || Ok(()))
                    .is_err()
            );
            assert_eq!(state.cells(), &before);
        }
    }

    #[test]
    fn late_cancellation_and_duplicate_id_refuse_the_complete_batch() {
        let mut state = state();
        let before = state.cells().clone();
        let updates = [exchange(&state, 1), exchange(&state, 2)];
        let mut calls = 0;
        assert_eq!(
            state.apply_solid_exchanges(&updates, || {
                calls += 1;
                if calls == 3 {
                    Err(DepositionError::Cancelled)
                } else {
                    Ok(())
                }
            }),
            Err(DepositionError::Cancelled)
        );
        assert_eq!(state.cells(), &before);
        assert!(
            state
                .apply_solid_exchanges(&[updates[0].clone(), updates[0].clone()], || Ok(()))
                .is_err()
        );
        assert_eq!(state.cells(), &before);
    }
}
