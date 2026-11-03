//! Thermal-owner admission for externally solved solid/gas exchange.
//!
//! The external owner supplies an admitted enthalpy interval and the complete
//! expected cell. This owner retains identities, caloric laws, reference
//! geometry, contacts, clock and schedule; only accepted mass and h can change.
//! Phase identity is caller-admitted: a temperature or an interval alone cannot
//! prove that an arbitrary material chart describes a pure solid. Mixed-phase
//! exchange additionally binds the original caloric identity; its physical
//! owner must admit the constituent inventories before submitting an endpoint.

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

/// A phase-partitioned endpoint on an explicitly identified original caloric law.
///
/// The external phase owner admits constituent masses, latent storage and exposed
/// solid area. This boundary only commits its complete bulk mass/enthalpy update;
/// it does not infer a phase fraction from temperature or add latent heat.
#[derive(Clone, Debug, PartialEq)]
pub struct CondensedCellExchange {
    pub expected: DepositedCell,
    pub interval_exposed_area_m2: f64,
    pub replacement: SolidSurfaceState,
    /// Inclusive original-source bulk enthalpy interval [J/kg].
    pub enthalpy_range_j_kg: [f64; 2],
    /// Exact identity of the original registered thermal law, including datum.
    pub caloric_identity: [u8; 32],
}

#[derive(Clone, Copy)]
struct Exchange<'a> {
    expected: &'a DepositedCell,
    area: f64,
    replacement: SolidSurfaceState,
    range: [f64; 2],
    caloric_identity: Option<[u8; 32]>,
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
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), DepositionError> {
        self.apply_phase_exchanges(
            exchanges.iter().map(|exchange| Exchange {
                expected: &exchange.expected,
                area: exchange.interval_exposed_area_m2,
                replacement: exchange.replacement,
                range: exchange.pure_solid_enthalpy_range_j_kg,
                caloric_identity: None,
            }),
            checkpoint,
        )
    }

    /// Atomically accept an original-source mixed-phase bulk endpoint.
    ///
    /// Shares the pure-solid path's stale-cell, storage, geometry, budget and
    /// cancellation checks. The caller must first resolve and admit every
    /// constituent on the SAME caloric identity. No melting/nucleation model
    /// or new caloric interpretation is introduced by this boundary.
    pub fn apply_condensed_exchanges(
        &mut self,
        exchanges: &[CondensedCellExchange],
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), DepositionError> {
        self.apply_phase_exchanges(
            exchanges.iter().map(|exchange| Exchange {
                expected: &exchange.expected,
                area: exchange.interval_exposed_area_m2,
                replacement: exchange.replacement,
                range: exchange.enthalpy_range_j_kg,
                caloric_identity: Some(exchange.caloric_identity),
            }),
            checkpoint,
        )
    }

    fn apply_phase_exchanges<'a>(
        &mut self,
        exchanges: impl ExactSizeIterator<Item = Exchange<'a>> + Clone,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), DepositionError> {
        checkpoint()?;
        if exchanges.len() > self.budget.max_cells {
            return Err(DepositionError::Budget("external solid exchange cells"));
        }
        let mut ids = BTreeSet::new();
        for (index, exchange) in exchanges.clone().enumerate() {
            if index % 256 == 0 {
                checkpoint()?;
            }
            let id = exchange.expected.id;
            let cell = self.cells.get(&id).ok_or(DepositionError::InvalidInput(
                "external solid exchange requires an active cell",
            ))?;
            if !ids.insert(id) || cell != exchange.expected {
                return Err(DepositionError::InvalidInput(
                    "duplicate or stale external solid exchange",
                ));
            }
            if self.material_storage_basis(cell.material_id)? != ThermalStorageBasis::Enthalpy {
                return Err(DepositionError::InvalidInput(
                    "external solid exchange requires enthalpy storage and an explicit open-mass model",
                ));
            }
            if let Some(identity) = exchange.caloric_identity {
                if self.material_caloric_identity(cell.material_id)? != Some(identity) {
                    return Err(DepositionError::InvalidInput(
                        "external condensed exchange caloric identity mismatch",
                    ));
                }
            }
            let [lo, hi] = exchange.range;
            let area = exchange.area;
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
        fn caloric_identity(&self) -> Option<[u8; 32]> {
            Some([17; 32])
        }
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
    fn condensed_batch_requires_original_source_identity_and_commits_atomically() {
        let mut state = state();
        let before = state.cells().clone();
        let update = |id| {
            let row = exchange(&state, id);
            CondensedCellExchange {
                expected: row.expected,
                interval_exposed_area_m2: row.interval_exposed_area_m2,
                replacement: row.replacement,
                enthalpy_range_j_kg: [-50_000.0, 30_000.0],
                caloric_identity: [17; 32],
            }
        };
        let mut rows = [update(1), update(2)];
        rows[1].caloric_identity = [18; 32];
        assert_eq!(
            state.apply_condensed_exchanges(&rows, || Ok(())),
            Err(DepositionError::InvalidInput(
                "external condensed exchange caloric identity mismatch"
            ))
        );
        assert_eq!(state.cells(), &before);
        rows[1].caloric_identity = [17; 32];
        let mut calls = 0;
        assert_eq!(
            state.apply_condensed_exchanges(&rows, || {
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
        state.apply_condensed_exchanges(&rows, || Ok(())).unwrap();
        for row in rows {
            let cell = &state.cells()[&row.expected.id];
            assert_eq!(cell.mass_kg, row.replacement.mass_kg);
            assert_eq!(cell.specific_h_j_kg, row.replacement.specific_enthalpy_j_kg);
            assert_eq!(cell.geometry, row.expected.geometry);
        }
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
