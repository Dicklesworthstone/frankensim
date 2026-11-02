//! Prepared analytic rheology on the existing sparse free-surface trajectory.
//!
//! The native owner still performs collision, forcing, gas reconstruction,
//! conservative mass exchange and phase conversion. This adapter supplies its
//! existing per-wet-cell BGK relaxation hook, using the same constitutive
//! stress inversion as the plain 2-D/3-D grids. No material card is synthesized.

use std::collections::BTreeMap;

use fs_exec::CancelGate;

use super::transport::MassLocation3;
use super::{Cell3, FreeSurface3, FreeSurfaceError3, TILE, TILE_CELLS, check_cancel};
use crate::d3q19::CollisionModel3;
use crate::rheology::{
    PreparedRheology, Rheology, RheologyField, RheologyFieldError, RheologyStats,
};

type Relaxations = BTreeMap<(u64, usize), f64>;

/// Native flow, field admission, or constitutive inversion refused publication.
#[derive(Debug, Clone, PartialEq)]
pub enum FreeSurfaceRheologyError3 {
    /// The native free-surface owner refused or observed cancellation.
    Flow(FreeSurfaceError3),
    /// A prepared field has the wrong full-domain shape.
    Field(RheologyFieldError),
    /// Uniform-law parameters or the requested collision family are invalid.
    Invalid(&'static str),
    /// A wet cell's local constitutive relation cannot be represented or solved.
    Constitutive {
        /// Stable global coordinates, independent of sparse slot ordering.
        coordinates: [usize; 3],
        /// Failed numerical condition.
        reason: &'static str,
    },
}
impl core::fmt::Display for FreeSurfaceRheologyError3 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Flow(error) => write!(f, "free-surface rheology: {error}"),
            Self::Field(error) => write!(f, "free-surface rheology: {error}"),
            Self::Invalid(reason) => write!(f, "free-surface rheology: {reason}"),
            Self::Constitutive {
                coordinates,
                reason,
            } => write!(f, "free-surface rheology at {coordinates:?}: {reason}"),
        }
    }
}
impl std::error::Error for FreeSurfaceRheologyError3 {}
impl From<FreeSurfaceError3> for FreeSurfaceRheologyError3 {
    fn from(error: FreeSurfaceError3) -> Self {
        Self::Flow(error)
    }
}
impl From<RheologyFieldError> for FreeSurfaceRheologyError3 {
    fn from(error: RheologyFieldError) -> Self {
        Self::Field(error)
    }
}

enum Assignment {
    Uniform(PreparedRheology),
    Field(RheologyField<3>),
}

/// Actual sparse VOF evolution with prepared Newtonian or non-Newtonian laws.
///
/// A uniform law needs no dense material array. A spatial field covers the
/// entire domain, including presently dry/inactive/solid cells, in x-fastest
/// [nx, ny, nz] order. It is a fixed Eulerian assignment: newly wet cells use
/// the law at their global coordinates, never the sparse slot that happens to
/// contain them. This does not advect material identities or infer mixture laws.
///
/// All wet collision cells are resolved, including zero-mass interface ghosts.
/// Gas and unlocated mass carry have no constitutive stress. The local law is
/// solved against the actual retained Guo force frame. A step advances a single
/// private native candidate, resolves its new wet cells, then publishes flow,
/// phase, mass, local relaxation and statistics together. A refusal or observed
/// cancellation leaves every accepted component unchanged.
///
/// Constitutive laws remain generalized Newtonian; finite yield regularization
/// permits creep and supplies no elastic, plastic, thixotropic or phase memory.
/// Use [MaterialScaling::bind_rheology](crate::material::MaterialScaling::bind_rheology)
/// for SI coefficients; retain that source/scale binding for physical results.
pub struct FreeSurfaceRheology3 {
    lattice: FreeSurface3,
    assignment: Assignment,
    tau: Relaxations,
    stats: RheologyStats,
}

impl FreeSurfaceRheology3 {
    /// Bind one law to all present and future wet cells without a dense field.
    pub fn uniform(
        lattice: FreeSurface3,
        law: Rheology,
        gate: &CancelGate,
    ) -> Result<Self, FreeSurfaceRheologyError3> {
        check_cancel(gate)?;
        let prepared = law
            .try_prepare()
            .map_err(FreeSurfaceRheologyError3::Invalid)?;
        Self::new(lattice, Assignment::Uniform(prepared), gate)
    }

    /// Bind a reusable prepared full-domain field. Its dimensions must exactly
    /// match the sparse domain, including dry and permanently solid cells.
    pub fn with_field(
        lattice: FreeSurface3,
        field: RheologyField<3>,
        gate: &CancelGate,
    ) -> Result<Self, FreeSurfaceRheologyError3> {
        Self::new(lattice, Assignment::Field(field), gate)
    }

    fn new(
        lattice: FreeSurface3,
        assignment: Assignment,
        gate: &CancelGate,
    ) -> Result<Self, FreeSurfaceRheologyError3> {
        let (tau, stats) = resolve(&lattice, &assignment, &Relaxations::new(), gate)?;
        Ok(Self {
            lattice,
            assignment,
            tau,
            stats,
        })
    }

    /// The actual native mass/phase/population owner, exposed read-only.
    #[must_use]
    pub const fn lattice(&self) -> &FreeSurface3 {
        &self.lattice
    }

    /// Return the native owner, ending this constitutive binding.
    #[must_use]
    pub fn into_lattice(self) -> FreeSurface3 {
        self.lattice
    }

    /// Prepared Eulerian field, when this owner uses a spatial assignment.
    #[must_use]
    pub fn field(&self) -> Option<&RheologyField<3>> {
        match &self.assignment {
            Assignment::Uniform(_) => None,
            Assignment::Field(field) => Some(field),
        }
    }

    /// Relaxation prepared for the next collision at a current wet cell.
    /// Gas, wall and out-of-domain coordinates return None.
    #[must_use]
    pub fn relaxation_at(&self, coordinates: [usize; 3]) -> Option<f64> {
        let MassLocation3::Cell(key, lane) = FreeSurface3::transport_location(coordinates)? else {
            return None;
        };
        self.tau.get(&(key, lane)).copied()
    }

    /// Constitutive ledger at the accepted endpoint, ready for the next collision.
    /// Counts include wet interface ghosts; this is not a cumulative ledger.
    #[must_use]
    pub const fn stats(&self) -> RheologyStats {
        self.stats
    }

    /// Atomically replace the uniform constitutive assignment at the current
    /// stress state. No populations, force frames, phase or mass are changed.
    pub fn set_uniform_law(
        &mut self,
        law: Rheology,
        gate: &CancelGate,
    ) -> Result<(), FreeSurfaceRheologyError3> {
        check_cancel(gate)?;
        let prepared = law
            .try_prepare()
            .map_err(FreeSurfaceRheologyError3::Invalid)?;
        self.rebind(Assignment::Uniform(prepared), gate)
    }

    /// Atomically replace the prepared field, for example after selecting new
    /// temperature-dependent SI parameters. Full-domain admission and all current
    /// wet-cell solves finish before the old assignment or relaxation is changed.
    pub fn set_field(
        &mut self,
        field: RheologyField<3>,
        gate: &CancelGate,
    ) -> Result<(), FreeSurfaceRheologyError3> {
        self.rebind(Assignment::Field(field), gate)
    }

    fn rebind(
        &mut self,
        assignment: Assignment,
        gate: &CancelGate,
    ) -> Result<(), FreeSurfaceRheologyError3> {
        let (tau, stats) = resolve(&self.lattice, &assignment, &self.tau, gate)?;
        check_cancel(gate)?;
        self.assignment = assignment;
        self.tau = tau;
        self.stats = stats;
        Ok(())
    }

    /// One actual free-surface step. Returns the constitutive ledger used by
    /// this collision; [Self::stats] then describes the accepted next state.
    pub fn step(&mut self) -> Result<RheologyStats, FreeSurfaceRheologyError3> {
        self.step_with_gate(&CancelGate::new_clock_free())
    }

    /// The same step with cancellation before work, between wet tiles, and
    /// immediately before publication. The native hook mutates only this one
    /// private candidate; it does not make another full-grid clone.
    pub fn step_with_gate(
        &mut self,
        gate: &CancelGate,
    ) -> Result<RheologyStats, FreeSurfaceRheologyError3> {
        check_cancel(gate)?;
        let applied = self.stats;
        let mut candidate = self.lattice.clone();
        candidate.step_with_material_relaxation_and_pressure(gate, &self.tau, None)?;
        let (tau, stats) = resolve(&candidate, &self.assignment, &self.tau, gate)?;
        check_cancel(gate)?;
        self.lattice = candidate;
        self.tau = tau;
        self.stats = stats;
        Ok(applied)
    }
}

fn coordinates(key: u64, lane: usize) -> [usize; 3] {
    let (x, y, z) = FreeSurface3::global_of(key, lane);
    [x as usize, y as usize, z as usize]
}

fn resolve(
    lattice: &FreeSurface3,
    assignment: &Assignment,
    warm: &Relaxations,
    gate: &CancelGate,
) -> Result<(Relaxations, RheologyStats), FreeSurfaceRheologyError3> {
    check_cancel(gate)?;
    let CollisionModel3::Bgk { tau: initial_tau } = lattice.grid.collision().0 else {
        return Err(FreeSurfaceRheologyError3::Invalid(
            "analytic free-surface rheology requires the native BGK collision",
        ));
    };
    // Static dispatch preserves the uniform-law path without dense field
    // indexing or repeated parameter validation inside the cell loop.
    match assignment {
        Assignment::Uniform(law) => solve_cells(lattice, |_, _| law, warm, initial_tau, gate),
        Assignment::Field(field) => {
            let (x, y, z) = lattice.grid.tile_dims();
            let shape = [x * TILE, y * TILE, z * TILE];
            field.admit_grid(shape, &[])?;
            solve_cells(
                lattice,
                |key, lane| {
                    let [x, y, z] = coordinates(key, lane);
                    field.prepared((z * shape[1] + y) * shape[0] + x)
                },
                warm,
                initial_tau,
                gate,
            )
        }
    }
}

fn solve_cells<'a>(
    lattice: &FreeSurface3,
    law_at: impl Fn(u64, usize) -> &'a PreparedRheology,
    warm: &Relaxations,
    initial_tau: f64,
    gate: &CancelGate,
) -> Result<(Relaxations, RheologyStats), FreeSurfaceRheologyError3> {
    let mut tau = Relaxations::new();
    let mut stats = RheologyStats {
        floored: 0,
        capped: 0,
        max_shear: 0.0,
        tau_range: (f64::INFINITY, 0.0),
    };
    for (&key, tile) in &lattice.tiles {
        check_cancel(gate)?;
        for lane in 0..TILE_CELLS {
            if tile.cells[lane] == Cell3::Gas {
                continue;
            }
            let numerator = lattice.transport_shear_rate(MassLocation3::Cell(key, lane), 1.0)?;
            let local = law_at(key, lane)
                .try_solve(
                    numerator,
                    warm.get(&(key, lane)).copied().unwrap_or(initial_tau),
                )
                .map_err(|reason| FreeSurfaceRheologyError3::Constitutive {
                    coordinates: coordinates(key, lane),
                    reason,
                })?;
            tau.insert((key, lane), local.tau);
            stats.floored = stats.floored.checked_add(u32::from(local.floored)).ok_or(
                FreeSurfaceRheologyError3::Invalid("rheology floor count overflows"),
            )?;
            stats.capped = stats.capped.checked_add(u32::from(local.capped)).ok_or(
                FreeSurfaceRheologyError3::Invalid("rheology cap count overflows"),
            )?;
            stats.max_shear = stats.max_shear.max(local.shear_rate);
            stats.tau_range.0 = stats.tau_range.0.min(local.tau);
            stats.tau_range.1 = stats.tau_range.1.max(local.tau);
        }
    }
    if stats.tau_range.0 == f64::INFINITY {
        stats.tau_range = (0.0, 0.0);
    }
    check_cancel(gate)?;
    Ok((tau, stats))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::d3q19::freesurface3::{ContactModel3, ConversionRedistribution3};
    use crate::d3q19::sparse::morton3;
    use crate::d3q19::{E3, W3, equilibrium3};

    fn laws() -> [Rheology; 4] {
        [
            Rheology::Newtonian { nu: 0.06 },
            Rheology::PowerLaw { k: 0.02, n: 0.7 },
            Rheology::Carreau {
                nu0: 0.2,
                nu_inf: 0.01,
                lambda: 100.0,
                n: 0.5,
            },
            Rheology::HerschelBulkley {
                k: 0.02,
                n: 1.0,
                yield_stress_over_density: 0.0002,
                regularization_time: 100.0,
                minimum_shear_rate: 1e-8,
            },
        ]
    }

    fn assert_same_state(a: &FreeSurface3, b: &FreeSurface3) {
        assert_eq!(a.grid.state_bits(), b.grid.state_bits());
        assert_eq!(a.mass_snapshot(), b.mass_snapshot());
        assert_eq!(a.conversions(), b.conversions());
        assert_eq!(a.geometry_epoch(), b.geometry_epoch());
        assert_eq!(
            a.worst_step_violation().to_bits(),
            b.worst_step_violation().to_bits()
        );
        assert_eq!(a.local_force, b.local_force);
        for (&key, tile) in &a.tiles {
            assert_eq!(tile.cells, b.tiles[&key].cells);
            assert_eq!(tile.mass, b.tiles[&key].mass);
        }
    }

    #[test]
    fn mixed_wet_cell_roots_and_actual_forced_collision_match_physical_stress() {
        let gate = CancelGate::new_clock_free();
        let force = [0.0002, -0.0001, 0.0003];
        let mut lattice = FreeSurface3::new(
            4,
            4,
            4,
            1.7,
            force,
            0.0,
            ContactModel3::Neutral,
            |_, _, _| true,
        )
        .unwrap();
        let families = laws();
        let mut field = vec![families[0]; 64];
        let mut expected = vec![0.0; 64];
        let mut max_rate = 0.0_f64;
        for z in 0..4 {
            for y in 0..4 {
                for x in 0..4 {
                    let index = (z * 4 + y) * 4 + x;
                    field[index] = families[(x + y + z) % families.len()];
                    let rho = 0.6 + 0.1 * x as f64;
                    let u = [0.03, -0.02, 0.01];
                    let strain = [
                        [0.002 + 0.0002 * y as f64, -0.001, 0.002],
                        [-0.001, -0.003, 0.001],
                        [0.002, 0.001, 0.002],
                    ];
                    let rate = (2.0 * strain.iter().flatten().map(|v| v * v).sum::<f64>()).sqrt();
                    max_rate = max_rate.max(rate);
                    let tau = 0.5 + field[index].viscosity(rate) / crate::CS2;
                    expected[index] = tau;
                    let mut f = equilibrium3(rho, u);
                    for (q, e) in E3.iter().enumerate() {
                        let e = [f64::from(e.0), f64::from(e.1), f64::from(e.2)];
                        let eu: f64 = (0..3).map(|a| e[a] * u[a]).sum();
                        let mut contraction = 0.0;
                        for a in 0..3 {
                            for b in 0..3 {
                                let hermite = e[a] * e[b] - if a == b { crate::CS2 } else { 0.0 };
                                contraction +=
                                    hermite * (-2.0 * rho * crate::CS2 * tau * strain[a][b]);
                            }
                        }
                        let source: f64 = (0..3)
                            .map(|a| (3.0 * (e[a] - u[a]) + 9.0 * eu * e[a]) * force[a])
                            .sum();
                        f[q] +=
                            W3[q] * (contraction / (2.0 * crate::CS2 * crate::CS2) - 0.5 * source);
                    }
                    let (slot, lane) = lattice
                        .grid
                        .resolve_source(x as i64, y as i64, z as i64)
                        .unwrap();
                    lattice.grid.set_populations(slot, lane, f);
                }
            }
        }
        let mass = lattice.ledger_mass();
        let mut run = FreeSurfaceRheology3::with_field(
            lattice,
            RheologyField::new([4; 3], &field).unwrap(),
            &gate,
        )
        .unwrap();
        for (i, expected) in expected.iter().enumerate() {
            assert!(
                (run.relaxation_at([i % 4, (i / 4) % 4, i / 16]).unwrap() - expected).abs() < 4e-12
            );
        }
        assert!((run.stats().max_shear - max_rate).abs() < 1e-13);
        assert_eq!((run.stats().floored, run.stats().capped), (0, 0));
        // Independent BGK+Guo collision, then pull from each interior source.
        // This proves that the resolved local law reaches actual VOF motion.
        let mut streamed = [0.0; 19];
        for (q, e) in E3.iter().enumerate() {
            let p = [2 - e.0, 2 - e.1, 2 - e.2];
            let sample = run
                .lattice
                .sample_at(p[0] as usize, p[1] as usize, p[2] as usize)
                .unwrap();
            let (slot, lane) = run
                .lattice
                .grid
                .resolve_source(i64::from(p[0]), i64::from(p[1]), i64::from(p[2]))
                .unwrap();
            let f = run.lattice.grid.populations(slot, lane);
            let tau = run.relaxation_at(p.map(|v| v as usize)).unwrap();
            let equilibrium = equilibrium3(sample.density, sample.velocity);
            let e = [f64::from(e.0), f64::from(e.1), f64::from(e.2)];
            let eu: f64 = (0..3).map(|a| e[a] * sample.velocity[a]).sum();
            let source: f64 = (0..3)
                .map(|a| (3.0 * (e[a] - sample.velocity[a]) + 9.0 * eu * e[a]) * force[a])
                .sum();
            streamed[q] = f[q] - (f[q] - equilibrium[q]) / tau + (1.0 - 0.5 / tau) * W3[q] * source;
        }
        run.step_with_gate(&gate).unwrap();
        let (slot, lane) = run.lattice.grid.resolve_source(2, 2, 2).unwrap();
        let actual = run.lattice.grid.populations(slot, lane);
        for q in 0..19 {
            assert!((actual[q] - streamed[q]).abs() < 3e-14);
        }
        assert!((run.lattice.ledger_mass() - mass).abs() < 2e-11);
    }

    #[test]
    fn uniform_and_field_paths_match_native_newtonian_and_each_other() {
        let gate = CancelGate::new_clock_free();
        let baseline = FreeSurface3::new(
            8,
            8,
            8,
            1.0,
            [1e-6, 0.0, -2e-6],
            0.0,
            ContactModel3::Neutral,
            |_, _, z| z < 4,
        )
        .unwrap()
        .with_conversion_redistribution(ConversionRedistribution3::ConservativeLocal);
        let law = Rheology::Newtonian {
            nu: crate::CS2 * 0.5,
        };
        let mut raw = baseline.clone();
        let mut prepared = FreeSurfaceRheology3::uniform(baseline.clone(), law, &gate).unwrap();
        for _ in 0..8 {
            raw.step_with_gate(&gate).unwrap();
            prepared.step_with_gate(&gate).unwrap();
            assert_same_state(&raw, prepared.lattice());
        }
        for law in laws() {
            let mut uniform = FreeSurfaceRheology3::uniform(baseline.clone(), law, &gate).unwrap();
            let mut spatial = FreeSurfaceRheology3::with_field(
                baseline.clone(),
                RheologyField::new([8; 3], &vec![law; 512]).unwrap(),
                &gate,
            )
            .unwrap();
            let mass = uniform.lattice.ledger_mass();
            for _ in 0..8 {
                let a = uniform.step_with_gate(&gate).unwrap();
                let b = spatial.step_with_gate(&gate).unwrap();
                assert_same_state(uniform.lattice(), spatial.lattice());
                assert_eq!(uniform.tau, spatial.tau);
                assert_eq!(
                    (a.floored, a.capped, a.max_shear, a.tau_range),
                    (b.floored, b.capped, b.max_shear, b.tau_range)
                );
            }
            assert!((uniform.lattice.ledger_mass() - mass).abs() < 2e-10);
        }
    }

    #[test]
    fn new_wet_cells_use_global_field_after_sparse_slot_reshuffle() {
        let gate = CancelGate::new_clock_free();
        let mut lattice = FreeSurface3::new(
            12,
            4,
            4,
            1.0,
            [0.0; 3],
            0.0,
            ContactModel3::Neutral,
            |x, _, _| x >= 8,
        )
        .unwrap();
        assert_eq!(lattice.retire_gas_tiles().unwrap(), 1);
        let original_key = morton3(2, 0, 0);
        let original_slot = lattice.grid.slot_of(original_key).unwrap();
        let (slot, lane) = lattice.grid.resolve_source(8, 1, 1).unwrap();
        let rho = lattice.grid.cell_macros(slot, lane).0;
        assert_eq!(lattice.tiles[&original_key].cells[lane], Cell3::Interface);
        // A bounded transport overshoot triggers the native conversion cascade.
        lattice.tiles.get_mut(&original_key).unwrap().mass[lane] = 1.05 * rho;
        let mass = lattice.ledger_mass();
        let field: Vec<_> = (0..192)
            .map(|i| Rheology::Newtonian {
                nu: 0.02 + 0.01 * (i % 12) as f64,
            })
            .collect();
        let mut run = FreeSurfaceRheology3::with_field(
            lattice,
            RheologyField::new([12, 4, 4], &field).unwrap(),
            &gate,
        )
        .unwrap();
        assert!(run.relaxation_at([7, 1, 1]).is_none());
        run.step_with_gate(&gate).unwrap();
        assert!(run.lattice.conversions().tiles_activated > 0);
        assert!(run.lattice.conversions().gas_to_interface > 0);
        assert_ne!(
            run.lattice.grid.slot_of(original_key).unwrap(),
            original_slot
        );
        assert_eq!(
            run.lattice.sample_at(7, 1, 1).unwrap().cell,
            Cell3::Interface
        );
        assert_eq!(
            run.relaxation_at([7, 1, 1]).unwrap(),
            0.5 + field[7].viscosity(0.0) / crate::CS2
        );
        for (&key, tile) in &run.lattice.tiles {
            for lane in 0..TILE_CELLS {
                let p = coordinates(key, lane);
                if tile.cells[lane] == Cell3::Gas {
                    assert!(run.relaxation_at(p).is_none());
                } else {
                    assert_eq!(
                        run.relaxation_at(p).unwrap(),
                        0.5 + field[p[0]].viscosity(0.0) / crate::CS2
                    );
                }
            }
        }
        assert!((run.lattice.ledger_mass() - mass).abs() < 2e-11);
    }

    #[test]
    fn rebind_cancellation_and_late_collision_refusal_preserve_accepted_state() {
        let gate = CancelGate::new_clock_free();
        let lattice = FreeSurface3::new(
            4,
            4,
            4,
            1.0,
            [0.0; 3],
            0.0,
            ContactModel3::Neutral,
            |_, _, _| true,
        )
        .unwrap();
        let mut run = FreeSurfaceRheology3::uniform(lattice, laws()[2], &gate).unwrap();
        let before = run.lattice.clone();
        let tau = run.tau.clone();
        assert!(
            run.set_uniform_law(Rheology::Newtonian { nu: f64::NAN }, &gate)
                .is_err()
        );
        let wrong = RheologyField::new([8, 2, 4], &vec![laws()[0]; 64]).unwrap();
        assert!(matches!(
            run.set_field(wrong, &gate),
            Err(FreeSurfaceRheologyError3::Field(_))
        ));
        assert_same_state(&before, &run.lattice);
        assert_eq!(run.tau, tau);
        let cancelled = CancelGate::new_clock_free();
        cancelled.request();
        assert!(matches!(
            run.step_with_gate(&cancelled),
            Err(FreeSurfaceRheologyError3::Flow(
                FreeSurfaceError3::Cancelled
            ))
        ));
        assert_same_state(&before, &run.lattice);
        assert_eq!(run.tau, tau);
        // Fault injection at the last wet lane forces the native collision
        // hook to refuse after earlier lanes have produced private post state.
        let last = *run.tau.last_key_value().unwrap().0;
        run.tau.remove(&last);
        let incomplete = run.tau.clone();
        assert!(matches!(
            run.step_with_gate(&gate),
            Err(FreeSurfaceRheologyError3::Flow(
                FreeSurfaceError3::Transport(_)
            ))
        ));
        assert_same_state(&before, &run.lattice);
        assert_eq!(run.tau, incomplete);
    }
}
