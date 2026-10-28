//! Held local forces in the original Guo population frame.

use super::*;
use crate::d3q19::{W3, guo_projection3};

pub(super) fn source(q: usize, velocity: [f64; 3], force: [f64; 3]) -> f64 {
    let e = [f64::from(E3[q].0), f64::from(E3[q].1), f64::from(E3[q].2)];
    let eu = e.iter().zip(velocity).map(|(e, u)| e * u).sum();
    W3[q] * guo_projection3(e, velocity, eu, force)
}

// Preserve the actual canonical population reduction, not only its analytic
// zeroth moment. Even a one-ulp density change can turn a full interface into
// a fictitious partial cell. The rest correction changes no momentum/stress;
// try its adjacent representations before closing the final ordered add.
fn preserve_density(populations: &mut [f64; Q3], rho: f64) -> Result<(), FreeSurfaceError3> {
    populations[0] += rho - populations.iter().sum::<f64>();
    let corrected_rest = populations[0];
    let actual = populations.iter().sum::<f64>();
    if actual == rho && rho.is_finite() {
        return Ok(());
    }
    let below = actual < rho;
    for _ in 0..8 {
        populations[0] = if below {
            populations[0].next_up()
        } else {
            populations[0].next_down()
        };
        let actual = populations.iter().sum::<f64>();
        if actual == rho && rho.is_finite() {
            return Ok(());
        }
        if !actual.is_finite() || (actual < rho) != below {
            break;
        }
    }
    // Rounded intermediate additions can skip the target for every nearby
    // rest value. Restore that correction and close the LAST add instead.
    // For physical positive populations prefix is near rho, so Sterbenz's
    // subtraction is exact; verify the resulting reduction even outside that
    // case. This moving-population correction affects momentum/stress
    // only at the explicitly bounded population-reduction roundoff scale.
    populations[0] = corrected_rest;
    let prefix = populations[..Q3 - 1].iter().sum::<f64>();
    let last = rho - prefix;
    let magnitude = populations.iter().map(|value| value.abs()).sum::<f64>();
    if !last.is_finite()
        || !magnitude.is_finite()
        || (last - populations[Q3 - 1]).abs() > 64.0 * f64::EPSILON * magnitude
    {
        return Err(FreeSurfaceError3::Transport(
            "force-frame density closure exceeds population roundoff",
        ));
    }
    populations[Q3 - 1] = last;
    if populations.iter().sum::<f64>() != rho {
        return Err(FreeSurfaceError3::Transport(
            "force-frame density is not exactly representable",
        ));
    }
    Ok(())
}

pub(super) fn reframe(
    mut populations: [f64; Q3],
    physical_velocity: [f64; 3],
    previous: [f64; 3],
    next: [f64; 3],
) -> Result<[f64; Q3], FreeSurfaceError3> {
    if previous == next {
        return Ok(populations);
    }
    let rho = populations.iter().sum::<f64>();
    let delta = core::array::from_fn(|d| previous[d] - next[d]);
    for (q, value) in populations.iter_mut().enumerate() {
        *value += 0.5 * source(q, physical_velocity, delta);
    }
    // The analytic source has zero mass. Close the same ordered reduction
    // used by native density and geometry, never changing VOF mass or heat.
    preserve_density(&mut populations, rho)?;
    if populations.iter().any(|v| !v.is_finite()) {
        return Err(FreeSurfaceError3::Transport(
            "unrepresentable force-frame change",
        ));
    }
    Ok(populations)
}

impl FreeSurface3 {
    pub(crate) fn reference_force(&self) -> [f64; 3] {
        self.grid.collision().1
    }

    pub(super) fn force_at_lane(&self, key: u64, lane: usize) -> [f64; 3] {
        self.local_force
            .as_ref()
            .map_or_else(|| self.grid.collision().1, |forces| forces[&(key, lane)])
    }

    /// Actual force density in lattice units at an accepted wet cell. Its
    /// half-source correction is already included by `sample_at`; adding it
    /// again would count the same force twice. Gas/wall/outside returns None.
    #[must_use]
    pub fn force_density_at(&self, x: usize, y: usize, z: usize) -> Option<[f64; 3]> {
        let (slot, lane) = self.grid.resolve_source(
            i64::try_from(x).ok()?,
            i64::try_from(y).ok()?,
            i64::try_from(z).ok()?,
        )?;
        let key = self.grid.active_keys()[slot];
        (self.tiles.get(&key)?.cells[lane] != Cell3::Gas).then(|| self.force_at_lane(key, lane))
    }

    /// Change only the source frame, preserving canonical density exactly and
    /// physical velocity/force-corrected second moment to population roundoff.
    /// The caller owns the enclosing candidate.
    /// No elapsed time, impulse, heat, mass route or mechanical work is added.
    pub(crate) fn rebind_force_field(
        &mut self,
        forces: BTreeMap<(u64, usize), [f64; 3]>,
        gate: &CancelGate,
    ) -> Result<(), FreeSurfaceError3> {
        check_cancel(gate)?;
        let wet = self
            .tiles
            .values()
            .map(|tile| {
                tile.cells
                    .iter()
                    .filter(|&&cell| cell != Cell3::Gas)
                    .count()
            })
            .sum::<usize>();
        if forces.len() != wet {
            return Err(FreeSurfaceError3::Transport(
                "local force must cover each wet cell once",
            ));
        }
        let uniform = self.grid.collision().1;
        let mut replacements = Vec::new();
        for (&key, tile) in &self.tiles {
            check_cancel(gate)?;
            let slot = self.grid.slot_of(key).expect("active");
            for lane in 0..TILE_CELLS {
                if tile.cells[lane] == Cell3::Gas {
                    continue;
                }
                let next = *forces
                    .get(&(key, lane))
                    .ok_or(FreeSurfaceError3::Transport(
                        "missing local force at a wet cell",
                    ))?;
                if next.iter().any(|f| !f.is_finite()) {
                    return Err(FreeSurfaceError3::Transport("nonfinite local force"));
                }
                let previous = self.force_at_lane(key, lane);
                // Preserve the original population bits when the force does
                // not change, including the constant-density gravity limit.
                if previous == next {
                    continue;
                }
                let populations = self.grid.populations(slot, lane);
                let (rho, raw) = self.grid.cell_macros(slot, lane);
                if !rho.is_finite() || rho <= 0.0 {
                    return Err(FreeSurfaceError3::Transport(
                        "invalid density in force-frame change",
                    ));
                }
                let velocity = core::array::from_fn(|d| raw[d] + 0.5 * previous[d] / rho);
                let populations = reframe(populations, velocity, previous, next)?;
                replacements.push((slot, lane, populations));
            }
        }
        check_cancel(gate)?;
        for (slot, lane, populations) in replacements {
            self.grid.set_populations(slot, lane, populations);
        }
        self.local_force = if forces.values().all(|&force| force == uniform) {
            None
        } else {
            Some(forces)
        };
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Guo's independently known first and second moments: F and uF + Fu.
    // The physical moments therefore add half that source to stored f.
    fn physical_moments(
        populations: &[f64; Q3],
        velocity: [f64; 3],
        force: [f64; 3],
    ) -> ([f64; 3], [[f64; 3]; 3]) {
        let mut momentum = force.map(|v| 0.5 * v);
        let mut stress = core::array::from_fn(|i| {
            core::array::from_fn(|j| 0.5 * (velocity[i] * force[j] + velocity[j] * force[i]))
        });
        for (q, &population) in populations.iter().enumerate() {
            let e = [E3[q].0, E3[q].1, E3[q].2].map(f64::from);
            for i in 0..3 {
                momentum[i] += population * e[i];
                for j in 0..3 {
                    stress[i][j] += population * e[i] * e[j];
                }
            }
        }
        (momentum, stress)
    }

    #[test]
    fn g1_force_reframe_preserves_exact_density_and_physical_moments() {
        let mut seed = 0x6a09_e667_f3bc_c909_u64;
        let mut uniform = || {
            seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
            (seed >> 11) as f64 / (1_u64 << 53) as f64
        };
        let mut tail_closures = 0;
        for _ in 0..512 {
            let nominal_density = 0.9 + 0.2 * uniform();
            let velocity = core::array::from_fn(|_| 0.1 * (uniform() - 0.5));
            let previous = core::array::from_fn(|_| 0.002 * (uniform() - 0.5));
            let next = core::array::from_fn(|_| 0.002 * (uniform() - 0.5));
            let initial = equilibrium3(nominal_density, velocity);
            let rho = initial.iter().sum::<f64>();
            let first = reframe(initial, velocity, [0.0; 3], previous).unwrap();
            let after = reframe(first, velocity, previous, next).unwrap();
            assert!(first.iter().chain(&after).all(|&value| value > 0.0));
            assert_eq!(first.iter().sum::<f64>().to_bits(), rho.to_bits());
            assert_eq!(after.iter().sum::<f64>().to_bits(), rho.to_bits());
            let raw_last = first[Q3 - 1]
                + 0.5
                    * source(
                        Q3 - 1,
                        velocity,
                        core::array::from_fn(|d| previous[d] - next[d]),
                    );
            tail_closures += usize::from(after[Q3 - 1].to_bits() != raw_last.to_bits());
            let (initial_momentum, initial_stress) = physical_moments(&initial, velocity, [0.0; 3]);
            let (momentum, stress) = physical_moments(&after, velocity, next);
            for i in 0..3 {
                assert!((momentum[i] - initial_momentum[i]).abs() <= 8.0 * f64::EPSILON * rho);
                for j in 0..3 {
                    assert!(
                        (stress[i][j] - initial_stress[i][j]).abs() <= 8.0 * f64::EPSILON * rho
                    );
                }
            }
            // A true full interface remains full, and every actual partial
            // volume keeps the same bit pattern when its stored mass is held.
            for fraction in [0.0, 1e-8, 0.2, 0.999999, 1.0] {
                let mass = fraction * rho;
                assert_eq!(
                    (mass / rho).to_bits(),
                    (mass / after.iter().sum::<f64>()).to_bits()
                );
            }
            let unchanged = reframe(after, velocity, next, next).unwrap();
            assert_eq!(unchanged.map(f64::to_bits), after.map(f64::to_bits));
        }
        assert!(
            tail_closures > 0,
            "exercise unattainable rest-only reductions"
        );
        eprintln!(
            "force-frame density: 512 positive low-Mach states, {tail_closures} final-add closures"
        );
    }
}
