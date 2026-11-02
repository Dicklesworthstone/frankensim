//! Capacity-bounded, momentum-consistent local conversion of signed VOF mass.

use super::*;

/// Balance of the most recent explicitly conservative conversion, in lattice
/// units. Collision, streaming, pressure traction and wall work are excluded.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ConversionMomentumReport3 {
    /// Sum of absolute actual mass transfers on adjacent interface routes.
    pub redistributed_mass: f64,
    /// Sum of absolute source/recipient mass-representation defects.
    pub absolute_mass_defect: f64,
    /// Sum of absolute local defects in each physical momentum component.
    pub absolute_momentum_defect: [f64; 3],
    /// Signed resolved kinetic change caused by numerical conversion alone.
    /// This is remap error, NEVER an authored physical heat source.
    pub kinetic_energy_change: f64,
    /// Number of nonzero represented adjacent-interface routes.
    pub routes: usize,
}

struct Source3 {
    velocity: [f64; 3],
    density: f64,
}

struct Recipient3 {
    mass: f64,
    velocity: [f64; 3],
    relative_momentum: [f64; 3],
    carried_momentum: [f64; 3],
    transferred_mass: f64,
}

pub(super) struct ConversionMomentum3 {
    sources: BTreeMap<(u64, usize), Source3>,
    recipients: BTreeMap<(u64, usize), Recipient3>,
    report: ConversionMomentumReport3,
}

fn refused(reason: &'static str) -> FreeSurfaceError3 {
    FreeSurfaceError3::Transport(reason)
}

fn state(
    surface: &FreeSurface3,
    key: u64,
    lane: usize,
) -> Result<(f64, [f64; 3]), FreeSurfaceError3> {
    let slot = surface.grid.slot_of(key).expect("active conversion cell");
    let (density, raw) = surface.grid.cell_macros(slot, lane);
    if !density.is_finite() || density <= 0.0 {
        return Err(refused(
            "conservative conversion needs positive finite density",
        ));
    }
    let force = surface.force_at_lane(key, lane);
    let velocity = core::array::from_fn(|d| raw[d] + 0.5 * force[d] / density);
    admit_velocity(velocity)?;
    Ok((density, velocity))
}

fn admit_velocity(velocity: [f64; 3]) -> Result<(), FreeSurfaceError3> {
    let speed = velocity[0].hypot(velocity[1]).hypot(velocity[2]);
    if !speed.is_finite() || speed / CS2_3.sqrt() >= crate::MACH_LIMIT {
        return Err(refused(
            "conservative conversion exceeds the lattice Mach domain; refine",
        ));
    }
    Ok(())
}

impl ConversionMomentum3 {
    pub(super) fn new(surface: &FreeSurface3) -> Result<Self, FreeSurfaceError3> {
        if surface.carry != 0.0 {
            return Err(refused(
                "conservative conversion cannot locate the existing mass carry",
            ));
        }
        Ok(Self {
            sources: BTreeMap::new(),
            recipients: BTreeMap::new(),
            report: ConversionMomentumReport3::default(),
        })
    }

    /// Capture BEFORE class/force removal. A negative old cell contributes
    /// signed m*u; its debt therefore carries THIS velocity when corrected,
    /// even though the scalar transport route points toward the removed cell.
    pub(super) fn capture_source(
        &mut self,
        surface: &FreeSurface3,
        key: u64,
        lane: usize,
        amount: f64,
    ) -> Result<(), FreeSurfaceError3> {
        let (density, velocity) = state(surface, key, lane)?;
        if !amount.is_finite() {
            return Err(refused("nonfinite conservative conversion mass"));
        }
        self.sources
            .insert((key, lane), Source3 { velocity, density });
        self.report.kinetic_energy_change -=
            0.5 * amount * velocity.iter().map(|u| u * u).sum::<f64>();
        Ok(())
    }

    pub(super) fn redistribute(
        &mut self,
        surface: &mut FreeSurface3,
        key: u64,
        lane: usize,
        amount: f64,
        mut transport: Option<&mut MassTransport3>,
    ) -> Result<(), FreeSurfaceError3> {
        if amount == 0.0 {
            return Ok(());
        }
        let source = &self.sources[&(key, lane)];
        let velocity = source.velocity;
        let mut neighbors = [(0, 0, 0.0, 0.0); Q3 - 1];
        let mut count = 0;
        let mut capacity = 0.0;
        let mut mass_scale = source.density + amount.abs();
        for q in 1..Q3 {
            if let Neighbor::Active(nb_key, nb_lane, Cell3::Interface) =
                surface.neighbor_state(key, lane, q)
            {
                let slot = surface.grid.slot_of(nb_key).expect("active recipient");
                let density = surface.grid.populations(slot, nb_lane).iter().sum::<f64>();
                let mass = surface.tiles[&nb_key].mass[nb_lane];
                if !density.is_finite()
                    || density <= 0.0
                    || !mass.is_finite()
                    || mass < 0.0
                    || mass > density
                {
                    return Err(refused(
                        "conservative conversion needs bounded adjacent inventories; refine",
                    ));
                }
                let available = if amount > 0.0 { density - mass } else { mass };
                if available > 0.0 {
                    neighbors[count] = (nb_key, nb_lane, density, available);
                    count += 1;
                    capacity += available;
                    mass_scale += density;
                }
            }
        }
        let mass_roundoff = 128.0 * f64::EPSILON * mass_scale;
        if count == 0
            || !capacity.is_finite()
            || !mass_roundoff.is_finite()
            || amount.abs() > capacity + mass_roundoff
        {
            return Err(refused(
                "conservative conversion lacks adjacent mass/capacity; refine or remap",
            ));
        }

        // Distribute in the existing direction order, proportionally to actual
        // capacity. A dry ghost cannot pay debt; a full cell cannot overfill.
        // The last recipient closes the remaining representable transfer.
        let mut remaining = amount;
        let mut remaining_capacity = capacity;
        for (i, &(nb_key, nb_lane, density, available)) in neighbors[..count].iter().enumerate() {
            let magnitude = if amount > 0.0 {
                remaining.max(0.0)
            } else {
                (-remaining).max(0.0)
            };
            let wanted = if i + 1 == count {
                magnitude
            } else {
                magnitude * (available / remaining_capacity.max(available))
            };
            let selected = wanted.min(available);
            let old_mass = surface.tiles[&nb_key].mass[nb_lane];
            let next_mass = if selected == available {
                if amount > 0.0 { density } else { 0.0 }
            } else if amount > 0.0 {
                old_mass + selected
            } else {
                old_mass - selected
            };
            let actual = next_mass - old_mass;
            if !next_mass.is_finite()
                || next_mass < 0.0
                || next_mass > density
                || (amount > 0.0 && actual < 0.0)
                || (amount < 0.0 && actual > 0.0)
            {
                return Err(refused(
                    "conservative conversion mass cannot be represented",
                ));
            }
            if actual != 0.0 {
                let recipient = match self.recipients.entry((nb_key, nb_lane)) {
                    std::collections::btree_map::Entry::Vacant(entry) => {
                        let (_, old_velocity) = state(surface, nb_key, nb_lane)?;
                        entry.insert(Recipient3 {
                            mass: old_mass,
                            velocity: old_velocity,
                            relative_momentum: [0.0; 3],
                            carried_momentum: [0.0; 3],
                            transferred_mass: 0.0,
                        })
                    }
                    std::collections::btree_map::Entry::Occupied(entry) => entry.into_mut(),
                };
                if let Some(trace) = transport.as_deref_mut() {
                    trace.exchange(
                        MassLocation3::Cell(key, lane),
                        MassLocation3::Cell(nb_key, nb_lane),
                        actual,
                    )?;
                }
                for d in 0..3 {
                    recipient.relative_momentum[d] +=
                        actual * (velocity[d] - recipient.velocity[d]);
                    recipient.carried_momentum[d] += actual * velocity[d];
                }
                recipient.transferred_mass += actual.abs();
                surface
                    .tiles
                    .get_mut(&nb_key)
                    .expect("active recipient")
                    .mass[nb_lane] = next_mass;
                self.report.redistributed_mass += actual.abs();
                self.report.routes += 1;
            }
            remaining -= actual;
            remaining_capacity -= available;
        }
        if !remaining.is_finite() || remaining.abs() > mass_roundoff {
            return Err(refused(
                "conservative conversion mass defect exceeds local roundoff; refine",
            ));
        }
        self.report.absolute_mass_defect += remaining.abs();
        for d in 0..3 {
            self.report.absolute_momentum_defect[d] += (remaining * velocity[d]).abs();
        }
        Ok(())
    }

    pub(super) fn finish(
        mut self,
        surface: &mut FreeSurface3,
    ) -> Result<ConversionMomentumReport3, FreeSurfaceError3> {
        // One population update per recipient, after ALL sources have routed.
        // Centering on its old velocity avoids subtracting large momenta and
        // preserves the co-moving state without a needless equilibrium reset.
        for ((key, lane), recipient) in self.recipients {
            let mass = surface.tiles[&key].mass[lane];
            let slot = surface.grid.slot_of(key).expect("active recipient");
            if mass > 0.0 {
                let next = core::array::from_fn(|d| {
                    recipient.velocity[d] + recipient.relative_momentum[d] / mass
                });
                admit_velocity(next)?;
                if next != recipient.velocity {
                    let populations = forcing::retarget_velocity(
                        surface.grid.populations(slot, lane),
                        recipient.velocity,
                        next,
                        surface.force_at_lane(key, lane),
                    )?;
                    surface.grid.set_populations(slot, lane, populations);
                }
            }
            let (_, actual) = state(surface, key, lane)?;
            let dm = mass - recipient.mass;
            let du: [f64; 3] = core::array::from_fn(|d| actual[d] - recipient.velocity[d]);
            let momentum_roundoff =
                256.0 * f64::EPSILON * (mass + recipient.mass + recipient.transferred_mass);
            for d in 0..3 {
                let change = dm * recipient.velocity[d] + mass * du[d];
                let defect = (change - recipient.carried_momentum[d]).abs();
                if !defect.is_finite() || defect > momentum_roundoff {
                    return Err(refused(
                        "conservative conversion momentum cannot fit positive bounded cells; refine",
                    ));
                }
                self.report.absolute_momentum_defect[d] += defect;
            }
            self.report.kinetic_energy_change +=
                0.5 * dm * recipient.velocity.iter().map(|u| u * u).sum::<f64>()
                    + mass
                        * (0..3)
                            .map(|d| (recipient.velocity[d] + 0.5 * du[d]) * du[d])
                            .sum::<f64>();
        }
        if [
            self.report.redistributed_mass,
            self.report.absolute_mass_defect,
            self.report.kinetic_energy_change,
        ]
        .iter()
        .chain(&self.report.absolute_momentum_defect)
        .any(|v| !v.is_finite())
        {
            return Err(refused("conservative conversion balance overflowed"));
        }
        Ok(self.report)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cell(coordinates: [usize; 3], density: f64, velocity: [f64; 3]) -> FreeSurfaceInitialCell3 {
        FreeSurfaceInitialCell3 {
            coordinates,
            fill_fraction: 0.4,
            density,
            velocity,
        }
    }

    fn surface(cells: &[FreeSurfaceInitialCell3], force: [f64; 3]) -> FreeSurface3 {
        let gate = CancelGate::new_clock_free();
        fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
            let cx = fs_exec::Cx::new(
                &gate,
                arena,
                fs_exec::StreamKey {
                    seed: 79,
                    kernel_id: 31,
                    tile: 0,
                    iteration: 0,
                },
                fs_exec::Budget::INFINITE,
                fs_exec::ExecMode::Deterministic,
            );
            FreeSurface3::from_cells(
                [16, 8, 8],
                0.8,
                force,
                0.0,
                ContactModel3::Neutral,
                &[],
                cells,
                &cx,
            )
            .unwrap()
            .with_conversion_redistribution(ConversionRedistribution3::ConservativeLocal)
        })
    }

    fn set_mass(surface: &mut FreeSurface3, coordinates: [usize; 3], mass: f64) {
        let [x, y, z] = coordinates.map(|p| p as i64);
        let (slot, lane) = surface.grid.resolve_source(x, y, z).unwrap();
        let key = surface.grid.active_keys()[slot];
        let tile = surface.tiles.get_mut(&key).unwrap();
        assert_eq!(tile.cells[lane], Cell3::Interface);
        tile.mass[lane] = mass;
    }

    // Independently reduce actual signed tracked mass and public physical
    // velocity, including the post-stream overshoot/debt before conversion.
    fn inventory(surface: &FreeSurface3, region: impl Fn(i64) -> bool) -> (f64, [f64; 3], f64) {
        let mut mass = 0.0;
        let mut momentum = [0.0; 3];
        let mut kinetic = 0.0;
        for (&key, tile) in &surface.tiles {
            for lane in 0..TILE_CELLS {
                let (x, y, z) = FreeSurface3::global_of(key, lane);
                if tile.cells[lane] == Cell3::Gas || !region(x) {
                    continue;
                }
                let sample = surface
                    .sample_at(x as usize, y as usize, z as usize)
                    .unwrap();
                let m = if tile.cells[lane] == Cell3::Fluid {
                    sample.density
                } else {
                    tile.mass[lane]
                };
                mass += m;
                for d in 0..3 {
                    momentum[d] += m * sample.velocity[d];
                }
                kinetic += 0.5 * m * sample.velocity.iter().map(|u| u * u).sum::<f64>();
            }
        }
        (mass, momentum, kinetic)
    }

    fn density_bits(surface: &FreeSurface3) -> BTreeMap<(u64, usize), u64> {
        let mut values = BTreeMap::new();
        for (&key, tile) in &surface.tiles {
            let slot = surface.grid.slot_of(key).unwrap();
            for lane in 0..TILE_CELLS {
                if tile.cells[lane] != Cell3::Gas {
                    values.insert(
                        (key, lane),
                        surface
                            .grid
                            .populations(slot, lane)
                            .iter()
                            .sum::<f64>()
                            .to_bits(),
                    );
                }
            }
        }
        values
    }

    fn check_balance(surface: &FreeSurface3, before: (f64, [f64; 3], f64)) {
        let after = inventory(surface, |_| true);
        let scale = before.0.abs().max(after.0.abs()).max(1.0);
        let report = surface.conversion_momentum_report().unwrap();
        assert!((after.0 - before.0).abs() < 2e-13 * scale);
        for d in 0..3 {
            assert!((after.1[d] - before.1[d]).abs() < 2e-14 * scale);
        }
        assert!((after.2 - before.2 - report.kinetic_energy_change).abs() < 2e-15 * scale);
        assert!(report.absolute_mass_defect < 2e-13 * scale);
        assert!(
            report
                .absolute_momentum_defect
                .iter()
                .all(|p| *p < 2e-14 * scale)
        );
        assert_eq!(surface.carry, 0.0);
        for (&key, tile) in &surface.tiles {
            let slot = surface.grid.slot_of(key).unwrap();
            for lane in 0..TILE_CELLS {
                if tile.cells[lane] == Cell3::Interface {
                    let rho = surface.grid.populations(slot, lane).iter().sum::<f64>();
                    assert!(tile.mass[lane] >= 0.0 && tile.mass[lane] <= rho);
                }
            }
        }
    }

    fn check_routes(trace: &MassTransport3, surface: &FreeSurface3) {
        let mut represented = trace.before.clone();
        for route in &trace.routes {
            assert!(route.mass > 0.0);
            *represented.entry(route.from).or_default() -= route.mass;
            *represented.entry(route.to).or_default() += route.mass;
        }
        let remap_defect = surface
            .conversion_momentum_report()
            .unwrap()
            .absolute_mass_defect;
        for location in represented.keys().chain(trace.after.keys()) {
            let routed = represented.get(location).copied().unwrap_or(0.0);
            let actual = trace.after.get(location).copied().unwrap_or(0.0);
            assert!(
                (routed - actual).abs()
                    <= remap_defect
                        + 128.0 * f64::EPSILON * routed.abs().max(actual.abs()).max(1.0),
                "{location:?}: routed {routed}, actual {actual}"
            );
        }
    }

    #[test]
    fn g1_advancing_front_conserves_density_weighted_momentum_and_reports_only_numerical_kinetic_change()
     {
        for scale in [0.4, 1.0, 2.5] {
            for boost in [[0.0; 3], [0.02, -0.01, 0.01]] {
                let incoming = core::array::from_fn(|d| [0.06, 0.01, -0.02][d] + boost[d]);
                let receiving = core::array::from_fn(|d| [-0.025, 0.03, 0.015][d] + boost[d]);
                let mut surface = surface(
                    &[
                        cell([3, 3, 3], 0.8 * scale, incoming),
                        cell([4, 3, 3], 1.4 * scale, receiving),
                    ],
                    [0.0002 * scale, -0.0001 * scale, 0.00005 * scale],
                );
                let rho = surface.sample_at(3, 3, 3).unwrap().density;
                set_mass(&mut surface, [3, 3, 3], 1.2 * rho);
                let before = inventory(&surface, |_| true);
                let densities = density_bits(&surface);
                let mut trace = MassTransport3::new(surface.mass_snapshot(), 1024);
                surface.apply_conversions(Some(&mut trace)).unwrap();
                trace.after = surface.mass_snapshot();
                check_routes(&trace, &surface);
                check_balance(&surface, before);
                let report = surface.conversion_momentum_report().unwrap();
                assert!(report.routes > 1 && report.kinetic_energy_change < 0.0);
                assert_eq!(trace.routes.len(), report.routes);
                assert!(
                    trace
                        .routes
                        .iter()
                        .all(|route| route.from != MassLocation3::Carry
                            && route.to != MassLocation3::Carry)
                );
                assert_eq!(
                    trace.routes.iter().map(|route| route.mass).sum::<f64>(),
                    report.redistributed_mass
                );
                assert_eq!(surface.sample_at(3, 3, 3).unwrap().cell, Cell3::Fluid);
                for (location, bits) in densities {
                    assert_eq!(density_bits(&surface)[&location], bits);
                }
                assert!(
                    (surface.sample_at(4, 3, 3).unwrap().velocity[0] - receiving[0]).abs() > 1e-6
                );
            }
        }
    }

    #[test]
    fn g1_receding_front_carries_the_signed_debt_velocity_and_never_borrows_from_dry_ghosts() {
        let mut surface = surface(
            &[
                cell([3, 3, 3], 0.8, [-0.06, 0.02, 0.0]),
                cell([4, 3, 3], 1.7, [0.02, -0.01, 0.03]),
                cell([3, 4, 3], 0.5, [-0.01, 0.04, 0.02]),
            ],
            [0.0002, -0.0001, 0.00005],
        );
        set_mass(&mut surface, [3, 3, 3], -0.02);
        set_mass(&mut surface, [3, 4, 3], 0.0);
        let before = inventory(&surface, |_| true);
        let old = surface.sample_at(4, 3, 3).unwrap();
        let donor = surface.sample_at(3, 3, 3).unwrap();
        let mut trace = MassTransport3::new(surface.mass_snapshot(), 1024);
        surface.apply_conversions(Some(&mut trace)).unwrap();
        trace.after = surface.mass_snapshot();
        check_routes(&trace, &surface);
        check_balance(&surface, before);
        assert!(surface.sample_at(3, 3, 3).is_none());
        assert_eq!(surface.sample_at(3, 4, 3).unwrap().fill_fraction, 0.0);
        let after = surface.sample_at(4, 3, 3).unwrap();
        let old_mass = 0.4 * 1.7;
        let mass = old_mass - 0.02;
        for d in 0..3 {
            assert!(
                (mass * after.velocity[d] - old_mass * old.velocity[d] + 0.02 * donor.velocity[d])
                    .abs()
                    < 2e-15
            );
        }
        assert!(
            surface
                .conversion_momentum_report()
                .unwrap()
                .kinetic_energy_change
                > 0.0
        );
    }

    #[test]
    fn g1_multiple_sources_preserve_each_disconnected_front_momentum() {
        let mut surface = surface(
            &[
                cell([3, 3, 3], 0.8, [0.05, -0.01, 0.02]),
                cell([4, 3, 3], 1.2, [-0.02, 0.03, 0.01]),
                cell([5, 3, 3], 1.6, [-0.04, -0.02, 0.0]),
                cell([11, 3, 3], 0.6, [0.01, 0.0, -0.03]),
                cell([12, 3, 3], 1.1, [0.0, 0.02, 0.04]),
            ],
            [0.0001, -0.0001, 0.00005],
        );
        for (coordinates, amount) in [([3, 3, 3], 0.12), ([5, 3, 3], 0.18), ([11, 3, 3], 0.1)] {
            let rho = surface
                .sample_at(coordinates[0], coordinates[1], coordinates[2])
                .unwrap()
                .density;
            set_mass(&mut surface, coordinates, rho + amount);
        }
        let before = inventory(&surface, |_| true);
        let left = inventory(&surface, |x| x < 8);
        let right = inventory(&surface, |x| x >= 8);
        surface.apply_conversions(None).unwrap();
        check_balance(&surface, before);
        for (actual, expected) in [
            (inventory(&surface, |x| x < 8), left),
            (inventory(&surface, |x| x >= 8), right),
        ] {
            assert!((actual.0 - expected.0).abs() < 2e-14);
            for d in 0..3 {
                assert!((actual.1[d] - expected.1[d]).abs() < 2e-14);
            }
        }
    }

    #[test]
    fn g3_comoving_conversion_keeps_existing_populations_and_zero_numerical_work() {
        let mut surface = surface(
            &[
                cell([3, 3, 3], 0.8, [0.0; 3]),
                cell([4, 3, 3], 1.6, [0.0; 3]),
            ],
            [0.0; 3],
        );
        set_mass(&mut surface, [3, 3, 3], 0.9);
        let before = surface.clone();
        surface.apply_conversions(None).unwrap();
        check_balance(&surface, inventory(&before, |_| true));
        for (&key, tile) in &before.tiles {
            let old_slot = before.grid.slot_of(key).unwrap();
            let slot = surface.grid.slot_of(key).unwrap();
            for lane in 0..TILE_CELLS {
                if tile.cells[lane] != Cell3::Gas {
                    assert_eq!(
                        surface.grid.populations(slot, lane).map(f64::to_bits),
                        before.grid.populations(old_slot, lane).map(f64::to_bits)
                    );
                }
            }
        }
        assert_eq!(
            surface
                .conversion_momentum_report()
                .unwrap()
                .kinetic_energy_change,
            0.0
        );
    }

    #[test]
    fn g4_debt_with_unrepresentable_momentum_or_no_local_capacity_refuses() {
        for (debt, available) in [(0.1, 0.1), (0.1, 0.101), (0.2, 0.1)] {
            let mut surface = surface(
                &[
                    cell([3, 3, 3], 1.0, [-0.15, 0.0, 0.0]),
                    cell([4, 3, 3], 1.0, [0.15, 0.0, 0.0]),
                ],
                [0.0; 3],
            );
            set_mass(&mut surface, [3, 3, 3], -debt);
            set_mass(&mut surface, [4, 3, 3], available);
            assert!(surface.apply_conversions(None).is_err());
        }
    }

    fn assert_unchanged(before: &FreeSurface3, after: &FreeSurface3) {
        assert_eq!(before.grid.active_keys(), after.grid.active_keys());
        for (&key, tile) in &before.tiles {
            let old_slot = before.grid.slot_of(key).unwrap();
            let slot = after.grid.slot_of(key).unwrap();
            assert_eq!(tile.cells, after.tiles[&key].cells);
            assert_eq!(
                tile.mass.map(f64::to_bits),
                after.tiles[&key].mass.map(f64::to_bits)
            );
            for lane in 0..TILE_CELLS {
                assert_eq!(
                    before.grid.populations(old_slot, lane).map(f64::to_bits),
                    after.grid.populations(slot, lane).map(f64::to_bits)
                );
            }
        }
        assert_eq!(before.carry.to_bits(), after.carry.to_bits());
        assert_eq!(before.conversions, after.conversions);
        assert_eq!(before.conversion_momentum, after.conversion_momentum);
        assert_eq!(before.local_force, after.local_force);
        assert_eq!(
            before.worst_step_violation.to_bits(),
            after.worst_step_violation.to_bits()
        );
    }

    #[test]
    fn g1_actual_native_timestep_advances_a_filled_front_and_zero_pressure_replays_it() {
        let mut cells = Vec::new();
        for z in 0..8 {
            for y in 0..8 {
                for x in 0..4 {
                    let mut initial = cell([x, y, z], 1.0, [0.04, 0.0, 0.0]);
                    initial.fill_fraction = if x == 3 { 0.99 } else { 1.0 };
                    cells.push(initial);
                }
            }
        }
        let mut normal = surface(&cells, [0.0; 3]);
        let mut pressure = normal.clone();
        let mut transported = normal.clone();
        let mass = normal.ledger_mass();
        let gate = CancelGate::new_clock_free();
        normal.step_with_gate(&gate).unwrap();
        pressure
            .step_with_gas_pressure(&gate, 0.2, |_| Ok(0.0))
            .unwrap();
        let mut trace = MassTransport3::new(transported.mass_snapshot(), 8192);
        transported
            .step_inner(&gate, Some(&mut trace), None, None)
            .unwrap();
        trace.after = transported.mass_snapshot();
        check_routes(&trace, &transported);
        assert_unchanged(&normal, &pressure);
        assert_unchanged(&normal, &transported);
        assert_eq!(normal.sample_at(3, 4, 4).unwrap().cell, Cell3::Fluid);
        assert!(normal.sample_at(4, 4, 4).unwrap().fill_fraction > 0.0);
        let report = normal.conversion_momentum_report().unwrap();
        assert!(report.routes > 0 && report.redistributed_mass > 0.0);
        assert!((normal.ledger_mass() - mass).abs() < 2e-12 * mass);
        assert!(report.absolute_momentum_defect.iter().all(|p| *p < 1e-12));
    }

    #[test]
    fn g4_public_normal_and_pressure_steps_rollback_late_conversion_and_cancellation() {
        let mut cells = Vec::new();
        for z in 0..8 {
            for y in 0..8 {
                for x in 0..8 {
                    let mut initial = cell([x, y, z], 1.0, [0.0; 3]);
                    initial.fill_fraction = if [x, y, z] == [3, 3, 3] { 0.4 } else { 1.0 };
                    cells.push(initial);
                }
            }
        }
        let mut surface = surface(&cells, [0.0; 3]);
        // This post-stream overshoot is surrounded entirely by full cells.
        // It needs a wider remap; borrowing from the distant free face is not
        // an adjacent-interface conversion. The error occurs after streaming.
        set_mass(&mut surface, [3, 3, 3], 1.1);
        let before = surface.clone();
        let gate = CancelGate::new_clock_free();
        assert!(surface.step_with_gate(&gate).is_err());
        assert_unchanged(&before, &surface);
        assert!(
            surface
                .step_with_gas_pressure(&gate, 0.2, |_| Ok(0.0))
                .is_err()
        );
        assert_unchanged(&before, &surface);
        surface.carry = 0.01;
        let before = surface.clone();
        assert!(surface.step_with_gate(&gate).is_err());
        assert_unchanged(&before, &surface);
        gate.request();
        assert!(surface.step_with_gate(&gate).is_err());
        assert!(
            surface
                .step_with_gas_pressure(&gate, 0.2, |_| Ok(0.0))
                .is_err()
        );
        assert_unchanged(&before, &surface);
    }
}
