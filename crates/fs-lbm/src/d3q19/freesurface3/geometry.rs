//! Conservative completion of native interface fills for resolved geometry.
//!
//! The ordinary VOF trajectory deliberately retains conversion hysteresis.
//! A physical surface reconstruction instead needs every actual interface
//! mass in [0,rho]. This explicit material-candidate operation moves the real
//! excess or debt on native local mass routes; no geometric clipping, mass
//! floor, heat source, or unlocated carry redistribution is introduced.

use super::*;

impl FreeSurface3 {
    fn interface_mass_and_density(
        &self,
        key: u64,
        lane: usize,
    ) -> Result<(f64, f64), FreeSurfaceError3> {
        let slot = self.grid.slot_of(key).ok_or(FreeSurfaceError3::Transport(
            "resolved interface requires an active cell",
        ))?;
        let density = self.grid.populations(slot, lane).iter().sum::<f64>();
        let mass = self.tiles[&key].mass[lane];
        if !density.is_finite() || density <= 0.0 || !mass.is_finite() {
            return Err(FreeSurfaceError3::Transport(
                "resolved interface needs finite mass and positive density",
            ));
        }
        Ok((mass, density))
    }

    /// Bring a private candidate's true interface masses into [0,rho].
    ///
    /// First transfer excess/debt to adjacent interfaces with available
    /// capacity. Recipients cannot be driven outside their physical bounds.
    /// If that one-ring storage is insufficient, the original native closure
    /// cascade exposes the next interface layer using an exact conversion
    /// threshold; subsequent local transfers finish its residuals. Every mass
    /// movement uses the existing companion-field route budget, and each
    /// unfinished round must consume a route. Unlocated carry refuses.
    ///
    /// The original trace.before remains untouched. The caller must advect
    /// enthalpy with the amended routes and trace.after, then admit source and
    /// geometry at the endpoint before publishing the private candidate. As
    /// with step_with_mass_transport, a refusal can leave THIS private clone
    /// partially advanced; accepted material state is owned by the caller.
    pub(crate) fn resolve_interface_fills(
        &mut self,
        trace: &mut MassTransport3,
        gate: &CancelGate,
    ) -> Result<(), FreeSurfaceError3> {
        check_cancel(gate)?;
        if self.redistribution != ConversionRedistribution3::Local || self.carry != 0.0 {
            return Err(FreeSurfaceError3::Transport(
                "resolved interface geometry requires local routes and no unlocated carry",
            ));
        }
        loop {
            let mut defects = Vec::new();
            for (&key, tile) in &self.tiles {
                check_cancel(gate)?;
                for (lane, &cell) in tile.cells.iter().enumerate() {
                    if cell == Cell3::Interface {
                        let (mass, density) = self.interface_mass_and_density(key, lane)?;
                        if mass < 0.0 || mass > density {
                            defects.push((key, lane));
                        }
                    }
                }
            }
            if defects.is_empty() {
                check_cancel(gate)?;
                trace.after = self.mass_snapshot();
                return Ok(());
            }
            let routes_before = trace.routes.len();
            for &(key, lane) in &defects {
                check_cancel(gate)?;
                let (mass, density) = self.interface_mass_and_density(key, lane)?;
                let excess = mass > density;
                let need = if excess { mass - density } else { -mass };
                if need <= 0.0 {
                    continue;
                }
                let mut recipients = [(0, 0, 0.0, 0.0, 0.0); Q3 - 1];
                let mut count = 0;
                let mut available = 0.0;
                for q in 1..Q3 {
                    let Neighbor::Active(other_key, other_lane, Cell3::Interface) =
                        self.neighbor_state(key, lane, q)
                    else {
                        continue;
                    };
                    let (other_mass, other_density) =
                        self.interface_mass_and_density(other_key, other_lane)?;
                    let capacity = if excess {
                        other_density - other_mass
                    } else {
                        other_mass
                    };
                    if capacity <= 0.0 {
                        continue;
                    }
                    recipients[count] =
                        (other_key, other_lane, other_mass, other_density, capacity);
                    count += 1;
                    available += capacity;
                }
                if !available.is_finite() {
                    return Err(FreeSurfaceError3::Transport(
                        "resolved interface recipient capacity is not representable",
                    ));
                }
                if count == 0 {
                    continue;
                }
                // Share by available storage, rather than selecting the first
                // lattice direction. Fixed arrays and a residual on the final
                // recipient retain bounded work and deterministic conservation.
                let amount = need.min(available);
                let mut remaining = amount;
                for (index, &(other_key, other_lane, other_mass, other_density, capacity)) in
                    recipients[..count].iter().enumerate()
                {
                    let transfer = if amount == available {
                        capacity
                    } else if index + 1 == count {
                        remaining.min(capacity)
                    } else {
                        (capacity * (amount / available))
                            .min(capacity)
                            .min(remaining)
                    };
                    let donor = MassLocation3::Cell(key, lane);
                    let neighbor = MassLocation3::Cell(other_key, other_lane);
                    if excess {
                        trace.exchange(donor, neighbor, transfer)?;
                    } else {
                        trace.exchange(neighbor, donor, transfer)?;
                    }
                    self.tiles.get_mut(&key).expect("active donor").mass[lane] +=
                        if excess { -transfer } else { transfer };
                    self.tiles
                        .get_mut(&other_key)
                        .expect("active recipient")
                        .mass[other_lane] = if transfer == capacity {
                        if excess { other_density } else { 0.0 }
                    } else if excess {
                        other_mass + transfer
                    } else {
                        other_mass - transfer
                    };
                    remaining -= transfer;
                }
                if amount == available {
                    // Each recipient took its entire representable capacity.
                    remaining = 0.0;
                }
                if amount == need && remaining == 0.0 {
                    // The complete defect is carried by the recorded edges.
                    // Assign its exact bound to avoid roundoff-only cut films;
                    // no mass residual is discarded by a geometry-only clamp.
                    self.tiles.get_mut(&key).expect("active donor").mass[lane] =
                        if excess { density } else { 0.0 };
                }
            }
            let mut unfinished = false;
            for &(key, lane) in &defects {
                let (mass, density) = self.interface_mass_and_density(key, lane)?;
                unfinished |= mass < 0.0 || mass > density;
            }
            if unfinished {
                check_cancel(gate)?;
                self.apply_conversions_with_tolerance(Some(trace), 0.0)?;
                if self.carry != 0.0 {
                    return Err(FreeSurfaceError3::Transport(
                        "resolved interface conversion has no local mass recipient",
                    ));
                }
            }
            if trace.routes.len() == routes_before {
                return Err(FreeSurfaceError3::Transport(
                    "resolved interface mass projection made no progress",
                ));
            }
        }
    }
}

#[cfg(test)]
impl FreeSurface3 {
    /// Manufactured native-route endpoint shared with the real enthalpy test.
    /// A finite flux oversends one interface cell into its adjacent peer;
    /// the other detached component must never participate in its correction.
    pub(crate) fn fill_projection_fixture(
        route_budget: usize,
    ) -> (Self, MassTransport3, MassLocation3, MassLocation3) {
        let mut surface = Self::new(
            24,
            12,
            12,
            0.8,
            [0.0; 3],
            0.002,
            ContactModel3::Neutral,
            |x, y, z| {
                (2..7).contains(&y)
                    && (2..7).contains(&z)
                    && ((2..7).contains(&x) || (17..22).contains(&x))
            },
        )
        .unwrap();
        let locate = |p: [i64; 3]| {
            let (slot, lane) = surface.grid.resolve_source(p[0], p[1], p[2]).unwrap();
            (surface.grid.active_keys()[slot], lane)
        };
        let (a_key, a_lane) = locate([6, 4, 4]);
        let (b_key, b_lane) = locate([6, 5, 4]);
        let a = MassLocation3::Cell(a_key, a_lane);
        let b = MassLocation3::Cell(b_key, b_lane);
        // Other neighboring interfaces are full: this manufactured endpoint
        // has exactly one local recipient, giving an independent two-node
        // analytic enthalpy solution after the actual native correction.
        for q in 1..Q3 {
            if let Neighbor::Active(key, lane, Cell3::Interface) =
                surface.neighbor_state(a_key, a_lane, q)
            {
                let density = surface.interface_mass_and_density(key, lane).unwrap().1;
                surface.tiles.get_mut(&key).unwrap().mass[lane] = density;
            }
        }
        let density = surface.interface_mass_and_density(a_key, a_lane).unwrap().1;
        surface.tiles.get_mut(&a_key).unwrap().mass[a_lane] = 0.5 * density;
        surface.tiles.get_mut(&b_key).unwrap().mass[b_lane] = 0.5 * density;
        let mut trace = MassTransport3::new(surface.mass_snapshot(), route_budget);
        // A dyadic residual isolates transport physics from subtraction ulps.
        let transfer = surface.tiles[&b_key].mass[b_lane] + 1.0 / 4096.0;
        trace.exchange(b, a, transfer).unwrap();
        surface.tiles.get_mut(&a_key).unwrap().mass[a_lane] += transfer;
        surface.tiles.get_mut(&b_key).unwrap().mass[b_lane] -= transfer;
        trace.after = surface.mass_snapshot();
        (surface, trace, a, b)
    }
}

#[cfg(test)]
mod tests {
    use super::super::thermocapillary::HeightGraphSurface3;
    use super::*;

    fn assert_balance(trace: &MassTransport3) {
        let mut balance = trace.before.clone();
        for route in &trace.routes {
            *balance.entry(route.from).or_default() -= route.mass;
            *balance.entry(route.to).or_default() += route.mass;
        }
        for (&key, &mass) in &trace.after {
            *balance.entry(key).or_default() -= mass;
        }
        for (key, residual) in balance {
            assert!(residual.abs() < 4e-13, "{key:?}: {residual}");
        }
    }

    fn assert_resolved(surface: &FreeSurface3) {
        assert_eq!(surface.carry, 0.0);
        for (&key, tile) in &surface.tiles {
            for (lane, &cell) in tile.cells.iter().enumerate() {
                if cell == Cell3::Interface {
                    let (mass, density) = surface.interface_mass_and_density(key, lane).unwrap();
                    assert!(
                        mass >= 0.0 && mass <= density,
                        "{key}/{lane}: {mass}/{density}"
                    );
                }
            }
        }
    }

    #[test]
    fn g1_resolved_fills_transport_real_excess_and_debt_without_cross_drop_mixing() {
        let (mut surface, mut trace, a, b) = FreeSurface3::fill_projection_fixture(1000);
        let before = trace.before.clone();
        let populations = surface.grid.state_bits();
        let mass = surface.ledger_mass();
        assert!(trace.after[&a] > 1.0 && trace.after[&b] < 0.0);
        surface
            .resolve_interface_fills(&mut trace, &CancelGate::new_clock_free())
            .unwrap();
        assert_resolved(&surface);
        assert_balance(&trace);
        assert_eq!(trace.before, before);
        assert_eq!(surface.grid.state_bits(), populations);
        assert!((surface.ledger_mass() - mass).abs() < 1e-12);
        assert_eq!(trace.after[&b], 0.0);
        assert_eq!(trace.routes.len(), 2);
        for route in &trace.routes {
            let p = FreeSurface3::transport_coordinates(route.from).unwrap();
            let q = FreeSurface3::transport_coordinates(route.to).unwrap();
            assert!(p[0] < 14 && q[0] < 14);
            assert!((0..3).all(|d| p[d].abs_diff(q[d]) <= 1));
        }
        for (&location, &mass) in &before {
            if FreeSurface3::transport_coordinates(location).is_some_and(|p| p[0] >= 14) {
                assert_eq!(trace.after[&location].to_bits(), mass.to_bits());
            }
        }
        let resolved = trace.after.clone();
        surface
            .resolve_interface_fills(&mut trace, &CancelGate::new_clock_free())
            .unwrap();
        assert_eq!(trace.routes.len(), 2);
        assert_eq!(trace.after, resolved);
    }

    #[test]
    fn g1_resolved_fills_reuse_native_front_conversion_when_local_capacity_is_absent() {
        let mut surface = FreeSurface3::new(
            12,
            12,
            12,
            0.8,
            [0.0; 3],
            0.002,
            ContactModel3::Neutral,
            |x, y, z| [x, y, z] == [6, 6, 6],
        )
        .unwrap();
        let (slot, lane) = surface.grid.resolve_source(6, 6, 6).unwrap();
        let key = surface.grid.active_keys()[slot];
        let density = surface.grid.populations(slot, lane).iter().sum::<f64>();
        surface.tiles.get_mut(&key).unwrap().mass[lane] = density + EPS / 4.0;
        let mut trace = MassTransport3::new(surface.mass_snapshot(), 1000);
        surface
            .resolve_interface_fills(&mut trace, &CancelGate::new_clock_free())
            .unwrap();
        assert_resolved(&surface);
        assert_balance(&trace);
        assert_eq!(surface.tiles[&key].cells[lane], Cell3::Fluid);
        assert_eq!(surface.conversions.to_fluid, 1);
        assert!(!trace.routes.is_empty());
        assert!(
            trace
                .routes
                .iter()
                .all(|route| route.from == MassLocation3::Cell(key, lane))
        );
    }

    #[test]
    fn g0_resolved_fills_cancel_or_exhaust_budget_before_partial_publication() {
        let (mut surface, mut trace, _, _) = FreeSurface3::fill_projection_fixture(1);
        let before = surface.mass_snapshot();
        let gate = CancelGate::new_clock_free();
        assert!(surface.resolve_interface_fills(&mut trace, &gate).is_err());
        assert_eq!(surface.mass_snapshot(), before);
        assert_eq!(trace.routes.len(), 1);
        let (mut surface, mut trace, _, _) = FreeSurface3::fill_projection_fixture(1000);
        let before = surface.mass_snapshot();
        gate.request();
        assert!(matches!(
            surface.resolve_interface_fills(&mut trace, &gate),
            Err(FreeSurfaceError3::Cancelled)
        ));
        assert_eq!(surface.mass_snapshot(), before);
        assert_eq!(trace.routes.len(), 1);
    }

    fn curved_surface(sign: f64) -> FreeSurface3 {
        let mut surface = FreeSurface3::new(
            12,
            12,
            12,
            0.8,
            [0.0; 3],
            0.002,
            ContactModel3::Neutral,
            |_, _, z| z < 6,
        )
        .unwrap();
        for (&key, tile) in &mut surface.tiles {
            let slot = surface.grid.slot_of(key).unwrap();
            for lane in 0..TILE_CELLS {
                let (x, y, z) = FreeSurface3::global_of(key, lane);
                let mut fill = 0.0;
                for a in 0..8 {
                    for b in 0..8 {
                        let u = x as f64 + (a as f64 + 0.5) / 8.0 - 6.0;
                        let v = y as f64 + (b as f64 + 0.5) / 8.0 - 6.0;
                        let height = 5.6 + sign * 0.02 * (u * u + v * v);
                        fill += (height - z as f64).clamp(0.0, 1.0) / 64.0;
                    }
                }
                let density = surface.grid.populations(slot, lane).iter().sum::<f64>();
                tile.mass[lane] = density * fill;
                tile.cells[lane] = if fill == 0.0 {
                    Cell3::Gas
                } else if fill == 1.0 {
                    Cell3::Fluid
                } else {
                    Cell3::Interface
                };
            }
        }
        let mut closure = Vec::new();
        for (&key, tile) in &surface.tiles {
            for (lane, &cell) in tile.cells.iter().enumerate() {
                if cell == Cell3::Fluid
                    && (1..Q3).any(|q| {
                        matches!(
                            surface.neighbor_state(key, lane, q),
                            Neighbor::Active(_, _, Cell3::Gas)
                        )
                    })
                {
                    closure.push((key, lane));
                }
            }
        }
        for (key, lane) in closure {
            surface.tiles.get_mut(&key).unwrap().cells[lane] = Cell3::Interface;
        }
        surface
    }

    fn tensions(surface: &FreeSurface3) -> BTreeMap<(u64, usize), f64> {
        surface
            .mass_snapshot()
            .keys()
            .filter_map(|&location| match location {
                MassLocation3::Cell(key, lane) => Some(((key, lane), 0.002)),
                MassLocation3::Carry => None,
            })
            .collect()
    }

    #[test]
    fn g1_concave_and_convex_capillary_graphs_keep_admissible_raw_geometry_after_motion() {
        let gate = CancelGate::new_clock_free();
        let graph = HeightGraphSurface3 {
            axis: 2,
            outward_sign: 1,
            maximum_slope: 0.5,
            minimum_fill: 0.001,
        };
        let mut repaired_steps = 0;
        for sign in [-1.0, 1.0] {
            let mut surface = curved_surface(sign);
            let initial_mass = surface.ledger_mass();
            let mut raw_refusals = 0;
            let mut worst_mass_defect = 0.0_f64;
            for _ in 0..6 {
                let mut sigma = tensions(&surface);
                surface
                    .height_graph_marangoni_forces(graph, &mut sigma, &gate)
                    .unwrap();
                surface.rebind_surface_tensions(sigma, &gate).unwrap();
                let tau = tensions(&surface)
                    .into_keys()
                    .map(|location| (location, 0.8))
                    .collect();
                let mut trace = surface
                    .step_with_mass_transport(&gate, &tau, None, 100_000)
                    .unwrap();
                let routes_before = trace.routes.len();
                let mut raw_sigma = tensions(&surface);
                raw_refusals += usize::from(
                    surface
                        .height_graph_marangoni_forces(graph, &mut raw_sigma, &gate)
                        .is_err(),
                );
                surface.resolve_interface_fills(&mut trace, &gate).unwrap();
                repaired_steps += usize::from(trace.routes.len() > routes_before);
                assert_resolved(&surface);
                assert_balance(&trace);
                let mut sigma = tensions(&surface);
                surface
                    .height_graph_marangoni_forces(graph, &mut sigma, &gate)
                    .unwrap();
                let mass_defect = (surface.ledger_mass() / initial_mass - 1.0).abs();
                worst_mass_defect = worst_mass_defect.max(mass_defect);
                assert!(mass_defect < 3e-13);
            }
            eprintln!(
                "curvature sign {sign}: {raw_refusals} raw geometry refusals resolved; relative mass defect {worst_mass_defect:e}"
            );
            assert!(
                raw_refusals > 0,
                "each curved graph must reproduce the strict-fill blockage"
            );
        }
        assert!(
            repaired_steps > 0,
            "fixture must exercise actual post-streaming fill correction"
        );
    }
}
