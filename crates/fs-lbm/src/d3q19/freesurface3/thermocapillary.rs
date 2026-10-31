//! Resolved surface traction on the original VOF/Guo owner.
//!
//! The planar reduction keeps one Cartesian normal and projected area. The
//! height-graph reduction reconstructs both tangent directions and actual
//! surface area from complete VOF columns, including multi-cell cut bands.
//! Local detached interfaces use complete height stencils, volume-matching cut
//! planes and connected wet-source gradients. All paths use the existing Guo
//! force owner; the normal Laplace load remains exclusively in the original
//! pressure boundary. Surfactant transport and surface-energy conversion are
//! outside these geometric stress reconstructions.

use super::*;

mod graph;
pub(crate) mod local;
pub(crate) use local::LocalSurface3;

/// All lengths are in lattice cell widths; material ingress supplies SI units.
#[derive(Clone, Copy, Debug)]
pub(crate) struct PlanarSurface3 {
    pub axis: usize,
    pub outward_sign: i8,
    pub position: f64,
    pub maximum_departure: f64,
    pub minimum_fill: f64,
}

fn invalid(reason: &'static str) -> FreeSurfaceError3 {
    FreeSurfaceError3::Transport(reason)
}

impl FreeSurface3 {
    pub(super) fn admit_surface_tensions(
        &self,
        field: &BTreeMap<(u64, usize), f64>,
        gate: &CancelGate,
    ) -> Result<(), FreeSurfaceError3> {
        let mut count = 0;
        for (&key, tile) in &self.tiles {
            check_cancel(gate)?;
            for (lane, cell) in tile.cells.iter().enumerate() {
                if *cell == Cell3::Gas {
                    continue;
                }
                let sigma = *field
                    .get(&(key, lane))
                    .ok_or(invalid("surface tension must cover every wet cell"))?;
                if !sigma.is_finite()
                    || sigma <= 0.0
                    || sigma * Self::CURVATURE_LIMIT / CS2_3 >= 1.0
                {
                    return Err(invalid(
                        "local surface tension exceeds native pressure-boundary domain",
                    ));
                }
                if !crate::admit_explicit_capillary_coefficient(sigma) {
                    return Err(invalid(
                        "local surface tension exceeds conservative explicit capillary-wave bound",
                    ));
                }
                count += 1;
            }
        }
        if count != field.len() {
            return Err(invalid("surface tension contains a foreign cell"));
        }
        Ok(())
    }

    /// Rebind the pressure coefficient without advancing momentum or heat.
    /// The enclosing material owner commits this with the rest of its candidate.
    pub(crate) fn rebind_surface_tensions(
        &mut self,
        field: BTreeMap<(u64, usize), f64>,
        gate: &CancelGate,
    ) -> Result<(), FreeSurfaceError3> {
        self.admit_surface_tensions(&field, gate)?;
        check_cancel(gate)?;
        self.local_sigma = if field.values().all(|&sigma| sigma == self.sigma) {
            None
        } else {
            Some(field)
        };
        Ok(())
    }

    /// Actual accepted capillary coefficient in lattice units at a wet cell.
    /// The physical source and SI conversion remain with the material adapter.
    #[must_use]
    pub fn surface_tension_at(&self, x: usize, y: usize, z: usize) -> Option<f64> {
        let (slot, lane) = self.grid.resolve_source(
            i64::try_from(x).ok()?,
            i64::try_from(y).ok()?,
            i64::try_from(z).ok()?,
        )?;
        let key = self.grid.active_keys()[slot];
        if self.tiles.get(&key)?.cells[lane] == Cell3::Gas {
            return None;
        }
        self.local_sigma
            .as_ref()
            .map_or(Some(self.sigma), |s| s.get(&(key, lane)).copied())
    }

    /// Sparse extra Guo force, with no normal component and no mutation.
    /// Source sigma is evaluated at each actual surface-cell temperature;
    /// interior normal gradients cannot become tangential traction.
    pub(crate) fn planar_marangoni_forces(
        &self,
        plane: PlanarSurface3,
        sigma: &BTreeMap<(u64, usize), f64>,
        gate: &CancelGate,
    ) -> Result<BTreeMap<(u64, usize), [f64; 3]>, FreeSurfaceError3> {
        self.admit_surface_tensions(sigma, gate)?;
        if plane.axis >= 3
            || ![-1, 1].contains(&plane.outward_sign)
            || !plane.position.is_finite()
            || !plane.maximum_departure.is_finite()
            || !(0.0..0.5).contains(&plane.maximum_departure)
            || !plane.minimum_fill.is_finite()
            || plane.minimum_fill <= 0.0
            || plane.minimum_fill > 1.0
            || !self.solid.is_empty()
            || self.planar_contact.is_some()
            || self.planar_contact_field.is_some()
            || self.contact != ContactModel3::Neutral
            || self.carry != 0.0
        {
            return Err(invalid(
                "planar thermocapillarity needs bounded Cartesian tank geometry and neutral walls",
            ));
        }
        let (nx, ny, nz) = self.grid.tile_dims();
        let size = [nx * TILE, ny * TILE, nz * TILE];
        if plane.position <= 1.0 || plane.position >= (size[plane.axis] - 1) as f64 {
            return Err(invalid(
                "planar surface needs resolved liquid depth and gas clearance",
            ));
        }
        let surface_layer = if plane.outward_sign == 1 {
            plane.position.ceil() as usize - 1
        } else {
            plane.position.floor() as usize
        };
        let tangent: Vec<usize> = (0..3).filter(|&d| d != plane.axis).collect();
        let mut surface = BTreeMap::new();
        for b in 0..size[tangent[1]] {
            check_cancel(gate)?;
            for a in 0..size[tangent[0]] {
                check_cancel(gate)?;
                let mut top: Option<((u64, usize), usize, f64)> = None;
                let mut in_gas = false;
                for depth in 0..size[plane.axis] {
                    if depth % TILE_CELLS == 0 {
                        check_cancel(gate)?;
                    }
                    let axis_position = if plane.outward_sign == 1 {
                        depth
                    } else {
                        size[plane.axis] - 1 - depth
                    };
                    let mut p = [0; 3];
                    p[plane.axis] = axis_position;
                    p[tangent[0]] = a;
                    p[tangent[1]] = b;
                    let location = self
                        .grid
                        .resolve_source(p[0] as i64, p[1] as i64, p[2] as i64)
                        .map(|(slot, lane)| (self.grid.active_keys()[slot], slot, lane));
                    // Geometry uses actual tracked fill, not the native
                    // curvature diagnostic's clamped rendering fraction.
                    let fill = match location {
                        Some((key, slot, lane)) => match self.tiles[&key].cells[lane] {
                            Cell3::Gas => 0.0,
                            Cell3::Fluid => 1.0,
                            Cell3::Interface => {
                                let rho = self.grid.populations(slot, lane).iter().sum::<f64>();
                                if !rho.is_finite() || rho <= 0.0 {
                                    return Err(invalid(
                                        "invalid density in planar surface geometry",
                                    ));
                                }
                                self.tiles[&key].mass[lane] / rho
                            }
                        },
                        None => 0.0,
                    };
                    if !fill.is_finite() || fill < 0.0 || fill > 1.0 + 64.0 * f64::EPSILON {
                        return Err(invalid(
                            "planar surface requires finite nonnegative geometric fill without overfill",
                        ));
                    }
                    if fill > 0.0 {
                        if in_gas {
                            return Err(invalid(
                                "planar surface has an overhang, overfill or disconnected liquid",
                            ));
                        }
                        if let Some((_, _, previous_fill)) = top {
                            if (previous_fill - 1.0_f64).abs() > 64.0 * f64::EPSILON {
                                return Err(invalid(
                                    "planar surface has more than one partial cell per column",
                                ));
                            }
                        }
                        let (key, _, lane) = location.expect("positive geometric fill");
                        top = Some(((key, lane), axis_position, fill));
                    } else {
                        in_gas = true;
                    }
                }
                let (location, position, fill) = top.ok_or(invalid(
                    "planar surface requires liquid in every tank column",
                ))?;
                let height = if plane.outward_sign == 1 {
                    position as f64 + fill
                } else {
                    position as f64 + 1.0 - fill
                };
                if !in_gas
                    || position != surface_layer
                    || fill < plane.minimum_fill
                    || (height - plane.position).abs() > plane.maximum_departure
                {
                    return Err(invalid(
                        "free surface left its admitted planar height or fill envelope",
                    ));
                }
                surface.insert([a, b], (location, fill, sigma[&location]));
            }
        }
        let mut forces = BTreeMap::new();
        for (&p, &(location, fill, value)) in &surface {
            check_cancel(gate)?;
            let mut force = [0.0; 3];
            for d in 0..2 {
                let previous = p[d].checked_sub(1).map(|index| {
                    let mut q = p;
                    q[d] = index;
                    surface[&q].2
                });
                let next = (p[d] + 1 < size[tangent[d]]).then(|| {
                    let mut q = p;
                    q[d] += 1;
                    surface[&q].2
                });
                // Centered interior difference; one-sided at the actual tank
                // walls. No unavailable gas temperature or ghost sigma is used.
                let gradient = match (previous, next) {
                    (Some(m), Some(p)) => (p - m) / 2.0,
                    (None, Some(p)) => p - value,
                    (Some(m), None) => value - m,
                    (None, None) => return Err(invalid("unresolved planar surface tangent")),
                };
                force[tangent[d]] = gradient / fill;
                if !force[tangent[d]].is_finite() {
                    return Err(invalid("planar surface traction is not representable"));
                }
            }
            if force != [0.0; 3] {
                forces.insert(location, force);
            }
        }
        Ok(forces)
    }
}

/// Resolved single-valued tank interface, with no small-slope approximation.
/// The maximum slope bounds |grad h|, not the departure from one lattice layer.
#[derive(Clone, Copy, Debug)]
pub(crate) struct HeightGraphSurface3 {
    pub axis: usize,
    pub outward_sign: i8,
    pub maximum_slope: f64,
    pub minimum_fill: f64,
}

#[derive(Clone, Copy)]
struct SurfaceColumn {
    height: f64,
    volume: f64,
    sigma: f64,
    force: [f64; 3],
}

// Centered differences inside the tank and second-order one-sided differences
// at its edges. Difference first, so a large common datum cancels exactly.
fn graph_derivative(
    p: [usize; 2],
    axis: usize,
    size: [usize; 2],
    sample: impl Fn([usize; 2]) -> f64,
) -> f64 {
    let at = |offset: isize| {
        let mut q = p;
        q[axis] = (q[axis] as isize + offset) as usize;
        sample(q)
    };
    if p[axis] == 0 {
        2.0 * (at(1) - at(0)) - 0.5 * (at(2) - at(0))
    } else if p[axis] + 1 == size[axis] {
        2.0 * (at(0) - at(-1)) - 0.5 * (at(0) - at(-2))
    } else {
        0.5 * (at(1) - at(-1))
    }
}

impl FreeSurface3 {
    /// Height-column Marangoni traction on a resolved graph q=h(u,v).
    ///
    /// Column sigma is the liquid-volume-weighted value in its interface band
    /// (Seric, Afkhami & Kondic, JCP 352 (2018), sec. 3.1,
    /// https://doi.org/10.1016/j.jcp.2017.10.008). With a=(h_u,h_v),
    /// b=(sigma_u,sigma_v), the graph metric gives
    /// grad_s sigma = (b-a*(a.b)/(1+a.a), (a.b)/(1+a.a)).
    /// Its traction is integrated over area sqrt(1+a.a), then distributed over
    /// the band's occupied volume. Thus sum(fill*force) is the patch traction
    /// and force.normal=0, without duplicating the Laplace pressure load.
    ///
    /// The derived column sigma also replaces the interface pressure
    /// coefficient, including zero-mass ghosts. Interior material temperatures
    /// are unchanged. The caller owns this temporary field and publishes it
    /// only with the complete thermal/momentum candidate.
    pub(crate) fn height_graph_marangoni_forces(
        &self,
        graph: HeightGraphSurface3,
        sigma: &mut BTreeMap<(u64, usize), f64>,
        gate: &CancelGate,
    ) -> Result<BTreeMap<(u64, usize), [f64; 3]>, FreeSurfaceError3> {
        self.admit_surface_tensions(sigma, gate)?;
        let patches = self.height_graph_geometry(graph, gate)?;
        let (nx, ny, nz) = self.grid.tile_dims();
        let size = [nx * TILE, ny * TILE, nz * TILE];
        let tangent = match graph.axis {
            0 => [1, 2],
            1 => [0, 2],
            _ => [0, 1],
        };
        let shape = [size[tangent[0]], size[tangent[1]]];
        let mut columns: Vec<_> = patches
            .iter()
            .map(|patch| SurfaceColumn {
                height: patch.centroid[graph.axis],
                volume: patch.volume,
                sigma: 0.0,
                force: [0.0; 3],
            })
            .collect();
        let mut references = vec![None; columns.len()];
        // Preserve the material band's occupied-volume average. Difference
        // first so uniform sigma produces exactly zero tangential force.
        for (&key, tile) in &self.tiles {
            check_cancel(gate)?;
            let slot = self.grid.slot_of(key).expect("active interface tile");
            for (lane, &cell) in tile.cells.iter().enumerate() {
                if cell != Cell3::Interface || tile.mass[lane] <= 0.0 {
                    continue;
                }
                let p = Self::transport_coordinates(MassLocation3::Cell(key, lane))
                    .ok_or(invalid("unlocated height-graph interface"))?;
                let index = p[tangent[0]] + shape[0] * p[tangent[1]];
                let fill = (tile.mass[lane]
                    / self.grid.populations(slot, lane).iter().sum::<f64>())
                .min(1.0);
                let value = sigma[&(key, lane)];
                let reference = *references[index].get_or_insert(value);
                columns[index].sigma += fill * (value - reference);
            }
        }
        for (column, reference) in columns.iter_mut().zip(references) {
            column.sigma = reference.ok_or(invalid("height graph has no liquid interface"))?
                + column.sigma / column.volume;
            if !column.sigma.is_finite() || column.sigma <= 0.0 {
                return Err(invalid("height-graph surface tension is not representable"));
            }
        }
        for b in 0..shape[1] {
            for a in 0..shape[0] {
                check_cancel(gate)?;
                let p = [a, b];
                let h = |q: [usize; 2]| columns[q[0] + shape[0] * q[1]].height;
                let s = |q: [usize; 2]| columns[q[0] + shape[0] * q[1]].sigma;
                let slope = [
                    graph_derivative(p, 0, shape, h),
                    graph_derivative(p, 1, shape, h),
                ];
                let gradient = [
                    graph_derivative(p, 0, shape, s),
                    graph_derivative(p, 1, shape, s),
                ];
                let metric = 1.0 + slope[0] * slope[0] + slope[1] * slope[1];
                let normal_component = (slope[0] * gradient[0] + slope[1] * gradient[1]) / metric;
                let column = &mut columns[a + shape[0] * b];
                let scale = patches[a + shape[0] * b].area / column.volume;
                let mut force = [0.0; 3];
                force[tangent[0]] = (gradient[0] - slope[0] * normal_component) * scale;
                force[tangent[1]] = (gradient[1] - slope[1] * normal_component) * scale;
                force[graph.axis] = normal_component * scale;
                if force.iter().any(|f| !f.is_finite()) {
                    return Err(invalid("height-graph traction is not representable"));
                }
                column.force = force;
            }
        }
        let mut forces = BTreeMap::new();
        for (&key, tile) in &self.tiles {
            check_cancel(gate)?;
            for (lane, &cell) in tile.cells.iter().enumerate() {
                if cell != Cell3::Interface {
                    continue;
                }
                let p = Self::transport_coordinates(MassLocation3::Cell(key, lane))
                    .ok_or(invalid("unlocated height-graph interface"))?;
                let column = columns[p[tangent[0]] + shape[0] * p[tangent[1]]];
                sigma.insert((key, lane), column.sigma);
                if tile.mass[lane] > 0.0 && column.force != [0.0; 3] {
                    forces.insert((key, lane), column.force);
                }
            }
        }
        Ok(forces)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn field(
        surface: &FreeSurface3,
        value: impl Fn([usize; 3]) -> f64,
    ) -> BTreeMap<(u64, usize), f64> {
        surface
            .mass_snapshot()
            .keys()
            .filter_map(|location| {
                let MassLocation3::Cell(key, lane) = *location else {
                    return None;
                };
                let p = FreeSurface3::transport_coordinates(*location).unwrap();
                Some(((key, lane), value(p)))
            })
            .collect()
    }

    #[test]
    fn g1_planar_marangoni_sign_tangent_and_occupied_volume_are_geometric() {
        let gate = CancelGate::new_clock_free();
        for axis in 0..3 {
            for sign in [-1, 1] {
                let surface = FreeSurface3::new(
                    8,
                    8,
                    8,
                    0.8,
                    [0.0; 3],
                    0.01,
                    ContactModel3::Neutral,
                    |x, y, z| {
                        let p = [x, y, z][axis];
                        if sign == 1 { p < 4 } else { p >= 4 }
                    },
                )
                .unwrap();
                let tangent: Vec<usize> = (0..3).filter(|&d| d != axis).collect();
                let plane = PlanarSurface3 {
                    axis,
                    outward_sign: sign,
                    position: if sign == 1 { 3.5 } else { 4.5 },
                    maximum_departure: 1e-10,
                    minimum_fill: 0.25,
                };
                // The normal slope must not become tangential traction.
                let sigma = field(&surface, |p| {
                    0.01 + 0.0002 * p[tangent[0]] as f64 - 0.0001 * p[tangent[1]] as f64
                        + 0.0003 * p[axis] as f64
                });
                let forces = surface
                    .planar_marangoni_forces(plane, &sigma, &gate)
                    .unwrap();
                assert_eq!(forces.len(), 64);
                let mut total = [0.0; 3];
                for (&(key, lane), force) in &forces {
                    let fill = surface.fill(key, lane);
                    assert!((fill - 0.5).abs() < 1e-14);
                    assert_eq!(force[axis], 0.0);
                    assert!((force[tangent[0]] - 0.0002 / fill).abs() < 1e-16);
                    assert!((force[tangent[1]] + 0.0001 / fill).abs() < 1e-16);
                    for d in 0..3 {
                        total[d] += fill * force[d];
                    }
                }
                assert!((total[tangent[0]] - 64.0 * 0.0002).abs() < 1e-15);
                assert!((total[tangent[1]] + 64.0 * 0.0001).abs() < 1e-15);
                let only_normal = field(&surface, |p| 0.01 + 0.0003 * p[axis] as f64);
                assert!(
                    surface
                        .planar_marangoni_forces(plane, &only_normal, &gate)
                        .unwrap()
                        .is_empty()
                );
            }
        }
    }

    #[test]
    fn g1_local_sigma_uses_original_laplace_boundary_and_constant_bits() {
        let gate = CancelGate::new_clock_free();
        let original = FreeSurface3::new(
            8,
            8,
            8,
            0.8,
            [0.0; 3],
            0.002,
            ContactModel3::Neutral,
            |x, y, z| {
                let r = [x, y, z].map(|p| p as f64 - 3.5);
                r.iter().map(|v| v * v).sum::<f64>() < 8.0
            },
        )
        .unwrap();
        let mut reference = original.clone();
        reference.refresh_fill(&gate).unwrap();
        let mut varying = reference.clone();
        let values = field(&varying, |p| 0.003 + 0.0001 * p[0] as f64);
        varying
            .rebind_surface_tensions(values.clone(), &gate)
            .unwrap();
        let mut saw_curvature = false;
        for (&(key, lane), &sigma) in &values {
            if varying.tiles[&key].cells[lane] != Cell3::Interface {
                continue;
            }
            let kappa = varying.curvature(key, lane);
            saw_curvature |= kappa.abs() > 1e-5;
            assert_eq!(
                varying.reference_density(key, lane),
                sigma.mul_add(kappa / CS2_3, 1.0)
            );
            let p = FreeSurface3::transport_coordinates(MassLocation3::Cell(key, lane)).unwrap();
            assert_eq!(varying.surface_tension_at(p[0], p[1], p[2]), Some(sigma));
        }
        assert!(saw_curvature);
        let mut rebound = original.clone();
        rebound
            .rebind_surface_tensions(field(&original, |_| 0.002), &gate)
            .unwrap();
        assert!(rebound.local_sigma.is_none());
        let mut ordinary = original;
        ordinary.step_with_gate(&gate).unwrap();
        rebound.step_with_gate(&gate).unwrap();
        assert_eq!(ordinary.grid().state_bits(), rebound.grid().state_bits());
        assert_eq!(
            ordinary.ledger_mass().to_bits(),
            rebound.ledger_mass().to_bits()
        );
    }

    #[test]
    fn g0_planar_geometry_and_incomplete_source_refuse_without_mutation() {
        let gate = CancelGate::new_clock_free();
        let surface = FreeSurface3::new(
            8,
            8,
            8,
            0.8,
            [0.0; 3],
            0.002,
            ContactModel3::Neutral,
            |x, _, z| z < 3 + usize::from(x > 3),
        )
        .unwrap();
        let before = surface.grid().state_bits();
        let values = field(&surface, |_| 0.002);
        assert!(
            surface
                .planar_marangoni_forces(
                    PlanarSurface3 {
                        axis: 2,
                        outward_sign: 1,
                        position: 3.5,
                        maximum_departure: 0.1,
                        minimum_fill: 0.2,
                    },
                    &values,
                    &gate
                )
                .is_err()
        );
        assert_eq!(before, surface.grid().state_bits());
        let mut candidate = surface.clone();
        let mut missing = values;
        missing.pop_first();
        assert!(candidate.rebind_surface_tensions(missing, &gate).is_err());
        assert!(candidate.local_sigma.is_none());
        assert_eq!(before, candidate.grid().state_bits());
    }

    #[test]
    fn g0_local_capillary_timestep_refuses_a_late_source_change_before_rebinding() {
        let gate = CancelGate::new_clock_free();
        let mut surface = FreeSurface3::new(8, 8, 8, 0.8, [0.0; 3], 0.01,
            ContactModel3::Neutral, |_, _, z| z < 4).unwrap();
        let accepted = field(&surface, |_| 0.02);
        surface.rebind_surface_tensions(accepted.clone(), &gate).unwrap();
        let before = surface.grid().state_bits();
        let mass = surface.ledger_mass().to_bits();
        let mut changed = accepted.clone();
        *changed.last_entry().unwrap().get_mut() = 1.01 * crate::EXPLICIT_CAPILLARY_LATTICE_LIMIT;
        assert!(matches!(surface.rebind_surface_tensions(changed, &gate),
            Err(FreeSurfaceError3::Transport(reason)) if reason.contains("capillary-wave")));
        assert_eq!(surface.local_sigma.as_ref(), Some(&accepted));
        assert_eq!(surface.grid().state_bits(), before);
        assert_eq!(surface.ledger_mass().to_bits(), mass);
    }

    fn graph_surface(axis: usize, sign: i8, height: impl Fn(f64, f64) -> f64) -> FreeSurface3 {
        let mut surface = FreeSurface3::new(
            8,
            8,
            8,
            0.8,
            [0.0; 3],
            0.01,
            ContactModel3::Neutral,
            |_, _, z| z < 4,
        )
        .unwrap();
        let tangent = match axis {
            0 => [1, 2],
            1 => [0, 2],
            _ => [0, 1],
        };
        // Independent subcell volume integration. Summing a complete column
        // recovers the mean graph height, even across several partial cells.
        let fill = |p: [usize; 3]| {
            let mut sum = 0.0;
            for a in 0..8 {
                for b in 0..8 {
                    let u = p[tangent[0]] as f64 + (a as f64 + 0.5) / 8.0;
                    let v = p[tangent[1]] as f64 + (b as f64 + 0.5) / 8.0;
                    let below = (height(u, v) - p[axis] as f64).clamp(0.0, 1.0);
                    sum += if sign > 0 { below } else { 1.0 - below };
                }
            }
            sum / 64.0
        };
        for (&key, tile) in &mut surface.tiles {
            for lane in 0..TILE_CELLS {
                let p =
                    FreeSurface3::transport_coordinates(MassLocation3::Cell(key, lane)).unwrap();
                let f = fill(p);
                let mut next = p;
                let next_fill = if sign > 0 && p[axis] + 1 < 8 {
                    next[axis] += 1;
                    fill(next)
                } else if sign < 0 && p[axis] > 0 {
                    next[axis] -= 1;
                    fill(next)
                } else {
                    0.0
                };
                tile.mass[lane] = f;
                tile.cells[lane] = if f == 0.0 {
                    Cell3::Gas
                } else if f < 1.0 || next_fill == 0.0 {
                    Cell3::Interface
                } else {
                    Cell3::Fluid
                };
            }
        }
        surface
    }

    fn column_impulses(
        surface: &FreeSurface3,
        axis: usize,
        forces: &BTreeMap<(u64, usize), [f64; 3]>,
    ) -> BTreeMap<[usize; 2], [f64; 3]> {
        let tangent = match axis {
            0 => [1, 2],
            1 => [0, 2],
            _ => [0, 1],
        };
        let mut columns = BTreeMap::new();
        for (&(key, lane), force) in forces {
            let p = FreeSurface3::transport_coordinates(MassLocation3::Cell(key, lane)).unwrap();
            let total = columns
                .entry([p[tangent[0]], p[tangent[1]]])
                .or_insert([0.0; 3]);
            for d in 0..3 {
                total[d] += surface.fill(key, lane) * force[d];
            }
        }
        columns
    }

    #[test]
    fn g1_height_graph_traction_is_tangent_with_correct_area_in_every_orientation() {
        let gate = CancelGate::new_clock_free();
        for axis in 0..3 {
            let tangent = match axis {
                0 => [1, 2],
                1 => [0, 2],
                _ => [0, 1],
            };
            for sign in [-1, 1] {
                for (hu, hv) in [(0.0, 0.0), (0.2, -0.3), (-0.4, 0.1)] {
                    let surface =
                        graph_surface(axis, sign, |u, v| 3.7 + hu * (u - 4.0) + hv * (v - 4.0));
                    let mut sigma = field(&surface, |p| {
                        0.01 + 0.0002 * (p[tangent[0]] as f64 + 0.5)
                            - 0.0001 * (p[tangent[1]] as f64 + 0.5)
                    });
                    let forces = surface
                        .height_graph_marangoni_forces(
                            HeightGraphSurface3 {
                                axis,
                                outward_sign: sign,
                                maximum_slope: 0.6,
                                minimum_fill: 0.01,
                            },
                            &mut sigma,
                            &gate,
                        )
                        .unwrap();
                    // Independent projection of the ambient gradient onto
                    // the exact plane, then multiply by actual surface area.
                    let n = [-hu, -hv, 1.0];
                    let norm2 = n.iter().map(|x| x * x).sum::<f64>();
                    let dot = -hu * 0.0002 + hv * 0.0001;
                    let ambient = [0.0002, -0.0001, 0.0];
                    let expected: [f64; 3] =
                        core::array::from_fn(|d| (ambient[d] - dot * n[d] / norm2) * norm2.sqrt());
                    let totals = column_impulses(&surface, axis, &forces);
                    assert_eq!(totals.len(), 64);
                    for total in totals.values() {
                        assert!((total[tangent[0]] - expected[0]).abs() < 3e-15);
                        assert!((total[tangent[1]] - expected[1]).abs() < 3e-15);
                        assert!((total[axis] - expected[2]).abs() < 3e-15);
                        assert!(
                            (total[axis] - hu * total[tangent[0]] - hv * total[tangent[1]]).abs()
                                < 3e-15
                        );
                    }
                    assert!(forces.keys().all(|&(key,lane)| surface.tiles[&key].cells[lane] == Cell3::Interface));
                }
            }
        }
    }

    #[test]
    fn g1_height_graph_curved_mixed_metric_matches_manufactured_surface_gradient() {
        let gate = CancelGate::new_clock_free();
        let surface = graph_surface(2, 1, |u, v| {
            let (u, v) = (u - 4.0, v - 4.0);
            3.5 + 0.03 * u * u - 0.02 * v * v + 0.01 * u * v
        });
        let mut sigma = field(&surface, |p| {
            let (u, v) = (p[0] as f64 - 3.5, p[1] as f64 - 3.5);
            0.01 + 0.00002 * u * u + 0.00001 * v
        });
        let forces = surface
            .height_graph_marangoni_forces(
                HeightGraphSurface3 {
                    axis: 2,
                    outward_sign: 1,
                    maximum_slope: 0.5,
                    minimum_fill: 0.01,
                },
                &mut sigma,
                &gate,
            )
            .unwrap();
        let totals = column_impulses(&surface, 2, &forces);
        assert_eq!(totals.len(), 64);
        for (p, total) in totals {
            let (u, v) = (p[0] as f64 - 3.5, p[1] as f64 - 3.5);
            let (hu, hv) = (0.06 * u + 0.01 * v, -0.04 * v + 0.01 * u);
            let normal = [-hu, -hv, 1.0];
            let square = normal.iter().map(|x| x * x).sum::<f64>();
            let ambient = [0.00004 * u, 0.00001, 0.0];
            let dot = ambient.iter().zip(normal).map(|(a, n)| a * n).sum::<f64>();
            for d in 0..3 {
                let expected = (ambient[d] - normal[d] * dot / square) * square.sqrt();
                assert!((total[d] - expected).abs() < 5e-15, "{p:?}: {total:?}");
            }
        }
    }

    #[test]
    fn g0_height_graph_constant_sigma_and_invalid_geometry_do_not_invent_stress() {
        let gate = CancelGate::new_clock_free();
        let surface = graph_surface(2, 1, |u, v| 3.7 + 0.2 * (u - 4.0) - 0.3 * (v - 4.0));
        let before = surface.grid().state_bits();
        let mut sigma = field(&surface, |_| 0.01);
        let original = sigma.clone();
        let graph = HeightGraphSurface3 {
            axis: 2,
            outward_sign: 1,
            maximum_slope: 0.6,
            minimum_fill: 0.01,
        };
        assert!(
            surface
                .height_graph_marangoni_forces(graph, &mut sigma, &gate)
                .unwrap()
                .is_empty()
        );
        assert_eq!(sigma, original);
        assert!(
            surface
                .height_graph_marangoni_forces(
                    HeightGraphSurface3 {
                        maximum_slope: 0.1,
                        ..graph
                    },
                    &mut sigma,
                    &gate,
                )
                .is_err()
        );
        assert_eq!(sigma, original);
        assert_eq!(before, surface.grid().state_bits());
        let mut broken = surface.clone();
        let (slot, lane) = broken.grid.resolve_source(3, 3, 1).unwrap();
        let key = broken.grid.active_keys()[slot];
        broken.tiles.get_mut(&key).unwrap().cells[lane] = Cell3::Gas;
        broken.tiles.get_mut(&key).unwrap().mass[lane] = 0.0;
        let mut broken_sigma = field(&broken, |_| 0.01);
        assert!(
            broken
                .height_graph_marangoni_forces(graph, &mut broken_sigma, &gate)
                .is_err()
        );
    }

    #[test]
    fn g1_height_graph_pressure_and_traction_share_volume_weighted_interface_sigma() {
        let gate = CancelGate::new_clock_free();
        let surface = graph_surface(2, 1, |u, v| 3.7 + 0.2 * (u - 4.0) - 0.3 * (v - 4.0));
        let mut sigma = field(&surface, |p| 0.01 + 0.0002 * p[2] as f64);
        let original = sigma.clone();
        let mut expected = BTreeMap::new();
        for (&(key, lane), &value) in &original {
            if surface.tiles[&key].cells[lane] != Cell3::Interface {
                continue;
            }
            let p = FreeSurface3::transport_coordinates(MassLocation3::Cell(key, lane)).unwrap();
            let entry = expected.entry([p[0], p[1]]).or_insert([0.0; 2]);
            let fill = surface.fill(key, lane);
            entry[0] += fill * value;
            entry[1] += fill;
        }
        surface
            .height_graph_marangoni_forces(
                HeightGraphSurface3 {
                    axis: 2,
                    outward_sign: 1,
                    maximum_slope: 0.6,
                    minimum_fill: 0.01,
                },
                &mut sigma,
                &gate,
            )
            .unwrap();
        for (&(key, lane), &value) in &sigma {
            if surface.tiles[&key].cells[lane] == Cell3::Interface {
                let p =
                    FreeSurface3::transport_coordinates(MassLocation3::Cell(key, lane)).unwrap();
                let [weighted, volume] = expected[&[p[0], p[1]]];
                assert!((value - weighted / volume).abs() < 1e-17);
            } else {
                assert_eq!(value, original[&(key, lane)]);
            }
        }
    }
}
