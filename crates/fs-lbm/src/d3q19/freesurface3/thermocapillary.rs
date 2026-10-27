//! Checked flat-surface traction on the original VOF/Guo owner.
//!
//! The interface is a single-valued, nearly planar tank surface, perpendicular
//! to one Cartesian axis. Each surface column has area one and occupied volume
//! `fill` in lattice units, so its force density is `grad_t(sigma)/fill`.
//! This is a declared flat-interface stress approximation, not a diffuse delta
//! function or a reconstructed arbitrary curved surface. The normal Laplace
//! load remains exclusively in the original pressure boundary.

use super::*;

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
}
