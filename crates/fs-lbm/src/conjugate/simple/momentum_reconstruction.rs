//! Conservative deferred momentum reconstruction on the original MAC faces.
//!
//! The implicit operator remains donor-cell upwind. A shared physical face
//! carries one limited correction, with opposite signs in its two rows. No
//! reconstruction crosses a prescribed velocity, pressure opening or solid.

use super::{Kind, Solver};

/// Reconstruction of advected velocity on staggered momentum control volumes.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum MomentumReconstruction {
    /// Original convection operator without a deferred correction.
    #[default]
    DonorCell,
    /// Minmod-limited linear velocity with physical nonuniform-grid gradients.
    /// Requires the implicit upwind convection scheme. Smooth interior affine
    /// fields are reproduced exactly; extrema and unsupported boundary
    /// stencils fall back to donor-cell transport. This alone is not a global
    /// TVD or kinetic-energy guarantee for the coupled multidimensional solve.
    Minmod,
}

impl Solver<'_> {
    /// Outward deferred flux `F * (u_face - u_upwind)` for one original side.
    /// `neighbor` is adjacent along `d`, with `plus` declaring its orientation.
    pub(super) fn momentum_reconstruction_flux(
        &self,
        a: usize,
        d: usize,
        index: usize,
        neighbor: usize,
        plus: bool,
        flux: f64,
    ) -> f64 {
        if self.config.momentum_reconstruction != MomentumReconstruction::Minmod || flux == 0.0 {
            return 0.0;
        }
        // The adjacent row must use precisely the same donor, differences,
        // distances and product. Only the final outward orientation changes.
        let (lower, upper) = if plus {
            (index, neighbor)
        } else {
            (neighbor, index)
        };
        let forward_flux = if plus { flux } else { -flux };
        let donor = if forward_flux > 0.0 { lower } else { upper };
        let comp = &self.comps[a];
        let at = comp.coords(donor);
        let p = at[d];
        if p == 0 || p + 1 >= comp.dims[d] {
            return 0.0;
        }
        let mut left = at;
        left[d] -= 1;
        let mut right = at;
        right[d] += 1;
        let left = comp.index(left);
        let right = comp.index(right);
        if ![left, donor, right]
            .into_iter()
            .all(|i| matches!(comp.kind[i], Kind::Unknown(_)))
        {
            return 0.0;
        }
        let width = |i| self.domain.width(d, i);
        let (left_distance, right_distance, offset) = if d == a {
            // Component-normal nodes are cell edges; the shared momentum side
            // is the center of the cell between the two velocity nodes.
            (
                width(p - 1),
                width(p),
                if forward_flux > 0.0 {
                    0.5 * width(p)
                } else {
                    -0.5 * width(p - 1)
                },
            )
        } else {
            // Transverse nodes are cell centers; their shared side is the
            // original cell edge, half the DONOR width from its center.
            (
                0.5 * (width(p - 1) + width(p)),
                0.5 * (width(p) + width(p + 1)),
                if forward_flux > 0.0 {
                    0.5 * width(p)
                } else {
                    -0.5 * width(p)
                },
            )
        };
        let values = &self.vel[a];
        let left_slope = (values[donor] - values[left]) / left_distance;
        let right_slope = (values[right] - values[donor]) / right_distance;
        let slope = if left_slope > 0.0 && right_slope > 0.0 {
            left_slope.min(right_slope)
        } else if left_slope < 0.0 && right_slope < 0.0 {
            left_slope.max(right_slope)
        } else {
            0.0
        };
        // Reconstruct the difference directly: adding/subtracting an absolute
        // face velocity would needlessly lose small corrections under a boost.
        let correction = forward_flux * (offset * slope);
        if plus { correction } else { -correction }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::conjugate::simple::{ConvectionScheme, FvBoundary, SimpleConfig};
    use crate::conjugate::{FluidProperties, Voxel, VoxelDomain};
    use fs_exec::CancelGate;

    fn fluid() -> FluidProperties {
        FluidProperties {
            density_kg_m3: 1.0,
            specific_heat_j_kg_k: 2.0,
            conductivity_w_m_k: 0.01,
            kinematic_viscosity_m2_s: 0.001,
        }
    }

    fn config(mode: MomentumReconstruction) -> SimpleConfig {
        let mut config = SimpleConfig::new([FvBoundary::Symmetry; 6]);
        config.scheme = ConvectionScheme::Upwind;
        config.momentum_reconstruction = mode;
        config
    }

    fn graded() -> VoxelDomain {
        VoxelDomain::graded([
            vec![0.25, 0.5, 1.0, 2.0, 0.75, 1.5, 0.375, 0.625],
            vec![1.0, 0.25, 1.5, 0.5, 2.0, 0.75, 0.375, 0.625],
            vec![0.5, 1.5, 0.25, 1.0, 0.75, 2.0, 0.625, 0.375],
        ])
        .unwrap()
    }

    fn coordinate(solver: &Solver<'_>, a: usize, d: usize, index: usize) -> f64 {
        let p = solver.comps[a].coords(index)[d];
        solver.domain.face_coord(d, p)
            + if d == a {
                0.0
            } else {
                0.5 * solver.domain.width(d, p)
            }
    }

    fn affine(solver: &mut Solver<'_>, a: usize, d: usize, slope: f64, datum: f64) {
        for i in 0..solver.vel[a].len() {
            solver.vel[a][i] = datum + slope * coordinate(solver, a, d, i);
        }
    }

    fn pair(solver: &Solver<'_>, a: usize, d: usize, at: [usize; 3]) -> (usize, usize) {
        let mut next = at;
        next[d] += 1;
        (solver.comps[a].index(at), solver.comps[a].index(next))
    }

    #[test]
    fn limited_momentum_requires_its_original_implicit_upwind_operator() {
        let domain = graded();
        let mut config = config(MomentumReconstruction::Minmod);
        for scheme in [ConvectionScheme::PowerLaw, ConvectionScheme::Exponential] {
            config.scheme = scheme;
            assert!(matches!(
                super::super::admit(&domain, &fluid(), &config),
                Err(crate::conjugate::ChtError::InvalidInput {
                    field: "simple.momentum_reconstruction",
                    ..
                })
            ));
        }
        config.scheme = ConvectionScheme::Upwind;
        assert!(super::super::admit(&domain, &fluid(), &config).is_ok());
        config.scheme = ConvectionScheme::EnthalpyMinmod;
        assert!(matches!(
            super::super::admit(&domain, &fluid(), &config),
            Err(crate::conjugate::ChtError::InvalidInput {
                field: "simple.scheme",
                ..
            })
        ));
    }

    #[test]
    fn graded_affine_faces_recover_both_mac_geometries_and_cancel_between_rows() {
        let domain = graded();
        let config = config(MomentumReconstruction::Minmod);
        let mut solver = Solver::new(&domain, &fluid(), &config);
        for a in 0..3 {
            for d in 0..3 {
                let at = [3; 3];
                let (lower, upper) = pair(&solver, a, d, at);
                let face = if d == a {
                    0.5 * (domain.face_coord(d, at[d]) + domain.face_coord(d, at[d] + 1))
                } else {
                    domain.face_coord(d, at[d] + 1)
                };
                for slope in [-1.5, 1.5] {
                    affine(&mut solver, a, d, slope, 8.0);
                    for flux in [-2.0, 2.0] {
                        let correction =
                            solver.momentum_reconstruction_flux(a, d, lower, upper, true, flux);
                        let reverse =
                            solver.momentum_reconstruction_flux(a, d, upper, lower, false, -flux);
                        assert_eq!(reverse.to_bits(), (-correction).to_bits());
                        let donor = if flux > 0.0 { lower } else { upper };
                        let reconstructed = solver.vel[a][donor] + correction / flux;
                        assert!(
                            (reconstructed - (8.0 + slope * face)).abs() < 2e-14,
                            "component {a}, direction {d}, flux {flux}"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn limiter_preserves_extrema_and_corrections_ignore_a_constant_velocity_datum() {
        let domain = graded();
        let config = config(MomentumReconstruction::Minmod);
        let mut solver = Solver::new(&domain, &fluid(), &config);
        for a in 0..3 {
            for d in 0..3 {
                let (lower, upper) = pair(&solver, a, d, [3; 3]);
                for flux in [-2.0, 2.0] {
                    affine(&mut solver, a, d, 1.5, 8.0);
                    let base = solver.momentum_reconstruction_flux(a, d, lower, upper, true, flux);
                    affine(&mut solver, a, d, 1.5, 8.0 + 1_048_576.0);
                    let boosted =
                        solver.momentum_reconstruction_flux(a, d, lower, upper, true, flux);
                    assert_eq!(base.to_bits(), boosted.to_bits());
                    let donor = if flux > 0.0 { lower } else { upper };
                    let mut left = solver.comps[a].coords(donor);
                    left[d] -= 1;
                    let mut right = solver.comps[a].coords(donor);
                    right[d] += 1;
                    let left = solver.comps[a].index(left);
                    let right = solver.comps[a].index(right);
                    for sign in [-1.0, 1.0] {
                        solver.vel[a][left] = sign;
                        solver.vel[a][donor] = 3.0 * sign;
                        solver.vel[a][right] = 2.0 * sign;
                        assert_eq!(
                            solver.momentum_reconstruction_flux(a, d, lower, upper, true, flux),
                            0.0
                        );
                    }
                    affine(&mut solver, a, d, 0.0, 8.0);
                    assert_eq!(
                        solver.momentum_reconstruction_flux(a, d, lower, upper, true, flux),
                        0.0
                    );
                }
            }
        }
    }

    #[test]
    fn original_boundaries_solids_and_disabled_reconstruction_keep_donor_transport() {
        for boundary in [FvBoundary::wall(), FvBoundary::Outlet] {
            let domain = graded();
            let mut config = config(MomentumReconstruction::Minmod);
            config.faces[0] = boundary;
            let mut solver = Solver::new(&domain, &fluid(), &config);
            affine(&mut solver, 0, 0, 1.5, 8.0);
            let (lower, upper) = pair(&solver, 0, 0, [1, 3, 3]);
            assert_eq!(
                solver.momentum_reconstruction_flux(0, 0, lower, upper, true, 2.0),
                0.0
            );
        }
        let mut domain = graded();
        domain.set(1, 3, 3, Voxel::Solid(0));
        let limited = config(MomentumReconstruction::Minmod);
        let mut solver = Solver::new(&domain, &fluid(), &limited);
        affine(&mut solver, 0, 0, 1.5, 8.0);
        let (lower, upper) = pair(&solver, 0, 0, [3; 3]);
        assert!(matches!(solver.comps[0].kind[lower], Kind::Unknown(_)));
        assert!(matches!(solver.comps[0].kind[upper], Kind::Unknown(_)));
        assert_eq!(
            solver.momentum_reconstruction_flux(0, 0, lower, upper, true, 2.0),
            0.0
        );
        // No transverse neighbor is invented beyond the original face lattice.
        affine(&mut solver, 1, 2, 1.5, 8.0);
        let (lower, upper) = pair(&solver, 1, 2, [3, 3, 0]);
        assert_eq!(
            solver.momentum_reconstruction_flux(1, 2, lower, upper, true, 2.0),
            0.0
        );

        let domain = graded();
        let original = config(MomentumReconstruction::DonorCell);
        let mut solver = Solver::new(&domain, &fluid(), &original);
        affine(&mut solver, 0, 1, 1.5, 8.0);
        let (lower, upper) = pair(&solver, 0, 1, [3; 3]);
        assert_eq!(
            solver.momentum_reconstruction_flux(0, 1, lower, upper, true, 2.0),
            0.0
        );
        assert_eq!(
            solver.momentum_reconstruction_flux(0, 1, lower, upper, true, 0.0),
            0.0
        );
    }

    #[test]
    fn native_momentum_residual_improves_for_smooth_transported_shear_on_each_axis() {
        let gate = CancelGate::new();
        let fluid = fluid();
        let g = |y: f64| y * y * (1.0 - y) * (1.0 - y);
        for (a, d) in [(0, 1), (1, 2), (2, 0)] {
            let mut widths = [vec![1.0], vec![1.0], vec![1.0]];
            widths[a] = vec![1.0 / 3.0; 3];
            widths[d] = vec![1.0 / 32.0; 32];
            let domain = VoxelDomain::graded(widths).unwrap();
            let residuals = [
                MomentumReconstruction::DonorCell,
                MomentumReconstruction::Minmod,
            ]
            .map(|mode| {
                let mut config = config(mode);
                config.faces[2 * a] = FvBoundary::Outlet;
                config.faces[2 * a + 1] = FvBoundary::Outlet;
                let mut transport = [0.0; 3];
                transport[d] = 1.0;
                config.faces[2 * d] = FvBoundary::Inlet {
                    velocity: transport,
                };
                config.faces[2 * d + 1] = FvBoundary::Outlet;
                let mut solver = Solver::new(&domain, &fluid, &config);
                for component in 0..3 {
                    for i in 0..solver.vel[component].len() {
                        let y = coordinate(&solver, component, d, i);
                        solver.vel[component][i] = if component == a {
                            g(y)
                        } else {
                            transport[component]
                        };
                    }
                }
                // Exact steady NS forcing for u_a=g(y), u_d=1, p=0:
                // f_a = rho*g'(y) - mu*g''(y). Boundary g=g'=0 removes an
                // imposed tangential inlet/outlet mismatch. No solve or
                // alternate reconstruction is used to measure this residual.
                solver.set_force(
                    (0..domain.cell_count())
                        .map(|cell| {
                            let at = domain.coords(cell);
                            let y = domain.center(at[0], at[1], at[2])[d];
                            let mut force = [0.0; 3];
                            force[a] = 2.0 * y - 6.0 * y * y + 4.0 * y * y * y
                                - fluid.kinematic_viscosity_m2_s * (2.0 - 12.0 * y + 12.0 * y * y);
                            force
                        })
                        .collect(),
                );
                solver.momentum(a, false, &gate).unwrap().rms
            });
            assert!(residuals[0] > 0.0);
            assert!(
                residuals[1] < 0.9 * residuals[0],
                "component {a}, direction {d}: {residuals:?}"
            );
        }
    }
}
