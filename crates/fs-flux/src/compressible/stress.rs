//! Compatible molecular mechanics on Cartesian grids and prescribed walls.
//!
//! With cell viscosities, centered derivatives C and half-differences H of
//! neighboring face derivatives, wall mechanical power minus kinetic-energy rate is
//! `2 mu |dev sym C(u)|^2 + bulk div(Cu)^2 + mu |H(u)|^2
//!  + (bulk + mu/3) sum_j H_j(u_j)^2`, summed over cells, plus physical
//! Navier interface friction heat, summed over walls.
//! In two dimensions the deviatoric trace retains the physical third dimension.
//! Every term is nonnegative. Expanding the adjoints gives the same compact
//! direct-normal/centered-cross stencil as constant-coefficient Newtonian stress,
//! but each cross derivative must carry its own cell's coefficient BEFORE the
//! two cell contributions are averaged onto the face.
//!
//! A prescribed-velocity wall uses the same odd ghost `2 U_wall - u` as the
//! cell gradient. Its positive-axis traction is `T-q` on the lower side and
//! `T+q` on the upper side, where `T=mu(C+C^T)+(bulk-2mu/3)tr(C)I` and
//! `q_aj=mu H_j(u_a)+(bulk+mu/3)delta_aj H_j(u_j)`. Hence its energy flux
//! uses `-U_wall dot traction`; outward wall power is `U_wall dot traction*n`.
//! This local identity also covers corners and singleton axes, without averaging
//! different planes' prescribed velocities. Zero-length Navier slip is the same
//! Dirichlet case. Finite Navier slip eliminates the fluid trace from this
//! positive quadratic form plus physical interface friction. Natural free slip
//! uses the matching zero-traction elimination, including its cross derivative;
//! it is not an even tangential ghost. Each off-diagonal derivative pair is a
//! local two-by-two solve. Pure-free singleton axes choose zero rotation only
//! when their provisional right-hand sides vanish.
//!
//! Prescribed baths and modeled equilibrium wall laws keep their complete native
//! operator. This is a semidiscrete molecular identity, not a time-step or SGS
//! theorem. It does not extend periodic second-order accuracy to wall traction;
//! finite slip with cell-sampled wall viscosity may converge at first order.

mod slip;
pub(super) use slip::{GradientSample, Slip, velocity_gradient, wall_stress};

/// Source-evaluated cell viscosities [Pa s]; sixteen numerical bytes per cell.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub(super) struct Viscosity {
    pub shear: f64,
    pub bulk: f64,
}

impl Viscosity {
    pub fn valid(self) -> bool {
        [self.shear, self.bulk]
            .into_iter()
            .all(|v| v.is_finite() && v >= 0.0)
    }

    /// The existing full-tensor row-sum scale, before density and cell widths.
    pub fn scale<const D: usize>(self) -> f64 {
        (4.0 / 3.0) * self.shear
            + self.bulk
            + (D - 1) as f64 * (self.bulk - (2.0 / 3.0) * self.shear).abs()
            + (2 * (D - 1)) as f64 * self.shear
    }
}

/// Compatible mechanical traction in the positive coordinate direction.
///
/// The normal is the positive coordinate `axis`; returned components are in
/// global order. `gradient(cell, component, derivative)` reads a cached centered
/// velocity derivative, with cell 0=left and 1=right. Source scalar/enthalpy
/// fluxes remain unchanged. The owner uses `-traction` and its matching work
/// `-u_face dot traction` in its conservative molecular face flux. Evaluate this
/// directly: subtracting an arbitrarily large midpoint-law stress would lose
/// finite cell stress and heat. Callers preserve their native arithmetic when
/// both cell coefficients equal the original face coefficients exactly.
#[inline]
pub(super) fn traction<const D: usize>(
    axis: usize,
    distance: f64,
    velocity: [[f64; D]; 2],
    cell: [Viscosity; 2],
    gradient: impl Fn(usize, usize, usize) -> f64,
) -> [f64; D] {
    let mean_mu = 0.5 * cell[0].shear + 0.5 * cell[1].shear;
    let mean_bulk = 0.5 * cell[0].bulk + 0.5 * cell[1].bulk;
    std::array::from_fn(|component| {
        let normal = (velocity[1][component] - velocity[0][component]) / distance;
        if component == axis {
            let mut value = ((4.0 / 3.0) * mean_mu + mean_bulk) * normal;
            for side in 0..2 {
                let cross = cell[side].bulk - (2.0 / 3.0) * cell[side].shear;
                let divergence: f64 = (0..D)
                    .filter(|&t| t != axis)
                    .map(|t| gradient(side, t, t))
                    .sum();
                value += 0.5 * cross * divergence;
            }
            value
        } else {
            mean_mu * normal
                + 0.5 * cell[0].shear * gradient(0, axis, component)
                + 0.5 * cell[1].shear * gradient(1, axis, component)
        }
    })
}

/// Complete compatible traction at a prescribed physical-velocity wall.
///
/// `outward_sign` is -1 on a lower plane and +1 on an upper plane;
/// `wall_distance` is the positive center-to-wall distance (half a grid cell).
/// Return components are global `tau[a, axis]`, along the POSITIVE coordinate
/// axis, not outward traction. Inputs use the same cell source viscosities and
/// cached gradient as the interior operator; `gradient(component, derivative)`
/// must include the actual wall's odd velocity ghost. Geometry owners admit the
/// normal wall speed against their conservative stage geometry before calling
/// this helper. Stationary owners continue to require zero normal velocity.
///
/// Both cross terms belong to the discrete adjoint, even for constant viscosity.
/// Replacing them by derivatives of a constant boundary trace breaks the energy
/// balance. They approach that physical trace under refinement, but this helper
/// does not claim second-order wall accuracy. Molecular energy flux is
/// `-wall_velocity dot traction` plus the original source's thermal flux;
/// thermal coefficients and any separate SGS contribution keep their own laws.
#[inline]
pub(super) fn no_slip_traction<const D: usize>(
    axis: usize,
    outward_sign: f64,
    wall_distance: f64,
    velocity: [f64; D],
    wall_velocity: [f64; D],
    cell: Viscosity,
    gradient: impl Fn(usize, usize) -> f64,
) -> [f64; D] {
    let longitudinal = (4.0 / 3.0) * cell.shear + cell.bulk;
    let cross = cell.bulk - (2.0 / 3.0) * cell.shear;
    std::array::from_fn(|component| {
        let normal =
            outward_sign * (wall_velocity[component] - velocity[component]) / wall_distance;
        if component == axis {
            let divergence: f64 = (0..D).filter(|&t| t != axis).map(|t| gradient(t, t)).sum();
            longitudinal * normal + cross * divergence
        } else {
            cell.shear * (normal + gradient(axis, component))
        }
    })
}

#[cfg(test)]
mod tests;
