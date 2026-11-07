//! Compatible molecular mechanics on periodic, zero-flux and free-slip grids.
//!
//! With cell viscosities, centered derivatives C and half-differences H of
//! neighboring face derivatives, the negative kinetic-energy rate is
//! `2 mu |dev sym C(u)|^2 + bulk div(Cu)^2 + mu |H(u)|^2
//!  + (bulk + mu/3) sum_j H_j(u_j)^2`, summed over cells.
//! In two dimensions the deviatoric trace retains the physical third dimension.
//! Every term is nonnegative. Expanding the adjoints gives the same compact
//! direct-normal/centered-cross stencil as constant-coefficient Newtonian stress,
//! but each cross derivative must carry its own cell's coefficient BEFORE the
//! two cell contributions are averaged onto the face.
//!
//! Prescribed velocity, finite-slip and bath boundaries need their own discrete
//! boundary-work derivation. Callers retain their original full operator there.
//! This is a semidiscrete molecular identity, not a time-step or SGS theorem.

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

#[cfg(test)]
mod tests;
