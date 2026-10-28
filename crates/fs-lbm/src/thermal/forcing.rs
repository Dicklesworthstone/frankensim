//! Shared admission and Guo source frames for held thermal/physical forces.

// Guo, Zheng & Shi, Phys. Rev. E 65, 046308 (2002),
// doi:10.1103/PhysRevE.65.046308. Stored populations are f_phys - S(u,F)/2.
// Replacing F without replacing that half source adds an unphysical
// (F_next-F_previous)/2 impulse before collision. Reframe at fixed physical u
// so both momentum and Pi_neq + (u F + F u)/2 remain unchanged. This is a
// coordinate change, not elapsed time or a heat/work contribution.
fn reframe<const Q: usize>(populations: &mut [f64; Q], source: impl Fn(usize) -> f64) {
    let rho = populations.iter().sum::<f64>();
    for (q, value) in populations.iter_mut().enumerate() {
        *value += source(q);
    }
    // The analytic source has zero zeroth moment. Keep summation roundoff in
    // the rest population, whose lattice velocity is zero for both stencils.
    populations[0] += rho - populations.iter().sum::<f64>();
}

pub(crate) fn reframe2(
    populations: &mut [f64; crate::Q],
    velocity: [f64; 2],
    previous: [f64; 2],
    next: [f64; 2],
) {
    let half_delta = [0.5 * previous[0] - 0.5 * next[0], 0.5 * previous[1] - 0.5 * next[1]];
    reframe(populations, |q| {
        let (ex, ey) = (f64::from(crate::E[q].0), f64::from(crate::E[q].1));
        let eu = ex * velocity[0] + ey * velocity[1];
        let fx = (ex - velocity[0]) / crate::CS2 + eu * ex / (crate::CS2 * crate::CS2);
        let fy = (ey - velocity[1]) / crate::CS2 + eu * ey / (crate::CS2 * crate::CS2);
        crate::W[q] * (fx * half_delta[0] + fy * half_delta[1])
    });
}

pub(crate) fn reframe3(
    populations: &mut [f64; crate::d3q19::Q3],
    velocity: [f64; 3],
    previous: [f64; 3],
    next: [f64; 3],
) {
    use crate::d3q19::{E3, W3, guo_projection3};
    let half_delta = core::array::from_fn(|axis| 0.5 * previous[axis] - 0.5 * next[axis]);
    reframe(populations, |q| {
        let e = [f64::from(E3[q].0), f64::from(E3[q].1), f64::from(E3[q].2)];
        let eu = e.iter().zip(velocity).map(|(e, u)| e * u).sum();
        W3[q] * guo_projection3(e, velocity, eu, half_delta)
    });
}

/// A force field that cannot be applied to the existing thermal solver.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ThermalForcingError {
    /// The field must contain one vector per cell, including wall cells.
    Length {
        /// Required number of cell vectors.
        expected: usize,
        /// Supplied number of cell vectors.
        actual: usize,
    },
    /// An external component or its sum with buoyancy is not finite.
    NonFinite {
        /// Linear grid-cell index.
        cell: usize,
        /// Zero-based force-vector component.
        component: usize,
    },
    /// The declared wall plane cannot describe the temperature storage.
    Layout,
}
impl core::fmt::Display for ThermalForcingError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Length { expected, actual } => write!(f, "thermal force field has {actual} cells, expected {expected}"),
            Self::NonFinite { cell, component } => write!(f, "nonfinite thermal force at cell {cell}, component {component}"),
            Self::Layout => write!(f, "invalid thermal force/temperature layout"),
        }
    }
}
impl std::error::Error for ThermalForcingError {}

pub(crate) fn validate_external_forces<const D: usize, const Q: usize>(
    forces: &[[f64; D]],
    temperature_populations: &[[f64; Q]],
    wall_plane: usize,
    gbeta: f64,
    reference_temperature: f64,
) -> Result<(), ThermalForcingError> {
    let n = temperature_populations.len();
    if forces.len() != n {
        return Err(ThermalForcingError::Length { expected: n, actual: forces.len() });
    }
    if D == 0 || wall_plane == 0 || wall_plane > n / 2 {
        return Err(ThermalForcingError::Layout);
    }
    for (cell, force) in forces.iter().enumerate() {
        for (component, value) in force.iter().enumerate() {
            if !value.is_finite() {
                return Err(ThermalForcingError::NonFinite { cell, component });
            }
        }
        if cell >= wall_plane && cell < n - wall_plane {
            let temperature: f64 = temperature_populations[cell].iter().sum();
            let total = gbeta * (temperature - reference_temperature) + force[D - 1];
            if !total.is_finite() {
                return Err(ThermalForcingError::NonFinite { cell, component: D - 1 });
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn physical_moments<const D: usize, const Q: usize>(
        populations: &[f64; Q],
        directions: [[f64; D]; Q],
        force: [f64; D],
    ) -> (f64, [f64; D], [[f64; D]; D]) {
        let rho = populations.iter().sum::<f64>();
        let velocity = core::array::from_fn(|a| {
            (populations.iter().zip(directions).map(|(f, e)| f * e[a]).sum::<f64>()
                + 0.5 * force[a]) / rho
        });
        let second = core::array::from_fn(|a| core::array::from_fn(|b| {
            populations.iter().zip(directions).map(|(f, e)| f * e[a] * e[b]).sum::<f64>()
                + 0.5 * (velocity[a] * force[b] + force[a] * velocity[b])
        }));
        (rho, velocity, second)
    }

    fn compare<const D: usize>(
        before: (f64, [f64; D], [[f64; D]; D]),
        after: (f64, [f64; D], [[f64; D]; D]),
    ) {
        assert!((before.0 - after.0).abs() < 2e-15);
        for a in 0..D {
            assert!((before.1[a] - after.1[a]).abs() < 2e-15);
            for b in 0..D {
                assert!((before.2[a][b] - after.2[a][b]).abs() < 2e-15);
            }
        }
    }

    #[test]
    fn g0_changing_force_frames_preserves_mass_velocity_and_full_stress() {
        // Include a nonequilibrium stress, so replacing the state by an
        // equilibrium distribution would fail even though velocity survives.
        let mut f = crate::equilibrium(1.7, 0.02, -0.03);
        f[1] += 0.002;
        f[3] += 0.002;
        f[0] -= 0.004;
        let e = crate::E.map(|(x, y)| [f64::from(x), f64::from(y)]);
        let initial = physical_moments(&f, e, [0.0; 2]);
        let shear = crate::core2::force_corrected_shear_rate(&f, initial.0, initial.1, [0.0; 2], 0.9);
        let mut previous = [0.0; 2];
        for next in [[0.01, -0.02], [-0.03, 0.015], [0.0; 2]] {
            let velocity = physical_moments(&f, e, previous).1;
            reframe2(&mut f, velocity, previous, next);
            let actual = physical_moments(&f, e, next);
            compare(initial, actual);
            let rate = crate::core2::force_corrected_shear_rate(&f, actual.0, actual.1, next, 0.9);
            assert!((rate - shear).abs() < 2e-15);
            previous = next;
        }

        let mut f = crate::equilibrium3(1.3, [0.01, -0.02, 0.03]);
        f[7] += 0.002;
        f[8] += 0.002;
        f[0] -= 0.004;
        let e = crate::d3q19::E3.map(|(x, y, z)| [f64::from(x), f64::from(y), f64::from(z)]);
        let initial = physical_moments(&f, e, [0.0; 3]);
        let shear = crate::d3q19::force_corrected_shear_rate3(&f, initial.0, initial.1, [0.0; 3], 0.9);
        let mut previous = [0.0; 3];
        for next in [[0.01, -0.02, 0.03], [-0.03, 0.015, -0.005], [0.0; 3]] {
            let velocity = physical_moments(&f, e, previous).1;
            reframe3(&mut f, velocity, previous, next);
            let actual = physical_moments(&f, e, next);
            compare(initial, actual);
            let rate = crate::d3q19::force_corrected_shear_rate3(&f, actual.0, actual.1, next, 0.9);
            assert!((rate - shear).abs() < 2e-15);
            previous = next;
        }
    }
}
