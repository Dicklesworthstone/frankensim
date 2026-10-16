//! The global-buckling spot check (§15.2 step 2): the critical
//! compression member, re-analyzed as an fs-solid geometrically exact
//! Cosserat rod with a seeded mid-span imperfection, loaded axially to
//! a factor beyond design. If the rod statics converge with bounded
//! transverse growth, the member is stable through the check factor;
//! a stall or runaway deflection is the buckling refusal the LP-level
//! Euler floor cannot see (it assumed perfect pin-ended members).

use fs_solid::{Rod, RodSection, TipLoad};

/// Re-analyze one compression member: returns (stable, max transverse
/// deflection / length) at `factor` × the design force.
///
/// The member is modeled clamped-free with a `imperfection`·length
/// mid-span bow (conservative vs pinned-pinned for the same Euler
/// load convention used in sizing).
#[must_use]
pub fn rod_buckling_check(
    length: f64,
    area: f64,
    youngs: f64,
    design_force: f64,
    factor: f64,
    imperfection: f64,
) -> (bool, f64) {
    // Solid square: I = A²/12; shear/torsion stiff.
    let inertia = area * area / 12.0;
    let section = RodSection {
        ea: youngs * area,
        ga: youngs * area, // shear-stiff screening model
        gj: youngs * inertia,
        ei: youngs * inertia,
    };
    let mut rod = Rod::straight(length, 10, section);
    // Seed the imperfection: a small transverse bow.
    let n = rod.positions.len();
    for (i, p) in rod.positions.iter_mut().enumerate() {
        #[allow(clippy::cast_precision_loss)]
        let s = i as f64 / (n - 1) as f64;
        p[1] += imperfection * length * (std::f64::consts::PI * s).sin();
    }
    let load = TipLoad {
        force: [-design_force.abs() * factor, 0.0, 0.0],
        moment: [0.0; 3],
    };
    match rod.solve_static(&load, 10, 1e-8) {
        Ok(_) => {
            let max_w = rod
                .positions
                .iter()
                .map(|p| p[1].abs().max(p[2].abs()))
                .fold(0.0f64, f64::max)
                / length;
            // Stable if the bow stays bounded (an order above the
            // seeded imperfection is the runaway line).
            (max_w < 10.0 * imperfection + 0.05, max_w)
        }
        Err(_) => (false, f64::INFINITY),
    }
}
