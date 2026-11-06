//! The global-buckling spot check (§15.2 step 2): the critical
//! compression member, re-analyzed as an fs-solid geometrically exact
//! Cosserat rod with a seeded mid-span imperfection, loaded axially to
//! a factor beyond design. Acceptance requires a converged dead-force
//! equilibrium, a positive free-coordinate tangent and bounded transverse
//! growth. Newton can converge to a straight but unstable compressed state;
//! a small residual or bow alone is never a stability test. This checks the
//! endpoint equilibrium, not every load along a continuation path.

use fs_solid::{Rod, RodSection, TipLoad};

/// Re-analyze one compression member: returns (stable, max transverse
/// deflection / length) at `factor` × the design force.
///
/// The member is modeled clamped-free with an initial `imperfection`·length
/// mid-span bow. This perturbs the Newton iterate, not the stress-free reference
/// geometry. The clamped-free Euler limit is one quarter of the pinned-pinned
/// load used in sizing. Invalid inputs or a failed solve return `(false, inf)`;
/// a nonpositive or numerically unresolved tangent returns `false` with the
/// converged bow. Zero imperfection is admitted and cannot bypass stability.
#[must_use]
pub fn rod_buckling_check(
    length: f64,
    area: f64,
    youngs: f64,
    design_force: f64,
    factor: f64,
    imperfection: f64,
) -> (bool, f64) {
    if [length, area, youngs, factor]
        .iter()
        .any(|x| !x.is_finite() || *x <= 0.0)
        || !design_force.is_finite()
        || !imperfection.is_finite()
        || imperfection < 0.0
    {
        return (false, f64::INFINITY);
    }
    // Solid square: I = A²/12; shear/torsion stiff.
    let inertia = area * area / 12.0;
    let section = RodSection {
        ea: youngs * area,
        ga: youngs * area, // shear-stiff screening model
        gj: youngs * inertia,
        ei: youngs * inertia,
    };
    if [section.ea, section.ga, section.gj, section.ei]
        .iter()
        .any(|x| !x.is_finite() || *x <= 0.0)
    {
        return (false, f64::INFINITY);
    }
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
            let stable = rod.is_stable_under_dead_load(&load, 1e-8).unwrap_or(false);
            (stable && max_w < 10.0 * imperfection + 0.05, max_w)
        }
        Err(_) => (false, f64::INFINITY),
    }
}
