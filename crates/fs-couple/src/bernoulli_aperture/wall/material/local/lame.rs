//! Lame's thick-walled cylinder as the radial image of an isotropic wall that is too thick for the
//! membrane law.
//!
//! `LocalMaterialWall::base_patch` calls this only above thickness/radius = 0.05; at or under it the
//! original membrane expressions run unchanged. The wall is a frozen-state, quasi-static elastic
//! cylinder with inner radius `a` (the bore), outer radius `b = a + h`, a free outer face and the
//! gas pressure `p` on the inner face. Its radial field is `u(r) = A r + B / r` with
//! `sigma_r(b) = 0` and `sigma_r(a) = -p`. Textbook source: Timoshenko and Goodier, *Theory of
//! Elasticity*, the thick-walled cylinder under internal pressure. Every input is a resolved
//! material number, a printed thickness and the bore radius; nothing is fitted.
//!
//! * Stress-free axial restraint is plane stress: `u(a) = (p a / E) [ (a^2 + b^2) / (b^2 - a^2) + nu ]`.
//! * Strain-fixed axial restraint is plane strain: the same form with `E' = E / (1 - nu^2)` and
//!   `nu' = nu / (1 - nu)`, which equals the textbook plane-strain form
//!   `(1 + nu) (p a / E) [ (1 - 2 nu) a^2 + b^2 ] / (b^2 - a^2)`.
//! * The stiffness per unit inner area is `k = p / u(a)`. It is returned divided by `E`, so the
//!   Kelvin-Voigt resistance and every Maxwell arm scale exactly as the membrane's do (the
//!   correspondence principle at a constant Poisson ratio).
//! * The surface density is the kinetic energy of the SAME field per unit inner area,
//!   `m = rho (integral of u^2 r dr from a to b) / (a u(a)^2)`, in closed form. This is a Rayleigh
//!   estimate: it holds below the wall's through-thickness resonances and says nothing about them.
//!
//! As `h / a` tends to 0 the two outputs tend to the membrane's `(h / a) / a` (times the axial
//! factor) and `rho h`.
use super::super::AxialConstraint;

/// `(stiffness_per_area / E, surface_density)` of an isotropic wall of the given thickness on a bore
/// of the given radius. `poisson` must lie in `(-1, 0.5)`, the caller's isotropic admission.
pub(super) fn thick_cylinder(
    inner_radius_m: f64,
    thickness_m: f64,
    density: f64,
    poisson: f64,
    axial: AxialConstraint,
) -> (f64, f64) {
    let (a, h) = (inner_radius_m, thickness_m);
    let b = a + h;
    let (modulus_scale, nu) = match axial {
        AxialConstraint::StressFree => (1.0, poisson),
        AxialConstraint::StrainFixed => (1.0 / (1.0 - poisson * poisson), poisson / (1.0 - poisson)),
    };
    // b^2 - a^2 without cancellation.
    let annulus = h * (2.0 * a + h);
    let hoop = (a * a + b * b) / annulus;
    let stiffness_over_young = modulus_scale / (a * (hoop + nu));
    // u = A r + B / r with B = 1 (the scale cancels) and sigma_r(b) = 0.
    let big_a = (1.0 - nu) / ((1.0 + nu) * b * b);
    let inner = big_a * a + 1.0 / a;
    let kinetic = big_a * big_a * (b * b * b * b - a * a * a * a) / 4.0
        + big_a * annulus
        + (h / a).ln_1p();
    (stiffness_over_young, density * kinetic / (a * inner * inner))
}

#[cfg(test)]
mod tests;
