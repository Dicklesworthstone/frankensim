//! Named dimensionless groups from ROLE-tagged inputs: Reynolds, Weber,
//! Capillary, Ohnesorge, Bond, Froude, Mach, Strouhal, Deborah, reduced
//! frequency, plus the structural set (slenderness, damping ratio,
//! P-Delta index). Every formula is evaluated through QtyAny dimension
//! arithmetic and REFUSED if the result is not dimensionless — groups
//! are dimensionless by construction, not by trust.

use crate::RegimeError;
use fs_math::det;
use fs_qty::QtyAny;

/// The physical role a quantity plays in group formation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Role {
    /// Mass density ρ (kg/m³).
    Density,
    /// Characteristic speed U (m/s).
    Velocity,
    /// Characteristic length L (m).
    Length,
    /// Dynamic viscosity μ (Pa·s).
    DynViscosity,
    /// Surface tension σ (N/m).
    SurfaceTension,
    /// Gravitational acceleration g (m/s²).
    Gravity,
    /// Speed of sound c (m/s).
    SoundSpeed,
    /// Characteristic frequency f (Hz).
    Frequency,
    /// Material relaxation time λ (s).
    RelaxationTime,
    /// Radius of gyration r (m) — structural slenderness.
    GyrationRadius,
    /// Viscous damping coefficient c (N·s/m).
    Damping,
    /// Stiffness k (N/m).
    Stiffness,
    /// Mass m (kg).
    Mass,
    /// Axial (gravity) load P (N) — P-Delta.
    AxialLoad,
    /// Lateral stiffness times height V·h proxy: story shear stiffness
    /// force (N) — P-Delta denominator.
    StoryShear,
}

impl Role {
    /// Stable report tag.
    #[must_use]
    pub fn tag(self) -> &'static str {
        match self {
            Role::Density => "density",
            Role::Velocity => "velocity",
            Role::Length => "length",
            Role::DynViscosity => "viscosity",
            Role::SurfaceTension => "surface-tension",
            Role::Gravity => "gravity",
            Role::SoundSpeed => "sound-speed",
            Role::Frequency => "frequency",
            Role::RelaxationTime => "relaxation-time",
            Role::GyrationRadius => "gyration-radius",
            Role::Damping => "damping",
            Role::Stiffness => "stiffness",
            Role::Mass => "mass",
            Role::AxialLoad => "axial-load",
            Role::StoryShear => "story-shear",
        }
    }
}

/// A role-tagged input.
#[derive(Debug, Clone, PartialEq)]
pub struct RoleInput {
    /// The role.
    pub role: Role,
    /// The dimensioned value (SI).
    pub qty: QtyAny,
}

/// A named, evaluated dimensionless group.
#[derive(Debug, Clone, PartialEq)]
pub struct NamedGroup {
    /// Canonical name ("Re", "We", …).
    pub name: &'static str,
    /// The value.
    pub value: f64,
}

struct Ctx<'a> {
    inputs: &'a [RoleInput],
}

impl Ctx<'_> {
    fn get(&self, role: Role) -> Option<QtyAny> {
        self.inputs.iter().find(|i| i.role == role).map(|i| i.qty)
    }
}

fn dimensionless(name: &'static str, q: QtyAny) -> Result<NamedGroup, RegimeError> {
    if q.dims.is_none() {
        Ok(NamedGroup {
            name,
            value: q.value,
        })
    } else {
        Err(RegimeError::NotDimensionless {
            context: format!("group {name}"),
            residual: q.dims.0,
        })
    }
}

/// Evaluate every standard group whose roles are present. Missing roles
/// simply skip a group; a present-but-inconsistent combination REFUSES.
///
/// # Errors
/// [`RegimeError::NotDimensionless`] when a formula's dimensional
/// self-check fails (wrongly-tagged input).
pub fn standard_groups(inputs: &[RoleInput]) -> Result<Vec<NamedGroup>, RegimeError> {
    let ctx = Ctx { inputs };
    let mut out = Vec::new();
    let (rho, u, l) = (
        ctx.get(Role::Density),
        ctx.get(Role::Velocity),
        ctx.get(Role::Length),
    );
    if let (Some(rho), Some(u), Some(l), Some(mu)) = (rho, u, l, ctx.get(Role::DynViscosity)) {
        out.push(dimensionless("Re", rho * u * l / mu)?);
    }
    if let (Some(rho), Some(u), Some(l), Some(sigma)) = (rho, u, l, ctx.get(Role::SurfaceTension)) {
        out.push(dimensionless("We", rho * u * u * l / sigma)?);
    }
    if let (Some(mu), Some(u), Some(sigma)) = (
        ctx.get(Role::DynViscosity),
        u,
        ctx.get(Role::SurfaceTension),
    ) {
        out.push(dimensionless("Ca", mu * u / sigma)?);
    }
    if let (Some(mu), Some(rho), Some(sigma), Some(l)) = (
        ctx.get(Role::DynViscosity),
        rho,
        ctx.get(Role::SurfaceTension),
        l,
    ) {
        // Oh = μ / sqrt(ρ σ L): compute as sqrt(μ² / (ρσL)) to stay in
        // dimension-checked arithmetic.
        let oh_sq = (mu * mu) / (rho * sigma * l);
        let checked = dimensionless("Oh", oh_sq)?;
        out.push(NamedGroup {
            name: "Oh",
            value: det::sqrt(checked.value),
        });
    }
    if let (Some(rho), Some(g), Some(l), Some(sigma)) = (
        rho,
        ctx.get(Role::Gravity),
        l,
        ctx.get(Role::SurfaceTension),
    ) {
        out.push(dimensionless("Bo", rho * g * l * l / sigma)?);
    }
    if let (Some(u), Some(g), Some(l)) = (u, ctx.get(Role::Gravity), l) {
        let fr_sq = (u * u) / (g * l);
        let checked = dimensionless("Fr", fr_sq)?;
        out.push(NamedGroup {
            name: "Fr",
            value: det::sqrt(checked.value),
        });
    }
    if let (Some(u), Some(c)) = (u, ctx.get(Role::SoundSpeed)) {
        out.push(dimensionless("Ma", u / c)?);
    }
    if let (Some(f), Some(l), Some(u)) = (ctx.get(Role::Frequency), l, u) {
        out.push(dimensionless("St", f * l / u)?);
    }
    if let (Some(lam), Some(u), Some(l)) = (ctx.get(Role::RelaxationTime), u, l) {
        out.push(dimensionless("De", lam * u / l)?);
    }
    if let (Some(l), Some(r)) = (l, ctx.get(Role::GyrationRadius)) {
        out.push(dimensionless("slenderness", l / r)?);
    }
    if let (Some(c), Some(k), Some(m)) = (
        ctx.get(Role::Damping),
        ctx.get(Role::Stiffness),
        ctx.get(Role::Mass),
    ) {
        // ζ = c / (2√(km)): dimension-check ζ² = c²/(4km).
        let zeta_sq = (c * c) / (k * m);
        let checked = dimensionless("zeta", zeta_sq)?;
        out.push(NamedGroup {
            name: "zeta",
            value: det::sqrt(checked.value) / 2.0,
        });
    }
    if let (Some(p), Some(v)) = (ctx.get(Role::AxialLoad), ctx.get(Role::StoryShear)) {
        out.push(dimensionless("theta-pdelta", p / v)?);
    }
    Ok(out)
}
