//! One reference-configuration load law for nodal assembly and weak residuals.
//! Callbacks are pure, mesh-independent physical densities. No nodal force
//! transfer or implicit density scaling is used when rebuilding on another grid.
use super::*;

/// Optional surface contribution per reference area. Pressure uses the same
/// inward-positive convention as `pressure_load`; neither law follows motion.
#[derive(Clone, Copy)]
pub enum SurfaceForce3<'a> {
    /// Traction from reference position and outward unit normal.
    Traction(&'a dyn Fn([f64; 3], [f64; 3]) -> [f64; 3]),
    /// Signed pressure; positive acts into the solid.
    Pressure(&'a dyn Fn([f64; 3]) -> f64),
}
impl SurfaceForce3<'_> {
    pub(in crate::elastic3) fn value(self, p: [f64; 3], n: [f64; 3]) -> [f64; 3] {
        match self {
            Self::Traction(f) => f(p, n),
            Self::Pressure(f) => { let pressure = f(p); n.map(|v| -pressure*v) }
        }
    }
}

/// Sum of a volume force and an optional force on the retained implicit surface.
/// The same law can also define a linear observation `q(u)`. For an observation,
/// its callback units are chosen to give the desired functional, not necessarily
/// force units. `Default` is an explicit zero functional. No surface quadrature
/// is required when `surface` is None, even for a surface-capable operator.
#[derive(Clone, Copy, Default)]
pub struct ReferenceLoad3<'a> {
    /// Body force per reference volume, or volume-observation density.
    pub body: Option<&'a dyn Fn([f64; 3]) -> [f64; 3]>,
    /// Surface traction/pressure, or surface-observation density.
    pub surface: Option<SurfaceForce3<'a>>,
}
impl<'a> ReferenceLoad3<'a> {
    /// Preserve the original body-only path, without surface work.
    #[must_use]
    pub const fn body(body: &'a dyn Fn([f64; 3]) -> [f64; 3]) -> Self {
        Self { body: Some(body), surface: None }
    }
    /// A pure surface traction/observation.
    #[must_use]
    pub const fn traction(traction: &'a dyn Fn([f64; 3], [f64; 3]) -> [f64; 3]) -> Self {
        Self { body: None, surface: Some(SurfaceForce3::Traction(traction)) }
    }
    /// A pure inward-positive reference pressure.
    #[must_use]
    pub const fn pressure(pressure: &'a dyn Fn([f64; 3]) -> f64) -> Self {
        Self { body: None, surface: Some(SurfaceForce3::Pressure(pressure)) }
    }
    pub(in crate::elastic3) fn admit(self, op: &CutElasticity3,
        checkpoint: &mut impl FnMut() -> ControlFlow<()>) -> Result<(), ElasticityError3> {
        surface_poll(checkpoint)?;
        if self.surface.is_some() {
            for cell in &op.cells {
                surface_poll(checkpoint)?;
                if cell.rules.surface().is_none() {
                    return Err(ElasticityError3::Invalid("reference surface law requires build_with_surface"));
                }
            }
        }
        Ok(())
    }
}
impl CutElasticity3 {
    /// Assemble one mixed load using the existing body and surface integrators.
    /// Missing surface support is refused before invoking either callback.
    /// Geometry and stiffness are never changed; an error returns no partial RHS.
    pub fn reference_load(&self, law: ReferenceLoad3<'_>,
        mut checkpoint: impl FnMut() -> ControlFlow<()>) -> Result<Vec<f64>, ElasticityError3> {
        law.admit(self, &mut checkpoint)?;
        let body = match law.body {
            Some(f) => Some(self.body_load(f, &mut checkpoint)?), None => None,
        };
        let surface = match law.surface {
            Some(f) => Some(self.surface_load(&|p, n| f.value(p, n), &mut checkpoint)?.rhs), None => None,
        };
        let result = match (body, surface) {
            (Some(mut b), Some(s)) => {
                for (i, (b, s)) in b.iter_mut().zip(s).enumerate() {
                    if i % 192 == 0 { surface_poll(&mut checkpoint)?; }
                    *b += s;
                }
                b
            }
            (Some(b), None) => b,
            (None, Some(s)) => s,
            (None, None) => vec![0.0; self.n()],
        };
        if !result.iter().all(|v| v.is_finite()) { return Err(ElasticityError3::Invalid("combined reference load overflow")); }
        surface_poll(&mut checkpoint)?;
        Ok(result)
    }
}
