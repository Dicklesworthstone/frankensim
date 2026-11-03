//! Prescribed and finite-inertia film surfaces sharing one native exchange.

use super::*;

/// A single planar film port in the simultaneous sphere/fluid response.
///
/// Each finite variant owns one independent carrier's spatial inertia. Equal
/// body values in different slots do not identify a shared multi-face carrier.
/// Prescribed variants retain their explicit translation and optional axial
/// spin. These remain missing-film surfaces, not resolved fluid boundaries.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum LubricationSurface3 {
    /// Driven geometry/rate; received film work belongs to its actuator port.
    Prescribed(LubricationPlane3),
    /// Finite spatial inertia; received film work changes its kinetic state.
    Finite(FiniteLubricationPlane3),
}

impl BoundaryGrid3 {
    /// Exchange spheres, fluid and mixed prescribed/finite film surfaces once.
    ///
    /// Every active sphere pair and surface film enters the SAME original
    /// implicit midpoint resistance solve. A finite carrier's recoil therefore
    /// changes the simultaneous reactions on prescribed surfaces and other
    /// carriers. There is one native collision/stream, not one per surface.
    ///
    /// `surface_external` has exactly one center impulse followed by center
    /// couple per surface, in the same stable order. Loads on finite carriers
    /// are independent of their film wrench. Prescribed surfaces require exact
    /// zero loads here: their known motion and received film work already define
    /// the driving-reservoir port; accepting a second impulse would hide its
    /// momentum destination.
    ///
    /// Each returned plane reaction retains its input position. A finite
    /// response has `finite_body: Some(_)`, including its actual kinetic change
    /// and external work. A prescribed response has `finite_body: None`; its
    /// received work is supplied with opposite sign by its actuator. The
    /// aggregate `coupling.lubrication.plane_work` includes BOTH categories and
    /// must not be treated entirely as external reservoir work.
    ///
    /// All spheres, finite carrier velocities and fluid populations publish
    /// together. Prescribed descriptors are unchanged. An empty surface slice
    /// retains the original sphere-only result, and homogeneous slices retain
    /// their existing native arithmetic. One bounded temporary descriptor array
    /// is prepared; the shared solver retains its original film/mode budgets.
    ///
    /// Geometry and world inertia are frozen, as in the existing finite-plane
    /// owner. No pose drift, gyroscopic update, solid contact, resolved moving
    /// wall or conversion of film loss into heat is introduced.
    ///
    /// # Errors
    /// Refuses mismatched load counts, more than 64 surfaces, nonfinite loads,
    /// a nonzero prescribed-surface load, and all original native geometry,
    /// inertia, positive-gap, speed, mode-budget and cancellation failures.
    /// Refusal preserves every accepted sphere, surface and population.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_lubricated_partial_saturation_surfaces<const N: usize>(
        &mut self,
        fractions: [&[f64]; N],
        radii: [f64; N],
        bodies: &mut [PartialSaturationBody3; N],
        external: [[f64; 6]; N],
        model: SphereLubrication3,
        surfaces: &mut [LubricationSurface3],
        surface_external: &[[f64; 6]],
        origin: [f64; 3],
        gate: &CancelGate,
    ) -> Result<LubricatedPlanesExchange3<N>> {
        let native = prepare(surfaces, surface_external, gate)?;
        let (coupling, (), planes) = self.step_with_partial_saturation_resistance_checked(
            fractions,
            bodies,
            external,
            origin,
            gate,
            Some((model, radii, &native)),
            |_, _, _| Ok(()),
        )?;
        accept(surfaces, &planes);
        Ok(LubricatedPlanesExchange3 { coupling, planes })
    }
}

fn prepare(
    surfaces: &[LubricationSurface3],
    external: &[[f64; 6]],
    gate: &CancelGate,
) -> Result<Vec<LubricationPlane3>> {
    checkpoint(gate)?;
    if surfaces.len() > 64 || external.len() != surfaces.len() {
        return Err(invalid(
            "invalid mixed surface count or external impulse count",
        ));
    }
    let mut native = Vec::new();
    native
        .try_reserve_exact(surfaces.len())
        .map_err(|_| invalid("mixed surface descriptor allocation refused"))?;
    for (&surface, &impulse) in surfaces.iter().zip(external) {
        checkpoint(gate)?;
        native.push(match surface {
            LubricationSurface3::Prescribed(plane) => {
                if impulse.into_iter().any(|x| x != 0.0) {
                    return Err(invalid(
                        "prescribed surface cannot receive an authored impulse",
                    ));
                }
                plane
            }
            LubricationSurface3::Finite(plane) => plane.driven(impulse)?,
        });
    }
    checkpoint(gate)?;
    Ok(native)
}

// The checked common owner produces exactly one correctly tagged response for
// each descriptor, before publishing fluid/sphere state. This phase is infallible
// and performs only the finite state assignments; no new cancellation boundary.
fn accept(surfaces: &mut [LubricationSurface3], reports: &[PlaneLubricationExchange3]) {
    for (surface, report) in surfaces.iter_mut().zip(reports) {
        if let LubricationSurface3::Finite(plane) = surface {
            plane.accept(
                report
                    .finite_body
                    .expect("every finite descriptor has its checked response")
                    .body_after,
            );
        }
    }
}

#[cfg(test)]
mod tests;
