//! Region-correct pressure and particle velocity from the original observers.
use super::*;
use crate::near_field::{Geometry, Observation, VectorObservation};

pub(super) fn region_geometry<'a>(
    surfaces: [&'a SpherePanels; 2],
    region: CoreShellRegion,
    points: &[[f64; 3]],
    clearance: f64,
    gate: &CancelGate,
) -> Result<[Geometry<'a>; 2], HelmholtzError> {
    let outer = match region {
        CoreShellRegion::Exterior => {
            Geometry::new_iterative(surfaces[0], points, clearance, gate)?
        }
        CoreShellRegion::Shell | CoreShellRegion::Core => {
            Geometry::new_interior(surfaces[0], points, clearance)?
        }
    };
    poll(gate)?;
    let core = match region {
        CoreShellRegion::Core => {
            Geometry::new_interior(surfaces[1], points, clearance)?
        }
        CoreShellRegion::Exterior | CoreShellRegion::Shell => {
            Geometry::new_iterative(surfaces[1], points, clearance, gate)?
        }
    };
    poll(gate)?;
    Ok([outer, core])
}

/// Shared representation topology. Signs follow the ORIGINAL outward mesh
/// normals; eliminating a transparent surface never alters region admission.
pub(super) struct Contribution {
    pub boundary: usize,
    pub region: usize,
    pub sign: f64,
    pub scattered: bool,
}
pub(super) fn region_plan(
    region: CoreShellRegion,
    offsets: [Option<usize>; 2],
) -> ([Option<Contribution>; 2], bool) {
    let [outer, core] = offsets.map(|offset| offset.is_some());
    let term = |boundary, region, sign, scattered| Some(Contribution { boundary, region, sign, scattered });
    match region {
        CoreShellRegion::Exterior => {
            let boundary = if outer { 0 } else { 1 };
            ([if outer || core { term(boundary, 0, 1.0, true) } else { None }, None], true)
        }
        CoreShellRegion::Shell if outer => (
            [term(0, 1, -1.0, false), if core { term(1, 1, 1.0, false) } else { None }], false,
        ),
        CoreShellRegion::Shell => ([if core { term(1, 0, 1.0, true) } else { None }, None], true),
        CoreShellRegion::Core if core => ([term(1, 2, -1.0, false), None], false),
        CoreShellRegion::Core if outer => ([term(0, 2, -1.0, false), None], false),
        CoreShellRegion::Core => ([None, None], true),
    }
}

impl CoreShellScatteredField<'_, '_> {
    /// Observe the complete field in one declared physical region. Winding
    /// and face/edge distance enforce the region using the original meshes.
    /// Shell observations use BOTH original boundaries with opposite signs,
    /// the shell's actual k/rho, and one shared quadrature-work allowance.
    ///
    /// The existing adaptive triangle owner supplies pressure and analytic
    /// particle velocity together. Returned discrepancy estimates add absolute
    /// contributions from both surfaces; they do not bound mesh/model error
    /// or guarantee relative accuracy of a nearly cancelling total field.
    pub fn pressure_velocity_at(
        &self,
        region: CoreShellRegion,
        points_m: &[[f64; 3]],
        minimum_clearance_m: f64,
        options: near_field::Options,
        gate: &CancelGate,
    ) -> Result<VectorObservation, HelmholtzError> {
        poll(gate)?;
        options.validate()?;
        let owner = self.interface;
        let geometries = region_geometry(owner.surfaces, region, points_m, minimum_clearance_m, gate)?;
        let (contributions, add_incident) = region_plan(region, owner.offsets);
        let mut remaining = options.maximum_kernel_evaluations;
        let mut result = None;
        for contribution in contributions.into_iter().flatten() {
            let (pressure, velocity) = if contribution.scattered {
                (&self.scattered.pressure, &self.scattered.velocity)
            } else {
                let trace = self.traces[contribution.boundary].as_ref().expect("active interface");
                (&trace.pressure, &trace.velocity)
            };
            let observed = self.observe_geometry(
                &geometries[contribution.boundary],
                contribution.boundary,
                contribution.region,
                pressure,
                velocity,
                contribution.sign,
                options,
                &mut remaining,
                gate,
            )?;
            if let Some(total) = &mut result {
                accumulate(total, observed)?;
            } else {
                result = Some(observed);
            }
        }
        let mut result = result.unwrap_or_else(|| zero_observation(points_m.len()));
        if add_incident {
            for (index, &point) in points_m.iter().enumerate() {
                poll(gate)?;
                let sample = self.incident.pressure_velocity_at(
                    owner.k[0],
                    medium(owner.states[0]),
                    point,
                )?;
                result.scalar.pressure[index] = result.scalar.pressure[index] + sample.pressure_pa;
                for axis in 0..3 {
                    result.particle_velocity_m_s[index][axis] =
                        result.particle_velocity_m_s[index][axis] + sample.velocity_m_s[axis];
                }
            }
        }
        validate_observation(&result)?;
        poll(gate)?;
        Ok(result)
    }

    #[allow(clippy::too_many_arguments)]
    fn observe_geometry(
        &self,
        geometry: &Geometry<'_>,
        boundary: usize,
        region: usize,
        pressure: &[C64],
        velocity: &[C64],
        sign: f64,
        mut options: near_field::Options,
        remaining: &mut usize,
        gate: &CancelGate,
    ) -> Result<VectorObservation, HelmholtzError> {
        poll(gate)?;
        if *remaining < 80 {
            return Err(bad(
                "core-shell total receiver kernel-evaluation budget exhausted",
            ));
        }
        options.maximum_kernel_evaluations = *remaining;
        let owner = self.interface;
        let representation = owner.radiation(
            boundary,
            region,
            pressure.iter().map(|p| p.scale(sign)).collect(),
            velocity.iter().map(|v| v.scale(sign)).collect(),
        )?;
        let prepared = geometry.prepare_velocity_cancellable(
            owner.k[region],
            medium(owner.states[region]),
            options,
            gate,
        )?;
        *remaining -= prepared.kernel_evaluations;
        let result = prepared.evaluate_velocity(&representation)?;
        poll(gate)?;
        Ok(result)
    }
}

pub(super) fn zero_observation(count: usize) -> VectorObservation {
    VectorObservation {
        scalar: Observation {
            pressure: vec![C64::ZERO; count],
            quadrature_error_estimate_pa: vec![0.0; count],
        },
        particle_velocity_m_s: vec![[C64::ZERO; 3]; count],
        quadrature_error_estimate_m_s: vec![[0.0; 3]; count],
    }
}

pub(super) fn accumulate(
    target: &mut VectorObservation,
    other: VectorObservation,
) -> Result<(), HelmholtzError> {
    for index in 0..target.scalar.pressure.len() {
        target.scalar.pressure[index] =
            target.scalar.pressure[index] + other.scalar.pressure[index];
        target.scalar.quadrature_error_estimate_pa[index] +=
            other.scalar.quadrature_error_estimate_pa[index];
        for axis in 0..3 {
            target.particle_velocity_m_s[index][axis] = target.particle_velocity_m_s[index][axis]
                + other.particle_velocity_m_s[index][axis];
            target.quadrature_error_estimate_m_s[index][axis] +=
                other.quadrature_error_estimate_m_s[index][axis];
        }
    }
    validate_observation(target)
}

pub(super) fn validate_observation(field: &VectorObservation) -> Result<(), HelmholtzError> {
    if field
        .scalar
        .pressure
        .iter()
        .chain(field.particle_velocity_m_s.iter().flatten())
        .any(|&value| !finite(value))
        || field
            .scalar
            .quadrature_error_estimate_pa
            .iter()
            .chain(field.quadrature_error_estimate_m_s.iter().flatten())
            .any(|error| !error.is_finite() || *error < 0.0)
    {
        Err(bad("nonfinite core-shell receiver field or discrepancy"))
    } else {
        Ok(())
    }
}
