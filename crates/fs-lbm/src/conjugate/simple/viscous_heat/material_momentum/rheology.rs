//! Endpoint generalized-Newtonian coefficients on one retained native solver.

use super::*;

/// Original physical inventory and native actual-mass geometry. Constitutive
/// coefficients are resolved from the native velocity, never cell-P inversion.
pub(in crate::conjugate) struct NativeMaterialRheologyInput<'a> {
    pub field: &'a FlowField,
    pub mass_kg: &'a [f64],
    pub momentum_kg_m_s: &'a [[f64; 3]],
    pub face_mass_kg: [Vec<f64>; 3],
    pub boundaries: [FvBoundary; 6],
    pub pressure_mouth_axis: Option<usize>,
    pub duration_s: f64,
}

pub(in crate::conjugate) struct NativeMaterialRheologyStep {
    pub viscosity: NativeMaterialViscosityStep,
    /// Exact unit-quadrature rates at the accepted native face velocity.
    pub endpoint_shear_rate_per_s: Vec<f64>,
    /// Initial source resolution plus every post-sweep endpoint resolution.
    pub source_updates: usize,
}

/// One bounded Picard solve of native endpoint constitutive momentum rows.
///
/// Every source update uses the exact current native rate. The shared kick
/// keeps its original velocity time level and actual face inertia throughout.
/// Refreshing both coefficients BEFORE the common residual readback measures
/// the nonlinear equations, not convergence of a previously frozen matrix.
/// The same admitted endpoint coefficients then own physical heat and work.
pub(in crate::conjugate) fn kick_with_source(
    domain: &VoxelDomain,
    input: NativeMaterialRheologyInput<'_>,
    controls: MaterialViscosityControls,
    gate: &CancelGate,
    resolve: impl FnMut(&[f64]) -> Result<(Vec<f64>, Vec<f64>), ChtError>,
) -> Result<NativeMaterialRheologyStep, ChtError> {
    kick_with_source_mode(domain, input, None, controls, gate, resolve)
}

/// Endpoint source rheology and resolved-gradient SGS on the SAME native time
/// level. The source sees only physical unit-quadrature strain. Its shear and
/// independent bulk coefficients are installed before refreshing SGS and its
/// material caps; molecular heat and modeled transfer remain separate.
pub(in crate::conjugate) fn kick_with_source_and_sgs(
    domain: &VoxelDomain,
    input: NativeMaterialRheologyInput<'_>,
    sgs: crate::conjugate::FvSgsConfig,
    controls: MaterialViscosityControls,
    gate: &CancelGate,
    resolve: impl FnMut(&[f64]) -> Result<(Vec<f64>, Vec<f64>), ChtError>,
) -> Result<NativeMaterialRheologyStep, ChtError> {
    poll(gate)?;
    sgs.validate()?;
    if sgs.model.coefficient() == 0.0 {
        return kick_with_source(domain, input, controls, gate, resolve);
    }
    kick_with_source_mode(domain, input, Some(sgs), controls, gate, resolve)
}

/// Temporarily expose the pure molecular unit operator, including reservoir
/// half-edges. Both the cell field AND the cached edge closures must be absent:
/// a unit molecular coefficient alone would still produce SGS-weighted rates.
/// Moving the buffers out and back retains all storage and restores the exact
/// accepted coefficients even when observation refuses or is cancelled.
fn observe_molecular_rate(
    solver: &mut Solver<'_>,
    unit_viscosity: &mut Vec<f64>,
    gate: &CancelGate,
) -> Result<Vec<f64>, ChtError> {
    let sgs = solver.sgs.take();
    let eddy = std::mem::take(&mut solver.eddy_viscosity);
    std::mem::swap(&mut solver.molecular_viscosity, unit_viscosity);
    let observed = solver.molecular_shear_rates(gate);
    std::mem::swap(&mut solver.molecular_viscosity, unit_viscosity);
    solver.eddy_viscosity = eddy;
    solver.sgs = sgs;
    observed
}

fn kick_with_source_mode(
    domain: &VoxelDomain,
    input: NativeMaterialRheologyInput<'_>,
    sgs: Option<crate::conjugate::FvSgsConfig>,
    controls: MaterialViscosityControls,
    gate: &CancelGate,
    mut resolve: impl FnMut(&[f64]) -> Result<(Vec<f64>, Vec<f64>), ChtError>,
) -> Result<NativeMaterialRheologyStep, ChtError> {
    poll(gate)?;
    let n = domain.cell_count();
    if n > controls.max_cells
        || controls.max_sweeps == 0
        || input.mass_kg.len() != n
        || input.momentum_kg_m_s.len() != n
    {
        return Err(refuse("endpoint rheology requires bounded original cells"));
    }
    let first = (0..n)
        .find(|&cell| domain.is_fluid(cell))
        .ok_or_else(|| refuse("endpoint rheology requires mobile material"))?;
    let density = input.mass_kg[first] / domain.volume(first);
    finite_positive("flow.rheology_reference_density", density)?;
    // A unit kinematic coefficient initializes the native owner without an
    // artificial restriction on reference density. It NEVER enters a solve:
    // the first refresh below installs the actual original-rate coefficients.
    let seed_viscosity = vec![density; n];
    let mut unit_viscosity = vec![1.0; n];
    let mut endpoint_shear_rate_per_s = Vec::new();
    let mut source_updates = 0_usize;
    let mut sgs_initialized = false;
    let mass = input.mass_kg;
    let viscosity = kick_with_refresh(
        domain,
        NativeMaterialViscosityInput {
            field: input.field,
            mass_kg: input.mass_kg,
            momentum_kg_m_s: input.momentum_kg_m_s,
            viscosity_pa_s: &seed_viscosity,
            bulk_viscosity_pa_s: None,
            face_mass_kg: input.face_mass_kg,
            boundaries: input.boundaries,
            pressure_mouth_axis: input.pressure_mouth_axis,
            duration_s: input.duration_s,
        },
        controls,
        gate,
        |solver| {
            // Velocity, wall traces, original inertia and TimeLevels.previous
            // stay untouched throughout this kinematic source observation.
            let rates = observe_molecular_rate(solver, &mut unit_viscosity, gate)?;
            poll(gate)?;
            let resolved = resolve(&rates);
            poll(gate)?;
            let (mut shear, bulk) = resolved?;
            if shear.len() != n || bulk.len() != n {
                return Err(refuse(
                    "endpoint rheology sources must match the original grid",
                ));
            }
            // The underlying transport owner admits one positive slot per
            // voxel. Immobile placeholders have no source or dissipation;
            // arbitrary callback values there are deliberately ignored.
            let reference = shear[first];
            for (cell, value) in shear.iter_mut().enumerate() {
                if cell % 256 == 0 {
                    poll(gate)?;
                }
                if !domain.is_fluid(cell) {
                    *value = reference;
                }
            }
            solver.set_material_transport(&shear, None, gate)?;
            set_bulk_viscosity(solver, &shear, &bulk, gate)?;
            if let Some(sgs) = sgs {
                if !sgs_initialized {
                    // The initialization seed never enters SGS admission.
                    // The first exact original source coefficients are already
                    // installed, including independently authored bulk response.
                    solver.set_sgs_density(mass, gate)?;
                    solver.enable_sgs(sgs, gate)?;
                    sgs_initialized = true;
                }
                // Refresh all cell, interior-edge and reservoir-edge closures
                // at this SAME endpoint and its freshly resolved molecular mu,
                // before any common momentum residual or power acceptance.
                solver.update_sgs(gate)?;
            }
            source_updates = source_updates
                .checked_add(1)
                .ok_or_else(|| refuse("endpoint rheology source work count overflow"))?;
            endpoint_shear_rate_per_s = rates;
            Ok(())
        },
    )?;
    poll(gate)?;
    Ok(NativeMaterialRheologyStep {
        viscosity,
        endpoint_shear_rate_per_s,
        source_updates,
    })
}

#[cfg(test)]
mod tests;
