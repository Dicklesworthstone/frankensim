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
            // Only coefficient weighting changes during this observation.
            // Velocity, wall traces, original inertia and TimeLevels.previous
            // stay untouched. Restore even when observation is cancelled.
            std::mem::swap(&mut solver.molecular_viscosity, &mut unit_viscosity);
            let observed = solver.molecular_shear_rates(gate);
            std::mem::swap(&mut solver.molecular_viscosity, &mut unit_viscosity);
            let rates = observed?;
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
