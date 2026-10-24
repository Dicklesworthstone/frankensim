//! Original reciprocal P1/P0 source, under the SAME unchanged matrix admission.
use super::*;
use crate::reciprocal_radiation::{Budget, Prepared};

/// Numerical diagnostics from the actual variational source, not collocation
/// condition probes. Flux/residual checks concern this discrete operator.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ReciprocalLoadingDiagnostics {
    pub angular_refinement: u32,
    pub unknowns: usize,
    pub cache_bytes: usize,
    pub frequency_work: u64,
    pub sweep_work: u64,
    pub minimum_panels_per_wavelength: f64,
    pub maximum_backward_residual: f64,
    pub maximum_flux_balance: f64,
}

/// Complete unchanged impedance rows and the owner that computed their work.
#[derive(Clone, Debug)]
pub struct ReciprocalLoadingBake {
    pub loading: ModalLoadingBake,
    pub diagnostics: ReciprocalLoadingDiagnostics,
}

/// Compute every original signed body/aperture velocity column with one
/// symmetric variational source. The existing strict reciprocity/passivity
/// assessment is applied verbatim before publishing any row. No collocation
/// matrix is repaired and no causal/acoustic-accuracy claim follows from this.
/// The caller's per-frequency budget bounds preparation and the complete batch;
/// the returned sweep work additionally counts every declared frequency.
/// This is the original owner's bounded, below-resonance chart, with no fallback.
#[allow(clippy::too_many_arguments)]
pub fn modal_radiation_loading_reciprocal(
    surface: &SpherePanels, shapes: &[Vec<f64>], omegas: &[f64], medium: Medium,
    budget: Budget, provenance: ContentHash, maximum_matrix_entries: usize,
    mut cancelled: impl FnMut() -> bool,
) -> Result<ReciprocalLoadingBake, ModalLoadingError> {
    poll(&mut cancelled)?;
    let modes = shapes.len();
    if modes == 0 || modes > helmholtz::MAX_RADIATION_FIELDS_PER_BATCH
        || omegas.is_empty() || omegas.iter().any(|w|!w.is_finite() || *w <= 0.0)
        || omegas.windows(2).any(|w|w[0] >= w[1])
        || !medium.density.is_finite() || medium.density <= 0.0
        || !medium.sound_speed.is_finite() || medium.sound_speed <= 0.0
        || shapes.iter().any(|row|row.len() != surface.areas().len()
            || row.iter().any(|v|!v.is_finite()))
    { return Err(ModalLoadingError::Invalid("complete finite reciprocal source and increasing frequency grid required")); }
    let entries = modes.checked_mul(modes).and_then(|n|n.checked_mul(omegas.len()))
        .ok_or(ModalLoadingError::Invalid("reciprocal loading extent overflow"))?;
    if entries > maximum_matrix_entries {
        return Err(ModalLoadingError::Invalid("reciprocal loading exceeds matrix-entry allowance"));
    }
    fs_la::eigen::admit_jacobi_eigh(modes).map_err(|_|
        ModalLoadingError::Invalid("modal resistance eigenproblem exceeds work envelope"))?;
    let owner = Prepared::new(surface,omegas[omegas.len()-1]/medium.sound_speed,
        budget,&mut cancelled)?;
    let frequency_work = owner.velocity_port_work(modes)?;
    if frequency_work > budget.maximum_work {
        return Err(ModalLoadingError::Invalid("complete reciprocal source exceeds per-frequency allowance"));
    }
    let sweep_work = frequency_work.checked_mul(omegas.len() as u64)
        .ok_or(ModalLoadingError::Invalid("reciprocal sweep work overflow"))?;
    let mut hash = DomainHasher::new("fs-bem/complete-modal-loading/reciprocal-p1p0/v1");
    hash.update(provenance.as_bytes());
    for value in [surface.areas().len() as u64,modes as u64,omegas.len() as u64,
        maximum_matrix_entries as u64,u64::from(budget.angular_refinement)] {
        hash.update(&value.to_le_bytes());
    }
    let triangles = surface.triangles().ok_or(ModalLoadingError::Invalid(
        "reciprocal source requires its original triangles"))?;
    for value in [medium.density,medium.sound_speed].iter()
        .chain(omegas).chain(triangles.iter().flatten().flatten()).chain(shapes.iter().flatten()) {
        hash.update(&value.to_bits().to_le_bytes());
    }
    let refs: Vec<_> = shapes.iter().map(Vec::as_slice).collect();
    let mut diagnostics = ReciprocalLoadingDiagnostics {
        angular_refinement:budget.angular_refinement,unknowns:owner.unknowns(),
        cache_bytes:owner.cache_bytes(),frequency_work,sweep_work,
        minimum_panels_per_wavelength:f64::INFINITY,
        maximum_backward_residual:0.0,maximum_flux_balance:0.0,
    };
    let mut rows = reserved(omegas.len())?;
    for &omega in omegas {
        poll(&mut cancelled)?;
        let source = owner.solve_velocity_ports(omega/medium.sound_speed,medium,&refs,&mut cancelled)?;
        for solution in source.solutions() {
            diagnostics.minimum_panels_per_wavelength = diagnostics.minimum_panels_per_wavelength
                .min(solution.panels_per_wavelength);
            diagnostics.maximum_backward_residual = diagnostics.maximum_backward_residual
                .max(solution.backward_residual);
            diagnostics.maximum_flux_balance = diagnostics.maximum_flux_balance.max(solution.flux_balance);
        }
        rows.push(assess(omega,source.impedance().to_vec(),modes,None)?);
    }
    poll(&mut cancelled)?;
    Ok(ReciprocalLoadingBake {
        loading:ModalLoadingBake {modes,declaration_hash:hash.finalize(),rows},diagnostics,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn independent_localized_ports_keep_original_reciprocal_work_without_repair() {
        let sphere = SpherePanels::icosphere(0.1,0).unwrap();
        let surface = SpherePanels::from_triangles(sphere.triangles().unwrap().iter()
            .map(|t|t.map(|p|[1.17*p[0],0.81*p[1],p[2]])).collect()).unwrap();
        let mut shapes = vec![vec![0.0;20];2];
        shapes[0][0] = 1.0; shapes[1][7] = -0.7;
        let budget = Budget {maximum_unknowns:64,maximum_cache_bytes:1_000_000,
            maximum_work:100_000_000,angular_refinement:1};
        let omega = [300.0,900.0];
        let medium = Medium::air();
        let bake = modal_radiation_loading_reciprocal(&surface,&shapes,&omega,medium,budget,
            fs_blake3::hash_domain("reciprocal-loading-test",b"independent original faces"),8,||false).unwrap();
        let owner = Prepared::new(&surface,omega[1]/medium.sound_speed,budget,||false).unwrap();
        let refs: Vec<_> = shapes.iter().map(Vec::as_slice).collect();
        for (index,row) in bake.loading.rows().iter().enumerate() {
            let original = owner.solve_velocity_ports(omega[index]/medium.sound_speed,medium,&refs,||false).unwrap();
            assert_eq!(row.impedance(),original.impedance());
            assert!(row.reciprocity_relative_defect() <= DENSE_RECIPROCITY_BOUND);
            assert!(row.impedance()[1].abs() > 0.0);
        }
        assert_eq!(bake.diagnostics.sweep_work,2*bake.diagnostics.frequency_work);
        assert_eq!(bake.diagnostics.unknowns,owner.unknowns());
        assert_eq!(bake.diagnostics.cache_bytes,owner.cache_bytes());
    }
}
