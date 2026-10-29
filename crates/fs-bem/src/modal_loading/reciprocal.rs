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
///
/// Inside [`helmholtz::with_dense_workers`] the frequencies are solved on the
/// declared workers; the result, and the refusal returned, are bit-identical
/// for every count. Each frequency in flight holds its own dense operator and
/// factor, about 32 n^2 bytes for n [`ReciprocalLoadingDiagnostics::unknowns`].
/// The caller's cancellation is then polled about every 20 ms while they run.
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
    // One frequency's whole work: its solutions' checks, in order, and its row.
    let frequency = |omega:f64,cancelled:&mut dyn FnMut()->bool|
        -> Result<(Vec<[f64;3]>,ModalLoadingRow),ModalLoadingError> {
        let source = owner.solve_velocity_ports(omega/medium.sound_speed,medium,&refs,cancelled)?;
        let checks = source.solutions().iter().map(|solution|
            [solution.panels_per_wavelength,solution.backward_residual,solution.flux_balance]).collect();
        Ok((checks,assess(omega,source.impedance().to_vec(),modes,None)?))
    };
    let workers = helmholtz::dense_workers().min(omegas.len());
    let mut solved = if workers > 1 {
        on_workers(omegas,workers,&mut cancelled,&frequency)?
    } else { Vec::new() }.into_iter();
    for &omega in omegas {
        let (checks,row) = if workers > 1 {
            solved.next().expect("every frequency up to the first failure has its outcome")?
        } else {
            poll(&mut cancelled)?;
            frequency(omega,&mut cancelled)?
        };
        for [panels_per_wavelength,backward_residual,flux_balance] in checks {
            diagnostics.minimum_panels_per_wavelength = diagnostics.minimum_panels_per_wavelength
                .min(panels_per_wavelength);
            diagnostics.maximum_backward_residual = diagnostics.maximum_backward_residual
                .max(backward_residual);
            diagnostics.maximum_flux_balance = diagnostics.maximum_flux_balance.max(flux_balance);
        }
        rows.push(row);
    }
    poll(&mut cancelled)?;
    Ok(ReciprocalLoadingBake {
        loading:ModalLoadingBake {modes,declaration_hash:hash.finalize(),rows},diagnostics,
    })
}

/// `solve` every frequency on `workers` scoped threads, each taking the lowest
/// unclaimed frequency. Outcomes return in frequency order, up to and including
/// the first failure. That is the failure the serial loop returns: a failure
/// abandons only higher frequencies. Spare workers go to each frequency's
/// factorization. This thread polls the caller's cancellation while it waits;
/// a cancellation stops every frequency and returns `Cancelled`.
fn on_workers<T:Send>(omegas:&[f64],workers:usize,cancelled:&mut impl FnMut()->bool,
    solve:&(impl Fn(f64,&mut dyn FnMut()->bool)->Result<T,ModalLoadingError>+Sync))
    ->Result<Vec<Result<T,ModalLoadingError>>,ModalLoadingError> {
    use std::sync::atomic::{AtomicBool,AtomicUsize,Ordering::Relaxed};
    use std::sync::mpsc::{RecvTimeoutError,channel};
    let (next,first_failure,stop)=(AtomicUsize::new(0),AtomicUsize::new(usize::MAX),AtomicBool::new(false));
    let (next,first_failure,stop)=(&next,&first_failure,&stop);
    let inner=helmholtz::dense_workers()/workers;
    let mut outcomes:Vec<Option<Result<T,ModalLoadingError>>>=omegas.iter().map(|_|None).collect();
    std::thread::scope(|scope| {
        // Nothing is sent: the channel disconnects when every worker is done.
        let (running,finished)=channel::<()>();
        let handles:Vec<_>=(0..workers).map(|_| {
            let running=running.clone();
            scope.spawn(move || {
                let _running=running;
                let mut mine=Vec::new();
                loop {
                    let index=next.fetch_add(1,Relaxed);
                    if index>=omegas.len() || index>first_failure.load(Relaxed) {break mine;}
                    let outcome=helmholtz::with_dense_workers(inner,||
                        solve(omegas[index],&mut || stop.load(Relaxed) || index>first_failure.load(Relaxed)));
                    if outcome.is_err() {first_failure.fetch_min(index,Relaxed);}
                    mine.push((index,outcome));
                }
            })
        }).collect();
        drop(running);
        loop {
            if !stop.load(Relaxed) && cancelled() {stop.store(true,Relaxed);}
            if let Err(RecvTimeoutError::Disconnected)=finished.recv_timeout(std::time::Duration::from_millis(20)) {break;}
        }
        for handle in handles {
            let mine=handle.join().unwrap_or_else(|panic|std::panic::resume_unwind(panic));
            for (index,outcome) in mine {outcomes[index]=Some(outcome);}
        }
    });
    if stop.load(Relaxed) {return Err(HelmholtzError::Cancelled.into());}
    let mut ordered=Vec::with_capacity(omegas.len());
    for outcome in outcomes {
        let outcome=outcome.expect("every frequency below the first failure ran");
        let failed=outcome.is_err();
        ordered.push(outcome);
        if failed {break;}
    }
    Ok(ordered)
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

    #[test]
    fn declared_workers_solve_the_frequencies_with_the_same_bits_and_cancel_cleanly() {
        let surface = SpherePanels::icosphere(0.1,0).unwrap();
        let mut shapes = vec![vec![0.0;20];3];
        shapes[0][0] = 1.0; shapes[1][7] = -0.7; shapes[2][12] = 0.4;
        let budget = Budget {maximum_unknowns:64,maximum_cache_bytes:1_000_000,
            maximum_work:100_000_000,angular_refinement:1};
        let omega = [200.0,350.0,500.0,650.0,800.0];
        let medium = Medium::air();
        let bake = |cancelled:&mut dyn FnMut()->bool| modal_radiation_loading_reciprocal(&surface,&shapes,
            &omega,medium,budget,fs_blake3::hash_domain("reciprocal-loading-test",b"declared workers"),45,cancelled);
        let bits = |bake:&ReciprocalLoadingBake| (bake.loading.declaration_hash.as_bytes().to_vec(),
            bake.loading.rows().iter().map(|row|row.impedance().iter()
                .map(|z|(z.re.to_bits(),z.im.to_bits())).collect::<Vec<_>>()).collect::<Vec<_>>(),
            [bake.diagnostics.minimum_panels_per_wavelength,bake.diagnostics.maximum_backward_residual,
                bake.diagnostics.maximum_flux_balance].map(f64::to_bits));
        let serial = bits(&bake(&mut ||false).unwrap());
        for workers in [2,3,5,8] {
            let parallel = helmholtz::with_dense_workers(workers,|| bake(&mut ||false)).unwrap();
            assert_eq!(bits(&parallel),serial,"{workers} workers");
        }
        // Cancelled once the sweep has started: no bake, and the next one is unchanged.
        let mut preamble = 0;
        Prepared::new(&surface,omega[4]/medium.sound_speed,budget,|| {preamble+=1;false}).unwrap();
        let mut calls = 0;
        let cancelled = helmholtz::with_dense_workers(3,|| bake(&mut || {calls+=1;calls>1+preamble}));
        assert!(matches!(cancelled,Err(ModalLoadingError::Helmholtz(HelmholtzError::Cancelled))));
        assert_eq!(calls,2+preamble);
        assert_eq!(bits(&helmholtz::with_dense_workers(3,|| bake(&mut ||false)).unwrap()),serial);
    }
}
