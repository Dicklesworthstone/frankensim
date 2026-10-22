//! Restore the physical work account together with the constrained search.
//! The application reconstructs the same immutable model, then supplies its
//! identity and expected search declarations. No callback runs while loading.
use super::*;
pub use fs_dfo::differential_evolution::constrained::snapshot::DifferentialEvolutionRestoreLimits;

const MAGIC: &[u8; 8] = b"FSPSCP01";
const DOMAIN: &str = "fs-scan-fit/constrained-population-checkpoint-v1";
const HEADER: usize = 64;
const OVERHEAD: usize = HEADER + 32;

impl<F, O: ScanObservationSpec> PopulationScanStudy<F, O> {
    /// Capture a drained owner. Includes failed attempts, prior probes/replays,
    /// saved limits and successful parallel siblings. This is a byte transport,
    /// not a filesystem transaction, authentication or a global work ledger.
    /// The model identity must cover mesh, material, sources, observations,
    /// all constraint definitions/scales, time grid and numerical build.
    pub fn checkpoint_bytes(&self, cx: &Cx<'_>, model: [u8; 32], max_bytes: usize)
        -> Result<Vec<u8>, PopulationScanError>
    {
        check(cx)?;
        let limit = max_bytes.checked_sub(OVERHEAD)
            .ok_or_else(|| invalid("population checkpoint byte allowance is too small"))?;
        let search = self.search.checkpoint_bytes(cx, model, limit)?;
        let total = search.len().checked_add(OVERHEAD)
            .ok_or_else(|| invalid("population checkpoint size overflow"))?;
        let mut bytes = Vec::new();
        bytes.try_reserve_exact(total).map_err(|_| DifferentialEvolutionError::Allocation)?;
        bytes.extend_from_slice(MAGIC);
        bytes.extend_from_slice(&model);
        for n in [self.study.evaluations(), self.study.evaluation_limit(), search.len()] {
            bytes.extend_from_slice(&u64::try_from(n)
                .map_err(|_| invalid("population checkpoint count exceeds u64"))?.to_le_bytes());
        }
        check(cx)?;
        bytes.extend_from_slice(&search);
        let digest = fs_blake3::hash_domain(DOMAIN, &bytes);
        bytes.extend_from_slice(digest.as_bytes());
        check(cx)?;
        Ok(bytes)
    }

    /// Restore into an unevaluated, caller-reconstructed physical model. Saved
    /// search configuration, bounds, initial design, physical limit and policy
    /// must match the supplied declarations; changing a budget is an explicit
    /// extension AFTER restoration. All admission precedes counter publication.
    ///
    /// The identity is caller-owned: the codec cannot compare arbitrary physical
    /// closures. A previously used replay slot stays used; restoration may leave
    /// no room for replay or search until the allowances are explicitly extended.
    /// Work after the last capture cannot be recovered from an older snapshot.
    #[allow(clippy::too_many_arguments)]
    pub fn from_checkpoint_bytes(cx: &Cx<'_>, mut study: ValueConstrainedScanStudy<F, O>,
        bytes: &[u8], expected_model: [u8; 32], expected_bounds: &[[f64; 2]],
        expected_config: &DifferentialEvolutionConfig, expected_tolerance: f64,
        limits: DifferentialEvolutionRestoreLimits) -> Result<Self, PopulationScanError>
    {
        check(cx)?;
        if study.evaluations() != 0 || bytes.len() < OVERHEAD || bytes.len() > limits.max_bytes {
            return Err(invalid("population restore requires a fresh study and bounded complete bytes").into());
        }
        let end = bytes.len() - 32;
        if bytes.get(..8) != Some(MAGIC.as_slice())
            || bytes.get(8..40) != Some(expected_model.as_slice())
            || fs_blake3::hash_domain(DOMAIN, &bytes[..end]).as_bytes().as_slice() != &bytes[end..]
        { return Err(invalid("population checkpoint format, model or checksum mismatch").into()); }
        let spent = count(&bytes[40..48])?;
        let physical_limit = count(&bytes[48..56])?;
        let search_len = count(&bytes[56..64])?;
        if spent > physical_limit || physical_limit != study.evaluation_limit()
            || search_len != end - HEADER
        { return Err(invalid("population checkpoint physical allowance or payload length changed").into()); }
        let search = ConstrainedDifferentialEvolutionState::from_checkpoint_bytes(cx,
            &bytes[HEADER..end], expected_model, study.constraint_count(), expected_tolerance,
            DifferentialEvolutionRestoreLimits { max_bytes: limits.max_bytes - OVERHEAD, ..limits })?;
        let actual = search.configuration();
        let same_target = actual.target.map(f64::to_bits) == expected_config.target.map(f64::to_bits);
        if actual != expected_config || !same_target
            || actual.differential_weight.to_bits() != expected_config.differential_weight.to_bits()
            || actual.crossover_probability.to_bits() != expected_config.crossover_probability.to_bits()
            || search.bounds().len() != expected_bounds.len()
            || !search.bounds().iter().zip(expected_bounds).all(|(a,b)| same(a,b))
            || search.initial_coordinates().is_none_or(|x| !same(x, study.initial_coordinates()))
        { return Err(invalid("population restart changed the search configuration, box or initial design").into()); }
        let lower: Vec<_> = expected_bounds.iter().map(|b| b[0]).collect();
        let upper: Vec<_> = expected_bounds.iter().map(|b| b[1]).collect();
        study.sources_at(cx, &lower).map_err(DifferentialEvolutionRunError::Evaluation)?;
        study.sources_at(cx, &upper).map_err(DifferentialEvolutionRunError::Evaluation)?;
        check(cx)?;
        study.restore_evaluation_count(cx, spent).map_err(DifferentialEvolutionRunError::Evaluation)?;
        Ok(Self { study, search, target: expected_config.target })
    }
}

fn count(bytes: &[u8]) -> Result<usize, DifferentialEvolutionError> {
    let word: [u8; 8] = bytes.try_into().map_err(|_| invalid("truncated population count"))?;
    usize::try_from(u64::from_le_bytes(word)).map_err(|_| invalid("population count exceeds this platform"))
}
fn same(a: &[f64], b: &[f64]) -> bool {
    a.len() == b.len() && a.iter().zip(b).all(|(a,b)| a.to_bits() == b.to_bits())
}
