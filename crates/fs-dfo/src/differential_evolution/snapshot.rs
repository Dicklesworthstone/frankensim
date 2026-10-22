//! Restart a drained DE search without discarding its completed simulations.
//!
//! The byte format retains exact binary64 words, algorithm configuration, spent
//! attempts, complete parents, partial replacements and cached parallel outcomes.
//! A domain-separated checksum detects corruption, not malicious forgery. The
//! caller supplies the identity of the immutable objective/model/build; equal
//! IDs are an assertion by that caller, not proof that two callbacks are equal.
//! Resume requires the same numerical implementation and floating-point behavior.
//! No cross-ISA equivalence, durable filesystem write or global work ledger is
//! provided by this in-memory codec. Work after the last saved snapshot is lost.

use super::{DifferentialEvolutionConfig, DifferentialEvolutionError as Error,
    DifferentialEvolutionPoint as Point, DifferentialEvolutionState, Pending,
    DifferentialEvolutionScore, RankedDifferentialEvolutionState,
    DIFFERENTIAL_EVOLUTION_STREAM_KERNEL, checkpoint, finite, invalid, reserved};
use fs_exec::Cx;

const MAGIC: &[u8; 8] = b"FSDECP01";
/// Bump when changing stored fields or the candidate/selection semantics.
pub const DIFFERENTIAL_EVOLUTION_CHECKPOINT_VERSION: u32 = 1;
const DOMAIN: &str = "fs-dfo/differential-evolution-checkpoint-v1";
const DIGEST_BYTES: usize = 32;

/// Independent limits supplied by the restoring caller, before allocation.
/// These bound encoded bytes and stored shapes, not all allocator/solver memory.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DifferentialEvolutionRestoreLimits {
    pub max_bytes: usize,
    pub max_dimension: usize,
    /// Maximum population size times dimension, including an empty initial state.
    pub max_population_entries: usize,
}

// One state transport for scalar and feasibility-first searches. A codec owns
// only its score representation and immutable selection policy, never mutation,
// candidate generation, callback accounting or physical-model reconstruction.
pub(super) trait SnapshotCodec<S: DifferentialEvolutionScore> {
    const MAGIC: &'static [u8; 8];
    const VERSION: u32;
    const DOMAIN: &'static str;
    fn write_policy(&self, out: &mut Encoder) -> Result<(), Error>;
    fn read_policy(&self, input: &mut Decoder<'_>) -> Result<(), Error>;
    fn write_score(&self, out: &mut Encoder, score: S) -> Result<(), Error>;
    fn read_score(&self, input: &mut Decoder<'_>) -> Result<S, Error>;
    fn equivalent_best(&self, best: S, minimum: S) -> bool;
}

struct ScalarCodec;
impl SnapshotCodec<f64> for ScalarCodec {
    const MAGIC: &'static [u8; 8] = self::MAGIC;
    const VERSION: u32 = DIFFERENTIAL_EVOLUTION_CHECKPOINT_VERSION;
    const DOMAIN: &'static str = self::DOMAIN;
    fn write_policy(&self, _: &mut Encoder) -> Result<(), Error> { Ok(()) }
    fn read_policy(&self, _: &mut Decoder<'_>) -> Result<(), Error> { Ok(()) }
    fn write_score(&self, out: &mut Encoder, score: f64) -> Result<(), Error> { out.number(score) }
    fn read_score(&self, input: &mut Decoder<'_>) -> Result<f64, Error> { input.number() }
    // Scalar replacement admits equal numeric values including opposite signed
    // zeros; retain the historical snapshot validator's numeric comparison.
    fn equivalent_best(&self, best: f64, minimum: f64) -> bool { best == minimum }
}

impl DifferentialEvolutionState {
    /// Encode a quiescent search, including successful siblings of a failed wave.
    /// `&self` excludes an in-flight mutable advance. No callback is invoked and
    /// no work is spent. On cancellation/allocation/size refusal no bytes escape.
    /// Persist the returned bytes with an application-owned crash-safe protocol.
    /// The checksum is one indivisible operation, bounded by `max_bytes`.
    pub fn checkpoint_bytes(&self, cx: &Cx<'_>, model_identity: [u8; 32], max_bytes: usize)
        -> Result<Vec<u8>, Error>
    {
        self.encode_snapshot(model_identity, max_bytes, &mut || cx.checkpoint().is_ok())
    }

    /// Restore exact retained work without invoking the objective or recomputing
    /// scores. Unknown versions, mismatched model IDs, truncated/trailing bytes,
    /// checksum failures, invalid shapes/counters and foreign pending candidates
    /// refuse before publication. The supplied limits are not read from the file.
    /// A saved exhausted budget stays exhausted; extend it explicitly afterward.
    /// Use a fresh noncancelled Cx when saving state after a cancelled advance.
    pub fn from_checkpoint_bytes(cx: &Cx<'_>, bytes: &[u8], expected_model: [u8; 32],
        limits: DifferentialEvolutionRestoreLimits) -> Result<Self, Error>
    {
        Self::decode_snapshot(bytes, expected_model, limits, &mut || cx.checkpoint().is_ok())
    }

    fn encode_snapshot(&self, model: [u8; 32], max_bytes: usize,
        poll: &mut impl FnMut() -> bool) -> Result<Vec<u8>, Error>
    {
        self.encode_snapshot_with(model, max_bytes, &ScalarCodec, poll)
    }

    fn decode_snapshot(bytes: &[u8], model: [u8; 32], limits: DifferentialEvolutionRestoreLimits,
        poll: &mut impl FnMut() -> bool) -> Result<Self, Error>
    {
        Self::decode_snapshot_with(bytes, model, limits, &ScalarCodec, poll)
    }
}

impl<S: DifferentialEvolutionScore> RankedDifferentialEvolutionState<S> {
    pub(super) fn encode_snapshot_with<C: SnapshotCodec<S>>(&self, model: [u8; 32], max_bytes: usize,
        codec: &C, poll: &mut impl FnMut() -> bool) -> Result<Vec<u8>, Error>
    {
        checkpoint(poll)?;
        let mut out = Encoder { bytes: Vec::new(), limit: max_bytes.checked_sub(DIGEST_BYTES)
            .ok_or_else(|| invalid("DE checkpoint byte allowance is too small"))? };
        out.put(C::MAGIC)?;
        out.put(&C::VERSION.to_le_bytes())?;
        out.put(&fs_rand::STREAM_SEMANTICS_VERSION.to_le_bytes())?;
        out.put(&DIFFERENTIAL_EVOLUTION_STREAM_KERNEL.to_le_bytes())?;
        out.put(&model)?;
        codec.write_policy(&mut out)?;
        out.count(self.bounds.len())?;
        out.count(self.config.population_size)?;
        out.number(self.config.differential_weight)?;
        out.number(self.config.crossover_probability)?;
        out.put(&self.config.seed.to_le_bytes())?;
        out.count(self.config.max_evaluations)?;
        out.count(self.config.max_population_entries)?;
        out.optional(self.config.target)?;
        out.count(self.evaluations)?;
        out.count(self.generations)?;
        for (i, bound) in self.bounds.iter().enumerate() {
            if i % 256 == 0 { checkpoint(poll)?; }
            for &value in bound { out.number(value)?; }
        }
        out.flag(self.initial.is_some())?;
        if let Some(initial) = &self.initial { out.vector(initial, poll)?; }
        for points in [&self.population, &self.next] {
            out.count(points.len())?;
            for point in points { checkpoint(poll)?; out.point(point, codec, poll)?; }
        }
        out.flag(self.pending.is_some())?;
        if let Some(pending) = &self.pending { out.pending(pending, codec, poll)?; }
        out.count(self.queued.len())?;
        for pending in &self.queued { checkpoint(poll)?; out.pending(pending, codec, poll)?; }
        out.flag(self.best.is_some())?;
        if let Some(best) = &self.best { out.point(best, codec, poll)?; }
        checkpoint(poll)?;
        let digest = fs_blake3::hash_domain(C::DOMAIN, &out.bytes);
        out.limit = max_bytes;
        out.put(digest.as_bytes())?;
        checkpoint(poll)?;
        Ok(out.bytes)
    }

    pub(super) fn decode_snapshot_with<C: SnapshotCodec<S>>(bytes: &[u8], model: [u8; 32],
        limits: DifferentialEvolutionRestoreLimits, codec: &C,
        poll: &mut impl FnMut() -> bool) -> Result<Self, Error>
    {
        checkpoint(poll)?;
        if bytes.len() > limits.max_bytes || bytes.len() < C::MAGIC.len() + 12 + 32 + DIGEST_BYTES {
            return Err(invalid("DE checkpoint size exceeds its allowance or is truncated"));
        }
        let (payload, checksum) = bytes.split_at(bytes.len() - DIGEST_BYTES);
        if fs_blake3::hash_domain(C::DOMAIN, payload).as_bytes().as_slice() != checksum {
            return Err(invalid("DE checkpoint checksum mismatch"));
        }
        checkpoint(poll)?;
        let mut input = Decoder { bytes: payload, offset: 0 };
        if input.take(8)? != C::MAGIC || input.u32()? != C::VERSION
            || input.u32()? != fs_rand::STREAM_SEMANTICS_VERSION
            || input.u32()? != DIFFERENTIAL_EVOLUTION_STREAM_KERNEL
        { return Err(invalid("DE checkpoint format, algorithm or random-stream version is unsupported")); }
        if input.take(32)? != model.as_slice() { return Err(invalid("DE checkpoint belongs to a different objective/model")); }
        codec.read_policy(&mut input)?;
        let dimension = input.count()?;
        let population_size = input.count()?;
        if dimension == 0 || dimension > limits.max_dimension || population_size < 4
            || population_size > u32::MAX as usize
            || dimension.checked_mul(population_size).is_none_or(|n| n > limits.max_population_entries)
        { return Err(invalid("DE checkpoint population shape exceeds the restoring caller's allowance")); }
        let config = DifferentialEvolutionConfig {
            population_size, differential_weight: input.number()?, crossover_probability: input.number()?,
            seed: input.u64()?, max_evaluations: input.count()?, max_population_entries: input.count()?,
            target: input.optional()?,
        };
        let evaluations = input.count()?;
        let generations = input.count()?;
        input.need(dimension.checked_mul(16).ok_or_else(|| invalid("DE checkpoint shape overflow"))?)?;
        let mut bounds = reserved(dimension)?;
        for i in 0..dimension {
            if i % 256 == 0 { checkpoint(poll)?; }
            bounds.push([input.number()?, input.number()?]);
        }
        let initial = if input.flag()? { Some(input.vector(dimension, poll)?) } else { None };
        // Reuse the original algorithm's configuration, shape and box admission.
        let mut state = Self::new_with(&bounds, initial.as_deref(), config, poll)?;
        state.evaluations = evaluations;
        state.generations = generations;
        let count = input.count()?;
        if count > population_size { return Err(invalid("DE checkpoint parent count is invalid")); }
        for _ in 0..count { checkpoint(poll)?; state.population.push(input.point(dimension, codec, poll)?); }
        let count = input.count()?;
        if count >= population_size { return Err(invalid("DE checkpoint replacement count is invalid")); }
        for _ in 0..count { checkpoint(poll)?; state.next.push(input.point(dimension, codec, poll)?); }
        if input.flag()? { state.pending = Some(input.pending(dimension, codec, poll)?); }
        let count = input.count()?;
        if count > population_size { return Err(invalid("DE checkpoint queued count is invalid")); }
        state.queued.try_reserve_exact(count).map_err(|_| Error::Allocation)?;
        for _ in 0..count { checkpoint(poll)?; state.queued.push_back(input.pending(dimension, codec, poll)?); }
        if input.flag()? { state.best = Some(input.point(dimension, codec, poll)?); }
        if input.offset != payload.len() { return Err(invalid("DE checkpoint has trailing payload")); }
        state.validate_snapshot(codec, poll)?;
        checkpoint(poll)?;
        Ok(state)
    }

    fn validate_snapshot<C: SnapshotCodec<S>>(&self, codec: &C,
        poll: &mut impl FnMut() -> bool) -> Result<(), Error> {
        let size = self.config.population_size;
        let initializing = self.population.len() < size;
        if self.evaluations > self.config.max_evaluations
            || (initializing && (self.generations != 0 || !self.next.is_empty()))
            || self.population.is_empty() != self.best.is_none()
        { return Err(invalid("DE checkpoint counters or initialization state are inconsistent")); }
        let member = if initializing { self.population.len() } else { self.next.len() };
        let retained = self.queued.len().checked_add(usize::from(self.pending.is_some()))
            .ok_or_else(|| invalid("DE checkpoint pending count overflow"))?;
        if retained > size - member { return Err(invalid("DE checkpoint wave crosses its frozen generation")); }
        let consumed = if initializing { Some(member) } else {
            self.generations.checked_add(1).and_then(|g| g.checked_mul(size)).and_then(|n| n.checked_add(member))
        }.ok_or_else(|| invalid("DE checkpoint consumed-work count overflow"))?;
        let cached = self.pending.iter().chain(&self.queued).filter(|p| p.value.is_some()).count();
        if consumed.checked_add(cached).is_none_or(|n| n > self.evaluations) {
            return Err(invalid("DE checkpoint does not account for its completed callback attempts"));
        }
        for point in self.population.iter().chain(&self.next).chain(self.best.iter()) {
            checkpoint(poll)?;
            point.value.checked()?;
            for (i, (&x, &[lo, hi])) in point.coordinates.iter().zip(&self.bounds).enumerate() {
                if i % 256 == 0 { checkpoint(poll)?; }
                if !x.is_finite() || x < lo || x > hi { return Err(invalid("DE checkpoint point is outside its box")); }
            }
        }
        if let Some(best) = &self.best {
            let minimum = self.population.iter().chain(&self.next)
                .map(|p| p.value).min_by(|a, b| a.rank(*b))
                .ok_or_else(|| invalid("DE checkpoint lacks scored parents"))?;
            if !codec.equivalent_best(best.value, minimum) {
                return Err(invalid("DE checkpoint best score disagrees with retained population"));
            }
        }
        for (offset, pending) in self.pending.iter().chain(&self.queued).enumerate() {
            checkpoint(poll)?;
            let expected = self.candidate_at(member + offset, poll)?;
            if !same_words(&expected, &pending.coordinates) {
                return Err(invalid("DE checkpoint pending coordinates do not match its original candidate stream"));
            }
        }
        Ok(())
    }
}

fn same_words(a: &[f64], b: &[f64]) -> bool {
    a.len() == b.len() && a.iter().zip(b).all(|(a, b)| a.to_bits() == b.to_bits())
}
pub(super) struct Encoder { bytes: Vec<u8>, limit: usize }
impl Encoder {
    fn put(&mut self, value: &[u8]) -> Result<(), Error> {
        if self.bytes.len().checked_add(value.len()).is_none_or(|n| n > self.limit) {
            return Err(invalid("DE checkpoint exceeds its byte allowance"));
        }
        self.bytes.try_reserve(value.len()).map_err(|_| Error::Allocation)?;
        self.bytes.extend_from_slice(value);
        Ok(())
    }
    pub(super) fn count(&mut self, n: usize) -> Result<(), Error> {
        self.put(&u64::try_from(n).map_err(|_| invalid("DE checkpoint count exceeds u64"))?.to_le_bytes())
    }
    pub(super) fn number(&mut self, value: f64) -> Result<(), Error> { self.put(&finite(value)?.to_bits().to_le_bytes()) }
    fn flag(&mut self, value: bool) -> Result<(), Error> { self.put(&[u8::from(value)]) }
    fn optional(&mut self, value: Option<f64>) -> Result<(), Error> {
        self.flag(value.is_some())?;
        if let Some(value) = value { self.number(value)?; }
        Ok(())
    }
    fn vector(&mut self, values: &[f64], poll: &mut impl FnMut() -> bool) -> Result<(), Error> {
        for (i, &value) in values.iter().enumerate() {
            if i % 256 == 0 { checkpoint(poll)?; }
            self.number(value)?;
        }
        Ok(())
    }
    fn point<S: DifferentialEvolutionScore>(&mut self, point: &Point<S>,
        codec: &impl SnapshotCodec<S>, poll: &mut impl FnMut() -> bool) -> Result<(), Error> {
        self.vector(&point.coordinates, poll)?;
        codec.write_score(self, point.value)
    }
    fn pending<S: DifferentialEvolutionScore>(&mut self, pending: &Pending<S>,
        codec: &impl SnapshotCodec<S>, poll: &mut impl FnMut() -> bool) -> Result<(), Error> {
        self.vector(&pending.coordinates, poll)?;
        self.flag(pending.value.is_some())?;
        if let Some(score) = pending.value { codec.write_score(self, score)?; }
        Ok(())
    }
}
pub(super) struct Decoder<'a> { bytes: &'a [u8], offset: usize }
impl<'a> Decoder<'a> {
    fn need(&self, n: usize) -> Result<(), Error> {
        if n > self.bytes.len() - self.offset { Err(invalid("DE checkpoint is truncated")) } else { Ok(()) }
    }
    fn take(&mut self, n: usize) -> Result<&'a [u8], Error> {
        self.need(n)?;
        let start = self.offset;
        self.offset += n;
        Ok(&self.bytes[start..self.offset])
    }
    fn u32(&mut self) -> Result<u32, Error> {
        Ok(u32::from_le_bytes(self.take(4)?.try_into().expect("checked word")))
    }
    fn u64(&mut self) -> Result<u64, Error> {
        Ok(u64::from_le_bytes(self.take(8)?.try_into().expect("checked word")))
    }
    pub(super) fn count(&mut self) -> Result<usize, Error> {
        usize::try_from(self.u64()?).map_err(|_| invalid("DE checkpoint count exceeds this platform"))
    }
    pub(super) fn number(&mut self) -> Result<f64, Error> { finite(f64::from_bits(self.u64()?)) }
    fn flag(&mut self) -> Result<bool, Error> {
        match self.take(1)?[0] { 0 => Ok(false), 1 => Ok(true), _ => Err(invalid("DE checkpoint has an invalid optional-value tag")) }
    }
    fn optional(&mut self) -> Result<Option<f64>, Error> {
        if self.flag()? { Ok(Some(self.number()?)) } else { Ok(None) }
    }
    fn vector(&mut self, n: usize, poll: &mut impl FnMut() -> bool) -> Result<Vec<f64>, Error> {
        self.need(n.checked_mul(8).ok_or_else(|| invalid("DE checkpoint vector size overflow"))?)?;
        let mut values = reserved(n)?;
        for i in 0..n { if i % 256 == 0 { checkpoint(poll)?; } values.push(self.number()?); }
        Ok(values)
    }
    fn point<S: DifferentialEvolutionScore>(&mut self, n: usize,
        codec: &impl SnapshotCodec<S>, poll: &mut impl FnMut() -> bool) -> Result<Point<S>, Error> {
        Ok(Point { coordinates: self.vector(n, poll)?, value: codec.read_score(self)? })
    }
    fn pending<S: DifferentialEvolutionScore>(&mut self, n: usize,
        codec: &impl SnapshotCodec<S>, poll: &mut impl FnMut() -> bool) -> Result<Pending<S>, Error> {
        let coordinates = self.vector(n, poll)?;
        let value = if self.flag()? { Some(codec.read_score(self)?) } else { None };
        Ok(Pending { coordinates, value })
    }
}

#[cfg(test)]
mod tests;
