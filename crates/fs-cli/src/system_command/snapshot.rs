//! Durable accepted scenes, including their sampling and adaptive driver memory.

use std::io::Read as _;
use std::sync::OnceLock;

use fs_blake3::identity::{
    Admitted, AuthorityAdmitter, AuthorityRef, AuthorityVerifier, CanonicalLimits, CanonicalSchema,
    ContentId, ExternalAnchorRef, FieldSpec, KeyPolicyId, Presented, StrongIdentity, Verified,
    VerifierId, WireType,
};
use fs_exec::solver::SolverStateV2;
use fs_exec::solver::snapshot_v2::{
    self, ExpectedResumeContextV2, PausedSnapshotBoundaryV2, SnapshotAlgorithmIdV2,
    SnapshotAuthoritySubjectIdV2, SnapshotBudgetStateIdV2, SnapshotDecoderV2,
    SnapshotDeterminismV2, SnapshotEncoderV2, SnapshotExecutionFingerprintIdV2, SnapshotLimitsV2,
    SnapshotPauseRequestIdV2, SnapshotProblemIdV2, SnapshotProvenanceIdV2, SnapshotRngCounterIdV2,
    SnapshotStateCodecIdV2, SnapshotStateSchemaIdV2, SnapshotStateTypeIdV2, SnapshotV2Error,
    StateIdentityCharterV2,
};
use fs_exec::{CancelGate, DrainTracker, RunId};
use fs_system::SystemImage;

use super::{MAX_SAMPLES, MAX_STEPS, SystemScene, SystemSceneError, exit};

/// The file adapter applies this same cap before reading candidate bytes.
pub(super) const MAX_CHECKPOINT_BYTES: usize = 64 * 1024 * 1024;

const CHARTER: StateIdentityCharterV2 = StateIdentityCharterV2 {
    owner: "fs-cli::SystemScene",
    state_family: "accepted-scene",
    schema_grammar: "system: fs-system::SystemImage with its type, schema, codec and version; iterations u64; shortenings u64; next_sample_s f64; sample_width u64; sample_times [f64]; sample_values [f64] in report coordinate order; adaptive: optional attempts u64, trial_steps u64, last_error_norm f64, max_error_norm f64, suggested_dt_s f64",
    codec_grammar: "v2 little-endian primitives; system ids as four u64 each then u32 version and original SystemImage codec; samples use two f64 slices; adaptive uses a f64 slice of zero or five u64 bit words",
    codec_version: 1,
};

#[derive(Clone)]
struct SceneImage {
    system: SystemImage,
    iterations: u64,
    shortenings: u64,
    next_sample_s: f64,
    sample_width: u64,
    sample_times: Vec<f64>,
    sample_values: Vec<f64>,
    adaptive: Vec<f64>,
}

fn limits() -> SnapshotLimitsV2 {
    SnapshotLimitsV2::new(
        (MAX_CHECKPOINT_BYTES - snapshot_v2::HEADER_LEN_V2) as u64,
        65_536,
        CanonicalLimits::new(16_384, 4096, 32, 32, 64),
        (MAX_CHECKPOINT_BYTES / 8) as u64,
        MAX_CHECKPOINT_BYTES as u64,
        65_536,
    )
}

fn refused(message: impl std::fmt::Display) -> SystemSceneError {
    SystemSceneError {
        exit_code: exit::REFUSED,
        code: "system-checkpoint",
        message: message.to_string(),
    }
}

fn hash(domain: &str, bytes: &[u8]) -> [u8; 32] {
    *fs_blake3::hash_domain(domain, bytes).as_bytes()
}

fn executable_identity() -> Result<[u8; 32], SystemSceneError> {
    static IDENTITY: OnceLock<Result<[u8; 32], String>> = OnceLock::new();
    IDENTITY
        .get_or_init(|| {
            let path = std::env::current_exe().map_err(|error| error.to_string())?;
            let mut file = std::fs::File::open(path).map_err(|error| error.to_string())?;
            let mut identity =
                fs_blake3::DomainHasher::new("org.frankensim.fs-cli.system-scene.executable.v1");
            let mut buffer = [0; 65_536];
            loop {
                let count = file.read(&mut buffer).map_err(|error| error.to_string())?;
                if count == 0 {
                    break;
                }
                identity.update(&buffer[..count]);
            }
            Ok(*identity.finalize().as_bytes())
        })
        .clone()
        .map_err(|error| refused(format!("cannot bind checkpoint executable: {error}")))
}

/// This synchronous owner has exclusive access between driver calls. The
/// short-lived worker records that quiescent serialization operation; no
/// active numerical worker or caller-owned cancellation gate is impersonated.
fn context<S: SolverStateV2>(
    problem: SnapshotProblemIdV2,
) -> Result<ExpectedResumeContextV2, SystemSceneError> {
    let gate = CancelGate::new_clock_free();
    let run = RunId(u64::from_le_bytes(
        problem.as_bytes()[..8].try_into().unwrap(),
    ));
    let tracker = DrainTracker::new(run, &gate);
    let worker = tracker.register_worker().map_err(refused)?;
    gate.request();
    drop(worker);
    let pause = PausedSnapshotBoundaryV2::from_drain_report(
        tracker.finalize().map_err(refused)?,
        SnapshotPauseRequestIdV2::from_bytes(hash(
            "org.frankensim.fs-cli.system-scene.pause.v1",
            problem.as_bytes(),
        )),
        1,
    );
    Ok(ExpectedResumeContextV2::for_paused_state::<S>(
        SnapshotAlgorithmIdV2::from_bytes(hash(
            "org.frankensim.fs-cli.system-scene.algorithm.v1",
            b"accepted-scene",
        )),
        1,
        problem,
        SnapshotRngCounterIdV2::from_bytes([0; 32]),
        SnapshotDeterminismV2::Deterministic,
        SnapshotExecutionFingerprintIdV2::from_bytes(executable_identity()?),
        SnapshotBudgetStateIdV2::from_bytes(hash(
            "org.frankensim.fs-cli.system-scene.limits.v1",
            &(MAX_CHECKPOINT_BYTES as u64).to_le_bytes(),
        )),
        SnapshotProvenanceIdV2::from_bytes(*problem.as_bytes()),
        pause,
    ))
}

impl SolverStateV2 for SceneImage {
    const STATE_TYPE_ID_V2: SnapshotStateTypeIdV2 = SnapshotStateTypeIdV2::from_bytes([
        0x31, 0xa9, 0xc2, 0x46, 0x72, 0x42, 0x8a, 0x67, 0x64, 0x65, 0xa6, 0x06, 0x9b, 0xef, 0xbe,
        0x69, 0x44, 0x25, 0x53, 0xae, 0x5c, 0x4d, 0x8f, 0x2f, 0x20, 0xf8, 0xfc, 0x4e, 0x9c, 0x8b,
        0x77, 0xb9,
    ]);
    const STATE_SCHEMA_ID_V2: SnapshotStateSchemaIdV2 = SnapshotStateSchemaIdV2::from_bytes([
        0xfc, 0xa7, 0xc8, 0xd6, 0x69, 0xdd, 0x29, 0xde, 0xc7, 0x74, 0xdf, 0x93, 0x44, 0x58, 0x1f,
        0x5f, 0x26, 0x55, 0x31, 0xae, 0x91, 0x75, 0x88, 0x12, 0x1f, 0x5b, 0xfd, 0x8a, 0xda, 0x17,
        0x11, 0x4f,
    ]);
    const STATE_CODEC_ID_V2: SnapshotStateCodecIdV2 = SnapshotStateCodecIdV2::from_bytes([
        0xd0, 0xb8, 0x77, 0x6b, 0x70, 0x79, 0xaa, 0x03, 0xea, 0x29, 0xb5, 0xd9, 0x20, 0x89, 0x2c,
        0xec, 0x11, 0x5e, 0x99, 0xf6, 0xf5, 0x56, 0xa9, 0xa2, 0x70, 0xa0, 0x03, 0x31, 0x84, 0xd2,
        0xe7, 0xd0,
    ]);
    const STATE_CODEC_VERSION_V2: u32 = 1;

    fn charter() -> Option<&'static StateIdentityCharterV2> {
        Some(&CHARTER)
    }

    fn encode_v2(&self, encoder: &mut SnapshotEncoderV2<'_>) -> Result<(), SnapshotV2Error> {
        for id in [
            SystemImage::STATE_TYPE_ID_V2.as_bytes(),
            SystemImage::STATE_SCHEMA_ID_V2.as_bytes(),
            SystemImage::STATE_CODEC_ID_V2.as_bytes(),
        ] {
            for word in id.as_chunks::<8>().0 {
                encoder.put_u64(u64::from_le_bytes(*word))?;
            }
        }
        encoder.put_u32(SystemImage::STATE_CODEC_VERSION_V2)?;
        self.system.encode_v2(encoder)?;
        encoder.put_u64(self.iterations)?;
        encoder.put_u64(self.shortenings)?;
        encoder.put_f64(self.next_sample_s)?;
        encoder.put_u64(self.sample_width)?;
        encoder.put_f64_slice(&self.sample_times)?;
        encoder.put_f64_slice(&self.sample_values)?;
        encoder.put_f64_slice(&self.adaptive)
    }

    fn decode_v2(decoder: &mut SnapshotDecoderV2<'_, '_>) -> Result<Self, SnapshotV2Error> {
        for id in [
            SystemImage::STATE_TYPE_ID_V2.as_bytes(),
            SystemImage::STATE_SCHEMA_ID_V2.as_bytes(),
            SystemImage::STATE_CODEC_ID_V2.as_bytes(),
        ] {
            for word in id.as_chunks::<8>().0 {
                if decoder.get_u64()? != u64::from_le_bytes(*word) {
                    return Err(SnapshotV2Error::WrongStateSchema);
                }
            }
        }
        if decoder.get_u32()? != SystemImage::STATE_CODEC_VERSION_V2 {
            return Err(SnapshotV2Error::WrongStateSchema);
        }
        Ok(Self {
            system: SystemImage::decode_v2(decoder)?,
            iterations: decoder.get_u64()?,
            shortenings: decoder.get_u64()?,
            next_sample_s: decoder.get_f64()?,
            sample_width: decoder.get_u64()?,
            sample_times: decoder.get_f64_vec()?,
            sample_values: decoder.get_f64_vec()?,
            adaptive: decoder.get_f64_vec()?,
        })
    }
}

// Authority is the root the caller retained at sealing, supplied separately
// from candidate bytes. Inspection alone cannot manufacture that authority.
enum KeptRootSchema {}
impl CanonicalSchema for KeptRootSchema {
    const DOMAIN: &'static str = "org.frankensim.fs-cli.system-scene.kept-root.v1";
    const NAME: &'static str = "scene-kept-root";
    const VERSION: u32 = 1;
    const CONTEXT: &'static str = "compare with the caller's independently retained seal root";
    const FIELDS: &'static [FieldSpec] = &[FieldSpec::required("version", WireType::U64)];
}
enum ResumePolicySchema {}
impl CanonicalSchema for ResumePolicySchema {
    const DOMAIN: &'static str = "org.frankensim.fs-cli.system-scene.resume-policy.v1";
    const NAME: &'static str = "scene-resume-policy";
    const VERSION: u32 = 1;
    const CONTEXT: &'static str = "resume only the exact caller-pinned scene checkpoint";
    const FIELDS: &'static [FieldSpec] = &[FieldSpec::required("version", WireType::U64)];
}
type Subject<S> = AuthorityRef<SnapshotAuthoritySubjectIdV2, KeptRootSchema, ResumePolicySchema, S>;
struct KeptRoot([u8; 32]);
impl AuthorityVerifier<SnapshotAuthoritySubjectIdV2, KeptRootSchema, ResumePolicySchema>
    for KeptRoot
{
    type Error = &'static str;
    fn verify(&self, candidate: &Subject<Presented>) -> Result<(), Self::Error> {
        if candidate.receipt().id().as_bytes() == &self.0 {
            Ok(())
        } else {
            Err("checkpoint does not match the independently retained resume root")
        }
    }
}
struct ResumePolicy;
impl AuthorityAdmitter<SnapshotAuthoritySubjectIdV2, KeptRootSchema, ResumePolicySchema>
    for ResumePolicy
{
    type Error = &'static str;
    fn admit(&self, _: &Subject<Verified>) -> Result<(), Self::Error> {
        Ok(())
    }
}

fn authority(bytes: &[u8], root: [u8; 32]) -> Result<Subject<Admitted>, SystemSceneError> {
    let receipt = snapshot_v2::inspect(bytes, limits(), || false)
        .map_err(refused)?
        .authority_subject_receipt();
    Subject::<Presented>::present(
        receipt,
        ExternalAnchorRef::presented(ContentId::of_bytes(&root)),
        VerifierId::<KeptRootSchema>::parse_slice(&hash(KeptRootSchema::DOMAIN, b"version 1"))
            .ok_or_else(|| refused("checkpoint verifier identity is not 32 bytes"))?,
        KeyPolicyId::<ResumePolicySchema>::parse_slice(&hash(
            ResumePolicySchema::DOMAIN,
            b"version 1",
        ))
        .ok_or_else(|| refused("checkpoint resume policy identity is not 32 bytes"))?,
    )
    .verify(&KeptRoot(root))
    .map_err(refused)?
    .admit(&ResumePolicy)
    .map_err(refused)
}

impl SystemScene {
    fn checkpoint_problem(&self) -> SnapshotProblemIdV2 {
        let mut identity =
            fs_blake3::DomainHasher::new("org.frankensim.fs-cli.system-scene.problem.v1");
        identity.update(&self.source_identity);
        identity.update(self.device.system.problem_id().as_bytes());
        for (name, material) in &self.device.materials {
            for text in [name.as_str(), material.report_json()] {
                identity.update(&(text.len() as u64).to_le_bytes());
                identity.update(text.as_bytes());
            }
        }
        SnapshotProblemIdV2::from_bytes(*identity.finalize().as_bytes())
    }

    /// Capture a durable accepted prefix and the 32-byte root needed to resume it.
    ///
    /// Keep the returned root independently from the candidate checkpoint file.
    /// The checkpoint uses the existing v2 envelope and native component codecs;
    /// an unsupported component is refused by name. Both the file and restored
    /// sample allocation are capped at 64 MiB. Capturing does not change the scene.
    pub fn checkpoint(&self) -> Result<(Vec<u8>, [u8; 32]), SystemSceneError> {
        let image = SceneImage::capture(self)?;
        let sealed = image
            .seal_v2(
                &context::<SceneImage>(self.checkpoint_problem())?,
                limits(),
                || false,
            )
            .map_err(refused)?;
        let root = *sealed.authority_subject_receipt().id().as_bytes();
        let (bytes, _evidence) = sealed.into_parts();
        Ok((bytes, root))
    }

    /// Resume a checkpoint under its independently retained root.
    ///
    /// The exact admitted source bytes, resolved material identities, component
    /// graph, and executable must match. The source label may differ. Every
    /// sample and adaptive-memory field is admitted before the kernel's atomic
    /// restore; a refusal changes no state, ledger, samples or driver counters.
    pub fn resume_checkpoint(
        &mut self,
        bytes: &[u8],
        trusted_root: [u8; 32],
    ) -> Result<(), SystemSceneError> {
        if bytes.len() > MAX_CHECKPOINT_BYTES {
            return Err(refused("checkpoint exceeds 64 MiB"));
        }
        let opened = SceneImage::unseal_v2_authorized(
            bytes,
            &authority(bytes, trusted_root)?,
            &context::<SceneImage>(self.checkpoint_problem())?,
            limits(),
            || false,
        )
        .map_err(refused)?;
        let image = opened.state();
        let keys = image.validate(self)?;
        let samples: Vec<_> = image
            .sample_times
            .iter()
            .enumerate()
            .map(|(row, &time)| {
                let values = keys
                    .iter()
                    .enumerate()
                    .map(|(column, key)| {
                        (key.clone(), image.sample_values[row * keys.len() + column])
                    })
                    .collect();
                (time, values)
            })
            .collect();

        // The enclosing envelope has already admitted this exact image. Reuse
        // the System owner's normal seal/open/restore path instead of exposing
        // a raw-image activation bypass or implementing another System codec.
        let system_context = context::<SystemImage>(self.device.system.problem_id())?;
        let system_sealed = image
            .system
            .seal_v2(&system_context, limits(), || false)
            .map_err(refused)?;
        let system_opened = SystemImage::unseal_v2_expected(
            system_sealed.bytes(),
            &system_sealed.expectation(),
            limits(),
            || false,
        )
        .map_err(refused)?;
        self.device.system.resume(&system_opened).map_err(refused)?;
        self.device.samples = samples;
        self.device.next_sample_s = image.next_sample_s;
        self.device.iterations = image.iterations;
        self.device.shortenings = image.shortenings;
        if let Some(adaptive) = &mut self.device.adaptive {
            adaptive.attempts = image.adaptive[0].to_bits();
            adaptive.trial_steps = image.adaptive[1].to_bits();
            adaptive.last_error_norm = image.adaptive[2];
            adaptive.max_error_norm = image.adaptive[3];
            adaptive.suggested_dt_s = image.adaptive[4];
        }
        Ok(())
    }
}

/// Bound the actual restored Vec/String allocation, not only wire floats.
fn sample_storage(keys: &[String], rows: usize) -> Result<(), SystemSceneError> {
    let row_bytes = keys
        .iter()
        .try_fold(
            std::mem::size_of::<(f64, Vec<(String, f64)>)>(),
            |bytes, key| {
                bytes
                    .checked_add(std::mem::size_of::<(String, f64)>())?
                    .checked_add(key.len())
            },
        )
        .ok_or_else(|| refused("checkpoint sample allocation overflow"))?;
    if row_bytes
        .checked_mul(rows)
        .is_none_or(|bytes| bytes > MAX_CHECKPOINT_BYTES)
    {
        return Err(refused("checkpoint restored samples exceed 64 MiB"));
    }
    Ok(())
}

impl SceneImage {
    fn capture(scene: &SystemScene) -> Result<Self, SystemSceneError> {
        let keys: Vec<_> = scene
            .device
            .samples
            .first()
            .into_iter()
            .flat_map(|(_, row)| row.iter().map(|(key, _)| key.clone()))
            .collect();
        sample_storage(&keys, scene.device.samples.len())?;
        let system = &scene.device.system;
        let sealed = system
            .snapshot(
                &context::<SystemImage>(system.problem_id())?,
                limits(),
                || false,
            )
            .map_err(refused)?;
        let opened = SystemImage::unseal_v2_expected(
            sealed.bytes(),
            &sealed.expectation(),
            limits(),
            || false,
        )
        .map_err(refused)?;
        let image = Self {
            system: opened.state().clone(),
            iterations: scene.device.iterations,
            shortenings: scene.device.shortenings,
            next_sample_s: scene.device.next_sample_s,
            sample_width: keys.len() as u64,
            sample_times: scene.device.samples.iter().map(|(time, _)| *time).collect(),
            sample_values: scene
                .device
                .samples
                .iter()
                .flat_map(|(_, row)| row.iter().map(|(_, value)| *value))
                .collect(),
            adaptive: scene.device.adaptive.as_ref().map_or_else(Vec::new, |a| {
                vec![
                    f64::from_bits(a.attempts),
                    f64::from_bits(a.trial_steps),
                    a.last_error_norm,
                    a.max_error_norm,
                    a.suggested_dt_s,
                ]
            }),
        };
        image.validate(scene)?;
        Ok(image)
    }

    fn validate(&self, scene: &SystemScene) -> Result<Vec<String>, SystemSceneError> {
        let device = &scene.device;
        let time = self.system.time_s();
        let steps = self.system.steps();
        let slack = 1.0e-12 * device.duration_s;
        if !time.is_finite()
            || time < 0.0
            || time - device.duration_s > slack
            || steps > 2 * MAX_STEPS as u64
            || (steps == 0) != (time == 0.0)
        {
            return Err(refused(
                "checkpoint clock or step count is outside the scene horizon",
            ));
        }
        let config = device.system.checkpoint().map_err(refused)?.config();
        if self.system.config() != config {
            return Err(refused(
                "checkpoint solver configuration differs from the scene",
            ));
        }

        let trial_steps = match (&device.adaptive, self.adaptive.as_slice()) {
            (None, []) => steps,
            (Some(settings), [attempts, trials, last, maximum, suggested]) => {
                let attempts = attempts.to_bits();
                let trials = trials.to_bits();
                let intervals = steps / 2;
                if steps % 2 != 0
                    || attempts < intervals
                    || attempts > intervals.saturating_mul(u64::from(settings.config.max_attempts))
                    || trials < 3 * intervals
                    || trials > attempts.saturating_mul(3)
                    || !suggested.is_finite()
                    || *suggested <= 0.0
                    || *suggested > 2.0 * device.step_s
                    || !maximum.is_finite()
                    || !(0.0..=1.0).contains(maximum)
                    || (steps == 0
                        && (!last.is_nan()
                            || *maximum != 0.0
                            || suggested.to_bits() != device.step_s.to_bits()))
                    || (steps != 0 && (!last.is_finite() || !(0.0..=*maximum).contains(last)))
                {
                    return Err(refused(
                        "checkpoint adaptive memory is inconsistent with the accepted steps",
                    ));
                }
                trials
            }
            _ => return Err(refused("checkpoint adaptive state has the wrong shape")),
        };
        let iteration_limit = trial_steps
            .saturating_mul(u64::from(config.max_iterations))
            .saturating_mul(device.elements.len() as u64);
        if self.iterations > iteration_limit
            || self.shortenings > trial_steps.saturating_mul(u64::from(config.max_shortenings))
        {
            return Err(refused(
                "checkpoint driver work counters exceed the admitted step budget",
            ));
        }

        let keys: Vec<_> = device
            .samples
            .first()
            .into_iter()
            .flat_map(|(_, row)| row.iter().map(|(key, _)| key.clone()))
            .collect();
        let rows = self.sample_times.len();
        if rows > MAX_SAMPLES + 1
            || self.sample_width != keys.len() as u64
            || rows.checked_mul(keys.len()) != Some(self.sample_values.len())
            || self.sample_values.iter().any(|value| value.is_infinite())
        {
            return Err(refused(
                "checkpoint sample dimensions or values are invalid",
            ));
        }
        sample_storage(&keys, rows)?;
        if let Some(every) = device.sample_every_s {
            if self
                .sample_times
                .first()
                .is_none_or(|time| time.to_bits() != 0.0_f64.to_bits())
            {
                return Err(refused("checkpoint is missing the initial sample"));
            }
            if !self
                .sample_values
                .iter()
                .take(keys.len())
                .zip(&device.samples[0].1)
                .all(|(value, (_, initial))| value.to_bits() == initial.to_bits())
            {
                return Err(refused("checkpoint initial sample differs from the scene"));
            }
            let mut next = every;
            let mut previous = 0.0;
            let mut scheduled = 0;
            for (index, &sample_time) in self.sample_times.iter().enumerate().skip(1) {
                if !sample_time.is_finite()
                    || sample_time <= previous
                    || sample_time > time
                    || (sample_time < next - slack
                        && !(index + 1 == rows && device.duration_s - sample_time <= slack))
                {
                    return Err(refused(
                        "checkpoint sample times are outside the sampling schedule",
                    ));
                }
                while next <= sample_time + slack {
                    scheduled += 1;
                    if scheduled > MAX_SAMPLES + 1 || next + every <= next {
                        return Err(refused(
                            "checkpoint sample schedule exceeds its admitted budget",
                        ));
                    }
                    next += every;
                }
                previous = sample_time;
            }
            if next.to_bits() != self.next_sample_s.to_bits()
                || (time >= next - slack && device.duration_s - time > slack)
                || (device.duration_s - time <= slack
                    && self.sample_times.last().copied() != Some(time))
            {
                return Err(refused(
                    "checkpoint next sample time is inconsistent with its retained samples",
                ));
            }
        } else if rows != 0 || self.next_sample_s != f64::INFINITY {
            return Err(refused(
                "checkpoint contains sampling state for an unsampled scene",
            ));
        }
        Ok(keys)
    }
}

#[cfg(test)]
mod tests;
