//! Durable snapshots of a system's accepted generation (bead rc-1109 .5.4,
//! item 4): written through fs-exec's [`SolverStateV2`] envelope, so that a
//! run resumes in another process and forks.
//!
//! The snapshot's state is one family owned by this crate, [`SystemImage`].
//! Its members are carried inside it as the plain words of each component's
//! own codec (`Component::encode_state`), because the envelope's codec has
//! no nested framing and its decoder cannot know a member's type. Which
//! system it belongs to is the snapshot's problem identity
//! ([`System::problem_id`]), from the graph and from each member's codec
//! charter and parameters.

use std::fmt::Write as _;

use fs_blake3::identity::CancellationProbe;
use fs_exec::solver::SolverStateV2;
use fs_exec::solver::snapshot_v2::{
    ExpectedResumeContextV2, OpenedSnapshotV2, SealedSnapshotV2, SnapshotDecoderV2,
    SnapshotEncoderV2, SnapshotLimitsV2, SnapshotProblemIdV2, SnapshotStateCodecIdV2,
    SnapshotStateSchemaIdV2, SnapshotStateTypeIdV2, SnapshotV2Error, StateIdentityCharterV2,
};

use super::{CheckpointError, SolveConfig, System};
use crate::ledger::Sum;

/// Domain of the problem identity of a system.
const PROBLEM_DOMAIN: &str = "org.frankensim.fs-system.problem.v1";

/// The charter of [`SystemImage`]. Its derived identities are pinned below
/// and checked against it by every seal and by a test.
const IMAGE_CHARTER: StateIdentityCharterV2 = StateIdentityCharterV2 {
    owner: "fs-system::System",
    state_family: "accepted-generation",
    schema_grammar: "problem: [u8; 32]; members: per component in index order, \
        {state type id [u8; 32], state schema id [u8; 32], state codec id [u8; 32], \
        codec version u32, words: [u64]}; accepted inputs: [f64]; accepted outputs: [f64]; \
        clock: sum and compensation as two u64 bits; steps: u64; config: relative_tolerance f64, \
        absolute_tolerance f64, max_iterations u32, max_shortenings u32, \
        interface_energy_tolerance f64; ledger: [u64]",
    codec_grammar: "v2 framed little-endian: [u8; 32] as four u64; u64 count of members; \
        [u64] as an f64 slice of the words' bits; [f64] as an f64 slice",
    codec_version: 1,
};

const IMAGE_TYPE_ID: SnapshotStateTypeIdV2 = SnapshotStateTypeIdV2::from_bytes(PINNED_TYPE);
const IMAGE_SCHEMA_ID: SnapshotStateSchemaIdV2 = SnapshotStateSchemaIdV2::from_bytes(PINNED_SCHEMA);
const IMAGE_CODEC_ID: SnapshotStateCodecIdV2 = SnapshotStateCodecIdV2::from_bytes(PINNED_CODEC);

// Derived from IMAGE_CHARTER (BLAKE3 derive-key under fs-exec's charter
// domains); `seal_v2` refuses if they drift, and so does a test.
const PINNED_TYPE: [u8; 32] = [
    0x48, 0xda, 0x3c, 0x2f, 0x41, 0x72, 0xca, 0x57, 0x0c, 0x03, 0xff, 0xf0, 0x27, 0x10, 0xdd, 0xf7,
    0x57, 0xba, 0x20, 0x74, 0xa0, 0x3c, 0x78, 0xd8, 0x26, 0x56, 0x92, 0xce, 0xb2, 0x2c, 0x76, 0xf3,
];
const PINNED_SCHEMA: [u8; 32] = [
    0xc6, 0xba, 0xbd, 0x96, 0x28, 0x33, 0x71, 0xb5, 0x5c, 0x6e, 0xbb, 0x51, 0xd1, 0xf0, 0x92, 0x77,
    0xf3, 0xbd, 0x70, 0x36, 0x68, 0x18, 0x67, 0x1f, 0xa3, 0xe5, 0xe6, 0x9d, 0x52, 0x2a, 0xb4, 0x99,
];
const PINNED_CODEC: [u8; 32] = [
    0xd1, 0xdb, 0x26, 0x47, 0xa7, 0x01, 0xd3, 0xf9, 0x51, 0xfa, 0x8e, 0x8a, 0xf0, 0x95, 0xa9, 0x2d,
    0x14, 0x04, 0xbd, 0xef, 0x94, 0x0a, 0xd9, 0x75, 0xea, 0x27, 0x75, 0x98, 0xa3, 0x85, 0xd4, 0xe3,
];

/// The accepted generation of a [`System`] as plain data: what a durable
/// snapshot of it holds.
///
/// It is made by [`System::snapshot`] and read by [`System::resume`]; its
/// contents are private. Each member's state is the words of that
/// component's own codec, with the identities of the codec's charter.
#[derive(Debug, Clone)]
pub struct SystemImage {
    problem: [u8; 32],
    members: Vec<MemberImage>,
    inputs: Vec<f64>,
    outputs: Vec<f64>,
    time: Sum,
    steps: u64,
    config: SolveConfig,
    ledger: Vec<u64>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct MemberImage {
    /// State type, schema and codec identities of the member's charter.
    charter: [[u8; 32]; 3],
    codec_version: u32,
    words: Vec<u64>,
}

impl MemberImage {
    fn charter_of(charter: &StateIdentityCharterV2) -> [[u8; 32]; 3] {
        [
            *charter.state_type_id().as_bytes(),
            *charter.state_schema_id().as_bytes(),
            *charter.state_codec_id().as_bytes(),
        ]
    }
}

fn put_id(encoder: &mut SnapshotEncoderV2<'_>, id: &[u8; 32]) -> Result<(), SnapshotV2Error> {
    for word in id.as_chunks::<8>().0 {
        encoder.put_u64(u64::from_le_bytes(*word))?;
    }
    Ok(())
}

fn get_id(decoder: &mut SnapshotDecoderV2<'_, '_>) -> Result<[u8; 32], SnapshotV2Error> {
    let mut id = [0_u8; 32];
    for word in id.as_chunks_mut::<8>().0 {
        *word = decoder.get_u64()?.to_le_bytes();
    }
    Ok(id)
}

/// Words travel as the bits of an `f64` slice, the envelope's one admitted
/// collection: the bits are moved, never computed with.
fn put_words(encoder: &mut SnapshotEncoderV2<'_>, words: &[u64]) -> Result<(), SnapshotV2Error> {
    let bits: Vec<f64> = words.iter().map(|&word| f64::from_bits(word)).collect();
    encoder.put_f64_slice(&bits)
}

fn get_words(decoder: &mut SnapshotDecoderV2<'_, '_>) -> Result<Vec<u64>, SnapshotV2Error> {
    Ok(decoder
        .get_f64_vec()?
        .into_iter()
        .map(f64::to_bits)
        .collect())
}

impl SolverStateV2 for SystemImage {
    const STATE_TYPE_ID_V2: SnapshotStateTypeIdV2 = IMAGE_TYPE_ID;
    const STATE_SCHEMA_ID_V2: SnapshotStateSchemaIdV2 = IMAGE_SCHEMA_ID;
    const STATE_CODEC_ID_V2: SnapshotStateCodecIdV2 = IMAGE_CODEC_ID;
    const STATE_CODEC_VERSION_V2: u32 = IMAGE_CHARTER.codec_version;

    fn charter() -> Option<&'static StateIdentityCharterV2> {
        Some(&IMAGE_CHARTER)
    }

    fn encode_v2(&self, encoder: &mut SnapshotEncoderV2<'_>) -> Result<(), SnapshotV2Error> {
        put_id(encoder, &self.problem)?;
        encoder.put_u64(self.members.len() as u64)?;
        for member in &self.members {
            for id in &member.charter {
                put_id(encoder, id)?;
            }
            encoder.put_u32(member.codec_version)?;
            put_words(encoder, &member.words)?;
        }
        encoder.put_f64_slice(&self.inputs)?;
        encoder.put_f64_slice(&self.outputs)?;
        for word in self.time.to_words() {
            encoder.put_u64(word)?;
        }
        encoder.put_u64(self.steps)?;
        encoder.put_f64(self.config.relative_tolerance)?;
        encoder.put_f64(self.config.absolute_tolerance)?;
        encoder.put_u32(self.config.max_iterations)?;
        encoder.put_u32(self.config.max_shortenings)?;
        encoder.put_f64(self.config.interface_energy_tolerance)?;
        put_words(encoder, &self.ledger)
    }

    fn decode_v2(decoder: &mut SnapshotDecoderV2<'_, '_>) -> Result<Self, SnapshotV2Error> {
        let problem = get_id(decoder)?;
        let count = decoder.get_u64()?;
        // The count is not trusted for an allocation: members are pushed as
        // they decode, and a count beyond the payload ends in its refusal.
        let mut members = Vec::new();
        for _ in 0..count {
            let charter = [get_id(decoder)?, get_id(decoder)?, get_id(decoder)?];
            members.push(MemberImage {
                charter,
                codec_version: decoder.get_u32()?,
                words: get_words(decoder)?,
            });
        }
        let inputs = decoder.get_f64_vec()?;
        let outputs = decoder.get_f64_vec()?;
        let time = Sum::from_words([decoder.get_u64()?, decoder.get_u64()?]);
        Ok(Self {
            problem,
            members,
            inputs,
            outputs,
            time,
            steps: decoder.get_u64()?,
            config: SolveConfig {
                relative_tolerance: decoder.get_f64()?,
                absolute_tolerance: decoder.get_f64()?,
                max_iterations: decoder.get_u32()?,
                max_shortenings: decoder.get_u32()?,
                interface_energy_tolerance: decoder.get_f64()?,
            },
            ledger: get_words(decoder)?,
        })
    }
}

/// One length-prefixed field of a problem identity's preimage.
fn field(preimage: &mut Vec<u8>, bytes: &[u8]) {
    preimage.extend_from_slice(&(bytes.len() as u64).to_le_bytes());
    preimage.extend_from_slice(bytes);
}

impl System {
    /// The identity of the problem this system is: its graph and what each
    /// member is.
    ///
    /// It hashes the description sealing made the system from (names,
    /// ports, classes and tiers, load lines, connections with their ratios,
    /// maps and heat routes, the driven ports and the heat sinks, as their
    /// `Debug` text) and, per member, the identities of its state codec's
    /// charter and the words of `Component::identity_words`. Boundary
    /// functions are not in it: a resumed run follows the same trajectory
    /// only under the same boundary responses. A toolchain that prints a
    /// description differently gives another identity, so it refuses older
    /// snapshots; it never makes two descriptions one.
    #[must_use]
    pub fn problem_id(&self) -> SnapshotProblemIdV2 {
        let mut preimage = Vec::new();
        let mut text = String::new();
        for count in [
            self.nodes.len(),
            self.connections.len(),
            self.boundaries.len(),
            self.ledger.sinks.len(),
        ] {
            field(&mut preimage, &(count as u64).to_le_bytes());
        }
        for (index, node) in self.nodes.iter().enumerate() {
            text.clear();
            let _ = write!(
                text,
                "{:?} {:?} {:?} {} {:?} {:?}",
                self.names[index],
                node.ports(),
                node.class(),
                node.has_tangent(),
                node.line_ports(),
                node.line_ends(),
            );
            field(&mut preimage, text.as_bytes());
            match node.state_charter() {
                None => field(&mut preimage, &[]),
                Some(charter) => {
                    for id in MemberImage::charter_of(&charter) {
                        field(&mut preimage, &id);
                    }
                    field(&mut preimage, &charter.codec_version.to_le_bytes());
                }
            }
            let mut words = Vec::new();
            node.identity_words(&mut words);
            let bytes: Vec<u8> = words.iter().flat_map(|word| word.to_le_bytes()).collect();
            field(&mut preimage, &bytes);
        }
        for connection in &self.connections {
            text.clear();
            let _ = write!(text, "{connection:?}");
            field(&mut preimage, text.as_bytes());
        }
        for (port, _) in &self.boundaries {
            text.clear();
            let _ = write!(text, "boundary {port:?}");
            field(&mut preimage, text.as_bytes());
        }
        for sink in &self.ledger.sinks {
            field(&mut preimage, sink.name.as_bytes());
        }
        SnapshotProblemIdV2::from_bytes(
            *fs_blake3::hash_domain(PROBLEM_DOMAIN, &preimage).as_bytes(),
        )
    }

    /// Seal the accepted generation as a durable snapshot under `context`.
    ///
    /// The caller makes `context` for [`SystemImage`] with this system's
    /// [`Self::problem_id`] and the pause boundary of its executor's drain.
    /// Capture also works after a refused or cancelled step: trial
    /// candidates and scratch values are not in it. The system is unchanged.
    ///
    /// # Errors
    /// A context naming another problem identity; a member without a state
    /// codec; or the envelope's refusal. Nothing is sealed then.
    pub fn snapshot<C: CancellationProbe>(
        &self,
        context: &ExpectedResumeContextV2,
        limits: SnapshotLimitsV2,
        cancellation: C,
    ) -> Result<SealedSnapshotV2, CheckpointError> {
        if context.context().problem() != self.problem_id() {
            return Err(CheckpointError::ForeignProblem);
        }
        self.image()?
            .seal_v2(context, limits, cancellation)
            .map_err(CheckpointError::Envelope)
    }

    fn image(&self) -> Result<SystemImage, CheckpointError> {
        let mut members = Vec::with_capacity(self.nodes.len());
        for (index, (node, state)) in self.nodes.iter().zip(&self.accepted).enumerate() {
            let no_codec = || CheckpointError::NoCodec {
                node: self.names[index].clone(),
            };
            let charter = node.state_charter().ok_or_else(no_codec)?;
            let mut words = Vec::new();
            let written = node
                .encode_state(state.as_ref(), &mut words)
                .map_err(|source| CheckpointError::State {
                    node: self.names[index].clone(),
                    source,
                })?;
            if !written {
                return Err(no_codec());
            }
            members.push(MemberImage {
                charter: MemberImage::charter_of(&charter),
                codec_version: charter.codec_version,
                words,
            });
        }
        let mut ledger = Vec::new();
        self.ledger.push_words(&mut ledger);
        Ok(SystemImage {
            problem: *self.problem_id().as_bytes(),
            members,
            inputs: self.accepted_inputs.clone(),
            outputs: self.accepted_outputs.clone(),
            time: self.time,
            steps: self.steps,
            config: self.config,
            ledger,
        })
    }

    /// Make an opened snapshot this system's accepted generation.
    ///
    /// The snapshot must be of this problem ([`Self::problem_id`]), so of a
    /// system sealed from the same description with members of the same
    /// codecs and parameters. Every member's state is decoded before
    /// anything changes, and a refusal changes nothing. As after
    /// [`Self::restore`], the allowance of the last step is zero until a
    /// step is taken. One opened snapshot may be resumed into several
    /// systems, which then go on independently.
    ///
    /// # Errors
    /// A snapshot of another problem; a member without a codec, or one that
    /// refuses its words; or an image that does not have this system's
    /// shape.
    pub fn resume(
        &mut self,
        opened: &OpenedSnapshotV2<SystemImage>,
    ) -> Result<(), CheckpointError> {
        let problem = self.problem_id();
        let image = opened.state();
        if opened.context().problem() != problem || image.problem != *problem.as_bytes() {
            return Err(CheckpointError::ForeignProblem);
        }
        if image.members.len() != self.nodes.len()
            || image.inputs.len() != self.accepted_inputs.len()
            || image.outputs.len() != self.accepted_outputs.len()
        {
            return Err(CheckpointError::Image);
        }
        let mut accepted = Vec::with_capacity(self.nodes.len());
        for (index, (node, member)) in self.nodes.iter().zip(&image.members).enumerate() {
            let charter = node
                .state_charter()
                .ok_or_else(|| CheckpointError::NoCodec {
                    node: self.names[index].clone(),
                })?;
            if member.charter != MemberImage::charter_of(&charter)
                || member.codec_version != charter.codec_version
            {
                return Err(CheckpointError::Image);
            }
            accepted.push(node.decode_state(&member.words).ok_or_else(|| {
                CheckpointError::Decode {
                    node: self.names[index].clone(),
                }
            })?);
        }
        let ledger = self
            .ledger
            .read(&image.ledger)
            .ok_or(CheckpointError::Image)?;

        self.accepted = accepted;
        self.accepted_inputs.copy_from_slice(&image.inputs);
        self.accepted_outputs.copy_from_slice(&image.outputs);
        self.accepted_allowance.fill(0.0);
        self.inputs.copy_from_slice(&self.accepted_inputs);
        self.outputs.copy_from_slice(&self.accepted_outputs);
        self.time = image.time;
        self.steps = image.steps;
        self.ledger = ledger;
        self.config = image.config;
        Ok(())
    }
}
