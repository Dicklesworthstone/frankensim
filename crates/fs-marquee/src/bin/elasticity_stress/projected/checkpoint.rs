//! Bounded, bit-exact checkpoints for the single-load projected stress command.
//! The digest detects corruption, not authorship. Restoration re-solves the
//! saved field and refuses changed mechanics before any new study is published.
use super::*;
use fs_topols::OptimizeCheckpoint;
use fs_topols::projected_stress::ProjectedStressSetupStage;
use std::fs::File;
use std::io::{self, Read};

const MAGIC: &[u8] = b"fs-marquee-projected-stress-checkpoint-v1\n";
const MAX_BYTES: u64 = 1024 * 1024;

#[derive(Clone)]
pub(super) struct Policy {
    pub fixed: Vec<(usize, f64)>,
    pub projection: VolumeProjectionSettings,
    pub controls: ProjectedSettings,
}

#[derive(Default)]
pub(super) struct Options {
    pub enabled: bool,
    pub pause_after: Option<usize>,
}

pub(super) fn options(args: &[String]) -> Result<(Vec<String>, Options), Box<dyn Error>> {
    let mut options = Options::default();
    let mut checkpoint_seen = false;
    let mut positional = Vec::new();
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--checkpoint" if !checkpoint_seen => {
                checkpoint_seen = true;
                options.enabled = true;
            }
            "--pause-after" if options.pause_after.is_none() => {
                i += 1;
                let count = args.get(i).ok_or("--pause-after requires an update count")?.parse::<usize>()?;
                if count > 200 { return Err("--pause-after must be in 0..=200".into()); }
                options.pause_after = Some(count);
                options.enabled = true;
            }
            flag if flag.starts_with("--") => return Err("unknown or repeated checkpoint option".into()),
            _ => positional.push(args[i].clone()),
        }
        i += 1;
    }
    Ok((positional, options))
}

// A fixed-block, domain-separated digest chain avoids loading the executable
// into memory or introducing another hashing dependency. The read loop assembles
// full blocks even when the underlying reader returns short reads.
fn executable_digest(mut input: impl Read) -> io::Result<String> {
    let mut digest = fs_ledger::hash_bytes(b"fs-marquee-stress-executable-chain-v1").to_string();
    let mut block = vec![0u8; 64 * 1024];
    let mut total = 0u64;
    loop {
        let mut filled = 0;
        while filled < block.len() {
            match input.read(&mut block[filled..]) {
                Ok(0) => break,
                Ok(count) => filled += count,
                Err(error) if error.kind() == io::ErrorKind::Interrupted => continue,
                Err(error) => return Err(error),
            }
        }
        if filled == 0 { break; }
        total += filled as u64;
        if total > 1024 * 1024 * 1024 {
            return Err(io::Error::other("executable exceeds fingerprint byte budget"));
        }
        let mut preimage = digest.into_bytes();
        preimage.extend_from_slice(&(filled as u64).to_le_bytes());
        preimage.extend_from_slice(&block[..filled]);
        digest = fs_ledger::hash_bytes(&preimage).to_string();
    }
    if total == 0 { return Err(io::Error::other("empty executable")); }
    Ok(digest)
}

pub(super) fn executable() -> Result<String, Box<dyn Error>> {
    Ok(executable_digest(File::open(std::env::current_exe()?)?)?)
}

fn word(out: &mut Vec<u8>, value: u64) { out.extend_from_slice(&value.to_le_bytes()); }
fn real(out: &mut Vec<u8>, value: f64) { word(out, value.to_bits()); }

fn encode(optimizer: &ProjectedStressOptimizer, policy: &Policy, executable: &str) -> Vec<u8> {
    let checkpoint = optimizer.checkpoint();
    let settings = checkpoint.settings();
    let fixture = checkpoint.fixture();
    let projection = policy.projection;
    let controls = policy.controls;
    let limit = optimizer.limit();
    let state = optimizer.current();
    let mut out = MAGIC.to_vec();
    out.extend_from_slice(executable.as_bytes());
    for value in [u64::from(settings.level), settings.iterations as u64,
        checkpoint.next_iteration() as u64, settings.nucleation_period as u64,
        projection.max_evaluations as u64, controls.max_candidates as u64,
        controls.poll_iters as u64, policy.fixed.len() as u64,
        checkpoint.geometry().nodes().len() as u64]
    { word(&mut out, value); }
    for value in [settings.volfrac, settings.band_cells, settings.move_cells,
        settings.ell0, settings.mu_al, settings.sobolev_alpha, settings.hole_radius_cells,
        settings.youngs, settings.poisson, fixture.load, fixture.band,
        projection.target, projection.tolerance, projection.max_shift,
        controls.contraction, controls.min_relative_improvement,
        limit.max_von_mises, limit.absolute_tolerance, checkpoint.ell(),
        state.compliance, state.volume, state.sampled_max_von_mises,
        state.max_location[0], state.max_location[1]]
    { real(&mut out, value); }
    word(&mut out, state.sample_count as u64);
    word(&mut out, state.snapshot);
    for &(index, value) in &policy.fixed {
        word(&mut out, index as u64);
        real(&mut out, value);
    }
    for &value in checkpoint.geometry().nodes() { real(&mut out, value); }
    out
}

pub(super) fn save(
    path: &Path, optimizer: &ProjectedStressOptimizer, policy: &Policy, executable: &str,
) -> Result<(), Box<dyn Error>> {
    if executable.len() != 64 { return Err("invalid executable fingerprint".into()); }
    let bytes = encode(optimizer, policy, executable);
    if bytes.len() as u64 + 65 > MAX_BYTES { return Err("checkpoint byte budget exceeded".into()); }
    let mut file = writer(path)?;
    writeln!(file, "{}", fs_ledger::hash_bytes(&bytes))?;
    file.write_all(&bytes)?;
    file.flush()?;
    file.get_ref().sync_all()?;
    Ok(())
}

struct Reader<'a>(&'a [u8]);
impl Reader<'_> {
    fn word(&mut self) -> Result<u64, Box<dyn Error>> {
        let bytes: [u8; 8] = self.0.get(..8).ok_or("truncated checkpoint")?.try_into()?;
        self.0 = &self.0[8..];
        Ok(u64::from_le_bytes(bytes))
    }
    fn count(&mut self, max: usize) -> Result<usize, Box<dyn Error>> {
        let value = usize::try_from(self.word()?)?;
        if value > max { return Err("checkpoint count exceeds admitted budget".into()); }
        Ok(value)
    }
    fn real(&mut self) -> Result<f64, Box<dyn Error>> {
        let value = f64::from_bits(self.word()?);
        if !value.is_finite() { return Err("checkpoint contains non-finite data".into()); }
        Ok(value)
    }
}

struct Saved {
    checkpoint: OptimizeCheckpoint,
    policy: Policy,
    limit: SampledStressLimit,
    expected: SampledStressEvaluation,
}

fn decode(bytes: &[u8], executable: &str) -> Result<Saved, Box<dyn Error>> {
    if bytes.len() as u64 > MAX_BYTES { return Err("checkpoint byte budget exceeded".into()); }
    let bytes = bytes.strip_prefix(MAGIC).ok_or("unsupported stress checkpoint schema")?;
    if bytes.get(..64) != Some(executable.as_bytes()) {
        return Err("checkpoint requires the original executable; cross-build resume refused".into());
    }
    let mut reader = Reader(&bytes[64..]);
    let level = reader.count(7)? as u32;
    if level < 2 { return Err("checkpoint level must be in 2..=7".into()); }
    let n = 1usize << level;
    let iterations = reader.count(200)?;
    let ordinal = reader.count(iterations)?;
    let nucleation_period = reader.count(200)?;
    let max_evaluations = reader.count(64)?;
    let max_candidates = reader.count(16)?;
    let poll_iters = reader.count(60_000)?;
    let fixed_count = reader.count(2 * (n + 1))?;
    let node_count = reader.count((n + 1) * (n + 1))?;
    if fixed_count != 2 * (n + 1) || node_count != (n + 1) * (n + 1) {
        return Err("checkpoint must retain the complete lattice and both fixed boundary traces".into());
    }
    let settings = OptimizeSettings {
        level, iterations, nucleation_period,
        volfrac: reader.real()?, band_cells: reader.real()?, move_cells: reader.real()?,
        ell0: reader.real()?, mu_al: reader.real()?, sobolev_alpha: reader.real()?,
        hole_radius_cells: reader.real()?, youngs: reader.real()?, poisson: reader.real()?,
    };
    let fixture = Cantilever { load: reader.real()?, band: reader.real()? };
    let projection = VolumeProjectionSettings {
        target: reader.real()?, tolerance: reader.real()?, max_shift: reader.real()?, max_evaluations,
    };
    let controls = ProjectedSettings {
        max_candidates, poll_iters, contraction: reader.real()?, min_relative_improvement: reader.real()?,
    };
    let limit = SampledStressLimit::new(reader.real()?, reader.real()?)?;
    let ell = reader.real()?;
    let expected = SampledStressEvaluation {
        compliance: reader.real()?, volume: reader.real()?, sampled_max_von_mises: reader.real()?,
        max_location: [reader.real()?, reader.real()?],
        sample_count: reader.count(100_000_000)?, snapshot: reader.word()?,
    };
    if reader.0.len() != fixed_count * 16 + node_count * 8 {
        return Err("checkpoint has missing or trailing field data".into());
    }
    // Resource admission above precedes lattice allocation and all PDE work.
    let mut fixed = Vec::with_capacity(fixed_count);
    for j in 0..=n {
        for side in [0, n] {
            let index = reader.count(node_count - 1)?;
            if index != j * (n + 1) + side {
                return Err("checkpoint fixed traces are missing, duplicated or reordered".into());
            }
            fixed.push((index, reader.real()?));
        }
    }
    let mut geometry = GridSdf::from_fn(n, &|_, _| 0.0);
    for value in geometry.nodes_mut() { *value = reader.real()?; }
    for &(index, value) in &fixed {
        if value.to_bits() != geometry.nodes()[index].to_bits() {
            return Err("checkpoint geometry violates its retained fixed boundary".into());
        }
    }
    let checkpoint = OptimizeCheckpoint::restore(geometry, fixture, settings, ordinal, ell)?;
    Ok(Saved { checkpoint, policy: Policy { fixed, projection, controls }, limit, expected })
}

fn same_state(a: &SampledStressEvaluation, b: &SampledStressEvaluation) -> bool {
    a.snapshot == b.snapshot && a.sample_count == b.sample_count
        && [a.compliance, a.volume, a.sampled_max_von_mises, a.max_location[0], a.max_location[1]]
            .iter().zip([b.compliance, b.volume, b.sampled_max_von_mises, b.max_location[0], b.max_location[1]])
            .all(|(left, right)| left.to_bits() == right.to_bits())
}

pub(super) fn load_controlled<B>(
    path: &Path, executable: &str,
    control: impl FnMut(ProjectedStressSetupStage) -> ControlFlow<B>,
) -> Result<ControlFlow<B, (ProjectedStressOptimizer, Policy)>, Box<dyn Error>> {
    let mut bytes = Vec::new();
    File::open(path)?.take(MAX_BYTES + 1).read_to_end(&mut bytes)?;
    if bytes.len() as u64 > MAX_BYTES || bytes.len() < 65 || bytes[64] != b'\n' {
        return Err("oversized or truncated stress checkpoint envelope".into());
    }
    if fs_ledger::hash_bytes(&bytes[65..]).to_string().as_bytes() != &bytes[..64] {
        return Err("stress checkpoint content hash mismatch".into());
    }
    let saved = decode(&bytes[65..], executable)?;
    let policy = saved.policy;
    // Neither serialized metrics nor the content hash confer feasibility. Two
    // canonical solves reconstruct mechanics and stress on the retained field.
    let optimizer = match ProjectedStressOptimizer::from_checkpoint_controlled(&saved.checkpoint,
        policy.fixed.clone(), policy.projection, policy.controls, saved.limit, control,
    )? {
        ControlFlow::Continue(optimizer) => optimizer,
        ControlFlow::Break(reason) => return Ok(ControlFlow::Break(reason)),
    };
    if !same_state(optimizer.current(), &saved.expected) {
        return Err("restored mechanics/stress differ from the retained accepted state".into());
    }
    Ok(ControlFlow::Continue((optimizer, policy)))
}

pub(super) fn resume(args: &[String], started: Instant) -> Result<u8, Box<dyn Error>> {
    if args.len() < 4 || args.len() > 6 {
        return Err("usage: --projected --resume CHECKPOINT NEW_OUTPUT_DIR --wall-seconds N [--pause-after N]".into());
    }
    let mut wall = None;
    let mut pause = None;
    let mut i = 2;
    while i < args.len() {
        let value = args.get(i + 1).ok_or("resume option requires a value")?.parse::<u64>()?;
        match args[i].as_str() {
            "--wall-seconds" if wall.is_none() && (1..=3600).contains(&value) => wall = Some(value),
            "--pause-after" if pause.is_none() && value <= 200 => pause = Some(value as usize),
            _ => return Err("unknown, repeated or invalid resume option; study constraints are immutable".into()),
        }
        i += 2;
    }
    let wall = wall.ok_or("resume requires an explicit --wall-seconds budget")?;
    let output = Path::new(&args[1]);
    if output.try_exists()? { return Err("output directory already exists; refusing to overwrite it".into()); }
    if poll_wall(started, wall).is_break() {
        return Ok(stopped_before_study("checkpoint recovery"));
    }
    let executable = executable()?;
    let (optimizer, policy) = match load_controlled(Path::new(&args[0]), &executable,
        |_| poll_wall(started, wall),
    )? {
        ControlFlow::Continue(state) => state,
        ControlFlow::Break(()) => return Ok(stopped_before_study("checkpoint recovery")),
    };
    run_optimizer(output, optimizer, policy, Options { enabled: true, pause_after: pause },
        Some(executable), started, wall, true)
}
