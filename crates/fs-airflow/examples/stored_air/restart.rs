//! Persistent recovery for the actual stored-air FEM example, not a new solver.
//! Checkpoints bind the exact executable, mesh/order, forcing and estimator.
//! Each invocation writes a NEW directory; previous generations are never replaced.
//! A synchronized pending file is renamed only in that writer-exclusive directory.
//! Interrupted pending writes are not selected by directory resume. Unix directory
//! sync is required; power-loss semantics of a filesystem remain platform-specific.
use std::fs::{self, File, OpenOptions};
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use fs_blake3::DomainHasher;
use fs_couple::iqn_ils::driver::march::adaptive::{AdaptiveEvolution, AdaptiveReport};
use super::{Inputs, Model, State, INITIAL_K, RHO_CP, HTC};

const LIMIT: usize = 1024 * 1024;
const EXECUTABLE_LIMIT: u64 = 2 * 1024 * 1024 * 1024;
const MAX_FILES: usize = 20_002;
const MAGIC: &[u8;8] = b"FSSA0001";

/// Cumulative work includes rejected attempts and carries across processes.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub(super) struct Progress { pub attempts: u64, pub evaluations: u64, pub rejected: u64 }
impl Progress {
    pub fn record(&mut self, report: &AdaptiveReport) -> Result<(), String> {
        let attempts = self.attempts.checked_add(report.attempts as u64).ok_or("attempt count overflow")?;
        let evaluations = self.evaluations.checked_add(report.evaluations as u64).ok_or("evaluation count overflow")?;
        let rejected = self.rejected.checked_add(report.rejected as u64).ok_or("rejection count overflow")?;
        *self = Self { attempts, evaluations, rejected }; Ok(())
    }
}

pub(super) fn executable_identity() -> Result<[u8;32], String> {
    let mut file = File::open(std::env::current_exe().map_err(|e| e.to_string())?)
        .map_err(|e| format!("cannot bind the current executable: {e}"))?;
    if file.metadata().map_err(|e| e.to_string())?.len() > EXECUTABLE_LIMIT {
        return Err("executable exceeds the 2 GiB restart-identity read cap".into());
    }
    let mut h = DomainHasher::new("org.frankensim.stored-air.executable.v1");
    let mut block = [0_u8; 65536]; let mut total = 0_u64;
    loop {
        let n = file.read(&mut block).map_err(|e| e.to_string())?;
        if n == 0 { break; }
        total += n as u64;
        if total > EXECUTABLE_LIMIT { return Err("executable changed beyond read cap".into()); }
        h.update(&block[..n]);
    }
    Ok(*h.finalize().as_bytes())
}

pub(super) fn binding(model: &Model<'_>, inputs: Inputs, executable: &[u8;32]) -> [u8;32] {
    let mut h = DomainHasher::new("org.frankensim.stored-air.model-restart.v1");
    h.update(executable);
    // Runtime policy/endpoint identity is separately bound by AdaptiveEvolution.
    // --attempts is an explicit ADDITIONAL invocation budget, not a model change.
    for value in [RHO_CP, INITIAL_K, HTC, 1.5, 2.0e6, model.capacity_j_k,
        model.ventilation_w_k, inputs.tolerance_k] { h.update(&value.to_bits().to_le_bytes()); }
    h.update(&(model.mesh.vertex_count() as u64).to_le_bytes());
    for p in model.mesh.positions() { for value in p { h.update(&value.to_bits().to_le_bytes()); } }
    h.update(&(model.mesh.element_count() as u64).to_le_bytes());
    for tet in &model.mesh.complex().tets { for v in tet { h.update(&v.to_le_bytes()); } }
    h.update(&(model.duty.segments().len() as u64).to_le_bytes());
    for segment in model.duty.segments() {
        for value in [segment.duration_s(), segment.start_scale(), segment.end_scale()] {
            h.update(&value.to_bits().to_le_bytes());
        }
        h.update(&[match segment.interpolation() {
            fs_conduction::duty::SegmentInterpolation::Constant => 0,
            fs_conduction::duty::SegmentInterpolation::Linear => 1,
        }]);
    }
    *h.finalize().as_bytes()
}

fn encode(state: &State, progress: Progress) -> Result<Vec<u8>, String> {
    let length = state.solid_k.len().checked_mul(8).and_then(|n| n.checked_add(64))
        .filter(|&n| n < LIMIT).ok_or("physical checkpoint payload too large")?;
    let mut out = Vec::new(); out.try_reserve_exact(length).map_err(|_| "checkpoint allocation")?;
    out.extend_from_slice(MAGIC);
    for n in [progress.attempts, progress.evaluations, progress.rejected, state.solid_k.len() as u64] {
        out.extend_from_slice(&n.to_le_bytes());
    }
    for value in state.solid_k.iter().copied().chain([state.air_k, state.net_input_j, state.source_input_j]) {
        out.extend_from_slice(&value.to_bits().to_le_bytes());
    }
    Ok(out)
}
fn decode(bytes: &[u8], vertices: usize) -> Result<(State, Progress), String> {
    let expected = vertices.checked_mul(8).and_then(|n| n.checked_add(64)).ok_or("field too large")?;
    if bytes.len() != expected || bytes.len() < 64 || &bytes[..8] != MAGIC {
        return Err("wrong physical checkpoint version, length or field shape".into());
    }
    let word = |at: usize| -> [u8;8] { bytes[at..at+8].try_into().expect("checked fixed payload layout") };
    if u64::from_le_bytes(word(32)) != vertices as u64 { return Err("checkpoint mesh vertex count differs".into()); }
    let progress = Progress { attempts: u64::from_le_bytes(word(8)),
        evaluations: u64::from_le_bytes(word(16)), rejected: u64::from_le_bytes(word(24)) };
    let mut values = Vec::with_capacity(vertices+3);
    for at in (40..bytes.len()).step_by(8) { values.push(f64::from_le_bytes(word(at))); }
    if values.iter().any(|v| !v.is_finite()) || values[..vertices+1].iter().any(|v| *v <= 0.0)
        || values[vertices+2] < 0.0
    { return Err("checkpoint temperatures/energy must be finite and physically admissible".into()); }
    let state = State { solid_k: values[..vertices].to_vec(), air_k: values[vertices],
        net_input_j: values[vertices+1], source_input_j: values[vertices+2] };
    Ok((state, progress))
}

pub(super) fn bytes(run: &AdaptiveEvolution<State>, progress: Progress, identity: &[u8;32])
    -> Result<Vec<u8>, String>
{
    run.checkpoint_bytes(identity, &encode(run.state(), progress)?, LIMIT).map_err(|e| e.to_string())
}

pub(super) fn restore(run: &mut AdaptiveEvolution<State>, bytes: &[u8], identity: &[u8;32])
    -> Result<Progress, String>
{
    let vertices = run.state().solid_k.len();
    let mut candidate = run.clone(); let mut progress = None;
    candidate.restore_checkpoint(bytes, identity, LIMIT, |payload| {
        let (state, work) = decode(payload, vertices)?; progress = Some(work); Ok(state)
    }).map_err(|e| e.to_string())?;
    let work = progress.ok_or("missing checkpoint work counts")?;
    let accepted = candidate.accepted_steps() as u64;
    let maximum = work.attempts.checked_mul(96).ok_or("checkpoint work range overflow")?;
    if accepted > work.attempts || work.rejected > work.attempts - accepted || work.evaluations > maximum {
        return Err("checkpoint work counts contradict its accepted prefix".into());
    }
    *run = candidate; Ok(work)
}

/// A file is explicit; a directory resolves its lexically last completed generation.
pub(super) fn read(path: &Path) -> Result<Vec<u8>, String> {
    let path = if path.is_dir() {
        let mut latest = None; let mut count = 0; let mut scanned = 0;
        for entry in fs::read_dir(path).map_err(|e| e.to_string())? {
            scanned += 1;
            if scanned > 2 * MAX_FILES { return Err("checkpoint directory scan cap exceeded".into()); }
            let entry = entry.map_err(|e| e.to_string())?;
            let name = entry.file_name(); let Some(name) = name.to_str() else { continue; };
            if name.starts_with("checkpoint-") && name.ends_with(".fscp") {
                count += 1; if count > MAX_FILES { return Err("checkpoint directory exceeds generation cap".into()); }
                let number = name.strip_prefix("checkpoint-").unwrap().strip_suffix(".fscp").unwrap();
                if number.len() != 8 || !number.bytes().all(|c| c.is_ascii_digit()) {
                    return Err("malformed checkpoint generation name".into());
                }
                let item = (name.to_owned(), entry.path());
                if latest.as_ref().is_none_or(|old: &(String, PathBuf)| item.0 > old.0) { latest = Some(item); }
            }
        }
        latest.ok_or("directory contains no completed checkpoint")?.1
    } else { path.to_owned() };
    let file = File::open(&path).map_err(|e| format!("cannot read checkpoint {}: {e}", path.display()))?;
    if file.metadata().map_err(|e| e.to_string())?.len() > LIMIT as u64 { return Err("checkpoint exceeds 1 MiB".into()); }
    let mut bytes = Vec::new(); file.take((LIMIT+1) as u64).read_to_end(&mut bytes).map_err(|e| e.to_string())?;
    if bytes.len() > LIMIT { return Err("checkpoint exceeds 1 MiB".into()); } Ok(bytes)
}

pub(super) struct Writer { directory: PathBuf, generation: usize }
impl Writer {
    pub fn new(directory: &Path) -> Result<Self, String> {
        fs::create_dir(directory).map_err(|e|
            format!("checkpoint output must be a NEW directory ({}): {e}", directory.display()))?;
        Ok(Self { directory: directory.to_owned(), generation: 0 })
    }
    pub fn publish(&mut self, run: &AdaptiveEvolution<State>, progress: Progress, identity: &[u8;32])
        -> Result<PathBuf, String>
    {
        if self.generation >= MAX_FILES { return Err("checkpoint generation cap reached".into()); }
        let bytes = bytes(run, progress, identity)?;
        let name = format!("checkpoint-{:08}", self.generation);
        let pending = self.directory.join(format!("{name}.pending"));
        let final_path = self.directory.join(format!("{name}.fscp"));
        if final_path.try_exists().map_err(|e| e.to_string())? { return Err("checkpoint generation already exists".into()); }
        let mut file = OpenOptions::new().create_new(true).write(true).open(&pending).map_err(|e| e.to_string())?;
        file.write_all(&bytes).and_then(|()| file.sync_all()).map_err(|e| e.to_string())?;
        drop(file);
        fs::rename(&pending, &final_path).map_err(|e| e.to_string())?;
        #[cfg(unix)]
        File::open(&self.directory).and_then(|d| d.sync_all()).map_err(|e| e.to_string())?;
        self.generation += 1; Ok(final_path)
    }
}

#[cfg(test)]
#[path = "restart/tests.rs"]
mod tests;
