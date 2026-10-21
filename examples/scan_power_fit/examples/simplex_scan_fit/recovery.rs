//! File recovery for this explicitly fixed synthetic physical model.
//! Files are create-only and synced. A crash during a write can leave the newest
//! file incomplete; its checksum then refuses, and earlier snapshots are intact.
//! This is not filesystem transactionality or an accounting ledger for work
//! completed after the last saved search state (including final validation).

use super::{Options, Result, DifferentialEvolutionState};
use fs_dfo::differential_evolution::snapshot::DifferentialEvolutionRestoreLimits;
use fs_exec::Cx;
use std::fs::{File, OpenOptions};
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

const MAX_BYTES: usize = 1 << 20;
const USAGE: &str = "population recovery: --checkpoint-dir DIR [--checkpoint-every N] [--pause-after N] [--resume FILE]";
static FILE_SEQUENCE: AtomicU64 = AtomicU64::new(0);

#[derive(Default)]
pub struct Recovery {
    directory: Option<PathBuf>,
    resume: Option<PathBuf>,
    interval: Option<usize>,
    pub pause_after: Option<usize>,
    identity: [u8; 32],
    pub last_saved: Option<PathBuf>,
}
impl Recovery {
    /// Remove only recovery flags; the existing physical command parser still
    /// validates every other argument. Repeated or incomplete options refuse.
    pub fn parse(args: &[String]) -> Result<(Vec<String>, Self)> {
        let mut out = Self::default();
        let mut remaining = Vec::new();
        let mut i = 0;
        while i < args.len() {
            match args[i].as_str() {
                "--checkpoint-dir" | "--resume" | "--checkpoint-every" | "--pause-after" => {
                    let flag = args[i].as_str();
                    i += 1;
                    let value = args.get(i).filter(|s| !s.is_empty() && !s.starts_with("--")).ok_or(USAGE)?;
                    match flag {
                        "--checkpoint-dir" if out.directory.is_none() => out.directory = Some(PathBuf::from(value.as_str())),
                        "--resume" if out.resume.is_none() => out.resume = Some(PathBuf::from(value.as_str())),
                        "--checkpoint-every" if out.interval.is_none() => out.interval = Some(value.parse()?),
                        "--pause-after" if out.pause_after.is_none() => out.pause_after = Some(value.parse()?),
                        _ => return Err(USAGE.into()),
                    }
                }
                _ => remaining.push(args[i].clone()),
            }
            i += 1;
        }
        if out.interval.into_iter().chain(out.pause_after).any(|n| n == 0 || n > 10_000)
            || ((out.pause_after.is_some() || out.interval.is_some()) && out.directory.is_none())
        { return Err(USAGE.into()); }
        Ok((remaining, out))
    }

    pub fn active(&self) -> bool { self.directory.is_some() || self.resume.is_some() }
    pub fn quantum(&self) -> usize { if self.active() { self.interval.unwrap_or(128) } else { usize::MAX } }

    /// Bind the exact executable plus physical inputs. Worker count, pauses and
    /// budget extension are intentionally not model changes. A rebuilt binary
    /// refuses old files rather than guessing numerical compatibility.
    pub fn bind(&mut self, cx: &Cx<'_>, input: Options) -> Result<()> {
        if !self.active() { return Ok(()); }
        if !input.population { return Err("checkpoint/resume requires --population".into()); }
        let binary = executable_identity(cx)?;
        let mut data = binary.to_vec();
        for value in input.targets.into_iter().chain([input.initial_h, input.maximum_w, super::TARGET]) {
            data.extend_from_slice(&value.to_bits().to_le_bytes());
        }
        data.extend_from_slice(&input.seed.to_le_bytes());
        self.identity = *fs_blake3::hash_domain("simplex-scan-fit/recovery-model-v1", &data).as_bytes();
        Ok(())
    }

    pub fn restore(&self, cx: &Cx<'_>, input: Options) -> Result<Option<DifferentialEvolutionState>> {
        let Some(path) = &self.resume else { return Ok(None); };
        let bytes = read_bounded(path, MAX_BYTES)?;
        let mut state = DifferentialEvolutionState::from_checkpoint_bytes(cx, &bytes, self.identity,
            DifferentialEvolutionRestoreLimits { max_bytes: MAX_BYTES,
                max_dimension: 2, max_population_entries: 48 })?;
        let ceiling = input.evaluations - 1;
        if ceiling < state.evaluation_limit() {
            return Err("resume cannot reduce the saved lifetime budget; use --pause-after for a shorter invocation".into());
        }
        if ceiling > state.evaluation_limit() { state.extend_evaluation_budget(ceiling)?; }
        Ok(Some(state))
    }

    pub fn save(&mut self, cx: &Cx<'_>, state: &DifferentialEvolutionState) -> Result<()> {
        let Some(directory) = &self.directory else { return Ok(()); };
        let bytes = state.checkpoint_bytes(cx, self.identity, MAX_BYTES)?;
        std::fs::create_dir_all(directory)?;
        // Never truncate or overwrite a previous checkpoint, including a file
        // left incomplete by an interrupted write. Each retry gets a new name.
        for _ in 0..32 {
            let serial = FILE_SEQUENCE.fetch_add(1, Ordering::Relaxed);
            let path = directory.join(format!("de-{:08}-{}-{serial:08}.fsde",
                state.evaluations(), std::process::id()));
            let mut file = match OpenOptions::new().write(true).create_new(true).open(&path) {
                Ok(file) => file,
                Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => continue,
                Err(e) => return Err(e.into()),
            };
            file.write_all(&bytes)?;
            file.sync_all()?;
            #[cfg(unix)]
            File::open(directory)?.sync_all()?;
            self.last_saved = Some(path.clone());
            eprintln!("# checkpoint={} search_attempts={}", path.display(), state.evaluations());
            return Ok(());
        }
        Err("checkpoint filenames collided repeatedly; choose another output directory".into())
    }

    pub fn save_after_failure(&mut self, state: &DifferentialEvolutionState) {
        if self.directory.is_none() { return; }
        // Saving a drained search is not additional physical work. It needs a
        // fresh context when the simulation's own context has been cancelled.
        if let Err(error) = super::with_cx(|_, cx| self.save(cx, state)) {
            eprintln!("# failed to save retained search after error: {error}");
        }
    }
}

fn read_bounded(path: &Path, maximum: usize) -> Result<Vec<u8>> {
    let file = File::open(path)?;
    if file.metadata()?.len() > maximum as u64 { return Err("checkpoint file is too large".into()); }
    let mut bytes = Vec::new();
    file.take(maximum as u64 + 1).read_to_end(&mut bytes)?;
    if bytes.len() > maximum { return Err("checkpoint file grew beyond its allowance".into()); }
    Ok(bytes)
}

fn executable_identity(cx: &Cx<'_>) -> Result<[u8; 32]> {
    let mut file = File::open(std::env::current_exe()?)?;
    let mut block = [0_u8; 64 * 1024];
    let mut digest = [0_u8; 32];
    let mut total = 0_u64;
    loop {
        cx.checkpoint().map_err(|_| "checkpoint identity cancelled")?;
        let mut count = 0;
        while count < block.len() {
            let n = match file.read(&mut block[count..]) {
                Ok(n) => n,
                Err(e) if e.kind() == std::io::ErrorKind::Interrupted => continue,
                Err(e) => return Err(e.into()),
            };
            if n == 0 { break; }
            count += n;
        }
        if count == 0 { break; }
        total += count as u64;
        if total > 1 << 30 { return Err("executable exceeds the 1 GiB identity-read allowance".into()); }
        let mut data = Vec::with_capacity(40 + count);
        data.extend_from_slice(&digest);
        data.extend_from_slice(&(count as u64).to_le_bytes());
        data.extend_from_slice(&block[..count]);
        digest = *fs_blake3::hash_domain("simplex-scan-fit/executable-block-v1", &data).as_bytes();
    }
    if total == 0 { return Err("cannot identify an empty executable".into()); }
    cx.checkpoint().map_err(|_| "checkpoint identity cancelled")?;
    Ok(digest)
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::{FitResult, fit, fit_with_recovery, options, with_cx};

    fn directory() -> PathBuf {
        std::env::temp_dir().join(format!("frankensim-de-recovery-test-{}-{}",
            std::process::id(), FILE_SEQUENCE.fetch_add(1, Ordering::Relaxed)))
    }

    #[test]
    fn g5_file_pause_and_new_context_resume_match_uninterrupted_physics() {
        let dir = directory();
        let args: Vec<String> = ["--population", "--workers", "1", "--initial-enthalpy", "80",
            "0.1", "0.3", "50", "65"].map(str::to_string).into();
        let input = options(&args).unwrap();
        let mut persistence = args.clone();
        persistence.extend(["--checkpoint-dir".into(), dir.to_string_lossy().into_owned(),
            "--checkpoint-every".into(), "5".into(), "--pause-after".into(), "13".into()]);
        let (_, mut first) = Recovery::parse(&persistence).unwrap();
        let saved = with_cx(|_, cx| {
            match fit_with_recovery(cx, input, false, &mut first)? {
                FitResult::Paused { search, checkpoint } => {
                    assert_eq!(search.evaluations(), 13);
                    Ok(checkpoint)
                }
                FitResult::Finished(_) => Err("expected deliberate checkpoint pause".into()),
            }
        }).unwrap();
        let original = read_bounded(&saved, MAX_BYTES).unwrap();
        drop(first);
        let (_, mut second) = Recovery::parse(&["--resume".into(), saved.to_string_lossy().into_owned()]).unwrap();
        let many = Options { workers: Some(4), ..input };
        let resumed = with_cx(|_, cx| match fit_with_recovery(cx, many, false, &mut second)? {
            FitResult::Finished(result) => Ok(result),
            FitResult::Paused { .. } => Err("unrequested second pause".into()),
        }).unwrap();
        let uninterrupted = with_cx(|_, cx| fit(cx, many, false)).unwrap();
        assert_eq!(resumed.search, uninterrupted.search);
        assert_eq!(resumed.sources, uninterrupted.sources);
        assert_eq!(resumed.replay.physical.predictions, uninterrupted.replay.physical.predictions);
        assert_eq!(resumed.spent, uninterrupted.spent);
        assert_eq!(read_bounded(&saved, MAX_BYTES).unwrap(), original);
    }

    #[test]
    fn g0_resume_refuses_changed_physics_corruption_and_budget_reduction() {
        let dir = directory();
        let input = options(&["--population", "0.1", "0.3", "50", "65"].map(str::to_string)).unwrap();
        let (_, mut saving) = Recovery::parse(&["--checkpoint-dir".into(), dir.to_string_lossy().into_owned(),
            "--pause-after".into(), "3".into()]).unwrap();
        let path = with_cx(|_, cx| match fit_with_recovery(cx, input, false, &mut saving)? {
            FitResult::Paused { checkpoint, .. } => Ok(checkpoint),
            _ => Err("expected pause".into()),
        }).unwrap();
        let (_, mut loading) = Recovery::parse(&["--resume".into(), path.to_string_lossy().into_owned()]).unwrap();
        with_cx(|_, cx| {
            for changed in [Options { targets: [0.2, 0.3], ..input },
                Options { initial_h: 80.0, ..input }, Options { maximum_w: 40.0, ..input },
                Options { seed: 18, ..input }, Options { evaluations: 64, ..input }]
            {
                loading.bind(cx, changed)?;
                assert!(loading.restore(cx, changed).is_err());
            }
            loading.bind(cx, input)?;
            assert_eq!(loading.restore(cx, Options { evaluations: 100, ..input })?
                .unwrap().evaluation_limit(), 99);
            let mut bad = read_bounded(&path, MAX_BYTES)?;
            let last = bad.len()-1;
            bad[last] ^= 1;
            let corrupt = dir.join("corrupt.fsde");
            OpenOptions::new().write(true).create_new(true).open(&corrupt)?.write_all(&bad)?;
            loading.resume = Some(corrupt);
            assert!(loading.restore(cx, input).is_err());
            Ok(())
        }).unwrap();
    }

    #[test]
    fn g0_recovery_options_are_explicit_and_cancelled_identity_does_no_save() {
        for words in [vec!["--resume"], vec!["--checkpoint-every", "5"],
            vec!["--pause-after", "1"], vec!["--checkpoint-dir", "x", "--pause-after", "0"],
            vec!["--resume", "x", "--resume", "y"]]
        { assert!(Recovery::parse(&words.into_iter().map(str::to_string).collect::<Vec<_>>()).is_err()); }
        let (_, mut recovery) = Recovery::parse(&["--checkpoint-dir".into(), directory().to_string_lossy().into_owned()]).unwrap();
        with_cx(|gate, cx| {
            assert!(recovery.bind(cx, options(&[])?).is_err());
            let input = options(&["--population".into()])?;
            gate.request();
            assert!(recovery.bind(cx, input).is_err());
            assert!(recovery.last_saved.is_none());
            Ok(())
        }).unwrap();
    }
}
