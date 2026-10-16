//! Explicit fixed-grid restart. Files are immutable, bounded and checked before
//! physics. The library owns optimizer restoration; this module only binds its
//! bytes to this executable, problem and original work allowances.
use super::*;
use std::{fs::{File, OpenOptions}, io::{Read, Write}, path::PathBuf};
use fs_blake3::{ContentHash, DomainHasher, hash_domain};
use fs_topopt::SolveWork;
use fs_cutfem::quad3::QuadratureWork3;
use fs_ascent::projected_al::ProjectedAlStop;

type Error = Box<dyn std::error::Error>;
const MAGIC: &[u8; 8] = b"FSMMRUN1";
const MAX_BYTES: usize = 64 * 1024 * 1024;
const FILE_DOMAIN: &str = "org.frankensim.example.minimum-material.restart.v1";

pub(super) struct Config { directory: PathBuf, source: Option<PathBuf> }
pub(super) fn parse(args: &mut Vec<String>) -> Result<Option<Config>, Error> {
    if !args.iter().any(|a| matches!(a.as_str(), "--checkpoint-dir" | "--resume")) { return Ok(None); }
    let mut directory = None; let mut source = None; let mut positional = Vec::new(); let mut i = 0;
    while i < args.len() {
        let arg = &args[i];
        let dest = match arg.as_str() {
            "--checkpoint-dir" => Some(&mut directory), "--resume" => Some(&mut source),
            _ => None,
        };
        if let Some(dest) = dest {
            if dest.is_some() { return Err(format!("duplicate {arg}").into()); }
            let path = args.get(i+1).filter(|s| !s.is_empty() && !s.starts_with("--"))
                .ok_or_else(|| format!("missing path after {arg}"))?;
            *dest = Some(PathBuf::from(path)); i += 2;
        } else if arg.starts_with("--") {
            return Err("restart is fixed-grid: do not combine it with binary/adaptive/estimate modes".into());
        } else { positional.push(arg.clone()); i += 1; }
    }
    let directory = directory.ok_or("--resume requires a NEW --checkpoint-dir to retain further spent work")?;
    *args = positional; Ok(Some(Config { directory, source }))
}
fn add(a: usize, b: usize) -> Result<usize, Error> { a.checked_add(b).ok_or_else(|| "restart work counter overflow".into()) }
fn total(a: SolveWork, b: SolveWork) -> Result<SolveWork, Error> {
    Ok(SolveWork { linear_solves: add(a.linear_solves, b.linear_solves)?,
        linear_iterations: add(a.linear_iterations, b.linear_iterations)?,
        preconditioner_operator_applications: add(a.preconditioner_operator_applications, b.preconditioner_operator_applications)?,
        preconditioner_galerkin_products: add(a.preconditioner_galerkin_products, b.preconditioner_galerkin_products)? })
}
fn executable() -> Result<ContentHash, Error> {
    let mut file = File::open(std::env::current_exe()?)?;
    let mut h = DomainHasher::new("org.frankensim.example.minimum-material.executable.v1");
    let mut block = [0u8; 65536];
    loop { let n = file.read(&mut block)?; if n == 0 { break; } h.update(&block[..n]); }
    Ok(h.finalize())
}
struct Saved { checkpoint: MinimumMaterialCheckpoint3, work: SolveWork, geometry: QuadratureWork3 }
struct Binding { executable: ContentHash, budget: usize, lo: f64, hi: f64 }
fn number(bytes: &mut &[u8]) -> Result<u64, Error> {
    let head = bytes.get(..8).ok_or("truncated restart header")?;
    let value = u64::from_le_bytes(head.try_into()?); *bytes = &bytes[8..]; Ok(value)
}
fn count(bytes: &mut &[u8]) -> Result<usize, Error> { Ok(usize::try_from(number(bytes)?)?) }
fn decode(bytes: &[u8], binding: &Binding) -> Result<Saved, Error> {
    if bytes.len() > MAX_BYTES || bytes.len() < 8+32+10*8+32 { return Err("restart file size refused".into()); }
    let (payload, digest) = bytes.split_at(bytes.len()-32);
    if hash_domain(FILE_DOMAIN, payload).as_bytes().as_slice() != digest { return Err("restart content checksum mismatch".into()); }
    if &payload[..8] != MAGIC || &payload[8..40] != binding.executable.as_bytes() {
        return Err("restart version or executable changed".into());
    }
    let mut r = &payload[40..];
    if count(&mut r)? != binding.budget || number(&mut r)? != binding.lo.to_bits() || number(&mut r)? != binding.hi.to_bits() {
        return Err("restart changed its original Krylov budget or physical bounds".into());
    }
    let work = SolveWork { linear_solves: count(&mut r)?, linear_iterations: count(&mut r)?,
        preconditioner_operator_applications: count(&mut r)?, preconditioner_galerkin_products: count(&mut r)? };
    let geometry = QuadratureWork3 { boxes: count(&mut r)?, field_evaluations: count(&mut r)?, points: count(&mut r)? };
    if work.linear_iterations > binding.budget { return Err("restart work exceeds original budget".into()); }
    let checkpoint = MinimumMaterialCheckpoint3::from_bytes(r, Default::default())?;
    Ok(Saved { checkpoint, work, geometry })
}
fn encode(saved: &Saved, binding: &Binding) -> Result<Vec<u8>, Error> {
    let body = saved.checkpoint.to_bytes(MinimumMaterialCodecBudget3 { max_bytes: MAX_BYTES-152, ..Default::default() })?;
    let mut bytes = Vec::new(); bytes.extend(MAGIC); bytes.extend(binding.executable.as_bytes());
    bytes.extend(u64::try_from(binding.budget)?.to_le_bytes()); bytes.extend(binding.lo.to_bits().to_le_bytes()); bytes.extend(binding.hi.to_bits().to_le_bytes());
    for n in [saved.work.linear_solves, saved.work.linear_iterations, saved.work.preconditioner_operator_applications,
        saved.work.preconditioner_galerkin_products, saved.geometry.boxes, saved.geometry.field_evaluations, saved.geometry.points] {
        bytes.extend(u64::try_from(n)?.to_le_bytes());
    }
    bytes.extend(body); let digest = hash_domain(FILE_DOMAIN, &bytes); bytes.extend(digest.as_bytes()); Ok(bytes)
}
struct Output { directory: PathBuf, next: usize, binding: Binding, prior: SolveWork, geometry: QuadratureWork3 }
impl Output {
    fn save(&mut self, checkpoint: &MinimumMaterialCheckpoint3, work: SolveWork) -> Result<(), Error> {
        let saved = Saved { checkpoint: checkpoint.clone(), work: total(self.prior, work)?, geometry: self.geometry };
        let bytes = encode(&saved, &self.binding)?;
        let path = self.directory.join(format!("{:08}.fsmcp", self.next));
        let mut file = OpenOptions::new().write(true).create_new(true).open(&path)?;
        file.write_all(&bytes)?; file.sync_all()?;
        self.next = self.next.checked_add(1).ok_or("checkpoint file counter overflow")?;
        eprintln!("checkpoint={}; cumulative_work={:?}", path.display(), saved.work); Ok(())
    }
}

// Returning from this scope releases the solver/control borrows even when
// restoration fails. The outer scope persists the UPDATED mutable reservation
// and exact consumed solver work rather than overwriting them with a new budget.
fn fit(study: &mut CutDensityStudy3<AdaptiveSolveSpace3>, cases: &[ResponseLimitCase3<'_>], steps: usize,
    checkpoint: &mut Option<MinimumMaterialCheckpoint3>, output: &mut Output, control: &mut SolveControl<'_>) -> Result<(), Error> {
    let options = MinimumMaterialOptions3::default();
    let initial = vec![0.5; study.cells()];
    let mut session = match checkpoint {
        Some(saved) => MinimumMaterialStudy3::restore(study, cases, saved, options, control)?,
        None => MinimumMaterialStudy3::new(study, cases, &initial, options, control)?,
    };
    *checkpoint = Some(session.checkpoint()); output.save(checkpoint.as_ref().unwrap(), session.work())?;
    let mut outcome = session.run(0);
    for _ in 0..steps {
        if !outcome.as_ref().is_ok_and(|r| r.stop == ProjectedAlStop::IterationLimit) { break; }
        outcome = session.run(1);
        *checkpoint = Some(session.checkpoint()); output.save(checkpoint.as_ref().unwrap(), session.work())?;
    }
    *checkpoint = Some(session.checkpoint());
    println!("iteration,volume_fraction,smooth_constraint,worst_original_residual");
    for row in session.history() { println!("{},{:.17e},{:.17e},{:.17e}", row.iteration, row.volume_fraction, row.constraint, row.worst_residual); }
    let best = session.best_feasible().unwrap_or(session.accepted());
    println!("case,limit,value,bound,original_residual");
    for (i, case) in cases.iter().enumerate() { for (j, limit) in case.limits.iter().enumerate() {
        println!("{i},{j},{:.17e},{:.17e},{:.17e}", best.measurements[i][j], limit.bound, best.residuals[i][j]);
    } }
    println!("cell,raw_density,projected_density");
    for (i, (r, p)) in best.rho.iter().zip(&best.projected_rho).enumerate() { println!("{i},{r:.17e},{p:.17e}"); }
    eprintln!("selected_latest={}; feasible={}; optimizer_evaluations={}; restoration_reservations={}; stop={:?}; continuum_certified=false",
        best.rho == session.accepted().rho, best.feasible(options.optimizer.tolerance),
        session.optimizer_work().evaluations, checkpoint.as_ref().unwrap().restoration_evaluations, outcome);
    outcome?;
    if !best.feasible(options.optimizer.tolerance) { return Err("no accepted response-limit-feasible endpoint; checkpoint remains resumable".into()); }
    Ok(())
}

pub(super) fn run(steps: usize, budget: usize, lo: f64, hi: f64, config: Config) -> Result<(), Error> {
    let binding = Binding { executable: executable()?, budget, lo, hi };
    let prior = if let Some(path) = &config.source {
        let mut bytes = Vec::new(); File::open(path)?.take((MAX_BYTES+1) as u64).read_to_end(&mut bytes)?;
        Some(decode(&bytes, &binding)?)
    } else { None };
    let old_geometry = prior.as_ref().map_or(QuadratureWork3::default(), |s| s.geometry);
    let old_work = prior.as_ref().map_or(SolveWork::default(), |s| s.work);
    let mut checkpoint = prior.map(|s| s.checkpoint);
    let defaults = QuadratureOptions3::default();
    let geometry_options = QuadratureOptions3 {
        max_boxes: defaults.max_boxes.checked_sub(old_geometry.boxes).ok_or("original geometry box budget exhausted")?,
        max_points: defaults.max_points.checked_sub(old_geometry.points).ok_or("original geometry point budget exhausted")?, ..defaults };
    let remaining = budget.checked_sub(old_work.linear_iterations).ok_or("original Krylov budget exhausted")?;
    // Never overwrite an earlier run, including the source checkpoint directory.
    std::fs::create_dir(&config.directory)?;
    let mut output = Output { directory: config.directory, next: 0, binding, prior: old_work, geometry: old_geometry };
    let mut gp = |_| ControlFlow::Continue(());
    let mut geometry = QuadratureControl3::new(geometry_options, &mut gp)?;
    let tree = Octree3::uniform(1, 4, 4096)?;
    let built = AdaptiveElasticity3::build_with_embedded_dirichlet(HexCell::try_new([0.0; 3], [1.0; 3])?,
        &tree, &Slab, &IsotropicElastic::new(1.0, 0.0, 1.0)?, &|_| false, &|_, _| true,
        ElasticityOptions3::default(), Default::default(), Default::default(), &mut geometry);
    let spent = geometry.work();
    output.geometry = QuadratureWork3 { boxes: add(old_geometry.boxes, spent.boxes)?,
        field_evaluations: add(old_geometry.field_evaluations, spent.field_evaluations)?, points: add(old_geometry.points, spent.points)? };
    let op = match built { Ok(op) => op, Err(error) => {
        if let Some(saved) = &checkpoint { output.save(saved, SolveWork::default())?; }
        return Err(error.into());
    } };
    let f = op.body_load(&|_| [0.0; 3], || ControlFlow::Continue(()))?;
    let q = op.body_load(&|_| [1.0, 0.0, 0.0], || ControlFlow::Continue(()))?;
    let mut study = CutDensityStudy3::new(AdaptiveSolveSpace3::jacobi(op, 100_000_000), 0.15, SimpParams::default());
    let tension = [
        ResponseLimit3 { response: BoundedResponse3::Reaction(&right), sense: ResponseLimitSense3::AtLeast, bound: lo, scale: 0.002 },
        ResponseLimit3 { response: BoundedResponse3::Reaction(&right), sense: ResponseLimitSense3::AtMost, bound: hi, scale: 0.002 },
        ResponseLimit3 { response: BoundedResponse3::Displacement(&q), sense: ResponseLimitSense3::AtMost, bound: 0.006, scale: 0.005 } ];
    let negative = [ResponseLimit3 { sense: ResponseLimitSense3::AtMost, bound: -lo, ..tension[0] },
        ResponseLimit3 { sense: ResponseLimitSense3::AtLeast, bound: -hi, ..tension[1] },
        ResponseLimit3 { sense: ResponseLimitSense3::AtLeast, bound: -0.006, ..tension[2] }];
    let cases = [ResponseLimitCase3 { force: &f, prescribed: Some(&motion), limits: &tension },
        ResponseLimitCase3 { force: &f, prescribed: Some(&compression), limits: &negative }];
    let mut cp = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget { total_iterations: remaining, ..Default::default() }, &mut cp);
    let result = fit(&mut study, &cases, steps, &mut checkpoint, &mut output, &mut control);
    if let Some(saved) = &checkpoint { output.save(saved, control.work())?; }
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::{AtomicUsize, Ordering};
    static NEXT: AtomicUsize = AtomicUsize::new(0);
    fn path(label: &str) -> PathBuf {
        std::env::temp_dir().join(format!("fs-minimum-restart-{}-{}-{label}", std::process::id(), NEXT.fetch_add(1, Ordering::Relaxed)))
    }
    fn last(directory: &std::path::Path) -> PathBuf {
        let mut entries: Vec<_> = std::fs::read_dir(directory).unwrap().map(|e| e.unwrap().path()).collect();
        entries.sort(); entries.pop().unwrap()
    }
    fn args(s: &str) -> Vec<String> { s.split_whitespace().map(str::to_owned).collect() }

    #[test]
    fn restart_arguments_are_explicit_and_legacy_arguments_are_unchanged() {
        let mut a = args("2 250000 --resume prior.fsmcp --checkpoint-dir next-run");
        let cfg = parse(&mut a).unwrap().unwrap(); assert_eq!(a, args("2 250000"));
        assert_eq!(cfg.source.unwrap(), PathBuf::from("prior.fsmcp"));
        for text in ["--resume prior.fsmcp", "--checkpoint-dir", "--checkpoint-dir next --refit",
            "--checkpoint-dir next --binary", "--checkpoint-dir next --checkpoint-dir other"] {
            let mut a = args(text); let before = a.clone(); assert!(parse(&mut a).is_err()); assert_eq!(a, before);
        }
        for text in ["80 250000", "80 --binary", "--adapt-rounds 2"] {
            let mut a = args(text); let before = a.clone(); assert!(parse(&mut a).unwrap().is_none()); assert_eq!(a, before);
        }
    }

    #[test]
    fn actual_disk_restart_matches_uninterrupted_search_and_retains_original_work() {
        let one = path("first"); let two = path("second"); let all = path("all");
        run(1, 250000, 0.0012, 0.003, Config { directory: one.clone(), source: None }).unwrap();
        let first_path = last(&one); let bytes_before = std::fs::read(&first_path).unwrap();
        run(1, 250000, 0.0012, 0.003, Config { directory: two.clone(), source: Some(first_path.clone()) }).unwrap();
        run(2, 250000, 0.0012, 0.003, Config { directory: all.clone(), source: None }).unwrap();
        let binding = Binding { executable: executable().unwrap(), budget: 250000, lo: 0.0012, hi: 0.003 };
        let first = decode(&bytes_before, &binding).unwrap();
        let mut resumed = decode(&std::fs::read(last(&two)).unwrap(), &binding).unwrap();
        let full = decode(&std::fs::read(last(&all)).unwrap(), &binding).unwrap();
        assert_eq!(first.checkpoint.optimizer.work.iterations, 1);
        assert_eq!(full.checkpoint.optimizer.work.iterations, 2);
        assert_eq!(resumed.checkpoint.accepted, full.checkpoint.accepted);
        assert_eq!(resumed.checkpoint.best_feasible, full.checkpoint.best_feasible);
        assert_eq!(resumed.checkpoint.history, full.checkpoint.history);
        let cost = resumed.checkpoint.restoration_evaluations;
        assert!(cost > 0); resumed.checkpoint.optimizer.work.evaluations -= cost;
        assert_eq!(resumed.checkpoint.optimizer, full.checkpoint.optimizer);
        assert!(resumed.work.linear_iterations > full.work.linear_iterations);
        assert_eq!(resumed.geometry.boxes, 2*full.geometry.boxes);
        assert_eq!(std::fs::read(first_path).unwrap(), bytes_before);
        assert!(run(0, 250000, 0.0012, 0.003, Config { directory: one, source: None }).is_err());
    }

    #[test]
    fn changed_problem_build_checksum_or_truncation_refuses_before_new_output() {
        let dir = path("source");
        run(0, 250000, 0.0012, 0.003, Config { directory: dir.clone(), source: None }).unwrap();
        let source = last(&dir); let bytes = std::fs::read(&source).unwrap();
        let mut binding = Binding { executable: executable().unwrap(), budget: 250000, lo: 0.0012, hi: 0.003 };
        assert!(decode(&bytes[..bytes.len()-1], &binding).is_err());
        let mut corrupt = bytes.clone(); corrupt[bytes.len()-40] ^= 1; assert!(decode(&corrupt, &binding).is_err());
        binding.executable.0[0] ^= 1; assert!(decode(&bytes, &binding).is_err());
        let rejected = path("rejected");
        assert!(run(0, 250001, 0.0012, 0.003, Config { directory: rejected.clone(), source: Some(source) }).is_err());
        assert!(!rejected.exists());
        assert_eq!(std::fs::read(last(&dir)).unwrap(), bytes);
    }
}
