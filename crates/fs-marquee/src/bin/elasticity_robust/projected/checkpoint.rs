//! Opt-in durable accepted-state checkpoints for the existing projected runner.
use super::*;

#[derive(Default)]
pub(super) struct Options {
    pub enabled: bool,
    pub pause_after: Option<usize>,
    pub recovery_solves: usize,
}

pub(super) fn options(args: &[String]) -> Result<(Vec<String>, Options), Box<dyn Error>> {
    let mut out = Options::default();
    let mut checkpoint_seen = false;
    let mut positional = Vec::new();
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--checkpoint" => {
                if checkpoint_seen { return Err("duplicate --checkpoint".into()); }
                checkpoint_seen = true;
                out.enabled = true;
            }
            "--pause-after" => {
                if out.pause_after.is_some() { return Err("duplicate --pause-after".into()); }
                i += 1;
                let value = args.get(i).ok_or("--pause-after requires an accepted-update count")?.parse::<usize>()?;
                if value > 200 { return Err("--pause-after must lie in 0..=200".into()); }
                out.pause_after = Some(value);
                out.enabled = true;
            }
            _ => positional.push(args[i].clone()),
        }
        i += 1;
    }
    Ok((positional, out))
}

pub(super) fn save(path: &Path, optimizer: &MultiLoadProjectedOptimizer) -> Result<(), Box<dyn Error>> {
    let bytes = optimizer.checkpoint_bytes();
    let digest = fs_ledger::hash_bytes(&bytes).to_string();
    if digest.len() != 64 { return Err("unexpected checkpoint digest format".into()); }
    let mut file = writer(path)?;
    writeln!(file, "{digest}")?;
    file.write_all(&bytes)?;
    file.flush()?;
    file.get_ref().sync_all()?;
    Ok(())
}

pub(super) fn load(path: &Path) -> Result<Vec<u8>, Box<dyn Error>> {
    const MAX_BYTES: u64 = 4 * 1024 * 1024 + 65;
    let mut bytes = Vec::new();
    File::open(path)?.take(MAX_BYTES + 1).read_to_end(&mut bytes)?;
    if bytes.len() as u64 > MAX_BYTES || bytes.len() < 65 || bytes[64] != b'\n' {
        return Err("oversized or malformed projected checkpoint envelope".into());
    }
    if fs_ledger::hash_bytes(&bytes[65..]).to_string().as_bytes() != &bytes[..64] {
        return Err("projected checkpoint content hash mismatch".into());
    }
    Ok(bytes.split_off(65))
}

pub(super) fn resume(args: &[String]) -> Result<u8, Box<dyn Error>> {
    if args.len() < 4 || args.len() > 6 {
        return Err("usage: --projected --resume CHECKPOINT.fscp NEW_OUTPUT_DIR --recovery-solves N [--pause-after N]".into());
    }
    let mut recovery = None;
    let mut pause = None;
    let mut i = 2;
    while i < args.len() {
        let value = args.get(i + 1).ok_or("resume option requires a value")?.parse::<usize>()?;
        match args[i].as_str() {
            "--recovery-solves" if recovery.is_none() && (1..=128).contains(&value) => recovery = Some(value),
            "--pause-after" if pause.is_none() && value <= 200 => pause = Some(value),
            _ => return Err("unknown, repeated or out-of-range resume option; study policy is immutable".into()),
        }
        i += 2;
    }
    let budget = recovery.ok_or("resume requires explicit --recovery-solves (two solves per case)")?;
    let output = Path::new(&args[1]);
    if output.try_exists()? { return Err("output directory already exists; refusing to overwrite it".into()); }
    let bytes = load(Path::new(&args[0]))?;
    let mut remaining = budget;
    let optimizer = MultiLoadProjectedOptimizer::restore_checkpoint(&bytes, &mut remaining)
        .map_err(|error| format!("{error}; recovery_solves_started={}", budget - remaining))?;
    let input = optimizer.geometry().clone();
    run_optimizer(output, optimizer, &input, Options {
        enabled: true, pause_after: pause, recovery_solves: budget - remaining,
    }, None)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pause_is_bounded_opt_in_and_does_not_consume_stress_options() {
        let args = ["out", "loads.csv", "--stress-limit", "100", "--pause-after", "0"]
            .map(str::to_string);
        let (rest, parsed) = options(&args).unwrap();
        assert!(parsed.enabled);
        assert_eq!(parsed.pause_after, Some(0));
        assert_eq!(rest, ["out", "loads.csv", "--stress-limit", "100"]);
        for bad in [vec!["--pause-after"], vec!["--pause-after", "201"],
            vec!["--checkpoint", "--checkpoint"], vec!["--pause-after", "0", "--pause-after", "1"]]
        {
            assert!(options(&bad.into_iter().map(str::to_string).collect::<Vec<_>>()).is_err());
        }
    }
}
