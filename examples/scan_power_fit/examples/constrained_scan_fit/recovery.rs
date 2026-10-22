//! Bounded checkpoint files for the constrained example, not a new optimizer.
use super::{Options, Result};
use fs_exec::Cx;
use std::{fs::{self, File, OpenOptions}, io::{Read, Write}, path::{Path, PathBuf}};
use std::sync::atomic::{AtomicU64, Ordering};

pub const MAX_BYTES: usize = 1 << 20;
#[derive(Default)]
pub struct Recovery {
    pub checkpoint: Option<PathBuf>,
    pub resume: Option<PathBuf>,
    pub pause_after: Option<usize>,
}
impl Recovery {
    pub fn parse(args: Vec<String>) -> Result<(Vec<String>, Self)> {
        let mut out = Self::default();
        let mut rest = Vec::new();
        let mut args = args.into_iter();
        while let Some(arg) = args.next() {
            match arg.as_str() {
                "--checkpoint" if out.checkpoint.is_none() => out.checkpoint = Some(args.next().ok_or("missing checkpoint path")?.into()),
                "--resume" if out.resume.is_none() => out.resume = Some(args.next().ok_or("missing resume path")?.into()),
                "--pause-after" if out.pause_after.is_none() => out.pause_after = Some(args.next().ok_or("missing pause count")?.parse()?),
                "--checkpoint" | "--resume" | "--pause-after" => return Err("duplicate recovery option".into()),
                _ => rest.push(arg),
            }
        }
        if out.pause_after.is_some_and(|n| n == 0) || (out.pause_after.is_some() && out.checkpoint.is_none()) {
            return Err("a positive --pause-after requires --checkpoint FILE".into());
        }
        Ok((rest, out))
    }
    pub fn active(&self) -> bool { self.checkpoint.is_some() || self.resume.is_some() }
    pub fn load(&self) -> Result<Option<Vec<u8>>> {
        self.resume.as_deref().map(|p| {
            let mut bytes = Vec::new();
            File::open(p)?.take((MAX_BYTES + 1) as u64).read_to_end(&mut bytes)?;
            if bytes.len() > MAX_BYTES { return Err("checkpoint exceeds 1 MiB".into()); }
            Ok(bytes)
        }).transpose()
    }
    /// Bind actual executable bytes, physical limits and time discretization. Fixed
    /// model, solver and target settings are compiled into that executable.
    /// Fixed-size blocks make the identity independent of short filesystem reads.
    pub fn model(&self, cx: &Cx<'_>, input: Options) -> Result<[u8; 32]> {
        if !self.active() { return Ok([0;32]); }
        let mut file = File::open(std::env::current_exe()?)?;
        let mut block = vec![0u8; 64 << 10];
        let mut roots = Vec::new();
        let mut size = 0u64;
        loop {
            cx.checkpoint().map_err(|_| "cancelled while binding executable")?;
            let mut n = 0;
            while n < block.len() {
                match file.read(&mut block[n..]) {
                    Ok(0) => break,
                    Ok(k) => n += k,
                    Err(e) if e.kind() == std::io::ErrorKind::Interrupted => continue,
                    Err(e) => return Err(e.into()),
                }
            }
            if n == 0 { break; }
            size += n as u64;
            if size > (1u64 << 30) { return Err("executable fingerprint exceeds its 1 GiB work allowance".into()); }
            roots.extend_from_slice(fs_blake3::hash_domain("scan-fit-executable-block-v1", &block[..n]).as_bytes());
            if n < block.len() { break; }
        }
        roots.extend_from_slice(&size.to_le_bytes());
        roots.extend_from_slice(&input.energy_j.to_bits().to_le_bytes());
        roots.extend_from_slice(&input.peak_fraction.to_bits().to_le_bytes());
        roots.push(u8::from(input.power_slew.is_some()));
        if let Some(rates) = input.power_slew {
            for rate in rates { roots.extend_from_slice(&rate.to_bits().to_le_bytes()); }
            roots.extend_from_slice(&(input.ramp_steps as u64).to_le_bytes());
        }
        Ok(*fs_blake3::hash_domain("constrained-scan-fit-model-v1", &roots).as_bytes())
    }
    pub fn save(&self, bytes: &[u8]) -> Result<()> {
        if let Some(path) = &self.checkpoint { atomic_write(path, bytes)?; }
        Ok(())
    }
}

fn atomic_write(path: &Path, bytes: &[u8]) -> Result<()> {
    static SEQUENCE: AtomicU64 = AtomicU64::new(0);
    if bytes.len() > MAX_BYTES { return Err("checkpoint exceeds 1 MiB".into()); }
    let parent = path.parent().filter(|p| !p.as_os_str().is_empty()).unwrap_or(Path::new("."));
    let name = path.file_name().ok_or("checkpoint must name a file")?;
    let mut temporary_name = name.to_os_string();
    temporary_name.push(format!(".tmp-{}-{}", std::process::id(), SEQUENCE.fetch_add(1, Ordering::Relaxed)));
    let temporary = parent.join(temporary_name);
    let mut file = OpenOptions::new().write(true).create_new(true).open(&temporary)?;
    let result = (|| -> Result<()> {
        file.write_all(bytes)?;
        file.sync_all()?;
        drop(file);
        fs::rename(&temporary, path)?;
        #[cfg(unix)]
        File::open(parent)?.sync_all()?;
        Ok(())
    })();
    if result.is_err() { let _ = fs::remove_file(&temporary); }
    result
}
