//! fs-substrate — hardware capability probes, machine fingerprints, topology
//! facts, and one-shot SIMD dispatch resolution (plan §5.1, patch Rev Q).
//!
//! The founding rule: hardware facts are MEASURED VALUES, never static
//! claims. Spec-sheet numbers are design inputs at best; the probe is the
//! source of truth, its results persist to the ledger `capability_probes`
//! table, and every performance target is keyed by the resulting
//! fingerprint. GPU / Neural-Engine facts are recorded as SEPARATE optional
//! entries and are NOT consumed by the default pure-CPU build.
//!
//! Fingerprint discipline: the fingerprint hashes STABLE topology facts only
//! (ISA, brand, core counts, cache/page geometry, memory size, feature set).
//! Measured bandwidth lives BESIDE the fingerprint with run-to-run jitter
//! expected — hashing it would make every reboot a "new machine".
//!
//! Probing I/O: on macOS the probe shells out to `sysctl` (no FFI — P1);
//! on Linux it reads `/sys` and `/proc`. Probing is a once-per-startup,
//! latency-lane operation; determinism of COMPUTE is never affected.

use std::sync::OnceLock;
use std::sync::atomic::{AtomicU32, Ordering};

pub mod bandwidth;

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// Cache-line size as a compile-time constant, resolved per target (plan
/// §5.1 consequence 1: 128 B on Apple aarch64, 64 B elsewhere). Getting this
/// wrong silently halves effective bandwidth on contended structures; the
/// probe VERIFIES it against the OS-reported value at runtime.
#[cfg(all(target_arch = "aarch64", target_vendor = "apple"))]
pub const CACHE_LINE: usize = 128;
/// Cache-line size (non-Apple targets).
#[cfg(not(all(target_arch = "aarch64", target_vendor = "apple")))]
pub const CACHE_LINE: usize = 64;

/// Instruction-set families FrankenSim schedules for.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Isa {
    /// Apple Silicon (NEON committed; matrix units probe-gated, never assumed).
    Aarch64Apple,
    /// Other aarch64 (NEON baseline).
    Aarch64Other,
    /// x86-64 (AVX2/AVX-512 runtime-detected).
    X86_64,
    /// Anything else: scalar tier only.
    Other,
}

impl Isa {
    fn name(self) -> &'static str {
        match self {
            Isa::Aarch64Apple => "aarch64-apple",
            Isa::Aarch64Other => "aarch64",
            Isa::X86_64 => "x86_64",
            Isa::Other => "other",
        }
    }
}

/// SIMD execution tier, resolved ONCE at startup into dispatch decisions —
/// no per-call branching in hot paths (plan §5.1 consequence 5).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SimdTier {
    /// Portable scalar (the correctness reference; always available).
    Scalar,
    /// aarch64 NEON (baseline on all aarch64).
    Neon,
    /// x86-64 AVX2.
    Avx2,
    /// x86-64 AVX-512F.
    Avx512,
}

impl SimdTier {
    /// Stable lowercase name (ledger rows, tune-table keys, log fields).
    #[must_use]
    pub fn name(self) -> &'static str {
        match self {
            SimdTier::Scalar => "scalar",
            SimdTier::Neon => "neon",
            SimdTier::Avx2 => "avx2",
            SimdTier::Avx512 => "avx512",
        }
    }
}

/// Measured (jittery) facts — stored BESIDE the fingerprint, never hashed.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Measured {
    /// Single-thread sustained triad bandwidth, GB/s.
    pub single_thread_gbs: f64,
    /// All-logical-core sustained triad bandwidth, GB/s.
    pub all_core_gbs: f64,
}

/// The capability probe: stable topology facts + measured bandwidth +
/// separate optional accelerator facts.
#[derive(Debug, Clone, PartialEq)]
pub struct CapabilityProbe {
    /// ISA family.
    pub isa: Isa,
    /// CPU brand string (e.g. "Apple M4 Pro").
    pub cpu_brand: String,
    /// Detected SIMD/vector features (sorted, deduplicated).
    pub features: Vec<String>,
    /// OS-reported cache-line size in bytes.
    pub cache_line: u32,
    /// OS-reported base page size in bytes (16 KiB on Apple — plan §5.1).
    pub page_size: u32,
    /// Physical memory in bytes.
    pub mem_bytes: u64,
    /// Logical CPU count.
    pub logical_cpus: u32,
    /// Apple performance-core count (perflevel0), when exposed.
    pub perf_cores: Option<u32>,
    /// Apple efficiency-core count (perflevel1), when exposed.
    pub eff_cores: Option<u32>,
    /// Per-cluster L2 size in bytes, when exposed (Apple).
    pub l2_per_cluster: Option<u64>,
    /// NUMA node count (Linux), when exposed.
    pub numa_nodes: Option<u32>,
    /// GPU description — a SEPARATE fact class, never consumed by the
    /// default CPU build (patch Rev Q: no single "accelerator" bucket).
    pub gpu: Option<String>,
    /// Measured bandwidth (excluded from the fingerprint).
    pub measured: Measured,
}

impl CapabilityProbe {
    /// Run the full probe (topology reads + bandwidth measurement).
    /// Bandwidth measurement takes ~100–300 ms by design.
    #[must_use]
    pub fn run() -> CapabilityProbe {
        let mut p = Self::topology_only();
        p.measured = bandwidth::measure(p.logical_cpus as usize);
        p
    }

    /// Topology facts only (fast; no bandwidth measurement).
    #[must_use]
    pub fn topology_only() -> CapabilityProbe {
        let isa = detect_isa();
        let os = OsFacts::read();
        CapabilityProbe {
            isa,
            cpu_brand: os.brand,
            features: detect_features(),
            cache_line: os.cache_line.unwrap_or(CACHE_LINE as u32),
            page_size: os.page_size.unwrap_or(4096),
            mem_bytes: os.mem_bytes.unwrap_or(0),
            logical_cpus: std::thread::available_parallelism().map_or(1, |n| n.get() as u32),
            perf_cores: os.perf_cores,
            eff_cores: os.eff_cores,
            l2_per_cluster: os.l2_per_cluster,
            numa_nodes: os.numa_nodes,
            gpu: os.gpu,
            measured: Measured {
                single_thread_gbs: 0.0,
                all_core_gbs: 0.0,
            },
        }
    }

    /// Stable machine fingerprint: FNV-1a 64 over the STABLE facts only
    /// (measured bandwidth deliberately excluded — see module docs).
    #[must_use]
    pub fn fingerprint(&self) -> u64 {
        let stable = format!(
            "{}|{}|{}|{}|{}|{}|{:?}|{:?}|{:?}|{:?}|{}",
            self.isa.name(),
            self.cpu_brand,
            self.cache_line,
            self.page_size,
            self.mem_bytes,
            self.logical_cpus,
            self.perf_cores,
            self.eff_cores,
            self.l2_per_cluster,
            self.numa_nodes,
            self.features.join(",")
        );
        fs_obs::fnv1a64(stable.as_bytes())
    }

    /// Serialize to the ledger `capability_probes` row shape (canonical JSON,
    /// in-house per P1). `probe_time` is the CALLER's to supply when
    /// inserting — this crate stays clock-free (determinism discipline).
    #[must_use]
    pub fn to_json(&self) -> String {
        use std::fmt::Write as _;
        let mut s = String::with_capacity(512);
        let esc = |t: &str| t.replace('\\', "\\\\").replace('"', "\\\"");
        let _ = write!(
            s,
            "{{\"machine\":\"{:016x}\",\"cpu\":{{\"isa\":\"{}\",\"brand\":\"{}\",\"features\":[",
            self.fingerprint(),
            self.isa.name(),
            esc(&self.cpu_brand)
        );
        for (i, f) in self.features.iter().enumerate() {
            let _ = write!(s, "{}\"{}\"", if i > 0 { "," } else { "" }, esc(f));
        }
        let _ = write!(
            s,
            "]}},\"topology\":{{\"logical_cpus\":{},\"perf_cores\":{},\"eff_cores\":{},\"numa_nodes\":{}}}",
            self.logical_cpus,
            opt_u32(self.perf_cores),
            opt_u32(self.eff_cores),
            opt_u32(self.numa_nodes)
        );
        let _ = write!(
            s,
            ",\"memory\":{{\"bytes\":{},\"page_size\":{},\"cache_line\":{},\"l2_per_cluster\":{}}}",
            self.mem_bytes,
            self.page_size,
            self.cache_line,
            self.l2_per_cluster
                .map_or_else(|| "null".to_string(), |v| v.to_string())
        );
        let _ = write!(
            s,
            ",\"gpu\":{},\"math\":{{\"single_thread_gbs\":{:.3},\"all_core_gbs\":{:.3}}}}}",
            self.gpu
                .as_ref()
                .map_or_else(|| "null".to_string(), |g| format!("\"{}\"", esc(g))),
            self.measured.single_thread_gbs,
            self.measured.all_core_gbs
        );
        s
    }
}

fn opt_u32(v: Option<u32>) -> String {
    v.map_or_else(|| "null".to_string(), |x| x.to_string())
}

fn detect_isa() -> Isa {
    #[cfg(all(target_arch = "aarch64", target_vendor = "apple"))]
    {
        Isa::Aarch64Apple
    }
    #[cfg(all(target_arch = "aarch64", not(target_vendor = "apple")))]
    {
        Isa::Aarch64Other
    }
    #[cfg(target_arch = "x86_64")]
    {
        Isa::X86_64
    }
    #[cfg(not(any(target_arch = "aarch64", target_arch = "x86_64")))]
    {
        Isa::Other
    }
}

/// Runtime feature detection (sorted for fingerprint stability).
fn detect_features() -> Vec<String> {
    let mut f: Vec<String> = Vec::new();
    #[cfg(target_arch = "aarch64")]
    {
        // NEON is architecturally guaranteed on aarch64.
        f.push("neon".to_string());
        if std::arch::is_aarch64_feature_detected!("fp16") {
            f.push("fp16".to_string());
        }
        if std::arch::is_aarch64_feature_detected!("dotprod") {
            f.push("dotprod".to_string());
        }
    }
    #[cfg(target_arch = "x86_64")]
    {
        if std::arch::is_x86_feature_detected!("avx2") {
            f.push("avx2".to_string());
        }
        if std::arch::is_x86_feature_detected!("fma") {
            f.push("fma".to_string());
        }
        if std::arch::is_x86_feature_detected!("avx512f") {
            f.push("avx512f".to_string());
        }
    }
    f.sort();
    f.dedup();
    f
}

// ---------------------------------------------------------------------------
// OS fact reading (macOS: sysctl subprocess; Linux: /sys and /proc).
// ---------------------------------------------------------------------------

#[derive(Debug, Default)]
struct OsFacts {
    brand: String,
    cache_line: Option<u32>,
    page_size: Option<u32>,
    mem_bytes: Option<u64>,
    perf_cores: Option<u32>,
    eff_cores: Option<u32>,
    l2_per_cluster: Option<u64>,
    numa_nodes: Option<u32>,
    gpu: Option<String>,
}

impl OsFacts {
    #[cfg(target_os = "macos")]
    fn read() -> OsFacts {
        let get = |key: &str| -> Option<String> {
            let out = std::process::Command::new("sysctl")
                .arg("-n")
                .arg(key)
                .output()
                .ok()?;
            if !out.status.success() {
                return None;
            }
            let s = String::from_utf8_lossy(&out.stdout).trim().to_string();
            if s.is_empty() { None } else { Some(s) }
        };
        let num = |key: &str| get(key).and_then(|s| s.parse::<u64>().ok());
        OsFacts {
            brand: get("machdep.cpu.brand_string").unwrap_or_else(|| "unknown".to_string()),
            cache_line: num("hw.cachelinesize").map(|v| v as u32),
            page_size: num("hw.pagesize").map(|v| v as u32),
            mem_bytes: num("hw.memsize"),
            perf_cores: num("hw.perflevel0.physicalcpu").map(|v| v as u32),
            eff_cores: num("hw.perflevel1.physicalcpu").map(|v| v as u32),
            l2_per_cluster: num("hw.perflevel0.l2cachesize"),
            numa_nodes: None,
            // Apple GPUs share the SoC brand; recorded as a separate fact,
            // never consumed by the CPU build (Rev Q separation).
            gpu: get("machdep.cpu.brand_string").map(|b| format!("{b} (integrated)")),
        }
    }

    #[cfg(target_os = "linux")]
    fn read() -> OsFacts {
        let read = |p: &str| {
            std::fs::read_to_string(p)
                .ok()
                .map(|s| s.trim().to_string())
        };
        let brand = read("/proc/cpuinfo")
            .and_then(|t| {
                t.lines()
                    .find(|l| l.starts_with("model name"))
                    .and_then(|l| l.split(':').nth(1))
                    .map(|s| s.trim().to_string())
            })
            .unwrap_or_else(|| "unknown".to_string());
        let cache_line = read("/sys/devices/system/cpu/cpu0/cache/index0/coherency_line_size")
            .and_then(|s| s.parse().ok());
        let mem_bytes = read("/proc/meminfo").and_then(|t| {
            t.lines().find(|l| l.starts_with("MemTotal")).and_then(|l| {
                l.split_whitespace()
                    .nth(1)
                    .and_then(|kb| kb.parse::<u64>().ok())
                    .map(|kb| kb * 1024)
            })
        });
        let numa_nodes = std::fs::read_dir("/sys/devices/system/node")
            .ok()
            .map(|rd| {
                rd.flatten()
                    .filter(|e| e.file_name().to_string_lossy().starts_with("node"))
                    .count() as u32
            });
        OsFacts {
            brand,
            cache_line,
            page_size: Some(4096),
            mem_bytes,
            perf_cores: None,
            eff_cores: None,
            l2_per_cluster: None,
            numa_nodes,
            gpu: None,
        }
    }

    #[cfg(not(any(target_os = "macos", target_os = "linux")))]
    fn read() -> OsFacts {
        OsFacts {
            brand: "unknown".to_string(),
            ..OsFacts::default()
        }
    }
}

// ---------------------------------------------------------------------------
// One-shot dispatch resolution.
// ---------------------------------------------------------------------------

static DISPATCH: OnceLock<SimdTier> = OnceLock::new();
static RESOLUTIONS: AtomicU32 = AtomicU32::new(0);

/// The process-wide SIMD tier, resolved exactly once (function-table
/// consumers key off this; no hot-path branching).
pub fn dispatch_tier() -> SimdTier {
    *DISPATCH.get_or_init(|| {
        RESOLUTIONS.fetch_add(1, Ordering::Relaxed);
        resolve_tier()
    })
}

/// Test hook: how many times resolution actually ran (must be ≤ 1).
#[must_use]
pub fn dispatch_resolution_count() -> u32 {
    RESOLUTIONS.load(Ordering::Relaxed)
}

fn resolve_tier() -> SimdTier {
    #[cfg(target_arch = "aarch64")]
    {
        SimdTier::Neon
    }
    #[cfg(target_arch = "x86_64")]
    {
        if std::arch::is_x86_feature_detected!("avx512f") {
            SimdTier::Avx512
        } else if std::arch::is_x86_feature_detected!("avx2") {
            SimdTier::Avx2
        } else {
            SimdTier::Scalar
        }
    }
    #[cfg(not(any(target_arch = "aarch64", target_arch = "x86_64")))]
    {
        SimdTier::Scalar
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn probe_reports_plausible_stable_facts() {
        let p = CapabilityProbe::topology_only();
        assert!(
            matches!(p.cache_line, 64 | 128),
            "cache line {}",
            p.cache_line
        );
        assert!(p.logical_cpus >= 1);
        assert!(
            p.mem_bytes > 1 << 30,
            "at least 1 GiB expected, got {}",
            p.mem_bytes
        );
        assert!(
            p.page_size == 4096 || p.page_size == 16384,
            "page {}",
            p.page_size
        );
        #[cfg(all(target_arch = "aarch64", target_vendor = "apple"))]
        {
            assert_eq!(p.isa, Isa::Aarch64Apple);
            assert!(p.features.iter().any(|f| f == "neon"));
            assert!(p.perf_cores.is_some() && p.eff_cores.is_some());
        }
    }

    #[test]
    fn compile_time_cache_line_matches_os_report() {
        let p = CapabilityProbe::topology_only();
        assert_eq!(
            p.cache_line as usize, CACHE_LINE,
            "the CACHE_LINE const must agree with the OS on this target — if this fires, \
             padding geometry is wrong for this machine (plan §5.1)"
        );
    }

    #[test]
    fn fingerprint_is_stable_and_excludes_measured_bandwidth() {
        let a = CapabilityProbe::topology_only();
        let b = CapabilityProbe::topology_only();
        assert_eq!(
            a.fingerprint(),
            b.fingerprint(),
            "same machine, same fingerprint"
        );
        let mut c = a.clone();
        c.measured = Measured {
            single_thread_gbs: 123.4,
            all_core_gbs: 567.8,
        };
        assert_eq!(
            a.fingerprint(),
            c.fingerprint(),
            "measured jitter must not change identity"
        );
        let mut d = a.clone();
        d.logical_cpus += 1;
        assert_ne!(
            a.fingerprint(),
            d.fingerprint(),
            "topology change must change identity"
        );
    }

    #[test]
    fn dispatch_resolves_exactly_once() {
        let t1 = dispatch_tier();
        let t2 = dispatch_tier();
        assert_eq!(t1, t2);
        assert!(dispatch_resolution_count() <= 1);
        #[cfg(target_arch = "aarch64")]
        assert_eq!(t1, SimdTier::Neon);
        let _ = t1.name();
    }

    #[test]
    fn json_row_has_the_capability_probes_shape() {
        let p = CapabilityProbe::topology_only();
        let j = p.to_json();
        for key in [
            "\"machine\":",
            "\"cpu\":",
            "\"topology\":",
            "\"memory\":",
            "\"gpu\":",
            "\"math\":",
        ] {
            assert!(j.contains(key), "missing {key} in {j}");
        }
        // Balanced braces (single-line, writer-canonical).
        let depth = j
            .chars()
            .fold(0i32, |d, c| d + i32::from(c == '{') - i32::from(c == '}'));
        assert_eq!(depth, 0, "unbalanced JSON: {j}");
        assert!(!j.contains('\n'));
    }
}
