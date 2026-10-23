//! The CUDA backend: uploads a `Plan` and a batch of `SimState`s and runs kernel.cu, many steps per launch.
//!
//! Optional twice over: built only with the `cuda` feature, and at run time the driver and NVRTC are loaded
//! dynamically (cudarc), so a binary built with CUDA support still runs, on the CPU, on a machine without an
//! NVIDIA driver. The kernel is compiled at start-up for the GPU that is present, straight to its machine
//! code (CUBIN), so a runtime compiler newer than the driver is fine.

use crate::plan::{Plan, SimState, TPB};
use cudarc::driver::{CudaContext, CudaSlice, CudaStream, DevicePtr, DeviceRepr, LaunchConfig, PushKernelArg, sys};
use cudarc::nvrtc::{Ptx, result, sys as nsys};
use std::sync::Arc;

/// Names and C types of the kernel's parameter struct, in order. Every field is 8 bytes on both sides
/// (pointers, `long long`, `double`), so the Rust struct below and the generated C struct share one layout.
macro_rules! params {
    ($($name:ident : $cty:literal),* $(,)?) => {
        #[repr(C)]
        #[derive(Clone, Copy, Default)]
        #[allow(non_snake_case)]
        pub struct Params { $(pub $name: u64,)* }
        const PARAMS_C: &[(&str, &str)] = &[$((stringify!($name), $cty),)*];
    };
}

params! {
    B: "ll", NT: "ll", NB: "ll", NC: "ll", NI: "ll", NCI: "ll", NMD: "ll", NM: "ll", NCON: "ll", NCB: "ll",
    NFR: "ll", NFB: "ll", NCAV: "ll", NRAD: "ll", NLINK: "ll", NDRV: "ll", NVK: "ll", RINGTOT: "ll",
    NROWS: "ll", NCN: "ll", NSEG: "ll", KL: "ll", KS: "ll", KG: "ll", steps: "ll", n0: "ll", gstep0: "ll", k: "double", rate: "double",
    gas_density: "double", gas_bulk_modulus: "double",
    bias: "double*", ia: "double*", cw: "double*", cwp: "double*", cs1: "double*", ms0: "double*", mass: "double*", z0: "double*",
    lapE_col: "int*", lapE_val: "double*", stfE_col: "int*", stfE_val: "double*", gtE_core: "int*", gtE_slot: "int*",
    visc_ptr: "int*", visc_col: "int*", visc_val: "double*", vroot_ptr: "int*", vroot_col: "int*", vroot_val: "double*", bviscous: "int*",
    bstiff: "int*",
    bh: "double*", bih2: "double*", bs1: "double*", bten: "double*", bber: "double*", bairm: "double*", barea: "double*", beh: "double*",
    bvk: "int*", bcoff: "int*", bnc: "int*", bitem0: "int*", bitem1: "int*", bcitem0: "int*", bcitem1: "int*", bcav: "int*", brad: "int*",
    cnb: "int*", crest: "double*", vk_body: "int*", airy_off: "ll*", airy64: "double*", airy32: "float*",
    item_body: "int*", item_start: "int*", item_len: "int*", citem_body: "int*", citem_start: "int*", citem_len: "int*",
    W0: "double*", W1: "double*", W2: "double*", LW0: "double*", LW1: "double*", F: "double*", G: "double*", TMP: "double*",
    Q: "double*", QF: "float*", CXX: "double*", CYY: "double*", CXY: "double*", AF: "double*",
    psiA: "double*", rprevA: "double*", budgetA: "double*", psi1: "double*", rprev1: "double*", savf: "double*",
    viscous_dissipated: "double*",
    p1: "double*", p3: "double*", pc3: "double*", p4: "double*", p7: "double*",
    cav_hptr: "int*", cav_head: "int*", cav_vol: "double*", cav_port: "int*", port_area: "double*", port_len: "double*",
    port_res: "double*", port_x: "double*", port_u: "double*", port_du: "double*",
    rad_ptr: "int*", rad_body: "int*", rad_res: "double*",
    lk_a: "int*", lk_b: "int*", lk_k: "double*", lk_c: "double*", lnk_ptr: "int*", lnk_idx: "int*", LF: "double*",
    dv_ptr: "int*", dv_t: "double*", dv_v: "double*", drv_ptr: "int*", drv_idx: "int*",
    cnode: "int*", cn_node: "int*", cn_body: "int*",
    ct_a0: "int*", ct_b0: "int*", ct_end: "int*", ct_node: "int*", ct_w: "double*", ct_gap: "double*", ct_kk: "double*",
    ct_alpha: "double*", ct_mu: "double*", ct_comp: "double*", ct_fx: "int*", cb_ptr: "int*", cb_idx: "int*",
    ct_active: "int*", ct_force: "double*", ct_dnow: "double*", ct_peak: "double*",
    bc_ptr: "int*", bc_idx: "int*",
    fr_x0: "int*", fr_y0: "int*", fr_end: "int*", fr_node: "int*", fr_w: "double*", fr_contact: "int*", fr_mu: "double*",
    fr_cxx: "double*", fr_cyy: "double*", fr_cxy: "double*", fb_ptr: "int*", fb_idx: "int*", fr_force: "double*",
    st_ptr: "int*", st_body: "int*", st_contact: "int*", st_time: "double*", st_speed: "double*", st_clear: "double*", st_cursor: "int*",
    seg_q0: "int*", seg_part: "double*", md_seg: "int*", md_node: "int*", md_w: "double*", md_mic: "int*", md_delay: "int*",
    mic_roff: "int*", mic_rmask: "int*",
    ring: "double*", pt_ptr: "int*", pt_cav: "int*", pt_delay: "int*", pt_lo: "double*", pt_hi: "double*", mic_tau: "double*",
    mic_hp: "double*", out: "double*", bar: "unsigned int*", prof: "ll*",
}

unsafe impl DeviceRepr for Params {}

/// Simulations per pass over a dense Airy inverse (more with f32 storage, where the pass is bound by
/// reading the matrix, not by FP64 arithmetic); the batch is padded to a multiple of it.
fn gb(airy_f32: bool) -> usize {
    if airy_f32 { 16 } else { 4 }
}

fn source(airy_f32: bool) -> String {
    let mut s = String::from("typedef long long ll;\nstruct P {\n");
    for (name, ty) in PARAMS_C {
        s += &format!("  {ty} {name};\n");
    }
    s += "};\n";
    s += &format!(
        "#define GB {}\n#define AIRY_F32 {}\n#define SPIN_NS {}\n",
        gb(airy_f32), airy_f32 as i32,
        std::env::var("RENDER_GPU_SPIN_NS").ok().and_then(|v| v.parse::<u32>().ok()).unwrap_or(32)
    );
    s += include_str!("kernel.cu");
    s
}

/// The compiled kernel for this source and GPU architecture, from the on-disk cache when present
/// (compiling takes seconds on an idle host and over a minute on a loaded one).
fn compile(src: &str, arch: &str) -> Result<Ptx, String> {
    use std::hash::{Hash, Hasher};
    let opts = [format!("--gpu-architecture={arch}"), "-std=c++17".into(), "--fmad=true".into()];
    let mut h = std::collections::hash_map::DefaultHasher::new();
    (src, &opts).hash(&mut h);
    let dir = std::env::var_os("XDG_CACHE_HOME")
        .map(std::path::PathBuf::from)
        .or_else(|| std::env::var_os("HOME").map(|h| std::path::PathBuf::from(h).join(".cache")))
        .map(|d| d.join("fs-render-gpu"));
    let file = dir.as_ref().map(|d| d.join(format!("{:016x}-{arch}.cubin", h.finish())));
    if let Some(Ok(bin)) = file.as_ref().map(std::fs::read) {
        return Ok(Ptx::from_binary(bin));
    }
    let bin = nvrtc_cubin(src, &opts)?;
    if let (Some(d), Some(f)) = (dir, file) {
        // best effort (a read-only home only costs the next run a compile), and atomic: written aside and
        // renamed, so a second process (`--device all` compiles once per GPU) never loads half a file
        let tmp = d.join(format!(".{}.{}", f.file_name().map_or(String::new(), |n| n.to_string_lossy().into_owned()), std::process::id()));
        let _ = std::fs::create_dir_all(&d).and_then(|_| std::fs::write(&tmp, &bin)).and_then(|_| std::fs::rename(&tmp, &f));
    }
    Ok(Ptx::from_binary(bin))
}

/// Compile to the device's own machine code (CUBIN); the driver never JIT-compiles PTX.
fn nvrtc_cubin(src: &str, opts: &[String]) -> Result<Vec<u8>, String> {
    let c = std::ffi::CString::new(src).map_err(|e| e.to_string())?;
    let prog = result::create_program(&c, None).map_err(|e| format!("nvrtc: {e:?}"))?;
    if let Err(e) = unsafe { result::compile_program(prog, opts) } {
        let log = unsafe { result::get_program_log(prog) }.map_err(|e| format!("{e:?}"))?;
        let log = unsafe { std::ffi::CStr::from_ptr(log.as_ptr()) }.to_string_lossy().into_owned();
        return Err(format!("nvrtc {e:?}:\n{log}"));
    }
    let mut size = 0usize;
    unsafe { nsys::nvrtcGetCUBINSize(prog, &mut size) }.result().map_err(|e| format!("{e:?}"))?;
    let mut buf = vec![0u8; size];
    unsafe { nsys::nvrtcGetCUBIN(prog, buf.as_mut_ptr().cast()) }.result().map_err(|e| format!("{e:?}"))?;
    unsafe { result::destroy_program(prog) }.map_err(|e| format!("{e:?}"))?;
    Ok(buf)
}

/// Compile both kernel variants (f64 and f32 Airy storage) for `arch` into the on-disk cache, without a
/// GPU: NVRTC needs none. For hosts where a compile is slow (a busy build machine), compile elsewhere and
/// copy ~/.cache/fs-render-gpu over.
pub fn precompile(arch: &str) -> Result<(), String> {
    let arch: &'static str = Box::leak(arch.to_string().into_boxed_str());
    for f32_airy in [false, true] {
        let t = std::time::Instant::now();
        compile(&source(f32_airy), arch)?;
        eprintln!("{arch} airy {}: {:.1} s", if f32_airy { "f32" } else { "f64" }, t.elapsed().as_secs_f64());
    }
    Ok(())
}

/// Number of CUDA devices (0 without a driver).
pub fn device_count() -> usize {
    std::panic::catch_unwind(|| CudaContext::device_count().ok()).ok().flatten().map_or(0, |n| n.max(0) as usize)
}

/// The first CUDA device's name, if a driver and a device are present. Never panics: a missing driver
/// library is reported as `None`.
pub fn probe(device: usize) -> Option<String> {
    std::panic::catch_unwind(|| {
        let ctx = CudaContext::new(device).ok()?;
        ctx.name().ok()
    })
    .ok()
    .flatten()
}

fn e<T: std::fmt::Debug>(x: T) -> String {
    format!("{x:?}")
}

/// Device memory for one run; every allocation stays alive while the kernel can use it.
struct Mem {
    stream: Arc<CudaStream>,
    f64s: Vec<CudaSlice<f64>>,
    f32s: Vec<CudaSlice<f32>>,
    i32s: Vec<CudaSlice<i32>>,
    i64s: Vec<CudaSlice<i64>>,
    u32s: Vec<CudaSlice<u32>>,
}

impl Mem {
    fn ptr<T: DeviceRepr>(&self, s: &CudaSlice<T>) -> u64 {
        let (p, _sync) = s.device_ptr(&self.stream);
        p
    }
    fn f64(&mut self, v: &[f64]) -> Result<u64, String> {
        let s = if v.is_empty() { self.stream.alloc_zeros::<f64>(1) } else { self.stream.clone_htod(v) }.map_err(e)?;
        let p = self.ptr(&s);
        self.f64s.push(s);
        Ok(p)
    }
    fn f32(&mut self, v: &[f32]) -> Result<u64, String> {
        let s = if v.is_empty() { self.stream.alloc_zeros::<f32>(1) } else { self.stream.clone_htod(v) }.map_err(e)?;
        let p = self.ptr(&s);
        self.f32s.push(s);
        Ok(p)
    }
    fn i32(&mut self, v: &[i32]) -> Result<u64, String> {
        let s = if v.is_empty() { self.stream.alloc_zeros::<i32>(1) } else { self.stream.clone_htod(v) }.map_err(e)?;
        let p = self.ptr(&s);
        self.i32s.push(s);
        Ok(p)
    }
    fn i64(&mut self, v: &[i64]) -> Result<u64, String> {
        let s = if v.is_empty() { self.stream.alloc_zeros::<i64>(1) } else { self.stream.clone_htod(v) }.map_err(e)?;
        let p = self.ptr(&s);
        self.i64s.push(s);
        Ok(p)
    }
    fn zeros(&mut self, n: usize) -> Result<u64, String> {
        let s = self.stream.alloc_zeros::<f64>(n.max(1)).map_err(e)?;
        let p = self.ptr(&s);
        self.f64s.push(s);
        Ok(p)
    }
}

pub struct Run {
    pub device_name: String,
    pub blocks: u32,
    pub compile_ms: f64,
    /// Seconds spent stepping (the launches and their output copies), setup excluded.
    pub step_s: f64,
    /// Accumulated material curvature-rate loss [J], laid out [sim][body].
    /// Includes each input state's retained counter; other losses are separate.
    pub bending_loss_j: Vec<f64>,
}

/// Step every simulation of the batch `total` times; `on_block(first_step, samples)` receives the
/// microphone samples of each launch, laid out [step][sim][mic].
pub fn run(
    plan: &Plan,
    sims: &[SimState],
    total: usize,
    airy_f32: bool,
    device: usize,
    steps_per_launch: usize,
    mut on_block: impl FnMut(usize, &[f64]),
) -> Result<Run, String> {
    let b = sims.len();
    let ctx = CudaContext::new(device).map_err(e)?;
    let stream = ctx.default_stream();
    let name = ctx.name().map_err(e)?;
    let major = ctx.attribute(sys::CUdevice_attribute::CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MAJOR).map_err(e)?;
    let minor = ctx.attribute(sys::CUdevice_attribute::CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MINOR).map_err(e)?;
    let sms = ctx.attribute(sys::CUdevice_attribute::CU_DEVICE_ATTRIBUTE_MULTIPROCESSOR_COUNT).map_err(e)? as u32;
    let t0 = std::time::Instant::now();
    let arch: &'static str = Box::leak(format!("sm_{major}{minor}").into_boxed_str());
    let module = ctx.load_module(compile(&source(airy_f32), arch)?).map_err(e)?;
    let func = module.load_function("vk_step").map_err(e)?;
    let compile_ms = t0.elapsed().as_secs_f64() * 1e3;
    let mut m = Mem { stream: stream.clone(), f64s: vec![], f32s: vec![], i32s: vec![], i64s: vec![], u32s: vec![] };
    let cat = |f: &dyn Fn(&SimState) -> &Vec<f64>| sims.iter().flat_map(|s| f(s).iter().copied()).collect::<Vec<f64>>();
    let cati = |f: &dyn Fn(&SimState) -> &Vec<i32>| sims.iter().flat_map(|s| f(s).iter().copied()).collect::<Vec<i32>>();
    let (nt, nc, nb, ni, nci) = (plan.nt, plan.nc, plan.nb, plan.item_body.len(), plan.citem_body.len());
    let ncon = plan.ct_gap.len();
    let mut p = Params {
        B: b as u64, NT: nt as u64, NB: nb as u64, NC: nc as u64, NI: ni as u64, NCI: nci as u64,
        NMD: plan.md_mic.len() as u64, NM: plan.nm as u64, NCON: ncon as u64, NCB: (plan.cb_ptr.len().max(1) - 1) as u64,
        NFR: plan.fr_mu.len() as u64, NFB: (plan.fb_ptr.len().max(1) - 1) as u64, NCAV: plan.cav_vol.len() as u64,
        NRAD: plan.rad_res.len() as u64, NLINK: plan.lk_k.len() as u64, NDRV: plan.dv_node.len() as u64,
        NVK: plan.vk_body.len() as u64, RINGTOT: plan.ringtot as u64,
        k: plan.k.to_bits(), rate: plan.rate.to_bits(),
        gas_density: plan.gas_density.to_bits(), gas_bulk_modulus: plan.gas_bulk_modulus.to_bits(),
        ..Default::default()
    };
    p.bias = m.f64(&plan.bias)?; p.ia = m.f64(&plan.ia)?; p.cw = m.f64(&plan.cw)?; p.cwp = m.f64(&plan.cwp)?;
    p.cs1 = m.f64(&plan.cs1)?; p.ms0 = m.f64(&plan.ms0)?; p.mass = m.f64(&plan.mass)?; p.z0 = m.f64(&plan.z0)?;
    p.lapE_col = m.i32(&plan.lap_ecol)?; p.lapE_val = m.f64(&plan.lap_eval)?;
    p.stfE_col = m.i32(&plan.stf_ecol)?; p.stfE_val = m.f64(&plan.stf_eval)?;
    p.visc_ptr = m.i32(&plan.visc_ptr)?; p.visc_col = m.i32(&plan.visc_col)?; p.visc_val = m.f64(&plan.visc_val)?;
    p.vroot_ptr = m.i32(&plan.vroot_ptr)?; p.vroot_col = m.i32(&plan.vroot_col)?; p.vroot_val = m.f64(&plan.vroot_val)?;
    p.bviscous = m.i32(&plan.bviscous)?;
    p.gtE_core = m.i32(&plan.gt_ecore)?; p.gtE_slot = m.i32(&plan.gt_eslot)?; p.bstiff = m.i32(&plan.bstiff)?;
    p.NROWS = plan.bnc.iter().map(|&n| n as u64).sum();
    (p.KL, p.KS, p.KG) = (plan.kl as u64, plan.ks as u64, plan.kg as u64);
    p.bh = m.f64(&plan.bh)?;
    let bih2: Vec<f64> = plan.bh.iter().map(|&h| 1.0 / (h * h)).collect();
    p.bih2 = m.f64(&bih2)?; p.bs1 = m.f64(&plan.bs1)?; p.bten = m.f64(&plan.bten)?; p.bber = m.f64(&plan.bber)?;
    p.bairm = m.f64(&plan.bairm)?; p.barea = m.f64(&plan.barea)?; p.beh = m.f64(&plan.beh)?;
    p.bvk = m.i32(&plan.bvk)?; p.bcoff = m.i32(&plan.bcoff)?; p.bnc = m.i32(&plan.bnc)?;
    p.bitem0 = m.i32(&plan.bitem0)?; p.bitem1 = m.i32(&plan.bitem1)?; p.bcitem0 = m.i32(&plan.bcitem0)?;
    p.bcitem1 = m.i32(&plan.bcitem1)?; p.bcav = m.i32(&plan.bcav)?; p.brad = m.i32(&plan.brad)?;
    p.cnb = m.i32(&plan.cnb)?; p.crest = m.f64(&plan.crest)?; p.vk_body = m.i32(&plan.vk_body)?; p.airy_off = m.i64(&plan.airy_off)?;
    if airy_f32 {
        let a32: Vec<f32> = plan.airy.iter().map(|&x| x as f32).collect();
        p.airy32 = m.f32(&a32)?;
        p.airy64 = m.f64(&[])?;
    } else {
        p.airy64 = m.f64(&plan.airy)?;
        p.airy32 = m.f32(&[])?;
    }
    p.item_body = m.i32(&plan.item_body)?; p.item_start = m.i32(&plan.item_start)?; p.item_len = m.i32(&plan.item_len)?;
    p.citem_body = m.i32(&plan.citem_body)?; p.citem_start = m.i32(&plan.citem_start)?; p.citem_len = m.i32(&plan.citem_len)?;
    // state: W[n%3] = wp, W[(n+1)%3] = w at the first step (n = 0); L wp lives in LW[1]
    p.W0 = m.f64(&cat(&|s| &s.wp))?; p.W1 = m.f64(&cat(&|s| &s.w))?; p.W2 = m.zeros(b * nt)?;
    p.LW0 = m.zeros(b * nt)?; p.LW1 = m.f64(&cat(&|s| &s.lwp))?;
    p.F = m.zeros(b * nt)?; p.G = m.zeros(b * nt)?; p.TMP = m.zeros(b * nt)?;
    let bp = b.next_multiple_of(gb(airy_f32));
    p.Q = m.zeros(bp * nc)?; p.QF = m.f32(&vec![0.0f32; (bp * nc).max(1)])?;
    p.CXX = m.zeros(b * nc)?; p.CYY = m.zeros(b * nc)?; p.CXY = m.zeros(b * nc)?; p.AF = m.zeros(bp * nc)?;
    p.psiA = m.f64(&cat(&|s| &s.psi))?; p.rprevA = m.f64(&cat(&|s| &s.rprev))?; p.budgetA = m.f64(&cat(&|s| &s.budget))?;
    let viscous_counter_index = m.f64s.len();
    p.viscous_dissipated = m.f64(&cat(&|s| &s.viscous_dissipated))?;
    p.psi1 = m.zeros(b * nb)?; p.rprev1 = m.zeros(b * nb)?; p.savf = m.zeros(b * nb)?;
    p.p1 = m.zeros(b * ni * 4)?; p.p3 = m.zeros(b * ni)?; p.pc3 = m.zeros(b * nci)?; p.p4 = m.zeros(b * ni * 5)?; p.p7 = m.zeros(b * ni * 4)?;
    p.cav_hptr = m.i32(&plan.cav_hptr)?; p.cav_head = m.i32(&plan.cav_head)?; p.cav_vol = m.f64(&plan.cav_vol)?;
    p.cav_port = m.i32(&plan.cav_port)?; p.port_area = m.f64(&plan.port_area)?; p.port_len = m.f64(&plan.port_len)?;
    // port state is double-buffered by step parity: [0] holds the state the first step reads
    let twice = |v: Vec<f64>| v.iter().copied().chain(v.iter().copied()).collect::<Vec<f64>>();
    p.port_res = m.f64(&plan.port_res)?; p.port_x = m.f64(&twice(cat(&|s| &s.port_x)))?;
    p.port_u = m.f64(&twice(cat(&|s| &s.port_u)))?; p.port_du = m.f64(&twice(cat(&|s| &s.port_du)))?;
    p.rad_ptr = m.i32(&plan.rad_ptr)?; p.rad_body = m.i32(&plan.rad_body)?; p.rad_res = m.f64(&plan.rad_res)?;
    p.lk_a = m.i32(&plan.lk_a)?; p.lk_b = m.i32(&plan.lk_b)?; p.lk_k = m.f64(&plan.lk_k)?; p.lk_c = m.f64(&plan.lk_c)?;
    p.lnk_ptr = m.i32(&plan.lnk_ptr)?; p.lnk_idx = m.i32(&plan.lnk_idx)?; p.LF = m.zeros(b * plan.lk_k.len())?;
    p.dv_ptr = m.i32(&plan.dv_ptr)?; p.dv_t = m.f64(&plan.dv_t)?; p.dv_v = m.f64(&plan.dv_v)?;
    p.drv_ptr = m.i32(&plan.drv_ptr)?; p.drv_idx = m.i32(&plan.drv_idx)?;
    p.cnode = m.i32(&plan.cnode)?; p.cn_node = m.i32(&plan.cn_node)?; p.cn_body = m.i32(&plan.cn_body)?;
    p.NCN = plan.cn_node.len() as u64;
    p.ct_a0 = m.i32(&plan.ct_a0)?; p.ct_b0 = m.i32(&plan.ct_b0)?; p.ct_end = m.i32(&plan.ct_end)?; p.ct_node = m.i32(&plan.ct_node)?;
    p.ct_w = m.f64(&plan.ct_w)?; p.ct_gap = m.f64(&plan.ct_gap)?; p.ct_kk = m.f64(&plan.ct_kk)?; p.ct_alpha = m.f64(&plan.ct_alpha)?;
    p.ct_mu = m.f64(&plan.ct_mu)?; p.ct_comp = m.f64(&plan.ct_comp)?; p.ct_fx = m.i32(&plan.ct_fx)?;
    p.cb_ptr = m.i32(&plan.cb_ptr)?; p.cb_idx = m.i32(&plan.cb_idx)?;
    p.ct_active = m.i32(&cati(&|s| &s.ct_active))?; p.ct_force = m.f64(&cat(&|s| &s.ct_force))?;
    p.ct_dnow = m.f64(&cat(&|s| &s.ct_dnow))?; p.ct_peak = m.f64(&cat(&|s| &s.ct_peak))?;
    p.bc_ptr = m.i32(&plan.bc_ptr)?; p.bc_idx = m.i32(&plan.bc_idx)?;
    p.fr_x0 = m.i32(&plan.fr_x0)?; p.fr_y0 = m.i32(&plan.fr_y0)?; p.fr_end = m.i32(&plan.fr_end)?; p.fr_node = m.i32(&plan.fr_node)?;
    p.fr_w = m.f64(&plan.fr_w)?; p.fr_contact = m.i32(&plan.fr_contact)?; p.fr_mu = m.f64(&plan.fr_mu)?;
    p.fr_cxx = m.f64(&plan.fr_cxx)?; p.fr_cyy = m.f64(&plan.fr_cyy)?; p.fr_cxy = m.f64(&plan.fr_cxy)?;
    p.fb_ptr = m.i32(&plan.fb_ptr)?; p.fb_idx = m.i32(&plan.fb_idx)?; p.fr_force = m.zeros(b * plan.fr_mu.len())?;
    let mut st_ptr = vec![0i32];
    let (mut st_body, mut st_contact, mut st_time, mut st_speed, mut st_clear) = (vec![], vec![], vec![], vec![], vec![]);
    for s in sims {
        for &(t, body, c, v, clear) in &s.strikes {
            st_time.push(t);
            st_body.push(body);
            st_contact.push(c);
            st_speed.push(v);
            st_clear.push(clear);
        }
        st_ptr.push(st_time.len() as i32);
    }
    let cursor: Vec<i32> = st_ptr[..b].to_vec();
    p.st_ptr = m.i32(&st_ptr)?; p.st_body = m.i32(&st_body)?; p.st_contact = m.i32(&st_contact)?;
    p.st_time = m.f64(&st_time)?; p.st_speed = m.f64(&st_speed)?; p.st_clear = m.f64(&st_clear)?; p.st_cursor = m.i32(&cursor)?;
    p.NSEG = (plan.seg_q0.len() - 1) as u64;
    p.seg_q0 = m.i32(&plan.seg_q0)?; p.seg_part = m.zeros(b * (plan.seg_q0.len() - 1))?; p.md_seg = m.i32(&plan.md_seg)?;
    p.md_node = m.i32(&plan.md_node)?; p.md_w = m.f64(&plan.md_w)?;
    p.md_mic = m.i32(&plan.md_mic)?; p.md_delay = m.i32(&plan.md_delay)?; p.mic_roff = m.i32(&plan.mic_roff)?;
    p.mic_rmask = m.i32(&plan.mic_rmask)?; p.ring = m.zeros(b * plan.ringtot)?;
    p.pt_ptr = m.i32(&plan.pt_ptr)?; p.pt_cav = m.i32(&plan.pt_cav)?; p.pt_delay = m.i32(&plan.pt_delay)?;
    p.pt_lo = m.f64(&plan.pt_lo)?; p.pt_hi = m.f64(&plan.pt_hi)?; p.mic_tau = m.f64(&plan.mic_tau)?;
    p.mic_hp = m.f64(&cat(&|s| &s.mic_hp))?;
    // a GPU that also drives a display kills kernels that run for seconds: keep each launch short there
    let watchdog = ctx.attribute(sys::CUdevice_attribute::CU_DEVICE_ATTRIBUTE_KERNEL_EXEC_TIMEOUT).map_err(e)? != 0;
    let spl = if watchdog { steps_per_launch.clamp(1, 2000) } else { steps_per_launch.max(1) };
    let out = stream.alloc_zeros::<f64>(spl * b * plan.nm.max(1)).map_err(e)?;
    p.out = m.ptr(&out);
    let bar = stream.alloc_zeros::<u32>(64).map_err(e)?;
    p.bar = m.ptr(&bar);
    m.u32s.push(bar);
    let prof = stream.alloc_zeros::<i64>(16).map_err(e)?;
    p.prof = m.ptr(&prof);
    // one block per multiprocessor (measured on an RTX 4090: two per multiprocessor cost more in barriers
    // than they hid in latency)
    let per_sm = std::env::var("RENDER_GPU_BLOCKS_PER_SM").ok().and_then(|v| v.parse::<u32>().ok()).unwrap_or(1);
    let blocks = sms * per_sm;
    let cfg = LaunchConfig { grid_dim: (blocks, 1, 1), block_dim: (TPB as u32, 1, 1), shared_mem_bytes: 0 };
    let mut done = 0usize;
    let t_step = std::time::Instant::now();
    let mut host = vec![0.0f64; spl * b * plan.nm.max(1)];
    while done < total {
        let steps = spl.min(total - done);
        p.steps = steps as u64;
        p.n0 = done as u64;
        p.gstep0 = done as u64;
        let mut l = stream.launch_builder(&func);
        l.arg(&p);
        if let Err(err) = unsafe { l.launch_cooperative(cfg) } {
            return Err(format!("cooperative launch of {blocks} blocks ({per_sm} per multiprocessor): {err:?}"));
        }
        stream.memcpy_dtoh(&out, &mut host).map_err(e)?;
        stream.synchronize().map_err(e)?;
        on_block(done, &host[..steps * b * plan.nm]);
        done += steps;
    }
    let pr: Vec<i64> = stream.clone_dtoh(&prof).map_err(e)?;
    if std::env::var("RENDER_GPU_PROF").is_ok() {
        let tot: i64 = pr[..8].iter().sum();
        eprintln!("cycles per step by phase (block 0): {:?} total {}", pr.iter().take(7).map(|c| c / total.max(1) as i64).collect::<Vec<_>>(), tot / total.max(1) as i64);
    }
    let mut bending_loss_j = stream.clone_dtoh(&m.f64s[viscous_counter_index]).map_err(e)?;
    bending_loss_j.truncate(b * nb);
    drop(m);
    Ok(Run { device_name: name, blocks, compile_ms, step_s: t_step.elapsed().as_secs_f64(), bending_loss_j })
}
