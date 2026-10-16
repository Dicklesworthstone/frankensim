//! The kernel on an NVIDIA GPU (CUDA through cudarc: the driver and NVRTC are loaded at run time, so a
//! machine without them renders on the CPU). One cooperative launch per stretch: the driver guarantees
//! every block is resident, so the grid barrier is safe. `real` is `double` by default; FS_PIANO_GPU_REAL=
//! df64 selects double-float pairs on the FP32 units (GeForce cards run FP64 at 1/64 rate). Either way a
//! real is 8 bytes, so host buffers are u64 words in the chosen encoding. The kernel is compiled for the
//! GPU found, straight to its machine code (CUBIN), and cached in ~/.cache/fs-render-gpu-piano.
use crate::accel::{Backend, SEGMENT};
use crate::pack::{self, Packed};
use crate::plan::Program;
use cudarc::driver::{
    CudaContext, CudaFunction, CudaSlice, CudaStream, LaunchConfig, PushKernelArg, sys,
};
use cudarc::nvrtc::{Ptx, result, sys as nsys};
use std::sync::Arc;

const PRELUDE_COMMON: &str = r#"
#define PC const
#define CG const
#define CI const
#define SG
#define OG
#define TGM
#define TH
#define UNROLL _Pragma("unroll 4")
#define UNROLL_ALL _Pragma("unroll")
#define BLOCK_SYNC() __syncthreads()
#define BLOCK_SYNC_TG() __syncthreads()
"#;

const PRELUDE_DOUBLE: &str = r#"
typedef double real;
#define R(x) (x)
__device__ __forceinline__ real shfl(real x, int i) { return __shfl_sync(0xffffffffu, x, i); }
__device__ __forceinline__ real warp_sum(real x) {
  for (int m = 16; m >= 1; m >>= 1) x += __shfl_xor_sync(0xffffffffu, x, m);
  return x;
}
"#;

const PRELUDE_DF64: &str = r#"
struct df { float hi; float lo; };
__device__ __forceinline__ df mk(float h, float l) { df x; x.hi = h; x.lo = l; return x; }
__device__ __forceinline__ df quick_two_sum(float a, float b) { float s = a + b; return mk(s, b - (s - a)); }
__device__ __forceinline__ df two_sum(float a, float b) { float s = a + b; float v = s - a; return mk(s, (a - (s - v)) + (b - v)); }
__device__ __forceinline__ df operator+(df a, df b) { df s = two_sum(a.hi, b.hi); return quick_two_sum(s.hi, s.lo + (a.lo + b.lo)); }
__device__ __forceinline__ df operator-(df a) { return mk(-a.hi, -a.lo); }
__device__ __forceinline__ df operator-(df a, df b) { return a + (-b); }
__device__ __forceinline__ df operator*(df a, df b) {
  float p = a.hi * b.hi;
  float e = fmaf(a.hi, b.hi, -p);
  e = fmaf(a.hi, b.lo, e);
  e = fmaf(a.lo, b.hi, e);
  return quick_two_sum(p, e);
}
typedef df real;
#define R(x) mk(x##f, 0.0f)
__device__ __forceinline__ real shfl(real x, int i) { return mk(__shfl_sync(0xffffffffu, x.hi, i), __shfl_sync(0xffffffffu, x.lo, i)); }
__device__ __forceinline__ real warp_sum(real x) {
  for (int m = 16; m >= 1; m >>= 1) x = x + mk(__shfl_xor_sync(0xffffffffu, x.hi, m), __shfl_xor_sync(0xffffffffu, x.lo, m));
  return x;
}
"#;

const PRELUDE_TAIL: &str = r#"
__device__ __forceinline__ void block_sum(real* v, int nv, real* sh, int tid) {
  const int lane = tid & 31, warp = tid >> 5;
  for (int i = 0; i < nv; i++) { real w = warp_sum(v[i]); if (lane == 0) sh[i * NW + warp] = w; }
  __syncthreads();
  for (int i = 0; i < nv; i++) { real s = sh[i * NW]; for (int w = 1; w < NW; w++) s = s + sh[i * NW + w]; v[i] = s; }
  __syncthreads();
}
// render-gpu's grid barrier: the last arrival resets the counter and bumps the generation (fenced)
__device__ __forceinline__ void gsync(unsigned int* bar, unsigned int nb) {
  __syncthreads();
  if (threadIdx.x == 0) {
    volatile unsigned int* gen = bar + 32;
    unsigned int g = *gen;
    __threadfence();
    if (atomicAdd(bar, 1u) == nb - 1) { atomicExch(bar, 0u); __threadfence(); atomicAdd(bar + 32, 1u); }
    else {
      while (*gen == g) {
#if __CUDA_ARCH__ >= 700
        __nanosleep(32);
#endif
      }
    }
    __threadfence();
  }
  __syncthreads();
}
#define GRID_SYNC() gsync(bar, gridDim.x)
__device__ __forceinline__ void store_out(real* o, int i, real v) { o[i] = v; }
#define KERNEL_SIGNATURE extern "C" __global__ void __launch_bounds__(NW * 32) piano(const int* P, const real* C, \
  const real* D, const int* I, const int* ACT, real* S, real* TRACE, real* POS, real* LOSS, unsigned int* bar)
#define KERNEL_PROLOGUE extern __shared__ real tmem[]; const int tid = threadIdx.x; const int blk = blockIdx.x;
"#;

fn df64() -> bool {
    std::env::var("FS_PIANO_GPU_REAL").as_deref() == Ok("df64")
}

fn source(df64: bool) -> String {
    format!(
        "#define NW {}\n{}{PRELUDE_COMMON}{}{PRELUDE_TAIL}{}",
        pack::warps(),
        pack::defines(),
        if df64 { PRELUDE_DF64 } else { PRELUDE_DOUBLE },
        include_str!("kernel.h")
    )
}

/// f64 values as 8-byte device words: the double's bits, or a (hi, lo) float pair.
fn encode(x: &[f64], df64: bool) -> Vec<u64> {
    if !df64 {
        return x.iter().map(|v| v.to_bits()).collect();
    }
    x.iter()
        .map(|&v| {
            let hi = v as f32;
            let lo = (v - f64::from(hi)) as f32;
            u64::from(hi.to_bits()) | (u64::from(lo.to_bits()) << 32)
        })
        .collect()
}
fn decode(x: &[u64], df64: bool, out: &mut [f64]) {
    for (o, &w) in out.iter_mut().zip(x) {
        *o = if df64 {
            f64::from(f32::from_bits(w as u32)) + f64::from(f32::from_bits((w >> 32) as u32))
        } else {
            f64::from_bits(w)
        };
    }
}

fn e<T: std::fmt::Debug>(x: T) -> String {
    format!("{x:?}")
}

/// The compiled kernel for this source and architecture, from the on-disk cache when present.
fn compile(src: &str, arch: &str) -> Result<Ptx, String> {
    use std::hash::{Hash, Hasher};
    let opts = [
        format!("--gpu-architecture={arch}"),
        "-std=c++17".into(),
        "--fmad=true".into(),
    ];
    let mut h = std::collections::hash_map::DefaultHasher::new();
    (src, &opts).hash(&mut h);
    let dir = std::env::var_os("XDG_CACHE_HOME")
        .map(std::path::PathBuf::from)
        .or_else(|| std::env::var_os("HOME").map(|h| std::path::PathBuf::from(h).join(".cache")))
        .map(|d| d.join("fs-render-gpu-piano"));
    let file = dir
        .as_ref()
        .map(|d| d.join(format!("{:016x}-{arch}.cubin", h.finish())));
    if let Some(Ok(bin)) = file.as_ref().map(std::fs::read) {
        return Ok(Ptx::from_binary(bin));
    }
    let bin = nvrtc_cubin(src, &opts)?;
    if let (Some(d), Some(f)) = (dir, file) {
        // best effort and atomic (written aside, then renamed)
        let tmp = d.join(format!(
            ".{}.{}",
            f.file_name()
                .map_or(String::new(), |n| n.to_string_lossy().into_owned()),
            std::process::id()
        ));
        let _ = std::fs::create_dir_all(&d)
            .and_then(|_| std::fs::write(&tmp, &bin))
            .and_then(|_| std::fs::rename(&tmp, &f));
    }
    Ok(Ptx::from_binary(bin))
}

fn nvrtc_cubin(src: &str, opts: &[String]) -> Result<Vec<u8>, String> {
    let c = std::ffi::CString::new(src).map_err(|e| e.to_string())?;
    let prog = result::create_program(&c, None).map_err(|e| format!("nvrtc: {e:?}"))?;
    if let Err(err) = unsafe { result::compile_program(prog, opts) } {
        let log = unsafe { result::get_program_log(prog) }.map_err(e)?;
        let log = unsafe { std::ffi::CStr::from_ptr(log.as_ptr()) }
            .to_string_lossy()
            .into_owned();
        return Err(format!("nvrtc {err:?}:\n{log}"));
    }
    let mut size = 0usize;
    unsafe { nsys::nvrtcGetCUBINSize(prog, &mut size) }
        .result()
        .map_err(e)?;
    let mut buf = vec![0u8; size];
    unsafe { nsys::nvrtcGetCUBIN(prog, buf.as_mut_ptr().cast()) }
        .result()
        .map_err(e)?;
    unsafe { result::destroy_program(prog) }.map_err(e)?;
    Ok(buf)
}

/// Compile both kernel variants for `arch` into the cache without a GPU (NVRTC needs none): on a busy
/// build host a compile can take many minutes, so compile elsewhere and copy the cache directory.
pub fn precompile(arch: &str) -> Result<(), String> {
    for df in [false, true] {
        let t = std::time::Instant::now();
        compile(&source(df), arch)?;
        eprintln!(
            "{arch} {}: {:.1} s",
            if df { "df64" } else { "double" },
            t.elapsed().as_secs_f64()
        );
    }
    Ok(())
}

pub struct Cuda {
    stream: Arc<CudaStream>,
    func: CudaFunction,
    packed: Packed,
    df64: bool,
    name: String,
    params: CudaSlice<i32>,
    consts: CudaSlice<u64>,
    damper: CudaSlice<u64>,
    ints: CudaSlice<i32>,
    act: CudaSlice<i32>,
    state: CudaSlice<u64>,
    trace: CudaSlice<u64>,
    pos: CudaSlice<u64>,
    loss: CudaSlice<u64>,
    bar: CudaSlice<u32>,
    blocks: usize,
    shared: usize,
}

impl Cuda {
    pub fn new(p: &Program) -> Result<Self, String> {
        let device: usize = std::env::var("FS_PIANO_GPU_DEVICE")
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or(0);
        let ctx = std::panic::catch_unwind(|| CudaContext::new(device))
            .map_err(|_| "no CUDA driver".to_string())?
            .map_err(e)?;
        let stream = ctx.default_stream();
        let major = ctx
            .attribute(sys::CUdevice_attribute::CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MAJOR)
            .map_err(e)?;
        let minor = ctx
            .attribute(sys::CUdevice_attribute::CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MINOR)
            .map_err(e)?;
        let sms = ctx
            .attribute(sys::CUdevice_attribute::CU_DEVICE_ATTRIBUTE_MULTIPROCESSOR_COUNT)
            .map_err(e)? as usize;
        let df64 = df64();
        let module = ctx
            .load_module(compile(&source(df64), &format!("sm_{major}{minor}"))?)
            .map_err(e)?;
        let func = module.load_function("piano").map_err(e)?;
        // a cooperative launch needs every block resident: one per multiprocessor is always possible
        let budget = std::env::var("FS_PIANO_GPU_BLOCKS")
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or(sms);
        let packed = pack::pack(p, budget)
            .map_err(|err| format!("the GPU kernel cannot step this piano: {err}"))?;
        let blocks = packed.group_blocks + packed.board_blocks;
        let shared = pack::threadgroup_reals(&packed) * 8;
        if shared > 48 * 1024 {
            return Err(format!("a block needs {shared} B of shared memory"));
        }
        let up = |x: &[f64]| stream.clone_htod(&encode(x, df64)).map_err(e);
        let params = stream.clone_htod(&packed.params).map_err(e)?;
        let consts = up(&packed.consts)?;
        let damper = up(&packed.damper)?;
        let ints = stream.clone_htod(&packed.ints).map_err(e)?;
        let act = stream.clone_htod(&pack::active_flags(p)).map_err(e)?;
        let state = stream.alloc_zeros::<u64>(packed.state_len).map_err(e)?;
        let trace = stream.alloc_zeros::<u64>(SEGMENT * p.r).map_err(e)?;
        let pos = stream
            .alloc_zeros::<u64>(SEGMENT * p.contacts.len().max(1))
            .map_err(e)?;
        let loss = stream.alloc_zeros::<u64>(2 * blocks).map_err(e)?;
        let bar = stream.alloc_zeros::<u32>(64).map_err(e)?;
        let name = format!(
            "{} (CUDA, {}, {blocks} blocks x {})",
            ctx.name().map_err(e)?,
            if df64 { "double-float" } else { "double" },
            pack::threads()
        );
        Ok(Self {
            stream,
            func,
            packed,
            df64,
            name,
            params,
            consts,
            damper,
            ints,
            act,
            state,
            trace,
            pos,
            loss,
            bar,
            blocks,
            shared,
        })
    }
}

impl Backend for Cuda {
    fn name(&self) -> String {
        self.name.clone()
    }

    fn load_dampers(&mut self, p: &Program) -> Result<(), String> {
        let old = self.packed.damper.len();
        pack::pack_dampers(p, &mut self.packed);
        if self.packed.damper.len() != old {
            return Err("damper data changed size".into());
        }
        self.stream
            .memcpy_htod(&encode(&self.packed.damper, self.df64), &mut self.damper)
            .map_err(e)?;
        self.stream
            .memcpy_htod(&pack::active_flags(p), &mut self.act)
            .map_err(e)?;
        Ok(())
    }

    fn run(
        &mut self,
        p: &Program,
        q: &mut [f64],
        v: &mut [f64],
        steps: usize,
        trace: &mut [f64],
        positions: &mut [f64],
    ) -> Result<(f64, f64), String> {
        if steps == 0 || steps > SEGMENT {
            return Err("stretch length out of range".into());
        }
        let (n, r, nc) = (p.n, p.r, p.contacts.len());
        self.packed.params[pack::slot("STEPS")] = steps as i32;
        self.packed.params[pack::slot("SKIP")] = std::env::var("FS_PIANO_GPU_SKIP")
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or(0);
        let at = |s: &str| self.packed.params[pack::slot(s)] as usize;
        let (sq, sv) = (at("SQ"), at("SV"));
        self.stream
            .memcpy_htod(&self.packed.params, &mut self.params)
            .map_err(e)?;
        self.stream
            .memcpy_htod(
                &encode(q, self.df64),
                &mut self.state.slice_mut(sq..sq + n + r),
            )
            .map_err(e)?;
        self.stream
            .memcpy_htod(
                &encode(v, self.df64),
                &mut self.state.slice_mut(sv..sv + n + r),
            )
            .map_err(e)?;
        self.stream.memset_zeros(&mut self.bar).map_err(e)?;
        let cfg = LaunchConfig {
            grid_dim: (self.blocks as u32, 1, 1),
            block_dim: (pack::threads() as u32, 1, 1),
            shared_mem_bytes: self.shared as u32,
        };
        let mut l = self.stream.launch_builder(&self.func);
        l.arg(&self.params)
            .arg(&self.consts)
            .arg(&self.damper)
            .arg(&self.ints)
            .arg(&self.act)
            .arg(&mut self.state)
            .arg(&mut self.trace)
            .arg(&mut self.pos)
            .arg(&mut self.loss)
            .arg(&mut self.bar);
        unsafe { l.launch_cooperative(cfg) }
            .map_err(|err| format!("cooperative launch of {} blocks: {err:?}", self.blocks))?;
        let mut words = vec![0u64; n + r];
        self.stream
            .memcpy_dtoh(&self.state.slice(sq..sq + n + r), &mut words)
            .map_err(e)?;
        decode(&words, self.df64, q);
        self.stream
            .memcpy_dtoh(&self.state.slice(sv..sv + n + r), &mut words)
            .map_err(e)?;
        decode(&words, self.df64, v);
        let mut t = vec![0u64; steps * r];
        self.stream
            .memcpy_dtoh(&self.trace.slice(0..steps * r), &mut t)
            .map_err(e)?;
        decode(&t, self.df64, &mut trace[..steps * r]);
        let mut c = vec![0u64; steps * nc];
        if nc > 0 {
            self.stream
                .memcpy_dtoh(&self.pos.slice(0..steps * nc), &mut c)
                .map_err(e)?;
            decode(&c, self.df64, &mut positions[..steps * nc]);
        }
        let mut lw = vec![0u64; 2 * self.blocks];
        self.stream.memcpy_dtoh(&self.loss, &mut lw).map_err(e)?;
        self.stream.synchronize().map_err(e)?;
        let mut losses = vec![0.0; 2 * self.blocks];
        decode(&lw, self.df64, &mut losses);
        Ok((
            losses.iter().step_by(2).sum(),
            losses.iter().skip(1).step_by(2).sum(),
        ))
    }
}
