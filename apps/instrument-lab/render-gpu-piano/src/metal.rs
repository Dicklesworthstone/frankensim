//! The kernel on Apple silicon (Metal). Apple GPUs have no FP64, so `real` is a double-float pair
//! (hi + lo, about 48 significant bits, error-free transformations on fused multiply-adds; fast math
//! off so the compiler keeps them). All threadgroups stay resident and meet at an atomic grid barrier;
//! the state buffer is `coherent(device)`, without which some threadgroups read stale values across the
//! barrier. Buffers are in shared (unified) memory, so uploads and read-backs are plain copies.
use crate::accel::{Backend, SEGMENT};
use crate::pack::{self, Packed};
use crate::plan::Program;
use metal::*;

const PRELUDE: &str = r#"
#include <metal_stdlib>
using namespace metal;
struct df { float hi; float lo; };
static inline df mk(float h, float l) { df x; x.hi = h; x.lo = l; return x; }
static inline df quick_two_sum(float a, float b) { float s = a + b; return mk(s, b - (s - a)); }
static inline df two_sum(float a, float b) { float s = a + b; float v = s - a; return mk(s, (a - (s - v)) + (b - v)); }
#ifdef SLOPPY_ADD
static inline df operator+(df a, df b) {
  df s = two_sum(a.hi, b.hi);
  return quick_two_sum(s.hi, s.lo + (a.lo + b.lo));
}
#else
static inline df operator+(df a, df b) {
  df s = two_sum(a.hi, b.hi); df t = two_sum(a.lo, b.lo);
  s = quick_two_sum(s.hi, s.lo + t.hi);
  return quick_two_sum(s.hi, s.lo + t.lo);
}
#endif
static inline df operator-(df a) { return mk(-a.hi, -a.lo); }
static inline df operator-(df a, df b) { return a + (-b); }
static inline df operator*(df a, df b) {
  float p = a.hi * b.hi;
  float e = fma(a.hi, b.hi, -p);
  e = fma(a.hi, b.lo, e);
  e = fma(a.lo, b.hi, e);
  return quick_two_sum(p, e);
}
typedef df real;
#define R(x) mk(x##f, 0.0f)
#define PC constant
#define CG const device
#define CI const device
#define SG coherent(device) device
#define OG device
#define TH thread
#define TGM threadgroup
static inline real shfl(real x, int i) { return mk(simd_shuffle(x.hi, (ushort)i), simd_shuffle(x.lo, (ushort)i)); }
static inline real warp_sum(real x) {
  for (ushort m = 16; m >= 1; m >>= 1) x = x + mk(simd_shuffle_xor(x.hi, m), simd_shuffle_xor(x.lo, m));
  return x;
}
#define UNROLL _Pragma("unroll 4")
#define UNROLL_ALL _Pragma("unroll")
#define BLOCK_SYNC() threadgroup_barrier(mem_flags::mem_device | mem_flags::mem_threadgroup)
#define BLOCK_SYNC_TG() threadgroup_barrier(mem_flags::mem_threadgroup)
static inline void block_sum(TH real* v, int nv, TGM real* sh, int tid) {
  const int lane = tid & 31, warp = tid >> 5;
  for (int i = 0; i < nv; i++) { real w = warp_sum(v[i]); if (lane == 0) sh[i * NW + warp] = w; }
  BLOCK_SYNC_TG();
  for (int i = 0; i < nv; i++) { real s = sh[i * NW]; for (int w = 1; w < NW; w++) s = s + sh[i * NW + w]; v[i] = s; }
  BLOCK_SYNC_TG();
}
// Two-level grid barrier: threadgroups arrive in groups of 8 (one counter per 256-byte span each); the
// last of each group arrives at the root; the last at the root bumps the generation everyone waits on.
// Counters only grow (the host zeroes them per launch) and a round's last arrival is recognised modulo
// the group size, so no counter is reset while another threadgroup may already use it again.
// Metal does not promise that all threadgroups of a dispatch run concurrently, and very rarely one is
// not scheduled while the others spin (about once per few hundred launches here). So a wait that
// lasts far longer than any phase raises the abort word, every threadgroup returns, and the host
// reruns the stretch from its saved state.
#define BAR_ABORT (33 * 64)
#define SPIN_LIMIT (1u << 19)
static inline bool grid_sync(device atomic_uint* bar, uint nblk, uint tg, int tid, threadgroup int* flag) {
  threadgroup_barrier(mem_flags::mem_device | mem_flags::mem_threadgroup);
  if (tid == 0) {
    atomic_thread_fence(mem_flags::mem_device, memory_order_seq_cst, thread_scope_device);
    device atomic_uint* gw = &bar[32 * 64];
    const uint gen = atomic_load_explicit(gw, memory_order_relaxed);
    const uint c = tg / 8, groups = (nblk + 7) / 8, members = min(8u, nblk - c * 8);
    bool released = false;
    int ok = 1;
    if (atomic_fetch_add_explicit(&bar[(1 + c) * 64], 1u, memory_order_relaxed) % members == members - 1) {
      if (atomic_fetch_add_explicit(&bar[0], 1u, memory_order_relaxed) % groups == groups - 1) {
        atomic_fetch_add_explicit(gw, 1u, memory_order_relaxed);
        released = true;
      }
    }
    if (!released) {
      uint spins = 0;
      while (atomic_load_explicit(gw, memory_order_relaxed) == gen) {
        spins++;
        if ((spins & 1023u) == 0u && atomic_load_explicit(&bar[BAR_ABORT], memory_order_relaxed) != 0u) { ok = 0; break; }
        if (spins > SPIN_LIMIT) { atomic_store_explicit(&bar[BAR_ABORT], 1u, memory_order_relaxed); ok = 0; break; }
      }
    }
    atomic_thread_fence(mem_flags::mem_device, memory_order_seq_cst, thread_scope_device);
    *flag = ok;
  }
  threadgroup_barrier(mem_flags::mem_device | mem_flags::mem_threadgroup);
  return *flag != 0;
}
#define GRID_SYNC() do { if (!grid_sync(bar, nblk_, blk_, tid, &gflag)) return; } while (0)
static inline void store_out(OG real* o, int i, real v) { o[i] = v; }
#define KERNEL_SIGNATURE kernel void piano(PC int* P [[buffer(0)]], CG real* C [[buffer(1)]], \
  CG real* D [[buffer(2)]], CI int* I [[buffer(3)]], CI int* ACT [[buffer(4)]], SG real* S [[buffer(5)]], \
  OG real* TRACE [[buffer(6)]], OG real* POS [[buffer(7)]], OG real* LOSS [[buffer(8)]], \
  device atomic_uint* bar [[buffer(9)]], threadgroup real* tmem [[threadgroup(0)]], \
  uint tid_ [[thread_position_in_threadgroup]], uint blk_ [[threadgroup_position_in_grid]], \
  uint nblk_ [[threadgroups_per_grid]])
#define KERNEL_PROLOGUE threadgroup int gflag; const int tid = (int)tid_; const int blk = (int)blk_;
"#;

pub const KERNEL: &str = include_str!("kernel.h");
/// The grid barrier's words: the root counter, one counter per group of 8 threadgroups, the generation,
/// the abort word (BAR_ABORT).
const BAR_WORDS: usize = 34 * 64;
const BAR_ABORT: usize = 33 * 64;
/// Launches of one stretch before giving up (each abort costs well under a second).
const ATTEMPTS: usize = 4;

/// How many threadgroups of the kernel can be resident at once: one per GPU core (more relied on
/// time-slicing and hung a grid barrier on long launches). FS_PIANO_GPU_BLOCKS overrides.
fn co_resident_blocks() -> Result<usize, String> {
    if let Some(b) = std::env::var("FS_PIANO_GPU_BLOCKS")
        .ok()
        .and_then(|s| s.parse().ok())
    {
        return Ok(b);
    }
    let out = std::process::Command::new("/usr/sbin/ioreg")
        .args(["-rc", "AGXAccelerator"])
        .output()
        .map_err(|e| format!("cannot read the GPU core count: {e}"))?;
    let text = String::from_utf8_lossy(&out.stdout);
    text.lines()
        .find_map(|l| {
            l.split_once("\"gpu-core-count\" = ")
                .and_then(|(_, v)| v.trim().parse().ok())
        })
        .filter(|&c: &usize| c >= 8 && c <= 256)
        .ok_or_else(|| {
            "cannot read the GPU core count (ioreg gpu-core-count); set FS_PIANO_GPU_BLOCKS".into()
        })
}

pub struct Metal {
    queue: CommandQueue,
    pipeline: ComputePipelineState,
    packed: Packed,
    name: String,
    params: Buffer,
    consts: Buffer,
    damper: Buffer,
    ints: Buffer,
    act: Buffer,
    state: Buffer,
    trace: Buffer,
    pos: Buffer,
    loss: Buffer,
    bar: Buffer,
    blocks: usize,
    tgmem: usize,
    aborts: usize,
    device: Device,
}

fn buffer<T: Copy>(device: &Device, data: &[T]) -> Buffer {
    let bytes = std::mem::size_of_val(data).max(16);
    let b = device.new_buffer(bytes as u64, MTLResourceOptions::StorageModeShared);
    unsafe {
        std::ptr::copy_nonoverlapping(
            data.as_ptr().cast::<u8>(),
            b.contents().cast::<u8>(),
            std::mem::size_of_val(data),
        );
    }
    b
}
fn write<T: Copy>(b: &Buffer, offset: usize, data: &[T]) {
    assert!(std::mem::size_of_val(data) + offset * std::mem::size_of::<T>() <= b.length() as usize);
    unsafe {
        std::ptr::copy_nonoverlapping(
            data.as_ptr(),
            b.contents().cast::<T>().add(offset),
            data.len(),
        );
    }
}
fn read<T: Copy>(b: &Buffer, offset: usize, len: usize) -> &[T] {
    assert!((offset + len) * std::mem::size_of::<T>() <= b.length() as usize);
    unsafe { std::slice::from_raw_parts(b.contents().cast::<T>().add(offset), len) }
}

impl Metal {
    pub fn new(p: &Program) -> Result<Self, String> {
        let device = Device::system_default().ok_or("no Metal device")?;
        // the double-float add: "sloppy" (default; error relative to the operands, like rounding the
        // inputs) or "ieee" (error relative to the result; about 10% slower here)
        let sloppy = if std::env::var("FS_PIANO_GPU_ADD").as_deref() == Ok("ieee") {
            ""
        } else {
            "#define SLOPPY_ADD\n"
        };
        let source = format!(
            "#define NW {}\n{sloppy}{}{PRELUDE}{KERNEL}",
            pack::warps(),
            pack::defines()
        );
        let options = CompileOptions::new();
        options.set_fast_math_enabled(false);
        let library = device
            .new_library_with_source(&source, &options)
            .map_err(|e| format!("Metal kernel compile: {e}"))?;
        let function = library
            .get_function("piano", None)
            .map_err(|e| format!("Metal kernel: {e}"))?;
        let pipeline = device
            .new_compute_pipeline_state_with_function(&function)
            .map_err(|e| format!("Metal pipeline: {e}"))?;
        if (pipeline.max_total_threads_per_threadgroup() as usize) < pack::threads() {
            return Err(format!(
                "the Metal kernel allows {} threads per threadgroup, needs {}",
                pipeline.max_total_threads_per_threadgroup(),
                pack::threads()
            ));
        }
        let budget = co_resident_blocks()?;
        let packed = pack::pack(p, budget)
            .map_err(|e| format!("the GPU kernel cannot step this piano: {e}"))?;
        let blocks = packed.group_blocks + packed.board_blocks;
        let df = |x: &[f64]| pack::to_df64(x);
        let params = buffer(&device, &packed.params);
        let consts = buffer(&device, &df(&packed.consts));
        let damper = buffer(&device, &df(&packed.damper));
        let ints = buffer(&device, &packed.ints);
        let act = buffer(&device, &pack::active_flags(p));
        let state = buffer(&device, &vec![[0.0f32; 2]; packed.state_len]);
        let trace = buffer(&device, &vec![[0.0f32; 2]; SEGMENT * p.r]);
        let pos = buffer(
            &device,
            &vec![[0.0f32; 2]; SEGMENT * p.contacts.len().max(1)],
        );
        let loss = buffer(&device, &vec![[0.0f32; 2]; 2 * blocks]);
        let bar = buffer(&device, &[0u32; BAR_WORDS]);
        let tgmem = pack::threadgroup_reals(&packed) * 8;
        if tgmem > device.max_threadgroup_memory_length() as usize {
            return Err(format!(
                "a bridge group needs {tgmem} B of threadgroup memory; this GPU has {}",
                device.max_threadgroup_memory_length()
            ));
        }
        let name = format!(
            "{} (Metal, double-float, {blocks} threadgroups x {})",
            device.name(),
            pack::threads()
        );
        Ok(Self {
            queue: device.new_command_queue(),
            pipeline,
            packed,
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
            tgmem,
            aborts: 0,
            device,
        })
    }
}

impl Backend for Metal {
    fn name(&self) -> String {
        self.name.clone()
    }
    fn stalls(&self) -> usize {
        self.aborts
    }

    fn load_dampers(&mut self, p: &Program) -> Result<(), String> {
        let old = self.packed.damper.len();
        pack::pack_dampers(p, &mut self.packed);
        if self.packed.damper.len() != old {
            return Err("damper data changed size".into());
        }
        write(&self.damper, 0, &pack::to_df64(&self.packed.damper));
        write(&self.act, 0, &pack::active_flags(p));
        write(&self.params, 0, &self.packed.params);
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
        // a timing diagnostic only (kernel.h P_SKIP)
        self.packed.params[pack::slot("SKIP")] = std::env::var("FS_PIANO_GPU_SKIP")
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or(0);
        let params = self.packed.params.clone();
        let slot = |s: &str| params[pack::slot(s)] as usize;
        write(&self.params, 0, &self.packed.params);
        write(&self.state, slot("SQ"), &pack::to_df64(q));
        write(&self.state, slot("SV"), &pack::to_df64(v));
        let mut attempt = 0;
        loop {
            attempt += 1;
            write(&self.bar, 0, &[0u32; BAR_WORDS]);
            let cb = self.queue.new_command_buffer();
            let enc = cb.new_compute_command_encoder();
            enc.set_compute_pipeline_state(&self.pipeline);
            for (i, b) in [
                &self.params,
                &self.consts,
                &self.damper,
                &self.ints,
                &self.act,
                &self.state,
                &self.trace,
                &self.pos,
                &self.loss,
                &self.bar,
            ]
            .into_iter()
            .enumerate()
            {
                enc.set_buffer(i as u64, Some(b), 0);
            }
            enc.set_threadgroup_memory_length(0, self.tgmem as u64);
            enc.dispatch_thread_groups(
                MTLSize::new(self.blocks as u64, 1, 1),
                MTLSize::new(pack::threads() as u64, 1, 1),
            );
            enc.end_encoding();
            cb.commit();
            cb.wait_until_completed();
            if cb.status() != MTLCommandBufferStatus::Completed {
                return Err(format!("Metal command buffer ended {:?}", cb.status()));
            }
            if read::<u32>(&self.bar, BAR_ABORT, 1)[0] == 0 {
                break;
            }
            self.aborts += 1;
            if attempt >= ATTEMPTS {
                return Err(format!(
                    "the GPU stretch stalled {attempt} times (a threadgroup was not scheduled)"
                ));
            }
            // rerun from the same start
            write(&self.state, slot("SQ"), &pack::to_df64(q));
            write(&self.state, slot("SV"), &pack::to_df64(v));
        }
        pack::from_df64(read(&self.state, slot("SQ"), n + r), q);
        pack::from_df64(read(&self.state, slot("SV"), n + r), v);
        pack::from_df64(read(&self.trace, 0, steps * r), &mut trace[..steps * r]);
        pack::from_df64(read(&self.pos, 0, steps * nc), &mut positions[..steps * nc]);
        let mut losses = vec![0.0; 2 * self.blocks];
        pack::from_df64(read(&self.loss, 0, 2 * self.blocks), &mut losses);
        let modal = losses.iter().step_by(2).sum();
        let damper = losses.iter().skip(1).step_by(2).sum();
        let _ = &self.device;
        Ok((modal, damper))
    }
}
