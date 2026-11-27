//! Wasm32 `wasm-bindgen` surface over [`fs_runtime::AssemblyRuntime`].

#![forbid(unsafe_code)]
#![warn(missing_docs)]

use fs_runtime::{AssemblyRuntime, MAX_PACKS, RuntimeRefusal};

#[cfg(target_arch = "wasm32")]
use wasm_bindgen::prelude::*;

/// Browser / Node / AudioWorklet wrapper around [`AssemblyRuntime`].
#[cfg_attr(target_arch = "wasm32", wasm_bindgen)]
pub struct WasmRuntime {
    inner: AssemblyRuntime,
}

fn framing_err(msg: impl Into<String>) -> String {
    RuntimeRefusal::new("invalid-pack-framing", "load", msg).to_json()
}

fn slice_packs<'a>(pack_bytes: &'a [u8], pack_lengths: &[u32]) -> Result<Vec<&'a [u8]>, String> {
    if pack_lengths.len() > MAX_PACKS {
        return Err(framing_err(format!(
            "too many packs: {} exceeds limit {MAX_PACKS}",
            pack_lengths.len()
        )));
    }
    let mut packs = Vec::with_capacity(pack_lengths.len());
    let mut pos = 0usize;
    for (i, &len_u32) in pack_lengths.iter().enumerate() {
        let len = usize::try_from(len_u32)
            .map_err(|_| framing_err(format!("pack {i} length exceeds address space")))?;
        let end = pos
            .checked_add(len)
            .ok_or_else(|| framing_err("pack length overflow"))?;
        let slice = pack_bytes.get(pos..end).ok_or_else(|| {
            framing_err(format!(
                "pack {i} range {pos}..{end} exceeds {} bytes",
                pack_bytes.len()
            ))
        })?;
        packs.push(slice);
        pos = end;
    }
    if pos != pack_bytes.len() {
        return Err(framing_err(format!(
            "trailing pack bytes: consumed {pos} of {}",
            pack_bytes.len()
        )));
    }
    Ok(packs)
}

#[cfg_attr(target_arch = "wasm32", wasm_bindgen)]
impl WasmRuntime {
    /// Load a `.fsim` assembly and concatenated `FSMCDPK\0`/`FSINTPK\0` packs.
    #[cfg_attr(target_arch = "wasm32", wasm_bindgen(js_name = load))]
    pub fn load(
        fsim_bytes: &[u8],
        pack_bytes: &[u8],
        pack_lengths: &[u32],
        seed_override: Option<u64>,
    ) -> Result<Self, String> {
        let packs = slice_packs(pack_bytes, pack_lengths)?;
        let inner =
            AssemblyRuntime::load(fsim_bytes, &packs, seed_override).map_err(|r| r.to_json())?;
        Ok(Self { inner })
    }

    /// Canonical JSON array of interactive input descriptors.
    #[must_use]
    pub fn inputs(&self) -> String {
        self.inner.inputs_json()
    }

    /// Set one interactive input by index.
    pub fn set_input(&self, index: usize, value: f64) -> Result<(), String> {
        self.inner.set_input(index, value).map_err(|r| r.to_json())
    }

    /// Set one interactive input by declared name.
    pub fn set_input_by_name(&self, name: &str, value: f64) -> Result<(), String> {
        self.inner
            .set_input_by_name(name, value)
            .map_err(|r| r.to_json())
    }

    /// Canonical JSON snapshot descriptor (`SnapshotDescriptor`).
    #[must_use]
    pub fn probes(&self) -> String {
        self.inner.probes_json()
    }

    /// Required `Float64Array` length for [`Self::read_probes`].
    #[must_use]
    pub fn snapshot_len(&self) -> usize {
        self.inner.snapshot_len()
    }

    /// Advance by `dt_s` seconds and return total accepted steps.
    pub fn step(&mut self, dt_s: f64) -> Result<u64, String> {
        self.inner
            .step(dt_s)
            .map(|_| self.inner.steps())
            .map_err(|r| r.to_json())
    }

    /// Fill `out[..snapshot_len()]` with current flat probe values.
    pub fn read_probes(&self, out: &mut [f64]) -> Result<usize, String> {
        self.inner.read_probes(out).map_err(|r| r.to_json())
    }

    /// AudioWorklet block render into `out` (`Float32Array`), halting with silence on refusal.
    pub fn render_f32(
        &mut self,
        probe_index: usize,
        dt_s: f64,
        substeps: usize,
        out: &mut [f32],
    ) -> usize {
        self.inner.render_f32(probe_index, dt_s, substeps, out)
    }

    /// Canonical JSON export of the kernel window ledger.
    #[must_use]
    pub fn ledger(&self) -> String {
        self.inner.ledger_json()
    }

    /// Poll and clear the refusal that halted the runtime, reported at most once.
    pub fn take_refusal(&mut self) -> Option<String> {
        self.inner.take_refusal().map(|r| r.to_json())
    }

    /// Whether the runtime has halted after a refusal.
    #[must_use]
    pub fn is_halted(&self) -> bool {
        self.inner.is_halted()
    }

    /// Current accepted simulation time [s].
    #[must_use]
    pub fn time_s(&self) -> f64 {
        self.inner.time_s()
    }

    /// Total accepted kernel steps.
    #[must_use]
    pub fn steps(&self) -> u64 {
        self.inner.steps()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn forced_oscillator_native_matches_wasm_parity_bits() {
        let fsim = include_bytes!("../../../examples/assemblies/forced_oscillator.fsim");
        let steel: &[u8] = include_bytes!(
            "../../../examples/assemblies/packs/steel-music-wire-heptonstall-ligo-p1400229.fsmcdpk"
        );
        let len = u32::try_from(steel.len()).unwrap();
        let mut rt = WasmRuntime::load(fsim, steel, &[len], None).unwrap();
        assert_eq!(rt.snapshot_len(), 2);
        assert!(rt.inputs().contains("\"name\":\"load\""));
        assert!(rt.probes().contains("\"name\":\"mass_pos\""));
        for _ in 0..50 {
            rt.step(1.0e-3).unwrap();
        }
        let mut snap = [0.0_f64; 2];
        rt.read_probes(&mut snap).unwrap();
        rt.set_input(0, -5.0).unwrap();
        for _ in 0..25 {
            rt.step(1.0e-3).unwrap();
        }
        let mut snap2 = [0.0_f64; 2];
        rt.read_probes(&mut snap2).unwrap();
        let mut block = [0.0_f32; 8];
        assert_eq!(rt.render_f32(0, 1.0e-3, 1, &mut block), 8);
        assert_eq!(snap[0].to_bits(), 0x3f78_d8c3_5603_6c4c);
        assert_eq!(snap[1].to_bits(), 0x3fcd_e0a2_9744_0181);
        assert_eq!(snap2[0].to_bits(), 0x3f85_f471_fa9d_6701);
        assert_eq!(snap2[1].to_bits(), 0x3fc1_a15b_df0f_7a33);
        assert_eq!(block[7].to_bits(), 0x3c3f_9dbb);
        let bad = WasmRuntime::load(fsim, steel, &[len - 1], None);
        assert!(bad.err().unwrap().contains("invalid-pack-framing"));
        let too_many = vec![0u32; MAX_PACKS + 1];
        let over = WasmRuntime::load(fsim, &[], &too_many, None);
        assert!(over.err().unwrap().contains("too many packs"));
    }
}
