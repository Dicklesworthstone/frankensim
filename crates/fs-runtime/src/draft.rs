//! Lowering `.fsim` sources and admitted [`DraftAssembly`] graphs into [`AssemblyRuntime`].

use fs_fsim::{PackCatalog, admit_assembly_source};
use fs_scenario::assembly::draft::DraftAssembly;
use fs_system::{CardStore, LoweringRegistry};

use crate::descriptor::RuntimeRefusal;
use crate::probes::{from_lowered, load_err};
use crate::runtime::AssemblyRuntime;

impl AssemblyRuntime {
    /// Parse and admit a `.fsim` source with canonical `fs-matdb` card packs
    /// (`FSMCDPK\0` / `FSINTPK\0`), then lower into a runtime.
    pub fn load(
        fsim_source: &[u8],
        packs: &[&[u8]],
        seed_override: Option<u64>,
    ) -> Result<Self, RuntimeRefusal> {
        let cat = PackCatalog::from_raw_bytes(packs).map_err(|r| load_err(r.code, r.what))?;
        let adm = admit_assembly_source(fsim_source, &cat).map_err(|r| load_err(r.code, r.what))?;
        let seed = seed_override.unwrap_or(adm.seed);
        Self::from_draft(&adm.draft, &adm.card_store, seed)
    }

    /// Lower an admitted [`DraftAssembly`] through the shared kernel registry.
    pub fn from_draft(
        draft: &DraftAssembly,
        cards: &CardStore,
        seed: u64,
    ) -> Result<Self, RuntimeRefusal> {
        let lowered = LoweringRegistry::with_kernel_defaults()
            .lower(draft, cards)
            .map_err(|e| load_err("lowering-refused", e.to_string()))?;
        from_lowered(draft, lowered, Vec::new(), Vec::new(), seed)
    }
}
