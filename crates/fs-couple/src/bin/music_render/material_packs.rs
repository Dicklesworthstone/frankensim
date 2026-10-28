//! Shared bounded loading of explicitly named native material-card sources.
use fs_matdb::NormalizedMaterialCardPack;
use std::{io::Read, path::Path};

/// The existing native section input limit, shared by all callers.
const MAX_SECTION_PACK_BYTES: usize = 8 * 1024 * 1024;

/// Load already-admitted material-card bytes relative to their owning input.
/// Callers parse names with `MaterialRecord`; no material or state is inferred.
pub(super) fn section_pack_loader(parent: &Path)
    -> impl FnMut(&str) -> Result<NormalizedMaterialCardPack, String> + '_
{
    move |name| {
        let file = std::fs::File::open(parent.join(name))
            .map_err(|e| format!("material-card pack {name:?} open failed: {e}"))?;
        let mut bytes = Vec::new();
        file.take((MAX_SECTION_PACK_BYTES + 1) as u64).read_to_end(&mut bytes)
            .map_err(|e| format!("material-card pack {name:?} read failed: {e}"))?;
        if bytes.len() > MAX_SECTION_PACK_BYTES {
            return Err("material-card pack exceeds 8 MiB".into());
        }
        NormalizedMaterialCardPack::from_bytes(&bytes)
            .map_err(|e| format!("material-card pack {name:?}: {e}"))
    }
}
