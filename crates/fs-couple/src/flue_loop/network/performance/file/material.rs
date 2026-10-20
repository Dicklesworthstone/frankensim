//! One bounded immutable material-source cache for walls and moving boundaries.
use super::{FlueInputError, bad, poll};
use fs_exec::CancelGate;
use fs_matdb::NormalizedMaterialCardPack;
use std::{collections::BTreeMap, sync::Arc};

pub(super) const MAX_PACK_BYTES: usize = 8 * 1024 * 1024;
const MAX_TOTAL_BYTES: usize = 64 * 1024 * 1024;

pub(super) fn relative_path(path: &str, line: usize) -> Result<(), FlueInputError> {
    if path.len() > 512
        || path.contains('\\')
        || path.contains(':')
        || path
            .split('/')
            .any(|part| part.is_empty() || part == "." || part == "..")
    {
        return Err(bad(
            line,
            "material pack path must remain relative to the performance directory",
        ));
    }
    Ok(())
}

#[derive(Default)]
pub(super) struct Sources {
    packs: BTreeMap<String, Arc<NormalizedMaterialCardPack>>,
    bytes: usize,
}
impl Sources {
    pub fn load(
        &mut self,
        path: &str,
        line: usize,
        loader: &mut impl FnMut(&str) -> Result<NormalizedMaterialCardPack, String>,
        gate: &CancelGate,
    ) -> Result<Arc<NormalizedMaterialCardPack>, FlueInputError> {
        poll(gate)?;
        relative_path(path, line)?;
        if !self.packs.contains_key(path) {
            let pack =
                loader(path).map_err(|source| FlueInputError::MaterialSource { line, source })?;
            poll(gate)?;
            let bytes = pack.to_bytes().len();
            if bytes > MAX_PACK_BYTES || bytes > MAX_TOTAL_BYTES.saturating_sub(self.bytes) {
                return Err(bad(
                    line,
                    "material sources exceed 8 MiB per pack or 64 MiB combined",
                ));
            }
            self.bytes += bytes;
            self.packs.insert(path.into(), Arc::new(pack));
        }
        Ok(Arc::clone(&self.packs[path]))
    }
}
