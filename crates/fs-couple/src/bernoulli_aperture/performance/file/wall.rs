//! Authored DG wall sources share the existing material loader with the plate.
//! Geometry, source partitions and constitutive admission remain DG-owned.

use super::{DuctInput, PlateValveInputError, Reader, bad, checkpoint};
use crate::dg_port::input::wall::regions::{DgMaterialWallLayout, DgMaterialWallLayoutInput};
use fs_exec::CancelGate;
use fs_matdb::NormalizedMaterialCardPack;
use std::{collections::BTreeMap, sync::Arc};

/// Wall declarations have the shared count/byte caps; the original section
/// count bounds additional plate artifacts. Every path is loaded once.
pub(super) type MaterialPacks = BTreeMap<String, Arc<NormalizedMaterialCardPack>>;

pub(super) fn load_pack(
    cache: &mut MaterialPacks,
    path: &str,
    line: usize,
    load_material: &mut dyn FnMut(&str) -> Result<NormalizedMaterialCardPack, String>,
    gate: &CancelGate,
) -> Result<Arc<NormalizedMaterialCardPack>, PlateValveInputError> {
    checkpoint(gate)?;
    if let Some(pack) = cache.get(path) {
        return Ok(Arc::clone(pack));
    }
    let pack = Arc::new(
        load_material(path)
            .map_err(|source| PlateValveInputError::MaterialSource { line, source })?,
    );
    checkpoint(gate)?;
    cache.insert(path.to_owned(), Arc::clone(&pack));
    Ok(pack)
}

pub(super) fn read(
    rows: &mut Reader<'_>,
    duct: &DuctInput,
    cache: &mut MaterialPacks,
    load_material: &mut dyn FnMut(&str) -> Result<NormalizedMaterialCardPack, String>,
    gate: &CancelGate,
) -> Result<Option<DgMaterialWallLayout>, PlateValveInputError> {
    if rows
        .lines
        .clone()
        .next()
        .and_then(|line| line.split_ascii_whitespace().next())
        != Some("dg_wall")
    {
        return Ok(None);
    }
    let mut row = rows.row("dg_wall")?;
    if !matches!(duct, DuctInput::Dg(_)) {
        return Err(bad(
            row.line,
            "dg_wall requires an original dg_tube or dg_bore source",
        ));
    }
    let declaration = DgMaterialWallLayoutInput::read(row.line, &mut row.fields, &mut || {
        let next = rows.next()?;
        Ok::<_, PlateValveInputError>((next.line, next.fields))
    })?;
    row.finish()?;
    let layout = declaration.resolve(gate, &mut |path, line| {
        load_pack(cache, path, line, load_material, gate)
    })?;
    checkpoint(gate)?;
    Ok(Some(layout))
}

#[cfg(test)]
mod tests;
