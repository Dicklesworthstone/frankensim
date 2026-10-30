//! Physical axial assignments of the existing material wall source.
//! Regions are coordinates in the original bore, not a list tied to one mesh.
//! Every interface must coincide with an actual DG element face; no averaged
//! source or fractional element is introduced to make a declaration fit.

use super::{DgMaterialWallInput, DgMaterialWallRecord, MAX_DG_WALL_PACK_BYTES};
use super::super::{DgSourceError, Row, bad, checkpoint};
use crate::bernoulli_aperture::wall::material::RegionalMaterialWallDg;
use crate::dg_port::{DgBranch, DgPortWallSource};
use fs_duct::dg::{DgDuct, PortId};
use fs_exec::CancelGate;
use fs_matdb::NormalizedMaterialCardPack;
use std::{collections::BTreeMap, str::SplitAsciiWhitespace, sync::Arc};

/// Bounds source resolution independently of the later complete-state budget.
pub const MAX_DG_WALL_REGIONS: usize = 256;

struct RegionInput {
    interval_m: [f64; 2],
    source: DgMaterialWallInput,
}

enum LayoutInput {
    Rigid,
    Uniform(DgMaterialWallInput),
    Regions(Vec<RegionInput>),
}

/// One uniform source or a complete, physically positioned source partition.
pub struct DgMaterialWallLayoutInput {
    line: usize,
    layout: LayoutInput,
    chimneys: Vec<(PortId, DgMaterialWallLayoutInput)>,
}

/// Original material source and its authored physical support.
#[derive(Clone, Debug)]
pub struct DgMaterialWallRegion {
    interval_m: [f64; 2],
    source: DgMaterialWallRecord,
}

impl DgMaterialWallRegion {
    /// Axial start and end from the actual bore inlet [m].
    pub const fn interval_m(&self) -> [f64; 2] { self.interval_m }
    /// Full source artifact and resolved original constitutive law.
    pub const fn source(&self) -> &DgMaterialWallRecord { &self.source }
}

/// Frozen source assignment, retained independently of propagation resolution.
/// Empty/default means the original rigid wall, with no invented stiff solid.
#[derive(Clone, Debug, Default)]
pub struct DgMaterialWallLayout {
    line: usize,
    uniform: Option<DgMaterialWallRecord>,
    regions: Vec<DgMaterialWallRegion>,
    chimneys: Vec<DgMaterialWallChimney>,
}

/// Explicit wall source for a finite DG branch, addressed by its original port.
#[derive(Clone, Debug)]
pub struct DgMaterialWallChimney {
    port: PortId,
    layout: DgMaterialWallLayout,
}

impl DgMaterialWallChimney {
    /// Original side-port identity, independent of authored or numerical order.
    pub const fn port(&self) -> PortId { self.port }
    /// Uniform or regional wall source, in coordinates from this branch inlet.
    pub const fn layout(&self) -> &DgMaterialWallLayout { &self.layout }
}

impl DgMaterialWallLayoutInput {
    /// Read the existing rigid/isotropic row, or "regions N" followed by N rows:
    /// "dg_wall_region X_BEGIN X_END isotropic PACK T OMEGA H AXIAL FMIN FMAX LOSS".
    /// Regions cover the complete selected duct. An assembly names one main
    /// duct and independently supplied finite side chimneys:
    /// `assembly COUNT`, followed by `dg_wall_main LAYOUT` and
    /// `dg_wall_chimney PORT_ID LAYOUT`. Each LAYOUT uses this same rigid,
    /// isotropic, orthotropic or regional grammar; nested assemblies refuse.
    /// Omitted chimneys retain their existing explicit rigid-wall model.
    /// Partial rigid/material regional mixtures are not silently completed.
    pub fn read<'a, E: From<DgSourceError>>(
        line: usize,
        fields: &mut SplitAsciiWhitespace<'a>,
        next: &mut impl FnMut() -> Result<(usize, SplitAsciiWhitespace<'a>), E>,
    ) -> Result<Self, E> {
        if fields.clone().next() == Some("assembly") {
            let mut row = Row { line, fields };
            row.word()?;
            let count = row.count(MAX_DG_WALL_REGIONS + 1)?;
            row.finish()?;
            if count == 0 {
                return Err(bad(line, "DG wall assembly needs one main duct and explicit branch targets").into());
            }
            let mut main = None;
            let mut chimneys = Vec::new();
            for _ in 0..count {
                let (line, mut fields) = next()?;
                let mut row = Row { line, fields: &mut fields };
                let target = match row.word()? {
                    "dg_wall_main" => {
                        if main.is_some() {
                            return Err(bad(line, "DG wall assembly repeats its main duct").into());
                        }
                        None
                    }
                    "dg_wall_chimney" => {
                        let port = PortId(row.word()?.parse::<u32>()
                            .map_err(|_| bad(line, "DG wall chimney port must be a u32 identity"))?);
                        if chimneys.iter().any(|(previous, _)| *previous == port) {
                            return Err(bad(line, "DG wall assembly repeats a chimney port").into());
                        }
                        Some(port)
                    }
                    _ => return Err(bad(line, "expected dg_wall_main or dg_wall_chimney").into()),
                };
                if row.fields.clone().next() == Some("assembly") {
                    return Err(bad(line, "DG wall targets cannot contain nested assemblies").into());
                }
                let layout = Self::read(line, row.fields, next)?;
                match target {
                    Some(port) => chimneys.push((port, layout)),
                    None => main = Some(layout),
                }
            }
            let mut main = main.ok_or_else(|| bad(line, "DG wall assembly must declare its main duct"))?;
            main.chimneys = chimneys;
            return Ok(main);
        }
        if fields.clone().next() != Some("regions") {
            let layout = DgMaterialWallInput::read(line, fields)?
                .map_or(LayoutInput::Rigid, LayoutInput::Uniform);
            return Ok(Self { line, layout, chimneys: Vec::new() });
        }
        let mut row = Row { line, fields };
        row.word()?;
        let count = row.count(MAX_DG_WALL_REGIONS)?;
        row.finish()?;
        if count == 0 {
            return Err(bad(line, "DG wall regions require a nonempty complete source partition").into());
        }
        let mut regions: Vec<RegionInput> = Vec::with_capacity(count);
        for _ in 0..count {
            let (line, mut fields) = next()?;
            let mut row = Row { line, fields: &mut fields };
            if row.word()? != "dg_wall_region" {
                return Err(bad(line, "expected dg_wall_region with physical axial coordinates").into());
            }
            let interval_m = [row.scalar()?, row.positive()?];
            if interval_m[0] < 0.0 || interval_m[0] >= interval_m[1] {
                return Err(bad(line, "DG wall region needs finite ordered nonnegative coordinates").into());
            }
            // Adjacent intervals must name the same physical boundary. Mesh
            // roundoff is handled at binding, never by filling a source gap.
            let expected_start = regions.last().map_or(0.0, |region| region.interval_m[1]);
            if interval_m[0] != expected_start {
                return Err(bad(line, "DG wall regions must start at zero and be contiguous in authored order").into());
            }
            let source = DgMaterialWallInput::read(line, row.fields)?
                .ok_or_else(|| bad(line, "every DG wall region requires an original isotropic source"))?;
            row.finish()?;
            regions.push(RegionInput { interval_m, source });
        }
        Ok(Self { line, layout: LayoutInput::Regions(regions), chimneys: Vec::new() })
    }

    /// Resolve each original artifact once, sharing it across regions and with
    /// the caller's existing blade/material cache. Distinct local temperatures,
    /// thicknesses, restraints and loss spectra remain separate source queries.
    /// The total unique wall artifacts share the existing 8 MiB source budget.
    pub fn resolve<E: From<DgSourceError>>(
        self,
        gate: &CancelGate,
        load: &mut impl FnMut(&str, usize) -> Result<Arc<NormalizedMaterialCardPack>, E>,
    ) -> Result<DgMaterialWallLayout, E> {
        let mut cache = BTreeMap::<String, Arc<NormalizedMaterialCardPack>>::new();
        let mut bytes = 0usize;
        self.resolve_cached(gate, load, &mut cache, &mut bytes)
    }

    fn resolve_cached<E: From<DgSourceError>>(
        self, gate: &CancelGate,
        load: &mut impl FnMut(&str, usize) -> Result<Arc<NormalizedMaterialCardPack>, E>,
        cache: &mut BTreeMap<String, Arc<NormalizedMaterialCardPack>>, bytes: &mut usize,
    ) -> Result<DgMaterialWallLayout, E> {
        checkpoint(gate, self.line)?;
        let mut resolve = |input: DgMaterialWallInput| -> Result<DgMaterialWallRecord, E> {
            checkpoint(gate, input.line())?;
            let path = input.path().to_owned();
            let pack = if let Some(pack) = cache.get(&path) {
                Arc::clone(pack)
            } else {
                let pack = load(&path, input.line())?;
                *bytes = bytes.checked_add(pack.to_bytes().len())
                    .ok_or_else(|| bad(input.line(), "DG wall source byte count overflow"))?;
                if *bytes > MAX_DG_WALL_PACK_BYTES {
                    return Err(bad(input.line(), "complete DG wall sources exceed 8 MiB").into());
                }
                cache.insert(path, Arc::clone(&pack));
                pack
            };
            input.resolve(pack, gate).map_err(E::from)
        };
        let (uniform, regions) = match self.layout {
            LayoutInput::Rigid => (None, Vec::new()),
            LayoutInput::Uniform(input) => (Some(resolve(input)?), Vec::new()),
            LayoutInput::Regions(inputs) => {
                let mut regions = Vec::with_capacity(inputs.len());
                for input in inputs {
                    regions.push(DgMaterialWallRegion {
                        interval_m: input.interval_m,
                        source: resolve(input.source)?,
                    });
                }
                (None, regions)
            }
        };
        let mut chimneys = Vec::with_capacity(self.chimneys.len());
        for (port, input) in self.chimneys {
            chimneys.push(DgMaterialWallChimney {
                port, layout: input.resolve_cached(gate, load, cache, bytes)?,
            });
        }
        let result = DgMaterialWallLayout { line: self.line, uniform, regions, chimneys };
        if result.complete_band_hz().is_some_and(|band| band[0] >= band[1]) {
            return Err(bad(self.line, "DG wall source regions have no common positive frequency band").into());
        }
        checkpoint(gate, self.line)?;
        Ok(result)
    }
}

impl DgMaterialWallLayout {
    /// Original whole-bore source; absent for rigid or regional assignments.
    pub const fn uniform_record(&self) -> Option<&DgMaterialWallRecord> {
        self.uniform.as_ref()
    }
    /// Original physical regions in inlet-to-outlet order.
    pub fn regions(&self) -> &[DgMaterialWallRegion] { &self.regions }
    /// Explicit material/rigid assignments for original finite side branches.
    pub fn chimneys(&self) -> &[DgMaterialWallChimney] { &self.chimneys }
    /// Intersection of this duct's original applicability bands.
    pub fn band_hz(&self) -> Option<[f64; 2]> {
        self.uniform.iter().map(|source| source.wall().band_hz())
            .chain(self.regions.iter().map(|region| region.source.wall().band_hz()))
            .reduce(|band, next| [band[0].max(next[0]), band[1].min(next[1])])
    }
    /// Common source band across the main bore and all supplied finite branches.
    pub fn complete_band_hz(&self) -> Option<[f64; 2]> {
        self.band_hz().into_iter()
            .chain(self.chimneys.iter().filter_map(|chimney| chimney.layout.band_hz()))
            .reduce(|band, next| [band[0].max(next[0]), band[1].min(next[1])])
    }

    /// Bind each declared finite branch to the SAME local solid-wall owner.
    /// Its original geometry and local gas/diffusion state remain unchanged.
    /// Unknown targets refuse; no port is invented and no main-bore material
    /// is implicitly copied to a chimney. The caller includes all resulting
    /// gas, solid and hereditary states in its complete DG state allowance.
    pub fn bind_chimneys(&self, branches: Vec<DgBranch>, gate: &CancelGate)
        -> Result<Vec<DgBranch>, DgSourceError>
    {
        checkpoint(gate, self.line)?;
        for chimney in &self.chimneys {
            if !branches.iter().any(|branch| branch.port == chimney.port) {
                return Err(bad(self.line, "DG wall chimney target is not an actual finite branch"));
            }
        }
        branches.into_iter().map(|branch| {
            checkpoint(gate, self.line)?;
            let duct = match self.chimneys.iter().find(|source| source.port == branch.port) {
                Some(source) => source.layout.bind(branch.duct, gate)?,
                None => branch.duct,
            };
            Ok(DgBranch { port: branch.port, duct })
        }).collect()
    }

    /// Bind to the original mesh, retaining every wall and gas memory in its
    /// native operator. Refining an aligned mesh repeats source assignments;
    /// an interface that cuts an element refuses instead of mixing materials.
    pub fn bind(&self, duct: DgDuct, gate: &CancelGate) -> Result<DgDuct, DgSourceError> {
        self.bind_with_source(duct, gate).map(|(duct, _)| duct)
    }

    /// Retain the exact original main-bore sources alongside the one operator.
    /// Regional source order follows actual elements after physical admission.
    pub(crate) fn bind_with_source(
        &self, duct: DgDuct, gate: &CancelGate,
    ) -> Result<(DgDuct, Option<DgPortWallSource>), DgSourceError> {
        checkpoint(gate, self.line)?;
        if let Some(source) = &self.uniform {
            return source.wall().clone().bind_dg_controlled(duct, &mut || gate.is_requested())
                .map(|bound| (bound.into_duct(), Some(DgPortWallSource::Uniform(source.wall().clone()))))
                .map_err(DgSourceError::Physics);
        }
        if self.regions.is_empty() { return Ok((duct, None)); }
        let mesh = duct.mesh();
        if self.regions.len() > mesh.elements {
            return Err(bad(self.line, "each DG material region needs at least one original element"));
        }
        let nodes = duct.node_positions();
        let local = mesh.order + 1;
        let length = duct.total_length();
        let tolerance = 32.0 * f64::EPSILON * length.abs();
        let coincides = |a: f64, b: f64| (a - b).abs() <= tolerance;
        if !coincides(self.regions.last().expect("nonempty source regions").interval_m[1], length) {
            return Err(bad(self.line, "DG wall regions must cover the actual complete bore length"));
        }
        let mut walls = Vec::with_capacity(mesh.elements);
        let mut element = 0usize;
        for region in &self.regions {
            checkpoint(gate, self.line)?;
            if element == mesh.elements
                || !coincides(nodes[element * local], region.interval_m[0])
            {
                return Err(bad(self.line, "DG wall material interfaces must align with original element faces"));
            }
            let first = element;
            while element < mesh.elements {
                let right = nodes[(element + 1) * local - 1];
                if right > region.interval_m[1] + tolerance {
                    return Err(bad(self.line, "DG wall material interface cuts an element; refine the declared mesh"));
                }
                walls.push(region.source.wall().clone());
                element += 1;
                if coincides(right, region.interval_m[1]) { break; }
            }
            if element == first || !coincides(nodes[element * local - 1], region.interval_m[1]) {
                return Err(bad(self.line, "DG wall region has no matching physical end face"));
            }
        }
        if element != mesh.elements {
            return Err(bad(self.line, "DG wall regions leave original elements uncovered"));
        }
        RegionalMaterialWallDg::new_controlled(walls, duct, &mut || gate.is_requested())
            .map(|bound| {
                let source = DgPortWallSource::Regional(bound.walls().to_vec());
                (bound.into_duct(), Some(source))
            }).map_err(DgSourceError::Physics)
    }
}
