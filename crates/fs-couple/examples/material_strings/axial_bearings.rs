//! Native source ingress for the shared longitudinal/body bearing owner.
//! Coefficients use the same resolver as every sourced scalar bridge.
use super::*;
use fs_couple::bowed_string::runtime::reciprocal::bank::nonlinear::longitudinal::{
    AxialBearingGeometry, MaterialAxialBearing, bearings::MAX_AXIAL_BEARINGS,
};

#[derive(Clone, Debug)]
pub struct Declaration {
    pub geometry: AxialBearingGeometry,
    pub source: bridge_interfaces::Declaration,
}

#[derive(Clone, Default)]
pub struct Inputs {
    declarations: Vec<Declaration>,
    packs: BTreeMap<String, NormalizedInterfacePack>,
    bindings: Vec<MaterialAxialBearing>,
    loaded_state: Option<(u64, u64)>,
}
impl core::fmt::Debug for Inputs {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        let packs: Vec<_> = self.packs.iter().map(|(path, pack)| (path, pack.content_hash())).collect();
        f.debug_struct("MaterialAxialBearings").field("declarations", &self.declarations)
            .field("packs", &packs).field("bindings", &self.bindings).finish()
    }
}
impl Inputs {
    pub fn new(declarations: Vec<Declaration>) -> Result<Self> {
        if declarations.len() > MAX_AXIAL_BEARINGS {
            return Err("axial bearings exceed the shared owner connection capacity".into());
        }
        for declaration in &declarations {
            declaration.geometry.validate()?;
            declaration.source.validate()?;
            if declaration.source.string != declaration.geometry.string {
                return Err("axial bearing geometry and interface name different strings".into());
            }
        }
        Ok(Self { declarations, ..Self::default() })
    }

    pub fn len(&self) -> usize { self.declarations.len() }

    /// Resolve every supplied interface before publishing any source binding.
    pub fn load(&mut self, temperature: f64, omega: f64, gate: &CancelGate,
        mut loader: impl FnMut(&str) -> Result<NormalizedInterfacePack>) -> Result<()>
    {
        let query = QueryPoint::new().with("T", temperature)?.with("omega", omega)?;
        let mut packs = BTreeMap::new();
        let mut bindings = Vec::with_capacity(self.declarations.len());
        let mut remaining = input::MAX_TOTAL_PACK_BYTES;
        for declaration in &self.declarations {
            poll(gate)?;
            let source = &declaration.source;
            if !packs.contains_key(&source.pack) {
                let pack = loader(&source.pack)?;
                poll(gate)?;
                let size = pack.to_bytes().len();
                if size > input::MAX_PACK_BYTES || size > remaining {
                    return Err("axial bearing packs exceed the shared source byte budget".into());
                }
                remaining -= size;
                packs.insert(source.pack.clone(), pack);
            }
            let interface = bridge_interfaces::resolve(source, &packs[&source.pack], &query, omega)?;
            bindings.push(MaterialAxialBearing { geometry: declaration.geometry.clone(), interface });
        }
        poll(gate)?;
        self.packs = packs;
        self.bindings = bindings;
        self.loaded_state = Some((temperature.to_bits(), omega.to_bits()));
        Ok(())
    }

    pub fn pack_bytes(&self) -> Result<usize> {
        if self.bindings.len() != self.declarations.len() {
            return Err("load every declared axial bearing interface before constructing material strings".into());
        }
        self.packs.values().try_fold(0usize, |sum, pack| sum.checked_add(pack.to_bytes().len())
            .ok_or_else(|| "axial bearing source byte count overflow".into()))
    }

    pub fn bindings(&self, temperature: f64, omega: f64) -> Result<&[MaterialAxialBearing]> {
        self.pack_bytes()?;
        if !self.declarations.is_empty()
            && self.loaded_state != Some((temperature.to_bits(), omega.to_bits())) {
            return Err("axial bearings must be resolved at the current declared solid state".into());
        }
        Ok(&self.bindings)
    }
}
