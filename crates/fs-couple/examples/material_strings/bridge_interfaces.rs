//! Source ingress for the same linear interface used by material shell joints.
//! The bank validates the actual outer string and original body footprint;
//! this adapter only resolves source data and supplies the existing bridge law.
use super::*;
use fs_couple::bowed_string::runtime::reciprocal::bank::MaterialBridgeInterface;
use fs_material::state_point::{ResolvedLinearInterface, resolve_interface_state_point,
    interface_stiffness_per_area_dims, interface_viscosity_per_area_dims,
    LINEAR_INTERFACE_CONSTANT_LOSS_REFERENCE_OMEGA, constant_loss_reference_omega_dims};

#[derive(Clone, Debug)]
pub struct Declaration {
    pub string: usize,
    pub pack: String,
    pub area_m2: f64,
    pub stiffness_property: String,
    pub viscosity_property: String,
}

impl Declaration {
    pub fn validate(&self) -> Result<()> {
        if !self.area_m2.is_finite() || self.area_m2 <= 0.0
            || self.stiffness_property == self.viscosity_property
            || [&self.stiffness_property, &self.viscosity_property].iter()
                .any(|name| name.is_empty() || name.len() > 128) {
            return Err("bridge interfaces need distinct properties and positive physical area".into());
        }
        Ok(())
    }
}

/// One coefficient resolver for ordinary bridges and axial bearings alike.
pub(super) fn resolve(
    declaration: &Declaration, pack: &NormalizedInterfacePack, query: &QueryPoint, omega: f64,
) -> Result<ResolvedLinearInterface> {
    declaration.validate()?;
    let mut requirements = vec![
        Requirement::try_new(&declaration.stiffness_property, interface_stiffness_per_area_dims(), Admit::NonNegative)?,
        Requirement::try_new(&declaration.viscosity_property, interface_viscosity_per_area_dims(), Admit::NonNegative)?,
    ];
    if !pack.card().claims_for(LINEAR_INTERFACE_CONSTANT_LOSS_REFERENCE_OMEGA).is_empty() {
        requirements.push(Requirement::try_new(LINEAR_INTERFACE_CONSTANT_LOSS_REFERENCE_OMEGA,
            constant_loss_reference_omega_dims(), Admit::StrictlyPositive)?);
    }
    let state = resolve_interface_state_point(pack.card(), query, &requirements, Pick::SingleClaimOnly)?;
    let interface = ResolvedLinearInterface::new(state, &declaration.stiffness_property,
        &declaration.viscosity_property, declaration.area_m2)?;
    interface.require_angular_frequency_rad_s(omega)?;
    Ok(interface)
}

#[derive(Clone, Default)]
pub struct Inputs {
    declarations: Vec<Declaration>,
    packs: BTreeMap<String, NormalizedInterfacePack>,
    bindings: Vec<MaterialBridgeInterface>,
    loaded_state: Option<(u64, u64)>,
}
impl core::fmt::Debug for Inputs {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        let packs: Vec<_> = self.packs.iter().map(|(path, pack)|
            (path, pack.content_hash())).collect();
        f.debug_struct("MaterialBridgeInterfaces")
            .field("declarations", &self.declarations)
            .field("packs", &packs)
            .field("bindings", &self.bindings).finish()
    }
}
impl Inputs {
    pub fn new(declarations: Vec<Declaration>) -> Result<Self> {
        for (i, d) in declarations.iter().enumerate() {
            d.validate()?;
            if declarations[..i].iter().any(|p| p.string == d.string) {
                return Err("bridge interfaces need unique strings".into());
            }
        }
        Ok(Self { declarations, ..Self::default() })
    }

    /// All-or-nothing cold source loading at the performance's declared solid
    /// temperature and angular frequency. No source-dependent value is guessed.
    pub fn load(&mut self, temperature: f64, omega: f64, gate: &CancelGate,
        mut loader: impl FnMut(&str) -> Result<NormalizedInterfacePack>) -> Result<()>
    {
        let query = QueryPoint::new().with("T", temperature)?.with("omega", omega)?;
        let mut packs = BTreeMap::new();
        let mut bindings = Vec::with_capacity(self.declarations.len());
        let mut remaining = input::MAX_TOTAL_PACK_BYTES;
        for d in &self.declarations {
            poll(gate)?;
            if !packs.contains_key(&d.pack) {
                let pack = loader(&d.pack)?;
                poll(gate)?;
                let size = pack.to_bytes().len();
                if size > input::MAX_PACK_BYTES || size > remaining {
                    return Err("bridge interface packs exceed the shared source byte budget".into());
                }
                remaining -= size;
                packs.insert(d.pack.clone(), pack);
            }
            let interface = resolve(d, &packs[&d.pack], &query, omega)?;
            bindings.push(MaterialBridgeInterface { string: d.string, interface });
        }
        poll(gate)?;
        self.packs = packs;
        self.bindings = bindings;
        self.loaded_state = Some((temperature.to_bits(), omega.to_bits()));
        Ok(())
    }

    pub fn pack_bytes(&self) -> Result<usize> {
        if self.bindings.len() != self.declarations.len() {
            return Err("load every declared bridge interface before constructing material strings".into());
        }
        self.packs.values().try_fold(0usize, |sum, pack| sum.checked_add(pack.to_bytes().len())
            .ok_or_else(|| "bridge interface source byte count overflow".into()))
    }

    pub fn apply(&self, temperature: f64, omega: f64, strings: &mut [input::Wire]) -> Result<()> {
        self.pack_bytes()?;
        if !self.declarations.is_empty()
            && self.loaded_state != Some((temperature.to_bits(), omega.to_bits())) {
            return Err("bridge interfaces must be resolved at the current declared solid state".into());
        }
        // Prepare every assignment before mutating the candidate geometry.
        let bridges = self.bindings.iter().map(|binding| {
            let original = strings.get(binding.string).ok_or("bridge interface names an absent string")?.bridge;
            Ok((binding.string, SharedPlateBridge::from_interface(original.station_fraction,
                original.rest_extension_m, &binding.interface)?))
        }).collect::<Result<Vec<_>>>()?;
        for (string, bridge) in bridges { strings[string].bridge = bridge; }
        Ok(())
    }

    pub fn bindings(&self) -> &[MaterialBridgeInterface] { &self.bindings }
}
