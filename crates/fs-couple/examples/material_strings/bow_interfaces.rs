//! Ordered surface-data input for the existing simultaneous bow solver.
//! This adapter selects coefficients and their common applicability; it owns
//! no traction equation. The original packs and query receipts remain retained.
//! Original source-qualified environmental coordinates stay held. Coulomb
//! coefficients that also respond to contact temperature require the joint
//! finite-cell contact-heat owner; pressure and speed remain actual mechanical
//! applicability coordinates. No rosin, sorption or wear law is inferred.
use super::*;
use fs_contact::interface_binding::{bind_dry_interface_system_card,
    thermal_friction::{SourceThermalCoulomb, ThermalFrictionAxes, FrictionQueryAxis,
        FrictionQueryAxisKind}};
use fs_couple::render::schedule::force::coupled::contact::friction::interface::DryFrictionBinding;
use fs_material::state_point::{ResolvedInterfaceStatePoint, ResolvedScalarProperty,
    append_held_environment, query_with_source_conventions, resolve_interface_state_point};
use fs_qty::{QuantitySpec, semantic::{QuantityKind, SemanticType, ValueForm}};
use fs_tribo::{
    ApplicabilityRange, DryFrictionApplicability, DryInterfaceSystemCard, InputAuthority,
    InterfaceMedium, InterfaceSystemRef,
};

#[derive(Clone, Copy, Debug)]
pub enum Model {
    Coulomb,
    Stribeck,
    Hyperbolic,
}
impl Model {
    pub const fn property_count(self) -> usize {
        match self {
            Self::Coulomb => 2,
            Self::Stribeck => 4,
            Self::Hyperbolic => 3,
        }
    }
}

/// Property names are in law order: static mu, kinetic mu, then the physical
/// weakening speed [m/s] for Stribeck/Hyperbolic, and (Stribeck only) viscous
/// coefficient per speed [s/m]. A selected name never supplies a fallback value.
#[derive(Clone, Debug)]
pub struct Declaration {
    pub track: String,
    pub pack: String,
    pub temperature: f64,
    pub area: f64,
    pub pressure: f64,
    pub speed: f64,
    pub model: Model,
    pub names: Vec<String>,
}
#[derive(Clone)]
struct Bound {
    state: ResolvedInterfaceStatePoint,
    binding: DryFrictionBinding,
}
#[derive(Clone)]
pub struct Inputs {
    research: Option<FrictionLaw>,
    layers: Option<shear_layers::Inputs>,
    declarations: Vec<Declaration>,
    packs: BTreeMap<String, NormalizedInterfacePack>,
    environment: BTreeMap<String, Vec<(String, f64)>>,
    bound: BTreeMap<String, Bound>,
}
impl core::fmt::Debug for Inputs {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        let packs: Vec<_> = self
            .packs
            .iter()
            .map(|(path, pack)| (path, pack.content_hash()))
            .collect();
        let states: Vec<_> = self
            .bound
            .iter()
            .map(|(track, b)| (track, b.state.identity()))
            .collect();
        f.debug_struct("BowFriction")
            .field("research_override", &if self.layers.is_none() { self.research.as_ref() } else { None })
            .field("material_shear_layers", &self.layers)
            .field("declarations", &self.declarations)
            .field("packs", &packs)
            .field("environment", &self.environment)
            .field("states", &states)
            .finish()
    }
}
impl Inputs {
    pub fn research(law: FrictionLaw) -> Result<Self> {
        law.regularized_traction_1d(0.0, 0.0, 1.0)?;
        Ok(Self {
            research: Some(law),
            layers: None,
            declarations: Vec::new(),
            packs: BTreeMap::new(),
            environment: BTreeMap::new(),
            bound: BTreeMap::new(),
        })
    }
    /// The zero dry law constructs port topology only. Every layer declaration
    /// must resolve and bind before scheduling; no fallback renderer is exposed.
    pub fn layers(declarations: Vec<shear_layers::Declaration>) -> Result<Self> {
        let layers = shear_layers::Inputs::new(declarations)?;
        let mut inputs = Self::research(FrictionLaw::Coulomb { static_mu: 0.0, kinetic_mu: 0.0 })?;
        inputs.layers = Some(layers);
        Ok(inputs)
    }
    pub fn shear_layers(&self) -> Option<&shear_layers::Inputs> {
        self.layers.as_ref()
    }
    pub fn validate_thermal_mode(&self, dry_thermal: bool) -> Result<()> {
        if !dry_thermal && self.bound.values().any(|bound|
            bound.binding.temperature_source().is_some()) {
            return Err("temperature-dependent bow friction requires explicit finite-cell dry-contact heat routes".into());
        }
        Ok(())
    }
    pub fn interfaces(declarations: Vec<Declaration>) -> Result<Self> {
        for (i, d) in declarations.iter().enumerate() {
            if d.track.is_empty()
                || d.track.len() > 128
                || declarations[..i].iter().any(|p| p.track == d.track)
                || !d.speed.is_finite()
                || d.speed < 0.0
                || d.names.len() != d.model.property_count()
                || d.names
                    .iter()
                    .enumerate()
                    .any(|(i, n)| n.is_empty() || n.len() > 128 || d.names[..i].contains(n))
            {
                return Err("bow interfaces need distinct tracks/properties and nonnegative query slip speed".into());
            }
        }
        Ok(Self {
            research: None,
            layers: None,
            declarations,
            packs: BTreeMap::new(),
            environment: BTreeMap::new(),
            bound: BTreeMap::new(),
        })
    }
    pub fn validate_tracks(&self, schedule: &GestureSchedule) -> Result<()> {
        if let Some(layers) = &self.layers { return layers.validate_tracks(schedule); }
        if self.research.is_some() {
            return Ok(());
        }
        if self.declarations.len() != schedule.tracks().len()
            || schedule
                .tracks()
                .iter()
                .any(|t| !self.declarations.iter().any(|d| d.track == t.id))
        {
            return Err(
                "source-bound friction requires exactly one named interface for every bow track"
                    .into(),
            );
        }
        Ok(())
    }

    /// Cold, all-or-nothing loading; repeated paths decode once. An absent or
    /// incompatible source cannot silently fall back to research coefficients.
    pub fn load(
        &mut self,
        gate: &CancelGate,
        loader: impl FnMut(&str) -> Result<NormalizedInterfacePack>,
    ) -> Result<()> {
        self.load_with_environment(&BTreeMap::new(), gate, loader)
    }

    pub fn load_with_environment(
        &mut self,
        environment: &BTreeMap<String, Vec<(String, f64)>>,
        gate: &CancelGate,
        mut loader: impl FnMut(&str) -> Result<NormalizedInterfacePack>,
    ) -> Result<()> {
        let mut packs = BTreeMap::new();
        let mut retained_environment = BTreeMap::new();
        let mut bound = BTreeMap::new();
        let mut remaining = input::MAX_TOTAL_PACK_BYTES;
        for d in &self.declarations {
            poll(gate)?;
            if !packs.contains_key(&d.pack) {
                let pack = loader(&d.pack)?;
                poll(gate)?;
                let size = pack.to_bytes().len();
                if size > input::MAX_PACK_BYTES || size > remaining {
                    return Err("bow interface packs exceed the shared source byte budget".into());
                }
                remaining -= size;
                packs.insert(d.pack.clone(), pack);
            }
            let held = environment.get(&d.pack).map_or(&[][..], Vec::as_slice);
            let source = resolve(d, &packs[&d.pack], held)?;
            if let Some(held) = environment.get(&d.pack) {
                retained_environment.insert(d.pack.clone(), held.clone());
            }
            bound.insert(d.track.clone(), source);
        }
        poll(gate)?;
        self.packs = packs;
        self.environment = retained_environment;
        self.bound = bound;
        Ok(())
    }
    pub fn require_environment(&self, environment: &BTreeMap<String, Vec<(String, f64)>>) -> Result<()> {
        if self.declarations.iter().any(|d| self.environment.get(&d.pack) != environment.get(&d.pack)) {
            return Err("bow interfaces must be resolved at the current declared held environment".into());
        }
        Ok(())
    }
    pub fn has_source(&self, name: &str) -> bool { self.packs.contains_key(name) }
    pub fn pack_bytes(&self) -> Result<usize> {
        if self
            .declarations
            .iter()
            .any(|d| !self.bound.contains_key(&d.track))
        {
            return Err(
                "load every declared bow interface before constructing material strings".into(),
            );
        }
        self.packs.values().try_fold(0usize, |n, p| {
            n.checked_add(p.to_bytes().len())
                .ok_or_else(|| "bow interface byte count overflow".into())
        })
    }
    pub fn binding(&self, track: &str) -> Result<Option<&DryFrictionBinding>> {
        if self.research.is_some() {
            return Ok(None);
        }
        Ok(Some(
            &self
                .bound
                .get(track)
                .ok_or("bow interface is missing or not loaded")?
                .binding,
        ))
    }
    pub fn law(&self, track: &str) -> Result<&FrictionLaw> {
        if let Some(law) = &self.research {
            return Ok(law);
        }
        Ok(self
            .binding(track)?
            .ok_or("missing source-bound bow law")?
            .card()
            .friction_law())
    }
    /// Ordered source participants, retained before reduction to a scalar law.
    pub fn surface_materials(&self, track: &str) -> Result<&[MaterialStateId; 2]> {
        Ok(self.bound.get(track)
            .ok_or("dry contact heat requires an original source-bound interface")?
            .state.surface_materials())
    }
    pub fn validate_surfaces(&self, spec: &input::Spec, catalog: &Catalog) -> Result<()> {
        if self.research.is_some() {
            return Ok(());
        }
        for track in spec.schedule.tracks() {
            let member = fs_scenario::gesture::path::legacy_member(&track.target, "bow")
                .ok_or("interface binding requires a material-string bow track")?;
            let wire = spec
                .strings
                .get(member)
                .ok_or("bow interface names an absent string")?;
            let outer = match &wire.section {
                input::WireSection::Measured { .. } => return Err("whole-string measurements do not identify the outer contact material".into()),
                input::WireSection::Homogeneous => &wire.pack,
                input::WireSection::CloseWound { pack, .. } => pack,
                input::WireSection::Bonded { sleeves, .. } => {
                    &sleeves.last().ok_or("bonded surface is absent")?.pack
                }
            };
            let material = catalog
                .get(outer)
                .ok_or("outer string material has not been resolved")?
                .card
                .id();
            let source = self
                .bound
                .get(&track.id)
                .ok_or("bow interface is not loaded")?;
            if &source.state.surface_materials()[0] != material {
                return Err(format!(
                    "bow interface {:?} surface A does not match its string's outer material",
                    track.id
                )
                .into());
            }
        }
        Ok(())
    }
}

fn resolve(d: &Declaration, pack: &NormalizedInterfacePack, environment: &[(String, f64)]) -> Result<Bound> {
    if pack.card().medium() != "dry" {
        return Err("the dry bow law requires an explicitly dry ordered interface card".into());
    }
    let requirements = d
        .names
        .iter()
        .enumerate()
        .map(|(i, name)| {
            let (dims, admit) = match i {
                2 => (fs_qty::Velocity::DIMS, Admit::StrictlyPositive),
                3 => (
                    Dims::NONE.saturating_minus(fs_qty::Velocity::DIMS),
                    Admit::NonNegative,
                ),
                _ => (Dims::NONE, Admit::NonNegative),
            };
            Requirement::try_new(name, dims, admit)
        })
        .collect::<std::result::Result<Vec<_>, _>>()?;
    let claims = pack.claims_pack().claims();
    let static_quantity = |kind| QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static));
    let point = query_with_source_conventions(claims, &[
        ("T", static_quantity(QuantityKind::AbsoluteTemperature), d.temperature),
        ("p", static_quantity(QuantityKind::Pressure), d.pressure),
        ("v", QuantitySpec::dimensional(fs_qty::Velocity::DIMS), d.speed),
    ])?;
    let point = append_held_environment(claims, point, environment)?;
    let state =
        resolve_interface_state_point(pack.card(), &point, &requirements, Pick::SingleClaimOnly)?;
    // Temperature, pressure and actual slip belong to physical owners. Only
    // explicitly supplied, source-qualified extra coordinates may be frozen.
    let held_response = |property: &ResolvedScalarProperty| {
        if property.source_axes().is_empty() {
            property.answer().receipt.decision == fs_matdb::EvaluationDecision::ConstantWithinValidity
        } else {
            property.source_axes().iter().all(|(axis, _)|
                environment.iter().any(|(name, _)| name == axis))
        }
    };
    let thermal = matches!(d.model, Model::Coulomb)
        && state.properties().iter().any(|property| !held_response(property));
    let mut ranges = [(0.0f64, f64::MAX); 3];
    for property in state.properties() {
        if !thermal && !held_response(property) {
            return Err(
                "evolving bow coefficients require the sourced Coulomb thermal owner; only declared environmental coordinates may be held"
                    .into(),
            );
        }
        let validity = &property.answer().evidence.model.validity;
        for (i, name) in ["T", "p", "v"].iter().enumerate() {
            let (lo, hi) = validity
                .bound(name)
                .ok_or("every bow coefficient needs explicit T/p/v applicability")?;
            if !lo.is_finite()
                || !hi.is_finite()
                || hi < lo
                || (i == 0 && lo <= 0.0)
                || (i != 0 && lo < 0.0)
            {
                return Err("bow applicability needs positive temperature and nonnegative finite pressure/speed bounds".into());
            }
            ranges[i].0 = ranges[i].0.max(lo);
            ranges[i].1 = ranges[i].1.min(hi);
        }
    }
    let value = |i: usize| {
        state
            .property(&d.names[i])
            .map(|p| p.value_si())
            .ok_or("missing resolved friction coefficient")
    };
    let law = match d.model {
        Model::Coulomb => FrictionLaw::Coulomb {
            static_mu: value(0)?,
            kinetic_mu: value(1)?,
        },
        Model::Stribeck => FrictionLaw::Stribeck {
            static_mu: value(0)?,
            kinetic_mu: value(1)?,
            characteristic_speed: value(2)?,
            viscous_per_speed: value(3)?,
        },
        Model::Hyperbolic => FrictionLaw::Hyperbolic {
            static_mu: value(0)?,
            kinetic_mu: value(1)?,
            characteristic_speed: value(2)?,
        },
    };
    let range = |i: usize| ApplicabilityRange::new(ranges[i].0, ranges[i].1);
    let domain = DryFrictionApplicability::new(range(0)?, range(1)?, range(2)?)?;
    if thermal {
        let point = point.with_quantity("v", fs_qty::QuantitySpec::dimensional(fs_qty::Velocity::DIMS), d.speed)?;
        let pins = requirements.iter().map(|requirement| {
            let property = state.property(requirement.name()).expect("resolved required property");
            (requirement.name().to_owned(), property.answer().receipt.selected)
        }).collect();
        let source = SourceThermalCoulomb::resolve_with_query_axes(
            bind_dry_interface_system_card(pack.card(), InputAuthority::CallerDeclared)?, point,
            requirements[0].clone(), requirements[1].clone(), Pick::PinnedByProperty(pins),
            ThermalFrictionAxes { temperature: "T".into(), nominal_pressure: Some("p".into()) },
            vec![FrictionQueryAxis { name: "v".into(), kind: FrictionQueryAxisKind::SlidingSpeed }],
        )?;
        let state = source.initial_source().clone();
        let binding = DryFrictionBinding::new_temperature_source(source, domain, d.area)?;
        return Ok(Bound { state, binding });
    }
    let reference = InterfaceSystemRef::new(
        format!("{}", pack.card().content_hash()),
        pack.card().history(),
        format!("pack:{}; state:{}", pack.content_hash(), state.identity()),
        InputAuthority::CallerDeclared,
        InterfaceMedium::Dry,
    )?;
    let card = DryInterfaceSystemCard::new(reference, law, domain)?;
    let binding = DryFrictionBinding::new(card, d.temperature, d.area)?;
    Ok(Bound { state, binding })
}
