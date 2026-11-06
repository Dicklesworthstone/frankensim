//! `friction rosin-thermal-sw2000 ROSIN_PACK CONTEXT RADIUS_M REF_N LAYER_M STRING_SURFACE_PACK CONTEXT
//! (HAIR_PACK CONTEXT | rosin-per-wg2025-eq5) CALIBRATION RIG_PACK CONTEXT`:
//! Smith & Woodhouse (2000) thermal rosin on the bow's contact, in place of the
//! steady `rosin-wg2025-eq13` row (README "Thermal rosin").
//!
//! `RADIUS_M` is the WHOLE contact's radius at the whole bow load `REF_N`, and
//! `LAYER_M` the actively deforming rosin layer. The rosin pack gives density
//! and specific heat. The string surface, stationary relative to the contact,
//! and the hair, which moves through it, give conductivity, density and
//! specific heat. The steady sliding curve and its rig's geometry come only
//! from the named `CALIBRATION` (`fs_tribo::rosin_thermal::RosinCalibration`),
//! never from numbers on the row; `RIG_PACK` is the solid that rig's heat
//! conducted into. Every bow contact cell is one cell of this one contact
//! (`fs_couple ... friction::multiple::rosin_thermal`).
//!
//! Each pack names its source context on the card, as a measured string section
//! does (MEASURED_STRINGS.md): `COUNT AXIS VALUE ...`. `0` reads the pack at the
//! performance's `state` (T, omega) and its `material-state` coordinates. A
//! nonzero count reads it at exactly the named source coordinates, without the
//! performance's T: `temperature_known 0` says the source stated no test
//! temperature. It is not zero kelvin, room temperature, or the performance's T.
//! Such a pack's values enter the law as the constants Smith & Woodhouse used.
//! `rosin-per-wg2025-eq5` in the hair's place puts the rosin pack's constants on
//! the hair side, which is Woodhouse & Galluzzo's (2025) own stated assumption in
//! their eq. (5); no printed keratin value is on file.
use super::*;
use fs_couple::material_source::{MaterialRecord, MAX_RECORD_AXES};
use fs_couple::render::schedule::force::coupled::contact::friction::multiple::rosin_thermal::RosinThermalPorts;
use fs_tribo::rosin_thermal::{RosinCalibration, RosinContact, SmithWoodhouse, SmithWoodhouseParams, ThermalSolid};

/// The calibrations on file, by name.
/// `wg2025-eq13`: Woodhouse & Galluzzo (2025) eq. (13) on Smith & Woodhouse's
/// (2000) steady-sliding rig, whose solid was perspex.
const CALIBRATIONS: [&str; 1] = ["wg2025-eq13"];

/// The literal that, in the hair's place, puts the rosin's constants on the hair.
pub const ROSIN_HAIR: &str = "rosin-per-wg2025-eq5";

fn on_file(name: &str) -> Result<()> {
    if CALIBRATIONS.contains(&name) { return Ok(()); }
    Err(format!("rosin-thermal-sw2000: no calibration on file named {name:?}; on file: {}", CALIBRATIONS.join(", ")).into())
}

/// The named curve and rig, with the rig's solid read from `RIG_PACK`.
fn calibration(name: &str, rig_solid: ThermalSolid) -> Result<RosinCalibration> {
    match name {
        "wg2025-eq13" => Ok(RosinCalibration::woodhouse_galluzzo_2025_eq13(rig_solid)?),
        _ => { on_file(name)?; Err("rosin-thermal-sw2000: a calibration on file has no constructor".into()) }
    }
}

/// One pack on the row and the source context it is read at, as written on the card.
#[derive(Clone, Debug, PartialEq)]
pub struct PackRead {
    pub path: String,
    /// Empty: the performance's state. Otherwise exactly these source coordinates.
    pub context: Vec<(String, f64)>,
}

impl PackRead {
    /// # Errors
    /// More than eight coordinates, or a repeated, nonfinite, `T` or `omega` axis.
    pub fn new(path: String, context: Vec<(String, f64)>) -> Result<Self> {
        if context.len() > MAX_RECORD_AXES {
            return Err(format!("rosin-thermal-sw2000: a source context takes at most {MAX_RECORD_AXES} coordinates").into());
        }
        for (index, (axis, value)) in context.iter().enumerate() {
            if axis == "T" || axis == "omega" || !value.is_finite()
                || context[..index].iter().any(|(prior, _)| prior == axis)
            {
                return Err("rosin-thermal-sw2000: a source context takes distinct finite coordinates; T and omega belong to the state row".into());
            }
        }
        Ok(Self { path, context })
    }

    /// The pack's properties at its own context, in the assembly's budgeted catalog.
    fn read(&self, spec: &input::Spec, catalog: &mut Catalog,
        loader: &mut impl FnMut(&str) -> Result<MaterialPack>, requirements: &[Requirement])
        -> Result<ResolvedMaterialStatePoint>
    {
        if self.context.is_empty() {
            return resolve_properties(&self.path, catalog, loader, spec, requirements);
        }
        if spec.material_states.contains_key(&self.path) {
            return Err(format!("rosin-thermal-sw2000: {} names its source context on the row; a material-state for it would be a second context", self.path).into());
        }
        load_pack(&self.path, catalog, loader)?;
        let card = &catalog[&self.path].card;
        let query = MaterialRecord { pack: self.path.clone(), state: card.id().process.clone(),
            axes: self.context.clone() }.query_point(card)?;
        Ok(resolve_material_state_point(card, &query, requirements, Pick::SingleClaimOnly)?)
    }

    fn describe(&self) -> String {
        if self.context.is_empty() {
            return format!("{:?} at the performance's state", self.path);
        }
        let axes = self.context.iter().map(|(axis, value)| format!("{axis}={value}")).collect::<Vec<_>>().join(",");
        format!("{:?} at source context {axes}, without the performance's T: its values enter as Smith & Woodhouse's constants", self.path)
    }
}

/// The hair side: a pack, or the rosin's constants as WG2025 eq. (5) assumes.
#[derive(Clone, Debug, PartialEq)]
pub enum Hair {
    Pack(PackRead),
    RosinPerWg2025Eq5,
}

/// The parsed row; its calibration name is checked on parse.
#[derive(Clone, Debug)]
pub struct Row {
    pub rosin: PackRead,
    pub radius_m: f64,
    pub reference_normal_n: f64,
    pub layer_m: f64,
    pub string: PackRead,
    pub hair: Hair,
    pub calibration_name: String,
    pub rig: PackRead,
}

impl Row {
    /// Check the calibration name first: an unknown name refuses before any pack is read.
    #[allow(clippy::too_many_arguments)]
    pub fn new(rosin: PackRead, radius_m: f64, reference_normal_n: f64, layer_m: f64,
        string: PackRead, hair: Hair, calibration_name: String, rig: PackRead) -> Result<Self> {
        on_file(&calibration_name)?;
        Ok(Self { rosin, radius_m, reference_normal_n, layer_m, string, hair, calibration_name, rig })
    }

    /// Read the packs into the assembly's budgeted catalog, each at its own
    /// context, and build the whole contact and the calibration.
    pub fn resolve(&self, spec: &input::Spec, catalog: &mut Catalog,
        loader: &mut impl FnMut(&str) -> Result<MaterialPack>) -> Result<Resolved> {
        let density = Requirement::try_new("density", Density::DIMS, Admit::StrictlyPositive)?;
        let heat = Requirement::try_new(SPECIFIC_HEAT_CAPACITY_PROPERTY, SPECIFIC_HEAT_CAPACITY_DIMS, Admit::StrictlyPositive)?;
        let conduct = Requirement::try_new(THERMAL_CONDUCTIVITY_PROPERTY, THERMAL_CONDUCTIVITY_DIMS, Admit::StrictlyPositive)?;
        let value = |state: &ResolvedMaterialStatePoint, pack: &str, name: &str| -> Result<f64> {
            Ok(state.property(name).ok_or_else(|| format!("rosin-thermal-sw2000: {pack} resolved no {name}"))?.value_si())
        };
        let rosin = self.rosin.read(spec, catalog, loader, &[density.clone(), heat.clone()])?;
        let mut solid = |pack: &PackRead| -> Result<ThermalSolid> {
            let state = pack.read(spec, catalog, loader, &[conduct.clone(), density.clone(), heat.clone()])?;
            Ok(ThermalSolid {
                conductivity_w_per_m_k: value(&state, &pack.path, THERMAL_CONDUCTIVITY_PROPERTY)?,
                density_kg_per_m3: value(&state, &pack.path, "density")?,
                specific_heat_j_per_kg_k: value(&state, &pack.path, SPECIFIC_HEAT_CAPACITY_PROPERTY)?,
            })
        };
        let stationary = solid(&self.string)?;
        let moving = match &self.hair { Hair::Pack(pack) => solid(pack)?, Hair::RosinPerWg2025Eq5 => solid(&self.rosin)? };
        let rig_solid = solid(&self.rig)?;
        Ok(Resolved {
            row: self.clone(),
            calibration: calibration(&self.calibration_name, rig_solid)?,
            rosin_density_kg_per_m3: value(&rosin, &self.rosin.path, "density")?,
            rosin_specific_heat_j_per_kg_k: value(&rosin, &self.rosin.path, SPECIFIC_HEAT_CAPACITY_PROPERTY)?,
            contact: RosinContact { radius_m: self.radius_m, reference_normal_n: self.reference_normal_n,
                layer_thickness_m: self.layer_m, stationary, moving },
        })
    }
}

/// The row with its packs read.
#[derive(Clone, Debug)]
pub struct Resolved {
    row: Row,
    calibration: RosinCalibration,
    rosin_density_kg_per_m3: f64,
    rosin_specific_heat_j_per_kg_k: f64,
    contact: RosinContact,
}

impl Resolved {
    /// The law on the bank's own step, with one history step kept for every
    /// sample of the performance. Validated here, before any output exists.
    pub fn params(&self, step_s: f64, samples: u64) -> Result<SmithWoodhouseParams> {
        let c = self.calibration;
        let params = SmithWoodhouseParams {
            rosin_density_kg_per_m3: self.rosin_density_kg_per_m3,
            rosin_specific_heat_j_per_kg_k: self.rosin_specific_heat_j_per_kg_k,
            steady_fit: c.curve, calibration: c.rig, calibration_speed_max_mps: c.speed_max_mps,
            contact: self.contact, step_s,
            history_capacity: usize::try_from(samples)?,
        };
        SmithWoodhouse::new(params)?;
        Ok(params)
    }

    /// The setup receipt, printed before the render.
    pub fn describe(&self, params: &SmithWoodhouseParams, cells: usize) -> Result<()> {
        let law = SmithWoodhouse::new(*params)?;
        let c = &params.contact;
        let initial_mu = core::f64::consts::PI * c.radius_m * c.radius_m * law.yield_stress_pa(0.0)? / c.reference_normal_n;
        let hair = match &self.row.hair {
            Hair::Pack(pack) => format!("hair_pack={} (moving)", pack.describe()),
            Hair::RosinPerWg2025Eq5 => format!("hair={ROSIN_HAIR}: the rosin pack's constants on the hair side (moving), \
                Woodhouse & Galluzzo (2025) eq. (5)'s own assumption; no printed keratin value is on file"),
        };
        eprintln!("rosin_thermal=smith-woodhouse-2000; calibration={:?}; curve={:?}; rig={:?}; rig_pack={}; calibration_speed_max_m_s={}; rosin_pack={}; string_surface_pack={} (stationary); {hair}; whole_radius_m={}; whole_reference_N={}; layer_m={}; cells={cells}; rosin_density_kg_m3={}; rosin_specific_heat_J_kg_K={}; string_surface_solid={:?}; hair_solid={:?}; temperature_limit_K={}; step_s={:e}; history_steps={}; initial_mu={initial_mu}; each cell takes the step-start Coulomb mu_i = pi a^2 k_y(T_i)/N_ref; its heat advances on accepted loaded steps only; each cell's history starts at its first loaded step (its t = 0), and after its lift its force is mu x 0 = 0 and its history is frozen; a second landing refuses, since off-string cooling is unsourced",
            self.row.calibration_name, params.steady_fit, params.calibration, self.row.rig.describe(), params.calibration_speed_max_mps,
            self.row.rosin.describe(), self.row.string.describe(),
            c.radius_m, c.reference_normal_n, c.layer_thickness_m,
            params.rosin_density_kg_per_m3, params.rosin_specific_heat_j_per_kg_k,
            c.stationary, c.moving, law.temperature_limit_k(), params.step_s, params.history_capacity);
        Ok(())
    }
}

/// The receipt after the render, or after a refusal part way through. Steps
/// are 0-based accepted samples; a lifted cell's temperature is frozen.
pub fn receipt(ports: &RosinThermalPorts) {
    eprintln!("rosin_thermal_accepted_steps={}; heat_steps={:?}; phases={:?}; temperature_rise_K={:?}; next_coefficients={:?}; heat_seconds={}; heat_seconds is the wall time spent on the heat histories, a receipt, not physics",
        ports.accepted_steps(), ports.heat_steps(), ports.phases(), ports.temperatures_k(), ports.coefficients(), ports.heat_seconds());
}
