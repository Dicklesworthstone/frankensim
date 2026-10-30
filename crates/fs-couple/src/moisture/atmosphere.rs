//! Isothermal moisture exchange with a finite, well-mixed atmosphere.
//!
//! One scalar implicit balance couples the ORIGINAL slab diffusion solves to
//! one atmosphere's actual water inventory. Sorption is a pinned material
//! property in kg water/kg dry solid as a function of relative humidity; RH,
//! dry-basis moisture and gas composition are never interchangeable inputs.
//! The existing gas owner resolves the accepted composition, pressure,
//! transport and acoustic medium. Carrier species and dry skeleton are fixed.
//!
//! Arbitrarily many admitted specimens may share this same inventory; every
//! trial starts from all their accepted fields and uses one common endpoint RH.
//! Each face may contact this atmosphere or be sealed. Temperature and volume
//! are prescribed; a thermostat maintains the isothermal condition. As in the
//! original moisture owner, bound-water caloric storage/sorption heat are NOT
//! represented, so the gas energy change below is not a total heat balance.
//! There is no condensation, hysteresis, gas momentum or moving interface.

use super::surface::{DryBasisSurfaceTransfer, ResolvedMoistureFace, admit_environment};
use super::{
    MaterialCard, MaterialPropertySelection, MoistureError, MoistureSlab, MoistureStepConfig,
    MoistureStepReport, QueryPoint, ResolvedMaterialStatePoint, ScalarAdmissibility,
    ScalarPropertyRequirement, poll, positive, quantity, resolve_material_state_point,
};
use fs_exec::Cx;
use fs_material::thermochemical_gas::{
    ThermochemicalGasModel,
    energy::{GasEnergyError, parcel::ClosedGasParcel},
};
use fs_qty::semantic::{MoistureBasis, QuantityKind};
use fs_thermochem::{Composition, CompositionBasis, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};

/// Source-backed equilibrium surface content on an explicit dry-mass basis.
/// This is an isotherm, not instantaneous equilibration of the slab interior.
#[derive(Clone, Debug, PartialEq)]
pub struct DryBasisSorption {
    card: MaterialCard,
    source: ResolvedMaterialStatePoint,
    property: String,
    surface_transfer: Option<DryBasisSurfaceTransfer>,
}
impl DryBasisSorption {
    /// Pin one property with semantic kg-water/kg-dry-solid output. The query
    /// must carry typed absolute temperature and relative humidity. Total gas
    /// pressure dependence and additional condition coordinates are excluded
    /// from this two-coordinate dilute-vapor activity model.
    pub fn new(
        card: MaterialCard,
        point: QueryPoint,
        property: &str,
        selection: MaterialPropertySelection,
    ) -> Result<Self, MoistureError> {
        let rh = quantity(QuantityKind::Moisture(MoistureBasis::RelativeHumidity));
        admit_environment(&point)?;
        let requirement = ScalarPropertyRequirement::try_with_quantity(
            property,
            quantity(QuantityKind::Moisture(MoistureBasis::DryMass)),
            ScalarAdmissibility::NonNegative,
        )?;
        let source = resolve_material_state_point(&card, &point, &[requirement], selection)?;
        if !source
            .property(property)
            .expect("required sorption property")
            .has_only_owned_and_fixed_external_conditions("relative_humidity", rh)
        {
            return Err(MoistureError::Invalid(
                "sorption needs static source support for the owned RH coordinate",
            ));
        }
        Ok(Self {
            card,
            source,
            property: property.into(),
            surface_transfer: None,
        })
    }

    /// Limit exchange by the original source-backed surface emission law.
    /// Both laws are re-evaluated at the actual environment RH on each trial.
    /// This does not change the equilibrium isotherm or create a surface store.
    pub fn with_surface_transfer(
        mut self,
        transfer: DryBasisSurfaceTransfer,
    ) -> Result<Self, MoistureError> {
        let point = self.source.query();
        transfer.at(
            point.axes()["temperature"],
            point.axes()["relative_humidity"],
        )?;
        self.surface_transfer = Some(transfer);
        Ok(self)
    }

    pub(super) fn resolve(
        &self,
        temperature: f64,
        rh: f64,
    ) -> Result<ResolvedMoistureFace, MoistureError> {
        let (equilibrium, sorption) = self.at(temperature, rh)?;
        let (coefficient_m_s, transfer) = self
            .surface_transfer
            .as_ref()
            .map(|law| law.at(temperature, rh))
            .transpose()?
            .unzip();
        Ok(ResolvedMoistureFace {
            equilibrium: Some(equilibrium),
            coefficient_m_s,
            sorption: Some(sorption),
            transfer,
        })
    }

    fn at(
        &self,
        temperature: f64,
        rh: f64,
    ) -> Result<(f64, ResolvedMaterialStatePoint), MoistureError> {
        if self.source.query().axes()["temperature"].to_bits() != temperature.to_bits()
            || !rh.is_finite()
            || !(0.0..=1.0).contains(&rh)
        {
            return Err(MoistureError::Invalid(
                "gas, sorption and slab need the same fixed temperature and RH in [0,1]",
            ));
        }
        let point = self
            .source
            .query()
            .with("relative_humidity", rh)
            .map_err(|_| MoistureError::Invalid("invalid relative humidity"))?;
        let resolved = self.source.at_query_point(&self.card, &point)?;
        let value = resolved
            .property(&self.property)
            .expect("pinned sorption property")
            .value_si();
        Ok((value, resolved))
    }
}

/// Original material, gas or diffusion refusal; no partial endpoint is published.
#[derive(Debug)]
pub enum AtmosphereMoistureError {
    Moisture(MoistureError),
    Gas(GasEnergyError),
    Invalid(&'static str),
    Unbracketed { residual_kg: [f64; 2] },
    NotConverged { residual_kg: f64, iterations: usize },
    MassBalance { residual_kg: f64, tolerance_kg: f64 },
}
impl core::fmt::Display for AtmosphereMoistureError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite moisture atmosphere: {self:?}")
    }
}
impl std::error::Error for AtmosphereMoistureError {}
impl From<MoistureError> for AtmosphereMoistureError {
    fn from(value: MoistureError) -> Self {
        Self::Moisture(value)
    }
}
impl From<GasEnergyError> for AtmosphereMoistureError {
    fn from(value: GasEnergyError) -> Self {
        Self::Gas(value)
    }
}
impl From<fs_material::MaterialError> for AtmosphereMoistureError {
    fn from(value: fs_material::MaterialError) -> Self {
        Self::Moisture(MoistureError::Constitutive(value))
    }
}
fn invalid(what: &'static str) -> AtmosphereMoistureError {
    AtmosphereMoistureError::Invalid(what)
}

/// Bounded scalar solve of endpoint vapor + endpoint bound water = initial water.
#[derive(Clone, Copy, Debug)]
pub struct AtmosphereMoistureSolve {
    /// Source-supported closed RH interval in [0,1], containing the initial RH.
    /// All proposed diffusion/source states in this bracket must be admissible.
    pub relative_humidity_bracket: [f64; 2],
    /// Maximum bisection iterations, in 1..=256; initial and two endpoint trials
    /// are extra. An already balanced initial RH needs only one native solve.
    pub maximum_iterations: usize,
    /// Absolute total-water, total-mass AND per-carrier mass tolerance [kg].
    pub mass_tolerance_kg: f64,
}

/// One native gas inventory and its retained complete source catalog. The catalog
/// retains vapor data even after exact depletion; it is not another species store.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteMoistureAtmosphere {
    gas: ClosedGasParcel,
    catalog: ThermochemicalGasModel,
    water: usize,
    saturation_mass_kg: f64,
}

/// Accepted finite-reservoir exchange; positive transfer enters the solid.
#[derive(Clone, Debug)]
pub struct AtmosphereMoistureStep {
    pub transport: MoistureStepReport,
    pub sorption: [Option<ResolvedMaterialStatePoint>; 2],
    pub relative_humidity: [f64; 2],
    pub vapor_mass_kg: [f64; 2],
    pub water_to_solid_kg: f64,
    pub water_balance_residual_kg: f64,
    pub maximum_carrier_mass_residual_kg: f64,
    /// Actual change of gas m*u at prescribed T. This excludes bound-water
    /// storage and sorption heat, and is NOT the thermostat's total heat duty.
    pub gas_internal_energy_change_j: f64,
    pub iterations: usize,
}

/// Gas and solid endpoints accepted together after source and mass checks.
#[derive(Clone, Debug)]
pub struct AtmosphereMoistureCandidate {
    atmosphere: FiniteMoistureAtmosphere,
    moisture: MoistureSlab,
    step: AtmosphereMoistureStep,
}
impl AtmosphereMoistureCandidate {
    pub const fn atmosphere(&self) -> &FiniteMoistureAtmosphere {
        &self.atmosphere
    }
    pub const fn moisture(&self) -> &MoistureSlab {
        &self.moisture
    }
    pub const fn step(&self) -> &AtmosphereMoistureStep {
        &self.step
    }
    pub fn into_parts(
        self,
    ) -> (
        FiniteMoistureAtmosphere,
        MoistureSlab,
        AtmosphereMoistureStep,
    ) {
        (self.atmosphere, self.moisture, self.step)
    }
}

/// One physical specimen exposed to a shared atmosphere. Stable, distinct IDs
/// name separate physical inventories, including geometrically identical ones.
/// Each specimen retains its own original material, diffusion and surface laws.
pub struct MoistureParticipant<'a> {
    pub id: u64,
    pub slab: &'a MoistureSlab,
    pub faces: &'a [Option<DryBasisSorption>; 2],
    pub transport: MoistureStepConfig,
}

/// One accepted specimen and its native local transfer/source report.
#[derive(Clone, Debug)]
pub struct SharedMoistureSolid {
    pub id: u64,
    pub moisture: MoistureSlab,
    pub transport: MoistureStepReport,
}

/// One common gas inventory debit, with all internal vapor/solid transfers
/// canceled. The gas energy account excludes sorption heat as in `advance`.
#[derive(Clone, Debug)]
pub struct SharedAtmosphereMoistureStep {
    pub interval_s: [f64; 2],
    pub relative_humidity: [f64; 2],
    pub vapor_mass_kg: [f64; 2],
    pub water_to_solids_kg: f64,
    pub water_balance_residual_kg: f64,
    pub total_mass_residual_kg: f64,
    pub maximum_carrier_mass_residual_kg: f64,
    /// Gas m*u change at prescribed T; not the thermostat's total heat duty.
    pub gas_internal_energy_change_j: f64,
    pub iterations: usize,
}

/// Original gas and every moisture field accepted together. Solids are in
/// increasing stable-ID order, independent of caller enumeration order.
#[derive(Clone, Debug)]
pub struct SharedAtmosphereMoistureCandidate {
    atmosphere: FiniteMoistureAtmosphere,
    solids: Vec<SharedMoistureSolid>,
    step: SharedAtmosphereMoistureStep,
}
impl SharedAtmosphereMoistureCandidate {
    pub const fn atmosphere(&self) -> &FiniteMoistureAtmosphere {
        &self.atmosphere
    }
    pub fn solids(&self) -> &[SharedMoistureSolid] {
        &self.solids
    }
    pub const fn step(&self) -> &SharedAtmosphereMoistureStep {
        &self.step
    }
    pub fn into_parts(
        self,
    ) -> (
        FiniteMoistureAtmosphere,
        Vec<SharedMoistureSolid>,
        SharedAtmosphereMoistureStep,
    ) {
        (self.atmosphere, self.solids, self.step)
    }
}

impl FiniteMoistureAtmosphere {
    /// Bind the original gas and explicitly identified water species. At least
    /// one noncondensing carrier is required so exact drying leaves a gas state.
    /// The existing Buck liquid-water saturation law supplies its stated domain.
    pub fn new(gas: ClosedGasParcel) -> Result<Self, AtmosphereMoistureError> {
        let catalog = gas.model().clone();
        let id = catalog.water_species().ok_or_else(|| {
            invalid("finite humidity needs explicitly identified vapor source data")
        })?;
        let water = catalog
            .thermodynamics()
            .components()
            .iter()
            .position(|s| s.species() == id)
            .expect("water in catalog");
        if catalog.thermodynamics().components().len() < 2 {
            return Err(invalid(
                "finite atmosphere requires a noncondensing carrier",
            ));
        }
        let temperature = gas.state().state().temperature;
        let saturation_mass_kg = fs_material::gas::saturation_pressure_water_pa(temperature)?
            * gas.volume_m3()
            * catalog.thermodynamics().components()[water]
                .molar_mass()
                .value()
            / (R * temperature);
        if !positive(saturation_mass_kg) {
            return Err(invalid("finite saturation inventory is unrepresentable"));
        }
        let result = Self {
            gas,
            catalog,
            water,
            saturation_mass_kg,
        };
        let rh = result.relative_humidity()?;
        if !(0.0..=1.0).contains(&rh) {
            return Err(invalid("supersaturation requires condensation physics"));
        }
        Ok(result)
    }
    pub const fn gas(&self) -> &ClosedGasParcel {
        &self.gas
    }
    pub fn relative_humidity(&self) -> Result<f64, AtmosphereMoistureError> {
        Ok(self.inventory()?[self.water] / self.saturation_mass_kg)
    }
    /// Actual component masses in the retained catalog's canonical order.
    fn inventory(&self) -> Result<Vec<f64>, AtmosphereMoistureError> {
        let active = self.gas.model().thermodynamics();
        let molar_mass = self
            .gas
            .state()
            .thermodynamics()
            .properties()
            .molar_mass()
            .value();
        let amount = self.gas.mass_kg() / molar_mass;
        if !positive(amount) {
            return Err(invalid("gas amount is unrepresentable"));
        }
        let mut result = vec![0.0; self.catalog.thermodynamics().components().len()];
        for (component, fraction) in active
            .components()
            .iter()
            .zip(active.mole_composition().fractions())
        {
            let index = self
                .catalog
                .thermodynamics()
                .components()
                .iter()
                .position(|s| s.species() == component.species())
                .ok_or_else(|| invalid("gas species is outside the retained source catalog"))?;
            let mass = amount * fraction * component.molar_mass().value();
            if !positive(mass) {
                return Err(invalid("positive species inventory was lost"));
            }
            result[index] = mass;
        }
        Ok(result)
    }

    fn endpoint(&self, water_mass: f64, old: &[f64]) -> Result<Self, AtmosphereMoistureError> {
        let mut inventory = old.to_vec();
        inventory[self.water] = water_mass;
        let amounts: Vec<f64> = inventory
            .iter()
            .zip(self.catalog.thermodynamics().components())
            .map(|(mass, source)| mass / source.molar_mass().value())
            .collect();
        let total: f64 = amounts.iter().sum();
        if !positive(total)
            || amounts
                .iter()
                .zip(&inventory)
                .any(|(n, m)| *m > 0.0 && !positive(*n))
        {
            return Err(invalid(
                "accepted vapor/carrier amounts are unrepresentable",
            ));
        }
        let mut fractions: Vec<f64> = amounts.iter().map(|n| n / total).collect();
        let sum: f64 = fractions.iter().sum();
        if !sum.is_finite() || (sum - 1.0).abs() > 256.0 * f64::EPSILON {
            return Err(invalid("composition normalization exceeds roundoff"));
        }
        if sum != 1.0 {
            // Only reconcile roundoff; never create a trace species or floor one.
            let last = amounts
                .iter()
                .rposition(|n| *n > 0.0)
                .expect("positive total");
            fractions[last] = 1.0 - fractions[..last].iter().sum::<f64>();
        }
        if fractions
            .iter()
            .zip(&amounts)
            .any(|(x, n)| *n > 0.0 && !positive(*x))
        {
            return Err(invalid("a positive atmospheric species fraction was lost"));
        }
        let composition = Composition::new(CompositionBasis::MoleFraction, fractions)
            .map_err(|_| invalid("invalid accepted atmosphere composition"))?;
        let model = self.catalog.with_composition(composition)?;
        let pressure = total * R * self.gas.state().state().temperature / self.gas.volume_m3();
        if !positive(pressure) {
            return Err(invalid("accepted gas pressure is unrepresentable"));
        }
        let point = self
            .gas
            .point()
            .clone()
            .with("pressure", pressure)
            .map_err(|_| invalid("invalid accepted pressure coordinate"))?;
        let gas = ClosedGasParcel::new(model, point, self.gas.volume_m3())?;
        Ok(Self {
            gas,
            ..self.clone()
        })
    }

    /// Solve one interval using endpoint RH in the existing backward-Euler slab
    /// boundary. Optional finite surface exchange is resolved at the same trial
    /// RH; nonlinear interior diffusivity keeps the original accepted-start
    /// lagging. No accepted state is changed during trials or endpoint checks.
    pub fn advance(
        &self,
        slab: &MoistureSlab,
        dt_s: f64,
        faces: &[Option<DryBasisSorption>; 2],
        transport: MoistureStepConfig,
        controls: AtmosphereMoistureSolve,
        cx: &Cx<'_>,
    ) -> Result<AtmosphereMoistureCandidate, AtmosphereMoistureError> {
        let participant = MoistureParticipant {
            id: 0,
            slab,
            faces,
            transport,
        };
        let candidate =
            self.advance_shared(&[participant], dt_s, slab.moisture().len(), controls, cx)?;
        let (atmosphere, mut solids, shared) = candidate.into_parts();
        let solid = solids.pop().expect("one admitted participant");
        let step = AtmosphereMoistureStep {
            sorption: solid.transport.surface_sorption.clone(),
            transport: solid.transport,
            relative_humidity: shared.relative_humidity,
            vapor_mass_kg: shared.vapor_mass_kg,
            water_to_solid_kg: shared.water_to_solids_kg,
            water_balance_residual_kg: shared.water_balance_residual_kg,
            maximum_carrier_mass_residual_kg: shared.maximum_carrier_mass_residual_kg,
            gas_internal_energy_change_j: shared.gas_internal_energy_change_j,
            iterations: shared.iterations,
        };
        Ok(AtmosphereMoistureCandidate {
            atmosphere,
            moisture: solid.moisture,
            step,
        })
    }

    /// Advance heterogeneous specimens against ONE finite vapor inventory.
    ///
    /// Every scalar trial uses the same endpoint RH for all exposed faces and
    /// restarts each native diffusion owner from its accepted field. Source
    /// sorption and surface-transfer laws remain specimen-local; the atmosphere
    /// is debited only once, after all fields and the gas endpoint pass. This
    /// permits simultaneous adsorption and desorption without ordering bias or
    /// treating the same vapor as a separate reservoir for every specimen.
    ///
    /// Stable IDs must be unique and accepted clocks and temperatures must agree.
    /// The total cell budget includes sealed specimens, whose internal diffusion
    /// still advances. Native iteration budgets apply to each scalar trial.
    /// The existing isothermal, fixed-volume, dilute bound-water assumptions
    /// remain: this does not add sorption heat, gas flow, or condensation.
    pub fn advance_shared(
        &self,
        participants: &[MoistureParticipant<'_>],
        dt_s: f64,
        maximum_total_cells: usize,
        controls: AtmosphereMoistureSolve,
        cx: &Cx<'_>,
    ) -> Result<SharedAtmosphereMoistureCandidate, AtmosphereMoistureError> {
        poll(cx)?;
        if participants.is_empty()
            || maximum_total_cells == 0
            || maximum_total_cells > super::MAX_MOISTURE_CELLS
            || participants.len() > maximum_total_cells
        {
            return Err(invalid(
                "a nonempty bounded collection of moisture specimens is required",
            ));
        }
        let mut ordered: Vec<_> = participants.iter().collect();
        ordered.sort_unstable_by_key(|body| body.id);
        if ordered.windows(2).any(|pair| pair[0].id == pair[1].id) {
            return Err(invalid("shared-atmosphere specimen IDs must be unique"));
        }
        let time_s = ordered[0].slab.time_s();
        let temperature = self.gas.state().state().temperature;
        let mut cells = 0usize;
        for body in &ordered {
            poll(cx)?;
            cells = cells
                .checked_add(body.slab.moisture().len())
                .filter(|count| *count <= maximum_total_cells)
                .ok_or_else(|| invalid("shared-atmosphere total cell budget exhausted"))?;
            if body.slab.time_s().to_bits() != time_s.to_bits()
                || body.slab.reference().query().axes()["temperature"].to_bits()
                    != temperature.to_bits()
            {
                return Err(invalid(
                    "shared moisture specimens need the same accepted clock and gas temperature",
                ));
            }
        }
        let old_inventory = self.inventory()?;
        let initial_rh = old_inventory[self.water] / self.saturation_mass_kg;
        let [mut lower, mut upper] = controls.relative_humidity_bracket;
        if !lower.is_finite()
            || !upper.is_finite()
            || lower < 0.0
            || upper > 1.0
            || lower >= upper
            || initial_rh < lower
            || initial_rh > upper
            || !positive(controls.mass_tolerance_kg)
            || controls.maximum_iterations == 0
            || controls.maximum_iterations > 256
        {
            return Err(invalid(
                "finite mass tolerance, bounded iterations and a supported RH bracket containing the start are required",
            ));
        }

        let evaluate = |rh| -> Result<_, AtmosphereMoistureError> {
            poll(cx)?;
            let mut solids = Vec::with_capacity(ordered.len());
            // Balance increments: large pre-existing bound-water inventories
            // must not provide an artificial tolerance scale for vapor exchange.
            let mut water_change = 0.0;
            for body in &ordered {
                poll(cx)?;
                let mut boundary = std::array::from_fn(|_| ResolvedMoistureFace::default());
                for side in 0..2 {
                    if let Some(law) = &body.faces[side] {
                        boundary[side] = law.resolve(temperature, rh)?;
                    }
                }
                let (moisture, transport) =
                    body.slab
                        .advance_resolved(dt_s, boundary, body.transport, cx)?;
                water_change += moisture.water_mass_kg() - body.slab.water_mass_kg();
                solids.push(SharedMoistureSolid {
                    id: body.id,
                    moisture,
                    transport,
                });
            }
            let residual = self.saturation_mass_kg * (rh - initial_rh) + water_change;
            if !residual.is_finite() {
                return Err(invalid("coupled water residual is unrepresentable"));
            }
            Ok((residual, solids, water_change))
        };
        let mut iterations = 0;
        let start = evaluate(initial_rh)?;
        let (rh, accepted) = if start.0.abs() <= controls.mass_tolerance_kg {
            (initial_rh, start)
        } else {
            if ordered
                .iter()
                .all(|body| body.faces.iter().all(Option::is_none))
            {
                return Err(AtmosphereMoistureError::MassBalance {
                    residual_kg: start.0,
                    tolerance_kg: controls.mass_tolerance_kg,
                });
            }
            let mut low = evaluate(lower)?;
            let high = evaluate(upper)?;
            if low.0.abs() <= controls.mass_tolerance_kg {
                (lower, low)
            } else if high.0.abs() <= controls.mass_tolerance_kg {
                (upper, high)
            } else {
                if low.0.is_sign_positive() == high.0.is_sign_positive() {
                    return Err(AtmosphereMoistureError::Unbracketed {
                        residual_kg: [low.0, high.0],
                    });
                }
                loop {
                    poll(cx)?;
                    let middle = lower + 0.5 * (upper - lower);
                    let trial = evaluate(middle)?;
                    iterations += 1;
                    if trial.0.abs() <= controls.mass_tolerance_kg {
                        break (middle, trial);
                    }
                    if iterations == controls.maximum_iterations
                        || middle == lower
                        || middle == upper
                    {
                        return Err(AtmosphereMoistureError::NotConverged {
                            residual_kg: trial.0,
                            iterations,
                        });
                    }
                    if trial.0.is_sign_positive() == low.0.is_sign_positive() {
                        lower = middle;
                        low = trial;
                    } else {
                        upper = middle;
                    }
                }
            }
        };
        let (_, solids, water_to_solids_kg) = accepted;
        let atmosphere = if rh.to_bits() == initial_rh.to_bits() {
            self.clone()
        } else {
            self.endpoint(self.saturation_mass_kg * rh, &old_inventory)?
        };
        let new_inventory = atmosphere.inventory()?;
        let residual = new_inventory[self.water] - old_inventory[self.water] + water_to_solids_kg;
        let carrier_residual = new_inventory
            .iter()
            .zip(&old_inventory)
            .enumerate()
            .filter(|(i, _)| *i != self.water)
            .map(|(_, (a, b))| (a - b).abs())
            .fold(0.0, f64::max);
        let mass_residual = atmosphere.gas.mass_kg() - self.gas.mass_kg() + water_to_solids_kg;
        let lost_species = new_inventory
            .iter()
            .zip(&old_inventory)
            .enumerate()
            .any(|(i, (a, b))| i != self.water && *b > 0.0 && !positive(*a));
        for error in [residual, carrier_residual, mass_residual] {
            if !error.is_finite() || error.abs() > controls.mass_tolerance_kg {
                return Err(AtmosphereMoistureError::MassBalance {
                    residual_kg: error,
                    tolerance_kg: controls.mass_tolerance_kg,
                });
            }
        }
        if lost_species {
            return Err(invalid("accepted exchange lost a noncondensing species"));
        }
        let gas_internal_energy_change_j = atmosphere.gas.mass_kg()
            * atmosphere.gas.state().specific_internal_energy_j_kg()
            - self.gas.mass_kg() * self.gas.state().specific_internal_energy_j_kg();
        if !gas_internal_energy_change_j.is_finite() {
            return Err(invalid("gas storage change is unrepresentable"));
        }
        let step = SharedAtmosphereMoistureStep {
            interval_s: [time_s, solids[0].moisture.time_s()],
            relative_humidity: [initial_rh, atmosphere.relative_humidity()?],
            vapor_mass_kg: [old_inventory[self.water], new_inventory[self.water]],
            water_to_solids_kg,
            water_balance_residual_kg: residual,
            total_mass_residual_kg: mass_residual,
            maximum_carrier_mass_residual_kg: carrier_residual,
            gas_internal_energy_change_j,
            iterations,
        };
        poll(cx)?;
        Ok(SharedAtmosphereMoistureCandidate {
            atmosphere,
            solids,
            step,
        })
    }
}
