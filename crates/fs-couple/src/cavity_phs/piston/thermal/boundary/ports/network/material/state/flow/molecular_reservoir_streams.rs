//! Directional stream accuracy over the same coarse and two fine intervals.
//!
//! Bath sensible enthalpy uses ONE fixed native chart for each bath/species.
//! It is additional to separately controlled species and reduced-heat streams;
//! it is not an accuracy claim about an absolute formation-energy datum. Native
//! per-step relative-energy values use moving gas references and are never
//! integrated here as a physical boundary-energy account.

use crate::cavity_phs::piston::thermal::boundary::ports::network::material::state::flow::{
    MaterialFlowError, diffusion::coupled::scaled_error,
};
use fs_exec::Cx;
use fs_gas::material::transport::parcel::reservoir::{
    ReservoirMolecularLinkStep, ReservoirMolecularStep,
};
use fs_material::thermochemical_gas::energy::parcel::transfer::{
    snapshot::MAX_SNAPSHOT_TRANSFERS, species::network::boundary::SpeciesNetworkBoundaryFlux,
};
use fs_thermochem::SpeciesId;
use std::collections::BTreeMap;

/// Positive absolute LOCAL error allowances for separately exposed streams.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MolecularReservoirStreamScales {
    /// Each face/species/direction, including steady throughflow [kg].
    pub species_mass_kg: f64,
    /// Each directional additional-heat or fixed-chart sensible-energy stream [J].
    pub heat_j: f64,
}
impl MolecularReservoirStreamScales {
    /// Admit accuracy scales before spending any complete physical trial.
    pub fn validate(self) -> Result<(), MaterialFlowError> {
        if [self.species_mass_kg, self.heat_j]
            .into_iter()
            .all(|value| value.is_finite() && value > 0.0)
        {
            Ok(())
        } else {
            Err(MaterialFlowError::Invalid(
                "positive finite molecular species and heat stream scales required",
            ))
        }
    }
}

/// Dimensionless directional coarse/fine discrepancies; each must be <= 1.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct MolecularReservoirStreamEstimate {
    pub species_mass: f64,
    pub heat: f64,
    /// Bath sensible enthalpy in a fixed native species chart, plus reduced heat.
    /// The independent heat/species errors prevent cancellation from hiding either.
    pub thermal_energy: f64,
}
impl MolecularReservoirStreamEstimate {
    pub fn maximum(self) -> f64 {
        self.species_mass.max(self.heat).max(self.thermal_energy)
    }
    pub fn accumulate(&mut self, other: Self) {
        self.species_mass = self.species_mass.max(other.species_mass);
        self.heat = self.heat.max(other.heat);
        self.thermal_energy = self.thermal_energy.max(other.thermal_energy);
    }
}

type Error = MaterialFlowError;
fn check(cx: &Cx<'_>) -> Result<(), Error> {
    cx.checkpoint().map_err(|_| Error::Cancelled)
}
fn finite(value: f64) -> Result<f64, Error> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(Error::Invalid("unrepresentable molecular stream account"))
    }
}

/// Including owners keeps an ID reassignment from masquerading as a zero error.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
enum Face {
    Internal(u64, [u64; 2]),
    Reservoir(u64, u64),
}
type SpeciesKey<'a> = (Face, &'a SpeciesId);

struct View<'a> {
    step: &'a ReservoirMolecularStep,
    faces: BTreeMap<u64, Face>,
    baths: BTreeMap<Face, &'a ReservoirMolecularLinkStep>,
    boundaries: BTreeMap<SpeciesKey<'a>, &'a SpeciesNetworkBoundaryFlux>,
}
impl<'a> View<'a> {
    fn new(step: &'a ReservoirMolecularStep, cx: &Cx<'_>) -> Result<Self, Error> {
        check(cx)?;
        let count = step.species.len();
        if count == 0
            || count > MAX_SNAPSHOT_TRANSFERS
            || !step.species.windows(2).all(|pair| pair[0] < pair[1])
            || !step.duration_s.is_finite()
            || step.duration_s <= 0.0
            || step.exchange.boundaries.len() > MAX_SNAPSHOT_TRANSFERS
            || step
                .links
                .len()
                .checked_mul(count + 2)
                .and_then(|bath| {
                    step.internal_links
                        .len()
                        .checked_mul(count + 1)
                        .and_then(|internal| bath.checked_add(internal))
                })
                .is_none_or(|work| work > MAX_SNAPSHOT_TRANSFERS)
        {
            return Err(Error::Invalid(
                "bounded canonical molecular stream graph required",
            ));
        }
        let mut view = Self {
            step,
            faces: BTreeMap::new(),
            baths: BTreeMap::new(),
            boundaries: BTreeMap::new(),
        };
        for link in &step.internal_links {
            check(cx)?;
            if link.parcels[0] == link.parcels[1]
                || link.mass_left_to_right_kg.len() != count
                || !link
                    .mass_left_to_right_kg
                    .iter()
                    .copied()
                    .all(f64::is_finite)
                || !link.heat_correction_left_to_right_j.is_finite()
                || view
                    .faces
                    .insert(link.id, Face::Internal(link.id, link.parcels))
                    .is_some()
            {
                return Err(Error::Invalid(
                    "invalid internal molecular stream identity or axis",
                ));
            }
        }
        for link in &step.links {
            check(cx)?;
            let face = Face::Reservoir(link.id, link.parcel);
            if link.mass_into_gas_kg.len() != count
                || !link.mass_into_gas_kg.iter().copied().all(f64::is_finite)
                || !link.reduced_heat_into_gas_j.is_finite()
                || view.faces.insert(link.id, face).is_some()
            {
                return Err(Error::Invalid("invalid reservoir stream identity or axis"));
            }
            view.baths.insert(face, link);
        }
        for row in &step.exchange.boundaries {
            check(cx)?;
            let face = Face::Reservoir(row.id, row.parcel);
            let link = view
                .baths
                .get(&face)
                .ok_or(Error::Invalid("native stream has no declared held bath"))?;
            let index = step
                .species
                .binary_search(&row.transfer.species)
                .map_err(|_| Error::Invalid("native bath species is outside the original axis"))?;
            let mass = link.mass_into_gas_kg[index];
            let transfer = &row.transfer;
            if transfer.net_mass_into_gas_kg != mass
                || transfer.incoming_mass_kg != mass.max(0.0)
                || transfer.outgoing_mass_kg != (-mass).max(0.0)
                || !transfer.relative_energy_into_gas_j.is_finite()
                || (mass == 0.0 && transfer.relative_energy_into_gas_j != 0.0)
                || view
                    .boundaries
                    .insert((face, &transfer.species), row)
                    .is_some()
            {
                return Err(Error::Invalid(
                    "native directional bath receipt does not match its face",
                ));
            }
        }
        for (&face, link) in &view.baths {
            for (species, &mass) in step.species.iter().zip(&link.mass_into_gas_kg) {
                check(cx)?;
                if mass != 0.0 && !view.boundaries.contains_key(&(face, species)) {
                    return Err(Error::Invalid(
                        "nonzero bath stream omitted its native caloric chart",
                    ));
                }
            }
        }
        Ok(view)
    }

    fn admit_same_baths(&self, other: &Self, cx: &Cx<'_>) -> Result<(), Error> {
        if self.step.species != other.step.species || self.faces != other.faces {
            return Err(Error::Invalid(
                "molecular comparison changed its species or exposed graph",
            ));
        }
        for (face, link) in &self.baths {
            check(cx)?;
            if other
                .baths
                .get(face)
                .is_none_or(|other| link.reservoir != other.reservoir)
            {
                return Err(Error::Invalid(
                    "molecular comparison changed a held bath source or state",
                ));
            }
        }
        Ok(())
    }
}

#[derive(Default)]
struct Sum {
    value: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, value: f64) -> Result<(), Error> {
        let value = finite(value)?;
        let next = finite(self.value + value)?;
        self.correction = finite(
            self.correction
                + if self.value.abs() >= value.abs() {
                    (self.value - next) + value
                } else {
                    (value - next) + self.value
                },
        )?;
        self.value = next;
        Ok(())
    }
    fn total(&self) -> Result<f64, Error> {
        finite(self.value + self.correction)
    }
}
#[derive(Default)]
struct Directional([Sum; 2]);
impl Directional {
    fn add(&mut self, signed: f64) -> Result<(), Error> {
        self.0[usize::from(signed < 0.0)].add(signed.abs())
    }
}
#[derive(Default)]
struct Streams<'a> {
    species: BTreeMap<SpeciesKey<'a>, Directional>,
    heat: BTreeMap<Face, Directional>,
    thermal: BTreeMap<Face, Directional>,
}
impl<'a> Streams<'a> {
    fn record(
        &mut self,
        view: &View<'a>,
        charts: &BTreeMap<SpeciesKey<'a>, f64>,
        cx: &Cx<'_>,
    ) -> Result<(), Error> {
        for link in &view.step.internal_links {
            check(cx)?;
            let face = Face::Internal(link.id, link.parcels);
            for (species, &mass) in view.step.species.iter().zip(&link.mass_left_to_right_kg) {
                check(cx)?;
                self.species.entry((face, species)).or_default().add(mass)?;
            }
            self.heat
                .entry(face)
                .or_default()
                .add(link.heat_correction_left_to_right_j)?;
        }
        for (&face, link) in &view.baths {
            check(cx)?;
            let thermal = self.thermal.entry(face).or_default();
            for (species, &mass) in view.step.species.iter().zip(&link.mass_into_gas_kg) {
                check(cx)?;
                self.species.entry((face, species)).or_default().add(mass)?;
                if mass != 0.0 {
                    let coefficient = charts.get(&(face, species)).ok_or(Error::Invalid(
                        "nonzero bath stream has no fixed native chart",
                    ))?;
                    thermal.add(finite(mass * coefficient)?)?;
                }
            }
            self.heat
                .entry(face)
                .or_default()
                .add(link.reduced_heat_into_gas_j)?;
            thermal.add(link.reduced_heat_into_gas_j)?;
        }
        Ok(())
    }
}

fn error<K: Ord>(
    coarse: &BTreeMap<K, Directional>,
    fine: &BTreeMap<K, Directional>,
    absolute: f64,
    relative: f64,
    factor: f64,
    cx: &Cx<'_>,
) -> Result<f64, Error> {
    if !coarse.keys().eq(fine.keys()) {
        return Err(Error::Invalid(
            "molecular stream comparison omitted a physical coordinate",
        ));
    }
    let mut maximum = 0.0_f64;
    for (coarse, fine) in coarse.values().zip(fine.values()) {
        check(cx)?;
        for direction in 0..2 {
            let ratio = scaled_error(
                coarse.0[direction].total()?,
                fine.0[direction].total()?,
                absolute,
                relative,
            )
            .ok_or(Error::Invalid("invalid molecular stream error scale"))?;
            // Infinity must reject; clamping or treating it as malformed zero
            // would accept an arbitrarily inaccurate physical stream.
            maximum = maximum.max(ratio * factor);
        }
    }
    Ok(maximum)
}

/// Compare one complete interval with its two accepted-order fine candidates.
///
/// Each held bath/species receives one native sensible-enthalpy coefficient
/// from the first nonzero receipt in coarse/first/second order. Applying that
/// SAME coefficient to all three mass accounts puts them in a common chart.
/// Dividing an already reference-cancelled native increment avoids subtracting
/// absolute enthalpies; no caloric source is reconstructed or reevaluated.
/// If all three streams vanish, that species needs no energy-chart coefficient.
/// Species, additional heat and sensible energy retain separate gross signs.
#[allow(clippy::too_many_arguments)]
pub(crate) fn compare(
    coarse: &ReservoirMolecularStep,
    first: &ReservoirMolecularStep,
    second: &ReservoirMolecularStep,
    scales: MolecularReservoirStreamScales,
    relative_tolerance: f64,
    fine_error_factor: f64,
    cx: &Cx<'_>,
) -> Result<MolecularReservoirStreamEstimate, Error> {
    check(cx)?;
    scales.validate()?;
    if !relative_tolerance.is_finite()
        || !(0.0..1.0).contains(&relative_tolerance)
        || !fine_error_factor.is_finite()
        || fine_error_factor < 1.0
    {
        return Err(Error::Invalid(
            "finite molecular relative accuracy and fine-step factor required",
        ));
    }
    let coarse = View::new(coarse, cx)?;
    let first = View::new(first, cx)?;
    let second = View::new(second, cx)?;
    coarse.admit_same_baths(&first, cx)?;
    coarse.admit_same_baths(&second, cx)?;
    let mut charts = BTreeMap::new();
    for view in [&coarse, &first, &second] {
        for (&key, row) in &view.boundaries {
            check(cx)?;
            let mass = row.transfer.net_mass_into_gas_kg;
            if mass != 0.0 {
                if let std::collections::btree_map::Entry::Vacant(entry) = charts.entry(key) {
                    entry.insert(finite(row.transfer.relative_energy_into_gas_j / mass)?);
                }
            }
        }
    }
    let mut a = Streams::default();
    let mut b = Streams::default();
    a.record(&coarse, &charts, cx)?;
    b.record(&first, &charts, cx)?;
    b.record(&second, &charts, cx)?;
    let estimate = MolecularReservoirStreamEstimate {
        species_mass: error(
            &a.species,
            &b.species,
            scales.species_mass_kg,
            relative_tolerance,
            fine_error_factor,
            cx,
        )?,
        heat: error(
            &a.heat,
            &b.heat,
            scales.heat_j,
            relative_tolerance,
            fine_error_factor,
            cx,
        )?,
        thermal_energy: error(
            &a.thermal,
            &b.thermal,
            scales.heat_j,
            relative_tolerance,
            fine_error_factor,
            cx,
        )?,
    };
    check(cx)?;
    Ok(estimate)
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_evidence::ValidityDomain;
    use fs_exec::CancelGate;
    use fs_gas::material::transport::parcel::ParcelMolecularLinkStep;
    use fs_gas::multicomponent::dissipation::TransportFaceReport;
    use fs_matdb::*;
    use fs_material::thermochemical_gas::energy::parcel::transfer::species::network::{
        SpeciesNetworkStep, boundary::BoundarySpeciesTransfer, transport::GasReservoirWindowAudit,
    };
    use fs_material::thermochemical_gas::{
        ResolvedThermochemicalGasState, SpeciesTransportMaterial, ThermochemicalGasModel,
    };
    use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
    use fs_qty::{DynViscosity, MolarMass, Pressure, QuantitySpec};
    use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;
    use fs_thermochem::transport::kinetic::ThermalConductivityQuantityV1;
    use fs_thermochem::*;

    fn with_cx(run: impl FnOnce(&Cx<'_>, &CancelGate)) {
        let gate = CancelGate::new();
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(
                &gate,
                arena,
                fs_exec::StreamKey {
                    seed: 31,
                    kernel_id: 71,
                    tile: 0,
                    iteration: 0,
                },
                fs_exec::Budget::INFINITE,
                fs_exec::ExecMode::Deterministic,
            );
            run(&cx, &gate);
        });
    }
    fn close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() <= 2e-14 * expected.abs().max(1.0),
            "{actual:e} != {expected:e}"
        );
    }

    // Real source receipts; synthetic integrated streams below exercise only
    // the estimator, independently of any transport or caloric time trajectory.
    fn held_state(formation_j_kg: f64) -> ResolvedThermochemicalGasState {
        let make = |name: &str| {
            let species = SpeciesId::new(name).unwrap();
            let mut coefficients = [0.0; 9];
            coefficients[2] = 3.5;
            coefficients[7] = formation_j_kg * 0.028 / UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K;
            let mut parameters: BTreeMap<_, _> = required_nasa9_parameter_dims()
                .into_iter()
                .take(9)
                .zip(coefficients)
                .map(|((name, dims), value)| (name.to_string(), LawParameter { value, dims }))
                .collect();
            parameters.insert(
                "reference_pressure".into(),
                LawParameter {
                    value: 1e5,
                    dims: Pressure::DIMS,
                },
            );
            let provenance = Provenance {
                source: "synthetic molecular stream error fixture".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            };
            let region = Nasa9RegionV1::from_card(ConstitutiveModelCard {
                law: LawId(NASA9_LAW_ID_V1.into()),
                law_version: NASA9_LAW_VERSION_V1,
                parameters,
                state_schema_version: NASA9_STATE_SCHEMA_VERSION_V1,
                initial_state: InitialStatePolicy::ZeroInternalState,
                validity: ValidityDomain::unconstrained().with("T", 200.0, 2000.0),
                sources: vec![],
                provenance: provenance.clone(),
            })
            .unwrap();
            let caloric = Nasa9StandardStateModelV1::new(
                species.clone(),
                MolarMass::new(0.028),
                StandardStateConventionV1::new(
                    StandardStatePhaseV1::Gas,
                    ReferenceEquationOfStateV1::IdealGas,
                    Pressure::new(1e5),
                    ElementalReferenceIdV1::new("synthetic").unwrap(),
                )
                .unwrap(),
                vec![region],
            )
            .unwrap();
            let mut claims = ClaimSet::new();
            for (name, dims, value) in [
                ("dynamic_viscosity", DynViscosity::DIMS, 2e-5),
                (
                    "thermal_conductivity",
                    ThermalConductivityQuantityV1::DIMS,
                    0.03,
                ),
            ] {
                claims
                    .insert_claim(PropertyClaim {
                        key: PropertyKey::new(name, dims),
                        value: PropertyValue::Scalar { value, dims },
                        validity: ValidityDomain::unconstrained(),
                        uncertainty: UncertaintyModel::Unstated,
                        interpolation: InterpolationPolicy::ConstantWithinValidity,
                        observations: vec![],
                        provenance: provenance.clone(),
                    })
                    .unwrap();
            }
            let material = MaterialCard::assemble(
                MaterialStateId {
                    chemistry: name.into(),
                    phase: "gas".into(),
                    process: "synthetic stream fixture".into(),
                    revision: 0,
                },
                claims,
                vec![],
            )
            .unwrap();
            (
                caloric,
                SpeciesTransportMaterial {
                    species,
                    material,
                    selection: fs_material::state_point::MaterialPropertySelection::SingleClaimOnly,
                },
            )
        };
        let (a, ta) = make("A");
        let (b, tb) = make("B");
        let mixture = FrozenIdealGasMixtureModelV1::new(
            vec![a, b],
            Composition::new(CompositionBasis::MoleFraction, vec![0.5, 0.5]).unwrap(),
        )
        .unwrap();
        let bank = ThermochemicalGasModel::new(mixture, vec![ta, tb], None).unwrap();
        let quantity = |kind| QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static));
        let point = QueryPoint::new()
            .with_quantity(
                "temperature",
                quantity(QuantityKind::AbsoluteTemperature),
                1600.0,
            )
            .unwrap()
            .with_quantity("pressure", quantity(QuantityKind::Pressure), 1e5)
            .unwrap();
        bank.evaluate(&point).unwrap()
    }

    // (bath ID, parcel, signed A mass, reduced heat, native relative h per kg).
    // B carries the opposite mass so each face retains the mass frame.
    fn step(
        bath: &ResolvedThermochemicalGasState,
        rows: &[(u64, u64, f64, f64, f64)],
        duration_s: f64,
    ) -> ReservoirMolecularStep {
        let species = vec![SpeciesId::new("A").unwrap(), SpeciesId::new("B").unwrap()];
        let mut boundaries = Vec::new();
        let links = rows
            .iter()
            .map(|&(id, parcel, mass, heat, relative_h)| {
                let masses = vec![mass, -mass];
                for (species, &mass) in species.iter().zip(&masses) {
                    if mass != 0.0 {
                        let relative_energy = mass * relative_h;
                        let original_h = bath.specific_enthalpy_j_kg() - relative_h;
                        boundaries.push(SpeciesNetworkBoundaryFlux {
                            id,
                            parcel,
                            transfer: BoundarySpeciesTransfer {
                                species: species.clone(),
                                incoming_mass_kg: mass.max(0.0),
                                outgoing_mass_kg: (-mass).max(0.0),
                                net_mass_into_gas_kg: mass,
                                relative_energy_into_gas_j: relative_energy,
                            },
                            original_partial_enthalpy_j_kg: original_h,
                            energy_into_gas_j: mass * original_h + relative_energy,
                        });
                    }
                }
                ReservoirMolecularLinkStep {
                    id,
                    parcel,
                    mass_into_gas_kg: masses,
                    reduced_heat_into_gas_j: heat,
                    reservoir: bath.clone(),
                    face: TransportFaceReport {
                        enthalpy_temperature_k: bath.state().temperature,
                        ..TransportFaceReport::default()
                    },
                }
            })
            .collect();
        ReservoirMolecularStep {
            duration_s,
            species,
            internal_links: vec![],
            links,
            exchange: SpeciesNetworkStep {
                fluxes: vec![],
                heat: vec![],
                external_heat: vec![],
                boundaries,
                parcels: vec![],
                total_mass_residual_kg: 0.0,
                maximum_species_mass_residual_kg: 0.0,
                total_energy_residual_j: 0.0,
            },
            audit: GasReservoirWindowAudit {
                total_mass_residual_kg: 0.0,
                maximum_species_mass_residual_kg: 0.0,
                energy_residual_j: 0.0,
                entropy_production_j_k: 0.0,
            },
        }
    }

    #[test]
    fn opposed_and_reversed_baths_cannot_hide_gross_stream_error_in_zero_net_state() {
        with_cx(|cx, _| {
            let bath = held_state(0.0);
            let mass = 2.0_f64.powi(-20);
            let scales = MolecularReservoirStreamScales {
                species_mass_kg: 0.25 * mass,
                heat_j: 0.25,
            };
            let coarse = step(
                &bath,
                &[
                    (7, 10, 2.0 * mass, 2.0, 0.0),
                    (8, 10, -2.0 * mass, -2.0, 0.0),
                ],
                2.0,
            );
            let first = step(
                &bath,
                &[(7, 10, mass, 1.0, 0.0), (8, 10, -mass, -1.0, 0.0)],
                1.0,
            );
            // Steady throughflow becomes stronger, or both baths reverse.
            // The receiving gas's aggregate species and heat change stays zero.
            for second_rate in [2.0, -1.0] {
                let second = step(
                    &bath,
                    &[
                        (7, 10, second_rate * mass, second_rate, 0.0),
                        (8, 10, -second_rate * mass, -second_rate, 0.0),
                    ],
                    1.0,
                );
                for row in [&coarse, &first, &second] {
                    assert_eq!(
                        row.links
                            .iter()
                            .map(|link| link.mass_into_gas_kg[0])
                            .sum::<f64>(),
                        0.0
                    );
                    assert_eq!(
                        row.links
                            .iter()
                            .map(|link| link.reduced_heat_into_gas_j)
                            .sum::<f64>(),
                        0.0
                    );
                }
                let estimate = compare(&coarse, &first, &second, scales, 0.0, 1.0, cx).unwrap();
                close(estimate.species_mass, 4.0);
                close(estimate.heat, 4.0);
                close(estimate.thermal_energy, 4.0);
                assert!(estimate.maximum() > 1.0);
            }
        });
    }

    #[test]
    fn changing_substep_caloric_references_do_not_become_a_boundary_energy_error() {
        with_cx(|cx, _| {
            let bath = held_state(1e24);
            let mass = 2.0_f64.powi(-20);
            let coefficient = 2.0_f64.powi(20);
            let coarse = step(&bath, &[(7, 10, 2.0 * mass, 4.0, coefficient)], 2.0);
            let first = step(&bath, &[(7, 10, mass, 2.0, coefficient)], 1.0);
            let second = step(&bath, &[(7, 10, mass, 2.0, 0.5 * coefficient)], 1.0);
            // The second native coefficient has a different original gas T.
            // Both paths still carry exactly the same species and held-bath Q.
            let estimate = compare(
                &coarse,
                &first,
                &second,
                MolecularReservoirStreamScales {
                    species_mass_kg: 1e-15,
                    heat_j: 1e-12,
                },
                0.0,
                1.0,
                cx,
            )
            .unwrap();
            assert_eq!(estimate, MolecularReservoirStreamEstimate::default());
        });
    }

    #[test]
    fn common_chart_thermal_error_survives_large_formation_datums_and_fine_time_scaling() {
        with_cx(|cx, _| {
            let mass = 2.0_f64.powi(-20);
            let coefficient = 2.0_f64.powi(20);
            let mut original = None;
            for datum in [0.0, 1e24] {
                let bath = held_state(datum);
                let coarse = step(&bath, &[(7, 10, 2.0 * mass, 0.25, coefficient)], 2.0);
                let first = step(&bath, &[(7, 10, mass, 0.125, coefficient)], 1.0);
                let second = step(&bath, &[(7, 10, 2.0 * mass, 0.125, coefficient / 2.0)], 1.0);
                let scales = MolecularReservoirStreamScales {
                    species_mass_kg: 1.0,
                    heat_j: 0.25,
                };
                let estimate = compare(&coarse, &first, &second, scales, 0.0, 1.25, cx).unwrap();
                assert!(estimate.species_mass < 1.0);
                assert_eq!(estimate.heat, 0.0);
                close(estimate.thermal_energy, 5.0);
                if let Some(original) = original {
                    assert_eq!(estimate, original);
                }
                original = Some(estimate);
                // Relative scale also uses only this common sensible chart.
                let relative = compare(&coarse, &first, &second, scales, 0.25, 1.0, cx).unwrap();
                close(relative.thermal_energy, 1.0);
            }
        });
    }

    #[test]
    fn newly_nonzero_bath_species_require_their_native_chart_and_internal_streams_stay_separate() {
        with_cx(|cx, _| {
            let bath = held_state(0.0);
            let mass = 2.0_f64.powi(-20);
            let scales = MolecularReservoirStreamScales {
                species_mass_kg: 1.0,
                heat_j: 0.25,
            };
            let coarse = step(&bath, &[(7, 10, 0.0, 0.0, 0.0)], 2.0);
            let first = step(&bath, &[(7, 10, 0.0, 0.0, 0.0)], 1.0);
            let mut second = step(&bath, &[(7, 10, mass, 0.0, 1.0 / mass)], 1.0);
            let estimate = compare(&coarse, &first, &second, scales, 0.0, 1.0, cx).unwrap();
            close(estimate.thermal_energy, 4.0);
            second.exchange.boundaries.clear();
            assert!(matches!(
                compare(&coarse, &first, &second, scales, 0.0, 1.0, cx),
                Err(Error::Invalid(_))
            ));
            let internal = |multiplier: f64, duration| {
                let mut report = step(&bath, &[], duration);
                report.internal_links.push(ParcelMolecularLinkStep {
                    id: 19,
                    parcels: [10, 20],
                    mass_left_to_right_kg: vec![multiplier * mass, -multiplier * mass],
                    heat_correction_left_to_right_j: multiplier,
                    face: TransportFaceReport::default(),
                });
                report
            };
            let internal = compare(
                &internal(2.0, 2.0),
                &internal(1.0, 1.0),
                &internal(2.0, 1.0),
                MolecularReservoirStreamScales {
                    species_mass_kg: 0.25 * mass,
                    heat_j: 0.25,
                },
                0.0,
                1.0,
                cx,
            )
            .unwrap();
            close(internal.species_mass, 4.0);
            close(internal.heat, 4.0);
            assert_eq!(internal.thermal_energy, 0.0);
        });
    }

    #[test]
    fn malformed_coordinates_changed_baths_and_nonfinite_accuracy_refuse() {
        with_cx(|cx, _| {
            let bath = held_state(0.0);
            let coarse = step(&bath, &[(7, 10, 2e-6, 2.0, 0.0)], 2.0);
            let first = step(&bath, &[(7, 10, 1e-6, 1.0, 0.0)], 1.0);
            let scales = MolecularReservoirStreamScales {
                species_mass_kg: 1e-6,
                heat_j: 1.0,
            };
            let mut wrong = first.clone();
            wrong.links[0].parcel = 20;
            wrong
                .exchange
                .boundaries
                .iter_mut()
                .for_each(|row| row.parcel = 20);
            assert!(matches!(
                compare(&coarse, &first, &wrong, scales, 0.0, 1.0, cx),
                Err(Error::Invalid(_))
            ));
            let mut wrong = first.clone();
            wrong.links[0].reservoir = held_state(1e24);
            assert!(matches!(
                compare(&coarse, &first, &wrong, scales, 0.0, 1.0, cx),
                Err(Error::Invalid(_))
            ));
            let mut wrong = first.clone();
            wrong.links[0].mass_into_gas_kg.pop();
            assert!(matches!(
                compare(&coarse, &first, &wrong, scales, 0.0, 1.0, cx),
                Err(Error::Invalid(_))
            ));
            for (relative, factor) in [
                (1.0, 1.0),
                (f64::NAN, 1.0),
                (0.0, 0.0),
                (0.0, f64::INFINITY),
            ] {
                assert!(matches!(
                    compare(&coarse, &first, &first, scales, relative, factor, cx),
                    Err(Error::Invalid(_))
                ));
            }
            assert!(
                MolecularReservoirStreamScales {
                    species_mass_kg: 0.0,
                    ..scales
                }
                .validate()
                .is_err()
            );
            assert!(
                MolecularReservoirStreamScales {
                    heat_j: f64::NAN,
                    ..scales
                }
                .validate()
                .is_err()
            );
        });
    }

    #[test]
    fn overflowing_error_ratios_reject_and_cancellation_is_not_masked() {
        with_cx(|cx, gate| {
            let bath = held_state(0.0);
            let coarse = step(&bath, &[(7, 10, 0.0, 2e300, 0.0)], 2.0);
            let first = step(&bath, &[(7, 10, 0.0, 1e300, 0.0)], 1.0);
            let second = step(&bath, &[(7, 10, 0.0, 2e300, 0.0)], 1.0);
            let scales = MolecularReservoirStreamScales {
                species_mass_kg: 1.0,
                heat_j: f64::MIN_POSITIVE,
            };
            let estimate = compare(&coarse, &first, &second, scales, 0.0, 1.0, cx).unwrap();
            assert!(estimate.heat.is_infinite());
            assert!(estimate.maximum() > 1.0);
            let mut accepted = MolecularReservoirStreamEstimate::default();
            accepted.accumulate(estimate);
            assert_eq!(accepted, estimate);
            gate.request();
            assert!(matches!(
                compare(&coarse, &first, &second, scales, 0.0, 1.0, cx),
                Err(Error::Cancelled)
            ));
        });
    }
}
