//! Source-consistent elementary ideal-gas kinetics, without time integration.
//!
//! Each caller-authored reaction uses integer elementary orders and
//! an explicitly declared Arrhenius, third-body or falloff forward law. Reverse
//! rates of reversible reactions follow the original NASA-9 Gibbs energies,
//! including the pressure-to-concentration standard-state conversion. No empirical mechanism,
//! ignition criterion or falloff law is inferred. An explicitly irreversible
//! source disables its reverse rate and exposes its signed entropy production;
//! such a source does not acquire a detailed-balance theorem. L3 consumers own
//! state evolution and numerical step acceptance.
//!
//! Concentrations are mol/m3. The rate scale is mol/(m3 s) at the reference
//! temperature and at each participating concentration c*=1 mol/m3:
//! q_f = q_ref (T/T_ref)^b exp[-Ea/R (1/T-1/T_ref)] product(c_i/c*)^alpha_i.
//! Consequently ln K_c = -Delta g^0/(RT) + Delta nu ln[p0/(RT c*)]. This
//! normalization keeps the rate-scale quantity dimensionally explicit even
//! when elementary reactions have different total orders. See Cantera 3.2,
//! <https://cantera.org/3.2/reference/kinetics/reaction-rates.html> and
//! <https://cantera.org/3.2/reference/kinetics/rate-constants.html>.

pub use crate::equilibrium::ReactionProgressRateQuantityV1;
use crate::equilibrium::composition::IdealGasEquilibriumV1;
use crate::{
    ChemistryError, ConservationCertificate, MolarEnergyQuantityV1, ReactionId, SpeciesId,
    StoichiometricMatrix, ThermochemErrorV1, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R,
    verify_conservation,
};
use fs_qty::{Pressure, Temperature};

mod caloric;
pub use caloric::admit_caloric_interval;
/// Complete pinned GRI-Mech 3.0 source using these same kinetic owners.
pub mod gri30;
/// Pinned hydrogen/oxygen source mechanism using these same kinetic owners.
pub mod hydrogen;
mod pressure;
pub use pressure::{
    FalloffBroadeningV1, PressureDependenceV1, ThirdBodyEfficienciesV1, TroeFalloffV1,
};

/// Maximum admitted elementary reactions in one mechanism, including the
/// 325-channel GRI-Mech source. Workspace allocation remains explicitly bounded.
pub const MAX_KINETIC_REACTIONS_V1: usize = 512;
/// Maximum sum of elementary exponents on either side of one reaction.
pub const MAX_ELEMENTARY_ORDER_V1: u32 = 32;
/// Fixed operation-tree version of the concentration-based kinetic closure.
pub const IDEAL_GAS_KINETICS_VERSION_V1: u32 = 3;

/// Invalid source declaration, input domain, arithmetic range or cancellation.
#[derive(Clone, Debug, PartialEq)]
pub enum KineticsErrorV1 {
    /// Invalid dimensions, identities, controls or physical input.
    Invalid(&'static str),
    /// Exact elemental/charge conservation refused the mechanism.
    Conservation(ChemistryError),
    /// One admitted reaction disagrees with source molecular masses.
    MolecularMassBalance {
        /// Canonical reaction whose declared masses disagree.
        reaction: ReactionId,
        /// Absolute net mass divided by the sum of absolute mass contributions.
        relative_residual: f64,
    },
    /// Original NASA-9 source refused the queried temperature.
    SpeciesEvaluation {
        /// Original species identity.
        species: SpeciesId,
        /// Exact original source refusal.
        source: ThermochemErrorV1,
    },
    /// A nonzero physical value or derivative is outside binary64 range.
    Arithmetic(&'static str),
    /// Bounded scratch reservation failed.
    AllocationRefused,
    /// Cancellation before a complete evaluation was available.
    Cancelled,
}
impl core::fmt::Display for KineticsErrorV1 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "ideal-gas kinetics: {self:?}")
    }
}
impl std::error::Error for KineticsErrorV1 {}

/// Modified Arrhenius scale referenced to one explicit temperature.
///
/// A conventional k(T)=A T^b exp(-Ea/RT) becomes
/// q_ref=k(T_ref)*(1 mol/m3)^forward_order. The numerical coherent-SI value
/// is k(T_ref), but its quantity is always a volumetric progress rate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ArrheniusRateV1 {
    reference_temperature: Temperature,
    reference_rate: ReactionProgressRateQuantityV1,
    temperature_exponent: f64,
    activation_energy: MolarEnergyQuantityV1,
}
impl ArrheniusRateV1 {
    /// Retain a finite, nonnegative authored rate. Zero disables both directions.
    /// Negative activation energies are permitted as explicit source data.
    pub fn new(
        reference_temperature: Temperature,
        reference_rate: ReactionProgressRateQuantityV1,
        temperature_exponent: f64,
        activation_energy: MolarEnergyQuantityV1,
    ) -> Result<Self, KineticsErrorV1> {
        if !positive(reference_temperature.value())
            || !reference_rate.value().is_finite()
            || reference_rate.value() < 0.0
            || !temperature_exponent.is_finite()
            || !activation_energy.value().is_finite()
        {
            return Err(KineticsErrorV1::Invalid("invalid Arrhenius reference law"));
        }
        Ok(Self {
            reference_temperature,
            reference_rate,
            temperature_exponent,
            activation_energy,
        })
    }
    /// Original rate scale at the declared reference temperature, mol/(m3 s).
    pub const fn reference_rate(&self) -> ReactionProgressRateQuantityV1 {
        self.reference_rate
    }
    /// Original reference temperature.
    pub const fn reference_temperature(&self) -> Temperature {
        self.reference_temperature
    }
    /// Original dimensionless temperature exponent.
    pub const fn temperature_exponent(&self) -> f64 {
        self.temperature_exponent
    }
    /// Original molar activation energy, J/mol.
    pub const fn activation_energy(&self) -> MolarEnergyQuantityV1 {
        self.activation_energy
    }

    fn log_rate(self, t: f64) -> Result<(f64, f64), KineticsErrorV1> {
        if self.reference_rate.value() == 0.0 {
            return Ok((f64::NEG_INFINITY, 0.0));
        }
        let reference = self.reference_temperature.value();
        let ea_r = self.activation_energy.value() / R;
        let log = fs_math::det::ln(self.reference_rate.value())
            + self.temperature_exponent * crate::entropy::log_ratio(t, reference)
            + ea_r * ((t - reference) / t) / reference;
        let derivative = (self.temperature_exponent + ea_r / t) / t;
        if !log.is_finite() || !derivative.is_finite() {
            return Err(KineticsErrorV1::Arithmetic(
                "Arrhenius logarithm or derivative",
            ));
        }
        Ok((log, derivative))
    }
}

/// Original caller-declared elementary reaction, reversible by default.
///
/// Separate reactant/product sides retain explicit catalysts appearing on both
/// sides. An optional authored pressure dependence multiplies both directions.
/// A source may explicitly disable the reverse direction. Fractional/global
/// empirical orders and ions are outside this thermodynamic closure.
#[derive(Clone, Debug, PartialEq)]
pub struct ElementaryReactionV1 {
    reaction: ReactionId,
    reactants: Vec<(SpeciesId, u32)>,
    products: Vec<(SpeciesId, u32)>,
    forward_rate: ArrheniusRateV1,
    pressure_dependence: PressureDependenceV1,
    reversible: bool,
}
impl ElementaryReactionV1 {
    /// Canonicalize both sides and admit positive bounded integer orders.
    pub fn new(
        reaction: ReactionId,
        mut reactants: Vec<(SpeciesId, u32)>,
        mut products: Vec<(SpeciesId, u32)>,
        forward_rate: ArrheniusRateV1,
    ) -> Result<Self, KineticsErrorV1> {
        for side in [&mut reactants, &mut products] {
            if side.is_empty() || side.len() > MAX_ELEMENTARY_ORDER_V1 as usize {
                return Err(KineticsErrorV1::Invalid("elementary side has invalid size"));
            }
            side.sort_by(|a, b| a.0.cmp(&b.0));
            if side.windows(2).any(|pair| pair[0].0 == pair[1].0)
                || side
                    .iter()
                    .any(|(_, order)| *order == 0 || *order > MAX_ELEMENTARY_ORDER_V1)
                || side.iter().map(|(_, order)| *order as u64).sum::<u64>()
                    > MAX_ELEMENTARY_ORDER_V1 as u64
            {
                return Err(KineticsErrorV1::Invalid(
                    "duplicate species or invalid elementary order",
                ));
            }
        }
        if reactants == products {
            return Err(KineticsErrorV1::Invalid(
                "reaction has no net stoichiometry",
            ));
        }
        Ok(Self {
            reaction,
            reactants,
            products,
            forward_rate,
            pressure_dependence: PressureDependenceV1::Independent,
            reversible: true,
        })
    }
    /// Attach an explicit collider/falloff law. For falloff, `forward_rate` is
    /// the high-pressure law and the attached law supplies the low-pressure
    /// scale. No collider is inferred from either stoichiometric side.
    pub fn with_pressure_dependence(mut self, pressure_dependence: PressureDependenceV1) -> Self {
        self.pressure_dependence = pressure_dependence;
        self
    }
    /// Retain an original irreversible source channel without inventing a
    /// reverse reaction. Its chemical entropy production can be negative:
    /// callers must retain their physical-domain and finite-step entropy gates.
    pub fn irreversible(mut self) -> Self {
        self.reversible = false;
        self
    }
    /// Whether the reverse rate is derived from the original standard states.
    pub const fn is_reversible(&self) -> bool {
        self.reversible
    }
    /// Original reaction identity.
    pub const fn reaction(&self) -> &ReactionId {
        &self.reaction
    }
    /// Canonical elementary reactant orders.
    pub fn reactants(&self) -> &[(SpeciesId, u32)] {
        &self.reactants
    }
    /// Canonical elementary product orders.
    pub fn products(&self) -> &[(SpeciesId, u32)] {
        &self.products
    }
    /// Original forward Arrhenius law; a retained reverse is thermodynamically derived.
    pub const fn forward_rate(&self) -> ArrheniusRateV1 {
        self.forward_rate
    }
    /// Exact caller-authored pressure law, including collider efficiencies.
    pub const fn pressure_dependence(&self) -> &PressureDependenceV1 {
        &self.pressure_dependence
    }
}

#[derive(Clone, Debug, PartialEq)]
struct BoundReaction {
    reactants: Vec<(usize, u32)>,
    products: Vec<(usize, u32)>,
    net: Vec<(usize, f64)>,
    delta_molecules: f64,
    /// Canonical ln(epsilon_i); -infinity means an exactly inactive collider.
    /// Empty for pressure-independent rates.
    collider_logs: Vec<f64>,
}

/// Immutable neutral ideal-gas mechanism sharing one original NASA-9 bank.
#[derive(Clone, Debug, PartialEq)]
pub struct IdealGasKineticsV1 {
    thermodynamics: IdealGasEquilibriumV1,
    mechanism_id: String,
    mechanism_version: u32,
    reactions: Vec<ElementaryReactionV1>,
    bound: Vec<BoundReaction>,
    stoichiometry: StoichiometricMatrix,
    conservation: ConservationCertificate,
    relative_mass_tolerance: f64,
}
impl IdealGasKineticsV1 {
    /// Validate exact A*N=0 and z^T*N=0 using the original fs-qty owner.
    /// Source molecular-mass disagreement is independently bounded by the
    /// caller's tolerance (0 through1e-6 relative), never repaired numerically.
    pub fn new(
        thermodynamics: IdealGasEquilibriumV1,
        mechanism_id: impl Into<String>,
        mechanism_version: u32,
        mut reactions: Vec<ElementaryReactionV1>,
        relative_mass_tolerance: f64,
    ) -> Result<Self, KineticsErrorV1> {
        let mechanism_id = mechanism_id.into();
        if mechanism_id.is_empty()
            || mechanism_id.len() > 128
            || !mechanism_id
                .bytes()
                .all(|x| x.is_ascii_alphanumeric() || b"-_.:/".contains(&x))
            || mechanism_version == 0
            || reactions.is_empty()
            || reactions.len() > MAX_KINETIC_REACTIONS_V1
            || !relative_mass_tolerance.is_finite()
            || !(0.0..=1e-6).contains(&relative_mass_tolerance)
        {
            return Err(KineticsErrorV1::Invalid(
                "invalid mechanism identity, dimensions or mass tolerance",
            ));
        }
        reactions.sort_by(|a, b| a.reaction.cmp(&b.reaction));
        if reactions.windows(2).any(|p| p[0].reaction == p[1].reaction) {
            return Err(KineticsErrorV1::Invalid("duplicate reaction identity"));
        }
        let species = thermodynamics.elemental_matrix().species();
        let mut rows = vec![vec![0_i128; reactions.len()]; species.len()];
        let mut bound = Vec::new();
        bound
            .try_reserve_exact(reactions.len())
            .map_err(|_| KineticsErrorV1::AllocationRefused)?;
        for (column, reaction) in reactions.iter().enumerate() {
            let bind = |side: &[(SpeciesId, u32)]| -> Result<Vec<(usize, u32)>, KineticsErrorV1> {
                side.iter()
                    .map(|(id, order)| {
                        species.binary_search(id).map(|i| (i, *order)).map_err(|_| {
                            KineticsErrorV1::Invalid(
                                "reaction species is absent from the source bank",
                            )
                        })
                    })
                    .collect()
            };
            let reactants = bind(&reaction.reactants)?;
            let products = bind(&reaction.products)?;
            for &(i, order) in &reactants {
                rows[i][column] -= i128::from(order);
            }
            for &(i, order) in &products {
                rows[i][column] += i128::from(order);
            }
            let net: Vec<_> = rows
                .iter()
                .enumerate()
                .filter_map(|(i, row)| (row[column] != 0).then_some((i, row[column] as f64)))
                .collect();
            let mut dm = Sum::default();
            let mut mass_scale = Sum::default();
            for &(i, coefficient) in &net {
                let term = coefficient * thermodynamics.models()[i].molar_mass().value();
                dm.add(term);
                mass_scale.add(term.abs());
            }
            let residual = dm.value().abs() / mass_scale.value();
            if !residual.is_finite() || residual > relative_mass_tolerance {
                return Err(KineticsErrorV1::MolecularMassBalance {
                    reaction: reaction.reaction.clone(),
                    relative_residual: residual,
                });
            }
            let delta_molecules = net.iter().map(|(_, n)| *n).sum();
            bound.push(BoundReaction {
                reactants,
                products,
                net,
                delta_molecules,
                collider_logs: reaction.pressure_dependence.bind(species)?,
            });
        }
        let stoichiometry = StoichiometricMatrix::new(
            species.to_vec(),
            reactions.iter().map(|r| r.reaction.clone()).collect(),
            rows,
        )
        .map_err(KineticsErrorV1::Conservation)?;
        let conservation = verify_conservation(
            thermodynamics.elemental_matrix(),
            &stoichiometry,
            thermodynamics.charge_vector(),
        )
        .map_err(KineticsErrorV1::Conservation)?;
        Ok(Self {
            thermodynamics,
            mechanism_id,
            mechanism_version,
            reactions,
            bound,
            stoichiometry,
            conservation,
            relative_mass_tolerance,
        })
    }
    /// Original full species/element/charge/NASA-9 owner.
    pub const fn thermodynamics(&self) -> &IdealGasEquilibriumV1 {
        &self.thermodynamics
    }
    /// Caller-authored immutable mechanism identity (not source authentication).
    pub fn mechanism_id(&self) -> &str {
        &self.mechanism_id
    }
    /// Caller-authored immutable mechanism revision.
    pub const fn mechanism_version(&self) -> u32 {
        self.mechanism_version
    }
    /// Original complete reaction laws, in canonical reaction order.
    pub fn reactions(&self) -> &[ElementaryReactionV1] {
        &self.reactions
    }
    /// Exact source species-by-reaction matrix, including inert rows.
    pub const fn stoichiometry(&self) -> &StoichiometricMatrix {
        &self.stoichiometry
    }
    /// Original exact fs-qty conservation certificate.
    pub const fn conservation(&self) -> ConservationCertificate {
        self.conservation
    }
    /// Admitted relative molecular-mass source discrepancy per reaction.
    pub const fn relative_mass_tolerance(&self) -> f64 {
        self.relative_mass_tolerance
    }

    /// Require support of an entire temperature interval for every possible
    /// species. Crossing NASA regions requires explicit caloric reconciliation.
    pub fn admit_temperature_interval(
        &self,
        lower: Temperature,
        upper: Temperature,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), KineticsErrorV1> {
        if !positive(lower.value()) || !positive(upper.value()) || lower.value() > upper.value() {
            return Err(KineticsErrorV1::Invalid(
                "invalid kinetic source temperature interval",
            ));
        }
        for model in self.thermodynamics.models() {
            admit_caloric_interval(model, lower, upper, &mut cancelled)?;
        }
        Ok(())
    }

    /// Allocation-free trial thermodynamics, including exact species absence.
    /// Concentrations follow the canonical source species order, in mol/m3.
    /// Scratch contents are invalid after any refusal; immutable inputs remain.
    pub fn evaluate_thermodynamics(
        &self,
        temperature: Temperature,
        concentrations: &[f64],
        work: &mut KineticWorkspaceV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<KineticThermodynamicsV1, KineticsErrorV1> {
        poll(&mut cancelled)?;
        self.admit_work(concentrations, work)?;
        let t = temperature.value();
        if !positive(t) {
            return Err(KineticsErrorV1::Invalid("invalid kinetic temperature"));
        }
        let mut n = Sum::default();
        let mut mass = Sum::default();
        let mut energy = Sum::default();
        let mut cv = Sum::default();
        let mut entropy = Sum::default();
        let log_standard_c = fs_math::det::ln(
            self.thermodynamics
                .convention()
                .reference_pressure()
                .value(),
        ) - fs_math::det::ln(R)
            - fs_math::det::ln(t);
        for (i, model) in self.thermodynamics.models().iter().enumerate() {
            poll(&mut cancelled)?;
            let p = model.evaluate_properties(temperature).map_err(|source| {
                KineticsErrorV1::SpeciesEvaluation {
                    species: model.species().clone(),
                    source,
                }
            })?;
            let cvi = p.cp().value() - R;
            if !positive(cvi) {
                return Err(KineticsErrorV1::Invalid(
                    "kinetic source requires positive species Cv",
                ));
            }
            work.u[i] = p.u().value();
            work.cv[i] = cvi;
            work.g[i] = p.g().value();
            work.log_c[i] = if concentrations[i] == 0.0 {
                f64::NEG_INFINITY
            } else {
                fs_math::det::ln(concentrations[i])
            };
            n.add(concentrations[i]);
            mass.add(concentrations[i] * model.molar_mass().value());
            energy.add(concentrations[i] * p.u().value());
            cv.add(concentrations[i] * cvi);
            if concentrations[i] > 0.0 {
                entropy.add(
                    concentrations[i] * (p.s().value() - R * (work.log_c[i] - log_standard_c)),
                );
            }
        }
        let pressure = (R * t) * n.value();
        if !positive(pressure)
            || !positive(mass.value())
            || !positive(cv.value())
            || !energy.value().is_finite()
            || !entropy.value().is_finite()
        {
            return Err(KineticsErrorV1::Arithmetic("gas thermodynamic reduction"));
        }
        Ok(KineticThermodynamicsV1 {
            temperature,
            pressure: Pressure::new(pressure),
            total_concentration_mol_m3: n.value(),
            mass_density_kg_m3: mass.value(),
            internal_energy_density_j_m3: energy.value(),
            heat_capacity_v_density_j_m3_k: cv.value(),
            entropy_density_j_m3_k: entropy.value(),
        })
    }

    /// Coupled reaction/species rates and exact analytic concentration and
    /// temperature derivatives. Source evaluation occurs once per species,
    /// independent of the number of reactions. No successful-call allocation.
    pub fn evaluate(
        &self,
        temperature: Temperature,
        concentrations: &[f64],
        work: &mut KineticWorkspaceV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<KineticEvaluationV1, KineticsErrorV1> {
        let state =
            self.evaluate_thermodynamics(temperature, concentrations, work, &mut cancelled)?;
        let t = temperature.value();
        let ns = concentrations.len();
        work.species_rates.fill(0.0);
        work.species_rate_correction.fill(0.0);
        work.rate_jacobian.fill(0.0);
        let mut entropy_rate = Sum::default();
        let mut infinite_entropy_rate = false;
        let standard_c_log = fs_math::det::ln(
            self.thermodynamics
                .convention()
                .reference_pressure()
                .value(),
        ) - fs_math::det::ln(R)
            - fs_math::det::ln(t);
        for (j, (source, reaction)) in self.reactions.iter().zip(&self.bound).enumerate() {
            poll(&mut cancelled)?;
            let mut dg = Sum::default();
            let mut du = Sum::default();
            for &(i, coefficient) in &reaction.net {
                dg.add(coefficient * work.g[i]);
                du.add(coefficient * work.u[i]);
            }
            let log_kc = -dg.value() / (R * t) + reaction.delta_molecules * standard_c_log;
            let dlog_kc = (du.value() / (R * t)) / t;
            let scale = pressure::effective_rate(source, reaction, t, &work.log_c)?;
            let log_forward_scale = scale.log_forward;
            let dlog_forward_scale = scale.dlog_temperature;
            if !log_kc.is_finite() || !dlog_kc.is_finite() {
                return Err(KineticsErrorV1::Arithmetic(
                    "concentration equilibrium constant",
                ));
            }
            work.log_equilibrium_constants[j] = log_kc;
            let force = reaction_affinity(reaction, &work.log_c, log_kc)?;
            work.affinities[j] = force;
            if scale.inactive() {
                work.reaction_rates[j] = 0.0;
                work.temperature_derivatives[j] = 0.0;
                continue;
            }
            let log_reverse_scale = if source.reversible {
                log_forward_scale - log_kc
            } else {
                f64::NEG_INFINITY
            };
            let lf = side_log(&reaction.reactants, &work.log_c, None);
            let lr = if source.reversible {
                side_log(&reaction.products, &work.log_c, None)
            } else {
                f64::NEG_INFINITY
            };
            let forward = rate_from_log(log_forward_scale + lf)?;
            let reverse = rate_from_log(log_reverse_scale + lr)?;
            // An indeterminate both-sides-absent affinity has exact zero
            // progress. No physical affinity is manufactured for that state.
            let affinity = force.unwrap_or(0.0);
            let rate = net_rate(forward, reverse, affinity)?;
            let collider_response = scale.collider_response(lf, lr, log_kc, affinity)?;
            let dt = rate * dlog_forward_scale + reverse * dlog_kc;
            if !dt.is_finite() {
                return Err(KineticsErrorV1::Arithmetic(
                    "kinetic temperature derivative",
                ));
            }
            work.reaction_rates[j] = rate;
            work.temperature_derivatives[j] = dt;
            for i in 0..ns {
                let alpha = order(&reaction.reactants, i);
                let beta = if source.reversible {
                    order(&reaction.products, i)
                } else {
                    0
                };
                let mut derivative = if alpha == 0 && beta == 0 {
                    0.0
                } else if concentrations[i] > 0.0 {
                    let common = alpha.min(beta) as f64;
                    let extra = if alpha > beta {
                        forward * f64::from(alpha - beta)
                    } else {
                        -reverse * f64::from(beta - alpha)
                    };
                    (rate * common + extra) / concentrations[i]
                } else {
                    let df = if alpha == 1 {
                        rate_from_log(
                            log_forward_scale + side_log(&reaction.reactants, &work.log_c, Some(i)),
                        )?
                    } else {
                        0.0
                    };
                    let dr = if beta == 1 {
                        rate_from_log(
                            log_reverse_scale + side_log(&reaction.products, &work.log_c, Some(i)),
                        )?
                    } else {
                        0.0
                    };
                    // When an absent explicit catalyst occurs once on both
                    // sides, its derivative remains a near-equilibrium net
                    // rate and needs the same cancellation-safe difference.
                    if alpha == 1 && beta == 1 && force.is_some_and(f64::is_finite) {
                        net_rate(df, dr, affinity)?
                    } else {
                        df - dr
                    }
                };
                if let Some(&log_efficiency) = reaction.collider_logs.get(i) {
                    derivative += collider_response.scaled(log_efficiency)?;
                }
                if !derivative.is_finite() {
                    return Err(KineticsErrorV1::Arithmetic(
                        "kinetic concentration derivative",
                    ));
                }
                work.rate_jacobian[j * ns + i] = derivative;
            }
            for &(i, coefficient) in &reaction.net {
                add_compensated(
                    &mut work.species_rates[i],
                    &mut work.species_rate_correction[i],
                    coefficient * rate,
                );
            }
            if rate != 0.0 {
                let production = R * rate * affinity;
                if production.is_infinite() && production.is_sign_positive() {
                    infinite_entropy_rate = true;
                } else if !production.is_finite() || (source.reversible && production < 0.0) {
                    return Err(KineticsErrorV1::Arithmetic("chemical entropy production"));
                } else {
                    entropy_rate.add(production);
                }
            }
        }
        let mut energy_rate = Sum::default();
        for i in 0..ns {
            work.species_rates[i] += work.species_rate_correction[i];
            if !work.species_rates[i].is_finite() {
                return Err(KineticsErrorV1::Arithmetic("species production reduction"));
            }
            energy_rate.add(work.u[i] * work.species_rates[i]);
        }
        let temperature_rate = -energy_rate.value() / state.heat_capacity_v_density_j_m3_k;
        if !temperature_rate.is_finite() {
            return Err(KineticsErrorV1::Arithmetic("adiabatic temperature rate"));
        }
        poll(&mut cancelled)?;
        Ok(KineticEvaluationV1 {
            thermodynamics: state,
            adiabatic_temperature_rate_k_s: temperature_rate,
            entropy_production_w_m3_k: if infinite_entropy_rate {
                f64::INFINITY
            } else {
                entropy_rate.value()
            },
        })
    }

    fn admit_work(
        &self,
        concentrations: &[f64],
        work: &KineticWorkspaceV1,
    ) -> Result<(), KineticsErrorV1> {
        let ns = self.thermodynamics.models().len();
        if concentrations.len() != ns
            || work.u.len() != ns
            || work.reaction_rates.len() != self.reactions.len()
            || concentrations.iter().any(|c| !c.is_finite() || *c < 0.0)
        {
            return Err(KineticsErrorV1::Invalid(
                "invalid concentration vector or workspace shape",
            ));
        }
        Ok(())
    }
}

/// Reused bounded trial-state storage. Arrays follow the canonical mechanism.
/// Do not read scratch results after a failed/cancelled evaluation.
#[derive(Clone, Debug)]
pub struct KineticWorkspaceV1 {
    u: Vec<f64>,
    cv: Vec<f64>,
    g: Vec<f64>,
    log_c: Vec<f64>,
    reaction_rates: Vec<f64>,
    species_rates: Vec<f64>,
    species_rate_correction: Vec<f64>,
    rate_jacobian: Vec<f64>,
    temperature_derivatives: Vec<f64>,
    log_equilibrium_constants: Vec<f64>,
    affinities: Vec<Option<f64>>,
}
impl KineticWorkspaceV1 {
    /// Allocate once for the full admitted species/reaction dimensions.
    pub fn new(owner: &IdealGasKineticsV1) -> Result<Self, KineticsErrorV1> {
        let ns = owner.thermodynamics.models().len();
        let nr = owner.reactions.len();
        Ok(Self {
            u: zeros(ns)?,
            cv: zeros(ns)?,
            g: zeros(ns)?,
            log_c: zeros(ns)?,
            reaction_rates: zeros(nr)?,
            species_rates: zeros(ns)?,
            species_rate_correction: zeros(ns)?,
            rate_jacobian: zeros(nr * ns)?,
            temperature_derivatives: zeros(nr)?,
            log_equilibrium_constants: zeros(nr)?,
            affinities: filled(nr, None)?,
        })
    }
    /// Canonical signed net reaction progress, mol/(m3 s).
    pub fn reaction_rates(&self) -> &[f64] {
        &self.reaction_rates
    }
    /// Canonical signed species production, mol/(m3 s).
    pub fn species_rates(&self) -> &[f64] {
        &self.species_rates
    }
    /// Reaction-major d q_j / d c_i at fixed temperature, in1/s.
    pub fn rate_jacobian(&self) -> &[f64] {
        &self.rate_jacobian
    }
    /// d q_j / d T at fixed concentrations, mol/(m3 s K).
    pub fn temperature_derivatives(&self) -> &[f64] {
        &self.temperature_derivatives
    }
    /// Original species molar internal energies including formation references.
    pub fn species_internal_energies(&self) -> &[f64] {
        &self.u
    }
    /// Original species molar Cv, J/(mol K).
    pub fn species_heat_capacities_v(&self) -> &[f64] {
        &self.cv
    }
    /// ln Kc using the explicit1mol/m3 concentration normalization.
    pub fn log_equilibrium_constants(&self) -> &[f64] {
        &self.log_equilibrium_constants
    }
    /// Thermodynamic affinity divided by RT, independent of kinetic shutdown.
    /// One-sided exact absence gives infinity; competing absent net reactants
    /// and products give None. An absent explicit catalyst cancels from force.
    pub fn dimensionless_affinities(&self) -> &[Option<f64>] {
        &self.affinities
    }
}

/// Actual ideal-gas state functions per unit fixed volume.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct KineticThermodynamicsV1 {
    /// Evaluated source temperature.
    pub temperature: Temperature,
    /// Pressure including the current molecular count.
    pub pressure: Pressure,
    /// Sum of canonical concentrations, mol/m3.
    pub total_concentration_mol_m3: f64,
    /// Sum of c_i M_i, kg/m3.
    pub mass_density_kg_m3: f64,
    /// Sum of c_i u_i, J/m3; formation energy occurs once.
    pub internal_energy_density_j_m3: f64,
    /// Frozen composition sum of c_i Cv_i, J/(m3 K).
    pub heat_capacity_v_density_j_m3_k: f64,
    /// Ideal mixing/pressure entropy sum, J/(m3 K), with0ln0=0.
    pub entropy_density_j_m3_k: f64,
}

/// Constitutive rates at one immutable state, not a finite time-step result.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct KineticEvaluationV1 {
    /// Original source thermodynamics.
    pub thermodynamics: KineticThermodynamicsV1,
    /// dT/dt required by fixed-volume adiabatic conservation, K/s.
    pub adiabatic_temperature_rate_k_s: f64,
    /// Local R sum(q_j affinity_j), signed when original irreversible channels
    /// are present. Reversible channels retain nonnegative production; an
    /// irreversible source is not guaranteed to follow the thermodynamic force.
    /// Positive infinity denotes the exact ideal mixture boundary limit, not a
    /// finite-step entropy increment. This diagnostic never relaxes a consumer's
    /// total finite-step entropy acceptance.
    pub entropy_production_w_m3_k: f64,
}

fn side_log(side: &[(usize, u32)], log_c: &[f64], omit: Option<usize>) -> f64 {
    let mut sum = Sum::default();
    for &(i, exponent) in side {
        if omit == Some(i) {
            continue;
        }
        if log_c[i] == f64::NEG_INFINITY {
            return f64::NEG_INFINITY;
        }
        sum.add(f64::from(exponent) * log_c[i]);
    }
    sum.value()
}
fn order(side: &[(usize, u32)], i: usize) -> u32 {
    side.binary_search_by_key(&i, |&(index, _)| index)
        .ok()
        .map_or(0, |j| side[j].1)
}
fn reaction_affinity(
    reaction: &BoundReaction,
    log_c: &[f64],
    log_kc: f64,
) -> Result<Option<f64>, KineticsErrorV1> {
    let mut absent_reactant = false;
    let mut absent_product = false;
    let mut force = Sum::default();
    force.add(log_kc);
    for &(i, coefficient) in &reaction.net {
        if log_c[i] == f64::NEG_INFINITY {
            if coefficient < 0.0 {
                absent_reactant = true;
            } else {
                absent_product = true;
            }
        } else {
            force.add(-coefficient * log_c[i]);
        }
    }
    Ok(match (absent_reactant, absent_product) {
        (true, true) => None,
        (true, false) => Some(f64::NEG_INFINITY),
        (false, true) => Some(f64::INFINITY),
        (false, false) if force.value().is_finite() => Some(force.value()),
        _ => return Err(KineticsErrorV1::Arithmetic("reaction affinity")),
    })
}
fn rate_from_log(log: f64) -> Result<f64, KineticsErrorV1> {
    if log == f64::NEG_INFINITY {
        return Ok(0.0);
    }
    let value = fs_math::det::exp(log);
    if !log.is_finite() || !positive(value) {
        return Err(KineticsErrorV1::Arithmetic(
            "nonzero directional rate outside binary64 range",
        ));
    }
    Ok(value)
}
fn net_rate(forward: f64, reverse: f64, affinity: f64) -> Result<f64, KineticsErrorV1> {
    let value = if forward == 0.0 || reverse == 0.0 {
        forward - reverse
    } else if affinity > 0.0 {
        forward * -fs_math::det::expm1(-affinity)
    } else if affinity < 0.0 {
        reverse * fs_math::det::expm1(affinity)
    } else if affinity == 0.0 {
        0.0
    } else {
        return Err(KineticsErrorV1::Arithmetic("reaction affinity"));
    };
    if !value.is_finite() || (forward > 0.0 && reverse > 0.0 && affinity != 0.0 && value == 0.0) {
        return Err(KineticsErrorV1::Arithmetic(
            "nonzero net rate outside binary64 range",
        ));
    }
    Ok(value)
}
fn zeros(n: usize) -> Result<Vec<f64>, KineticsErrorV1> {
    filled(n, 0.0)
}
fn filled<T: Clone>(n: usize, value: T) -> Result<Vec<T>, KineticsErrorV1> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(n)
        .map_err(|_| KineticsErrorV1::AllocationRefused)?;
    values.resize(n, value);
    Ok(values)
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), KineticsErrorV1> {
    if cancelled() {
        Err(KineticsErrorV1::Cancelled)
    } else {
        Ok(())
    }
}
fn add_compensated(sum: &mut f64, correction: &mut f64, term: f64) {
    let next = *sum + term;
    *correction += if sum.abs() >= term.abs() {
        (*sum - next) + term
    } else {
        (term - next) + *sum
    };
    *sum = next;
}
#[derive(Default)]
struct Sum {
    sum: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, term: f64) {
        add_compensated(&mut self.sum, &mut self.correction, term);
    }
    fn value(&self) -> f64 {
        self.sum + self.correction
    }
}

#[cfg(test)]
mod tests;
