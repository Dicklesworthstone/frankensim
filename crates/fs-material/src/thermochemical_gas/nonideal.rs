//! Explicit Peng--Robinson binding over the existing thermochemical species bank.
//!
//! The ideal-gas `ThermochemicalGasModel::evaluate` remains unchanged. This
//! module exposes nonideal mechanics separately so consumers cannot combine PR
//! density with ideal-gas sound speed or caloric departures by accident.

use fs_matdb::QueryPoint;
use fs_qty::{Density, Pressure, Temperature};
use fs_thermochem::SpeciesId;
use fs_thermochem::eos::peng_robinson::{
    PengRobinsonComponentV1, PengRobinsonEvaluationV1, PengRobinsonMixtureV1,
    PengRobinsonPhaseRootV1, PengRobinsonPhaseV1,
    flash::{PengRobinsonFlashOptionsV1, PengRobinsonFlashV1},
};
use fs_thermochem::mixture::FrozenIdealGasMixtureEvaluationV1;

use super::{
    ResolvedThermochemicalGasState, ThermochemicalGasModel, coordinates, failure, invalid,
};
use crate::MaterialError;
use crate::fluid_transport::{
    FluidTransport, FluidTransportState, ThermalFluidTransport, ThermalFluidTransportState,
};

/// Conservative density/energy reconstruction on a homogeneous PR branch.
pub mod energy;

#[derive(Clone, Debug, PartialEq)]
pub struct PengRobinsonSpeciesBinding {
    pub species: SpeciesId,
    pub critical_temperature: Temperature,
    pub critical_pressure: Pressure,
    pub acentric_factor: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct ThermochemicalPengRobinsonModel {
    source: ThermochemicalGasModel,
    mechanics: PengRobinsonMixtureV1,
}

impl ThermochemicalPengRobinsonModel {
    /// Bind critical properties and pair interactions to the original species.
    ///
    /// A nonempty `binary_interactions` matrix is row-major in the supplied
    /// `bindings` order, on both axes. Species and both matrix axes are permuted
    /// together into the source model's canonical order. An empty matrix means
    /// zero pair interactions, as in the shared Peng--Robinson owner.
    pub fn new(
        source: ThermochemicalGasModel,
        bindings: Vec<PengRobinsonSpeciesBinding>,
        binary_interactions: Vec<f64>,
    ) -> Result<Self, MaterialError> {
        let components = source.thermodynamics().components();
        if bindings.len() != components.len() {
            return Err(invalid(
                "real-gas binding requires exactly one critical-property record per active species",
            ));
        }
        let mut order: Vec<usize> = (0..bindings.len()).collect();
        order.sort_by(|&a, &b| bindings[a].species.cmp(&bindings[b].species));
        if order.iter().zip(components).any(
            |(&index, component)| &bindings[index].species != component.species()
        ) {
            return Err(invalid(
                "real-gas critical properties contain duplicate, missing or foreign species",
            ));
        }
        let binary_interactions = if binary_interactions.is_empty() {
            binary_interactions
        } else {
            let n = bindings.len();
            if binary_interactions.len() != n * n {
                return Err(invalid(
                    "real-gas binary interactions require a square matrix in the supplied binding order",
                ));
            }
            let mut canonical = Vec::with_capacity(n * n);
            for &row in &order {
                for &column in &order {
                    canonical.push(binary_interactions[row * n + column]);
                }
            }
            canonical
        };
        let pr_components = order.iter().zip(components).map(
            |(&index, component)| PengRobinsonComponentV1 {
                critical_temperature: bindings[index].critical_temperature,
                critical_pressure: bindings[index].critical_pressure,
                acentric_factor: bindings[index].acentric_factor,
                molar_mass: component.molar_mass(),
            }
        ).collect();
        let mechanics = PengRobinsonMixtureV1::new(
            pr_components,
            source.thermodynamics().mole_composition().fractions().to_vec(),
            binary_interactions,
        ).map_err(|error| failure("Peng-Robinson binding", error))?;
        Ok(Self { source, mechanics })
    }

    #[must_use]
    pub const fn source(&self) -> &ThermochemicalGasModel { &self.source }

    #[must_use]
    pub const fn mechanics(&self) -> &PengRobinsonMixtureV1 { &self.mechanics }

    pub fn evaluate_mechanics(
        &self,
        point: &QueryPoint,
    ) -> Result<ResolvedThermochemicalPengRobinsonState, MaterialError> {
        let (temperature, pressure) = coordinates(point)?;
        let thermodynamics = self.source.thermodynamics().evaluate(temperature, pressure)
            .map_err(|error| failure("nonideal source thermodynamics", error))?;
        let mechanics = self.mechanics.evaluate_pt(temperature, pressure)
            .map_err(|error| failure("Peng-Robinson mechanics", error))?;
        Ok(ResolvedThermochemicalPengRobinsonState {
            temperature, pressure, thermodynamics, mechanics,
        })
    }


    /// Resolve one mechanically stable PR branch into a complete nonideal state.
    ///
    /// NASA supplies the ideal-reference calorics, PR supplies residual u/h/Cv
    /// plus exact pressure derivatives, and the original material cards supply
    /// viscosity/conductivity. The legacy ideal GasState is retained only as
    /// source provenance and is never reinterpreted as this dense state.
    pub fn evaluate_real_gas(
        &self,
        point: &QueryPoint,
        phase: PengRobinsonPhaseV1,
    ) -> Result<ResolvedThermochemicalRealGasState, MaterialError> {
        let source = self.source.evaluate(point)?;
        let temperature = Temperature::new(source.state().temperature);
        let pressure = Pressure::new(source.state().pressure);
        let mechanics = self.mechanics.evaluate_pt(temperature, pressure)
            .map_err(|error| failure("Peng-Robinson mechanics", error))?;
        let root = mechanics.select(phase).clone();
        Self::resolve_real_gas(source, mechanics, root)
    }

    /// Resolve the same constitutive state at the supplied mass density.
    ///
    /// The point supplies temperature and retained material coordinates; its
    /// pressure is replaced by the PR pressure at v=M/rho before source
    /// evaluation. The given density is retained exactly, with no selection of
    /// a different PT root. Requires positive pressure, (dp/dv)_T<0 and Cv>0.
    /// Mechanically stable metastable homogeneous states remain admissible;
    /// this does not perform an equilibrium flash or two-phase reconstruction.
    pub fn evaluate_density_temperature(
        &self, point: &QueryPoint, density_kg_m3: f64,
    ) -> Result<ResolvedThermochemicalRealGasState, MaterialError> {
        let (temperature, _) = coordinates(point)?;
        let direct = self.mechanics.evaluate_density_temperature(
            temperature, Density::new(density_kg_m3),
        ).map_err(|error| failure("Peng-Robinson density mechanics", error))?;
        let point = super::energy::updated_point(point, temperature.value(), direct.pressure().value())
            .map_err(|error| failure("density state coordinates", error))?;
        let source = self.source.evaluate(&point)?;
        // Keep the established PT root inventory available to existing
        // observers. Only this final material resolution solves a cubic;
        // caloric inversion uses the direct density evaluator throughout.
        let mechanics = self.mechanics.evaluate_pt(temperature, direct.pressure())
            .map_err(|error| failure("Peng-Robinson root inventory", error))?;
        if !mechanics.phase_roots().iter().any(|root| {
            (root.mass_density_kg_per_m3() / density_kg_m3 - 1.0).abs()
                <= 16_384.0 * f64::EPSILON
        }) {
            return Err(invalid("Peng-Robinson PT inventory cannot resolve the retained density branch"));
        }
        Self::resolve_real_gas(source, mechanics, direct.root().clone())
    }

    fn resolve_real_gas(
        source: ResolvedThermochemicalGasState,
        mechanics: PengRobinsonEvaluationV1,
        root: PengRobinsonPhaseRootV1,
    ) -> Result<ResolvedThermochemicalRealGasState, MaterialError> {
        let properties = source.thermodynamics().properties();
        let molar_mass = properties.molar_mass().value();
        let ideal_cv_molar = properties.mass_specific_cv().value() * molar_mass;
        let caloric = root.caloric_mechanics(ideal_cv_molar, molar_mass)
            .ok_or_else(|| invalid("selected Peng-Robinson root is mechanically or calorically unstable"))?;
        let cp = caloric.molar_cp_j_mol_k / molar_mass;
        let cv = caloric.molar_cv_j_mol_k / molar_mass;
        let density = root.mass_density_kg_per_m3();
        let sound_speed = caloric.sound_speed_m_s;
        let source_transport = source.transport_state()?;
        let transport = FluidTransportState::try_new(
            density,
            cp,
            source_transport.thermal_conductivity_w_m_k(),
            source_transport.dynamic_viscosity_pa_s(),
        )?;
        let specific_internal_energy_j_kg = source.specific_internal_energy_j_kg()
            + root.residual_internal_energy_j_mol() / molar_mass;
        let specific_enthalpy_j_kg = source.specific_enthalpy_j_kg()
            + root.residual_enthalpy_j_mol() / molar_mass;
        let gamma = cp / cv;
        let specific_gas_constant = fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / molar_mass;
        let characteristic_impedance = density * sound_speed;
        // At fixed composition, dp = p_T dT + p_v dv. The SAME selected
        // mechanically stable root therefore gives beta = (1/v)(dv/dT)_p
        // = -p_T/(v*p_v). Do not attach ideal-gas 1/T to a PR density.
        let volumetric_expansion_per_k =
            -(caloric.dp_dt_v_pa_k / caloric.dp_dv_t_pa_mol_m3) / root.molar_volume_m3_per_mol();
        if ![
            cp, cv, density, sound_speed, gamma, specific_gas_constant,
            characteristic_impedance, specific_internal_energy_j_kg,
            specific_enthalpy_j_kg, volumetric_expansion_per_k,
        ].into_iter().all(|value| value.is_finite())
            || cp <= 0.0 || cv <= 0.0 || density <= 0.0 || sound_speed <= 0.0
            || gamma <= 0.0 || specific_gas_constant <= 0.0 || characteristic_impedance <= 0.0
            || (caloric.dp_dt_v_pa_k != 0.0 && volumetric_expansion_per_k == 0.0)
        {
            return Err(invalid("resolved real-gas state is unrepresentable"));
        }
        Ok(ResolvedThermochemicalRealGasState {
            source,
            mechanics,
            root,
            transport,
            density_kg_m3: density,
            sound_speed_m_s: sound_speed,
            specific_heat_cp_j_kg_k: cp,
            specific_heat_cv_j_kg_k: cv,
            gamma,
            specific_gas_constant_j_kg_k: specific_gas_constant,
            characteristic_impedance_pa_s_m: characteristic_impedance,
            specific_internal_energy_j_kg,
            specific_enthalpy_j_kg,
            volumetric_expansion_per_k,
        })
    }

    pub fn two_phase_flash(
        &self,
        point: &QueryPoint,
        options: PengRobinsonFlashOptionsV1,
    ) -> Result<PengRobinsonFlashV1, MaterialError> {
        let (temperature, pressure) = coordinates(point)?;
        self.source.thermodynamics().evaluate(temperature, pressure)
            .map_err(|error| failure("nonideal source thermodynamics", error))?;
        self.mechanics.two_phase_flash_pt(temperature, pressure, options)
            .map_err(|error| failure("Peng-Robinson flash", error))
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedThermochemicalPengRobinsonState {
    temperature: Temperature,
    pressure: Pressure,
    thermodynamics: FrozenIdealGasMixtureEvaluationV1,
    mechanics: PengRobinsonEvaluationV1,
}
impl ResolvedThermochemicalPengRobinsonState {
    #[must_use] pub const fn temperature(&self) -> Temperature { self.temperature }
    #[must_use] pub const fn pressure(&self) -> Pressure { self.pressure }
    #[must_use] pub const fn thermodynamics(&self) -> &FrozenIdealGasMixtureEvaluationV1 {
        &self.thermodynamics
    }
    #[must_use] pub const fn mechanics(&self) -> &PengRobinsonEvaluationV1 { &self.mechanics }
}


#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedThermochemicalRealGasState {
    source: ResolvedThermochemicalGasState,
    mechanics: PengRobinsonEvaluationV1,
    root: PengRobinsonPhaseRootV1,
    transport: FluidTransportState,
    density_kg_m3: f64,
    sound_speed_m_s: f64,
    specific_heat_cp_j_kg_k: f64,
    specific_heat_cv_j_kg_k: f64,
    gamma: f64,
    specific_gas_constant_j_kg_k: f64,
    characteristic_impedance_pa_s_m: f64,
    specific_internal_energy_j_kg: f64,
    specific_enthalpy_j_kg: f64,
    volumetric_expansion_per_k: f64,
}
impl ResolvedThermochemicalRealGasState {
    /// Actual absolute temperature shared by reference and real-fluid views.
    #[must_use]pub const fn temperature_k(&self) -> f64 { self.source.state().temperature }
    /// Actual EOS pressure used for all original material source queries.
    #[must_use]pub const fn pressure_pa(&self) -> f64 { self.source.state().pressure }
    #[must_use]pub const fn source(&self) -> &ResolvedThermochemicalGasState { &self.source }
    #[must_use]pub const fn mechanics(&self) -> &PengRobinsonEvaluationV1 { &self.mechanics }
    #[must_use]pub const fn root(&self) -> &PengRobinsonPhaseRootV1 { &self.root }
    #[must_use]pub const fn density_kg_m3(&self) -> f64 { self.density_kg_m3 }
    #[must_use]pub const fn sound_speed_m_s(&self) -> f64 { self.sound_speed_m_s }
    #[must_use]pub const fn specific_heat_cp_j_kg_k(&self) -> f64 { self.specific_heat_cp_j_kg_k }
    #[must_use]pub const fn specific_heat_cv_j_kg_k(&self) -> f64 { self.specific_heat_cv_j_kg_k }
    #[must_use]pub const fn gamma(&self) -> f64 { self.gamma }
    #[must_use]pub const fn specific_gas_constant_j_kg_k(&self) -> f64 {
        self.specific_gas_constant_j_kg_k
    }
    #[must_use]pub const fn characteristic_impedance_pa_s_m(&self) -> f64 {
        self.characteristic_impedance_pa_s_m
    }
    #[must_use]pub const fn specific_internal_energy_j_kg(&self) -> f64 {
        self.specific_internal_energy_j_kg
    }
    #[must_use]pub const fn specific_enthalpy_j_kg(&self) -> f64 {
        self.specific_enthalpy_j_kg
    }
    /// Local isobaric volumetric expansion [1/K] of the same frozen-composition
    /// PR branch as density, Cp and sound speed. This is an EOS derivative,
    /// not an ideal-gas assumption, a phase-equilibrium derivative or a
    /// temperature interval over which a Boussinesq model is automatically valid.
    #[must_use]
    pub const fn volumetric_expansion_per_k(&self) -> f64 {
        self.volumetric_expansion_per_k
    }
}
impl FluidTransport for ResolvedThermochemicalRealGasState {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> { Ok(self.transport) }
}
impl ThermalFluidTransport for ResolvedThermochemicalRealGasState {
    fn thermal_transport_state(&self) -> Result<ThermalFluidTransportState, MaterialError> {
        ThermalFluidTransportState::try_new(
            self.transport,
            self.source.state().temperature,
            Some(self.volumetric_expansion_per_k),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::thermochemical_gas::tests::{model, point, species, transport};
    use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;
    use fs_thermochem::{Composition, CompositionBasis};

    fn bindings() -> Vec<PengRobinsonSpeciesBinding> {
        vec![
            PengRobinsonSpeciesBinding {
                species: SpeciesId::new("B").unwrap(),
                critical_temperature: Temperature::new(150.687),
                critical_pressure: Pressure::new(4.863e6),
                acentric_factor: -0.00219,
            },
            PengRobinsonSpeciesBinding {
                species: SpeciesId::new("A").unwrap(),
                critical_temperature: Temperature::new(190.564),
                critical_pressure: Pressure::new(4.5992e6),
                acentric_factor: 0.01142,
            },
        ]
    }

    #[test]
    fn canonical_species_bind_to_original_composition_and_molar_mass() {
        let source = model(vec![0.75, 0.25], 0.0);
        let bound = ThermochemicalPengRobinsonModel::new(
            source.clone(), bindings(), vec![]
        ).unwrap();
        let state = bound.evaluate_mechanics(&point(300.0, 5.0e6)).unwrap();
        assert_eq!(
            state.mechanics().molar_mass().value(),
            state.thermodynamics().properties().molar_mass().value()
        );
        assert_ne!(
            state.mechanics().phase_roots()[0].compressibility_factor().to_bits(),
            1.0f64.to_bits()
        );
        assert_eq!(bound.source(), &source);
    }

    fn ternary_source() -> ThermochemicalGasModel {
        // Manufactured caloric/transport data; this exercises species binding,
        // not empirical accuracy of a particular real mixture.
        ThermochemicalGasModel::new(
            FrozenIdealGasMixtureModelV1::new(
                vec![
                    species("A", 0.016_042_46, 3.5, 0.0),
                    species("B", 0.039_948, 2.5, 0.0),
                    species("C", 0.044_009_5, 4.5, 0.0),
                ],
                Composition::new(CompositionBasis::MoleFraction, vec![0.25, 0.25, 0.5])
                    .unwrap(),
            )
            .unwrap(),
            vec![
                transport("A", 1.1e-5, 0.03),
                transport("B", 2.2e-5, 0.02),
                transport("C", 1.5e-5, 0.025),
            ],
            None,
        )
        .unwrap()
    }

    fn ternary_bindings() -> Vec<PengRobinsonSpeciesBinding> {
        let mut values = bindings();
        values.swap(0, 1);
        values.push(PengRobinsonSpeciesBinding {
            species: SpeciesId::new("C").unwrap(),
            critical_temperature: Temperature::new(304.1282),
            critical_pressure: Pressure::new(7.3773e6),
            acentric_factor: 0.22394,
        });
        values
    }

    #[test]
    fn g3_ternary_pair_interactions_follow_species_through_every_permutation() {
        let source = ternary_source();
        let canonical_bindings = ternary_bindings();
        let interactions = vec![0.0, 0.02, 0.12, 0.02, 0.0, -0.03, 0.12, -0.03, 0.0];
        let canonical = ThermochemicalPengRobinsonModel::new(
            source.clone(), canonical_bindings.clone(), interactions.clone(),
        ).unwrap();
        for order in [[0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]] {
            let supplied = order.iter().map(|&i| canonical_bindings[i].clone()).collect();
            let mut supplied_interactions = Vec::new();
            for row in order {
                for column in order {
                    supplied_interactions.push(interactions[3 * row + column]);
                }
            }
            let permuted = ThermochemicalPengRobinsonModel::new(
                source.clone(), supplied, supplied_interactions,
            ).unwrap();
            assert_eq!(permuted.mechanics(), canonical.mechanics());
            for (temperature, pressure) in [(260.0, 1.0e6), (300.0, 5.0e6), (450.0, 8.0e6)] {
                let query = point(temperature, pressure);
                // The complete state includes density, calorics, sound speed,
                // expansion and transport, all of which must remain identical.
                assert_eq!(
                    permuted.evaluate_real_gas(&query, PengRobinsonPhaseV1::Vapor).unwrap(),
                    canonical.evaluate_real_gas(&query, PengRobinsonPhaseV1::Vapor).unwrap(),
                );
            }
        }
    }

    #[test]
    fn g1_ternary_density_matches_original_pair_assignment_in_eos_owner() {
        let source = ternary_source();
        let bindings = ternary_bindings();
        // Input order C,A,B: k_CA=0.12, k_CB=-0.03, k_AB=0.02.
        let wrapped = ThermochemicalPengRobinsonModel::new(
            source.clone(),
            vec![bindings[2].clone(), bindings[0].clone(), bindings[1].clone()],
            vec![0.0, 0.12, -0.03, 0.12, 0.0, 0.02, -0.03, 0.02, 0.0],
        ).unwrap();
        let direct = PengRobinsonMixtureV1::new(
            bindings.iter().zip(source.thermodynamics().components()).map(|(binding, species)| {
                PengRobinsonComponentV1 {
                    critical_temperature: binding.critical_temperature,
                    critical_pressure: binding.critical_pressure,
                    acentric_factor: binding.acentric_factor,
                    molar_mass: species.molar_mass(),
                }
            }).collect(),
            vec![0.25, 0.25, 0.5],
            vec![0.0, 0.02, 0.12, 0.02, 0.0, -0.03, 0.12, -0.03, 0.0],
        ).unwrap();
        let direct_state = direct.evaluate_pt(Temperature::new(300.0), Pressure::new(5.0e6)).unwrap();
        let wrapped_state = wrapped.evaluate_mechanics(&point(300.0, 5.0e6)).unwrap();
        assert_eq!(wrapped_state.mechanics(), &direct_state);
        // Independent scalar PR cubic solved with numpy.roots, then rho=M*p/(Z*R*T).
        // The old misplaced interactions produce 85.73957547172306 kg/m3.
        let density = wrapped_state.mechanics().select(PengRobinsonPhaseV1::Vapor)
            .mass_density_kg_per_m3();
        assert!((density - 85.102_022_789_600_11).abs() < 2.0e-10);
    }

    #[test]
    fn g0_permuted_interactions_refuse_bad_shape_and_invalid_pairs() {
        let source = ternary_source();
        let mut bindings = ternary_bindings();
        bindings.rotate_left(1);
        for bad in [
            vec![0.0; 8],
            vec![0.0; 10],
            vec![0.0, 0.1, 0.2, 0.3, 0.0, 0.4, 0.2, 0.4, 0.0],
            vec![0.0, f64::NAN, 0.2, f64::NAN, 0.0, 0.3, 0.2, 0.3, 0.0],
            vec![0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        ] {
            assert!(ThermochemicalPengRobinsonModel::new(source.clone(), bindings.clone(), bad).is_err());
        }
        let zero = ThermochemicalPengRobinsonModel::new(source, bindings, vec![]).unwrap();
        assert_eq!(zero.mechanics().binary_interactions(), &[0.0; 9]);
    }

    #[test]
    fn wrapped_mechanics_match_direct_pr_owner_exactly() {
        let bound = ThermochemicalPengRobinsonModel::new(
            model(vec![0.75, 0.25], 0.0),
            bindings(),
            vec![0.0, 0.03, 0.03, 0.0],
        ).unwrap();
        let query = point(280.0, 8.0e6);
        let wrapped = bound.evaluate_mechanics(&query).unwrap();
        let direct = bound.mechanics().evaluate_pt(
            wrapped.temperature(), wrapped.pressure()
        ).unwrap();
        assert_eq!(*wrapped.mechanics(), direct);
    }

    #[test]
    fn real_gas_state_closes_pr_calorics_acoustics_and_shared_transport() {
        let bound = ThermochemicalPengRobinsonModel::new(
            model(vec![0.75, 0.25], 0.0), bindings(), vec![]
        ).unwrap();
        let dense = bound.evaluate_real_gas(
            &point(300.0, 5.0e6), PengRobinsonPhaseV1::Stable
        ).unwrap();
        assert_eq!(dense.density_kg_m3(), dense.root().mass_density_kg_per_m3());
        assert_eq!(
            dense.transport_state().unwrap().density_kg_m3(),
            dense.density_kg_m3()
        );
        assert_eq!(
            dense.transport_state().unwrap().specific_heat_j_kg_k(),
            dense.specific_heat_cp_j_kg_k()
        );
        assert_eq!(
            dense.characteristic_impedance_pa_s_m(),
            dense.density_kg_m3() * dense.sound_speed_m_s()
        );
        assert!((dense.gamma()
            - dense.specific_heat_cp_j_kg_k() / dense.specific_heat_cv_j_kg_k()).abs()
            < 8.0 * f64::EPSILON);
        assert_ne!(
            dense.density_kg_m3().to_bits(),
            dense.source().state().density.to_bits()
        );
        assert_ne!(
            dense.sound_speed_m_s().to_bits(),
            dense.source().state().sound_speed.to_bits()
        );
        let acoustic =
            crate::fluid_acoustics::FluidAcoustics::acoustic_state(&dense).unwrap();
        assert_eq!(acoustic.density_kg_m3().to_bits(), dense.density_kg_m3().to_bits());
        assert_eq!(
            acoustic.sound_speed_m_s().to_bits(),
            dense.sound_speed_m_s().to_bits()
        );
        assert_eq!(
            acoustic.characteristic_impedance_pa_s_m().to_bits(),
            dense.characteristic_impedance_pa_s_m().to_bits()
        );
        assert!(dense.specific_internal_energy_j_kg().is_finite());
        assert!(dense.specific_enthalpy_j_kg().is_finite());
    }

    #[test]
    fn g1_real_gas_expansion_matches_independent_isobaric_density_differences() {
        for fractions in [vec![0.75, 0.25], vec![0.25, 0.75]] {
            let bound =
                ThermochemicalPengRobinsonModel::new(model(fractions, 0.001), bindings(), vec![])
                    .unwrap();
            for temperature in [260.0, 300.0, 450.0] {
                for pressure in [1.0e5, 5.0e6, 8.0e6] {
                    let evaluate = |t| {
                        bound
                            .evaluate_real_gas(&point(t, pressure), PengRobinsonPhaseV1::Stable)
                            .unwrap()
                    };
                    let state = evaluate(temperature);
                    let delta = 0.001;
                    let numerical = -(evaluate(temperature + delta).density_kg_m3()
                        - evaluate(temperature - delta).density_kg_m3())
                        / (2.0 * delta * state.density_kg_m3());
                    let beta = state.volumetric_expansion_per_k();
                    assert!(
                        (beta / numerical - 1.0).abs() < 2.0e-7,
                        "T={temperature}, p={pressure}: {beta} != {numerical}"
                    );
                    let thermal = state.thermal_transport_state().unwrap();
                    assert_eq!(thermal.transport(), state.transport_state().unwrap());
                    assert_eq!(thermal.temperature_k(), temperature);
                    assert_eq!(thermal.volumetric_expansion_per_k(), Some(beta));
                }
            }
        }
    }

    #[test]
    fn g3_real_gas_thermal_view_retains_dense_response_and_ideal_limit() {
        let source = crate::thermochemical_gas::tests::model_in_domains(
            vec![0.75, 0.25],
            0.0,
            [200.0, 2000.0],
            [200.0, 2000.0],
            [1.0e-3, 1.0e7],
        );
        let bound = ThermochemicalPengRobinsonModel::new(source, bindings(), vec![]).unwrap();
        let dilute = bound
            .evaluate_real_gas(&point(300.0, 0.01), PengRobinsonPhaseV1::Vapor)
            .unwrap();
        assert!((dilute.volumetric_expansion_per_k() * 300.0 - 1.0).abs() < 2.0e-8);
        let dense = bound
            .evaluate_real_gas(&point(300.0, 5.0e6), PengRobinsonPhaseV1::Stable)
            .unwrap();
        assert!((dense.volumetric_expansion_per_k() * 300.0 - 1.0).abs() > 0.01);
        // A heterogeneous consumer borrows the actual original source, not a
        // detached transport tuple or the retained ideal-reference GasState.
        let views: [&dyn ThermalFluidTransport; 2] = [&dilute, &dense];
        let first = views[0].thermal_transport_state().unwrap();
        let second = views[1].thermal_transport_state().unwrap();
        assert_ne!(first.transport(), second.transport());
        assert_ne!(
            first.volumetric_expansion_per_k(),
            second.volumetric_expansion_per_k()
        );
        assert_eq!(second.transport(), dense.transport_state().unwrap());
        assert_ne!(
            second.transport(),
            dense.source().transport_state().unwrap()
        );
    }

    #[test]
    fn foreign_species_refuses_instead_of_reassigning_critical_data() {
        let mut bad = bindings();
        bad[0].species = SpeciesId::new("FOREIGN").unwrap();
        assert!(ThermochemicalPengRobinsonModel::new(
            model(vec![0.75, 0.25], 0.0), bad, vec![]
        ).is_err());
    }
}
