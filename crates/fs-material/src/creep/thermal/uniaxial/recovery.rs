//! Exact thermal fiber reductions with separate permanent and recovery histories.
//!
//! The native axial owners eliminate free transverse stresses analytically.
//! Their compatible tensor strains retain the actual old and new temperatures;
//! the original thermal tensor owner audits all storage and heat. Recovering
//! resistance changes the flow kinetics without creating another energy law.

use super::add_free;
use crate::creep::recovery::RecoveringCreepState;
use crate::creep::recovery::plastic::RecoveringNortonJ2State;
use crate::creep::recovery::plastic::source::{
    RecoveringNortonJ2Source, ResolvedRecoveringNortonJ2,
};
use crate::creep::recovery::plastic::uniaxial::{
    UniaxialRecoveringNortonJ2, UniaxialRecoveringNortonJ2State,
};
use crate::creep::recovery::source::{RecoveringNortonSource, ResolvedRecoveringNorton};
use crate::creep::recovery::uniaxial::{UniaxialRecoveringNorton, UniaxialRecoveringNortonState};
use crate::creep::thermal::expansion::{ExpandingCreepSource, ThermalCreep, ThermalCreepState};
use crate::creep::{CreepError, CreepIncrement, CreepState, CreepTolerances, NortonCreep};
use crate::fiber::Uniaxial;
use crate::state_point::ResolvedMaterialStatePoint;
use crate::{MaterialAdmissibility, SmallStrainLaw, Voigt};
use fs_evidence::ModelCard;

// Pure and joint recovery use the same temperature transport, axial restriction
// and energy admission. Only their original mechanical history checks differ.
macro_rules! thermal_recovery {
    ($law:ident, $state:ident, $response:ident, $mechanical:ident, $history:ident,
     $resolved:ident, $source:ident) => {
        /// Accepted free-sided mechanical history and its actual source temperature.
        #[derive(Clone, Debug, PartialEq)]
        pub struct $state {
            pub mechanical: $history,
            temperature_k: f64,
        }
        impl $state {
            pub const fn temperature_k(&self) -> f64 {
                self.temperature_k
            }
        }

        /// One uncommitted axial candidate, including compatible transverse strain.
        #[derive(Clone, Debug, PartialEq)]
        pub struct $response {
            pub strain: Voigt,
            pub stress_pa: f64,
            pub tangent_pa: f64,
            pub state: $state,
        }

        /// Source-bound axial recovery with the native tensor thermal energy owner.
        #[derive(Clone, Debug)]
        pub struct $law {
            mechanical: $mechanical,
            tensor: ThermalCreep<$resolved>,
        }
        impl $law {
            /// Original kinetic coefficients and physical interval.
            pub fn kinetics(&self) -> &NortonCreep {
                self.mechanical.kinetics()
            }
            pub fn material(&self) -> Option<&ResolvedMaterialStatePoint> {
                self.mechanical.material()
            }
            pub fn expansion_material(&self) -> &ResolvedMaterialStatePoint {
                self.tensor.expansion_material()
            }
            pub fn tensor_law(&self) -> &ThermalCreep<$resolved> {
                &self.tensor
            }
            /// Original recovery coefficients and axial history owner for rebinding.
            pub fn mechanical_law(&self) -> &$mechanical {
                &self.mechanical
            }

            /// Reconstruct accepted total strain at its retained temperature,
            /// without consuming another flow or recovery interval.
            pub fn tensor_strain(
                &self,
                axial_strain: f64,
                state: &$state,
            ) -> Result<Voigt, CreepError> {
                let free = self
                    .tensor
                    .free_strain_at_temperature(state.temperature_k)?;
                let mechanical = self
                    .mechanical
                    .tensor_strain(axial_strain - free, &state.mechanical)?;
                Ok(add_free(mechanical, free))
            }

            pub fn response(
                &self,
                axial_strain: f64,
                old: &$state,
            ) -> Result<$response, CreepError> {
                self.tensor_strain(axial_strain, old)?;
                let free = self.tensor.current_free_strain()?;
                let response = self
                    .mechanical
                    .response(axial_strain - free, &old.mechanical)?;
                Ok($response {
                    strain: add_free(response.strain, free),
                    stress_pa: response.stress_pa,
                    tangent_pa: response.tangent_pa,
                    state: $state {
                        mechanical: response.state,
                        temperature_k: self.tensor.thermoelastic().temperature_k(),
                    },
                })
            }

            /// Admit through the original thermal tensor owner. Its full tensor
            /// tangent and energy are returned; `response` supplies the condensed
            /// axial tangent. No resistance energy or recovery heat is introduced.
            pub fn advance_with_energy(
                &self,
                old_strain: f64,
                old: &$state,
                strain: f64,
                tolerances: CreepTolerances,
            ) -> Result<CreepIncrement<$state>, CreepError> {
                let old_tensor = self.tensor_strain(old_strain, old)?;
                let response = self.response(strain, old)?;
                let tensor_state = ThermalCreepState::from_parts(
                    old.mechanical.tensor_history(),
                    old.temperature_k,
                );
                let checked = self.tensor.advance_with_energy(
                    &old_tensor,
                    &tensor_state,
                    &response.strain,
                    tolerances,
                )?;
                let expected = response.state.mechanical.tensor_history();
                let mut residual = self.history_residual(&checked.state.mechanical, &expected);
                residual = residual.max((checked.stress_pa[0] - response.stress_pa).abs());
                for stress in &checked.stress_pa[1..] {
                    residual = residual.max(stress.abs());
                }
                if checked.state.temperature_k() != response.state.temperature_k
                    || !residual.is_finite()
                    || residual > tolerances.stress_pa
                {
                    return Err(CreepError::StressResidual {
                        residual_pa: residual,
                        tolerance_pa: tolerances.stress_pa,
                    });
                }
                Ok(CreepIncrement {
                    strain: response.strain,
                    stress_pa: checked.stress_pa,
                    tangent_pa: checked.tangent_pa,
                    state: response.state,
                    energy: checked.energy,
                })
            }
        }

        impl Uniaxial for $law {
            type State = $state;
            fn initial_state(&self) -> Self::State {
                Self::State {
                    mechanical: self.mechanical.initial_state(),
                    temperature_k: self.tensor.thermoelastic().temperature_k(),
                }
            }
            fn stress(&self, strain: f64, state: &Self::State) -> f64 {
                self.response(strain, state)
                    .map(|v| v.stress_pa)
                    .unwrap_or(f64::NAN)
            }
            fn tangent(&self, strain: f64, state: &Self::State) -> f64 {
                self.response(strain, state)
                    .map(|v| v.tangent_pa)
                    .unwrap_or(f64::NAN)
            }
            fn update_state(&self, strain: f64, state: &Self::State) -> Self::State {
                self.response(strain, state)
                    .map(|v| v.state)
                    .unwrap_or(Self::State {
                        mechanical: state.mechanical.clone(),
                        temperature_k: f64::NAN,
                    })
            }
            fn admissibility(&self) -> MaterialAdmissibility {
                self.tensor.admissibility()
            }
            fn card(&self) -> ModelCard {
                self.tensor.card()
            }
        }

        impl<'a> ExpandingCreepSource<$source<'a>> {
            /// Resolve the tensor source once, retaining the same complete
            /// mechanical receipts in its exact free-sided axial reduction.
            pub fn at_uniaxial_temperature(
                &self,
                temperature_k: f64,
                duration_s: f64,
            ) -> Result<$law, CreepSourceError> {
                let tensor = self.at_temperature(temperature_k, duration_s)?;
                let mechanical = $mechanical::from_source(tensor.mechanical().clone());
                Ok($law { mechanical, tensor })
            }
        }
    };
}

use crate::creep::source::CreepSourceError;

thermal_recovery!(
    UniaxialThermalRecoveringNorton,
    UniaxialThermalRecoveringNortonState,
    UniaxialThermalRecoveringNortonResponse,
    UniaxialRecoveringNorton,
    UniaxialRecoveringNortonState,
    ResolvedRecoveringNorton,
    RecoveringNortonSource
);
thermal_recovery!(
    UniaxialThermalRecoveringNortonJ2,
    UniaxialThermalRecoveringNortonJ2State,
    UniaxialThermalRecoveringNortonJ2Response,
    UniaxialRecoveringNortonJ2,
    UniaxialRecoveringNortonJ2State,
    ResolvedRecoveringNortonJ2,
    RecoveringNortonJ2Source
);

fn creep_residual(
    e: f64,
    actual: &CreepState,
    expected: &CreepState,
    actual_resistance: f64,
    expected_resistance: f64,
) -> f64 {
    actual
        .creep_strain
        .iter()
        .copied()
        .zip(expected.creep_strain)
        .chain([
            (
                actual.equivalent_creep_strain,
                expected.equivalent_creep_strain,
            ),
            (actual_resistance, expected_resistance),
        ])
        .map(|(a, b)| e * (a - b).abs())
        .fold(0.0_f64, f64::max)
}

impl UniaxialThermalRecoveringNorton {
    fn history_residual(
        &self,
        actual: &RecoveringCreepState,
        expected: &RecoveringCreepState,
    ) -> f64 {
        creep_residual(
            self.kinetics().elastic().youngs,
            &actual.native,
            &expected.native,
            actual.resistance_strain,
            expected.resistance_strain,
        )
    }
}

impl UniaxialThermalRecoveringNortonJ2 {
    pub fn hardening(&self) -> f64 {
        self.mechanical.backbone().norton_j2().plastic().hardening
    }
    pub fn yield_stress(&self) -> f64 {
        self.mechanical
            .backbone()
            .norton_j2()
            .plastic()
            .yield_stress
    }

    fn history_residual(
        &self,
        actual: &RecoveringNortonJ2State,
        expected: &RecoveringNortonJ2State,
    ) -> f64 {
        let e = self.kinetics().elastic().youngs;
        let mut residual = creep_residual(
            e,
            &actual.native.creep,
            &expected.native.creep,
            actual.resistance_strain,
            expected.resistance_strain,
        );
        for (a, b) in actual
            .native
            .plastic
            .plastic_strain
            .iter()
            .zip(expected.native.plastic.plastic_strain)
        {
            residual = residual.max(e * (a - b).abs());
        }
        residual.max(
            e.max(self.hardening())
                * (actual.native.plastic.alpha - expected.native.plastic.alpha).abs(),
        )
    }
}
