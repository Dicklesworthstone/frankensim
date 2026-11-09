//! Mixed admission for the complete original finite kinetic owners.
//!
//! The shared adapter borrows each original law and its complete typed history.
//! Only volumetric assembly is selected; viscosity, the frozen interval, joint
//! flow, primary resistance and recovery remain in their native return maps.

use super::{FrozenFiniteInelasticRemainder, MixedFiniteInelasticMaterial};
use crate::plastic::finite::combined::FiniteNortonJ2;
use crate::plastic::finite::creep::FiniteNortonCreep;
use crate::plastic::finite::creep::recovery::{
    FiniteRecoveringCreepState, FiniteRecoveringNortonCreep, FiniteRecoveringNortonJ2,
};
use crate::plastic::finite::rate::{FiniteOverstressJ2, validate_viscous_committed_selected};
use crate::plastic::finite::{FiniteJ2, FiniteJ2Error, FiniteJ2Response, FiniteJ2State};

fn stored(backbone: &FiniteJ2, f: &[f64; 9], state: &FiniteJ2State) -> Result<f64, FiniteJ2Error> {
    let endpoint = backbone.elastic_endpoint(f, state)?;
    backbone.check_elastic_domain(endpoint.log_j, &endpoint.deviator)?;
    backbone.remainder_energy(&endpoint.deviator, state.alpha)
}

fn validate(
    backbone: &FiniteJ2,
    f: &[f64; 9],
    state: &FiniteJ2State,
    tolerance: f64,
    equilibrium_plastic_yield: bool,
) -> Result<(), FiniteJ2Error> {
    if equilibrium_plastic_yield {
        backbone.validate_committed_selected(f, state, tolerance, false)
    } else {
        validate_viscous_committed_selected(backbone, f, state, tolerance, false)
    }
}

fn norton_backbone(law: &FiniteNortonCreep) -> &FiniteJ2 {
    law.mixed_native().backbone()
}

fn recovering_norton_backbone(law: &FiniteRecoveringNortonCreep) -> &FiniteJ2 {
    norton_backbone(law.norton())
}

fn recovering_combined_backbone(law: &FiniteRecoveringNortonJ2) -> &FiniteJ2 {
    law.norton_j2().backbone()
}

// These owners already share the exact native FiniteJ2State and ReturnMap.
// Dispatch retains the original concrete law, never its inspection backbone.
macro_rules! native_kinetic {
    ($law:ty, $backbone:path, $equilibrium_yield:expr) => {
        impl MixedFiniteInelasticMaterial for $law {
            fn freeze_remainder<'a>(
                &'a self,
                accepted: &'a Self::State,
            ) -> Result<FrozenFiniteInelasticRemainder<'a, Self>, FiniteJ2Error> {
                self.validate_history(accepted)?;
                Ok(FrozenFiniteInelasticRemainder {
                    law: self,
                    accepted,
                    backbone: $backbone(self),
                    native_accepted: accepted,
                    prepare_return: |law, f, old| law.return_map_selected(f, old, false),
                })
            }

            fn remainder_response(
                &self,
                f: &[f64; 9],
                old: &Self::State,
            ) -> Result<FiniteJ2Response<Self::State>, FiniteJ2Error> {
                let returned = self.return_map_selected(f, old, false)?;
                $backbone(self).response_from_return(f, old, returned)
            }

            fn stored_remainder_energy(
                &self,
                f: &[f64; 9],
                state: &Self::State,
            ) -> Result<f64, FiniteJ2Error> {
                stored($backbone(self), f, state)
            }

            fn validate_remainder_committed(
                &self,
                f: &[f64; 9],
                state: &Self::State,
                yield_tolerance_pa: f64,
            ) -> Result<(), FiniteJ2Error> {
                validate(
                    $backbone(self),
                    f,
                    state,
                    yield_tolerance_pa,
                    $equilibrium_yield,
                )
            }
        }
    };
}

native_kinetic!(FiniteOverstressJ2, FiniteOverstressJ2::backbone, false);
native_kinetic!(FiniteNortonCreep, norton_backbone, false);
native_kinetic!(FiniteNortonJ2, FiniteNortonJ2::backbone, true);

// Recovery's independent resistance is validated and advanced by its original
// owner. Reconstructing the native tensor tangent does not advance it again.
macro_rules! recovering_kinetic {
    ($law:ty, $backbone:path, $equilibrium_yield:expr) => {
        impl MixedFiniteInelasticMaterial for $law {
            fn freeze_remainder<'a>(
                &'a self,
                accepted: &'a Self::State,
            ) -> Result<FrozenFiniteInelasticRemainder<'a, Self>, FiniteJ2Error> {
                self.validate_history(accepted)?;
                Ok(FrozenFiniteInelasticRemainder {
                    law: self,
                    accepted,
                    backbone: $backbone(self),
                    native_accepted: &accepted.native,
                    prepare_return: |law, f, old| {
                        law.return_map_selected(f, old, false)
                            .map(|(returned, _)| returned)
                    },
                })
            }

            fn remainder_response(
                &self,
                f: &[f64; 9],
                old: &Self::State,
            ) -> Result<FiniteJ2Response<Self::State>, FiniteJ2Error> {
                let (returned, resistance_strain) = self.return_map_selected(f, old, false)?;
                $backbone(self)
                    .response_from_return(f, &old.native, returned)
                    .map(|response| {
                        response.map_state(|native| FiniteRecoveringCreepState {
                            native,
                            resistance_strain,
                        })
                    })
            }

            fn stored_remainder_energy(
                &self,
                f: &[f64; 9],
                state: &Self::State,
            ) -> Result<f64, FiniteJ2Error> {
                self.validate_history(state)?;
                stored($backbone(self), f, &state.native)
            }

            fn validate_remainder_committed(
                &self,
                f: &[f64; 9],
                state: &Self::State,
                yield_tolerance_pa: f64,
            ) -> Result<(), FiniteJ2Error> {
                self.validate_history(state)?;
                validate(
                    $backbone(self),
                    f,
                    &state.native,
                    yield_tolerance_pa,
                    $equilibrium_yield,
                )
            }
        }
    };
}

recovering_kinetic!(
    FiniteRecoveringNortonCreep,
    recovering_norton_backbone,
    false
);
recovering_kinetic!(FiniteRecoveringNortonJ2, recovering_combined_backbone, true);
