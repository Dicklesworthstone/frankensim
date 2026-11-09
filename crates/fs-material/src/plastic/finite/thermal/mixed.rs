//! Original-volume thermal inelasticity with independent mixed pressure.
//!
//! These laws use Fm=F/theta(T), with ALL energy per original volume. For
//! q=theta^3, pressure conjugate to total J maps to native pm=q*p. The exact
//! dual is C(pm)+p*(q-1), with reconstructed J=q*Jm(pm), compliance q²*C''
//! and unscaled physical bulk storage U(Jm). A relaxed-volume reference
//! wrapper, which multiplies energy by q, would describe a different law.
//!
//! Original isochoric return maps, histories, hardening and kinetic clocks
//! supply the remainder. Bulk is omitted at its original assembly sites;
//! no large penalty is evaluated and subtracted. Inspection uses the accepted
//! temperature and performs no flow or recovery. Original bulk branches,
//! total/mechanical deformation bounds and typed-history checks still apply.
//!
//! Mixed entropy is the negative temperature derivative of the fixed-history
//! FREE-ENERGY saddle, not the incremental flow potential. Its bulk part is
//! 3*a*p*J_reconstructed-(K'/K)*U. Exact elastic incompressibility enforces
//! J=theta³ and has zero bulk storage, but retains 3*a*p*theta³ entropy.
//! Physical internal storage is free energy plus T*s. Reversible caloric
//! transfer, native inelastic heat, weak-volume work and signed endpoint
//! defects stay distinct. Supplied temperatures are not a heat solve; no
//! fixed-F heat-capacity bound is relabeled a mixed fixed-pressure derivative.

use super::recovery::{FrozenRecovery, FrozenThermalRecoveringCreep};
use super::{
    AffineThermalFiniteJ2, FiniteJ2Response, FiniteJ2State, FiniteJ2ThermalError,
    FiniteThermalInelasticMaterial, FiniteThermalInelasticPoint, FrozenNative,
    FrozenThermalFiniteJ2, finite,
};
use crate::MaterialError;
use crate::hyper::Volumetric;
use crate::hyper::mixed::{
    MixedFiniteStrainPotential, MixedPointResponse, MixedPointState, PressurePotential,
    RemainderResponse, VolumetricTreatment, cofactor, combine_response, combine_state,
    pressure_potential,
};
use crate::plastic::finite::creep::recovery::FiniteRecoveringCreepState;
use crate::plastic::finite::{FiniteJ2Error, ReturnMap, dot, log_one_plus};
use fs_math::det;

/// Fixed-history storage and entropy on the selected mixed pressure chart.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MixedThermalThermodynamics {
    pub first_piola_pa: [f64; 9],
    pub free_energy_j_m3: f64,
    pub entropy_j_m3_k: f64,
    pub internal_energy_j_m3: f64,
    pub bulk_free_energy_j_m3: f64,
    pub bulk_entropy_j_m3_k: f64,
    /// Actual total det(F), which need not satisfy volume pointwise.
    pub jacobian: f64,
    /// Total volume reconstructed from pressure and thermal stretch.
    pub bulk_volume_ratio: f64,
}

/// Native returned memory and independent mixed thermal work accounts.
#[derive(Clone, Debug)]
pub struct MixedThermalIncrement<S = FiniteJ2State> {
    /// BULK-FREE native response mapped to total F; state and physical heat
    /// retain their original meaning. Full pressure stress is in `mixed`.
    pub mechanical: FiniteJ2Response<S>,
    pub mixed: MixedPointState,
    pub old: MixedThermalThermodynamics,
    pub new: MixedThermalThermodynamics,
    pub coefficient_storage_change_j_m3: f64,
    pub thermoelastic_storage_change_j_m3: f64,
    /// Complete endpoint Piola work, including actual pressure work.
    pub mechanical_work_j_m3: f64,
    /// Remainder coordinate work plus bulk work excluding thermal-volume work.
    pub mechanical_coordinate_work_j_m3: f64,
    pub thermal_deformation_work_j_m3: f64,
    /// Signed endpoint and weak-volume defect; never material heat.
    pub endpoint_work_defect_j_m3: f64,
    pub reversible_heat_j_m3: f64,
    pub balance_error_j_m3: f64,
}

impl<S> MixedThermalIncrement<S> {
    pub fn caloric_heat_j_m3(&self) -> f64 {
        self.mechanical.inelastic_dissipation_j_m3() + self.reversible_heat_j_m3
    }

    pub fn internal_energy_change_j_m3(&self) -> f64 {
        self.new.free_energy_j_m3 - self.old.free_energy_j_m3
            + self.thermoelastic_storage_change_j_m3
    }
}

/// Complete endpoint-temperature point with an original-volume mixed split.
pub trait MixedThermalInelasticPoint: FiniteThermalInelasticPoint {
    fn freeze_mixed<'a>(
        &'a self,
        accepted: &'a Self::State,
    ) -> Result<FrozenMixedThermalInelastic<'a, Self>, FiniteJ2ThermalError>;

    fn mixed_thermodynamics(
        &self,
        f: &[f64; 9],
        state: &Self::State,
        pressure_pa: f64,
        treatment: VolumetricTreatment,
    ) -> Result<MixedThermalThermodynamics, FiniteJ2ThermalError>;

    fn validate_mixed_committed(
        &self,
        f: &[f64; 9],
        state: &Self::State,
        pressure_pa: f64,
        treatment: VolumetricTreatment,
        yield_tolerance_pa: f64,
    ) -> Result<(), FiniteJ2ThermalError>;

    /// `self` is the new frozen point; both views must belong to the same
    /// original family. The old view may be inspection-only; self must retain
    /// the physical duration whenever the selected law has kinetics.
    #[allow(clippy::too_many_arguments)]
    fn checked_mixed_increment_from(
        &self,
        old_point: &Self,
        old_f: &[f64; 9],
        old_state: &Self::State,
        old_pressure_pa: f64,
        new_f: &[f64; 9],
        new_pressure_pa: f64,
        treatment: VolumetricTreatment,
        yield_tolerance_pa: f64,
    ) -> Result<MixedThermalIncrement<Self::State>, FiniteJ2ThermalError>;
}

/// Borrowed endpoint law and complete accepted history for native mixed Newton.
pub struct FrozenMixedThermalInelastic<'a, P: FiniteThermalInelasticPoint + ?Sized> {
    point: &'a P,
    accepted: &'a P::State,
    thermal: &'a FrozenThermalFiniteJ2,
    native_accepted: &'a FiniteJ2State,
    prepare_return: fn(&P, &[f64; 9], &P::State) -> Result<ReturnMap, FiniteJ2ThermalError>,
}

impl<P: FiniteThermalInelasticPoint + ?Sized> core::fmt::Debug
    for FrozenMixedThermalInelastic<'_, P>
{
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("FrozenMixedThermalInelastic")
            .field("temperature_k", &self.thermal.temperature_k)
            .field("accepted", self.accepted)
            .finish_non_exhaustive()
    }
}

fn material_error(error: FiniteJ2ThermalError) -> MaterialError {
    MaterialError::State {
        what: format!("mixed finite thermal material: {error}"),
    }
}

struct Bulk {
    mixed: PressurePotential,
    mechanical: PressurePotential,
    thermal_volume: f64,
}

fn bulk(
    thermal: &FrozenThermalFiniteJ2,
    actual_log_volume: f64,
    elastic_deviator: &[f64; 9],
    pressure_pa: f64,
    treatment: VolumetricTreatment,
) -> Result<Bulk, FiniteJ2ThermalError> {
    if !pressure_pa.is_finite() {
        return Err(FiniteJ2Error::State("mixed thermal pressure must be finite").into());
    }
    let thermal_volume = finite(det::exp(3.0 * thermal.log_thermal_stretch))?;
    let thermal_increment = finite(det::expm1(3.0 * thermal.log_thermal_stretch))?;
    if thermal_volume <= 0.0 {
        return Err(FiniteJ2ThermalError::Numerical(
            "positive thermal volume required",
        ));
    }
    let backbone = thermal.backbone();
    let mechanical = match treatment {
        VolumetricTreatment::ExactIncompressible => PressurePotential {
            complementary_pa: 0.0,
            bulk_stored_pa: 0.0,
            volume_ratio: 1.0,
            volume_increment: 0.0,
            compliance_per_pa: 0.0,
        },
        VolumetricTreatment::OriginalMaterial => {
            if actual_log_volume >= 1.0 {
                return Err(FiniteJ2Error::Domain {
                    quantity: "mixed thermal actual ln(Jm), local convex branch",
                    value: actual_log_volume,
                    limit: 1.0,
                }
                .into());
            }
            pressure_potential(
                Volumetric::LogSquared {
                    kappa: backbone.parameters.bulk_modulus_pa,
                },
                finite(thermal_volume * pressure_pa)?,
            )
            .map_err(FiniteJ2ThermalError::MixedPressure)?
        }
    };
    if mechanical.volume_ratio <= backbone.limits.minimum_jacobian {
        return Err(FiniteJ2Error::Domain {
            quantity: "mixed thermal pressure-reconstructed Jm, strict lower bound",
            value: mechanical.volume_ratio,
            limit: backbone.limits.minimum_jacobian,
        }
        .into());
    }
    backbone.check_elastic_domain(log_one_plus(mechanical.volume_increment), elastic_deviator)?;
    let mixed = PressurePotential {
        complementary_pa: finite(mechanical.complementary_pa + pressure_pa * thermal_increment)?,
        bulk_stored_pa: mechanical.bulk_stored_pa,
        volume_ratio: finite(thermal_volume * mechanical.volume_ratio)?,
        volume_increment: finite(thermal_volume * mechanical.volume_increment + thermal_increment)?,
        compliance_per_pa: finite(
            thermal_volume * (thermal_volume * mechanical.compliance_per_pa),
        )?,
    };
    if mixed.volume_ratio <= 0.0
        || (treatment == VolumetricTreatment::OriginalMaterial && mixed.compliance_per_pa <= 0.0)
    {
        return Err(FiniteJ2ThermalError::Numerical(
            "mixed thermal volume or compliance underflow",
        ));
    }
    Ok(Bulk {
        mixed,
        mechanical,
        thermal_volume,
    })
}

fn inspect(
    thermal: &FrozenThermalFiniteJ2,
    f: &[f64; 9],
    state: &FiniteJ2State,
    pressure_pa: f64,
    treatment: VolumetricTreatment,
) -> Result<(MixedThermalThermodynamics, Bulk), FiniteJ2ThermalError> {
    let remainder = thermal.inspect_selected(f, state, false)?;
    let bulk = bulk(
        thermal,
        remainder.log_volume,
        &remainder.deviator,
        pressure_pa,
        treatment,
    )?;
    let coefficient_entropy = if treatment == VolumetricTreatment::OriginalMaterial {
        finite(
            (thermal.slopes.bulk_modulus_pa_k / thermal.backbone().parameters.bulk_modulus_pa)
                * bulk.mixed.bulk_stored_pa,
        )?
    } else {
        0.0
    };
    let bulk_entropy = finite(
        3.0 * thermal.expansion_per_k * pressure_pa * bulk.mixed.volume_ratio - coefficient_entropy,
    )?;
    let native = remainder.thermodynamics;
    let free = finite(native.free_energy_j_m3 + bulk.mixed.bulk_stored_pa)?;
    let entropy = finite(native.entropy_j_m3_k + bulk_entropy)?;
    let cof = cofactor(f);
    let first_piola_pa = core::array::from_fn(|i| native.first_piola_pa[i] + pressure_pa * cof[i]);
    if first_piola_pa.iter().any(|value| !value.is_finite()) {
        return Err(FiniteJ2ThermalError::Numerical(
            "mixed thermal Piola is not representable",
        ));
    }
    Ok((
        MixedThermalThermodynamics {
            first_piola_pa,
            free_energy_j_m3: free,
            entropy_j_m3_k: entropy,
            internal_energy_j_m3: finite(thermal.temperature_k.mul_add(entropy, free))?,
            bulk_free_energy_j_m3: bulk.mixed.bulk_stored_pa,
            bulk_entropy_j_m3_k: bulk_entropy,
            jacobian: native.jacobian,
            bulk_volume_ratio: bulk.mixed.volume_ratio,
        },
        bulk,
    ))
}

impl<P: FiniteThermalInelasticPoint + ?Sized> FrozenMixedThermalInelastic<'_, P> {
    fn prepare(&self, f: &[f64; 9]) -> Result<([f64; 9], ReturnMap), FiniteJ2ThermalError> {
        let (fm, _) = self.thermal.mechanical_deformation(f)?;
        Ok((fm, (self.prepare_return)(self.point, &fm, self.accepted)?))
    }

    fn pressure(
        &self,
        returned: &ReturnMap,
        p: f64,
        treatment: VolumetricTreatment,
    ) -> Result<PressurePotential, FiniteJ2ThermalError> {
        bulk(
            self.thermal,
            returned.trial.log_j,
            &returned.trial.deviator.map(|value| value * returned.scale),
            p,
            treatment,
        )
        .map(|value| value.mixed)
    }

    fn remainder(
        &self,
        fm: &[f64; 9],
        returned: &ReturnMap,
    ) -> Result<RemainderResponse, FiniteJ2ThermalError> {
        let scale = self.thermal.inverse_thermal_stretch;
        let tangent_scale = scale * scale;
        let tangent = self
            .thermal
            .backbone()
            .return_tangent(fm, self.native_accepted, returned)?;
        let result = RemainderResponse {
            potential_pa: returned.potential,
            first_piola_pa: returned.first_piola.map(|value| value * scale),
            tangent_pa: tangent.map(|row| row.map(|value| value * tangent_scale)),
        };
        if result
            .first_piola_pa
            .iter()
            .chain(result.tangent_pa.iter().flatten())
            .any(|v| !v.is_finite())
        {
            return Err(FiniteJ2ThermalError::Numerical(
                "mixed thermal tangent pullback is not representable",
            ));
        }
        Ok(result)
    }
}

impl<P: FiniteThermalInelasticPoint + ?Sized> MixedFiniteStrainPotential
    for FrozenMixedThermalInelastic<'_, P>
{
    fn validate_mixed(&self) -> Result<(), MaterialError> {
        self.thermal
            .validate_history(self.native_accepted)
            .map_err(material_error)
    }

    fn remainder_stress_and_potential(
        &self,
        f: &[f64; 9],
    ) -> Result<([f64; 9], f64), MaterialError> {
        let (_, returned) = self.prepare(f).map_err(material_error)?;
        let piola = returned
            .first_piola
            .map(|value| value * self.thermal.inverse_thermal_stretch);
        if piola.iter().any(|value| !value.is_finite()) {
            return Err(material_error(FiniteJ2ThermalError::Numerical(
                "mixed thermal Piola pullback is not representable",
            )));
        }
        Ok((piola, returned.potential))
    }

    fn remainder_response(&self, f: &[f64; 9]) -> Result<RemainderResponse, MaterialError> {
        let (fm, returned) = self.prepare(f).map_err(material_error)?;
        self.remainder(&fm, &returned).map_err(material_error)
    }

    fn mixed_pressure_potential(
        &self,
        f: &[f64; 9],
        p: f64,
        treatment: VolumetricTreatment,
    ) -> Result<PressurePotential, MaterialError> {
        let (_, returned) = self.prepare(f).map_err(material_error)?;
        self.pressure(&returned, p, treatment)
            .map_err(material_error)
    }

    fn mixed_state(
        &self,
        f: &[f64; 9],
        p: f64,
        treatment: VolumetricTreatment,
    ) -> Result<MixedPointState, MaterialError> {
        let (_, returned) = self.prepare(f).map_err(material_error)?;
        let bulk = self
            .pressure(&returned, p, treatment)
            .map_err(material_error)?;
        combine_state(
            f,
            p,
            (
                returned
                    .first_piola
                    .map(|value| value * self.thermal.inverse_thermal_stretch),
                returned.potential,
            ),
            bulk,
        )
    }

    fn mixed_response(
        &self,
        f: &[f64; 9],
        p: f64,
        treatment: VolumetricTreatment,
    ) -> Result<MixedPointResponse, MaterialError> {
        let (fm, returned) = self.prepare(f).map_err(material_error)?;
        let bulk = self
            .pressure(&returned, p, treatment)
            .map_err(material_error)?;
        combine_response(
            f,
            p,
            self.remainder(&fm, &returned).map_err(material_error)?,
            bulk,
        )
    }
}

fn native_return(
    point: &FrozenThermalFiniteJ2,
    fm: &[f64; 9],
    old: &FiniteJ2State,
) -> Result<ReturnMap, FiniteJ2ThermalError> {
    Ok(match &point.native {
        FrozenNative::RateIndependent(law) => law.return_map_selected(fm, old, None, false)?,
        FrozenNative::Overstress(law) => law.return_map_selected(fm, old, false)?,
        FrozenNative::Combined(law) => law.return_map_selected(fm, old, false)?,
        FrozenNative::Inspection(_) | FrozenNative::CombinedInspection(_) => {
            return Err(FiniteJ2ThermalError::DurationRequired);
        }
    })
}

fn recovery_return(
    point: &FrozenThermalRecoveringCreep,
    fm: &[f64; 9],
    old: &FiniteRecoveringCreepState,
) -> Result<(ReturnMap, f64), FiniteJ2ThermalError> {
    Ok(
        match point
            .native
            .as_ref()
            .ok_or(FiniteJ2ThermalError::DurationRequired)?
        {
            FrozenRecovery::Standalone(law) => law.return_map_selected(fm, old, false)?,
            FrozenRecovery::Combined(law) => law.return_map_selected(fm, old, false)?,
        },
    )
}

#[allow(clippy::too_many_arguments)]
fn increment<S>(
    old_point: &FrozenThermalFiniteJ2,
    new_point: &FrozenThermalFiniteJ2,
    old_f: &[f64; 9],
    old_state: &FiniteJ2State,
    old_pressure_pa: f64,
    new_f: &[f64; 9],
    new_pressure_pa: f64,
    treatment: VolumetricTreatment,
    yield_tolerance_pa: f64,
    response: impl FnOnce(&[f64; 9]) -> Result<FiniteJ2Response<S>, FiniteJ2ThermalError>,
    native_state: fn(&S) -> &FiniteJ2State,
) -> Result<MixedThermalIncrement<S>, FiniteJ2ThermalError> {
    let remainder = AffineThermalFiniteJ2::increment_with_response_selected(
        old_f,
        old_state,
        new_f,
        yield_tolerance_pa,
        old_point,
        new_point,
        response,
        native_state,
        false,
    )?;
    let (old, old_bulk) = inspect(old_point, old_f, old_state, old_pressure_pa, treatment)?;
    let (new, new_bulk) = inspect(
        new_point,
        new_f,
        native_state(&remainder.mechanical.state),
        new_pressure_pa,
        treatment,
    )?;
    let mixed = combine_state(
        new_f,
        new_pressure_pa,
        (
            remainder.mechanical.first_piola_pa,
            remainder.mechanical.incremental_potential_j_m3,
        ),
        new_bulk.mixed,
    )
    .map_err(FiniteJ2ThermalError::MixedPressure)?;
    let bulk_coefficient = if treatment == VolumetricTreatment::OriginalMaterial {
        let old_log_j = log_one_plus(old_bulk.mechanical.volume_increment);
        finite(
            0.5 * (new_point.backbone().parameters.bulk_modulus_pa
                - old_point.backbone().parameters.bulk_modulus_pa)
                * old_log_j
                * old_log_j,
        )?
    } else {
        0.0
    };
    let delta_temperature = new_point.temperature_k - old_point.temperature_k;
    let delta_log_theta = finite(
        delta_temperature * (0.5 * old_point.expansion_per_k + 0.5 * new_point.expansion_per_k),
    )?;
    let delta_thermal_volume = finite(old_bulk.thermal_volume * det::expm1(3.0 * delta_log_theta))?;
    let bulk_thermal_work =
        finite(new_pressure_pa * delta_thermal_volume * old_bulk.mechanical.volume_ratio)?;
    let delta_f = core::array::from_fn(|i| new_f[i] - old_f[i]);
    let bulk_work = finite(new_pressure_pa * dot(&cofactor(new_f), &delta_f))?;
    let bulk_ts = finite(delta_temperature.mul_add(
        old.bulk_entropy_j_m3_k,
        new_point.temperature_k * (new.bulk_entropy_j_m3_k - old.bulk_entropy_j_m3_k),
    ))?;
    let bulk_defect = finite(
        bulk_work
            - bulk_thermal_work
            - (new.bulk_free_energy_j_m3 - old.bulk_free_energy_j_m3 - bulk_coefficient),
    )?;
    let coefficient_storage_change_j_m3 =
        finite(remainder.coefficient_storage_change_j_m3 + bulk_coefficient)?;
    let thermoelastic_storage_change_j_m3 =
        finite(remainder.thermoelastic_storage_change_j_m3 + bulk_ts)?;
    let mechanical_work_j_m3 = finite(remainder.mechanical_work_j_m3 + bulk_work)?;
    let mechanical_coordinate_work_j_m3 =
        finite(remainder.mechanical_coordinate_work_j_m3 + bulk_work - bulk_thermal_work)?;
    let thermal_deformation_work_j_m3 =
        finite(remainder.thermal_deformation_work_j_m3 + bulk_thermal_work)?;
    let endpoint_work_defect_j_m3 = finite(remainder.endpoint_work_defect_j_m3 + bulk_defect)?;
    let reversible_heat_j_m3 =
        finite(remainder.reversible_heat_j_m3 + bulk_thermal_work - bulk_coefficient - bulk_ts)?;
    let balance_error_j_m3 = finite(
        mechanical_work_j_m3
            - (new.free_energy_j_m3 - old.free_energy_j_m3 + thermoelastic_storage_change_j_m3)
            - remainder.mechanical.inelastic_dissipation_j_m3()
            - reversible_heat_j_m3
            - endpoint_work_defect_j_m3,
    )?;
    finite(remainder.mechanical.inelastic_dissipation_j_m3() + reversible_heat_j_m3)?;
    Ok(MixedThermalIncrement {
        mechanical: remainder.mechanical,
        mixed,
        old,
        new,
        coefficient_storage_change_j_m3,
        thermoelastic_storage_change_j_m3,
        mechanical_work_j_m3,
        mechanical_coordinate_work_j_m3,
        thermal_deformation_work_j_m3,
        endpoint_work_defect_j_m3,
        reversible_heat_j_m3,
        balance_error_j_m3,
    })
}

impl MixedThermalInelasticPoint for FrozenThermalFiniteJ2 {
    fn freeze_mixed<'a>(
        &'a self,
        accepted: &'a Self::State,
    ) -> Result<FrozenMixedThermalInelastic<'a, Self>, FiniteJ2ThermalError> {
        self.validate_history(accepted)?;
        if matches!(
            self.native,
            FrozenNative::Inspection(_) | FrozenNative::CombinedInspection(_)
        ) {
            return Err(FiniteJ2ThermalError::DurationRequired);
        }
        Ok(FrozenMixedThermalInelastic {
            point: self,
            accepted,
            thermal: self,
            native_accepted: accepted,
            prepare_return: native_return,
        })
    }

    fn mixed_thermodynamics(
        &self,
        f: &[f64; 9],
        state: &Self::State,
        p: f64,
        treatment: VolumetricTreatment,
    ) -> Result<MixedThermalThermodynamics, FiniteJ2ThermalError> {
        inspect(self, f, state, p, treatment).map(|value| value.0)
    }

    fn validate_mixed_committed(
        &self,
        f: &[f64; 9],
        state: &Self::State,
        p: f64,
        treatment: VolumetricTreatment,
        tolerance: f64,
    ) -> Result<(), FiniteJ2ThermalError> {
        self.validate_committed_selected(f, state, tolerance, false)?;
        self.mixed_thermodynamics(f, state, p, treatment)?;
        Ok(())
    }

    fn checked_mixed_increment_from(
        &self,
        old_point: &Self,
        old_f: &[f64; 9],
        old_state: &Self::State,
        old_p: f64,
        new_f: &[f64; 9],
        new_p: f64,
        treatment: VolumetricTreatment,
        tolerance: f64,
    ) -> Result<MixedThermalIncrement<Self::State>, FiniteJ2ThermalError> {
        increment(
            old_point,
            self,
            old_f,
            old_state,
            old_p,
            new_f,
            new_p,
            treatment,
            tolerance,
            |fm| {
                let returned = native_return(self, fm, old_state)?;
                Ok(self
                    .backbone()
                    .response_from_return(fm, old_state, returned)?)
            },
            |state| state,
        )
    }
}

impl MixedThermalInelasticPoint for FrozenThermalRecoveringCreep {
    fn freeze_mixed<'a>(
        &'a self,
        accepted: &'a Self::State,
    ) -> Result<FrozenMixedThermalInelastic<'a, Self>, FiniteJ2ThermalError> {
        self.validate_history(accepted)?;
        if self.native.is_none() {
            return Err(FiniteJ2ThermalError::DurationRequired);
        }
        Ok(FrozenMixedThermalInelastic {
            point: self,
            accepted,
            thermal: &self.thermal,
            native_accepted: &accepted.native,
            prepare_return: |point, f, old| recovery_return(point, f, old).map(|value| value.0),
        })
    }

    fn mixed_thermodynamics(
        &self,
        f: &[f64; 9],
        state: &Self::State,
        p: f64,
        treatment: VolumetricTreatment,
    ) -> Result<MixedThermalThermodynamics, FiniteJ2ThermalError> {
        self.validate_history(state)?;
        inspect(&self.thermal, f, &state.native, p, treatment).map(|value| value.0)
    }

    fn validate_mixed_committed(
        &self,
        f: &[f64; 9],
        state: &Self::State,
        p: f64,
        treatment: VolumetricTreatment,
        tolerance: f64,
    ) -> Result<(), FiniteJ2ThermalError> {
        self.validate_history(state)?;
        self.thermal
            .validate_committed_selected(f, &state.native, tolerance, false)?;
        self.mixed_thermodynamics(f, state, p, treatment)?;
        Ok(())
    }

    fn checked_mixed_increment_from(
        &self,
        old_point: &Self,
        old_f: &[f64; 9],
        old_state: &Self::State,
        old_p: f64,
        new_f: &[f64; 9],
        new_p: f64,
        treatment: VolumetricTreatment,
        tolerance: f64,
    ) -> Result<MixedThermalIncrement<Self::State>, FiniteJ2ThermalError> {
        old_point.validate_history(old_state)?;
        let result = increment(
            &old_point.thermal,
            &self.thermal,
            old_f,
            &old_state.native,
            old_p,
            new_f,
            new_p,
            treatment,
            tolerance,
            |fm| {
                let (returned, resistance_strain) = recovery_return(self, fm, old_state)?;
                Ok(self
                    .thermal
                    .backbone()
                    .response_from_return(fm, &old_state.native, returned)?
                    .map_state(|native| FiniteRecoveringCreepState {
                        native,
                        resistance_strain,
                    }))
            },
            |state| &state.native,
        )?;
        self.validate_history(&result.mechanical.state)?;
        Ok(result)
    }
}

/// Reuse cached temperature views across quadrature points. The old view owns
/// old-temperature history admission; the new view owns one physical return.
#[allow(clippy::too_many_arguments)]
pub fn checked_mixed_increment<P: MixedThermalInelasticPoint>(
    old_point: &P,
    new_point: &P,
    old_f: &[f64; 9],
    old_state: &P::State,
    old_pressure_pa: f64,
    new_f: &[f64; 9],
    new_pressure_pa: f64,
    treatment: VolumetricTreatment,
    yield_tolerance_pa: f64,
) -> Result<MixedThermalIncrement<P::State>, FiniteJ2ThermalError> {
    new_point.checked_mixed_increment_from(
        old_point,
        old_f,
        old_state,
        old_pressure_pa,
        new_f,
        new_pressure_pa,
        treatment,
        yield_tolerance_pa,
    )
}

/// Source-preserving family convenience. Spatial owners should cache the two
/// temperature views and use [`checked_mixed_increment`] for each quadrature point.
#[allow(clippy::too_many_arguments)]
pub fn checked_mixed_increment_for_step<L>(
    family: &L,
    old_f: &[f64; 9],
    old_state: &L::State,
    old_temperature_k: f64,
    old_pressure_pa: f64,
    new_f: &[f64; 9],
    new_temperature_k: f64,
    new_pressure_pa: f64,
    treatment: VolumetricTreatment,
    yield_tolerance_pa: f64,
    duration_s: f64,
) -> Result<MixedThermalIncrement<L::State>, FiniteJ2ThermalError>
where
    L: FiniteThermalInelasticMaterial + ?Sized,
    L::Frozen: MixedThermalInelasticPoint,
{
    let old = family.at_temperature_for_inspection(old_temperature_k)?;
    let new = family.at_temperature_for_step(new_temperature_k, duration_s)?;
    checked_mixed_increment(
        &old,
        &new,
        old_f,
        old_state,
        old_pressure_pa,
        new_f,
        new_pressure_pa,
        treatment,
        yield_tolerance_pa,
    )
}
