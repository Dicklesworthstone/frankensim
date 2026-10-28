//! Equilibrium solid-liquid state resolved on a specific-enthalpy coordinate.
//!
//! Total enthalpy is the primary coordinate because it remains single valued
//! through an isothermal phase-change plateau. A caller supplies a monotone,
//! evidence-derived curve; this module neither recognizes material names nor
//! invents melting points, latent heats, densities, or phase fractions.
//!
//! This is a constitutive state primitive, not a heat-transfer or free-surface
//! solver. Thermal transport owns changes in specific enthalpy. Solid,
//! finite-strain, fluid, remeshing, acoustic, and optical consumers own their
//! respective responses to the returned phase state.

/// Declared stationary solid/liquid mixture mobility on this same h chart.
pub mod mushy;
/// Shared declared phase dependence of thermal conductivity.
pub mod conductivity;
/// Phase-first thermal/electrical sources and optional mobile fluid transport.
pub mod transport;
/// Source-backed solid/liquid Cp branches joined by one sharp fusion plateau.
pub mod melting;
/// Source-backed liquid/vapor coexistence and isobaric heat with expansion work.
pub mod vaporization;
/// Fixed-mass calorimetry with explicit constant-pressure boundary work.
pub mod isobaric {
    //! Fixed-mass calorimetry at a caller-declared constant absolute pressure.
    //!
    //! This is a quasistatic, spatially uniform reference with pressure-volume work
    //! as its only mechanical work. It uses the original equilibrium phase chart;
    //! it does not supply heat transfer, pressure evolution or moving geometry.

    use super::{EquilibriumEnthalpyPhaseCurve, EquilibriumPhaseState, PhaseStateError};

    /// One proposed isobaric heat step; the caller decides whether to accept it.
    #[derive(Clone, Copy, Debug, PartialEq)]
    pub struct IsobaricHeatStep {
        /// Original accepted state, including material and curve identities.
        pub initial: EquilibriumPhaseState,
        /// Proposed state on the same chart; no accepted state is mutated.
        pub state: EquilibriumPhaseState,
        /// Invariant specimen mass [kg].
        pub mass_kg: f64,
        /// Declared constant absolute pressure [Pa].
        pub pressure_pa: f64,
        /// Signed requested heat into the specimen [J].
        pub supplied_heat_j: f64,
        /// Initial volume from the fixed mass and original state density [m3].
        pub initial_volume_m3: f64,
        /// Final volume from that same mass and proposed state density [m3].
        pub final_volume_m3: f64,
        /// Signed volume change [m3], evaluated using relative density change to
        /// retain small expansion/contraction instead of subtracting large volumes.
        pub volume_change_m3: f64,
        /// Represented enthalpy change m*(h_final-h_initial) [J]. Sensible and
        /// latent energy already belong to this single caloric account.
        pub enthalpy_change_j: f64,
        /// Work delivered by the specimen to its surroundings, p*delta_V [J].
        /// Expansion is positive; contraction returns work to the specimen.
        pub boundary_work_j: f64,
        /// Internal-energy change delta_H - p*delta_V [J]. No absolute enthalpy
        /// datum or separate latent-energy term is added to this difference.
        pub internal_energy_change_j: f64,
        /// Requested heat minus the represented enthalpy change [J]. This exposes
        /// floating-point step quantization, not source or integration uncertainty.
        /// A caller must check it against the caller's energy-accuracy requirement.
        pub heat_residual_j: f64,
    }

    impl EquilibriumEnthalpyPhaseCurve {
        /// Advance a fixed mass by signed heat at a declared constant pressure.
        ///
        /// For this closed, quasistatic reference, heat equals the change of
        /// enthalpy. The original chart supplies phase and density, including its
        /// additive-specific-volume rule through coexistence. Pressure-volume work
        /// is positive outwards, so delta_U = delta_H - p*delta_V. Cooling follows
        /// the same admitted equilibrium chart; this does not add hysteresis,
        /// nucleation, stress locking or a resolidification mechanics law.
        ///
        /// The caller must supply the absolute pressure at which BOTH caloric and
        /// density data apply. The generic chart does not discover or verify that
        /// pressure from material names, query axes or unknown-condition flags.
        /// Source-backed callers retain and check their original source receipts.
        /// No kinetic, potential, elastic, chemical or other work account is added.
        ///
        /// Arithmetic uses enthalpy differences before multiplying by mass, avoiding
        /// cancellation or overflow from large arbitrary enthalpy datums. Relative
        /// density differences preserve small pressure work. The reported volumes
        /// are separately rounded; their direct subtraction can lose that change.
        /// A nonzero heat step completely lost at the chosen enthalpy datum refuses.
        /// Other representational heat residuals are returned explicitly for the
        /// caller's accuracy decision. This method neither chooses a time step nor
        /// claims a heat-transfer or material-uncertainty bound.
        ///
        /// No extrapolation, state mutation or allocation is performed. Invalid
        /// mass/pressure/heat, unrepresentable required products and nonfinite
        /// energy differences refuse before a complete proposal is returned.
        pub fn advance_isobaric_heat(
            &self,
            current: EquilibriumPhaseState,
            mass_kg: f64,
            pressure_pa: f64,
            supplied_heat_j: f64,
        ) -> Result<IsobaricHeatStep, PhaseStateError> {
            let invalid = |what| PhaseStateError::InvalidIsobaricHeatStep { what };
            if !mass_kg.is_finite() || mass_kg <= 0.0 {
                return Err(invalid("mass must be finite and positive"));
            }
            if !pressure_pa.is_finite() || pressure_pa < 0.0 {
                return Err(invalid("absolute pressure must be finite and nonnegative"));
            }
            if !supplied_heat_j.is_finite() {
                return Err(invalid("supplied heat must be finite"));
            }
            let dh_requested = supplied_heat_j / mass_kg;
            if !dh_requested.is_finite() || (supplied_heat_j != 0.0 && dh_requested == 0.0) {
                return Err(invalid("specific heat increment is not representable"));
            }
            let state = self.advance_specific_energy(current, dh_requested)?;
            let dh = state.specific_enthalpy_j_kg() - current.specific_enthalpy_j_kg();
            if !dh.is_finite() || (supplied_heat_j != 0.0 && dh == 0.0) {
                return Err(invalid(
                    "heat increment is lost at the chosen enthalpy datum",
                ));
            }
            let enthalpy_change_j = mass_kg * dh;
            if !enthalpy_change_j.is_finite() || (dh != 0.0 && enthalpy_change_j == 0.0) {
                return Err(invalid("represented enthalpy change is not representable"));
            }
            let initial_density = current.bulk_density_kg_m3();
            let final_density = state.bulk_density_kg_m3();
            let initial_volume_m3 = mass_kg / initial_density;
            let final_volume_m3 = mass_kg / final_density;
            if !initial_volume_m3.is_finite()
                || initial_volume_m3 <= 0.0
                || !final_volume_m3.is_finite()
                || final_volume_m3 <= 0.0
            {
                return Err(invalid("fixed-mass volumes must be finite and positive"));
            }
            // Delta V = max(V0,V1)*(rho0-rho1)/max(rho0,rho1). The relative
            // factor is bounded by one and avoids subtracting nearly equal volumes.
            let density_change = initial_density - final_density;
            let relative_change = density_change / initial_density.max(final_density);
            let volume_change_m3 = initial_volume_m3.max(final_volume_m3) * relative_change;
            if !volume_change_m3.is_finite() || (density_change != 0.0 && volume_change_m3 == 0.0) {
                return Err(invalid("volume change is not representable"));
            }
            let boundary_work_j = pressure_pa * volume_change_m3;
            if !boundary_work_j.is_finite()
                || (pressure_pa != 0.0 && volume_change_m3 != 0.0 && boundary_work_j == 0.0)
            {
                return Err(invalid("pressure-volume work is not representable"));
            }
            let internal_energy_change_j = enthalpy_change_j - boundary_work_j;
            let heat_residual_j = supplied_heat_j - enthalpy_change_j;
            if !internal_energy_change_j.is_finite() || !heat_residual_j.is_finite() {
                return Err(invalid("energy differences are not representable"));
            }
            Ok(IsobaricHeatStep {
                initial: current,
                state,
                mass_kg,
                pressure_pa,
                supplied_heat_j,
                initial_volume_m3,
                final_volume_m3,
                volume_change_m3,
                enthalpy_change_j,
                boundary_work_j,
                internal_energy_change_j,
                heat_residual_j,
            })
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use crate::phase::{EnthalpyPhaseKnot, SolidLiquidPhase};
        use fs_blake3::ContentHash;

        fn near(actual: f64, expected: f64) {
            assert!(
                (actual - expected).abs() <= 2e-14 * expected.abs().max(1.0),
                "{actual} != {expected}"
            );
        }

        fn curve(datum: f64) -> EquilibriumEnthalpyPhaseCurve {
            EquilibriumEnthalpyPhaseCurve::try_new(
                ContentHash([0x61; 32]),
                [
                    (0.0, 300.0, 0.0, 1000.0),
                    (131_072.0, 400.0, 0.0, 950.0),
                    (393_216.0, 400.0, 1.0, 800.0),
                    (655_360.0, 500.0, 1.0, 750.0),
                ]
                .map(
                    |(h, temperature_k, liquid_mass_fraction, bulk_density_kg_m3)| EnthalpyPhaseKnot {
                        specific_enthalpy_j_kg: h + datum,
                        temperature_k,
                        liquid_mass_fraction,
                        bulk_density_kg_m3,
                    },
                )
                .to_vec(),
            )
            .unwrap()
        }

        fn single(h: [f64; 2], rho: [f64; 2]) -> EquilibriumEnthalpyPhaseCurve {
            EquilibriumEnthalpyPhaseCurve::try_single_phase(
                ContentHash([0x62; 32]),
                SolidLiquidPhase::Solid,
                vec![
                    EnthalpyPhaseKnot {
                        specific_enthalpy_j_kg: h[0],
                        temperature_k: 300.0,
                        liquid_mass_fraction: 0.0,
                        bulk_density_kg_m3: rho[0],
                    },
                    EnthalpyPhaseKnot {
                        specific_enthalpy_j_kg: h[1],
                        temperature_k: 400.0,
                        liquid_mass_fraction: 0.0,
                        bulk_density_kg_m3: rho[1],
                    },
                ],
            )
            .unwrap()
        }

        #[test]
        fn g1_latent_heat_and_specific_volume_drive_one_isobaric_energy_account() {
            let chart = curve(0.0);
            let solid = chart.state_at_specific_enthalpy(131_072.0).unwrap();
            let half = chart
                .advance_isobaric_heat(solid, 2.0, 200_000.0, 262_144.0)
                .unwrap();
            assert_eq!(half.state.temperature_k(), 400.0);
            assert_eq!(half.state.liquid_mass_fraction(), 0.5);
            near(half.initial_volume_m3, 2.0 / 950.0);
            near(half.final_volume_m3, 1.0 / 950.0 + 1.0 / 800.0);
            near(half.volume_change_m3, 1.0 / 800.0 - 1.0 / 950.0);
            near(
                half.boundary_work_j,
                200_000.0 * (1.0 / 800.0 - 1.0 / 950.0),
            );
            near(half.enthalpy_change_j, 262_144.0);
            near(
                half.internal_energy_change_j + half.boundary_work_j,
                half.enthalpy_change_j,
            );
            assert_eq!(half.heat_residual_j, 0.0);
            assert_eq!(
                half.state.phase_curve_identity(),
                solid.phase_curve_identity()
            );
            assert_eq!(
                half.state.material_card_identity(),
                solid.material_card_identity()
            );
        }

        #[test]
        fn g3_melt_cool_cycle_preserves_mass_and_reverses_work_and_energy() {
            let chart = curve(0.0);
            let initial = chart.state_at_specific_enthalpy(0.0).unwrap();
            let hot = chart
                .advance_isobaric_heat(initial, 2.0, 100_000.0, 1_310_720.0)
                .unwrap();
            assert_eq!(hot.state.phase(), SolidLiquidPhase::Liquid);
            near(hot.final_volume_m3, 2.0 / 750.0);
            let cold = chart
                .advance_isobaric_heat(hot.state, 2.0, 100_000.0, -1_310_720.0)
                .unwrap();
            assert_eq!(cold.state, initial);
            assert_eq!(cold.enthalpy_change_j, -hot.enthalpy_change_j);
            assert_eq!(cold.volume_change_m3, -hot.volume_change_m3);
            assert_eq!(cold.boundary_work_j, -hot.boundary_work_j);
            assert_eq!(cold.internal_energy_change_j, -hot.internal_energy_change_j);
            let unchanged = chart
                .advance_isobaric_heat(initial, 2.0, 100_000.0, 0.0)
                .unwrap();
            assert_eq!(unchanged.state, initial);
            assert_eq!(unchanged.boundary_work_j, 0.0);
            assert_eq!(unchanged.internal_energy_change_j, 0.0);
            let vacuum = chart
                .advance_isobaric_heat(initial, 2.0, 0.0, 262_144.0)
                .unwrap();
            assert_eq!(vacuum.boundary_work_j, 0.0);
            assert_eq!(vacuum.internal_energy_change_j, vacuum.enthalpy_change_j);
        }

        #[test]
        fn g3_enthalpy_datum_shifts_do_not_change_heat_volume_or_pressure_work() {
            let ordinary = curve(0.0);
            let datum = 2.0_f64.powi(60);
            let shifted = curve(datum);
            let a = ordinary
                .advance_isobaric_heat(
                    ordinary.state_at_specific_enthalpy(0.0).unwrap(),
                    2.0,
                    1e5,
                    524_288.0,
                )
                .unwrap();
            let b = shifted
                .advance_isobaric_heat(
                    shifted.state_at_specific_enthalpy(datum).unwrap(),
                    2.0,
                    1e5,
                    524_288.0,
                )
                .unwrap();
            assert_eq!(a.state.temperature_k(), b.state.temperature_k());
            assert_eq!(
                a.state.liquid_mass_fraction(),
                b.state.liquid_mass_fraction()
            );
            assert_eq!(a.enthalpy_change_j, b.enthalpy_change_j);
            assert_eq!(a.boundary_work_j, b.boundary_work_j);
            assert_eq!(a.internal_energy_change_j, b.internal_energy_change_j);
            let h0 = 1e308;
            let h1 = h0 + 4e292;
            let chart = single([h0, h1], [1000.0; 2]);
            let heat = 2.0 * (h1 - h0);
            assert!(!(2.0 * h0).is_finite());
            let result = chart
                .advance_isobaric_heat(
                    chart.state_at_specific_enthalpy(h0).unwrap(),
                    2.0,
                    1e5,
                    heat,
                )
                .unwrap();
            assert_eq!(result.enthalpy_change_j, heat);
            assert_eq!(result.internal_energy_change_j, heat);
        }

        #[test]
        fn g1_small_density_changes_retain_pressure_work_beyond_rounded_volume_subtraction() {
            let epsilon = f64::EPSILON;
            let chart = single([0.0, 1.0], [1.0, 1.0 + epsilon]);
            let result = chart
                .advance_isobaric_heat(
                    chart.state_at_specific_enthalpy(0.0).unwrap(),
                    1e20,
                    1e5,
                    1e20,
                )
                .unwrap();
            let expected_delta_v = -1e20 * (epsilon / (1.0 + epsilon));
            near(result.volume_change_m3, expected_delta_v);
            near(result.boundary_work_j, 1e5 * expected_delta_v);
            assert!(
                (result.volume_change_m3 - (result.final_volume_m3 - result.initial_volume_m3)).abs()
                    > 5000.0
            );
            assert!(result.boundary_work_j < 0.0);
            assert!(result.internal_energy_change_j > result.enthalpy_change_j);
        }

        #[test]
        fn g0_heat_rounding_is_visible_and_fully_lost_steps_refuse() {
            let datum = 2.0_f64.powi(60);
            let chart = single([datum, datum + 4096.0], [1000.0; 2]);
            let initial = chart.state_at_specific_enthalpy(datum).unwrap();
            let result = chart
                .advance_isobaric_heat(initial, 1.0, 1e5, 300.0)
                .unwrap();
            assert_eq!(result.enthalpy_change_j, 256.0);
            assert_eq!(result.heat_residual_j, 44.0);
            assert_eq!(result.internal_energy_change_j, 256.0);
            assert!(chart.advance_isobaric_heat(initial, 1.0, 1e5, 1.0).is_err());
            assert!(
                chart
                    .advance_isobaric_heat(initial, f64::MAX, 1e5, f64::from_bits(1))
                    .is_err()
            );
        }

        #[test]
        fn g0_invalid_inputs_wrong_curve_and_outside_domain_refuse_transactionally() {
            let chart = curve(0.0);
            let initial = chart.state_at_specific_enthalpy(0.0).unwrap();
            for mass in [0.0, -1.0, f64::INFINITY, f64::NAN] {
                assert!(
                    chart
                        .advance_isobaric_heat(initial, mass, 1e5, 1.0)
                        .is_err()
                );
            }
            for pressure in [-1.0, f64::INFINITY, f64::NAN] {
                assert!(
                    chart
                        .advance_isobaric_heat(initial, 1.0, pressure, 1.0)
                        .is_err()
                );
            }
            for heat in [f64::INFINITY, f64::NAN, -1.0, 700_000.0] {
                assert!(
                    chart
                        .advance_isobaric_heat(initial, 1.0, 1e5, heat)
                        .is_err()
                );
            }
            assert_eq!(
                curve(1.0).advance_isobaric_heat(initial, 1.0, 1e5, 1.0),
                Err(PhaseStateError::StateCurveMismatch)
            );
            assert_eq!(initial, chart.state_at_specific_enthalpy(0.0).unwrap());
        }

        #[test]
        fn g0_unrepresentable_volumes_work_and_internal_energy_refuse() {
            let chart = single([0.0, 1.0], [f64::MIN_POSITIVE; 2]);
            let state = chart.state_at_specific_enthalpy(0.0).unwrap();
            assert!(
                chart
                    .advance_isobaric_heat(state, f64::MAX, 1.0, 0.0)
                    .is_err()
            );
            let chart = single([0.0, 1e-300], [1.0, 0.5]);
            let state = chart.state_at_specific_enthalpy(0.0).unwrap();
            assert!(
                chart
                    .advance_isobaric_heat(state, 1e300, 1e100, 1.0)
                    .is_err()
            );
            let energy = 0.75 * f64::MAX;
            let chart = single([0.0, energy], [0.5, 1.0]);
            let state = chart.state_at_specific_enthalpy(0.0).unwrap();
            assert!(
                chart
                    .advance_isobaric_heat(state, 1.0, energy, energy)
                    .is_err()
            );
        }
    }
}

// Relative equilibrium entropy uses this same admitted h -> T chart.
mod entropy;

use core::fmt;

use fs_blake3::{ContentHash, DomainHasher};
use fs_matdb::{
    ClaimSelection, ClaimSet, MaterialCard, PropertyKey, PropertyUsageReceipt,
    QueryPoint, SelectionPolicy,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, QuantitySpec};

const EQUILIBRIUM_ENTHALPY_PHASE_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-material.equilibrium-enthalpy-phase.v1";
const EQUILIBRIUM_PHASE_STATE_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-material.equilibrium-phase-state.v1";

/// Upper bound on generated enthalpy knots from one heat-capacity chart.
pub const MAX_HEAT_CAPACITY_ENTHALPY_KNOTS: usize = 4_096;

/// Proposed uniform-body thermal step shared by transport and coupled owners.
/// Geometry and mass belong to the specimen; this carrier does not justify
/// lumping, solve transport, or mutate the accepted state.
#[derive(Clone, Copy, Debug)]
pub struct UniformEnthalpyStepInput<'a> {
    /// Immutable equilibrium chart used by the accepted body.
    pub curve: &'a EquilibriumEnthalpyPhaseCurve,
    /// Accepted initial state on that chart.
    pub initial: EquilibriumPhaseState,
    /// Invariant specimen mass [kg].
    pub mass_kg: f64,
    /// Specimen volume [m3].
    pub volume_m3: f64,
    /// Whole exposed boundary area [m2].
    pub surface_area_m2: f64,
    /// Internally deposited heat over this step [J], counted exactly once.
    pub internal_heat_j: f64,
    /// Physical duration of this step [s].
    pub duration_s: f64,
}

/// Proposed thermal state and boundary transfer, with a solver-owned report.
/// The coupled owner must check the chart, energy balance and its own regime
/// before publishing this proposal. A transport callback must not debit a
/// mutable external reservoir before the coupled transaction accepts it.
#[derive(Clone, Debug, PartialEq)]
pub struct UniformEnthalpyStep<R> {
    /// Proposed state on the input equilibrium chart.
    pub state: EquilibriumPhaseState,
    /// Signed integrated boundary heat into the body [J].
    pub external_heat_j: f64,
    /// Declared transport solve residual allowance [J], separate from roundoff
    /// and from temporal or constitutive-model error.
    pub energy_residual_tolerance_j: f64,
    /// Actual transport result retained without erasing its concrete type.
    pub report: R,
}

/// One source-provided point on a solid-liquid equilibrium enthalpy curve.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EnthalpyPhaseKnot {
    /// Specific enthalpy [J/kg]. Knots must be strictly increasing in this
    /// coordinate.
    pub specific_enthalpy_j_kg: f64,
    /// Absolute temperature [K]. It must be positive and nondecreasing; equal
    /// temperatures are deliberately admitted for an isothermal latent-heat
    /// plateau.
    pub temperature_k: f64,
    /// Equilibrium liquid mass fraction in `[0, 1]`, nondecreasing with
    /// enthalpy. The solid mass fraction is exactly its complement.
    pub liquid_mass_fraction: f64,
    /// Equilibrium bulk density [kg/m3] at this state.
    pub bulk_density_kg_m3: f64,
}

/// One source point for a single-phase heat-capacity chart.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HeatCapacityKnot {
    /// Absolute temperature [K].
    pub temperature_k: f64,
    /// Specific heat capacity [J/(kg K)].
    pub specific_heat_capacity_j_kg_k: f64,
    /// Bulk density [kg/m3]. It is linearly interpolated in temperature only
    /// while producing the enthalpy chart; the resulting chart retains its
    /// existing specific-volume interpolation between generated knots.
    pub bulk_density_kg_m3: f64,
}

/// Explicit source keys and complete endpoints for one fixed-context thermal span.
#[derive(Clone, Copy, Debug)]
pub struct HeatCapacitySourceQuery<'a> {
    /// Specific heat capacity [J/(kg K)], not extensive heat capacity.
    pub heat_capacity: &'a PropertyKey,
    /// Bulk density [kg/m3].
    pub density: &'a PropertyKey,
    /// Absolute-temperature coordinate in kelvin.
    pub temperature_axis: &'a str,
    /// Lower-temperature endpoint with every source condition declared.
    pub lower: &'a QueryPoint,
    /// Upper endpoint; all non-temperature coordinates must remain fixed.
    pub upper: &'a QueryPoint,
    /// Source selection policy; sources cannot change within the span.
    pub selection: SelectionPolicy,
}

/// A single-phase enthalpy chart and its retained Cp/density source receipts.
/// Source conditions and receipt identities are bound into the chart identity.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedHeatCapacityCurve {
    curve: EquilibriumEnthalpyPhaseCurve,
    receipts: Vec<PropertyUsageReceipt>,
}

impl ResolvedHeatCapacityCurve {
    /// Resolve continuously supported Cp and density from one immutable card.
    /// Every source breakpoint along the fixed-context temperature line is
    /// retained before integrating piecewise-linear Cp, including crossed
    /// cells of an explicitly supported two-coordinate source surface.
    /// The reference enthalpy is declared at the lower endpoint. Phase is an
    /// explicit caller assumption: this does not discover phase boundaries.
    /// The interpolation-only error limit has the same meaning as
    /// [`EquilibriumEnthalpyPhaseCurve::try_from_heat_capacity`].
    pub fn try_new(
        card: &MaterialCard,
        query: HeatCapacitySourceQuery<'_>,
        phase: SolidLiquidPhase,
        reference_specific_enthalpy_j_kg: f64,
        maximum_interpolation_error_j_kg: f64,
    ) -> Result<Self, PhaseStateError> {
        let invalid = |what| PhaseStateError::InvalidCurve { what };
        single_phase_fraction(phase)?;
        let keys = [query.heat_capacity, query.density];
        for (key, dims) in keys
            .iter()
            .zip([Dims([2, 0, -2, -1, 0, 0]), Dims([-3, 1, 0, 0, 0, 0])])
        {
            if key.quantity() != QuantitySpec::dimensional(dims) {
                return Err(invalid(
                    "source keys must declare specific heat capacity and bulk density SI dimensions",
                ));
            }
        }
        let endpoints = heat_capacity_source_endpoints(query)?;
        let claims = card.claims();
        let temperatures = heat_capacity_source_temperatures(claims, query, endpoints)?;
        let mut knots = Vec::with_capacity(temperatures.len());
        let mut receipts = Vec::with_capacity(2 * temperatures.len());
        for temperature in temperatures {
            let point = query
                .lower
                .clone()
                .with(query.temperature_axis, temperature)
                .map_err(|error| PhaseStateError::SourceQuery {
                    property: query.temperature_axis.to_owned(),
                    reason: error.to_string(),
                })?;
            let mut values = [0.0; 2];
            for (index, key) in keys.iter().enumerate() {
                let answer = claims
                    .query_typed(key, &point, query.selection)
                    .map_err(|error| PhaseStateError::SourceQuery {
                        property: key.name().to_owned(),
                        reason: error.to_string(),
                    })?;
                values[index] = answer.evidence.value.value;
                receipts.push(answer.receipt);
            }
            knots.push(HeatCapacityKnot {
                temperature_k: temperature,
                specific_heat_capacity_j_kg_k: values[0],
                bulk_density_kg_m3: values[1],
            });
        }
        let mut curve = EquilibriumEnthalpyPhaseCurve::try_from_heat_capacity(
            card.content_hash(),
            phase,
            reference_specific_enthalpy_j_kg,
            &knots,
            maximum_interpolation_error_j_kg,
        )?;
        let mut identity =
            DomainHasher::new("org.frankensim.fs-material.sourced-heat-capacity-curve.v1");
        identity.update(curve.identity().as_bytes());
        identity.update(&(receipts.len() as u64).to_le_bytes());
        for receipt in &receipts {
            identity.update(receipt.content_hash().as_bytes());
        }
        curve.identity = identity.finalize();
        Ok(Self { curve, receipts })
    }

    /// The source-bound chart accepted by thermal and phase-state consumers.
    #[must_use]
    pub const fn curve(&self) -> &EquilibriumEnthalpyPhaseCurve {
        &self.curve
    }

    /// Receipts in increasing temperature order, Cp then density at each knot.
    #[must_use]
    pub fn receipts(&self) -> &[PropertyUsageReceipt] {
        &self.receipts
    }
}

fn heat_capacity_source_endpoints(
    query: HeatCapacitySourceQuery<'_>,
) -> Result<[f64; 2], PhaseStateError> {
    let invalid = |what| PhaseStateError::InvalidCurve { what };
    let absolute = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ));
    let mut endpoints = [0.0; 2];
    for (index, point) in [query.lower, query.upper].iter().enumerate() {
        let t = point
            .axes()
            .get(query.temperature_axis)
            .copied()
            .ok_or_else(|| invalid("source query lacks the temperature coordinate"))?;
        if !t.is_finite()
            || t <= 0.0
            || point
                .axis_quantities()
                .get(query.temperature_axis)
                .is_some_and(|q| {
                    *q != absolute && *q != QuantitySpec::dimensional(Dims([0, 0, 0, 1, 0, 0]))
                })
        {
            return Err(invalid(
                "source temperature must be positive absolute kelvin",
            ));
        }
        endpoints[index] = t;
    }
    if endpoints[0] >= endpoints[1]
        || query.lower.axis_quantities() != query.upper.axis_quantities()
        || !query
            .lower
            .axes()
            .iter()
            .filter(|(axis, _)| axis.as_str() != query.temperature_axis)
            .eq(query
                .upper
                .axes()
                .iter()
                .filter(|(axis, _)| axis.as_str() != query.temperature_axis))
    {
        return Err(invalid(
            "source temperatures must increase with fixed non-temperature context and descriptors",
        ));
    }
    Ok(endpoints)
}

fn heat_capacity_source_temperatures(
    claims: &ClaimSet,
    query: HeatCapacitySourceQuery<'_>,
    endpoints: [f64; 2],
) -> Result<Vec<f64>, PhaseStateError> {
    let mut temperatures = endpoints.to_vec();
    for key in [query.heat_capacity, query.density] {
        let slice = claims
            .query_line_typed(
                key,
                query.temperature_axis,
                query.lower,
                query.upper,
                ClaimSelection::Policy(query.selection),
            )
            .map_err(|error| PhaseStateError::SourceQuery {
                property: key.name().to_owned(),
                reason: format!("continuous support required: {error:?}"),
            })?;
        let selected = slice.knots[0].1.receipt.selected;
        let claim = claims
            .claim(selected)
            .expect("canonical source line selected an immutable claim");
        if claim.value.axes().any(|(axis, dims)| {
            axis == query.temperature_axis && dims != fs_qty::Temperature::DIMS
        }) {
            return Err(PhaseStateError::SourceQuery {
                property: key.name().to_owned(),
                reason: "source temperature coordinate must have temperature dimensions".into(),
            });
        }
        for (t, _) in slice.knots {
            if t > endpoints[0] && t < endpoints[1] {
                if temperatures.len() == MAX_HEAT_CAPACITY_ENTHALPY_KNOTS {
                    temperatures.sort_by(f64::total_cmp);
                    temperatures.dedup();
                    if temperatures
                        .binary_search_by(|value| value.total_cmp(&t))
                        .is_ok()
                    {
                        continue;
                    }
                    if temperatures.len() == MAX_HEAT_CAPACITY_ENTHALPY_KNOTS {
                        return Err(PhaseStateError::InvalidCurve {
                            what: "combined source temperature grid exceeds the knot budget",
                        });
                    }
                }
                temperatures.push(t);
            }
        }
    }
    temperatures.sort_by(f64::total_cmp);
    temperatures.dedup();
    Ok(temperatures)
}

/// Coarse phase topology selected from an equilibrium mass fraction.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SolidLiquidPhase {
    /// Exactly zero liquid mass fraction.
    Solid,
    /// Coexisting solid and liquid mass fractions.
    SolidLiquid,
    /// Exactly unit liquid mass fraction.
    Liquid,
}

/// State resolved from one admitted specific enthalpy.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EquilibriumPhaseState {
    /// Specific enthalpy [J/kg].
    specific_enthalpy_j_kg: f64,
    /// Absolute temperature [K].
    temperature_k: f64,
    /// Solid mass fraction in `[0, 1]`.
    solid_mass_fraction: f64,
    /// Liquid mass fraction in `[0, 1]`.
    liquid_mass_fraction: f64,
    /// Equilibrium bulk density [kg/m3].
    bulk_density_kg_m3: f64,
    /// Exact coarse phase topology.
    phase: SolidLiquidPhase,
    /// Material card that owns the equilibrium curve.
    material_card_identity: ContentHash,
    /// Identity of the complete admitted equilibrium curve.
    phase_curve_identity: ContentHash,
    /// Identity binding the curve and this resolved state.
    identity: ContentHash,
}

impl EquilibriumPhaseState {
    /// Specific enthalpy [J/kg].
    #[must_use]
    pub const fn specific_enthalpy_j_kg(self) -> f64 {
        self.specific_enthalpy_j_kg
    }

    /// Absolute temperature [K].
    #[must_use]
    pub const fn temperature_k(self) -> f64 {
        self.temperature_k
    }

    /// Solid mass fraction in `[0, 1]`.
    #[must_use]
    pub const fn solid_mass_fraction(self) -> f64 {
        self.solid_mass_fraction
    }

    /// Liquid mass fraction in `[0, 1]`.
    #[must_use]
    pub const fn liquid_mass_fraction(self) -> f64 {
        self.liquid_mass_fraction
    }

    /// Equilibrium bulk density [kg/m3].
    #[must_use]
    pub const fn bulk_density_kg_m3(self) -> f64 {
        self.bulk_density_kg_m3
    }

    /// Exact coarse phase topology.
    #[must_use]
    pub const fn phase(self) -> SolidLiquidPhase {
        self.phase
    }

    /// Material card that owns the equilibrium curve.
    #[must_use]
    pub const fn material_card_identity(self) -> ContentHash {
        self.material_card_identity
    }

    /// Identity of the complete admitted equilibrium curve.
    #[must_use]
    pub const fn phase_curve_identity(self) -> ContentHash {
        self.phase_curve_identity
    }

    /// Identity binding the curve and this resolved state.
    #[must_use]
    pub const fn identity(self) -> ContentHash {
        self.identity
    }

    /// Whether a solid-only constitutive/dynamics rung may consider this state.
    ///
    /// This is necessary but not sufficient admission: a solid consumer must
    /// still resolve its own temperature-, pressure-, rate-, and history-
    /// dependent properties at the same state.
    #[must_use]
    pub const fn is_fully_solid(self) -> bool {
        matches!(self.phase, SolidLiquidPhase::Solid)
    }
}

/// Immutable equilibrium phase curve bound to one material-card identity.
#[derive(Clone, Debug, PartialEq)]
pub struct EquilibriumEnthalpyPhaseCurve {
    material_card_identity: ContentHash,
    knots: Vec<EnthalpyPhaseKnot>,
    identity: ContentHash,
}

impl EquilibriumEnthalpyPhaseCurve {
    /// Admit a bounded, monotone solid-liquid enthalpy curve.
    ///
    /// At least two knots are required. Enthalpy must increase strictly;
    /// temperature and liquid fraction must not decrease. The first knot must
    /// be fully solid and the last fully liquid so that the curve owns both
    /// phase transitions rather than exposing a misleading partial taxonomy.
    pub fn try_new(
        material_card_identity: ContentHash,
        knots: Vec<EnthalpyPhaseKnot>,
    ) -> Result<Self, PhaseStateError> {
        validate_curve_common(material_card_identity, &knots)?;
        if knots[0].liquid_mass_fraction != 0.0
            || knots[knots.len() - 1]
                .liquid_mass_fraction
                .partial_cmp(&1.0)
                != Some(core::cmp::Ordering::Equal)
        {
            return Err(PhaseStateError::InvalidCurve {
                what: "the curve must begin fully solid and end fully liquid",
            });
        }
        for pair in knots.windows(2) {
            if pair[1].temperature_k < pair[0].temperature_k {
                return Err(PhaseStateError::InvalidCurve {
                    what: "temperature must be nondecreasing with specific enthalpy",
                });
            }
            if pair[1].liquid_mass_fraction < pair[0].liquid_mass_fraction {
                return Err(PhaseStateError::InvalidCurve {
                    what: "liquid mass fraction must be nondecreasing with specific enthalpy",
                });
            }
        }
        Ok(Self::from_valid_knots(material_card_identity, knots))
    }

    /// Admit a bounded equilibrium enthalpy curve entirely in one explicit
    /// phase. Solid curves carry liquid fraction zero; liquid curves carry
    /// liquid fraction one. Temperature must rise strictly, so this ingress
    /// cannot hide a latent plateau or mixed transition.
    pub fn try_single_phase(
        material_card_identity: ContentHash,
        phase: SolidLiquidPhase,
        knots: Vec<EnthalpyPhaseKnot>,
    ) -> Result<Self, PhaseStateError> {
        validate_curve_common(material_card_identity, &knots)?;
        let liquid_mass_fraction = single_phase_fraction(phase)?;
        if knots.iter().any(|knot| {
            knot.liquid_mass_fraction.partial_cmp(&liquid_mass_fraction)
                != Some(core::cmp::Ordering::Equal)
        }) {
            return Err(PhaseStateError::InvalidCurve {
                what: "every single-phase knot must carry the declared phase fraction",
            });
        }
        if knots
            .windows(2)
            .any(|pair| pair[1].temperature_k <= pair[0].temperature_k)
        {
            return Err(PhaseStateError::InvalidCurve {
                what: "single-phase temperature knots must be strictly increasing",
            });
        }
        Ok(Self::from_valid_knots(material_card_identity, knots))
    }

    /// Construct a bounded one-phase enthalpy chart from source heat-capacity
    /// knots. Heat capacity and density are linearly interpolated in
    /// temperature. Each segment's heat capacity is integrated exactly; it is
    /// subdivided until its straight enthalpy chord differs from that integral
    /// by no more than `maximum_interpolation_error_j_kg`.
    ///
    /// The requested bound covers only this interpolation arithmetic. It does
    /// not quantify source error, floating-point roundoff, or the density
    /// interpolation approximation.
    pub fn try_from_heat_capacity(
        material_card_identity: ContentHash,
        phase: SolidLiquidPhase,
        reference_specific_enthalpy_j_kg: f64,
        knots: &[HeatCapacityKnot],
        maximum_interpolation_error_j_kg: f64,
    ) -> Result<Self, PhaseStateError> {
        let liquid_mass_fraction = single_phase_fraction(phase)?;
        let generated_count = heat_capacity_generated_knot_count(
            material_card_identity,
            reference_specific_enthalpy_j_kg,
            knots,
            maximum_interpolation_error_j_kg,
        )?;
        let mut generated = Vec::with_capacity(generated_count);
        generated.push(EnthalpyPhaseKnot {
            specific_enthalpy_j_kg: reference_specific_enthalpy_j_kg,
            temperature_k: knots[0].temperature_k,
            liquid_mass_fraction,
            bulk_density_kg_m3: knots[0].bulk_density_kg_m3,
        });
        let mut enthalpy = reference_specific_enthalpy_j_kg;
        for pair in knots.windows(2) {
            let segment_start_enthalpy = enthalpy;
            let delta_temperature = pair[1].temperature_k - pair[0].temperature_k;
            let enthalpy_gain = (0.5 * pair[0].specific_heat_capacity_j_kg_k
                + 0.5 * pair[1].specific_heat_capacity_j_kg_k)
                * delta_temperature;
            let endpoint_enthalpy = enthalpy + enthalpy_gain;
            let curvature_error = (pair[1].specific_heat_capacity_j_kg_k
                - pair[0].specific_heat_capacity_j_kg_k)
                .abs()
                * delta_temperature
                / 8.0;
            let subdivisions = if curvature_error <= maximum_interpolation_error_j_kg {
                1
            } else {
                (curvature_error / maximum_interpolation_error_j_kg)
                    .sqrt()
                    .ceil() as usize
            };
            for step in 1..=subdivisions {
                let fraction = step as f64 / subdivisions as f64;
                let delta = delta_temperature * fraction;
                let temperature_k = if step == subdivisions {
                    pair[1].temperature_k
                } else {
                    pair[0].temperature_k + delta
                };
                let specific_enthalpy_j_kg = if step == subdivisions {
                    endpoint_enthalpy
                } else {
                    segment_start_enthalpy
                        + pair[0].specific_heat_capacity_j_kg_k * delta
                        + 0.5
                            * (pair[1].specific_heat_capacity_j_kg_k
                                - pair[0].specific_heat_capacity_j_kg_k)
                            / delta_temperature
                            * delta
                            * delta
                };
                let bulk_density_kg_m3 = if step == subdivisions {
                    pair[1].bulk_density_kg_m3
                } else {
                    pair[0].bulk_density_kg_m3
                        + (pair[1].bulk_density_kg_m3 - pair[0].bulk_density_kg_m3) * fraction
                };
                if !(temperature_k.is_finite()
                    && specific_enthalpy_j_kg.is_finite()
                    && specific_enthalpy_j_kg > enthalpy
                    && bulk_density_kg_m3.is_finite()
                    && bulk_density_kg_m3 > 0.0)
                {
                    return Err(PhaseStateError::InvalidCurve {
                        what: "generated heat-capacity enthalpy knot is invalid",
                    });
                }
                generated.push(EnthalpyPhaseKnot {
                    specific_enthalpy_j_kg,
                    temperature_k,
                    liquid_mass_fraction,
                    bulk_density_kg_m3,
                });
                enthalpy = specific_enthalpy_j_kg;
            }
        }
        Self::try_single_phase(material_card_identity, phase, generated)
    }

    fn from_valid_knots(
        material_card_identity: ContentHash,
        knots: Vec<EnthalpyPhaseKnot>,
    ) -> Self {
        Self {
            material_card_identity,
            identity: phase_curve_identity(material_card_identity, &knots),
            knots,
        }
    }

    /// Material card that owns the supplied equilibrium data.
    #[must_use]
    pub const fn material_card_identity(&self) -> ContentHash {
        self.material_card_identity
    }

    /// Complete phase-curve identity.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    /// Admitted curve knots in ascending specific enthalpy.
    #[must_use]
    pub fn knots(&self) -> &[EnthalpyPhaseKnot] {
        &self.knots
    }

    /// Piecewise-linear temperature derivative `dT/dh` [K/(J/kg)] on this
    /// specific-enthalpy chart, without extrapolation or phase smoothing.
    ///
    /// Inside a segment this is its exact chord slope, including exactly zero
    /// throughout an isothermal latent-heat plateau. At an interior knot the
    /// segment to the right (increasing enthalpy) is selected. The lower domain
    /// endpoint also uses its right derivative; the upper endpoint uses its
    /// left derivative. This explicit one-sided convention provides a
    /// generalized derivative for semismooth enthalpy solves; it does not claim
    /// classical differentiability at a slope-changing knot.
    ///
    /// For a chart generated from heat-capacity data, this differentiates the
    /// admitted piecewise-linear chart, not the original continuous Cp curve.
    /// No reciprocal or artificially large apparent heat capacity is formed
    /// on a latent plateau.
    ///
    /// # Errors
    /// Returns the same finite-input and domain refusals as
    /// [`Self::state_at_specific_enthalpy`]. A nonzero segment slope that
    /// overflows or underflows to zero returns
    /// [`PhaseStateError::UnrepresentableTemperatureDerivative`].
    pub fn temperature_derivative_at_specific_enthalpy(
        &self,
        specific_enthalpy_j_kg: f64,
    ) -> Result<f64, PhaseStateError> {
        let (lower_knot, upper_knot) = self.derivative_segment(specific_enthalpy_j_kg)?;
        let temperature_increment = upper_knot.temperature_k - lower_knot.temperature_k;
        let slope = temperature_increment
            / (upper_knot.specific_enthalpy_j_kg - lower_knot.specific_enthalpy_j_kg);
        if !slope.is_finite() || (temperature_increment > 0.0 && slope == 0.0) {
            return Err(PhaseStateError::UnrepresentableTemperatureDerivative);
        }
        Ok(slope)
    }

    /// Exact piecewise-linear derivative of liquid **mass** fraction with
    /// respect to specific enthalpy, in kg/J.
    ///
    /// Uses the same right derivative at interior knots and left derivative
    /// at the upper endpoint as [`Self::temperature_derivative_at_specific_enthalpy`].
    /// A latent plateau can have a positive fraction derivative while its
    /// temperature derivative is exactly zero. No phase smoothing, volume
    /// fraction conversion or extrapolation is performed.
    ///
    /// # Errors
    /// Refuses nonfinite or out-of-domain enthalpy. A nonzero segment slope
    /// that cannot be represented returns
    /// [`PhaseStateError::UnrepresentableLiquidMassFractionDerivative`].
    pub fn liquid_mass_fraction_derivative_at_specific_enthalpy(
        &self,
        specific_enthalpy_j_kg: f64,
    ) -> Result<f64, PhaseStateError> {
        let (lower_knot, upper_knot) = self.derivative_segment(specific_enthalpy_j_kg)?;
        let fraction_increment = upper_knot.liquid_mass_fraction - lower_knot.liquid_mass_fraction;
        let slope = fraction_increment
            / (upper_knot.specific_enthalpy_j_kg - lower_knot.specific_enthalpy_j_kg);
        if !slope.is_finite() || (fraction_increment > 0.0 && slope == 0.0) {
            return Err(PhaseStateError::UnrepresentableLiquidMassFractionDerivative);
        }
        Ok(slope)
    }

    fn derivative_segment(
        &self,
        specific_enthalpy_j_kg: f64,
    ) -> Result<(EnthalpyPhaseKnot, EnthalpyPhaseKnot), PhaseStateError> {
        if !specific_enthalpy_j_kg.is_finite() {
            return Err(PhaseStateError::NonFiniteSpecificEnthalpy);
        }
        let lower = self.knots[0].specific_enthalpy_j_kg;
        let upper = self.knots[self.knots.len() - 1].specific_enthalpy_j_kg;
        if specific_enthalpy_j_kg < lower || specific_enthalpy_j_kg > upper {
            return Err(PhaseStateError::OutsideEnthalpyDomain {
                specific_enthalpy_j_kg,
                lower_j_kg: lower,
                upper_j_kg: upper,
            });
        }
        let upper_index = self
            .knots
            .partition_point(|knot| knot.specific_enthalpy_j_kg <= specific_enthalpy_j_kg)
            .min(self.knots.len() - 1);
        Ok((self.knots[upper_index - 1], self.knots[upper_index]))
    }

    /// Resolve temperature in kelvin without constructing a phase-state identity.
    ///
    /// Uses the same chart evaluation and admission checks as
    /// [`Self::state_at_specific_enthalpy`], including latent plateaus and
    /// representability of interpolated density.
    ///
    /// # Errors
    /// Refuses nonfinite or out-of-domain enthalpy and an interpolated density
    /// whose specific-volume calculation is not finite and positive.
    pub fn temperature_at_specific_enthalpy(
        &self,
        specific_enthalpy_j_kg: f64,
    ) -> Result<f64, PhaseStateError> {
        self.values_at_specific_enthalpy(specific_enthalpy_j_kg)
            .map(|values| values.temperature_k)
    }

    /// Resolve equilibrium temperature, density, and phase fractions from
    /// specific enthalpy without extrapolation.
    ///
    /// # Errors
    /// Refuses nonfinite or out-of-domain enthalpy and an interpolated density
    /// whose specific-volume calculation is not finite and positive.
    pub fn state_at_specific_enthalpy(
        &self,
        specific_enthalpy_j_kg: f64,
    ) -> Result<EquilibriumPhaseState, PhaseStateError> {
        let values = self.values_at_specific_enthalpy(specific_enthalpy_j_kg)?;
        Ok(self.state_from_values(specific_enthalpy_j_kg, values))
    }

    fn values_at_specific_enthalpy(
        &self,
        specific_enthalpy_j_kg: f64,
    ) -> Result<EnthalpyPhaseKnot, PhaseStateError> {
        if !specific_enthalpy_j_kg.is_finite() {
            return Err(PhaseStateError::NonFiniteSpecificEnthalpy);
        }
        let lower = self.knots[0].specific_enthalpy_j_kg;
        let upper = self.knots[self.knots.len() - 1].specific_enthalpy_j_kg;
        if specific_enthalpy_j_kg < lower || specific_enthalpy_j_kg > upper {
            return Err(PhaseStateError::OutsideEnthalpyDomain {
                specific_enthalpy_j_kg,
                lower_j_kg: lower,
                upper_j_kg: upper,
            });
        }
        match self.knots.binary_search_by(|knot| {
            knot.specific_enthalpy_j_kg
                .partial_cmp(&specific_enthalpy_j_kg)
                .expect("admitted enthalpy knots and query are finite")
        }) {
            Ok(index) => Ok(self.knots[index]),
            Err(upper_index) => {
                let lower_knot = self.knots[upper_index - 1];
                let upper_knot = self.knots[upper_index];
                let alpha = (specific_enthalpy_j_kg - lower_knot.specific_enthalpy_j_kg)
                    / (upper_knot.specific_enthalpy_j_kg - lower_knot.specific_enthalpy_j_kg);
                let values = EnthalpyPhaseKnot {
                    specific_enthalpy_j_kg,
                    temperature_k: lerp(lower_knot.temperature_k, upper_knot.temperature_k, alpha),
                    liquid_mass_fraction: lerp(
                        lower_knot.liquid_mass_fraction,
                        upper_knot.liquid_mass_fraction,
                        alpha,
                    ),
                    // Density is mass per volume, so interpolation on a
                    // mass-specific enthalpy coordinate is performed in
                    // specific volume and inverted. Linear density would
                    // violate additive mixture volume through a two-phase
                    // interval.
                    // An admissible constant-density segment must retain
                    // its mass density bit for bit: reciprocal round trips may
                    // move even a constant value by 1 ULP. Specific volume
                    // must still be representable, as for every other segment.
                    bulk_density_kg_m3: if lower_knot.bulk_density_kg_m3.to_bits()
                        == upper_knot.bulk_density_kg_m3.to_bits()
                        && lower_knot.bulk_density_kg_m3.recip().is_finite()
                    {
                        lower_knot.bulk_density_kg_m3
                    } else {
                        lerp(
                            lower_knot.bulk_density_kg_m3.recip(),
                            upper_knot.bulk_density_kg_m3.recip(),
                            alpha,
                        )
                        .recip()
                    },
                };
                if !(values.bulk_density_kg_m3.is_finite() && values.bulk_density_kg_m3 > 0.0) {
                    return Err(PhaseStateError::InvalidCurve {
                        what: "interpolated bulk density is not finite and positive",
                    });
                }
                Ok(values)
            }
        }
    }

    /// Apply a signed specific-energy increment and resolve the resulting
    /// phase state. No heat source, boundary flux, or work term is invented;
    /// the caller owns that energy balance.
    pub fn advance_specific_energy(
        &self,
        current: EquilibriumPhaseState,
        net_specific_energy_j_kg: f64,
    ) -> Result<EquilibriumPhaseState, PhaseStateError> {
        if current.phase_curve_identity() != self.identity {
            return Err(PhaseStateError::StateCurveMismatch);
        }
        if !net_specific_energy_j_kg.is_finite() {
            return Err(PhaseStateError::NonFiniteSpecificEnergyIncrement);
        }
        let next = current.specific_enthalpy_j_kg() + net_specific_energy_j_kg;
        if !next.is_finite() {
            return Err(PhaseStateError::NonFiniteSpecificEnthalpy);
        }
        self.state_at_specific_enthalpy(next)
    }

    fn state_from_values(
        &self,
        specific_enthalpy_j_kg: f64,
        values: EnthalpyPhaseKnot,
    ) -> EquilibriumPhaseState {
        let liquid_mass_fraction = values.liquid_mass_fraction.clamp(0.0, 1.0);
        let solid_mass_fraction = 1.0 - liquid_mass_fraction;
        let phase = if liquid_mass_fraction == 0.0 {
            SolidLiquidPhase::Solid
        } else if liquid_mass_fraction.partial_cmp(&1.0) == Some(core::cmp::Ordering::Equal) {
            SolidLiquidPhase::Liquid
        } else {
            SolidLiquidPhase::SolidLiquid
        };
        let mut hasher = DomainHasher::new(EQUILIBRIUM_PHASE_STATE_IDENTITY_DOMAIN);
        hasher.update(self.identity.as_bytes());
        for value in [
            specific_enthalpy_j_kg,
            values.temperature_k,
            solid_mass_fraction,
            liquid_mass_fraction,
            values.bulk_density_kg_m3,
        ] {
            hasher.update(&value.to_bits().to_le_bytes());
        }
        hasher.update(&[match phase {
            SolidLiquidPhase::Solid => 0,
            SolidLiquidPhase::SolidLiquid => 1,
            SolidLiquidPhase::Liquid => 2,
        }]);
        EquilibriumPhaseState {
            specific_enthalpy_j_kg,
            temperature_k: values.temperature_k,
            solid_mass_fraction,
            liquid_mass_fraction,
            bulk_density_kg_m3: values.bulk_density_kg_m3,
            phase,
            material_card_identity: self.material_card_identity,
            phase_curve_identity: self.identity,
            identity: hasher.finalize(),
        }
    }
}

fn heat_capacity_generated_knot_count(
    material_card_identity: ContentHash,
    reference_specific_enthalpy_j_kg: f64,
    knots: &[HeatCapacityKnot],
    maximum_interpolation_error_j_kg: f64,
) -> Result<usize, PhaseStateError> {
    if material_card_identity == ContentHash([0; 32]) {
        return Err(PhaseStateError::InvalidCurve {
            what: "material-card identity must not be zero",
        });
    }
    if !reference_specific_enthalpy_j_kg.is_finite() {
        return Err(PhaseStateError::InvalidCurve {
            what: "reference specific enthalpy must be finite",
        });
    }
    if !(maximum_interpolation_error_j_kg.is_finite() && maximum_interpolation_error_j_kg > 0.0) {
        return Err(PhaseStateError::InvalidCurve {
            what: "maximum interpolation error must be finite and positive",
        });
    }
    if knots.len() < 2 {
        return Err(PhaseStateError::InvalidCurve {
            what: "a heat-capacity chart needs at least two knots",
        });
    }
    if knots.len() > MAX_HEAT_CAPACITY_ENTHALPY_KNOTS {
        return Err(PhaseStateError::InvalidCurve {
            what: "heat-capacity source knot count exceeds the generated-knot budget",
        });
    }
    let mut generated_count = 1_usize;
    let mut enthalpy = reference_specific_enthalpy_j_kg;
    for knot in knots {
        if !(knot.temperature_k.is_finite()
            && knot.temperature_k > 0.0
            && knot.specific_heat_capacity_j_kg_k.is_finite()
            && knot.specific_heat_capacity_j_kg_k > 0.0
            && knot.bulk_density_kg_m3.is_finite()
            && knot.bulk_density_kg_m3 > 0.0)
        {
            return Err(PhaseStateError::InvalidCurve {
                what: "heat-capacity knots need finite positive temperature, heat capacity, and density",
            });
        }
    }
    for pair in knots.windows(2) {
        let delta_temperature = pair[1].temperature_k - pair[0].temperature_k;
        if !(delta_temperature.is_finite() && delta_temperature > 0.0) {
            return Err(PhaseStateError::InvalidCurve {
                what: "heat-capacity temperatures must be strictly increasing",
            });
        }
        let curvature_error =
            (pair[1].specific_heat_capacity_j_kg_k - pair[0].specific_heat_capacity_j_kg_k).abs()
                * delta_temperature
                / 8.0;
        if !curvature_error.is_finite() {
            return Err(PhaseStateError::InvalidCurve {
                what: "heat-capacity interpolation error is not finite",
            });
        }
        let subdivisions = if curvature_error <= maximum_interpolation_error_j_kg {
            1
        } else {
            let required = (curvature_error / maximum_interpolation_error_j_kg)
                .sqrt()
                .ceil();
            if !(required.is_finite() && required <= MAX_HEAT_CAPACITY_ENTHALPY_KNOTS as f64) {
                return Err(PhaseStateError::InvalidCurve {
                    what: "heat-capacity interpolation exceeds the generated-knot budget",
                });
            }
            required as usize
        };
        generated_count =
            generated_count
                .checked_add(subdivisions)
                .ok_or(PhaseStateError::InvalidCurve {
                    what: "heat-capacity interpolation exceeds the generated-knot budget",
                })?;
        if generated_count > MAX_HEAT_CAPACITY_ENTHALPY_KNOTS {
            return Err(PhaseStateError::InvalidCurve {
                what: "heat-capacity interpolation exceeds the generated-knot budget",
            });
        }
        let enthalpy_gain = (0.5 * pair[0].specific_heat_capacity_j_kg_k
            + 0.5 * pair[1].specific_heat_capacity_j_kg_k)
            * delta_temperature;
        let next_enthalpy = enthalpy + enthalpy_gain;
        if !(enthalpy_gain.is_finite() && next_enthalpy.is_finite() && next_enthalpy > enthalpy) {
            return Err(PhaseStateError::InvalidCurve {
                what: "heat-capacity integration produced non-finite or collapsed enthalpy",
            });
        }
        enthalpy = next_enthalpy;
    }
    Ok(generated_count)
}

fn validate_curve_common(
    material_card_identity: ContentHash,
    knots: &[EnthalpyPhaseKnot],
) -> Result<(), PhaseStateError> {
    if material_card_identity == ContentHash([0; 32]) {
        return Err(PhaseStateError::InvalidCurve {
            what: "material-card identity must not be zero",
        });
    }
    if knots.len() < 2 {
        return Err(PhaseStateError::InvalidCurve {
            what: "an equilibrium phase curve needs at least two knots",
        });
    }
    for knot in knots {
        if !(knot.specific_enthalpy_j_kg.is_finite()
            && knot.temperature_k > 0.0
            && knot.temperature_k.is_finite()
            && (0.0..=1.0).contains(&knot.liquid_mass_fraction)
            && knot.bulk_density_kg_m3 > 0.0
            && knot.bulk_density_kg_m3.is_finite())
        {
            return Err(PhaseStateError::InvalidCurve {
                what: "every knot needs finite enthalpy, positive temperature/density, and liquid fraction in [0,1]",
            });
        }
    }
    for pair in knots.windows(2) {
        if pair[1].specific_enthalpy_j_kg <= pair[0].specific_enthalpy_j_kg
            || !(pair[1].specific_enthalpy_j_kg - pair[0].specific_enthalpy_j_kg).is_finite()
        {
            return Err(PhaseStateError::InvalidCurve {
                what: "specific enthalpy knots must be strictly increasing",
            });
        }
    }
    Ok(())
}

fn single_phase_fraction(phase: SolidLiquidPhase) -> Result<f64, PhaseStateError> {
    match phase {
        SolidLiquidPhase::Solid => Ok(0.0),
        SolidLiquidPhase::Liquid => Ok(1.0),
        SolidLiquidPhase::SolidLiquid => Err(PhaseStateError::InvalidCurve {
            what: "a single-phase curve must declare Solid or Liquid",
        }),
    }
}

/// Typed refusal from equilibrium phase-state admission or evaluation.
#[derive(Clone, Debug, PartialEq)]
pub enum PhaseStateError {
    /// A source query or its continuous support could not be admitted.
    SourceQuery {
        /// Requested source property or coordinate.
        property: String,
        /// Original query diagnosis.
        reason: String,
    },
    /// The source-provided curve violates a physical or ordering invariant.
    InvalidCurve {
        /// Failed curve invariant.
        what: &'static str,
    },
    /// A queried specific enthalpy was not finite.
    NonFiniteSpecificEnthalpy,
    /// A segment's positive `dT/dh` overflowed or underflowed to zero.
    /// A genuine isothermal plateau instead returns exactly zero successfully.
    UnrepresentableTemperatureDerivative,
    /// A positive liquid-mass-fraction slope overflowed or underflowed to zero.
    UnrepresentableLiquidMassFractionDerivative,
    /// A signed energy increment was not finite.
    NonFiniteSpecificEnergyIncrement,
    /// Evaluation would extrapolate beyond the source-provided curve.
    OutsideEnthalpyDomain {
        /// Refused query [J/kg].
        specific_enthalpy_j_kg: f64,
        /// Admitted lower endpoint [J/kg].
        lower_j_kg: f64,
        /// Admitted upper endpoint [J/kg].
        upper_j_kg: f64,
    },
    /// A proposed constant-pressure heat step violates its physical or numerical domain.
    InvalidIsobaricHeatStep {
        /// Failed step invariant.
        what: &'static str,
    },
    /// A state from a different phase curve was supplied to an update.
    StateCurveMismatch,
}

impl fmt::Display for PhaseStateError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{self:?}")
    }
}

impl std::error::Error for PhaseStateError {}

fn phase_curve_identity(
    material_card_identity: ContentHash,
    knots: &[EnthalpyPhaseKnot],
) -> ContentHash {
    let mut hasher = DomainHasher::new(EQUILIBRIUM_ENTHALPY_PHASE_IDENTITY_DOMAIN);
    hasher.update(material_card_identity.as_bytes());
    hasher.update(&u64::try_from(knots.len()).unwrap_or(u64::MAX).to_le_bytes());
    for knot in knots {
        for value in [
            knot.specific_enthalpy_j_kg,
            knot.temperature_k,
            knot.liquid_mass_fraction,
            knot.bulk_density_kg_m3,
        ] {
            hasher.update(&value.to_bits().to_le_bytes());
        }
    }
    hasher.finalize()
}

fn lerp(start: f64, end: f64, alpha: f64) -> f64 {
    (end - start).mul_add(alpha, start)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn plateau_curve() -> EquilibriumEnthalpyPhaseCurve {
        EquilibriumEnthalpyPhaseCurve::try_new(
            ContentHash([0x51; 32]),
            vec![
                EnthalpyPhaseKnot {
                    specific_enthalpy_j_kg: 0.0,
                    temperature_k: 300.0,
                    liquid_mass_fraction: 0.0,
                    bulk_density_kg_m3: 11_300.0,
                },
                EnthalpyPhaseKnot {
                    specific_enthalpy_j_kg: 30_000.0,
                    temperature_k: 600.0,
                    liquid_mass_fraction: 0.0,
                    bulk_density_kg_m3: 11_100.0,
                },
                EnthalpyPhaseKnot {
                    specific_enthalpy_j_kg: 55_000.0,
                    temperature_k: 600.0,
                    liquid_mass_fraction: 1.0,
                    bulk_density_kg_m3: 10_600.0,
                },
                EnthalpyPhaseKnot {
                    specific_enthalpy_j_kg: 95_000.0,
                    temperature_k: 800.0,
                    liquid_mass_fraction: 1.0,
                    bulk_density_kg_m3: 10_300.0,
                },
            ],
        )
        .unwrap()
    }

    #[test]
    fn temperature_queries_match_phase_states_across_sensible_and_latent_charts() {
        use fs_thermal::enthalpy::{DepositionError, ThermalLaw};

        let [cold, warm] = heat_capacity_knots(1_000.0, 2_000.0);
        let hot = HeatCapacityKnot {
            temperature_k: 600.0,
            specific_heat_capacity_j_kg_k: 500.0,
            bulk_density_kg_m3: 800.0,
        };
        let sensible = EquilibriumEnthalpyPhaseCurve::try_from_heat_capacity(
            ContentHash([0x75; 32]),
            SolidLiquidPhase::Solid,
            0.0,
            &[cold, warm, hot],
            1.0,
        )
        .unwrap();
        for curve in [plateau_curve(), sensible] {
            let check = |h| {
                let expected = curve.state_at_specific_enthalpy(h).unwrap().temperature_k();
                assert_eq!(
                    curve.temperature_at_specific_enthalpy(h).unwrap().to_bits(),
                    expected.to_bits()
                );
                assert_eq!(
                    ThermalLaw::temperature_k(&curve, h).unwrap().to_bits(),
                    expected.to_bits()
                );
            };
            for h in [-0.0, 0.0] {
                check(h);
            }
            for pair in curve.knots().windows(2) {
                for alpha in [0.0, 0.125, 0.5, 0.875, 1.0] {
                    check(lerp(
                        pair[0].specific_enthalpy_j_kg,
                        pair[1].specific_enthalpy_j_kg,
                        alpha,
                    ));
                }
            }
            let upper = curve.knots().last().unwrap().specific_enthalpy_j_kg;
            for h in [f64::NAN, f64::NEG_INFINITY, f64::INFINITY, -1.0, upper + 1.0] {
                assert_eq!(
                    curve.temperature_at_specific_enthalpy(h),
                    curve.state_at_specific_enthalpy(h).map(|state| state.temperature_k())
                );
                assert_eq!(
                    ThermalLaw::temperature_k(&curve, h),
                    Err(DepositionError::MaterialDomain)
                );
            }
        }
    }

    #[test]
    fn g0_enthalpy_crosses_an_isothermal_phase_plateau_without_name_switches() {
        let curve = plateau_curve();
        let solid = curve.state_at_specific_enthalpy(20_000.0).unwrap();
        assert_eq!(solid.phase(), SolidLiquidPhase::Solid);
        assert_eq!(solid.temperature_k().to_bits(), 500.0_f64.to_bits());

        let half_melted = curve.advance_specific_energy(solid, 22_500.0).unwrap();
        assert_eq!(half_melted.phase(), SolidLiquidPhase::SolidLiquid);
        assert_eq!(half_melted.temperature_k().to_bits(), 600.0_f64.to_bits());
        assert_eq!(
            half_melted.solid_mass_fraction().to_bits(),
            0.5_f64.to_bits()
        );
        assert_eq!(
            half_melted.liquid_mass_fraction().to_bits(),
            0.5_f64.to_bits()
        );
        assert!(!half_melted.is_fully_solid());

        let liquid = curve
            .advance_specific_energy(half_melted, 32_500.0)
            .unwrap();
        assert_eq!(liquid.phase(), SolidLiquidPhase::Liquid);
        assert_eq!(liquid.temperature_k().to_bits(), 700.0_f64.to_bits());
        assert_eq!(liquid.liquid_mass_fraction().to_bits(), 1.0_f64.to_bits());
    }

    #[test]
    fn g0_phase_curve_refuses_extrapolation_and_foreign_state() {
        let curve = plateau_curve();
        assert!(matches!(
            curve.state_at_specific_enthalpy(95_000.0 + f64::EPSILON * 95_000.0),
            Err(PhaseStateError::OutsideEnthalpyDomain { .. })
        ));
        let other =
            EquilibriumEnthalpyPhaseCurve::try_new(ContentHash([0x52; 32]), curve.knots().to_vec())
                .unwrap();
        let state = other.state_at_specific_enthalpy(40_000.0).unwrap();
        assert_eq!(
            curve.advance_specific_energy(state, 1.0),
            Err(PhaseStateError::StateCurveMismatch)
        );
    }

    #[test]
    fn g0_phase_curve_refuses_nonmonotone_or_partial_taxonomy() {
        let mut knots = plateau_curve().knots().to_vec();
        knots[3].liquid_mass_fraction = 0.25;
        assert!(matches!(
            EquilibriumEnthalpyPhaseCurve::try_new(ContentHash([0x53; 32]), knots),
            Err(PhaseStateError::InvalidCurve { .. })
        ));

        let mut knots = plateau_curve().knots().to_vec();
        knots[2].temperature_k = 599.0;
        assert!(matches!(
            EquilibriumEnthalpyPhaseCurve::try_new(ContentHash([0x54; 32]), knots),
            Err(PhaseStateError::InvalidCurve { .. })
        ));
    }

    fn single_phase_knots(liquid_mass_fraction: f64) -> Vec<EnthalpyPhaseKnot> {
        vec![
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 10.0,
                temperature_k: 300.0,
                liquid_mass_fraction,
                bulk_density_kg_m3: 1_000.0,
            },
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 110.0,
                temperature_k: 350.0,
                liquid_mass_fraction,
                bulk_density_kg_m3: 900.0,
            },
        ]
    }

    #[test]
    fn phase_state_refuses_nonrepresentable_density_between_finite_source_knots() {
        let tiny = f64::from_bits(1);
        for (lower_density, upper_density) in [(tiny, 1.0), (1.0, tiny), (tiny, tiny)] {
            let mut knots = single_phase_knots(0.0);
            knots[0].bulk_density_kg_m3 = lower_density;
            knots[1].bulk_density_kg_m3 = upper_density;
            let midpoint =
                0.5 * (knots[0].specific_enthalpy_j_kg + knots[1].specific_enthalpy_j_kg);
            let curve = EquilibriumEnthalpyPhaseCurve::try_single_phase(
                ContentHash([0x59; 32]),
                SolidLiquidPhase::Solid,
                knots,
            )
            .expect("finite positive source knots remain admissible");
            for knot in curve.knots() {
                let exact = curve
                    .state_at_specific_enthalpy(knot.specific_enthalpy_j_kg)
                    .expect("exact source-knot density is representable");
                assert_eq!(
                    exact.bulk_density_kg_m3().to_bits(),
                    knot.bulk_density_kg_m3.to_bits()
                );
                assert_eq!(
                    curve.temperature_at_specific_enthalpy(knot.specific_enthalpy_j_kg),
                    Ok(knot.temperature_k)
                );
            }
            assert_eq!(
                curve.state_at_specific_enthalpy(midpoint),
                Err(PhaseStateError::InvalidCurve {
                    what: "interpolated bulk density is not finite and positive",
                })
            );
            assert_eq!(
                curve.temperature_at_specific_enthalpy(midpoint),
                Err(PhaseStateError::InvalidCurve {
                    what: "interpolated bulk density is not finite and positive",
                })
            );
            assert_eq!(
                fs_thermal::enthalpy::ThermalLaw::temperature_k(&curve, midpoint),
                Err(fs_thermal::enthalpy::DepositionError::MaterialDomain)
            );
        }
    }

    fn heat_capacity_knots(first_cp: f64, second_cp: f64) -> [HeatCapacityKnot; 2] {
        [
            HeatCapacityKnot {
                temperature_k: 300.0,
                specific_heat_capacity_j_kg_k: first_cp,
                bulk_density_kg_m3: 1_000.0,
            },
            HeatCapacityKnot {
                temperature_k: 400.0,
                specific_heat_capacity_j_kg_k: second_cp,
                bulk_density_kg_m3: 900.0,
            },
        ]
    }

    #[test]
    fn g0_single_phase_enthalpy_curves_interpolate_without_phase_transition() {
        let solid = EquilibriumEnthalpyPhaseCurve::try_single_phase(
            ContentHash([0x61; 32]),
            SolidLiquidPhase::Solid,
            single_phase_knots(0.0),
        )
        .unwrap();
        let solid_midpoint = solid.state_at_specific_enthalpy(60.0).unwrap();
        assert_eq!(solid_midpoint.phase(), SolidLiquidPhase::Solid);
        assert_eq!(
            solid_midpoint.liquid_mass_fraction().to_bits(),
            0.0_f64.to_bits()
        );
        assert_eq!(
            solid_midpoint.solid_mass_fraction().to_bits(),
            1.0_f64.to_bits()
        );
        assert_eq!(
            solid_midpoint.temperature_k().to_bits(),
            325.0_f64.to_bits()
        );
        assert!((solid_midpoint.bulk_density_kg_m3() - 947.368_421_052_631_6).abs() < 1e-12);

        let liquid = EquilibriumEnthalpyPhaseCurve::try_single_phase(
            ContentHash([0x62; 32]),
            SolidLiquidPhase::Liquid,
            single_phase_knots(1.0),
        )
        .unwrap();
        let liquid_midpoint = liquid.state_at_specific_enthalpy(60.0).unwrap();
        assert_eq!(liquid_midpoint.phase(), SolidLiquidPhase::Liquid);
        assert_eq!(
            liquid_midpoint.liquid_mass_fraction().to_bits(),
            1.0_f64.to_bits()
        );
        assert_eq!(
            liquid_midpoint.solid_mass_fraction().to_bits(),
            0.0_f64.to_bits()
        );
        assert_eq!(
            liquid_midpoint.temperature_k().to_bits(),
            325.0_f64.to_bits()
        );
    }

    #[test]
    fn g0_single_phase_enthalpy_curves_refuse_transition_plateau_and_extrapolation() {
        let identity = ContentHash([0x63; 32]);
        let mut mixed_fraction = single_phase_knots(0.0);
        mixed_fraction[1].liquid_mass_fraction = 0.5;
        assert!(matches!(
            EquilibriumEnthalpyPhaseCurve::try_single_phase(
                identity,
                SolidLiquidPhase::Solid,
                mixed_fraction,
            ),
            Err(PhaseStateError::InvalidCurve { .. })
        ));
        assert!(matches!(
            EquilibriumEnthalpyPhaseCurve::try_single_phase(
                identity,
                SolidLiquidPhase::SolidLiquid,
                single_phase_knots(0.0),
            ),
            Err(PhaseStateError::InvalidCurve { .. })
        ));
        let mut plateau = single_phase_knots(0.0);
        plateau[1].temperature_k = plateau[0].temperature_k;
        assert!(matches!(
            EquilibriumEnthalpyPhaseCurve::try_single_phase(
                identity,
                SolidLiquidPhase::Solid,
                plateau
            ),
            Err(PhaseStateError::InvalidCurve { .. })
        ));
        let mut cooling = single_phase_knots(1.0);
        cooling[1].temperature_k = 299.0;
        assert!(matches!(
            EquilibriumEnthalpyPhaseCurve::try_single_phase(
                identity,
                SolidLiquidPhase::Liquid,
                cooling
            ),
            Err(PhaseStateError::InvalidCurve { .. })
        ));
        let curve = EquilibriumEnthalpyPhaseCurve::try_single_phase(
            identity,
            SolidLiquidPhase::Solid,
            single_phase_knots(0.0),
        )
        .unwrap();
        assert!(matches!(
            curve.state_at_specific_enthalpy(9.0),
            Err(PhaseStateError::OutsideEnthalpyDomain { .. })
        ));
    }

    #[test]
    fn g0_heat_capacity_curve_integrates_linear_capacity_and_preserves_endpoints() {
        let source = heat_capacity_knots(1_000.0, 2_000.0);
        let curve = EquilibriumEnthalpyPhaseCurve::try_from_heat_capacity(
            ContentHash([0x71; 32]),
            SolidLiquidPhase::Liquid,
            0.0,
            &source,
            1.0,
        )
        .unwrap();
        assert_eq!(
            curve.knots()[0].temperature_k.to_bits(),
            300.0_f64.to_bits()
        );
        assert_eq!(
            curve.knots()[0].specific_enthalpy_j_kg.to_bits(),
            0.0_f64.to_bits()
        );
        let last = curve.knots().last().unwrap();
        assert_eq!(last.temperature_k.to_bits(), 400.0_f64.to_bits());
        assert_eq!(last.bulk_density_kg_m3.to_bits(), 900.0_f64.to_bits());
        assert!((last.specific_enthalpy_j_kg - 150_000.0).abs() < 1e-10);

        // Exact h(350 K) = 1_000*50 + 0.5*10*50^2. The source chord error
        // is bounded by 1 J/kg, so inverse temperature error is at most that
        // divided by the minimum source Cp of 1_000 J/(kg K).
        let midpoint = curve.state_at_specific_enthalpy(62_500.0).unwrap();
        assert_eq!(midpoint.phase(), SolidLiquidPhase::Liquid);
        assert!((midpoint.temperature_k() - 350.0).abs() <= 0.001);
        // Include points between generated knots, not only the exact center.
        for temperature in [333.3, 347.25, 371.25] {
            let delta = temperature - 300.0;
            let exact_h = 1_000.0 * delta + 5.0 * delta * delta;
            let state = curve.state_at_specific_enthalpy(exact_h).unwrap();
            assert!((state.temperature_k() - temperature).abs() <= 0.001);
        }
    }

    #[test]
    fn g0_heat_capacity_curve_handles_decreasing_and_constant_capacity() {
        let decreasing = heat_capacity_knots(2_000.0, 1_000.0);
        let curve = EquilibriumEnthalpyPhaseCurve::try_from_heat_capacity(
            ContentHash([0x72; 32]),
            SolidLiquidPhase::Solid,
            20.0,
            &decreasing,
            1.0,
        )
        .unwrap();
        assert!((curve.knots().last().unwrap().specific_enthalpy_j_kg - 150_020.0).abs() < 1e-10);
        assert!(
            (curve
                .state_at_specific_enthalpy(87_520.0)
                .unwrap()
                .temperature_k()
                - 350.0)
                .abs()
                <= 0.001
        );

        let constant = heat_capacity_knots(1_000.0, 1_000.0);
        let constant_curve = EquilibriumEnthalpyPhaseCurve::try_from_heat_capacity(
            ContentHash([0x73; 32]),
            SolidLiquidPhase::Solid,
            10.0,
            &constant,
            1.0e-12,
        )
        .unwrap();
        assert_eq!(constant_curve.knots().len(), 2);
        assert_eq!(
            constant_curve.knots()[1].specific_enthalpy_j_kg.to_bits(),
            100_010.0_f64.to_bits()
        );
    }

    #[test]
    fn g0_heat_capacity_curve_refuses_invalid_source_phase_tolerance_and_budget() {
        let source = heat_capacity_knots(1_000.0, 2_000.0);
        let refuse = |phase, reference, knots: &[HeatCapacityKnot], tolerance| {
            EquilibriumEnthalpyPhaseCurve::try_from_heat_capacity(
                ContentHash([0x74; 32]),
                phase,
                reference,
                knots,
                tolerance,
            )
        };
        assert!(refuse(SolidLiquidPhase::SolidLiquid, 0.0, &source, 1.0).is_err());
        assert!(refuse(SolidLiquidPhase::Solid, f64::NAN, &source, 1.0).is_err());
        assert!(refuse(SolidLiquidPhase::Solid, 0.0, &source, 0.0).is_err());
        assert!(refuse(SolidLiquidPhase::Solid, 0.0, &source, f64::INFINITY).is_err());

        let mut nonpositive_cp = source;
        nonpositive_cp[1].specific_heat_capacity_j_kg_k = 0.0;
        assert!(refuse(SolidLiquidPhase::Solid, 0.0, &nonpositive_cp, 1.0).is_err());
        let mut unordered = source;
        unordered[1].temperature_k = unordered[0].temperature_k;
        assert!(refuse(SolidLiquidPhase::Solid, 0.0, &unordered, 1.0).is_err());
        let mut nonfinite = source;
        nonfinite[1].bulk_density_kg_m3 = f64::NAN;
        assert!(refuse(SolidLiquidPhase::Solid, 0.0, &nonfinite, 1.0).is_err());
        assert!(
            refuse(
                SolidLiquidPhase::Solid,
                1.0e300,
                &heat_capacity_knots(1.0, 1.0),
                1.0
            )
            .is_err()
        );

        let over_budget = [
            HeatCapacityKnot {
                temperature_k: 300.0,
                specific_heat_capacity_j_kg_k: 1.0,
                bulk_density_kg_m3: 1_000.0,
            },
            HeatCapacityKnot {
                temperature_k: 301.0,
                specific_heat_capacity_j_kg_k: 1.0e9,
                bulk_density_kg_m3: 999.0,
            },
        ];
        assert!(refuse(SolidLiquidPhase::Solid, 0.0, &over_budget, 1.0).is_err());

        let too_many = vec![
            HeatCapacityKnot {
                temperature_k: f64::NAN,
                specific_heat_capacity_j_kg_k: 1.0,
                bulk_density_kg_m3: 1.0,
            };
            MAX_HEAT_CAPACITY_ENTHALPY_KNOTS + 1
        ];
        assert!(refuse(SolidLiquidPhase::Solid, 0.0, &too_many, 1.0).is_err());
    }
}
