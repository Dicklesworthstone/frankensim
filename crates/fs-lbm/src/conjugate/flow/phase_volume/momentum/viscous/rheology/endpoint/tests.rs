//! Endpoint rates must belong to the accepted native MAC field, not scattered P.

use super::super::super::super::tests::{close, controls as momentum_controls};
use super::*;
use crate::conjugate::flow::CoexistenceVolumeInput;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::fluid_transport::FluidTransport;
use fs_material::liquid::{KINEMATIC_VISCOSITY_DIMS, KINEMATIC_VISCOSITY_PROPERTY as NU_KEY};
use fs_material::liquid_temperature::LiquidViscosityBasis;
use fs_material::phase::HeatCapacitySourceQuery;
use fs_material::phase::melting::{MeltingSourceQuery, ResolvedMeltingCurve};
use fs_material::phase::mushy::MushyPhaseLaw;
use fs_material::phase::transport::PhaseTransportSourceQuery;
use fs_material::phase::transport::resolve_phase_transport_state_with_shear_for_volume_change;
use fs_material::state_point::{
    MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS,
    SPECIFIC_HEAT_CAPACITY_PROPERTY as CP_KEY, THERMAL_CONDUCTIVITY_PROPERTY as K_KEY,
};
use fs_material::thermal::source::ThermalMaterialPoint;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Dims, Frequency, QuantitySpec, Temperature};
use std::cell::Cell;

fn absolute_temperature() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}

fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", absolute_temperature(), t)
        .unwrap()
        .with_quantity("rate", QuantitySpec::dimensional(Frequency::DIMS), 1.0)
        .unwrap()
}

#[derive(Clone, Copy)]
enum Rheology {
    Thinning,
    Thickening,
    ShearAndBulk,
}

/// Real phase/material resolution with finite rate support 0 <= s <= 5/4,
/// where s=sqrt(3)*gamma/2. rho=2, Cp=1 and T=300+h on the liquid branch.
fn fixture(rheology: Rheology) -> (MaterialCard, ResolvedMeltingCurve) {
    let cp = QuantitySpec::dimensional(SPECIFIC_HEAT_CAPACITY_DIMS);
    let rho = QuantitySpec::dimensional(Density::DIMS);
    let latent = QuantitySpec::dimensional(Dims([2, 0, -2, 0, 0, 0]));
    let nu = QuantitySpec::dimensional(KINEMATIC_VISCOSITY_DIMS);
    let k = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::ThermalConductivity,
        ValueForm::Static,
    ));
    let maximum_rate = 2.5 / fs_math::det::sqrt(3.0);
    let viscosity = match rheology {
        Rheology::Thickening => [1.0 / 4.0, 9.0 / 16.0], // mu=1/2+s/2.
        Rheology::Thinning => [7.0 / 16.0, 9.0 / 32.0],  // mu=7/8-s/4.
        Rheology::ShearAndBulk => [1.0 / 8.0, 9.0 / 32.0], // mu=1/4+s/4.
    }; // All original viscosity sources carry nu=mu/rho.

    let mut claims = ClaimSet::new();
    for (name, quantity, axis, range, values) in [
        ("solid_cp", cp, "T", [299.0, 300.0], [1.0; 2]),
        ("solid_density", rho, "T", [299.0, 300.0], [4.0; 2]),
        (CP_KEY, cp, "T", [300.0, 302.0], [1.0; 2]),
        ("density", rho, "T", [300.0, 302.0], [2.0; 2]),
        (NU_KEY, nu, "rate", [0.0, maximum_rate], viscosity),
        (K_KEY, k, "T", [299.0, 302.0], [1.0; 2]),
        (
            "melting_point",
            absolute_temperature(),
            "T",
            [299.0, 302.0],
            [300.0; 2],
        ),
        ("fusion", latent, "T", [299.0, 302.0], [1.0; 2]),
    ] {
        let temperatures = if axis == "T" { range } else { [299.0, 302.0] };
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::with_quantity(name, quantity),
                value: PropertyValue::Curve {
                    abscissa: axis.into(),
                    abscissa_dims: if axis == "T" {
                        Temperature::DIMS
                    } else {
                        Frequency::DIMS
                    },
                    knots: vec![(range[0], values[0]), (range[1], values[1])],
                    dims: quantity.dims(),
                },
                validity: ValidityDomain::unconstrained()
                    .with_quantity(
                        "T",
                        absolute_temperature(),
                        temperatures[0],
                        temperatures[1],
                    )
                    .with_quantity(
                        "rate",
                        QuantitySpec::dimensional(Frequency::DIMS),
                        0.0,
                        maximum_rate,
                    ),
                interpolation: InterpolationPolicy::LinearInside,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source:
                        "MANUFACTURED finite-support endpoint-rheology oracle, not measurements"
                            .into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let material = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured pure substance".into(),
            phase: "solid/liquid".into(),
            process: "native endpoint rheology regression".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let solid_cp = PropertyKey::with_quantity("solid_cp", cp);
    let solid_rho = PropertyKey::with_quantity("solid_density", rho);
    let liquid_cp = PropertyKey::with_quantity(CP_KEY, cp);
    let liquid_rho = PropertyKey::with_quantity("density", rho);
    let tm = PropertyKey::with_quantity("melting_point", absolute_temperature());
    let fusion = PropertyKey::with_quantity("fusion", latent);
    let points = [point(299.0), point(300.0), point(302.0)];
    let melting = ResolvedMeltingCurve::try_new(
        &material,
        MeltingSourceQuery {
            solid: HeatCapacitySourceQuery {
                heat_capacity: &solid_cp,
                density: &solid_rho,
                temperature_axis: "T",
                lower: &points[0],
                upper: &points[1],
                selection: SelectionPolicy::SingleClaimOnly,
            },
            liquid: HeatCapacitySourceQuery {
                heat_capacity: &liquid_cp,
                density: &liquid_rho,
                temperature_axis: "T",
                lower: &points[1],
                upper: &points[2],
                selection: SelectionPolicy::SingleClaimOnly,
            },
            melting_point: &tm,
            latent_heat_fusion: &fusion,
            transition: &points[1],
            transition_selection: SelectionPolicy::SingleClaimOnly,
            phase_selector: None,
        },
        -2.0,
        1e-12,
    )
    .unwrap();
    (material, melting)
}

fn resolve(
    material: &MaterialCard,
    law: &MushyPhaseLaw<'_>,
    h: f64,
    rate: f64,
    rheology: Rheology,
) -> Result<(ResolvedPhaseTransportState, f64), ChtError> {
    let source = resolve_phase_transport_state_with_shear_for_volume_change(
        law,
        h,
        PhaseTransportSourceQuery {
            source: ThermalMaterialPoint {
                material,
                context: &point(301.5),
                temperature_axis: "T",
                liquid_fraction_axis: None,
            },
            viscosity_basis: LiquidViscosityBasis::Kinematic,
            extra_requirements: &[],
            selection: MaterialPropertySelection::SingleClaimOnly,
        },
        "rate",
        rate,
    )
    .map_err(|error| ChtError::InvalidInput {
        field: "test.endpoint_rheology_source",
        reason: error.to_string(),
    })?;
    let bulk = match rheology {
        Rheology::ShearAndBulk => 0.25 + 0.25 * fs_math::det::sqrt(3.0) * rate,
        Rheology::Thinning | Rheology::Thickening => 0.0,
    };
    Ok((source, bulk)) // Explicit independent volumetric response.
}

fn controls() -> ThermalVolumeMomentumControls {
    let mut controls = momentum_controls();
    controls.mechanics.projection.duration_s = 1.0;
    controls.volume.thermal.volume.transport.projection = controls.mechanics.projection;
    controls
}

fn viscosity_controls() -> MaterialViscosityControls {
    MaterialViscosityControls {
        max_cells: 100,
        max_sweeps: 1000,
        relative_tolerance: 1e-12,
        maximum_impulse_defect_kg_m_s: 1e-9,
        maximum_angular_impulse_defect_kg_m2_s: 1e-9,
        maximum_energy_defect_j: 1e-8,
    }
}

fn input<'a>(
    mass: &'a [f64],
    phases: &'a [Option<EquilibriumPhaseState>],
    momentum: &'a [[f64; 3]],
    heat: &'a [f64],
    force: &'a [[f64; 3]],
) -> ThermalVolumeMomentumInput<'a> {
    ThermalVolumeMomentumInput {
        volume: CoexistenceVolumeInput {
            mass_kg: mass,
            phase_states: phases,
            net_enthalpy_increment_j: heat,
            boundary_inflow: [None; 6],
        },
        momentum_kg_m_s: momentum,
        acceleration_m_s2: [0.0; 3],
        force_density_n_m3: force,
        boundary_velocity_m_s: [None; 6],
    }
}

#[test]
fn real_shear_and_bulk_sources_close_on_native_endpoint_rate() {
    for rheology in [
        Rheology::Thinning,
        Rheology::Thickening,
        Rheology::ShearAndBulk,
    ] {
        let (material, melting) = fixture(rheology);
        let caloric = melting.curve();
        let law = MushyPhaseLaw::try_from_melting_curve(&melting, 1e-8).unwrap();
        let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
        let field = FlowField::quiescent(&domain);
        let phase = caloric.state_at_specific_enthalpy(1.0).unwrap();
        let mass = [2.0; 2];
        let phases = [Some(phase); 2];
        let momentum = [[2.0, 0.0, 0.0]; 2];
        let mut queries = vec![];
        let result = field
            .advance_piecewise_thermal_volume_momentum_with_endpoint_rheology(
                &domain,
                None,
                caloric,
                input(&mass, &phases, &momentum, &[0.0; 2], &[[0.0; 3]; 2]),
                [FvBoundary::Symmetry; 6],
                controls(),
                viscosity_controls(),
                &CancelGate::new_clock_free(),
                |cell, supplied_phase, rate| {
                    assert_eq!(supplied_phase, &phase);
                    queries.push((cell, rate));
                    resolve(
                        &material,
                        &law,
                        supplied_phase.specific_enthalpy_j_kg(),
                        rate,
                        rheology,
                    )
                },
            )
            .unwrap();
        // Native faces are [0,a,0], initially a0=1. Each actual face mass is 2.
        // gamma=2|a|/sqrt(3), so s=a for the admitted positive branch.
        // M*(a-a0)/dt +2*(4*mu/3+zeta)*a=0.
        // Thickening: 2*a^2+5*a-3=0. Thinning: 2*a^2-13*a+6=0.
        // Joint law: mu=1/4+s/4, zeta=1/4+s/2, hence 10*a^2+19*a-12=0.
        // All have the unique supported root a=1/2 and longitudinal coefficient 4*mu/3+zeta=1.
        // Frozen original coefficients give a=3/7, 6/11, or 12/29.
        let (expected_mu, expected_bulk) = match rheology {
            Rheology::ShearAndBulk => (3.0 / 8.0, 1.0 / 2.0),
            Rheology::Thinning | Rheology::Thickening => (3.0 / 4.0, 0.0),
        };
        let rate = 1.0 / fs_math::det::sqrt(3.0);
        let viscous = &result.evolution.viscosity;
        let volume = &result.evolution.evolution.thermal.volume;
        assert!(viscous.sweeps > 1);
        assert_eq!(result.source_updates, viscous.sweeps + 1);
        assert_eq!(result.source_evaluations, 2 * result.source_updates);
        assert_eq!(queries.len(), result.source_evaluations);
        for cell in 0..2 {
            assert_eq!(queries[cell].0, cell);
            close(queries[cell].1, 2.0 * rate);
            close(result.endpoint_shear_rate_per_s[cell], rate);
            let source = result.endpoint_sources[cell].as_ref().unwrap();
            assert_eq!(source.phase_state(), &phase);
            close(
                source.transport_state().unwrap().dynamic_viscosity_pa_s(),
                expected_mu,
            );
            close(
                source.material().property(NU_KEY).unwrap().value_si(),
                expected_mu / 2.0,
            );
            for property in source.material().properties() {
                let coordinates = &property.answer().receipt.query_point;
                close(
                    coordinates
                        .iter()
                        .find(|(axis, _)| axis == "rate")
                        .unwrap()
                        .1,
                    rate,
                );
                assert!(coordinates.contains(&("T".into(), 301.0)));
            }
            close(result.bulk_viscosity_pa_s[cell], expected_bulk);
            close(viscous.cell_molecular_heat_j[cell], 1.0 / 4.0);
            close(
                volume.phase_states[cell].unwrap().specific_enthalpy_j_kg(),
                9.0 / 8.0,
            );
            close(volume.transport.mass_kg[cell], 2.0);
        }
        close(viscous.molecular_heat_j, 1.0 / 2.0);
        close(viscous.native_kinetic_energy_change_j, -3.0 / 4.0);
        close(viscous.backward_euler_loss_j, 1.0 / 4.0);
        close(viscous.numerical_scatter_loss_j, 1.0 / 8.0);
        close(viscous.cell_kinetic_energy_change_j, -7.0 / 8.0);
        close(viscous.wall_work_j, 0.0);
        close(viscous.representation_work_j, 0.0);
        close(viscous.face_impulse_kg_m_s[0], -1.0);
        close(viscous.cell_momentum_change_kg_m_s[0], -1.0);
        close(viscous.energy_defect_j, 0.0);
        close(volume.report.material_enthalpy_change_j, 1.0 / 2.0);
        assert!(viscous.relative_residual <= viscosity_controls().relative_tolerance);
        assert!(viscous.absolute_impulse_defect_kg_m_s <= 1e-9);
        assert!(viscous.absolute_energy_defect_j <= 1e-8);
        // Native scatter would give P=3/2 per cell, hence reconstructed a=3/4.
        // That field would falsely report gamma=sqrt(3)/2, not this endpoint.
        assert!((result.endpoint_shear_rate_per_s[0] - 1.5 * rate).abs() > 0.1);
        assert_eq!(field, FlowField::quiescent(&domain));
        assert_eq!(phases, [Some(phase); 2]);
        assert_eq!(momentum, [[2.0, 0.0, 0.0]; 2]);
    }
}

#[test]
fn endpoint_work_source_and_cancellation_refusals_preserve_original_retry() {
    let (material, melting) = fixture(Rheology::Thickening);
    let caloric = melting.curve();
    let law = MushyPhaseLaw::try_from_melting_curve(&melting, 1e-8).unwrap();
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let field = FlowField::quiescent(&domain);
    let phase = caloric.state_at_specific_enthalpy(1.0).unwrap();
    let mass = [2.0; 2];
    let phases = [Some(phase); 2];
    let momentum = [[2.0, 0.0, 0.0]; 2];
    #[derive(Clone, Copy, PartialEq)]
    enum Fault {
        None,
        StalePhase,
        CancelEndpoint,
    }
    let run = |viscosity, momenta: &[[f64; 3]], fault, calls: &Cell<usize>, gate: &CancelGate| {
        field.advance_piecewise_thermal_volume_momentum_with_endpoint_rheology(
            &domain,
            None,
            caloric,
            input(&mass, &phases, momenta, &[0.0; 2], &[[0.0; 3]; 2]),
            [FvBoundary::Symmetry; 6],
            controls(),
            viscosity,
            gate,
            |_, supplied_phase, rate| {
                calls.set(calls.get() + 1);
                assert_eq!(supplied_phase, &phase);
                if calls.get() > 2 && fault == Fault::CancelEndpoint {
                    gate.request();
                    return Err(ChtError::InvalidInput {
                        field: "test.endpoint_cancel",
                        reason: "resolver failed after cancellation".into(),
                    });
                }
                let h = supplied_phase.specific_enthalpy_j_kg()
                    + if calls.get() > 2 && fault == Fault::StalePhase {
                        0.001
                    } else {
                        0.0
                    };
                resolve(&material, &law, h, rate, Rheology::Thickening)
            },
        )
    };
    let gate = CancelGate::new_clock_free();
    let calls = Cell::new(0);
    let expected = run(viscosity_controls(), &momentum, Fault::None, &calls, &gate).unwrap();
    calls.set(0);
    let exhausted = run(
        MaterialViscosityControls {
            max_sweeps: 1,
            ..viscosity_controls()
        },
        &momentum,
        Fault::None,
        &calls,
        &gate,
    );
    assert!(matches!(
        exhausted,
        Err(ChtError::SolverNotConverged { .. })
    ));
    assert_eq!(calls.get(), 4); // Initial and ONE final native field, two cells each.
    calls.set(0);
    assert!(
        run(
            viscosity_controls(),
            &momentum,
            Fault::StalePhase,
            &calls,
            &gate
        )
        .is_err()
    );
    assert_eq!(calls.get(), 3);
    // Genuine source-domain refusal: original gamma=4/sqrt(3) exceeds the
    // declared 5/(2*sqrt(3)) upper limit; no clipping or extrapolation.
    calls.set(0);
    assert!(
        run(
            viscosity_controls(),
            &[[4.0, 0.0, 0.0]; 2],
            Fault::None,
            &calls,
            &gate
        )
        .is_err()
    );
    assert_eq!(calls.get(), 1);

    let early_cancel = CancelGate::new_clock_free();
    early_cancel.request();
    calls.set(0);
    assert!(matches!(
        run(
            viscosity_controls(),
            &momentum,
            Fault::None,
            &calls,
            &early_cancel
        ),
        Err(ChtError::Cancelled)
    ));
    assert_eq!(calls.get(), 0);
    let endpoint_cancel = CancelGate::new_clock_free();
    calls.set(0);
    assert!(matches!(
        run(
            viscosity_controls(),
            &momentum,
            Fault::CancelEndpoint,
            &calls,
            &endpoint_cancel
        ),
        Err(ChtError::Cancelled)
    ));
    assert_eq!(calls.get(), 3);
    assert_eq!(field, FlowField::quiescent(&domain));
    assert_eq!(mass, [2.0; 2]);
    assert_eq!(phases, [Some(phase); 2]);
    assert_eq!(momentum, [[2.0, 0.0, 0.0]; 2]);
    calls.set(0);
    let retry = run(viscosity_controls(), &momentum, Fault::None, &calls, &gate).unwrap();
    assert_eq!(retry, expected);
    assert_eq!(calls.get(), retry.source_evaluations);
}

mod sgs;
