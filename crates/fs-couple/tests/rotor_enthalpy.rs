//! G1/G3/G4: source-owned rotor heat, thermal feedback and independent work.
#[path = "circuit_enthalpy/common.rs"]
mod common;

use fs_couple::circuit_enthalpy::{material::MaterialCell, rotor::*};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialCard, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, SurfaceSpec, SystemContext, UncertaintyModel,
};
use fs_material::{
    phase::ResolvedHeatCapacityCurve,
    state_point::{
        LINEAR_INTERFACE_CONSTANT_LOSS_REFERENCE_OMEGA, MaterialPropertySelection,
        ResolvedLinearInterface, ScalarAdmissibility, ScalarPropertyRequirement,
        constant_loss_reference_omega_dims, interface_stiffness_per_area_dims,
        interface_viscosity_per_area_dims, resolve_interface_state_point,
    },
    thermal::source::ThermalMaterialPoint,
};
use fs_mbd::rotordynamics::{RotorBearingModel, transient::RotorTransientState};
use fs_qty::{QuantitySpec, Temperature};
use fs_thermal::enthalpy::*;
use std::cell::Cell;

const A: CellId = CellId(10);
const B: CellId = CellId(20);

struct Fixture {
    card: MaterialCard,
    caloric: ResolvedHeatCapacityCurve,
    interface: InterfaceSystemCard,
    source: ResolvedLinearInterface,
}
fn fixture(c: f64, slope: f64, k: f64, upper: f64, anchor: bool) -> Fixture {
    let (card, caloric) = common::material_with_conductivity(1.0, 1.0, 0.0, 0.0);
    let mut properties = vec![
        ("k", interface_stiffness_per_area_dims(), [k, k]),
        (
            "c",
            interface_viscosity_per_area_dims(),
            [c - 0.5 * slope, c + slope * (upper - 300.0)],
        ),
    ];
    if anchor {
        properties.push((
            LINEAR_INTERFACE_CONSTANT_LOSS_REFERENCE_OMEGA,
            constant_loss_reference_omega_dims(),
            [1.0, 1.0],
        ));
    }
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (name, dims, values) in properties {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(299.5, values[0]), (upper, values[1])],
                    dims,
                },
                validity: ValidityDomain::unconstrained()
                    .with_quantity(
                        "T",
                        QuantitySpec::dimensional(Temperature::DIMS),
                        299.5,
                        upper,
                    )
                    .with("omega", 0.0, 1e4)
                    .with("lot", 7.0, 7.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic rotor thermal fixture, not measured bearing data".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
        requirements.push(
            ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::NonNegative)
                .unwrap(),
        );
    }
    let surface = SurfaceSpec {
        material: card.id().clone(),
        texture_frame: "rotor".into(),
    };
    let interface = InterfaceSystemCard::assemble(
        surface.clone(),
        surface,
        SystemContext {
            medium: "fixture".into(),
            third_body: None,
            environment: "fixed geometry".into(),
            history: "memoryless viscosity".into(),
        },
        claims,
        vec![],
    )
    .unwrap();
    let point = common::point(300.0).with("omega", 0.0).unwrap();
    let resolved = resolve_interface_state_point(
        &interface,
        &point,
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let source = ResolvedLinearInterface::new(resolved, "k", "c", 1.0).unwrap();
    Fixture {
        card,
        caloric,
        interface,
        source,
    }
}
impl Fixture {
    fn thermal(&self) -> DepositionSimulation<'_> {
        let events = [1.0, 3.0]
            .into_iter()
            .enumerate()
            .map(|(i, mass)| TimedEvent {
                time_s: 0.0,
                event: Event::Birth(DepositedCell {
                    id: [A, B][i],
                    segment_id: i as u64,
                    body_id: i as u64,
                    material_id: MaterialId(0),
                    mass_kg: mass,
                    specific_h_j_kg: 100.0,
                    geometry: CellGeometry {
                        start_m: [i as f64, 0.0, 0.0],
                        end_m: [i as f64 + mass, 0.0, 0.0],
                        volume_m3: mass,
                        surface_area_m2: 4.0 * mass + 2.0,
                    },
                    boundary: common::boundary(),
                }),
            })
            .collect();
        DepositionSimulation::new(
            events,
            &[(MaterialId(0), self.caloric.curve() as &dyn ThermalLaw)],
            DepositionBudget {
                max_cells: 2,
                max_links: 0,
                max_events: 2,
                max_substeps: 10,
            },
        )
        .unwrap()
        .advance_to(
            0.0,
            StepConfig {
                max_dt_s: 1.0,
                energy_tolerance_j: 1e-9,
            },
            || Ok(()),
        )
        .unwrap()
        .0
    }
    fn cells<'a>(&'a self, context: &'a QueryPoint) -> [MaterialCell<'a>; 2] {
        [A, B].map(|cell| MaterialCell {
            cell,
            material_id: MaterialId(0),
            curve: self.caloric.curve(),
            point: ThermalMaterialPoint {
                material: &self.card,
                context,
                temperature_axis: "T",
                liquid_fraction_axis: None,
            },
            conductivity: None,
        })
    }
    fn damper<'a>(&'a self, row: &'a [f64]) -> MaterialRotorDamper<'a> {
        MaterialRotorDamper {
            interface: &self.interface,
            source: &self.source,
            slip_row: row,
            cells: [A, B],
            temperature_side: 0,
            first_heat_fraction: 0.25,
        }
    }
}
fn network<'a>(
    model: &'a RotorBearingModel,
    cells: &'a [MaterialCell<'a>],
    dampers: &'a [MaterialRotorDamper<'a>],
) -> RotorThermalNetwork<'a> {
    RotorThermalNetwork {
        baseline: model,
        spin_rad_s: 0.0,
        cells,
        links: &[],
        dampers,
        qualified_angular_frequency_rad_s: 100.0,
    }
}
fn controls() -> RotorThermalControls {
    let mut solve = common::config();
    solve.difference_step_j = 1e-5;
    solve.newton.absolute_tolerance = 1e-11;
    solve.cell_energy_tolerance_j = 1e-8;
    solve.thermal_energy_tolerance_j = 1e-8;
    solve.total_energy_tolerance_j = 1e-8;
    RotorThermalControls {
        solve,
        relative_motion_tolerance: 1e-12,
        maximum_step_s: 1.0,
    }
}
fn scalar_model(c: f64) -> RotorBearingModel {
    RotorBearingModel::new(1, &[1.0], &[c], &[0.0], &[0.0], || false).unwrap()
}
fn state(v: f64) -> RotorTransientState {
    RotorTransientState {
        displacement: vec![0.0],
        velocity: vec![v],
    }
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 2e-7, "{a:e} != {b:e}");
}

#[test]
fn endpoint_material_feedback_matches_exact_discrete_solution_and_continues() {
    let f = fixture(0.25, 0.25, 0.0, 350.0, false);
    let thermal = f.thermal();
    let context = common::point(300.0);
    let cells = f.cells(&context);
    let model = scalar_model(0.0);
    let dampers = [f.damper(&[2.0])];
    let binding = network(&model, &cells, &dampers);
    let rotor = binding.bind(&thermal, state(3.0), || false).unwrap();
    let candidate = rotor
        .step_candidate(&thermal, 0.5, &[0.0], controls(), &[0.0; 2], || false)
        .unwrap();
    // c(T)=.25*(1+T-300), b=2, capacities1/3 and fractions1/4,3/4.
    // The exact coupled midpoint endpoint is c=.5, v1=1 and Q=4 J.
    near(candidate.rotor().state().velocity[0], 1.0);
    near(candidate.rotor().sources()[0].damping_n_s_m(), 0.5);
    near(candidate.report().electrical.damper_heat_j[0], 4.0);
    near(candidate.report().deposited_joule_j, 4.0);
    near(candidate.thermal().temperature_k(A).unwrap(), 301.0);
    near(candidate.thermal().temperature_k(B).unwrap(), 301.0);
    near(
        candidate.report().electrical.motion.kinetic_energy_change_j,
        -4.0,
    );
    assert!(candidate.report().total_energy_residual_j.abs() < 1e-8);
    let (thermal, rotor, _) = candidate.into_parts();
    let next = rotor
        .step_candidate(&thermal, 0.6, &[0.0], controls(), &[0.0; 2], || false)
        .unwrap();
    assert!(next.rotor().state().velocity[0] < 1.0);
    let t = next.thermal().temperature_k(A).unwrap();
    assert!(t > 301.0);
    near(
        next.rotor().sources()[0].damping_n_s_m(),
        0.25 * (1.0 + t - 300.0),
    );
}

#[test]
fn active_original_operator_can_supply_heat_with_zero_net_damping() {
    let f = fixture(0.5, 0.0, 0.0, 350.0, false);
    let thermal = f.thermal();
    let context = common::point(300.0);
    let cells = f.cells(&context);
    let model = scalar_model(-2.0);
    let dampers = [f.damper(&[2.0])];
    let binding = network(&model, &cells, &dampers);
    let rotor = binding.bind(&thermal, state(2.0), || false).unwrap();
    let result = rotor
        .step_candidate(&thermal, 0.25, &[0.0], controls(), &[0.0; 2], || false)
        .unwrap();
    let work = &result.report().electrical;
    near(work.motion.state_after.velocity[0], 2.0);
    near(work.motion.velocity_force_work_j, 0.0);
    near(work.motion.kinetic_energy_change_j, 0.0);
    near(work.baseline_velocity_work_j, -2.0);
    near(work.damper_heat_j[0], 2.0);
    near(work.supplied_to_storage_and_heat_j, 2.0);
    near(result.thermal().temperature_k(A).unwrap(), 300.5);
    near(result.thermal().temperature_k(B).unwrap(), 300.5);
}

#[test]
fn full_mass_gyro_and_skew_stiffness_keep_their_original_work() {
    let f = fixture(0.5, 0.0, 0.0, 350.0, false);
    let thermal = f.thermal();
    let context = common::point(300.0);
    let cells = f.cells(&context);
    let model = RotorBearingModel::new(
        2,
        &[2.0, 0.5, 0.5, 1.5],
        &[-1.0, 2.0, -2.0, 1.0],
        &[2.0, 3.0, -1.0, 2.0],
        &[0.0, -1.0, 1.0, 0.0],
        || false,
    )
    .unwrap();
    let dampers = [f.damper(&[1.0, -2.0])];
    let mut binding = network(&model, &cells, &dampers);
    binding.spin_rad_s = 2.0;
    let initial = RotorTransientState {
        displacement: vec![1.0, -1.0],
        velocity: vec![1.0, 0.0],
    };
    let rotor = binding.bind(&thermal, initial, || false).unwrap();
    // Manufactured binary-exact midpoint: v1=[0,1], q1=[1.125,-.875].
    let result = rotor
        .step_candidate(
            &thermal,
            0.25,
            &[-7.4375, 2.0625],
            controls(),
            &[0.0; 2],
            || false,
        )
        .unwrap();
    let work = &result.report().electrical;
    near(work.motion.state_after.displacement[0], 1.125);
    near(work.motion.state_after.displacement[1], -0.875);
    near(work.motion.state_after.velocity[0], 0.0);
    near(work.motion.state_after.velocity[1], 1.0);
    near(work.motion.kinetic_energy_change_j, -0.25);
    near(work.elastic_energy_change_j, 0.046875);
    near(work.nonconservative_stiffness_work_j, -0.5);
    near(work.baseline_velocity_work_j, 0.0);
    near(work.gyroscopic_work_j, 0.0);
    near(work.damper_heat_j[0], 0.03125);
    near(work.supplied_to_storage_and_heat_j, -0.171875);
}

#[test]
fn independent_channels_share_receivers_without_losing_heat_partition() {
    let f = fixture(0.5, 0.0, 0.0, 350.0, false);
    let thermal = f.thermal();
    let context = common::point(300.0);
    let cells = f.cells(&context);
    let model = scalar_model(0.0);
    let mut second = f.damper(&[2.0]);
    second.first_heat_fraction = 0.75;
    second.temperature_side = 1;
    let dampers = [f.damper(&[1.0]), second];
    let binding = network(&model, &cells, &dampers);
    let rotor = binding.bind(&thermal, state(2.0), || false).unwrap();
    let result = rotor
        .step_candidate(&thermal, 0.2, &[0.0], controls(), &[0.0; 2], || false)
        .unwrap();
    near(result.rotor().state().velocity[0], 1.2);
    near(result.report().electrical.damper_heat_j[0], 0.256);
    near(result.report().electrical.damper_heat_j[1], 1.024);
    near(result.report().deposited_joule_j, 1.28);
    near(result.thermal().temperature_k(A).unwrap(), 300.832);
    near(
        result.thermal().temperature_k(B).unwrap(),
        300.0 + 0.448 / 3.0,
    );
}

#[test]
fn source_domain_late_cancellation_and_nonviscous_laws_publish_no_half_state() {
    let model = scalar_model(0.0);
    for (k, anchor) in [(1.0, false), (0.0, true)] {
        let f = fixture(0.5, 0.0, k, 350.0, anchor);
        let thermal = f.thermal();
        let context = common::point(300.0);
        let cells = f.cells(&context);
        let dampers = [f.damper(&[2.0])];
        let binding = network(&model, &cells, &dampers);
        assert!(binding.bind(&thermal, state(2.0), || false).is_err());
    }
    for upper in [300.05, 350.0] {
        let f = fixture(0.5, 0.0, 0.0, upper, false);
        let thermal = f.thermal();
        let saved_cells = thermal.cells().clone();
        let context = common::point(300.0);
        let cells = f.cells(&context);
        let dampers = [f.damper(&[2.0])];
        let binding = network(&model, &cells, &dampers);
        let rotor = binding.bind(&thermal, state(2.0), || false).unwrap();
        let saved_state = rotor.state().clone();
        let saved_sources = rotor.sources().to_vec();
        let calls = Cell::new(0_usize);
        let result = rotor.step_candidate(&thermal, 0.25, &[0.0], controls(), &[0.0; 2], || {
            calls.set(calls.get() + 1);
            false
        });
        if upper < 301.0 {
            assert!(result.is_err());
        } else {
            assert!(result.is_ok());
            let full = calls.get();
            calls.set(0);
            assert!(
                rotor
                    .step_candidate(&thermal, 0.25, &[0.0], controls(), &[0.0; 2], || {
                        calls.set(calls.get() + 1);
                        calls.get() >= full.saturating_sub(1)
                    })
                    .is_err()
            );
            let unsupported = RotorThermalNetwork {
                qualified_angular_frequency_rad_s: 2e4,
                ..binding
            };
            assert!(unsupported.bind(&thermal, state(2.0), || false).is_err());
        }
        assert_eq!(thermal.cells(), &saved_cells);
        assert_eq!(rotor.state(), &saved_state);
        assert_eq!(rotor.sources(), saved_sources.as_slice());
        assert_eq!(thermal.time_s(), 0.0);
    }
}
