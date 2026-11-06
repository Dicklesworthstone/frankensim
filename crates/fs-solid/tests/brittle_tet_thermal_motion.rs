//! Source-integrated thermal strain, physical inertia and native tensile damage.
//! Manufactured coefficients verify this fixed-plane small-strain model;
//! no measured material qualification or independent transverse cracks are implied.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::brittle::source::expansion::ResolvedThermalBrittle;
use fs_material::brittle::source::expansion::dynamics::resolve_thermal_brittle_dynamic_material;
use fs_material::brittle::source::{
    BRITTLE_BK_EXPONENT_PROPERTY, BRITTLE_MODE_I_FRACTURE_ENERGY_PROPERTY,
    BRITTLE_MODE_II_FRACTURE_ENERGY_PROPERTY, BRITTLE_SHEAR_STRENGTH_PROPERTY,
    BRITTLE_TENSILE_STRENGTH_PROPERTY,
};
use fs_material::cohesive::source::FRACTURE_ENERGY_PER_AREA_DIMS;
use fs_material::state_point::{
    DENSITY_PROPERTY, INVERSE_TEMPERATURE_DIMS, LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
    MaterialPropertySelection, POISSON_RATIO_PROPERTY, YOUNG_MODULUS_PROPERTY,
};
use fs_qty::{Density, Dims, Pressure, QuantitySpec, Temperature};
use fs_solid::brittle::TetCrackBandGeometry;
use fs_solid::brittle::solid::dynamics::{BrittleTetDynamicControls, BrittleTetDynamicError};
use fs_solid::brittle::solid::expansion::ThermalBrittleTetSolid;
use fs_solid::brittle::solid::expansion::dynamics::{
    ThermalBrittleTetDynamicEnergy, ThermalBrittleTetMotion,
};
use fs_solid::brittle::solid::{BrittleTetControls, BrittleTetError};
use fs_solid::linear3::TetAssemblyBudget;
use fs_solid::loads3::TetLoadGeometry3;
use fs_solid::plastic3::{IncrementSettings, SmallStrainTetError};

const NODES: [[f64; 3]; 4] = [[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]];
const TETS: [[usize; 4]; 1] = [[0, 1, 2, 3]];
const YOUNG: f64 = 1200.0;
const DENSITY: f64 = 6.0;
const STRENGTH: f64 = 1.2;
const GC: f64 = 0.0012;
const VOLUME: f64 = 1.0 / 6.0;

fn with_cx<R>(run: impl FnOnce(&Cx<'_>, &CancelGate) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 593,
                kernel_id: 593,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        run(&cx, &gate)
    })
}

fn budget() -> TetAssemblyBudget {
    TetAssemblyBudget {
        maximum_nodes: 4,
        maximum_tetrahedra: 1,
        maximum_free_dofs: 12,
        minimum_scaled_jacobian: 1e-10,
    }
}

fn settings() -> IncrementSettings {
    IncrementSettings {
        max_newton_iterations: 30,
        relative_tolerance: 1e-11,
        max_backtracks: 18,
        linear_tolerance: 1e-13,
        max_linear_iterations: 120,
    }
}

fn controls() -> BrittleTetDynamicControls {
    BrittleTetDynamicControls {
        work: BrittleTetControls {
            material_energy_tolerance_j: 1e-10,
            equilibrium_work_tolerance_j: 1e-10,
            maximum_endpoint_work_defect_j: 1.0,
        },
        kinetic_work_tolerance_j: 1e-10,
        maximum_numerical_kinetic_loss_j: 1.0,
    }
}

fn card(tag: &str, alpha_per_k: f64, varying_density: bool) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        (YOUNG_MODULUS_PROPERTY, Pressure::DIMS, YOUNG),
        (POISSON_RATIO_PROPERTY, Dims::NONE, 0.0),
        (DENSITY_PROPERTY, Density::DIMS, DENSITY),
        (BRITTLE_TENSILE_STRENGTH_PROPERTY, Pressure::DIMS, STRENGTH),
        (BRITTLE_SHEAR_STRENGTH_PROPERTY, Pressure::DIMS, 1.5),
        (
            BRITTLE_MODE_I_FRACTURE_ENERGY_PROPERTY,
            FRACTURE_ENERGY_PER_AREA_DIMS,
            GC,
        ),
        (
            BRITTLE_MODE_II_FRACTURE_ENERGY_PROPERTY,
            FRACTURE_ENERGY_PER_AREA_DIMS,
            2.0 * GC,
        ),
        (BRITTLE_BK_EXPONENT_PROPERTY, Dims::NONE, 1.5),
        (
            LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
            INVERSE_TEMPERATURE_DIMS,
            alpha_per_k,
        ),
    ] {
        let varied = varying_density && name == DENSITY_PROPERTY;
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: if varied {
                    PropertyValue::Curve {
                        abscissa: "T".into(),
                        abscissa_dims: Temperature::DIMS,
                        knots: vec![(300.0, DENSITY), (400.0, 0.5 * DENSITY)],
                        dims,
                    }
                } else {
                    PropertyValue::Scalar { value, dims }
                },
                validity: if varied {
                    ValidityDomain::unconstrained().with_quantity(
                        "T",
                        QuantitySpec::dimensional(Temperature::DIMS),
                        300.0,
                        400.0,
                    )
                } else {
                    ValidityDomain::unconstrained()
                },
                interpolation: if varied {
                    InterpolationPolicy::LinearInside
                } else {
                    InterpolationPolicy::ConstantWithinValidity
                },
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: format!("manufactured thermal tensile motion {tag}"),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured brittle solid".into(),
            phase: "solid".into(),
            process: tag.into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn temperature(kelvin: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), kelvin)
        .unwrap()
}

fn source(
    card: &MaterialCard,
    geometry: &TetLoadGeometry3<'_>,
    reference_k: f64,
    current_k: f64,
) -> ResolvedThermalBrittle {
    let band = TetCrackBandGeometry::from_element(geometry, 0, [1.0, 0.0, 0.0]).unwrap();
    resolve_thermal_brittle_dynamic_material(
        card,
        &temperature(reference_k),
        &temperature(current_k),
        MaterialPropertySelection::SingleClaimOnly,
        band.material_settings(0.02, 1e-9),
    )
    .unwrap()
}

fn axial_supports() -> Vec<(usize, f64)> {
    (0..12)
        .filter(|&dof| dof != 3)
        .map(|dof| (dof, 0.0))
        .collect()
}

fn affine_targets(diagonal: f64) -> Vec<(usize, f64)> {
    NODES
        .iter()
        .enumerate()
        .flat_map(|(node, xyz)| (0..3).map(move |axis| (3 * node + axis, diagonal * xyz[axis])))
        .collect()
}

fn uniform_x_force(total_n: f64) -> Vec<(usize, f64)> {
    (0..4).map(|node| (3 * node, total_n / 4.0)).collect()
}

fn near(actual: f64, expected: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= absolute,
        "{actual:e} != {expected:e}"
    );
}

fn quadrature_kinetic(velocity: &[[f64; 3]]) -> f64 {
    // Independent positive degree-two tetrahedral quadrature retains velocity
    // cross terms; replacing consistent mass by row-sum lumping fails this check.
    let a = 0.5854101966249685;
    let b = 0.1381966011250105;
    let mut energy = 0.0;
    for q in 0..4 {
        let v: [f64; 3] = core::array::from_fn(|axis| {
            (0..4)
                .map(|node| (if node == q { a } else { b }) * velocity[node][axis])
                .sum()
        });
        energy += (DENSITY * VOLUME / 8.0) * v.iter().map(|x| x * x).sum::<f64>();
    }
    energy
}

fn check_balances(energy: &ThermalBrittleTetDynamicEnergy) {
    for defect in [
        energy.strain_split_error_j,
        energy.material_balance_error_j,
        energy.equilibrium_work_residual_j,
        energy.inertia.balance_error_j,
        energy.total_balance_error_j,
    ] {
        assert!(defect.abs() < 1e-10, "{energy:?}");
    }
}

#[test]
fn cooling_causes_inertial_tensile_softening_with_independent_motion_and_work() {
    with_cx(|cx, _| {
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, budget(), 4, cx).unwrap();
        let card = card("cooling", 4.5e-5, false);
        let reference = [source(&card, &geometry, 400.0, 400.0)];
        let cold = [reference[0].at_temperature(&card, 300.0).unwrap()];
        let solid = ThermalBrittleTetSolid::try_new(&geometry, &reference, budget(), cx).unwrap();
        let old =
            ThermalBrittleTetMotion::try_from_solid(&solid, &[[0.0; 3]; 4], 0.0, 8, cx).unwrap();
        let dt = 0.01_f64;
        let supports = axial_supports();
        let step = old
            .advance(&cold, dt, &supports, &[], settings(), controls(), cx)
            .unwrap();

        // One free x DOF: m*(u/dt^2) + V*sigma(u-free) = 0.
        // Its consistent diagonal mass is rho*V/10, not the lumped rho*V/4.
        // Inertia makes this softening branch locally nonsingular.
        let free = -4.5e-5 * 100.0;
        let onset = STRENGTH / YOUNG;
        let failure = 2.0 * GC / (STRENGTH * (1.0 / 3.0));
        let slope = STRENGTH / (failure - onset);
        let mass = DENSITY * VOLUME / 10.0;
        let u = -VOLUME * slope * (failure + free) / (mass / dt.powi(2) - VOLUME * slope);
        let strain = u - free;
        let stress = slope * (failure - strain);
        let velocity = u / dt;
        assert!(strain > onset && strain < failure);
        near(u, -0.0000625, 1e-16);
        near(stress, 0.375, 1e-13);
        near(
            step.motion.solid().solution().displacement_m[1][0],
            u,
            2e-12,
        );
        near(step.motion.velocity_mps()[1][0], velocity, 2e-10);
        near(
            step.motion.solid().solution().element_strain[0][0],
            u,
            2e-12,
        );
        near(
            step.motion.solid().solution().element_states[0].strain()[0],
            strain,
            2e-12,
        );
        near(
            step.motion.solid().solution().element_stress[0][0],
            stress,
            2e-10,
        );
        for axis in 1..3 {
            near(
                step.motion.solid().solution().element_strain[0][axis],
                0.0,
                1e-15,
            );
            near(
                step.motion.solid().solution().element_states[0].strain()[axis],
                -free,
                1e-15,
            );
            near(
                step.motion.solid().solution().element_stress[0][axis],
                -YOUNG * free,
                2e-12,
            );
        }

        let opening_work = 0.5 * STRENGTH * onset + 0.5 * (STRENGTH + stress) * (strain - onset);
        let opening_storage = 0.5 * stress * strain;
        let transverse_storage = YOUNG * free * free;
        let stored = VOLUME * (opening_storage + transverse_storage);
        let fracture = VOLUME * (opening_work - opening_storage);
        let mean_x = opening_work / strain;
        let mechanical_work = VOLUME * mean_x * u;
        let eigenstrain_work = -VOLUME * free * (mean_x - YOUNG * free);
        let endpoint_work = VOLUME * stress * u;
        let kinetic = 0.5 * mass * velocity * velocity;
        near(
            step.motion.solid().solution().element_states[0].damage(),
            1.0 - stress / (YOUNG * strain),
            2e-11,
        );
        near(step.motion.solid().stored_energy_j(), stored, 2e-12);
        near(
            step.motion.solid().fracture_dissipation_j(),
            fracture,
            2e-12,
        );
        near(step.energy.fracture_dissipation_j, fracture, 2e-12);
        near(
            step.energy.constitutive_strain_work_j,
            VOLUME * (opening_work + transverse_storage),
            2e-12,
        );
        near(step.energy.mechanical_work_j, mechanical_work, 2e-12);
        near(step.energy.eigenstrain_work_j, eigenstrain_work, 2e-12);
        near(step.energy.internal_endpoint_work_j, endpoint_work, 2e-12);
        near(
            step.energy.constitutive_endpoint_defect_j,
            endpoint_work - mechanical_work,
            2e-12,
        );
        near(step.motion.kinetic_energy_j(), kinetic, 2e-12);
        near(step.energy.inertia.numerical_loss_j, kinetic, 2e-12);
        near(step.energy.inertia.inertial_work_j, 2.0 * kinetic, 2e-12);
        near(step.energy.parameter_work_j, 0.0, 2e-14);
        near(step.energy.parameter_fracture_dissipation_j, 0.0, 2e-14);
        assert_eq!(step.energy.dead_load_endpoint_work_j, 0.0);
        assert_eq!(step.energy.support_endpoint_work_j, 0.0);
        near(
            eigenstrain_work,
            stored + fracture + 2.0 * kinetic + step.energy.constitutive_endpoint_defect_j,
            2e-12,
        );
        check_balances(&step.energy);
        assert_eq!(step.motion.time_s(), dt);
        assert_eq!(step.motion.accepted_steps(), 1);
        near(step.motion.reference_mass_kg(), DENSITY * VOLUME, 2e-15);
        assert!(std::ptr::eq(
            step.motion.solid().materials().as_ptr(),
            cold.as_ptr()
        ));
        assert_eq!(old.solid().solution().element_states[0].damage(), 0.0);
        assert_eq!(old.solid().fracture_dissipation_j(), 0.0);
        assert_eq!(old.velocity_mps(), &[[0.0; 3]; 4]);
        assert_eq!(old.time_s(), 0.0);
    });
}

#[test]
fn temperature_dependent_queried_density_cannot_change_reference_mass_or_momentum() {
    with_cx(|cx, _| {
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, budget(), 4, cx).unwrap();
        let card = card("inventory", 0.0, true);
        let reference = [source(&card, &geometry, 300.0, 300.0)];
        let hot = [reference[0].at_temperature(&card, 400.0).unwrap()];
        near(
            hot[0].band().reference_density_kg_m3().unwrap(),
            0.5 * DENSITY,
            1e-14,
        );
        let solid = ThermalBrittleTetSolid::try_new(&geometry, &reference, budget(), cx).unwrap();
        let velocity = [[0.03, -0.02, 0.01]; 4];
        let old = ThermalBrittleTetMotion::try_from_solid(&solid, &velocity, 1.5, 4, cx).unwrap();
        let dt = 0.05;
        let force = 2.0;
        let step = old
            .advance(
                &hot,
                dt,
                &[],
                &uniform_x_force(force),
                settings(),
                controls(),
                cx,
            )
            .unwrap();
        let mass = DENSITY * VOLUME;
        let expected_v = [
            velocity[0][0] + force * dt / mass,
            velocity[0][1],
            velocity[0][2],
        ];
        for (u, v) in step
            .motion
            .solid()
            .solution()
            .displacement_m
            .iter()
            .zip(step.motion.velocity_mps())
        {
            for axis in 0..3 {
                near(v[axis], expected_v[axis], 2e-11);
                near(u[axis], dt * expected_v[axis], 2e-12);
            }
        }
        let kinetic_change = 0.5 * mass * (expected_v[0].powi(2) - velocity[0][0].powi(2));
        let numerical = 0.5 * mass * (expected_v[0] - velocity[0][0]).powi(2);
        near(step.motion.reference_mass_kg(), mass, 2e-15);
        assert_eq!(step.motion.reference_density_kg_m3(), &[DENSITY]);
        near(step.energy.inertia.kinetic_change_j, kinetic_change, 2e-12);
        near(step.energy.inertia.numerical_loss_j, numerical, 2e-12);
        near(
            step.energy.dead_load_endpoint_work_j,
            force * dt * expected_v[0],
            2e-12,
        );
        near(
            step.energy.dead_load_endpoint_work_j,
            kinetic_change + numerical,
            2e-12,
        );
        near(
            step.motion.kinetic_energy_j(),
            quadrature_kinetic(step.motion.velocity_mps()),
            2e-12,
        );
        near(step.energy.mechanical_work_j, 0.0, 2e-12);
        assert_eq!(step.energy.eigenstrain_work_j, 0.0);
        near(step.energy.parameter_work_j, 0.0, 2e-12);
        near(step.motion.solid().stored_energy_j(), 0.0, 2e-12);
        assert_eq!(step.motion.solid().fracture_dissipation_j(), 0.0);
        assert_eq!(step.energy.support_endpoint_work_j, 0.0);
        check_balances(&step.energy);

        // A second source refresh supplies neither an impulse nor a mass update.
        let coast = step
            .motion
            .advance(&reference, dt, &[], &[], settings(), controls(), cx)
            .unwrap();
        for v in coast.motion.velocity_mps() {
            for axis in 0..3 {
                near(v[axis], expected_v[axis], 2e-11);
            }
        }
        near(coast.motion.reference_mass_kg(), mass, 2e-15);
        assert_eq!(coast.motion.reference_density_kg_m3(), &[DENSITY]);
        near(coast.energy.inertia.kinetic_change_j, 0.0, 2e-12);
        near(coast.energy.inertia.numerical_loss_j, 0.0, 2e-12);
        check_balances(&coast.energy);
        assert_eq!(old.velocity_mps(), &velocity);
        assert_eq!(old.time_s(), 1.5);
        assert_eq!(old.accepted_steps(), 0);
    });
}

#[test]
fn matched_free_expansion_keeps_native_damage_zero_and_support_inertia_explicit() {
    with_cx(|cx, _| {
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, budget(), 4, cx).unwrap();
        let card = card("free-expansion", 1e-5, false);
        let reference = [source(&card, &geometry, 300.0, 300.0)];
        let hot = [reference[0].at_temperature(&card, 400.0).unwrap()];
        let solid = ThermalBrittleTetSolid::try_new(&geometry, &reference, budget(), cx).unwrap();
        let dt = 0.1;
        let free = 1e-5 * 100.0;
        let targets = affine_targets(free);
        let compatible: [[f64; 3]; 4] =
            core::array::from_fn(|node| core::array::from_fn(|axis| free * NODES[node][axis] / dt));
        let coasting =
            ThermalBrittleTetMotion::try_from_solid(&solid, &compatible, 0.0, 2, cx).unwrap();
        let accelerated =
            ThermalBrittleTetMotion::try_from_solid(&solid, &[[0.0; 3]; 4], 0.0, 2, cx).unwrap();
        let coast = coasting
            .advance(&hot, dt, &targets, &[], settings(), controls(), cx)
            .unwrap();
        let driven = accelerated
            .advance(&hot, dt, &targets, &[], settings(), controls(), cx)
            .unwrap();
        let kinetic = quadrature_kinetic(&compatible);
        assert!(kinetic > 0.0);
        for step in [&coast, &driven] {
            for (node, (u, v)) in step
                .motion
                .solid()
                .solution()
                .displacement_m
                .iter()
                .zip(step.motion.velocity_mps())
                .enumerate()
            {
                for axis in 0..3 {
                    near(u[axis], free * NODES[node][axis], 2e-15);
                    near(v[axis], compatible[node][axis], 2e-14);
                }
            }
            for component in 0..6 {
                near(
                    step.motion.solid().solution().element_strain[0][component],
                    if component < 3 { free } else { 0.0 },
                    2e-15,
                );
                near(
                    step.motion.solid().solution().element_states[0].strain()[component],
                    0.0,
                    2e-15,
                );
                near(
                    step.motion.solid().solution().element_stress[0][component],
                    0.0,
                    2e-12,
                );
            }
            assert_eq!(
                step.motion.solid().solution().element_states[0].damage(),
                0.0
            );
            assert_eq!(step.motion.solid().fracture_dissipation_j(), 0.0);
            near(step.energy.mechanical_work_j, 0.0, 2e-14);
            near(step.energy.eigenstrain_work_j, 0.0, 2e-14);
            near(step.energy.constitutive_strain_work_j, 0.0, 2e-14);
            near(step.energy.parameter_work_j, 0.0, 2e-14);
            near(step.motion.kinetic_energy_j(), kinetic, 2e-14);
            check_balances(&step.energy);
        }
        near(coast.energy.inertia.kinetic_change_j, 0.0, 2e-14);
        near(coast.energy.inertia.numerical_loss_j, 0.0, 2e-14);
        near(coast.energy.support_endpoint_work_j, 0.0, 2e-14);
        near(driven.energy.inertia.kinetic_change_j, kinetic, 2e-14);
        near(driven.energy.inertia.numerical_loss_j, kinetic, 2e-14);
        near(driven.energy.inertia.inertial_work_j, 2.0 * kinetic, 2e-14);
        near(driven.energy.support_endpoint_work_j, 2.0 * kinetic, 2e-14);
        assert_eq!(driven.energy.dead_load_endpoint_work_j, 0.0);
    });
}

#[test]
fn source_reference_work_solver_and_budget_refusals_preserve_damaged_motion() {
    with_cx(|cx, gate| {
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, budget(), 4, cx).unwrap();
        let card = card("rollback", 4.5e-5, false);
        let reference = [source(&card, &geometry, 400.0, 400.0)];
        let cold = [reference[0].at_temperature(&card, 300.0).unwrap()];
        let solid = ThermalBrittleTetSolid::try_new(&geometry, &reference, budget(), cx).unwrap();
        let initial =
            ThermalBrittleTetMotion::try_from_solid(&solid, &[[0.0; 3]; 4], 1.0, 2, cx).unwrap();
        let supports = axial_supports();
        let old = initial
            .advance(&cold, 0.01, &supports, &[], settings(), controls(), cx)
            .unwrap()
            .motion;
        assert!(old.solid().solution().element_states[0].damage() > 0.0);
        let snapshot = old.solid().solution().clone();
        let velocity = old.velocity_mps().to_vec();
        let kinetic = old.kinetic_energy_j();
        let loss = old.solid().fracture_dissipation_j();
        let reset_reference = [source(&card, &geometry, 350.0, 300.0)];
        let foreign_card = self::card("different-original-source", 4.5e-5, false);
        let foreign = [source(&foreign_card, &geometry, 400.0, 300.0)];
        for next in [&reset_reference, &foreign] {
            assert!(matches!(
                old.advance(next, 0.01, &supports, &[], settings(), controls(), cx),
                Err(BrittleTetDynamicError::Source { .. })
            ));
        }
        for invalid in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(matches!(
                old.advance(&cold, invalid, &supports, &[], settings(), controls(), cx),
                Err(BrittleTetDynamicError::Invalid(_))
            ));
        }
        let strict = BrittleTetDynamicControls {
            maximum_numerical_kinetic_loss_j: 0.0,
            ..controls()
        };
        assert!(matches!(
            old.advance(&cold, 0.01, &supports, &[], settings(), strict, cx),
            Err(BrittleTetDynamicError::Solid(BrittleTetError::EnergyGate { defect_j, .. }))
                if defect_j > 0.0
        ));
        let stalled = IncrementSettings {
            max_newton_iterations: 1,
            ..settings()
        };
        assert!(matches!(
            old.advance(&cold, 0.01, &supports, &[], stalled, controls(), cx),
            Err(BrittleTetDynamicError::Solid(BrittleTetError::Solve(
                SmallStrainTetError::NewtonStalled { .. }
            )))
        ));
        gate.request();
        assert!(matches!(
            old.advance(&cold, 0.01, &supports, &[], settings(), controls(), cx),
            Err(BrittleTetDynamicError::Solid(BrittleTetError::Cancelled))
        ));
        assert_eq!(old.solid().solution(), &snapshot);
        assert_eq!(old.velocity_mps(), velocity.as_slice());
        assert_eq!(old.kinetic_energy_j(), kinetic);
        assert_eq!(old.solid().fracture_dissipation_j(), loss);
        assert_eq!(old.time_s(), 1.01);
        assert_eq!(old.accepted_steps(), 1);
        assert_eq!(old.remaining_steps(), 1);
        assert_eq!(old.reference_density_kg_m3(), &[DENSITY]);

        let retry_gate = CancelGate::new_clock_free();
        let retry = Cx::new(
            &retry_gate,
            cx.arena(),
            cx.stream_key(),
            cx.budget(),
            cx.mode(),
        );
        let accepted = old
            .advance(&cold, 0.01, &supports, &[], settings(), controls(), &retry)
            .unwrap();
        let repeated = old
            .advance(&cold, 0.01, &supports, &[], settings(), controls(), &retry)
            .unwrap();
        assert_eq!(
            accepted.motion.solid().solution(),
            repeated.motion.solid().solution()
        );
        assert_eq!(
            accepted.motion.velocity_mps(),
            repeated.motion.velocity_mps()
        );
        assert_eq!(accepted.energy, repeated.energy);
        near(
            accepted.motion.solid().fracture_dissipation_j(),
            loss,
            2e-12,
        );
        assert_eq!(accepted.motion.accepted_steps(), 2);
        assert_eq!(accepted.motion.remaining_steps(), 0);
        check_balances(&accepted.energy);
        assert!(matches!(
            accepted
                .motion
                .advance(&cold, 0.01, &supports, &[], settings(), controls(), &retry),
            Err(BrittleTetDynamicError::StepBudget {
                accepted_steps: 2,
                maximum_steps: 2
            })
        ));
        assert_eq!(old.solid().solution(), &snapshot);
        assert_eq!(old.velocity_mps(), velocity.as_slice());
        assert_eq!(old.accepted_steps(), 1);
    });
}

#[path = "brittle_tet_thermal_motion/midpoint.rs"]
mod midpoint;

#[path = "brittle_tet_thermal_motion/modes.rs"]
mod modes;
