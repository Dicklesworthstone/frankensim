//! Manufactured Maxwell moments and finite-state conservation, not measured
//! surface coefficients. The existing example loader supplies real NASA-9 and
//! transport owners, including a temperature-dependent heat capacity.
#[path = "../../fs-material/examples/material_gas_spring/source.rs"]
mod source;

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;
use fs_blake3::hash_bytes;
use fs_couple::liquid_gas::evaporation::momentum::MovingGasParcel;
use fs_couple::molecular_shear::{
    self, MolecularShearControls, MolecularShearError, MolecularShearSurface, MolecularShearWall,
};
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};
use fs_material::sublimation::gas::{GasSurfaceBinding, momentum::MaxwellReflection};
use fs_material::thermochemical_gas::energy::{GasEnergySolve, parcel::ClosedGasParcel};
use fs_mbd::{MassProperties, Pose, RigidBodyState, Vec3};
use fs_thermochem::{SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};

fn chart() -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_new(
        hash_bytes(b"synthetic constant-volume latent wall"),
        [
            (0.0, 200.0, 0.0),
            (50_000.0, 300.0, 0.0),
            (150_000.0, 300.0, 1.0),
            (250_000.0, 500.0, 1.0),
        ]
        .map(|(h, t, liquid)| EnthalpyPhaseKnot {
            specific_enthalpy_j_kg: h,
            temperature_k: t,
            liquid_mass_fraction: liquid,
            bulk_density_kg_m3: 1000.0,
        })
        .to_vec(),
    )
    .unwrap()
}
fn gas_at(velocity: [f64; 3], temperature: f64) -> MovingGasParcel {
    let model = source::load(
        "A 0.75 0.028 3.5 0.001 2e-5 0.03\nB 0.25 0.004 2.5 0.0003 1.4e-5 0.021\n",
        "manufactured two-species molecular shear",
    )
    .unwrap();
    MovingGasParcel::new(
        ClosedGasParcel::new(model, source::point(temperature, 100_000.0).unwrap(), 0.01).unwrap(),
        velocity,
    )
    .unwrap()
}
fn gas(velocity: [f64; 3]) -> MovingGasParcel {
    gas_at(velocity, 300.0)
}
fn scattering(gas: &MovingGasParcel, chi: [f64; 2]) -> MaxwellReflection {
    MaxwellReflection::from_diffuse_fractions(
        gas.parcel().model(),
        &[
            (SpeciesId::new("A").unwrap(), chi[0]),
            (SpeciesId::new("B").unwrap(), chi[1]),
        ],
        [200.0, 1000.0],
    )
    .unwrap()
}
fn wall(
    chart: &EquilibriumEnthalpyPhaseCurve,
    h: f64,
    velocity: [f64; 3],
) -> MolecularShearWall<'_> {
    let mass = 0.02;
    let properties = MassProperties::new(mass, Vec3::ZERO, Vec3::new(0.001, 0.001, 0.001)).unwrap();
    let body = RigidBodyState::new(
        Pose::identity(),
        Vec3::new(mass * velocity[0], mass * velocity[1], mass * velocity[2]),
        Vec3::ZERO,
    )
    .unwrap();
    MolecularShearWall::new(
        body,
        properties,
        mass,
        chart,
        chart.state_at_specific_enthalpy(h).unwrap(),
    )
    .unwrap()
}
fn surface(area: f64) -> MolecularShearSurface {
    MolecularShearSurface {
        area_m2: area,
        outward_normal: [0.0, 0.0, 1.0],
        centroid_arm_body_m: [0.03, -0.01, 0.02],
    }
}
fn controls() -> MolecularShearControls {
    MolecularShearControls {
        gas_energy: GasEnergySolve {
            minimum_temperature_k: 200.0,
            maximum_temperature_k: 1000.0,
            energy_tolerance_j_kg: 1e-8,
            maximum_iterations: 128,
        },
        maximum_duration_s: 1.0,
        maximum_temperature_change_k: 100.0,
        maximum_species: 2,
        maximum_phase_knots: 8,
        momentum_absolute_n_s: 1e-13,
        energy_absolute_j: 1e-8,
        position_absolute_m: 1e-12,
        relative_tolerance: 1e-10,
    }
}
fn near(actual: f64, expected: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= absolute,
        "{actual:e} != {expected:e}, budget {absolute:e}"
    );
}
fn velocity(wall: &MolecularShearWall<'_>) -> [f64; 3] {
    let v = wall
        .body()
        .center_of_mass_velocity_world(wall.properties())
        .unwrap();
    [v.x, v.y, v.z]
}

#[test]
fn finite_momenta_follow_exact_relaxation_and_original_latent_and_nasa_stores() {
    let chart = chart();
    let gas = gas([40.0, -10.0, 0.0]);
    let wall = wall(&chart, 100_000.0, [-2.0, 3.0, 0.0]);
    let chi = [0.8, 0.25];
    let reflection = scattering(&gas, chi);
    let area = 0.001;
    let dt = 0.001;
    let accepted = molecular_shear::step_candidate(
        &gas,
        wall,
        surface(area),
        &reflection,
        dt,
        controls(),
        || false,
    )
    .unwrap();
    let step = accepted.step();
    let mg = gas.parcel().mass_kg();
    let mw = wall.thermal_mass_kg();
    let mu = mg * mw / (mg + mw);
    // Independent half-Maxwellian Gaussian mass integral for each species.
    let encounters: f64 = [(0.75, 0.028, chi[0]), (0.25, 0.004, chi[1])]
        .into_iter()
        .map(|(x, molar_mass, fraction)| {
            fraction * 100_000.0 * x * area * dt
                / (2.0 * core::f64::consts::PI * R / molar_mass * 300.0).sqrt()
        })
        .sum();
    let gamma = encounters / mu;
    let relaxed = -(-gamma).exp_m1();
    let before = [42.0, -13.0, 0.0];
    let loss = 0.5 * mu * -(-2.0 * gamma).exp_m1() * (42.0 * 42.0 + 13.0 * 13.0);
    near(step.diffuse_encounter_mass_kg, encounters, 1e-16);
    near(step.relaxed_slip_fraction, relaxed, 1e-14);
    near(step.dissipated_relative_kinetic_energy_j, loss, 1e-12);
    for d in 0..3 {
        let impulse = mu * relaxed * before[d];
        near(step.tangential_impulse_on_wall_n_s[d], impulse, 1e-14);
        near(
            accepted.gas().velocity_m_s()[d],
            gas.velocity_m_s()[d] - impulse / mg,
            1e-12,
        );
        near(
            velocity(accepted.wall())[d],
            velocity(&wall)[d] + impulse / mw,
            1e-12,
        );
        near(step.slip_after_m_s[d], before[d] * (-gamma).exp(), 1e-12);
        let distance = dt * (velocity(&wall)[d] + mu / mw * before[d] * (1.0 - relaxed / gamma));
        near(step.wall_displacement_m[d], distance, 1e-14);
        near(step.momentum_residual_n_s[d], 0.0, 1e-13);
    }
    assert_eq!(step.reflection.heat_into_gas_j, 0.0);
    near(step.gas_heat.stored_energy_change_j, 0.5 * loss, 1e-8);
    near(step.stored_wall_heat_j, 0.5 * loss, 1e-8);
    near(step.mechanical_energy_residual_j, 0.0, 1e-8);
    near(step.thermal_energy_residual_j, 0.0, 1e-8);
    near(step.total_energy_residual_j, 0.0, 1e-8);
    assert_eq!(accepted.wall().phase().temperature_k(), 300.0);
    assert!(accepted.wall().phase().liquid_mass_fraction() > wall.phase().liquid_mass_fraction());
    assert!(accepted.gas().parcel().state().state().temperature > 300.0);
    assert_eq!(accepted.gas().parcel().mass_kg(), mg);
    assert_eq!(
        accepted.gas().parcel().volume_m3(),
        gas.parcel().volume_m3()
    );
    assert_eq!(accepted.gas().parcel().model(), gas.parcel().model());
    // Continuation retains the actual warmed finite gas and latent wall state.
    let next = molecular_shear::step_candidate(
        accepted.gas(),
        *accepted.wall(),
        surface(area),
        &reflection,
        dt,
        controls(),
        || false,
    )
    .unwrap();
    assert!(next.step().slip_after_m_s[0] < step.slip_after_m_s[0]);
    assert!(
        next.wall().phase().specific_enthalpy_j_kg()
            > accepted.wall().phase().specific_enthalpy_j_kg()
    );
}

#[test]
fn complete_pressure_and_both_normal_supports_and_wall_guide_moment_are_retained() {
    let chart = chart();
    let gas = gas([20.0, 5.0, 0.0]);
    let wall = wall(&chart, 200_000.0, [-3.0, 2.0, 0.0]); // T_wall = 400 K.
    let chi = [1.0, 0.4];
    let reflection = scattering(&gas, chi);
    let surface = surface(0.0001);
    let dt = 0.0001;
    let accepted =
        molecular_shear::step_candidate(&gas, wall, surface, &reflection, dt, controls(), || false)
            .unwrap();
    let report = accepted.step();
    let expected_normal = 100_000.0
        * surface.area_m2
        * dt
        * (1.0 + 0.5 * (0.75 * chi[0] + 0.25 * chi[1]) * ((400.0_f64 / 300.0).sqrt() - 1.0));
    near(report.normal_molecular_impulse_n_s, expected_normal, 1e-14);
    near(
        report.molecular_impulse_on_wall_n_s[2],
        -expected_normal,
        1e-14,
    );
    near(
        report.normal_support_impulse_on_wall_n_s[2],
        expected_normal,
        1e-14,
    );
    for d in 0..3 {
        assert_eq!(
            report.molecular_impulse_on_gas_n_s[d],
            -report.molecular_impulse_on_wall_n_s[d]
        );
        assert_eq!(
            report.normal_support_impulse_on_gas_n_s[d],
            -report.normal_support_impulse_on_wall_n_s[d]
        );
        assert_eq!(
            report.guide_angular_impulse_on_wall_n_m_s[d],
            -report.molecular_angular_impulse_on_wall_n_m_s[d]
        );
    }
    let [x, y, z] = surface.centroid_arm_body_m;
    let [jx, jy, jz] = report.molecular_impulse_on_wall_n_s;
    for (actual, expected) in report
        .molecular_angular_impulse_on_wall_n_m_s
        .into_iter()
        .zip([y * jz - z * jy, z * jx - x * jz, x * jy - y * jx])
    {
        near(actual, expected, 1e-16);
    }
    assert_eq!(accepted.wall().body().angular_momentum_body(), Vec3::ZERO);
    assert_eq!(report.support_work_j, 0.0);
    assert_eq!(accepted.gas().velocity_m_s()[2], 0.0);
    assert_eq!(velocity(accepted.wall())[2], 0.0);
    assert!(report.reflection.heat_into_gas_j > 0.0);
    near(
        report.gas_heat.supplied_heat_j,
        report.reflection.heat_into_gas_j + 0.5 * report.dissipated_relative_kinetic_energy_j,
        1e-15,
    );
    near(
        report.requested_wall_heat_j,
        -report.reflection.heat_into_gas_j + 0.5 * report.dissipated_relative_kinetic_energy_j,
        1e-15,
    );
    // The new sealed seam is exactly the original zero-sticking source path.
    let binding = GasSurfaceBinding {
        species: SpeciesId::new("A").unwrap(),
        accommodation_coefficient: 0.0,
        solid_enthalpy_offset_j_kg: 0.0,
    };
    let original = reflection
        .thermal_accommodation()
        .evaluate(&binding, gas.parcel(), 400.0, surface.area_m2, dt, || false)
        .unwrap();
    assert_eq!(report.reflection, original);
}

#[test]
fn specular_reflection_retains_pressure_with_no_drag_or_heat() {
    let chart = chart();
    let gas = gas([40.0, -10.0, 0.0]);
    let wall = wall(&chart, 200_000.0, [-2.0, 3.0, 0.0]);
    let reflection = scattering(&gas, [0.0, 0.0]);
    let result = molecular_shear::step_candidate(
        &gas,
        wall,
        surface(0.001),
        &reflection,
        0.001,
        controls(),
        || false,
    )
    .unwrap();
    assert_eq!(result.gas(), &gas);
    assert_eq!(result.wall().phase(), wall.phase());
    assert_eq!(
        result.wall().body().linear_momentum_world(),
        wall.body().linear_momentum_world()
    );
    assert_eq!(result.step().diffuse_encounter_mass_kg, 0.0);
    assert_eq!(result.step().dissipated_relative_kinetic_energy_j, 0.0);
    assert_eq!(result.step().gas_heat.supplied_heat_j, 0.0);
    near(result.step().normal_molecular_impulse_n_s, 0.1, 1e-15);
    assert_eq!(result.step().wall_displacement_m, [-0.002, 0.003, 0.0]);
}

#[test]
fn stiff_finite_inertia_relaxation_approaches_common_velocity_without_slip_reversal() {
    let chart = chart();
    let gas = gas([40.0, 0.0, 0.0]);
    let wall = wall(&chart, 100_000.0, [-2.0, 0.0, 0.0]);
    let reflection = scattering(&gas, [1.0, 1.0]);
    let result = molecular_shear::step_candidate(
        &gas,
        wall,
        surface(100.0),
        &reflection,
        0.001,
        controls(),
        || false,
    )
    .unwrap();
    let mg = gas.parcel().mass_kg();
    let mw = wall.thermal_mass_kg();
    let common = (mg * 40.0 - mw * 2.0) / (mg + mw);
    assert_eq!(result.step().relaxed_slip_fraction, 1.0);
    near(result.gas().velocity_m_s()[0], common, 1e-12);
    near(velocity(result.wall())[0], common, 1e-12);
    near(result.step().slip_after_m_s[0], 0.0, 1e-12);
    let mu = mg * mw / (mg + mw);
    near(
        result.step().dissipated_relative_kinetic_energy_j,
        0.5 * mu * 42.0 * 42.0,
        1e-12,
    );
    assert!(result.step().normal_molecular_impulse_n_s > 1000.0);
    near(result.step().mechanical_energy_residual_j, 0.0, 1e-8);
}

#[test]
fn common_tangent_boost_does_not_change_physical_heat_or_relaxation() {
    let chart = chart();
    let run = |boost: f64| {
        let gas = gas([40.0 + boost, -10.0 + boost, 0.0]);
        let wall = wall(&chart, 100_000.0, [-2.0 + boost, 3.0 + boost, 0.0]);
        molecular_shear::step_candidate(
            &gas,
            wall,
            surface(0.001),
            &scattering(&gas, [0.8, 0.25]),
            0.001,
            controls(),
            || false,
        )
        .unwrap()
    };
    let original = run(0.0);
    let boosted = run(64.0);
    assert_eq!(original.step().reflection, boosted.step().reflection);
    near(
        original.step().dissipated_relative_kinetic_energy_j,
        boosted.step().dissipated_relative_kinetic_energy_j,
        1e-12,
    );
    near(
        original.gas().parcel().state().state().temperature,
        boosted.gas().parcel().state().state().temperature,
        1e-10,
    );
    assert_eq!(original.wall().phase(), boosted.wall().phase());
    near(boosted.step().total_energy_residual_j, 0.0, 1e-8);

    // A huge common momentum cannot hide an unrepresented small exchanged
    // impulse behind either the boost or the much larger normal pressure load.
    let gas = gas([1e16 + 4.0, 0.0, 0.0]);
    let wall = wall(&chart, 100_000.0, [1e16, 0.0, 0.0]);
    let mut tight = controls();
    tight.momentum_absolute_n_s = 1e-15;
    assert!(matches!(
        molecular_shear::step_candidate(
            &gas,
            wall,
            surface(1e-10),
            &scattering(&gas, [1.0, 1.0]),
            0.001,
            tight,
            || false
        ),
        Err(MolecularShearError::Balance {
            account: "gas impulse" | "wall impulse",
            ..
        })
    ));
}

#[test]
fn cancellation_source_work_and_physical_constraint_refusals_are_atomic() {
    let chart = chart();
    let gas = gas([40.0, 0.0, 0.0]);
    let original = gas.clone();
    let wall = wall(&chart, 100_000.0, [0.0; 3]);
    let reflection = scattering(&gas, [0.8, 0.25]);
    let mut polls = 0;
    molecular_shear::step_candidate(
        &gas,
        wall,
        surface(0.001),
        &reflection,
        0.001,
        controls(),
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    assert!(polls > 5);
    let mut late_polls = 0;
    assert!(matches!(
        molecular_shear::step_candidate(
            &gas,
            wall,
            surface(0.001),
            &reflection,
            0.001,
            controls(),
            || {
                late_polls += 1;
                late_polls >= polls
            }
        ),
        Err(MolecularShearError::Cancelled)
    ));
    for altered in [
        MolecularShearControls {
            maximum_species: 1,
            ..controls()
        },
        MolecularShearControls {
            maximum_phase_knots: 3,
            ..controls()
        },
        MolecularShearControls {
            maximum_duration_s: 1e-4,
            ..controls()
        },
        MolecularShearControls {
            maximum_temperature_change_k: 1e-10,
            ..controls()
        },
    ] {
        assert!(
            molecular_shear::step_candidate(
                &gas,
                wall,
                surface(0.001),
                &reflection,
                0.001,
                altered,
                || false
            )
            .is_err()
        );
    }
    let normal_gas = MovingGasParcel::new(gas.parcel().clone(), [40.0, 0.0, 1e-20]).unwrap();
    assert!(
        molecular_shear::step_candidate(
            &normal_gas,
            wall,
            surface(0.001),
            &reflection,
            0.001,
            controls(),
            || false
        )
        .is_err()
    );
    assert!(
        MolecularShearWall::new(wall.body(), wall.properties(), 0.03, &chart, wall.phase())
            .is_err()
    );
    let spin = RigidBodyState::new(Pose::identity(), Vec3::ZERO, Vec3::new(0.0, 0.0, 1.0)).unwrap();
    assert!(MolecularShearWall::new(spin, wall.properties(), 0.02, &chart, wall.phase()).is_err());
    assert_eq!(gas, original);
    assert_eq!(wall.phase().specific_enthalpy_j_kg(), 100_000.0);
    assert_eq!(wall.body().linear_momentum_world(), Vec3::ZERO);
}

#[test]
fn equal_endpoint_density_does_not_hide_expansion_inside_a_crossed_latent_plateau() {
    let chart = EquilibriumEnthalpyPhaseCurve::try_new(
        hash_bytes(b"synthetic hidden expansion"),
        [
            (0.0, 200.0, 0.0, 1000.0),
            (50_000.0, 300.0, 0.0, 1000.0),
            (100_000.0, 300.0, 0.5, 900.0),
            (150_000.0, 300.0, 1.0, 1000.0),
            (250_000.0, 500.0, 1.0, 1000.0),
        ]
        .map(|(h, t, f, rho)| EnthalpyPhaseKnot {
            specific_enthalpy_j_kg: h,
            temperature_k: t,
            liquid_mass_fraction: f,
            bulk_density_kg_m3: rho,
        })
        .to_vec(),
    )
    .unwrap();
    let initial = gas_at([0.0; 3], 250.0);
    let mg = initial.parcel().mass_kg();
    let mw = 0.02;
    let mu = mg * mw / (mg + mw);
    // D/2 = 3000 J would move h from 25k to 175k, both at rho=1000.
    let gas =
        MovingGasParcel::new(initial.parcel().clone(), [(12_000.0 / mu).sqrt(), 0.0, 0.0]).unwrap();
    let wall = wall(&chart, 25_000.0, [0.0; 3]);
    assert_eq!(
        chart
            .state_at_specific_enthalpy(175_000.0)
            .unwrap()
            .bulk_density_kg_m3(),
        wall.phase().bulk_density_kg_m3()
    );
    assert!(matches!(
        molecular_shear::step_candidate(
            &gas,
            wall,
            surface(100.0),
            &scattering(&gas, [1.0, 1.0]),
            0.001,
            controls(),
            || false
        ),
        Err(MolecularShearError::Invalid(
            "crossed wall source segment requires unsupported expansion work"
        ))
    ));
}
