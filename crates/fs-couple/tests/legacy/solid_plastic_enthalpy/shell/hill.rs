//! Directional yielding with original caloric feedback and retained local history.
//! These coefficients are manufactured; no rolled-sheet calibration is claimed.

use super::*;
use fs_material::plastic::hill::source::{
    HILL_COEFFICIENT_PROPERTIES, HILL_HARDENING_MODULUS_PROPERTY,
    HILL_REFERENCE_YIELD_STRESS_PROPERTY,
};
use fs_material::plastic::hill::thermal::{ThermalHillSource, ThermalHillSourceQuery};
use fs_qty::QuantitySpec;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};

const CP_HILL: f64 = 100.0;
const STRESS: f64 = 2.6e6;
const DT: f64 = 0.001;
const EXTERNAL_RISE_K: f64 = 4.0;

fn temperature_spec() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}

fn point(temperature_k: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", temperature_spec(), temperature_k)
        .unwrap()
}

fn fixture() -> Fixture {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        (YOUNG_MODULUS_PROPERTY, Pressure::DIMS, E),
        (POISSON_RATIO_PROPERTY, Dims::NONE, 0.2),
        (HILL_HARDENING_MODULUS_PROPERTY, Pressure::DIMS, H),
        (DENSITY_PROPERTY, Density::DIMS, RHO),
        (CP, CP_DIMS, CP_HILL),
    ]
    .into_iter()
    .chain(
        HILL_COEFFICIENT_PROPERTIES
            .into_iter()
            .zip([0.2, 1.2, 0.6, 1.8, 1.6, 1.1])
            .map(|(name, value)| (name, Dims::NONE, value)),
    ) {
        let mut claim = scalar(name, dims, value);
        claim.provenance.source = "manufactured anisotropic thermal shell".into();
        claim.validity =
            ValidityDomain::unconstrained().with_quantity("T", temperature_spec(), 280.0, 320.0);
        claims.insert_claim(claim).unwrap();
    }
    let mut yield_claim = scalar(HILL_REFERENCE_YIELD_STRESS_PROPERTY, Pressure::DIMS, YIELD);
    yield_claim.provenance.source = "manufactured anisotropic thermal shell".into();
    yield_claim.value = PropertyValue::Curve {
        abscissa: "T".into(),
        abscissa_dims: Temperature::DIMS,
        knots: vec![(280.0, 5e6), (300.0, YIELD), (320.0, 1e6)],
        dims: Pressure::DIMS,
    };
    yield_claim.interpolation = InterpolationPolicy::LinearInside;
    yield_claim.validity =
        ValidityDomain::unconstrained().with_quantity("T", temperature_spec(), 280.0, 320.0);
    claims.insert_claim(yield_claim).unwrap();
    // No isotropic J2 yield/hardening claim is supplied. The original Hill
    // reference radius and anisotropic flow metric own this material state.
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured anisotropic metal".into(),
            phase: "solid".into(),
            process: "test".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let caloric = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new(CP, CP_DIMS),
            density: &PropertyKey::new(DENSITY_PROPERTY, Density::DIMS),
            temperature_axis: "T",
            lower: &point(280.0),
            upper: &point(320.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        -20.0 * CP_HILL,
        1e-5,
    )
    .unwrap();
    Fixture {
        card,
        caloric,
        expansion: false,
    }
}

fn source(fixture: &Fixture, angle: f64) -> ThermalHillSource<'_> {
    let (s, c) = angle.sin_cos();
    ThermalHillSource::try_new(
        ThermalHillSourceQuery {
            card: &fixture.card,
            caloric: &fixture.caloric,
            context: point(300.0),
            temperature_axis: "T".into(),
            selection: MaterialPropertySelection::SingleClaimOnly,
            thermal_expansion: None,
            strain_limit: 0.04,
            material_to_spatial: [[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]],
        },
        0.0,
    )
    .unwrap()
}

fn endpoint(old_temperature: f64, old_alpha: f64, a: f64, stress: f64) -> (f64, f64) {
    // Under free-sided uniaxial stress q=a|sigma| and alpha=(q-Y(T))/H.
    // Only the increment above retained alpha produces new yield heat.
    let residual = |temperature: f64| {
        let yield_pa = YIELD - SOFTENING * (temperature - 300.0);
        let increment = ((a * stress - yield_pa) / H - old_alpha).max(0.0);
        RHO * CP_HILL * (temperature - old_temperature - EXTERNAL_RISE_K) - yield_pa * increment
    };
    let (mut lower, mut upper) = (old_temperature, old_temperature + 5.0);
    assert!(residual(lower) < 0.0 && residual(upper) > 0.0);
    for _ in 0..80 {
        let mid = 0.5 * (lower + upper);
        if residual(mid) > 0.0 {
            upper = mid;
        } else {
            lower = mid;
        }
    }
    let temperature = 0.5 * (lower + upper);
    let yield_pa = YIELD - SOFTENING * (temperature - 300.0);
    (temperature, old_alpha.max((a * stress - yield_pa) / H))
}

#[test]
fn g1_rolled_and_oblique_shells_heat_soften_and_release_with_original_hill_memory() {
    context(|cx| {
        let fixture = fixture();
        let mesh = mesh();
        let heat = RHO * VOLUME * CP_HILL * EXTERNAL_RISE_K;
        // Axial extension, lateral contraction, plastic shear and the three
        // drilling spins follow native equilibrium. Bending stays planar.
        let supports: Vec<_> = (0..18)
            .filter(|d| ![5, 6, 7, 11, 13, 17].contains(d))
            .map(|d| (d, 0.0))
            .collect();
        let mut permanent_extension = [0.0; 2];
        for (direction, (angle, a2, local_normal, spatial_normal)) in [
            (
                0.0,
                1.8_f64,
                [1.8, -0.6, -1.2, 0.0, 0.0, 0.0],
                [1.8, -0.6, -1.2, 0.0, 0.0, 0.0],
            ),
            (
                std::f64::consts::FRAC_PI_4,
                0.9,
                [0.6, 0.1, -0.7, -0.55, 0.0, 0.0],
                [0.9, -0.2, -0.7, 0.25, 0.0, 0.0],
            ),
        ]
        .into_iter()
        .enumerate()
        {
            let mut thermal = thermal(&fixture, &[heat / DT], 0.0);
            let layers = vec![vec![ShellThermoplasticLayer {
                cell: CellId(100),
                layer: PlaneStressLayer {
                    bottom_z_m: -THICKNESS / 2.0,
                    top_z_m: THICKNESS / 2.0,
                    subdivisions: 2,
                    law: source(&fixture, angle),
                },
            }]];
            let mut shell =
                ThermoplasticShell::try_new(&thermal, &mesh, layers, discretization(), cx).unwrap();
            let (mut old_temperature, mut old_alpha) = (300.0, 0.0);
            for (interval, stress) in [STRESS, STRESS, 0.0].into_iter().enumerate() {
                let target = (interval + 1) as f64 * DT;
                let old_solution = shell.solution().clone();
                let old_cells = thermal.cells().clone();
                let (expected_temperature, alpha) =
                    endpoint(old_temperature, old_alpha, a2.sqrt(), stress);
                let step = shell
                    .advance(
                        &thermal,
                        &(),
                        target,
                        &supports,
                        &[(6, stress * VOLUME / L)],
                        IncrementSettings::default(),
                        thermal_config(),
                        config(),
                        cx,
                    )
                    .unwrap();
                let temperature = step.thermal().temperature_k(CellId(100)).unwrap();
                close(temperature, expected_temperature, 5e-8);
                let plastic = spatial_normal.map(|normal| alpha * normal / a2.sqrt());
                let axial = stress / E + plastic[0];
                let transverse = -0.2 * stress / E + plastic[1];
                close(step.shell().solution().coordinates[1][0], L * axial, 3e-11);
                close(
                    step.shell().solution().coordinates[2][1],
                    L * transverse,
                    3e-11,
                );
                close(
                    step.shell().solution().coordinates[1][1],
                    2.0 * L * plastic[3],
                    3e-11,
                );
                for section in &step.shell().solution().section_states {
                    for point in &section.points {
                        close(point.history.alpha, alpha, 3e-9);
                        for (actual, normal) in
                            point.history.plastic_strain.iter().zip(local_normal)
                        {
                            close(*actual, alpha * normal / a2.sqrt(), 3e-9);
                        }
                        close(point.strain[2], -0.2 * stress / E + plastic[2], 3e-9);
                    }
                }
                let energy = &step.report().mechanical;
                let yield_pa = YIELD - SOFTENING * (expected_temperature - 300.0);
                let generated = VOLUME * yield_pa * (alpha - old_alpha);
                close(energy.physical_heat_j, generated, 3e-10);
                close(energy.creep_heat_j, 0.0, 0.0);
                close(energy.reversible_heat_j, 0.0, 1e-12);
                close(energy.coefficient_storage_change_j, 0.0, 1e-12);
                close(step.report().external_heat_j, heat, 3e-10);
                close(
                    step.report().cell_heat_j[&CellId(100)].total_j(),
                    generated,
                    3e-10,
                );
                close(
                    RHO * VOLUME * CP_HILL * (temperature - old_temperature),
                    heat + generated,
                    3e-10,
                );
                assert!(step.report().combined_residual_j.abs() <= config().energy_tolerance_j);
                assert_eq!(shell.solution(), &old_solution);
                assert_eq!(thermal.cells(), &old_cells);
                if interval == 1 {
                    assert!(alpha > old_alpha && generated > 0.0);
                } else if interval == 2 {
                    close(alpha, old_alpha, 0.0);
                    close(energy.physical_heat_j, 0.0, 1e-12);
                    permanent_extension[direction] = step.shell().solution().coordinates[1][0];
                }
                let (next_shell, next_thermal, _, _, _) = step.into_parts();
                shell = next_shell;
                thermal = next_thermal;
                old_temperature = expected_temperature;
                old_alpha = alpha;
            }
        }
        assert!(permanent_extension[0] > permanent_extension[1]);
        assert!(permanent_extension[1] > 0.0);
    });
}
