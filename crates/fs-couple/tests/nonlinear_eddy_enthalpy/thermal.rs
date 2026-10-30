//! Mixed material nonlinear field and physical heat against a scalar oracle.
use super::*;
#[path = "../../../fs-material/tests/nonlinear_magnetic_thermal/common.rs"]
mod source;
use fs_couple::circuit_enthalpy::magnetic::caloric::{
    MagneticConstitutiveThermalLaw, MagneticThermalLaw,
};
use fs_em::planar::{eddy::EddyConductor, nonlinear::IsotropicMagneticMaterial};
use fs_feec::weighted::WeightedAssemblyLimits;
use fs_material::{
    conductor::ResolvedConductorState, field_storage::thermal::magnetic::ThermalMagneticMaterial,
};

fn limits() -> WeightedAssemblyLimits {
    WeightedAssemblyLimits {
        max_cells: 1000,
        max_dofs: 3000,
        max_triplets: 10000,
    }
}
fn setup(
    slope: f64,
) -> (
    MaterialCard,
    ResolvedHeatCapacityCurve,
    PlanarNonlinearMagnetostatics,
    Vec<MagneticThermalLaw<MagneticConstitutiveThermalLaw>>,
) {
    let card = source::card(slope);
    let curve = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
            density: &PropertyKey::new("density", Density::DIMS),
            temperature_axis: "T",
            lower: &source::point(200.0),
            upper: &source::point(500.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        1e-9,
    )
    .unwrap();
    let nonlinear = source::law(&card);
    let linear =
        ThermalMagneticMaterial::new(&card, &source::bundle(&card, 300.0), "T", [290.0, 360.0])
            .unwrap();
    let mut materials = Vec::new();
    let mut laws = Vec::new();
    for triangle in 0..4 {
        let law = if triangle < 2 {
            materials.push(IsotropicMagneticMaterial::Reversible(
                nonlinear.at_temperature(300.0).unwrap().magnetic().clone(),
            ));
            MagneticConstitutiveThermalLaw::Separable(nonlinear.clone())
        } else {
            materials.push(IsotropicMagneticMaterial::Linear(
                linear.at_temperature(300.0).unwrap().material().clone(),
            ));
            MagneticConstitutiveThermalLaw::Linear(linear.clone())
        };
        laws.push(MagneticThermalLaw {
            triangle,
            law,
            maximum_flux_density_t: 4.0,
        });
    }
    let problem =
        PlanarNonlinearMagnetostatics::new(field::mesh(VOLUME), materials, limits(), &mut || false)
            .unwrap();
    (card, curve, problem, laws)
}
fn retained<'a>(
    card: &MaterialCard,
    problem: &'a PlanarNonlinearMagnetostatics,
) -> PlanarNonlinearEddy<'a> {
    let conductor = ResolvedConductorState::try_from_material_with_basis(
        source::bundle(card, 300.0),
        ElectricalPropertyBasis::Resistivity,
    )
    .unwrap();
    PlanarNonlinearEddy::new(
        &field::initial(problem),
        (0..4)
            .map(|triangle| EddyConductor {
                triangle,
                region: 7,
                conductor: conductor.clone(),
            })
            .collect(),
        &[],
        &field::boundary(0.0),
        0.0,
        field::controls(),
        &mut || false,
    )
    .unwrap()
}

// Independently eliminate the square's one free P1 A coordinate. Exactly half
// the volume is a saturating source-separable law; half is linear with mu=1.
// The total conductor covariance is V sigma/18. This oracle evaluates authored
// formulae and bracketed roots, never production source/field/thermal solvers.
fn oracle(a0: f64, slope: f64) -> (f64, f64, f64, f64, f64) {
    let scale = |t: f64| 1.0 + slope * (t - 300.0);
    let h = |b: f64, t: f64| 0.5 * (scale(t) * source::h(b) + b);
    let psi = |b: f64, t: f64| VOLUME * 0.5 * (scale(t) * source::psi(b) + 0.5 * b * b);
    let amplitude = |t: f64| {
        let c = VOLUME / (1.0 + 0.005 * (t - 300.0));
        let mut lo = 0.0;
        let mut hi = a0;
        for _ in 0..100 {
            let a = f64::midpoint(lo, hi);
            if c * (a - a0) / DT + 2.0 * VOLUME * h(2.0 * a, t) > 0.0 {
                hi = a;
            } else {
                lo = a;
            }
        }
        f64::midpoint(lo, hi)
    };
    let capacity = 6.0 * VOLUME;
    let mut lo = 290.0;
    let mut hi = 360.0;
    for _ in 0..100 {
        let t = f64::midpoint(lo, hi);
        let a = amplitude(t);
        let joule = VOLUME / (1.0 + 0.005 * (t - 300.0)) * (a - a0).powi(2) / DT;
        let reversible = 0.5 * VOLUME * slope * t * (source::psi(2.0 * a) - source::psi(2.0 * a0));
        if capacity * (t - 300.0) > joule + reversible {
            hi = t;
        } else {
            lo = t;
        }
    }
    let t = f64::midpoint(lo, hi);
    let a = amplitude(t);
    let joule = VOLUME / (1.0 + 0.005 * (t - 300.0)) * (a - a0).powi(2) / DT;
    let reversible = 0.5 * VOLUME * slope * t * (source::psi(2.0 * a) - source::psi(2.0 * a0));
    let numerical = VOLUME * h(2.0 * a, t) * 2.0 * (a - a0) - psi(2.0 * a, t) + psi(2.0 * a0, t);
    (t, a, joule, reversible, numerical)
}

#[test]
fn mixed_linear_and_nonlinear_sources_share_actual_field_entropy_and_finite_heat() {
    let mut outcomes = Vec::new();
    for slope in [-0.001, 0.0, 0.001] {
        let (card, curve, problem, laws) = setup(slope);
        let thermal = thermal(&curve);
        let state = retained(&card, &problem);
        let context = source::point(300.0);
        let cells = cells(&card, &curve, &context);
        let routes = routes();
        let network = MagneticThermalNetwork {
            cells: &cells,
            links: &[],
            magnetic_regions: &routes,
            volume_tolerance_m3: 1e-12,
        };
        let result = network
            .step_nonlinear_eddy_thermal_candidate(
                &thermal,
                DT,
                &state,
                &[0.0; 4],
                &laws,
                controls(),
                field::controls(),
                limits(),
                &[0.0],
                || false,
            )
            .unwrap();
        let (t, a, joule, reversible, numerical) = oracle(state.potential_wb_m()[4], slope);
        near(result.thermal().temperature_k(CELL).unwrap(), t);
        near(result.field().potential_wb_m()[4], a);
        let report = result.report();
        near(report.coupled.deposited_joule_j, joule);
        near(report.coupled.source.reversible_heat_j, reversible);
        near(report.coupled.electrical.numerical_dissipation_j, numerical);
        assert!(numerical > 0.0);
        near(
            report.coupled.thermal.stored_energy_change_j,
            joule + reversible,
        );
        near(report.total_energy_residual_j, 0.0);
        near(report.external_heat_j, 0.0);
        let initial_b = 2.0 * state.potential_wb_m()[4];
        for triangle in 0..4 {
            let field = &report.coupled.source.fields[triangle];
            near(field.temperature_k[1], t);
            near(field.flux_density_t[1], 2.0 * a);
            if triangle < 2 {
                near(
                    field.state[1].entropy_j_k,
                    -slope * source::psi(2.0 * a) * VOLUME / 4.0,
                );
                near(
                    field.parameter_work_j,
                    slope * (t - 300.0) * source::psi(initial_b) * VOLUME / 4.0,
                );
                near(field.reversible_heat_j, reversible / 2.0);
            } else {
                assert_eq!(field.state[1].entropy_j_k, 0.0);
                assert_eq!(field.reversible_heat_j, 0.0);
            }
            assert_eq!(
                report.coupled.source.material.magnetic[triangle]
                    .material
                    .properties()
                    .len(),
                7
            );
        }
        // Original source receipts and actual A/T histories continue together.
        let replay = network
            .step_nonlinear_eddy_thermal_candidate(
                &thermal,
                DT,
                &state,
                &[0.0; 4],
                &laws,
                controls(),
                field::controls(),
                limits(),
                &[0.0],
                || false,
            )
            .unwrap();
        assert_eq!(
            replay.field().potential_wb_m(),
            result.field().potential_wb_m()
        );
        assert_eq!(replay.thermal().cells(), result.thermal().cells());
        if slope == 0.001 {
            // The complete heat owner enters this same nonlinear residual;
            // empty interfaces reproduce the bare assembly exactly. Contact
            // and radiation interventions have their own full-heat tests.
            use fs_couple::circuit_enthalpy::{
                contacts::{heat::MaterialHeatNetwork, material::ElectricalContactNetwork},
                material::MaterialNetwork,
            };
            let bulk = MaterialNetwork {
                cells: &cells,
                resistors: &[],
                links: &[],
            };
            let electrical = ElectricalContactNetwork {
                bulk: &bulk,
                contacts: &[],
                electrical_bases: &[],
            };
            let heat = MaterialHeatNetwork {
                electrical: &electrical,
                thermal_contacts: &[],
                radiation: None,
            };
            let gate = fs_exec::CancelGate::new();
            let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
            pool.scope(|arena| {
                let cx = fs_exec::Cx::new(
                    &gate,
                    arena,
                    fs_exec::StreamKey {
                        seed: 59,
                        kernel_id: 73,
                        tile: 0,
                        iteration: 0,
                    },
                    fs_exec::Budget::INFINITE,
                    fs_exec::ExecMode::Deterministic,
                );
                let prepared = heat.prepare(&cx, &thermal).unwrap();
                let full = prepared
                    .step_nonlinear_eddy_thermal_candidate(
                        &cx,
                        &thermal,
                        DT,
                        &state,
                        &[0.0; 4],
                        &routes,
                        1e-12,
                        &laws,
                        controls(),
                        field::controls(),
                        limits(),
                        &[0.0],
                    )
                    .unwrap();
                assert_eq!(
                    full.field().potential_wb_m(),
                    result.field().potential_wb_m()
                );
                assert_eq!(full.thermal().cells(), result.thermal().cells());
                near(full.report().total_energy_residual_j, 0.0);
            });
        }
        let second = network
            .step_nonlinear_eddy_thermal_candidate(
                result.thermal(),
                2.0 * DT,
                result.field(),
                &[0.0; 4],
                &laws,
                controls(),
                field::controls(),
                limits(),
                &[0.0],
                || false,
            )
            .unwrap();
        near(second.report().total_energy_residual_j, 0.0);
        assert_eq!(second.field().time_s(), 2.0 * DT);
        if slope == 0.0 {
            let fixed = network
                .step_nonlinear_eddy_candidate(
                    &thermal,
                    DT,
                    &state,
                    &[0.0; 4],
                    controls(),
                    field::controls(),
                    &[0.0],
                    || false,
                )
                .unwrap();
            near(fixed.field().potential_wb_m()[4], a);
            near(fixed.thermal().temperature_k(CELL).unwrap(), t);
        } else {
            // A fixed-H heat consumer cannot drop nonzero field entropy, even
            // when its very first trial still has the accepted temperature.
            assert!(
                network
                    .step_nonlinear_eddy_candidate(
                        &thermal,
                        DT,
                        &state,
                        &[0.0; 4],
                        controls(),
                        field::controls(),
                        &[0.0],
                        || false
                    )
                    .is_err()
            );
            assert!(reversible * slope < 0.0);
        }
        outcomes.push((t, a));
    }
    assert!(outcomes[0].0 > outcomes[1].0 && outcomes[1].0 > outcomes[2].0);
    assert!(outcomes[0].1 < outcomes[1].1 && outcomes[1].1 < outcomes[2].1);
}

#[test]
fn nonlinear_entropy_candidates_refuse_source_mismatch_and_cancel_without_publishing() {
    let (card, curve, problem, laws) = setup(0.001);
    let thermal = thermal(&curve);
    let state = retained(&card, &problem);
    let context = source::point(300.0);
    let cells = cells(&card, &curve, &context);
    let routes = routes();
    let network = MagneticThermalNetwork {
        cells: &cells,
        links: &[],
        magnetic_regions: &routes,
        volume_tolerance_m3: 1e-12,
    };
    let run = |bindings: &[MagneticThermalLaw<MagneticConstitutiveThermalLaw>]| {
        network.step_nonlinear_eddy_thermal_candidate(
            &thermal,
            DT,
            &state,
            &[0.0; 4],
            bindings,
            controls(),
            field::controls(),
            limits(),
            &[0.0],
            || false,
        )
    };
    assert!(run(&laws[..3]).is_err());
    let mut wrong = laws.clone();
    wrong[0].law = wrong[2].law.clone();
    assert!(run(&wrong).is_err());
    let mut narrow = laws.clone();
    narrow[0].maximum_flux_density_t = 0.1;
    assert!(run(&narrow).is_err());
    let other = setup(-0.001).3;
    assert!(run(&other).is_err());
    let mut polls = 0usize;
    let accepted = network
        .step_nonlinear_eddy_thermal_candidate(
            &thermal,
            DT,
            &state,
            &[0.0; 4],
            &laws,
            controls(),
            field::controls(),
            limits(),
            &[0.0],
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let mut count = 0usize;
    let cancelled = network.step_nonlinear_eddy_thermal_candidate(
        &thermal,
        DT,
        &state,
        &[0.0; 4],
        &laws,
        controls(),
        field::controls(),
        limits(),
        &[0.0],
        || {
            count += 1;
            count >= polls - 1
        },
    );
    assert!(matches!(cancelled, Err(CircuitEnthalpyError::Cancelled)));
    assert_eq!(state.time_s(), 0.0);
    near(thermal.temperature_k(CELL).unwrap(), 300.0);
    let retried = run(&laws).unwrap();
    assert_eq!(
        retried.field().potential_wb_m(),
        accepted.field().potential_wb_m()
    );
    assert_eq!(retried.thermal().cells(), accepted.thermal().cells());
}
