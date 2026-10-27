//! Synthetic material data exercise the native VOF, not experimental accuracy.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_lbm::d3q19::freesurface3::ContactModel3;
use fs_lbm::d3q19::{E3, Q3, equilibrium3};
use fs_lbm::material::free_surface::rheology::{
    FreeSurfaceRheologyError, IsothermalRheologyLaw, MaterialFreeSurfaceRheology3, RheologySolve,
};
use fs_lbm::material::free_surface::{
    MaterialFreeSurfaceDomain3, MaterialFreeSurfaceInitialCell3, plan_free_surface_from_material,
};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::liquid::surface::{LiquidCapillarity, ResolvedLiquidSurfaceState};
use fs_material::liquid_temperature::LiquidViscosityBasis;
use fs_material::state_point::{
    MaterialPropertySelection, ResolvedMaterialStatePoint, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_material_state_point,
};
use fs_qty::{Dims, Frequency, QuantitySpec};

const CELL: [usize; 3] = [3; 3];
const DX: f64 = 0.001;
const DT: f64 = 0.001;

fn solve() -> RheologySolve {
    RheologySolve {
        maximum_iterations: 128,
        relaxation: 0.5,
        absolute_relaxation_tolerance: 1e-13,
        relative_relaxation_tolerance: 1e-13,
    }
}

fn card(slope: f64, max_rate: f64, sigma: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        ("density", Dims([-3, 1, 0, 0, 0, 0]), 1000.0),
        ("dynamic_viscosity", Dims([-1, 1, -1, 0, 0, 0]), 0.1),
        ("surface_tension", Dims([0, 1, -2, 0, 0, 0]), sigma),
    ] {
        let evolving = name == "dynamic_viscosity" && slope != 0.0;
        let mut validity = ValidityDomain::unconstrained()
            .with("T", 300.0, 500.0)
            .with("pressure", 100000.0, 100000.0);
        let value = if evolving {
            validity = validity.with_quantity(
                "rate",
                QuantitySpec::dimensional(Frequency::DIMS),
                0.0,
                max_rate,
            );
            PropertyValue::Curve {
                abscissa: "rate".into(),
                abscissa_dims: Frequency::DIMS,
                knots: vec![(0.0, value), (max_rate, value + slope * max_rate)],
                dims,
            }
        } else {
            PropertyValue::Scalar { value, dims }
        };
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value,
                validity,
                interpolation: if evolving {
                    InterpolationPolicy::LinearInside
                } else {
                    InterpolationPolicy::ConstantWithinValidity
                },
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic isothermal rheology".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic liquid".into(),
            phase: "liquid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn bundle(card: &MaterialCard) -> ResolvedMaterialStatePoint {
    let point = QueryPoint::new()
        .with("T", 400.0)
        .unwrap()
        .with("pressure", 100000.0)
        .unwrap()
        .with_quantity("rate", QuantitySpec::dimensional(Frequency::DIMS), 0.0)
        .unwrap();
    let requirements = [
        ("density", Dims([-3, 1, 0, 0, 0, 0])),
        ("dynamic_viscosity", Dims([-1, 1, -1, 0, 0, 0])),
        ("surface_tension", Dims([0, 1, -2, 0, 0, 0])),
    ]
    .map(|(name, dims)| {
        ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::StrictlyPositive)
            .unwrap()
    });
    resolve_material_state_point(
        card,
        &point,
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

fn law(card: &MaterialCard) -> IsothermalRheologyLaw<'_> {
    IsothermalRheologyLaw::try_new(
        card,
        bundle(card),
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Sourced,
        "T",
        "rate",
    )
    .unwrap()
}

fn domain(gravity: [f64; 3]) -> MaterialFreeSurfaceDomain3<'static> {
    MaterialFreeSurfaceDomain3 {
        cells: [8; 3],
        gravity_m_s2: gravity,
        contact: ContactModel3::Neutral,
        solid_tiles: &[],
        max_cells: 512,
    }
}

fn initial(full: bool) -> Vec<MaterialFreeSurfaceInitialCell3> {
    let mut cells = Vec::new();
    for z in 0..8 {
        for y in 0..8 {
            for x in 0..8 {
                if full || ((1..=3).contains(&x) && (1..=6).contains(&y) && (1..=6).contains(&z)) {
                    cells.push(MaterialFreeSurfaceInitialCell3 {
                        coordinates: [x, y, z],
                        fill_fraction: 1.0,
                        velocity_m_s: [0.001 * ((z as f64 - 3.5).powi(2) - 5.25), 0.0, 0.0],
                    });
                }
            }
        }
    }
    cells
}

fn with_cx<R>(gate: &CancelGate, run: impl FnOnce(&Cx<'_>) -> R) -> R {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        run(&Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 41,
                kernel_id: 27,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual}, expected={expected}"
    );
}

// Independent non-equilibrium second-moment reduction. The full fixture has
// all eight Morton-ordered tiles active; no production stress helper is called.
fn numerator(flow: &MaterialFreeSurfaceRheology3<'_, '_>, gravity: [f64; 3]) -> f64 {
    let bits = flow.lattice().grid().state_bits();
    let [x, y, z] = CELL;
    let slot = (x / 4) | ((y / 4) << 1) | ((z / 4) << 2);
    let lane = x % 4 + 4 * (y % 4 + 4 * (z % 4));
    let sample = flow.lattice().sample_at(x, y, z).unwrap();
    let eq = equilibrium3(sample.density, sample.velocity);
    let force = gravity.map(|g| g * DT * DT / DX);
    let mut tensor = [[0.0; 3]; 3];
    for q in 0..Q3 {
        let e = [f64::from(E3[q].0), f64::from(E3[q].1), f64::from(E3[q].2)];
        let neq = f64::from_bits(bits[(slot * Q3 + q) * 64 + lane]) - eq[q];
        for a in 0..3 {
            for b in 0..3 {
                tensor[a][b] += e[a] * e[b] * neq;
            }
        }
    }
    for a in 0..3 {
        for b in 0..3 {
            tensor[a][b] += 0.5 * (sample.velocity[a] * force[b] + sample.velocity[b] * force[a]);
        }
    }
    3.0 * (0.5 * tensor.iter().flatten().map(|v| v * v).sum::<f64>()).sqrt() / sample.density
}

#[test]
fn native_stress_closes_source_viscosity_and_newtonian_limit_has_exact_native_parity() {
    let gate = CancelGate::new_clock_free();
    let cells = initial(true);
    with_cx(&gate, |cx| {
        for gravity in [[0.0; 3], [0.03, -0.02, 0.01]] {
            let source = card(0.003, 1000.0, 0.001);
            let binding = law(&source);
            let mut flow = MaterialFreeSurfaceRheology3::from_cells(
                &binding,
                DX,
                DT,
                0.1,
                domain(gravity),
                &cells,
                solve(),
                cx,
                &gate,
            )
            .unwrap();
            flow.step_with_gate(&gate).unwrap();
            // mu = 0.1 + b*rate; tau = 0.8 + 3*b*rate. Solve this quadratic
            // analytically instead of reproducing the production iteration.
            let rhs = numerator(&flow, gravity) / DT;
            let expected = 2.0 * rhs / (0.8 + (0.64 + 12.0 * 0.003 * rhs).sqrt());
            let sample = flow.rheology_at(CELL).unwrap();
            near(sample.shear_rate_per_s, expected, 2e-10);
            near(sample.collision_tau, 0.8 + 3.0 * 0.003 * expected, 2e-11);
            let material = flow.material_at(CELL).unwrap();
            near(
                material.property("dynamic_viscosity").unwrap().value_si(),
                0.1 + 0.003 * expected,
                2e-11,
            );
            assert_eq!(material.card_identity(), source.content_hash());
            assert!(material.query_point().contains(&("T".into(), 400.0)));
            assert!(
                material
                    .query_point()
                    .contains(&("pressure".into(), 100000.0))
            );
        }
        let source = card(0.0, 1000.0, 0.001);
        let binding = law(&source);
        let frozen = ResolvedLiquidSurfaceState::try_from_material(
            bundle(&source),
            LiquidViscosityBasis::Dynamic,
            LiquidCapillarity::Sourced,
        )
        .unwrap();
        let mut native = plan_free_surface_from_material(&frozen, DX, DT, 0.1)
            .unwrap()
            .free_surface3_from_cells(domain([0.0; 3]), &cells, cx)
            .unwrap();
        let mut rheology = MaterialFreeSurfaceRheology3::from_cells(
            &binding,
            DX,
            DT,
            0.1,
            domain([0.0; 3]),
            &cells,
            solve(),
            cx,
            &gate,
        )
        .unwrap();
        for _ in 0..4 {
            native.step_with_gate(&gate).unwrap();
            rheology.step_with_gate(&gate).unwrap();
            assert_eq!(
                native.lattice().grid().state_bits(),
                rheology.lattice().grid().state_bits()
            );
            assert_eq!(
                native.ledger_mass_kg().unwrap(),
                rheology.ledger_mass_kg().unwrap()
            );
        }
    });
}

#[test]
fn source_only_viscosity_and_capillarity_changes_reach_native_motion_and_preserve_mass() {
    let gate = CancelGate::new_clock_free();
    let cells = initial(false);
    with_cx(&gate, |cx| {
        let mut results = Vec::new();
        for (slope, sigma) in [
            (0.0, 0.001),
            (0.003, 0.001),
            (-0.00005, 0.001),
            (0.003, 0.002),
        ] {
            let source = card(slope, 1000.0, sigma);
            let binding = law(&source);
            let mut flow = MaterialFreeSurfaceRheology3::from_cells(
                &binding,
                DX,
                DT,
                0.1,
                domain([0.0, 0.0, -0.02]),
                &cells,
                solve(),
                cx,
                &gate,
            )
            .unwrap();
            let mut replay = MaterialFreeSurfaceRheology3::from_cells(
                &binding,
                DX,
                DT,
                0.1,
                domain([0.0, 0.0, -0.02]),
                &cells,
                solve(),
                cx,
                &gate,
            )
            .unwrap();
            let mass = flow.ledger_mass_kg().unwrap();
            assert!(flow.capillary_density_excursion_bound() > 0.0);
            for _ in 0..4 {
                flow.step_with_gate(&gate).unwrap();
                replay.step_with_gate(&gate).unwrap();
                near(flow.ledger_mass_kg().unwrap(), mass, 1e-10 * mass);
            }
            assert_eq!(
                flow.lattice().grid().state_bits(),
                replay.lattice().grid().state_bits()
            );
            results.push((
                flow.lattice().grid().state_bits(),
                flow.sample_at_si(2, 3, 3).unwrap().unwrap().velocity_m_s,
            ));
        }
        for (a, b) in [(0, 1), (0, 2), (1, 3)] {
            assert_ne!(
                results[a].0, results[b].0,
                "source data must reach actual collision/boundary"
            );
            assert!(
                results[a]
                    .1
                    .iter()
                    .zip(results[b].1)
                    .any(|(x, y)| (x - y).abs() > 1e-10),
                "source data must change physical motion"
            );
        }
    });
}

#[test]
fn endpoint_source_and_convergence_refusals_and_cancellation_preserve_the_accepted_state() {
    let gate = CancelGate::new_clock_free();
    let cells = initial(true);
    with_cx(&gate, |cx| {
        for (max_rate, maximum_iterations) in [(1e-6, 128), (1000.0, 1)] {
            let source = card(0.003, max_rate, 0.001);
            let binding = law(&source);
            let mut flow = MaterialFreeSurfaceRheology3::from_cells(
                &binding,
                DX,
                DT,
                0.1,
                domain([0.0; 3]),
                &cells,
                RheologySolve {
                    maximum_iterations,
                    ..solve()
                },
                cx,
                &gate,
            )
            .unwrap();
            let bits = flow.lattice().grid().state_bits();
            let mass = flow.ledger_mass_kg().unwrap();
            let carry = flow.lattice().carry();
            let conversions = flow.lattice().conversions();
            let sample = *flow.rheology_at(CELL).unwrap();
            let receipt = flow.material_at(CELL).unwrap().clone();
            let refusal = flow.step_with_gate(&gate).unwrap_err();
            if maximum_iterations == 1 {
                assert!(matches!(
                    refusal,
                    FreeSurfaceRheologyError::Convergence { .. }
                ));
            } else {
                assert!(matches!(refusal, FreeSurfaceRheologyError::Material(_)));
            }
            assert_eq!(flow.lattice().grid().state_bits(), bits);
            assert_eq!(flow.ledger_mass_kg().unwrap(), mass);
            assert_eq!(flow.lattice().carry(), carry);
            assert_eq!(flow.lattice().conversions(), conversions);
            assert_eq!(*flow.rheology_at(CELL).unwrap(), sample);
            assert_eq!(*flow.material_at(CELL).unwrap(), receipt);
            let cancelled = CancelGate::new_clock_free();
            cancelled.request();
            assert!(flow.step_with_gate(&cancelled).is_err());
            assert_eq!(flow.lattice().grid().state_bits(), bits);
        }
        let source = card(0.003, 1000.0, 0.001);
        let binding = law(&source);
        // The shared capillary gate rejects nonpositive reconstructed density.
        assert!(
            MaterialFreeSurfaceRheology3::from_cells(
                &binding,
                DX,
                0.1,
                0.0,
                domain([0.0; 3]),
                &[],
                solve(),
                cx,
                &gate,
            )
            .is_err()
        );
        // The isothermal path refuses a changing storage or capillary law.
        for name in ["density", "surface_tension"] {
            let mut claims = ClaimSet::new();
            for (_, original) in source.claims().claims_ordered() {
                let mut claim = original.clone();
                if claim.key.name() == name {
                    let value = if name == "density" { 1000.0 } else { 0.001 };
                    claim.value = PropertyValue::Curve {
                        abscissa: "rate".into(),
                        abscissa_dims: Frequency::DIMS,
                        knots: vec![(0.0, value), (1000.0, value)],
                        dims: claim.key.dims(),
                    };
                    claim.interpolation = InterpolationPolicy::LinearInside;
                    claim.validity = claim.validity.with_quantity(
                        "rate",
                        QuantitySpec::dimensional(Frequency::DIMS),
                        0.0,
                        1000.0,
                    );
                }
                claims.insert_claim(claim).unwrap();
            }
            let changed = MaterialCard::assemble(source.id().clone(), claims, vec![]).unwrap();
            assert!(
                IsothermalRheologyLaw::try_new(
                    &changed,
                    bundle(&changed),
                    LiquidViscosityBasis::Dynamic,
                    LiquidCapillarity::Sourced,
                    "T",
                    "rate"
                )
                .is_err()
            );
        }
    });
}
