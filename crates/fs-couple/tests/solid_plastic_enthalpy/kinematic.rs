//! Sourced backstress softening changes native cyclic motion and reciprocal heat.
use super::*;
use fs_material::plastic::source::KINEMATIC_HARDENING_MODULUS_PROPERTY;

const HK: f64 = 2e8;
const CP_VALUE: f64 = 100.0;
const DURATION: f64 = 0.001;
const HEAT_PER_VOLUME: f64 = RHO * CP_VALUE * 4.0;

fn fixture(slope: f64) -> Fixture {
    let mut kinematic = scalar(KINEMATIC_HARDENING_MODULUS_PROPERTY, Pressure::DIMS, HK);
    kinematic.value = PropertyValue::Curve {
        abscissa: "T".into(),
        abscissa_dims: Temperature::DIMS,
        knots: [280.0, 300.0, 310.0, 320.0]
            .map(|t| (t, HK + slope * (t - 300.0)))
            .to_vec(),
        dims: Pressure::DIMS,
    };
    kinematic.interpolation = InterpolationPolicy::LinearInside;
    kinematic.validity = ValidityDomain::unconstrained().with("T", 280.0, 320.0);
    let mut claims = ClaimSet::new();
    for claim in [
        scalar(YOUNG_MODULUS_PROPERTY, Pressure::DIMS, E),
        scalar(POISSON_RATIO_PROPERTY, Dims([0; 6]), 0.2),
        scalar(YIELD_STRESS_PROPERTY, Pressure::DIMS, YIELD),
        scalar(BASIS.property_name(), Pressure::DIMS, H),
        scalar(DENSITY_PROPERTY, Density::DIMS, RHO),
        scalar(CP, CP_DIMS, CP_VALUE),
        kinematic,
    ] {
        claims.insert_claim(claim).unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic mixed-hardening thermal solid".into(),
            phase: "solid".into(),
            process: "manufactured cyclic test".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let lower = QueryPoint::new().with("T", 280.0).unwrap();
    let upper = QueryPoint::new().with("T", 320.0).unwrap();
    let caloric = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new(CP, CP_DIMS),
            density: &PropertyKey::new(DENSITY_PROPERTY, Density::DIMS),
            temperature_axis: "T",
            lower: &lower,
            upper: &upper,
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        -20.0 * CP_VALUE,
        1e-5,
    )
    .unwrap();
    Fixture {
        card,
        caloric,
        expansion: false,
    }
}

// Independent stress-controlled uniaxial reduction: beta = H_k p,
// |sigma - beta| = Y + H_iso alpha and dp = sign(sigma-beta) d(alpha).
#[derive(Clone, Copy)]
struct Reference {
    t: f64,
    sigma: f64,
    p: f64,
    accumulated: f64,
    hk: f64,
}
impl Reference {
    fn next(self, t: f64, sigma: f64, slope: f64) -> Self {
        let hk = HK + slope * (t - 300.0);
        let relative = sigma - hk * self.p;
        let increment = ((relative.abs() - YIELD - H * self.accumulated) / (H + hk)).max(0.0);
        Self {
            t,
            sigma,
            p: self.p + relative.signum() * increment,
            accumulated: self.accumulated + increment,
            hk,
        }
    }
    fn free_energy(self) -> f64 {
        0.5 * self.sigma * self.sigma / E
            + 0.5 * H * self.accumulated * self.accumulated
            + 0.5 * self.hk * self.p * self.p
    }
    fn entropy(self, slope: f64) -> f64 {
        -0.5 * slope * self.p * self.p
    }
    fn strain(self) -> f64 {
        self.sigma / E + self.p
    }
}
fn coefficient(old: Reference, new: Reference) -> f64 {
    0.5 * (new.hk - old.hk) * old.p * old.p
}
fn entropy_storage(old: Reference, new: Reference, slope: f64) -> f64 {
    new.t * new.entropy(slope) - old.t * old.entropy(slope)
}
fn reversible(old: Reference, new: Reference, slope: f64) -> f64 {
    -coefficient(old, new) - entropy_storage(old, new, slope)
}
fn endpoint(old: Reference, sigma: f64, slope: f64) -> Reference {
    let residual = |t| {
        let new = old.next(t, sigma, slope);
        RHO * CP_VALUE * (t - old.t)
            - HEAT_PER_VOLUME
            - YIELD * (new.accumulated - old.accumulated)
            - reversible(old, new, slope)
    };
    let (mut lower, mut upper) = (old.t, 320.0);
    assert!(residual(lower) < 0.0 && residual(upper) > 0.0);
    for _ in 0..80 {
        let mid = 0.5 * (lower + upper);
        if residual(mid) > 0.0 {
            upper = mid;
        } else {
            lower = mid;
        }
    }
    old.next(0.5 * (lower + upper), sigma, slope)
}

#[test]
fn sourced_backstress_softening_changes_native_cyclic_plastic_motion_and_reciprocal_heat() {
    context(|cx| {
        let mut responses = Vec::new();
        for slope in [0.0, -5e6] {
            let fixture = fixture(slope);
            let mut thermal = initial(&[&fixture], &[], &[HEAT_PER_VOLUME * V / DURATION]);
            let mut query = fixture.query();
            query.hardening_basis = HardeningPropertyBasis::MixedPlasticStrainModuli;
            let mut solid = ThermoplasticSolid::try_new(
                &thermal,
                &NODES,
                &TETS,
                vec![ThermoplasticCell {
                    cell: CellId(100),
                    query,
                }],
                TetAssemblyBudget::standard(),
                cx,
            )
            .unwrap();
            let mut old = Reference {
                t: 300.0,
                sigma: 0.0,
                p: 0.0,
                accumulated: 0.0,
                hk: HK,
            };
            for (i, sigma) in [3.4e6, -3.4e6, 3.4e6].into_iter().enumerate() {
                let original_solution = solid.solution().clone();
                let original_cells = thermal.cells().clone();
                let candidate = solid
                    .advance(
                        &thermal,
                        &(),
                        (i + 1) as f64 * DURATION,
                        &ROLLERS,
                        &[(3, sigma * V / L)],
                        IncrementSettings::default(),
                        thermal_config(),
                        config(),
                        cx,
                    )
                    .unwrap();
                let expected = endpoint(old, sigma, slope);
                close(
                    candidate.thermal().temperature_k(CellId(100)).unwrap(),
                    expected.t,
                    3e-8,
                );
                let solution = candidate.solid().solution();
                close(solution.element_states[0].alpha, expected.accumulated, 3e-9);
                for (axis, p) in [expected.p, -0.5 * expected.p, -0.5 * expected.p]
                    .into_iter()
                    .enumerate()
                {
                    close(solution.element_states[0].plastic_strain[axis], p, 3e-9);
                    close(
                        solution.element_stress[0][axis],
                        if axis == 0 { sigma } else { 0.0 },
                        3e-3,
                    );
                }
                close(solution.displacement_m[1][0], L * expected.strain(), 3e-11);
                let source = candidate.solid().materials()[0]
                    .kinematic_hardening_source()
                    .unwrap();
                close(source.value_si(), expected.hk, 0.3);
                fixture
                    .card
                    .claims()
                    .verify_receipt(&source.answer().receipt)
                    .unwrap();
                let report = candidate.report();
                close(
                    report.coefficient_storage_change_j,
                    V * coefficient(old, expected),
                    1e-8,
                );
                close(
                    report.thermoelastic_storage_change_j,
                    V * entropy_storage(old, expected, slope),
                    1e-8,
                );
                close(
                    report.reversible_heat_j,
                    V * reversible(old, expected, slope),
                    1e-8,
                );
                close(
                    report.stored_mechanical_change_j,
                    V * (expected.free_energy() - old.free_energy()),
                    1e-8,
                );
                close(
                    report.physical_heat_j,
                    V * YIELD * (expected.accumulated - old.accumulated),
                    1e-8,
                );
                close(
                    report.mechanical_work_j,
                    V * sigma * (expected.strain() - old.strain()),
                    1e-8,
                );
                let dp = expected.p - old.p;
                let da = expected.accumulated - old.accumulated;
                let ds = expected.sigma - old.sigma;
                close(
                    report.numerical_loss_j,
                    V * (0.5 * ds * ds / E + 0.5 * H * da * da + 0.5 * expected.hk * dp * dp),
                    1e-8,
                );
                close(report.combined_residual_j, 0.0, 1e-8);
                assert_eq!(candidate.thermal().cells()[&CellId(100)].mass_kg, RHO * V);
                assert_eq!(solid.solution(), &original_solution);
                assert_eq!(thermal.cells(), &original_cells);
                let (next_thermal, next_solid, _, _) = candidate.into_parts();
                thermal = next_thermal;
                solid = next_solid;
                old = expected;
            }
            responses.push(old);
        }
        assert!(responses[1].accumulated > responses[0].accumulated);
        assert!((responses[1].t - responses[0].t).abs() > 1e-3);
        assert!((responses[1].strain() - responses[0].strain()).abs() > 1e-5);
    });
}
