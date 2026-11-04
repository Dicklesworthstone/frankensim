//! Complete native damped-frame spectra against independent mechanical pencils.
//! Manufactured coefficients keep every expected pole analytic.

#[path = "frame2d_poles/coupled.rs"]
mod coupled;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::section::{ConcentricElasticSection, ConcentricMaterialLayer};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_material_state_point,
};
use fs_math::c64::C64;
use fs_qty::{Density, DynViscosity, Pressure};
use fs_solid::frame2d::{
    FrameMassFormulation, FrameViscousSection,
    poles::{FramePoleControls, FramePoleError},
};
use fs_solid::{Frame2d, Geometry, Rayleigh, SectionModel};

const RHO: f64 = 6.0;
const YOUNG: f64 = 8.0;
const RADIUS: f64 = 0.2;
const BAND: (f64, f64) = (0.5, 3.0);
const NO_RAYLEIGH: Rayleigh = Rayleigh { a0: 0.0, a1: 0.0 };

fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 91,
                kernel_id: 22,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

fn controls() -> FramePoleControls {
    FramePoleControls {
        maximum_coordinates: 12,
        maximum_payload_bytes: 1 << 20,
        maximum_qr_sweeps_per_window: 120,
        maximum_work: 1 << 36,
    }
}

fn source(eta_axial: f64, eta_bending: f64, band: (f64, f64)) -> FrameViscousSection {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (name, value, dims) in [
        ("density", RHO, Density::DIMS),
        ("young_modulus", YOUNG, Pressure::DIMS),
        ("axial_eta", eta_axial, DynViscosity::DIMS),
        ("bending_eta", eta_bending, DynViscosity::DIMS),
    ] {
        let mut validity = ValidityDomain::unconstrained().with("T", 300.0, 400.0);
        if name.ends_with("_eta") {
            validity = validity.with("omega", band.0, band.1);
        }
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity,
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured analytic frame-pole fixture".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
        requirements.push(
            ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::Finite).unwrap(),
        );
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured solid".into(),
            phase: "solid".into(),
            process: "linear viscosity fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let point = QueryPoint::new()
        .with("T", 300.0)
        .unwrap()
        .with("omega", 0.5 * (band.0 + band.1))
        .unwrap();
    let material = resolve_material_state_point(
        &card,
        &point,
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let section = ConcentricElasticSection::try_new(
        0.0,
        vec![ConcentricMaterialLayer {
            outer_radius_m: RADIUS,
            material,
        }],
    )
    .unwrap();
    FrameViscousSection::try_new(section, "axial_eta", "bending_eta").unwrap()
}

fn member(geometry: Geometry, formulation: FrameMassFormulation) -> Frame2d {
    let mut frame = Frame2d::with_mass_formulation(geometry, formulation);
    frame.add_node(0.0, 0.0);
    frame.add_node(1.0, 0.0);
    frame.fix(0, [true; 3]);
    frame
}

fn guided(eta: f64) -> Frame2d {
    let mut frame = member(Geometry::Linear, FrameMassFormulation::ConsistentReference);
    frame.fix(1, [false, true, true]);
    frame
        .add_material_viscous_element(0, 1, 5, source(eta, 0.0, BAND))
        .unwrap();
    frame
}

fn expect_rates(actual: &[C64], expected: &[C64], relative: f64, absolute: f64) {
    assert_eq!(
        actual.len(),
        expected.len(),
        "the complete 2n spectrum is required"
    );
    assert!(actual.windows(2).all(|pair| {
        pair[0]
            .re
            .total_cmp(&pair[1].re)
            .then_with(|| pair[0].im.total_cmp(&pair[1].im))
            .is_le()
    }));
    // QR may separate the real parts of conjugates by a few ulps. Match the
    // complete root set, without assuming those roundoff ties have one order.
    let mut unmatched = actual.to_vec();
    for &wanted in expected {
        let index = unmatched
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| (**a - wanted).abs().total_cmp(&(**b - wanted).abs()))
            .unwrap()
            .0;
        let found = unmatched.remove(index);
        let tolerance = absolute.max(relative * wanted.abs());
        assert!(
            (found - wanted).abs() <= tolerance,
            "actual={found:?}, expected={wanted:?}, tolerance={tolerance:.3e}"
        );
    }
}

#[test]
fn material_scalar_keeps_oscillatory_critical_and_both_overdamped_rates() {
    with_cx(&CancelGate::new(), |cx| {
        // Unit member, fixed base and guided tip: M=2A, K=8A, C=eta*A.
        // The second case adds Rayleigh C=A exactly once to material C=2A.
        let cases = [
            (
                0.0,
                NO_RAYLEIGH,
                [C64::new(0.0, -2.0), C64::new(0.0, 2.0)],
                0.0,
            ),
            (
                2.0,
                Rayleigh {
                    a0: 0.25,
                    a1: 0.0625,
                },
                [
                    C64::new(-0.75, -55.0_f64.sqrt() / 4.0),
                    C64::new(-0.75, 55.0_f64.sqrt() / 4.0),
                ],
                0.0,
            ),
            // Assembly roundoff perturbs a repeated root by O(sqrt(epsilon)).
            (8.0, NO_RAYLEIGH, [C64::new(-2.0, 0.0); 2], 3e-7),
            (
                20.0,
                NO_RAYLEIGH,
                [
                    C64::new(-5.0 - 21.0_f64.sqrt(), 0.0),
                    C64::new(-4.0 / (5.0 + 21.0_f64.sqrt()), 0.0),
                ],
                0.0,
            ),
            // Roots differ from these values by relative O(1e-24); a direct
            // subtraction of nearly equal discriminant terms loses the slow one.
            (
                8e12,
                NO_RAYLEIGH,
                [C64::new(-4e12, 0.0), C64::new(-1e-12, 0.0)],
                0.0,
            ),
        ];
        for (eta, rayleigh, expected, absolute) in cases {
            let mut frame = guided(eta);
            let mut force = vec![0.0; frame.ndof()];
            force[3] = YOUNG * core::f64::consts::PI * RADIUS.powi(2) * 1e-3;
            frame.static_load(&force, 1).unwrap();
            let accepted = format!("{frame:?}");
            let spectrum = frame.resonance_poles(rayleigh, controls(), cx).unwrap();
            assert_eq!(spectrum.free_dofs, [3]);
            // Retain the source band even for nonoscillatory and fast rates;
            // reporting the coefficients' poles does not extend material evidence.
            assert_eq!(spectrum.material_omega_band_rad_s, Some(BAND));
            expect_rates(&spectrum.rates_per_s, &expected, 2e-11, absolute);
            assert_eq!(format!("{frame:?}"), accepted);
        }
    });
}

#[test]
fn heterogeneous_member_losses_keep_the_full_nonproportional_pencil() {
    with_cx(&CancelGate::new(), |cx| {
        let mut frame = Frame2d::with_mass_formulation(
            Geometry::Linear,
            FrameMassFormulation::ConsistentReference,
        );
        for node in 0..3 {
            frame.add_node(node as f64, 0.0);
            frame.fix(
                node,
                if node == 0 {
                    [true; 3]
                } else {
                    [false, true, true]
                },
            );
        }
        frame
            .add_material_viscous_element(0, 1, 3, source(16.0, 0.0, (0.25, 4.0)))
            .unwrap();
        frame
            .add_material_viscous_element(1, 2, 7, source(1.0, 0.0, BAND))
            .unwrap();
        // Divide all matrices by common area A:
        // M=[[4,1],[1,2]], K=[[16,-8],[-8,8]], C=[[17,-1],[-1,1]].
        // Directly expanding the 2x2 determinant gives
        // 7s^4+40s^3+96s^2+136s+64
        // = (s^2+2s+4)(7s^2+26s+16).
        // Thus one physical loss distribution has both oscillatory and real
        // relaxation rates; diagonalizing C in an undamped basis changes them.
        let expected = [
            C64::new(-1.0, -3.0_f64.sqrt()),
            C64::new(-1.0, 3.0_f64.sqrt()),
            C64::new((-13.0 - 57.0_f64.sqrt()) / 7.0, 0.0),
            C64::new((-13.0 + 57.0_f64.sqrt()) / 7.0, 0.0),
        ];
        let spectrum = frame.resonance_poles(NO_RAYLEIGH, controls(), cx).unwrap();
        assert_eq!(spectrum.free_dofs, [3, 6]);
        assert_eq!(spectrum.material_omega_band_rad_s, Some(BAND));
        expect_rates(&spectrum.rates_per_s, &expected, 2e-10, 0.0);
    });
}

#[test]
fn bending_poles_retain_consistent_mutual_mass_and_material_section_spin() {
    with_cx(&CancelGate::new(), |cx| {
        let eta_bending = 0.16;
        let rayleigh = Rayleigh { a0: 0.25, a1: 0.01 };
        let mut frame = member(Geometry::Linear, FrameMassFormulation::ConsistentReference);
        frame.fix(1, [true, false, false]);
        frame
            .add_material_viscous_element(0, 1, 5, source(0.0, eta_bending, BAND))
            .unwrap();
        // Exact Hermite centerline plus section-spin integrals for the two
        // free tip coordinates [uy,theta], divided by area A; j=I/A.
        let j = RADIUS.powi(2) / 4.0;
        let m00 = RHO * (13.0 / 35.0 + 6.0 * j / 5.0);
        let m01 = RHO * (-11.0 / 210.0 - j / 10.0);
        let m11 = RHO * (1.0 / 105.0 + 2.0 * j / 15.0);
        let k00 = 12.0 * YOUNG * j;
        let k01 = -6.0 * YOUNG * j;
        let k11 = 4.0 * YOUNG * j;
        let det_m = m00 * m11 - m01 * m01;
        let det_k = k00 * k11 - k01 * k01;
        let mixed = k00 * m11 + k11 * m00 - 2.0 * k01 * m01;
        // det(K-lambda M)=0 is an independent scalar quadratic. Source
        // bending viscosity plus Rayleigh give C=a0 M+(a1+eta/E)K here.
        let high = (mixed + (mixed * mixed - 4.0 * det_m * det_k).sqrt()) / (2.0 * det_m);
        let low = det_k / (det_m * high);
        let mut expected = Vec::new();
        for lambda in [low, high] {
            let decay = 0.5 * (rayleigh.a0 + (rayleigh.a1 + eta_bending / YOUNG) * lambda);
            let oscillation = (lambda - decay * decay).sqrt();
            assert!(oscillation.is_finite() && oscillation > 0.0);
            expected.extend([
                C64::new(-decay, -oscillation),
                C64::new(-decay, oscillation),
            ]);
        }
        let spectrum = frame.resonance_poles(rayleigh, controls(), cx).unwrap();
        assert_eq!(spectrum.free_dofs, [4, 5]);
        assert_eq!(spectrum.material_omega_band_rad_s, Some(BAND));
        expect_rates(&spectrum.rates_per_s, &expected, 2e-10, 0.0);
    });
}

#[test]
fn free_point_mass_retains_rigid_position_and_velocity_relaxation() {
    with_cx(&CancelGate::new(), |cx| {
        let mut frame = Frame2d::new(Geometry::Linear);
        frame.add_node(0.0, 0.0);
        frame.add_mass(0, 2.0);
        frame.fix(0, [false, true, true]);
        for a0 in [0.0, 0.75] {
            let spectrum = frame
                .resonance_poles(Rayleigh { a0, a1: 0.0 }, controls(), cx)
                .unwrap();
            assert_eq!(spectrum.free_dofs, [0]);
            assert_eq!(spectrum.material_omega_band_rad_s, None);
            // 2*s*(s+a0)=0: no invented stiffness or removal of the zero rate.
            expect_rates(
                &spectrum.rates_per_s,
                &[C64::new(-a0, 0.0), C64::new(0.0, 0.0)],
                2e-14,
                0.0,
            );
        }
    });
}

#[test]
fn refusals_cover_full_native_budgets_nonlinearity_and_massless_rotation() {
    with_cx(&CancelGate::new(), |cx| {
        let frame = guided(2.0);
        let accepted = format!("{frame:?}");
        let mut cap = controls();
        // The scalar reduction still assembles six original coordinates.
        cap.maximum_coordinates = 1;
        assert!(matches!(
            frame.resonance_poles(NO_RAYLEIGH, cap, cx),
            Err(FramePoleError::Budget {
                requested: 6,
                cap: 1,
                ..
            })
        ));
        cap = controls();
        cap.maximum_payload_bytes = 0;
        assert!(matches!(
            frame.resonance_poles(NO_RAYLEIGH, cap, cx),
            Err(FramePoleError::Budget { .. })
        ));
        cap = controls();
        cap.maximum_work = 0;
        assert!(matches!(
            frame.resonance_poles(NO_RAYLEIGH, cap, cx),
            Err(FramePoleError::Budget { .. })
        ));
        cap = controls();
        cap.maximum_qr_sweeps_per_window = 0;
        assert!(matches!(
            frame.resonance_poles(NO_RAYLEIGH, cap, cx),
            Err(FramePoleError::Invalid(_))
        ));
        assert_eq!(format!("{frame:?}"), accepted);

        let elastic = || SectionModel::Elastic { ea: 8.0, ei: 2.0 };
        let mut nonlinear = member(
            Geometry::Corotational,
            FrameMassFormulation::ConsistentReference,
        );
        nonlinear.add_element(0, 1, 5, 6.0, &elastic).unwrap();
        let mut fiber = member(Geometry::Linear, FrameMassFormulation::ConsistentReference);
        fiber
            .add_element(0, 1, 5, 6.0, &|| {
                SectionModel::Fiber(fs_solid::fiber::rc_section(0.5, 0.5, 4, 1.5e-3))
            })
            .unwrap();
        let mut massless_rotation = member(Geometry::Linear, FrameMassFormulation::Lumped);
        massless_rotation
            .add_element(0, 1, 5, 6.0, &elastic)
            .unwrap();
        for unsupported in [nonlinear, fiber, massless_rotation] {
            assert!(matches!(
                unsupported.resonance_poles(NO_RAYLEIGH, controls(), cx),
                Err(FramePoleError::Invalid(_))
            ));
        }
    });
}

#[test]
fn cancelled_query_publishes_no_spectrum_or_source_change() {
    let frame = guided(2.0);
    let accepted = format!("{frame:?}");
    let gate = CancelGate::new();
    gate.request();
    with_cx(&gate, |cx| {
        assert!(matches!(
            frame.resonance_poles(NO_RAYLEIGH, controls(), cx),
            Err(FramePoleError::Cancelled)
        ));
    });
    assert_eq!(format!("{frame:?}"), accepted);
}
