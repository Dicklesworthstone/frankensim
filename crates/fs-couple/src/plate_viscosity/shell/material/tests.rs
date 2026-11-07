use super::*;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::plane_stress_viscosity::{
    MaterialPlane, PLANE_STRESS_VISCOSITY_PROPERTIES as KEYS,
};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_material_state_point,
};
use fs_plate::PlateSection;
use fs_plate::shell::{CstDktMass, ShellSupport, assemble_shell_with_mass};
use fs_qty::DynViscosity;

pub(super) fn source(values: [f64; 6], band: (f64, f64)) -> Arc<PlaneStressViscosity> {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (key, value) in KEYS.into_iter().zip(values) {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(key, DynViscosity::DIMS),
                value: PropertyValue::Scalar {
                    value,
                    dims: DynViscosity::DIMS,
                },
                validity: ValidityDomain::unconstrained()
                    .with("T", 290.0, 300.0)
                    .with("omega", band.0, band.1),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic shell viscosity".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
        requirements.push(
            ScalarPropertyRequirement::try_new(
                key,
                DynViscosity::DIMS,
                ScalarAdmissibility::Finite,
            )
            .unwrap(),
        );
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic".into(),
            phase: "solid".into(),
            process: "shell damping test".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let state = resolve_material_state_point(
        &card,
        &QueryPoint::new()
            .with("T", 295.0)
            .unwrap()
            .with("omega", band.0)
            .unwrap(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    Arc::new(PlaneStressViscosity::from_material(state, MaterialPlane::Axes12, KEYS).unwrap())
}

pub(super) fn fixture() -> (ShellMesh, ShellModel, Vec<ModePair>) {
    let mesh = ShellMesh::new(
        vec![[0.0, 0.0, 0.0], [0.2, 0.0, 0.0], [0.0, 0.1, 0.0]],
        vec![[0, 1, 2]],
    )
    .unwrap();
    let section = PlateSection::isotropic(1e6, 0.25, 0.005, 1000.0).unwrap();
    let model = assemble_shell_with_mass(
        &mesh,
        &section,
        &[0, 1],
        ShellSupport::Clamped,
        CstDktMass::Consistent,
    )
    .unwrap();
    let n = model.free;
    let mut k = vec![0.0; n * n];
    let mut m = k.clone();
    for i in 0..n {
        for j in 0..n {
            k[n * i + j] = model.k.get(i, j);
            m[n * i + j] = model.m.get(i, j);
        }
    }
    let modes = fs_modal::eigh_gen_dense(&k, &m, n).unwrap();
    (mesh, model, modes)
}

pub(super) fn limits() -> PlateViscosityLimits {
    PlateViscosityLimits {
        max_modes: 12,
        max_projection_terms: 2_000_000,
        reference_resistance_n_s_m: 1.0,
    }
}

fn quadratic(a: &[f64; 9], x: [f64; 3]) -> f64 {
    (0..3)
        .map(|i| (0..3).map(|j| x[i] * a[3 * i + j] * x[j]).sum::<f64>())
        .sum()
}

#[test]
fn source_thickness_and_orientation_control_independent_membrane_bending_power() {
    let law = source([4.0, 9.0, 16.0, 1.0, 0.5, -0.75], (0.0, 1e6));
    let h = 0.004;
    let angle = 0.37;
    let facet = ShellMaterialViscosity::new(law.clone(), h, angle).unwrap();
    let (a, b, d) = facet.section().matrices();
    assert_eq!(b, [0.0; 9]);
    assert_eq!(facet.source(), law.as_ref());
    assert_eq!(facet.thickness_m(), h);
    assert_eq!(facet.material_angle_rad(), angle);
    // Independently rotate a physical engineering strain into material axes.
    let (s, c) = angle.sin_cos();
    let local = |e: [f64; 3]| {
        [
            c * c * e[0] + s * s * e[1] + c * s * e[2],
            s * s * e[0] + c * c * e[1] - c * s * e[2],
            -2.0 * c * s * e[0] + 2.0 * c * s * e[1] + (c * c - s * s) * e[2],
        ]
    };
    let membrane = [0.2, -0.1, 0.07];
    let curvature = [0.3, 0.4, -0.2];
    let expected = h * law.stress_and_power(local(membrane)).unwrap().1
        + h * h * h / 12.0 * law.stress_and_power(local(curvature)).unwrap().1;
    let actual = quadratic(&a, membrane) + quadratic(&d, curvature);
    assert!((actual - expected).abs() < 3e-12 * expected);
    let thick = ShellMaterialViscosity::new(law, 2.0 * h, angle).unwrap();
    let (a2, _, d2) = thick.section().matrices();
    for i in 0..9 {
        assert_eq!(a2[i], 2.0 * a[i]);
        assert_eq!(d2[i], 8.0 * d[i]);
    }
}

#[test]
fn source_modal_shell_uses_original_full_projection_and_retains_material_records() {
    let (mesh, model, modes) = fixture();
    let law = source([400.0, 900.0, 1600.0, 100.0, 50.0, -75.0], (0.0, 1e6));
    let facet = ShellMaterialViscosity::new(law.clone(), 0.005, 0.31).unwrap();
    let gate = CancelGate::new_clock_free();
    let bound = ModalShellMaterialViscosity::project(
        &mesh,
        &[facet.clone()],
        &model,
        &modes,
        limits(),
        &gate,
    )
    .unwrap();
    let direct = ModalPlateViscosity::from_cst_dkt_shell(
        &mesh,
        &[facet.section().clone()],
        &model,
        &modes.iter().map(|m| m.phi.clone()).collect::<Vec<_>>(),
        limits(),
        &gate,
    )
    .unwrap();
    assert_eq!(bound.modal().matrix_per_s(), direct.matrix_per_s());
    assert_eq!(bound.omega_band_rad_s(), (0.0, 1e6));
    assert_eq!(bound.facets()[0].source(), law.as_ref());
    assert!(bound.modal().channel_count() > 0);
    for (omega, mode) in bound.angular_frequencies_rad_s().iter().zip(modes) {
        assert_eq!(*omega, fs_math::det::sqrt(mode.lambda));
    }
}

#[test]
fn source_support_applies_to_actual_modes_and_exact_zero_material() {
    let (mesh, model, modes) = fixture();
    let gate = CancelGate::new_clock_free();
    for values in [[0.0; 6], [4.0, 9.0, 16.0, 0.0, 0.0, 0.0]] {
        let law = source(values, (0.0, 1e-3));
        let facet = ShellMaterialViscosity::new(law, 0.005, 0.0).unwrap();
        assert!(
            ModalShellMaterialViscosity::project(&mesh, &[facet], &model, &modes, limits(), &gate)
                .is_err()
        );
    }
    let law = source([0.0; 6], (0.0, 1e6));
    let facet = ShellMaterialViscosity::new(law.clone(), 0.005, 0.0).unwrap();
    let bound = ModalShellMaterialViscosity::project(
        &mesh,
        &[facet.clone()],
        &model,
        &modes,
        limits(),
        &gate,
    )
    .unwrap();
    assert_eq!(bound.modal().channel_count(), 0);
    assert!(ShellMaterialViscosity::new(law.clone(), 0.0, 0.0).is_err());
    assert!(ShellMaterialViscosity::new(law, 0.005, f64::NAN).is_err());
    let exhausted = PlateViscosityLimits {
        max_projection_terms: 1,
        ..limits()
    };
    assert!(
        ModalShellMaterialViscosity::project(&mesh, &[facet], &model, &modes, exhausted, &gate)
            .is_err()
    );
}
