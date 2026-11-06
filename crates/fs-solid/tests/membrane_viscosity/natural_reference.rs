//! Shifted natural references use the original hyperelastic volume convention.
use super::*;
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter};
use fs_material::hyper::plane_stress::PlaneStressSettings;

const THICKNESS_M: f64 = 0.002;
const VALUES: [f64; 6] = [4.0, 3.0, 2.0, 1.0, 0.4, 0.2];

fn source(values: [f64; 6]) -> (ResolvedSheet, PlaneStressViscosity) {
    let provenance = Provenance {
        source: "synthetic natural-reference viscosity; no specimen authority".into(),
        license: "test".into(),
        artifact: None,
    };
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    let mut properties = vec![
        (AREAL_DENSITY, 6.0, Dims([-2, 1, 0, 0, 0, 0])),
        (THICKNESS, THICKNESS_M, Dims([1, 0, 0, 0, 0, 0])),
    ];
    properties.extend(
        KEYS.into_iter()
            .zip(values)
            .map(|(key, value)| (key, value, DynViscosity::DIMS)),
    );
    for (key, value, dims) in properties {
        // Reference inventory is static. Only the viscosity claims carry
        // frequency applicability; it cannot constrain hyperelastic mass.
        let mut validity = ValidityDomain::unconstrained().with("T", 290.0, 310.0);
        if KEYS.contains(&key) {
            validity = validity.with("omega", 0.0, 1000.0);
        }
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(key, dims),
                value: PropertyValue::Scalar { value, dims },
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                validity,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: provenance.clone(),
            })
            .unwrap();
        requirements.push(
            ScalarPropertyRequirement::try_new(key, dims, ScalarAdmissibility::Finite).unwrap(),
        );
    }
    let pa = Dims([-1, 1, -2, 0, 0, 0]);
    let model = ConstitutiveModelCard {
        law: LawId("neo-hookean-log".into()),
        law_version: 1,
        parameters: [
            ("mu", 3.0e5, pa),
            ("lambda", 8.0e5, pa),
            ("stretch_limit", 5.0, Dims::NONE),
        ]
        .into_iter()
        .map(|(name, value, dims)| (name.into(), LawParameter { value, dims }))
        .collect(),
        state_schema_version: 1,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with("T", 290.0, 310.0),
        sources: vec![],
        provenance,
    };
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic hyperelastic sheet".into(),
            phase: "solid".into(),
            process: "natural-reference damping verification".into(),
            revision: 0,
        },
        claims,
        vec![model],
    )
    .unwrap();
    let material = resolve_material_state_point(
        &card,
        &QueryPoint::new()
            .with("T", 300.0)
            .unwrap()
            .with("omega", 0.0)
            .unwrap(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let law =
        PlaneStressViscosity::from_material(material.clone(), MaterialPlane::Axes12, KEYS).unwrap();
    let sheet = ResolvedSheet::try_from_material(material)
        .unwrap()
        .with_hyperelastic_membrane(
            &card,
            0,
            PlaneStressSettings {
                relative_tolerance: 2e-13,
                ..PlaneStressSettings::default()
            },
        )
        .unwrap();
    (sheet, law)
}

fn from_sheet(sheet: ResolvedSheet, grain: f64) -> MembranePatch {
    MembranePatch::try_new(
        Mesh2 {
            nodes: vec![[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]],
            elems: vec![vec![0, 1, 2]],
            patches: vec![],
        },
        vec![sheet],
        vec![grain],
        MembraneDomain {
            minimum_area_ratio: 0.01,
            maximum_abs_strain: 0.5,
        },
        MembraneBudget {
            max_nodes: 3,
            max_triangles: 1,
        },
        || false,
    )
    .unwrap()
}

// Independent oracle: form Fe and Fedot from the current physical edges
// divided by their natural reference lengths. Integrate the source quadratic
// over actual natural area and thickness, without using a pulled viscosity
// tensor, the membrane's strain maps or its assembled force.
fn natural_power(
    u: &[[f64; 3]],
    velocity: &[[f64; 3]],
    stretches: [f64; 3],
    grain: f64,
    viscosity: &[f64; 9],
) -> f64 {
    let edges = [
        [
            1.0 + u[1][0] - u[0][0],
            u[1][1] - u[0][1],
            u[1][2] - u[0][2],
        ],
        [
            u[2][0] - u[0][0],
            1.0 + u[2][1] - u[0][1],
            u[2][2] - u[0][2],
        ],
    ];
    let edge_rates: [[f64; 3]; 2] =
        core::array::from_fn(|a| core::array::from_fn(|i| velocity[a + 1][i] - velocity[0][i]));
    let (s, c) = grain.sin_cos();
    let natural_columns = |edges: [[f64; 3]; 2]| -> [[f64; 3]; 2] {
        [
            core::array::from_fn(|i| (c * edges[0][i] + s * edges[1][i]) / stretches[0]),
            core::array::from_fn(|i| (-s * edges[0][i] + c * edges[1][i]) / stretches[1]),
        ]
    };
    let f = natural_columns(edges);
    let fdot = natural_columns(edge_rates);
    let dot = |a: [f64; 3], b: [f64; 3]| -> f64 { a.into_iter().zip(b).map(|(a, b)| a * b).sum() };
    let rate = [
        dot(f[0], fdot[0]),
        dot(f[1], fdot[1]),
        dot(f[0], fdot[1]) + dot(f[1], fdot[0]),
    ];
    let density: f64 = (0..3)
        .map(|i| {
            (0..3)
                .map(|j| rate[i] * viscosity[3 * i + j] * rate[j])
                .sum::<f64>()
        })
        .sum();
    let natural_area = 0.5 * stretches[0] * stretches[1];
    let natural_thickness = THICKNESS_M * stretches[2];
    natural_area * natural_thickness * density
}

#[test]
fn natural_volume_heat_and_rayleigh_force_agree_in_rotated_material_axes() {
    let stretches = [0.8, 0.9, 1.1];
    let (sheet, law) = source(VALUES);
    let matrix = *law.matrix_pa_s();
    let u = deformed();
    let velocity = [[0.2, -0.1, 0.3], [0.32, -0.14, 0.38], [0.18, 0.09, 0.47]];
    for grain in [0.0, 0.37] {
        let p = from_sheet(
            sheet
                .clone()
                .with_hyperelastic_reference_stretches(stretches)
                .unwrap(),
            grain,
        );
        let damping =
            MembraneViscosity::new(&p, vec![law.clone()], (0.0, 1000.0), || false).unwrap();
        assert_eq!(damping.laws()[0].matrix_pa_s(), &matrix);
        assert_eq!(damping.laws()[0].source(), p.sheets()[0].material());
        let response = damping.evaluate(&u, &velocity, || false).unwrap();
        let expected = natural_power(&u, &velocity, stretches, grain, &matrix);
        close(response.power_w, expected);
        close(response.element_power_w[0], expected);
        let work: f64 = response
            .resisting_force_n
            .iter()
            .flatten()
            .zip(velocity.iter().flatten())
            .map(|(f, v)| f * v)
            .sum();
        close(work, expected);
        // f = d(P/2)/dv checks every nodal force component independently.
        let step = 1e-4;
        for node in 0..3 {
            for axis in 0..3 {
                let mut plus = velocity;
                let mut minus = velocity;
                plus[node][axis] += step;
                minus[node][axis] -= step;
                let derivative = (natural_power(&u, &plus, stretches, grain, &matrix)
                    - natural_power(&u, &minus, stretches, grain, &matrix))
                    / (4.0 * step);
                close(response.resisting_force_n[node][axis], derivative);
            }
        }
        let unpulled = law.stress_and_power([0.14, 0.2219, 0.0111]).unwrap().1 * 0.5 * THICKNESS_M;
        if grain == 0.0 {
            assert!((expected - unpulled).abs() > 0.1 * expected);
        }
    }
}

#[test]
fn normal_reference_length_scales_loss_like_elastic_energy_without_changing_mass() {
    let (sheet, law) = source(VALUES);
    let make = |normal| {
        from_sheet(
            sheet
                .clone()
                .with_hyperelastic_reference_stretches([0.8, 0.9, normal])
                .unwrap(),
            0.0,
        )
    };
    let first = make(0.7);
    let twice = make(1.4);
    assert_eq!(first.reference_mass_kg(), twice.reference_mass_kg());
    assert_eq!(first.sheets()[0].thickness_m(), THICKNESS_M);
    assert_eq!(twice.sheets()[0].thickness_m(), THICKNESS_M);
    let u = deformed();
    close(
        twice.evaluate(&u, || false).unwrap().strain_energy_j,
        2.0 * first.evaluate(&u, || false).unwrap().strain_energy_j,
    );
    let velocity = [[0.0; 3], [0.0, 0.0, 1.3], [0.0, 0.0, -0.7]];
    let a = MembraneViscosity::new(&first, vec![law.clone()], (0.0, 1000.0), || false)
        .unwrap()
        .evaluate(&u, &velocity, || false)
        .unwrap();
    let b = MembraneViscosity::new(&twice, vec![law], (0.0, 1000.0), || false)
        .unwrap()
        .evaluate(&u, &velocity, || false)
        .unwrap();
    close(b.power_w, 2.0 * a.power_w);
    for (a, b) in a
        .resisting_force_n
        .iter()
        .flatten()
        .zip(b.resisting_force_n.iter().flatten())
    {
        close(*b, 2.0 * a);
    }
}

#[test]
fn identity_reference_preserves_original_force_heat_and_modal_arithmetic() {
    let (sheet, law) = source(VALUES);
    let original = from_sheet(sheet.clone(), 0.37);
    let identity = from_sheet(
        sheet
            .with_hyperelastic_reference_stretches([1.0; 3])
            .unwrap(),
        0.37,
    );
    let a = MembraneViscosity::new(&original, vec![law.clone()], (0.0, 1000.0), || false).unwrap();
    let b = MembraneViscosity::new(&identity, vec![law], (0.0, 1000.0), || false).unwrap();
    let u = deformed();
    let velocity = [[0.0; 3], [0.0, 0.0, 1.3], [0.0, 0.0, -0.7]];
    let left = a.evaluate(&u, &velocity, || false).unwrap();
    let right = b.evaluate(&u, &velocity, || false).unwrap();
    assert_eq!(left.resisting_force_n, right.resisting_force_n);
    assert_eq!(left.power_w, right.power_w);
    assert_eq!(left.element_power_w, right.element_power_w);
    let left_modes = modes(&original, &u, &[5, 8]);
    let right_modes = modes(&identity, &u, &[5, 8]);
    let left = a.project_modes(&left_modes, 100_000, || false).unwrap();
    let right = b.project_modes(&right_modes, 100_000, || false).unwrap();
    assert_eq!(left.matrix_per_s(), right.matrix_per_s());
    assert_eq!(left.dissipation_rows(), right.dissipation_rows());
}

#[test]
fn natural_reference_modal_gram_preserves_full_source_power() {
    let stretches = [0.8, 0.9, 1.1];
    let grain = 0.37;
    let (sheet, law) = source(VALUES);
    let matrix = *law.matrix_pa_s();
    let p = from_sheet(
        sheet
            .with_hyperelastic_reference_stretches(stretches)
            .unwrap(),
        grain,
    );
    let damping = MembraneViscosity::new(&p, vec![law], (0.0, 1000.0), || false).unwrap();
    let u = deformed();
    let basis = modes(&p, &u, &[5, 8]);
    let projected = damping.project_modes(&basis, 100_000, || false).unwrap();
    let c = projected.matrix_per_s();
    assert_eq!(c[1], c[2]);
    assert!(c[1].abs() > 1e-8);
    for qdot in [[1.0, 0.0], [0.0, 1.0], [0.7, -0.3], [0.7, 0.3]] {
        let velocity = basis.reconstruct(&qdot).unwrap();
        let response = damping.evaluate(&u, &velocity, || false).unwrap();
        let expected = natural_power(&u, &velocity, stretches, grain, &matrix);
        let modal = qdot[0] * (c[0] * qdot[0] + c[1] * qdot[1])
            + qdot[1] * (c[2] * qdot[0] + c[3] * qdot[1]);
        let rows: f64 = projected
            .dissipation_rows()
            .chunks_exact(2)
            .map(|row| (row[0] * qdot[0] + row[1] * qdot[1]).powi(2))
            .sum();
        close(response.power_w, expected);
        close(modal, expected);
        close(rows, expected);
    }
}

#[test]
fn natural_reference_keeps_rigid_motion_flat_transverse_and_directional_nulls() {
    let (sheet, law) = source(VALUES);
    let p = from_sheet(
        sheet
            .with_hyperelastic_reference_stretches([0.8, 0.9, 1.1])
            .unwrap(),
        0.37,
    );
    let damping = MembraneViscosity::new(&p, vec![law], (0.0, 1000.0), || false).unwrap();
    let u = deformed();
    let rigid_translation = damping
        .evaluate(&u, &[[0.7, -0.4, 0.2]; 3], || false)
        .unwrap();
    assert_eq!(rigid_translation.power_w, 0.0);
    assert!(
        rigid_translation
            .resisting_force_n
            .iter()
            .flatten()
            .all(|v| *v == 0.0)
    );
    let positions: Vec<[f64; 3]> = p
        .mesh()
        .nodes
        .iter()
        .zip(&u)
        .map(|(x, u)| [x[0] + u[0], x[1] + u[1], u[2]])
        .collect();
    let spin: Vec<_> = positions
        .iter()
        .map(|x| {
            [
                -1.1 * x[2] - 0.4 * x[1],
                0.4 * x[0] - 0.7 * x[2],
                0.7 * x[1] + 1.1 * x[0],
            ]
        })
        .collect();
    let rigid_spin = damping.evaluate(&u, &spin, || false).unwrap();
    assert!(rigid_spin.power_w < 1e-28);
    assert!(
        rigid_spin
            .resisting_force_n
            .iter()
            .flatten()
            .all(|f| f.abs() < 1e-14)
    );
    let velocity = [[0.2, -0.3, 0.1], [-0.4, 0.1, 0.8], [0.3, 0.6, -0.2]];
    let before = damping.evaluate(&u, &velocity, || false).unwrap();
    let rotate = |x: [f64; 3]| [x[2], x[0], x[1]];
    let rotated_u: Vec<_> = positions
        .iter()
        .zip(&p.mesh().nodes)
        .map(|(x, reference)| {
            let x = rotate(*x);
            [
                x[0] + 0.2 - reference[0],
                x[1] - 0.3 - reference[1],
                x[2] + 0.7,
            ]
        })
        .collect();
    let after = damping
        .evaluate(&rotated_u, &velocity.map(rotate), || false)
        .unwrap();
    close(after.power_w, before.power_w);
    for (actual, original) in after
        .resisting_force_n
        .iter()
        .zip(&before.resisting_force_n)
    {
        for (a, b) in actual.iter().zip(rotate(*original)) {
            close(*a, b);
        }
    }

    let flat = [[0.0; 3]; 3];
    let transverse = [[0.0; 3], [0.0, 0.0, 1.0], [0.0; 3]];
    assert!(p.tangent_apply(&flat, &transverse, || false).unwrap()[1][2] > 0.0);
    let flat_loss = damping.evaluate(&flat, &transverse, || false).unwrap();
    assert_eq!(flat_loss.power_w, 0.0);
    assert!(
        flat_loss
            .resisting_force_n
            .iter()
            .flatten()
            .all(|f| *f == 0.0)
    );
    let flat_modes = modes(&p, &flat, &[5]);
    assert_eq!(
        damping
            .project_modes(&flat_modes, 100_000, || false)
            .unwrap()
            .matrix_per_s(),
        &[0.0]
    );

    // Rank-one V = [1,2,0]^T [1,2,0]; natural rate [2,-1,0] is null.
    let (sheet, law) = source([1.0, 4.0, 0.0, 2.0, 0.0, 0.0]);
    let p = from_sheet(
        sheet
            .with_hyperelastic_reference_stretches([0.5, 1.0, 1.5])
            .unwrap(),
        0.0,
    );
    let damping = MembraneViscosity::new(&p, vec![law], (0.0, 1000.0), || false).unwrap();
    let null = [[0.0; 3], [0.5, 0.0, 0.0], [0.0, -1.0, 0.0]];
    let loss = damping.evaluate(&flat, &null, || false).unwrap();
    assert!(loss.power_w < 1e-28);
    assert!(
        loss.resisting_force_n
            .iter()
            .flatten()
            .all(|f| f.abs() < 1e-14)
    );
    let positive = [[0.0; 3], [1.0, 0.0, 0.0], [0.0; 3]];
    close(
        damping
            .evaluate(&flat, &positive, || false)
            .unwrap()
            .power_w,
        0.012,
    );
}
