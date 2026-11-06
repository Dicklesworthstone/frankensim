//! Independent cylindrical Ampere/energy fixtures, with synthetic source cards.
mod common;

use common::{config, limits, material};
use fs_em::planar::{MagneticError, axisymmetric::AxisymmetricMagnetostatics};
use fs_rep_mesh::{Metric2, TriComplex2};

fn rectangle(r0: f64, r1: f64, height: f64, nr: usize, nz: usize, angle: f64) -> TriComplex2 {
    let mut points = Vec::new();
    let mut faces = Vec::new();
    for j in 0..=nz {
        for i in 0..=nr {
            points.push([
                r0 + (r1 - r0) * i as f64 / nr as f64,
                height * j as f64 / nz as f64,
            ]);
        }
    }
    for j in 0..nz {
        for i in 0..nr {
            let a = (j * (nr + 1) + i) as u32;
            let b = a + (nr + 1) as u32;
            faces.push([a, a + 1, b + 1]);
            faces.push([a, b + 1, b]);
        }
    }
    TriComplex2::from_indexed_triangles(
        "axisymmetric magnetic fixture",
        points,
        faces,
        Metric2::axisymmetric(angle).unwrap(),
    )
    .unwrap()
}

fn uniform(mesh: TriComplex2, mu: f64) -> AxisymmetricMagnetostatics {
    let states = vec![material(mu); mesh.faces().len()];
    AxisymmetricMagnetostatics::new(mesh, states, limits(), &mut || false).unwrap()
}

fn outer(p: &AxisymmetricMagnetostatics, r: f64, value: f64) -> Vec<(usize, f64)> {
    p.mesh()
        .vertices()
        .iter()
        .enumerate()
        .filter(|(_, point)| point[0] == r)
        .map(|(v, _)| (v, value))
        .collect()
}

fn close(actual: f64, expected: f64, scale: f64) {
    assert!(
        (actual - expected).abs() <= 2e-8 * scale,
        "{actual} != {expected}; scale {scale}"
    );
}

#[test]
fn uniform_axial_field_is_exact_with_axis_regularization_and_natural_end_faces() {
    let mu = 2e-6;
    let b0 = 0.7;
    let radius = 0.03;
    let height = 0.08;
    let angle = core::f64::consts::TAU;
    let p = uniform(rectangle(0.0, radius, height, 8, 5, angle), mu);
    let f = p
        .solve(
            &vec![0.0; p.mesh().faces().len()],
            &outer(&p, radius, b0 / 2.0),
            config(),
            &mut || false,
        )
        .unwrap();
    for &u in f.potential_over_radius_t() {
        close(u, b0 / 2.0, b0);
    }
    for cell in 0..p.mesh().faces().len() {
        for i in 0..3 {
            let mut barycentric = [0.0; 3];
            barycentric[i] = 1.0;
            let point = p.mesh().vertices()[p.mesh().faces()[cell][i] as usize];
            let sample = f.sample(cell, barycentric).unwrap();
            close(sample.flux_density_t[0], 0.0, b0);
            close(sample.flux_density_t[1], b0, b0);
            close(sample.magnetic_field_a_m[1], b0 / mu, b0 / mu);
            close(
                sample.azimuthal_potential_wb_m,
                b0 * point[0] / 2.0,
                b0 * radius,
            );
            if point[0] == 0.0 {
                assert_eq!(sample.azimuthal_potential_wb_m, 0.0);
                assert_eq!(sample.flux_density_t[0], 0.0);
            }
        }
    }
    let energy = angle * radius * radius * height * b0 * b0 / (4.0 * mu);
    close(f.energy_j(), energy, energy);
    close(f.boundary_work_j(), 2.0 * energy, energy);
    close(f.impressed_current_work_j(), 0.0, energy);
    close(f.work_balance_error_j(), 0.0, energy);
    assert!(f.solve_report().unwrap().converged_euclidean());
    assert_eq!(p.materials().len(), p.mesh().faces().len());
}

#[test]
fn source_free_quadrupole_recovers_both_cylindrical_components() {
    // Maxwell quadrupole: Br=-beta*r, Bz=2c+2beta*z has div(B)=0 and
    // curl(H)=0. The exact regularized potential is the affine u=c+beta*z.
    let radius = 0.03_f64;
    let height = 0.08_f64;
    let c = 0.3;
    let beta = 2.0;
    let mu = 1.4e-6;
    let angle = core::f64::consts::TAU;
    let p = uniform(rectangle(0.0, radius, height, 8, 5, angle), mu);
    let boundary: Vec<_> = p
        .mesh()
        .vertices()
        .iter()
        .enumerate()
        .filter(|(_, p)| p[0] == radius || p[1] == 0.0 || p[1] == height)
        .map(|(v, p)| (v, c + beta * p[1]))
        .collect();
    let f = p
        .solve(
            &vec![0.0; p.mesh().faces().len()],
            &boundary,
            config(),
            &mut || false,
        )
        .unwrap();
    for (cell, face) in p.mesh().faces().iter().enumerate() {
        let r = face
            .iter()
            .map(|&v| p.mesh().vertices()[v as usize][0])
            .sum::<f64>()
            / 3.0;
        let z = face
            .iter()
            .map(|&v| p.mesh().vertices()[v as usize][1])
            .sum::<f64>()
            / 3.0;
        let s = f.sample(cell, [1.0 / 3.0; 3]).unwrap();
        close(s.flux_density_t[0], -beta * r, c);
        close(s.flux_density_t[1], 2.0 * (c + beta * z), c);
        close(s.azimuthal_potential_wb_m, r * (c + beta * z), radius * c);
    }
    let energy = angle / (2.0 * mu)
        * (beta * beta * radius.powi(4) * height / 4.0
            + 2.0 * radius * radius
                * (c * c * height
                    + c * beta * height * height
                    + beta * beta * height.powi(3) / 3.0));
    close(f.energy_j(), energy, energy);
    close(f.work_balance_error_j(), 0.0, energy);
}

#[test]
fn azimuthal_current_has_exact_cylindrical_curl_and_current_work() {
    // u=alpha(R-r): Bz=alpha(2R-3r), hence Jphi=-dHz/dr=3alpha/mu.
    // With u(R)=0, no prescribed boundary work is supplied. This finite-domain
    // flux closure is explicit; it is not the free-space field of a solenoid.
    let radius = 0.04;
    let height = 0.07;
    let alpha = 2.0;
    let mu = 1.5e-6;
    let j = 3.0 * alpha / mu;
    let angle = core::f64::consts::TAU;
    let p = uniform(rectangle(0.0, radius, height, 8, 4, angle), mu);
    let f = p
        .solve(
            &vec![j; p.mesh().faces().len()],
            &outer(&p, radius, 0.0),
            config(),
            &mut || false,
        )
        .unwrap();
    for (&[r, _], &u) in p.mesh().vertices().iter().zip(f.potential_over_radius_t()) {
        close(u, alpha * (radius - r), alpha * radius);
    }
    for (cell, face) in p.mesh().faces().iter().enumerate() {
        let r = face
            .iter()
            .map(|&v| p.mesh().vertices()[v as usize][0])
            .sum::<f64>()
            / 3.0;
        let s = f.sample(cell, [1.0 / 3.0; 3]).unwrap();
        close(s.flux_density_t[0], 0.0, alpha * radius);
        close(
            s.flux_density_t[1],
            alpha * (2.0 * radius - 3.0 * r),
            alpha * radius,
        );
    }
    let energy = angle * height * alpha * alpha * radius.powi(4) / (8.0 * mu);
    close(f.energy_j(), energy, energy);
    close(f.impressed_current_work_j(), 2.0 * energy, energy);
    close(f.boundary_work_j(), 0.0, energy);
    close(f.work_balance_error_j(), 0.0, energy);
    // Angle changes only extensive energy and work, not local fields.
    let sector = uniform(rectangle(0.0, radius, height, 8, 4, angle / 7.0), mu);
    let sf = sector
        .solve(
            &vec![j; sector.mesh().faces().len()],
            &outer(&sector, radius, 0.0),
            config(),
            &mut || false,
        )
        .unwrap();
    for (&u, &v) in f
        .potential_over_radius_t()
        .iter()
        .zip(sf.potential_over_radius_t())
    {
        close(u, v, alpha * radius);
    }
    close(sf.energy_j() * 7.0, f.energy_j(), energy);
    close(
        sf.impressed_current_work_j() * 7.0,
        f.impressed_current_work_j(),
        energy,
    );
}

#[test]
fn thick_solenoid_with_permeable_core_converges_to_independent_ampere_solution() {
    // Infinite-length solenoid slice: Jphi is constant in a<r<b. On natural
    // axial ends Br=0. H_z=J(b-a) in the core, J(b-r) in the winding, zero
    // outside. The exterior Aphi=F/r is fixed at the finite outer radius.
    let a = 0.02_f64;
    let b = 0.04_f64;
    let radius = 0.06_f64;
    let height = 0.08;
    let mu = 4.0e-7 * core::f64::consts::PI;
    let mu_core = 7.0 * mu;
    let j = 2e6;
    let angle = core::f64::consts::TAU;
    let flux_integral = mu_core * j * (b - a) * a * a / 2.0
        + mu * j * (b * (b * b - a * a) / 2.0 - (b.powi(3) - a.powi(3)) / 3.0);
    let integral = |r: f64| b * b * r * r / 2.0 - 2.0 * b * r.powi(3) / 3.0 + r.powi(4) / 4.0;
    let exact_energy = angle * height * j * j / 2.0
        * (mu_core * (b - a).powi(2) * a * a / 2.0 + mu * (integral(b) - integral(a)));
    let mut previous = f64::INFINITY;
    let mut previous_field_error = f64::INFINITY;
    for nr in [12, 24, 48, 96] {
        let mesh = rectangle(0.0, radius, height, nr, 4, angle);
        let mut currents = Vec::new();
        let mut states = Vec::new();
        for face in mesh.faces() {
            let r = face
                .iter()
                .map(|&v| mesh.vertices()[v as usize][0])
                .sum::<f64>()
                / 3.0;
            states.push(material(if r < a { mu_core } else { mu }));
            currents.push(if r > a && r < b { j } else { 0.0 });
        }
        let p = AxisymmetricMagnetostatics::new(mesh, states, limits(), &mut || false).unwrap();
        let f = p
            .solve(
                &currents,
                &outer(&p, radius, flux_integral / radius.powi(2)),
                config(),
                &mut || false,
            )
            .unwrap();
        let error = (f.energy_j() - exact_energy).abs() / exact_energy;
        assert!(
            error < 0.4 * previous,
            "nr={nr}: energy relative error {error}, previous {previous}"
        );
        previous = error;
        // Independent positive degree-four triangle quadrature integrates
        // the exact energy norm of the affine field error, including r dA.
        // Resolving the permeability interface gives first-order field and
        // second-order stored-energy convergence for this P1 formulation.
        let mut field_error_energy = 0.0;
        for (cell, face) in p.mesh().faces().iter().enumerate() {
            let radii = face.map(|v| p.mesh().vertices()[v as usize][0]);
            let mean_radius = radii.iter().sum::<f64>() / 3.0;
            let volume = p.mesh().face_measure(cell).unwrap();
            let permeability = p.permeability_h_m()[cell];
            for (coordinate, weight) in [
                (0.445948490915965, 0.223381589678011),
                (0.091576213509771, 0.109951743655322),
            ] {
                for barycentric in [
                    [coordinate, coordinate, 1.0 - 2.0 * coordinate],
                    [coordinate, 1.0 - 2.0 * coordinate, coordinate],
                    [1.0 - 2.0 * coordinate, coordinate, coordinate],
                ] {
                    let r = radii
                        .iter()
                        .zip(barycentric)
                        .map(|(&radius, shape)| radius * shape)
                        .sum::<f64>();
                    let expected_h = if r < a {
                        j * (b - a)
                    } else if r < b {
                        j * (b - r)
                    } else {
                        0.0
                    };
                    let sample = f.sample(cell, barycentric).unwrap();
                    let radial_error = sample.flux_density_t[0];
                    let axial_error = sample.flux_density_t[1] - permeability * expected_h;
                    field_error_energy += weight * volume * (r / mean_radius)
                        * (radial_error * radial_error + axial_error * axial_error)
                        / (2.0 * permeability);
                }
            }
        }
        let field_error = (field_error_energy / exact_energy).sqrt();
        assert!(
            field_error < 0.55 * previous_field_error,
            "nr={nr}: field energy-norm error {field_error}, previous {previous_field_error}"
        );
        previous_field_error = field_error;
        close(f.work_balance_error_j(), 0.0, exact_energy);
    }
    assert!(
        previous < 0.0003,
        "energy did not converge sufficiently: {previous}"
    );
    assert!(
        previous_field_error < 0.04,
        "field did not converge sufficiently: {previous_field_error}"
    );
}

#[test]
fn boundaries_budgets_sources_and_cancellation_are_admitted_before_publication() {
    let p = uniform(rectangle(0.0, 1.0, 1.0, 3, 3, 1.0), 1.0);
    let current = vec![0.0; p.mesh().faces().len()];
    // u on the axis is not A_phi; prescribing only the axis must not anchor it.
    let axis: Vec<_> = p
        .mesh()
        .vertices()
        .iter()
        .enumerate()
        .filter(|(_, p)| p[0] == 0.0)
        .map(|(v, _)| (v, 0.0))
        .collect();
    assert!(p.solve(&current, &axis, config(), &mut || false).is_err());
    assert!(p.solve(&current, &[], config(), &mut || false).is_err());
    assert!(
        p.solve(&current, &[(3, 0.0)], config(), &mut || false)
            .is_err()
    );
    assert!(
        p.solve(
            &current,
            &[(3, 0.0), (7, 0.0), (7, 0.0)],
            config(),
            &mut || false
        )
        .is_err()
    );
    assert!(
        p.solve(&[0.0], &outer(&p, 1.0, 1.0), config(), &mut || false)
            .is_err()
    );
    let mut invalid = current.clone();
    invalid[0] = f64::NAN;
    assert!(
        p.solve(&invalid, &outer(&p, 1.0, 1.0), config(), &mut || false)
            .is_err()
    );
    assert!(matches!(
        p.solve(&current, &outer(&p, 1.0, 1.0), config(), &mut || true),
        Err(MagneticError::Cancelled)
    ));
    let mesh = rectangle(0.0, 1.0, 1.0, 3, 3, 1.0);
    let count = mesh.faces().len();
    let mut budget = limits();
    budget.max_triplets = 9 * count - 1;
    assert!(
        AxisymmetricMagnetostatics::new(mesh, vec![material(1.0); count], budget, &mut || false)
            .is_err()
    );
    let mesh = rectangle(0.0, 1.0, 1.0, 3, 3, 1.0);
    let mut polls = 0;
    assert!(matches!(
        AxisymmetricMagnetostatics::new(mesh, vec![material(1.0); count], limits(), &mut || {
            polls += 1;
            polls > count + 10
        }),
        Err(MagneticError::Cancelled)
    ));
    let f = p
        .solve(&current, &outer(&p, 1.0, 1.0), config(), &mut || false)
        .unwrap();
    assert!(f.sample(0, [0.5, 0.5, 0.5]).is_err());
    assert!(f.sample(0, [-0.1, 0.6, 0.5]).is_err());
    assert!(f.sample(count, [1.0, 0.0, 0.0]).is_err());
}
