#![cfg(feature = "certified-speculation")]
//! Actual verifier regressions; analytical slab laws do not depend on a FEM solver.
use std::collections::BTreeMap;
use fs_verify::interval::Iv;
use fs_verify::tet::{AffineSourceTetProblem, BoundaryCondition as Bc, BoundaryFace,
    FluxBudget, TetError, TetProblem, affine_source_mean_bound, energy_bound, mean_bound};

struct Slabs {
    vertices: Vec<[f64; 3]>, tets: Vec<[usize; 4]>, k: Vec<f64>, source: Vec<f64>,
    boundary: Vec<BoundaryFace>, per: usize,
}
impl Slabs {
    fn problem(&self) -> TetProblem<'_> {
        TetProblem { vertices: &self.vertices, tets: &self.tets, conductivity: &self.k,
            source: &self.source, boundary: &self.boundary }
    }
    fn owner(&self, face: [usize; 3]) -> (usize, usize) {
        for (e, tet) in self.tets.iter().enumerate() {
            if face.iter().all(|v| tet.contains(v)) {
                return (e, tet.iter().position(|v| !face.contains(v)).unwrap());
            }
        }
        panic!("missing fixture face");
    }
}
fn slabs(n: usize, resistance: f64, heated: bool) -> Slabs {
    let per = (n+1)*(n+1)*(n+1);
    let mut vertices = Vec::new(); let mut tets = Vec::new();
    let mut k = Vec::new(); let mut source = Vec::new();
    for side in 0..2 {
        for z in 0..=n { for y in 0..=n { for x in 0..=n {
            vertices.push([side as f64+x as f64/n as f64, y as f64/n as f64, z as f64/n as f64]);
        }}}
        let index = |p: [usize; 3]| side*per+p[0]+(n+1)*(p[1]+(n+1)*p[2]);
        for z in 0..n { for y in 0..n { for x in 0..n {
            for axes in [[0,1,2],[0,2,1],[1,0,2],[1,2,0],[2,0,1],[2,1,0]] {
                let mut p = [x,y,z]; let mut tet = [index(p); 4];
                for (i, d) in axes.into_iter().enumerate() { p[d] += 1; tet[i+1] = index(p); }
                tets.push(tet); k.push(if side == 0 { 2.0 } else { 1.0 });
                source.push(if heated && side == 0 { 2.0 } else { 0.0 });
            }
        }}}
    }
    let mut incidence = BTreeMap::<[usize; 3], usize>::new();
    for tet in &tets { for opposite in 0..4 {
        let mut face = [0; 3]; let mut j = 0;
        for (i, &v) in tet.iter().enumerate() { if i != opposite { face[j] = v; j += 1; } }
        face.sort_unstable(); *incidence.entry(face).or_default() += 1;
    }}
    let boundary = incidence.into_iter().filter(|(_, count)| *count == 1).map(|(face, _)| {
        let at = |x| face.iter().all(|&v| vertices[v][0] == x);
        let condition = if at(1.0) {
            let partner = face.map(|v| if v < per { v+per-n } else { v-per+n });
            Bc::Contact { partner, resistance }
        } else if at(2.0) { Bc::Dirichlet([300.0; 3]) }
        else if at(0.0) && !heated { Bc::Dirichlet([400.0; 3]) }
        else { Bc::Neumann(0.0) };
        BoundaryFace { vertices: face, condition }
    }).collect();
    Slabs { vertices, tets, k, source, boundary, per }
}
fn contains(interval: Iv, truth: f64) {
    assert!(interval.lo <= truth && truth <= interval.hi, "{interval:?} excludes {truth}");
}

#[test]
fn exact_contact_patch_preserves_temperature_jump_materials_and_opposite_fluxes() {
    for r in [0.5, 2.5, 6.5] {
        let f = slabs(2, r, false); let q = 100.0/(1.5+r);
        let v: Vec<_> = f.vertices.iter().enumerate().map(|(i,p)|
            if i < f.per { 400.0-q*p[0]/2.0 } else { 300.0+q*(2.0-p[0]) }).collect();
        let energy = energy_bound(&f.problem(), &v, FluxBudget::default(), || true).unwrap();
        assert!(energy.energy_error_upper < 1e-6, "r={r}: {energy:?}");
        for flux in &energy.outward_flux_integrals {
            contains(flux.iter().copied().fold(Iv::zero(), Iv::add), 0.0);
        }
        // Every tet here has volume 1/48 and each contact triangle area 1/8.
        // With zero source, the conservative RT0 flux is constant in each
        // cell, as is the exact flux of this affine patch. Cauchy-Schwarz
        // therefore gives |F_h-F_exact| <= A*sqrt(k*eta_cell^2/V). The cell
        // majorant includes the volume defect (and possibly a contact term),
        // so this is an independent physical-flux allowance, evaluated
        // outward. An auxiliary flux box alone need not contain F_exact.
        let cell_volume = Iv::point(1.0).div_pos(Iv::point(48.0));
        let contact_area = Iv::point(0.125);
        for face in &f.boundary {
            if let Bc::Contact { partner, .. } = face.condition {
                let (a,i) = f.owner(face.vertices); let (b,j) = f.owner(partner);
                contains(energy.outward_flux_integrals[a][i].add(energy.outward_flux_integrals[b][j]), 0.0);
                let error = Iv::point(energy.cell_majorant_squared_upper[a])
                    .mul(Iv::point(f.k[a])).div_pos(cell_volume).sqrt().mul(contact_area);
                assert!(!error.is_unbounded());
                let physical_flux = energy.outward_flux_integrals[a][i]
                    .add(Iv { lo: -error.hi, hi: error.hi });
                contains(physical_flux, if face.vertices[0] < f.per { q/8.0 } else { -q/8.0 });
            }
        }
        let bounded = mean_bound(&f.problem(), &v, &vec![0.0; v.len()], FluxBudget::default(), || true).unwrap();
        contains(bounded.enclosure, 350.0+q/8.0);
        assert!(bounded.enclosure.hi-bounded.enclosure.lo < 1e-6);
    }
}

#[test]
fn contact_energy_cannot_be_dropped_or_replaced_by_temperature_continuity() {
    let f = slabs(1, 0.5, false);
    let v: Vec<_> = (0..f.vertices.len()).map(|i| if i < f.per { 400.0 } else { 300.0 }).collect();
    // Exact q=50. Volume error energy=3750; contact error energy=11250.
    // With q_h=0, the majorant is 20000. Omitting the jump term gives zero.
    let bounded = energy_bound(&f.problem(), &v, FluxBudget { max_iterations: 0, ..FluxBudget::default() }, || true).unwrap();
    assert!(bounded.energy_error_upper >= 15000.0_f64.sqrt());
    assert!(bounded.energy_error_upper < 142.0);
}

#[test]
fn a_heated_solid_can_be_anchored_only_through_finite_contact() {
    for n in [1,2,4] {
        let f = slabs(n, 0.5, true);
        let v: Vec<_> = f.vertices.iter().enumerate().map(|(i,p)|
            if i < f.per { 303.0+0.5*(1.0-p[0]*p[0]) } else { 300.0+2.0*(2.0-p[0]) }).collect();
        let z: Vec<_> = f.vertices.iter().enumerate().map(|(i,p)|
            if i < f.per { 2.0+0.25*(1.0-p[0]*p[0]) } else { 0.5*(4.0-p[0]*p[0]) }).collect();
        let bounded = mean_bound(&f.problem(), &v, &z, FluxBudget::default(), || true).unwrap();
        contains(bounded.enclosure, 302.0+1.0/6.0);
        assert!(bounded.integral.primal.energy_error_upper >= 1.0/(6.0_f64.sqrt()*n as f64));
    }
}

#[test]
fn affine_heating_uses_the_same_contact_operator_in_primal_and_dual_bounds() {
    let f = slabs(2, 0.5, true);
    let source: Vec<_> = f.tets.iter().map(|tet| tet.map(|i|
        if i < f.per { 6.0*f.vertices[i][0] } else { 0.0 })).collect();
    let k: Vec<_> = f.k.iter().map(|&x| [[x,0.0,0.0],[0.0,3.0,0.5],[0.0,0.5,2.0]]).collect();
    let v: Vec<_> = f.vertices.iter().enumerate().map(|(i,p)|
        if i < f.per { 304.5+0.5*(1.0-p[0]*p[0]*p[0]) } else { 300.0+3.0*(2.0-p[0]) }).collect();
    let z: Vec<_> = f.vertices.iter().enumerate().map(|(i,p)|
        if i < f.per { 2.0+0.25*(1.0-p[0]*p[0]) } else { 0.5*(4.0-p[0]*p[0]) }).collect();
    let p = AffineSourceTetProblem { vertices: &f.vertices, tets: &f.tets,
        conductivity: &k, source: &source, boundary: &f.boundary };
    let bounded = affine_source_mean_bound(&p, &v, &z, FluxBudget::default(), || true).unwrap();
    contains(bounded.enclosure, 303.1875);
    // The source, tensor and duplicate interface nodes are all unchanged.
    assert_eq!(p.source, source.as_slice());
}

#[test]
fn partner_order_and_declaration_order_do_not_change_the_conservation_graph() {
    let mut f = slabs(1, 0.5, false); let v = vec![300.0; f.vertices.len()];
    for row in &mut f.boundary { if let Bc::Dirichlet(ref mut values) = row.condition { *values = [300.0;3]; } }
    let first = mean_bound(&f.problem(), &v, &vec![0.0;v.len()], FluxBudget::default(), || true).unwrap();
    f.boundary.reverse();
    for row in &mut f.boundary {
        row.vertices.reverse();
        if let Bc::Contact { ref mut partner, .. } = row.condition { partner.reverse(); }
    }
    let second = mean_bound(&f.problem(), &v, &vec![0.0;v.len()], FluxBudget::default(), || true).unwrap();
    assert_eq!(first.enclosure, second.enclosure);
    assert_eq!(first.integral.primal.outward_flux_integrals, second.integral.primal.outward_flux_integrals);
}

#[test]
fn invalid_or_unanchored_contacts_and_cancellation_return_no_bound() {
    for mutation in 0..7 {
        let mut f = slabs(1, 0.5, true);
        let i = f.boundary.iter().position(|r| matches!(r.condition, Bc::Contact { .. })).unwrap();
        let Bc::Contact { partner, resistance } = f.boundary[i].condition else { unreachable!() };
        f.boundary[i].condition = match mutation {
            0 => Bc::Contact { partner, resistance: 0.0 },
            1 => Bc::Contact { partner, resistance: f64::NAN },
            2 => Bc::Contact { partner, resistance: resistance*2.0 },
            3 => Bc::Contact { partner: f.boundary[i].vertices, resistance },
            4 => Bc::Neumann(0.0),
            5 => Bc::Contact { partner: [usize::MAX;3], resistance },
            _ => Bc::Contact { partner, resistance },
        };
        if mutation == 6 { for face in &mut f.boundary {
            if matches!(face.condition, Bc::Dirichlet(_)) { face.condition = Bc::Neumann(0.0); }
        }}
        assert!(energy_bound(&f.problem(), &vec![300.0; f.vertices.len()], FluxBudget::default(), || true).is_err());
    }
    let f = slabs(2, 0.5, true); let v = vec![300.0; f.vertices.len()];
    for stop in [0,20,500] {
        let mut polls = 0;
        let result = mean_bound(&f.problem(), &v, &vec![0.0;v.len()], FluxBudget::default(), || { polls += 1; polls > 0 && polls <= stop });
        assert!(matches!(result, Err(TetError::Cancelled)));
    }
}
