//! The regional P1 volume functional on the unchanged native thermal mesh.
//! A node indicator is not this load: shared vertices carry only the volumes
//! of selected cells, and duplicated contact traces remain distinct vertices.
use fs_conduction::ConductionMesh;
use super::{Cx, SolveRefusal, bad, finite, poll, zeros};

pub(super) const OUTPUT: &str = "temperature-volume-mean-adjoint";
pub(super) const SCOPE: &str = "Estimated derivative of the P1 volume-mean temperature over the existing temperature-max requirement's region on the final accepted native mesh. The full assembly supplies the thermal Jacobian; selection changes the objective load only, not the physical domain or its boundaries. Each selected tetrahedron contributes V_e/(4 V_region) to its four original vertices, including prescribed nodes. Contact traces remain independent. This is a spatial mean, not a statistical expectation, a maximum-temperature derivative, or a continuum mean-temperature certificate.";

pub(super) struct Mean {
    pub(super) weights: Vec<f64>,
    pub(super) value_k: f64,
    pub(super) volume_m3: f64,
    pub(super) elements: usize,
}

// Compensate positive volume and weighted-temperature sums without changing
// source element order or calling a separate geometry/numerical producer.
#[derive(Default)]
struct Sum { value: f64, correction: f64 }
impl Sum {
    fn add(&mut self, value: f64) -> Result<(), SolveRefusal> {
        let y = finite(value - self.correction)?;
        let total = finite(self.value + y)?;
        self.correction = finite((total - self.value) - y)?;
        self.value = total;
        Ok(())
    }
}

pub(super) fn prepare(
    cx: &Cx<'_>, mesh: &ConductionMesh, labels: &[u32], region: u32, temperature: &[f64],
) -> Result<Mean, SolveRefusal> {
    poll(cx)?;
    if labels.len() != mesh.element_count() || temperature.len() != mesh.vertex_count() {
        return Err(bad("volume-mean objective does not match the retained mesh, labels or field"));
    }
    let mut volume = Sum::default();
    let mut elements = 0;
    for (e, &label) in labels.iter().enumerate() {
        if e % 512 == 0 { poll(cx)?; }
        if label != region { continue; }
        let v = mesh.element_volume(e);
        if !v.is_finite() || v <= 0.0 { return Err(bad("invalid selected element volume")); }
        volume.add(v)?;
        elements += 1;
    }
    if elements == 0 || volume.value <= 0.0 {
        return Err(bad("volume-mean objective requires a nonempty positive-volume region"));
    }
    let mut weights = zeros(mesh.vertex_count())?;
    for (e, tet) in mesh.complex().tets.iter().enumerate() {
        if e % 512 == 0 { poll(cx)?; }
        if labels[e] != region { continue; }
        // Divide BEFORE /4 so representable small element volumes are not
        // prematurely rounded to zero. Underflow is a refusal, not lost mass.
        let weight = finite((mesh.element_volume(e) / volume.value) / 4.0)?;
        if weight <= 0.0 { return Err(bad("volume-mean element weight underflow")); }
        for &vertex in tet {
            let value = weights.get_mut(vertex as usize)
                .ok_or_else(|| bad("volume-mean element refers to a missing vertex"))?;
            *value = finite(*value + weight)?;
        }
    }
    let mut mean = Sum::default();
    for (i, (&weight, &value)) in weights.iter().zip(temperature).enumerate() {
        if i % 512 == 0 { poll(cx)?; }
        // Validate the complete retained field, including unselected nodes.
        finite(value)?;
        mean.add(finite(weight * value)?)?;
    }
    poll(cx)?;
    Ok(Mean { weights, value_k: mean.value, volume_m3: volume.value, elements })
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_alloc::{ArenaConfig, ArenaPool};
    use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
    use fs_rep_mesh::TetComplex;

    fn with_cx(f: impl FnOnce(&Cx<'_>, &CancelGate)) {
        let gate = CancelGate::new_clock_free();
        ArenaPool::new(ArenaConfig::default()).scope(|arena| {
            f(&Cx::new(&gate, arena, StreamKey { seed: 7, kernel_id: 47, tile: 0, iteration: 0 },
                Budget::INFINITE, ExecMode::Deterministic), &gate)
        });
    }
    fn mesh() -> ConductionMesh {
        // Volumes 1/6 and 8/6, sharing ONE actual vertex, not overlapping cells.
        ConductionMesh::new(TetComplex::from_tets(7, vec![[0,1,2,3], [0,4,6,5]]),
            vec![[0.,0.,0.], [1.,0.,0.], [0.,1.,0.], [0.,0.,1.],
                [-2.,0.,0.], [0.,-2.,0.], [0.,0.,-2.]]).unwrap()
    }
    const TEMPERATURE: [f64;7] = [300.,302.,304.,306.,310.,320.,330.];
    fn near(a: f64, b: f64) { assert!((a-b).abs() < 2e-12, "{a} != {b}"); }

    #[test]
    fn regional_mean_uses_original_volumes_not_nodes_or_equal_cell_weights() {
        with_cx(|cx,_| {
            let mesh = mesh();
            let goal = prepare(cx,&mesh,&[9,9],9,&TEMPERATURE).unwrap();
            near(goal.value_k,(303.0+8.0*315.0)/9.0);
            near(goal.volume_m3,1.5); assert_eq!(goal.elements,2);
            near(goal.weights[0],0.25);
            for &v in &[1,2,3] { near(goal.weights[v],1.0/36.0); }
            for &v in &[4,5,6] { near(goal.weights[v],2.0/9.0); }
            assert!((goal.value_k-TEMPERATURE.iter().sum::<f64>()/7.0).abs()>1.0);
            assert!((goal.value_k-(303.0+315.0)/2.0).abs()>1.0);
            let reordered = ConductionMesh::new(TetComplex::from_tets(7,
                vec![[0,4,6,5],[0,1,2,3]]),mesh.positions().to_vec()).unwrap();
            let other = prepare(cx,&reordered,&[9,9],9,&TEMPERATURE).unwrap();
            near(other.value_k,goal.value_k);
            for (&a,&b) in other.weights.iter().zip(&goal.weights) { near(a,b); }
        });
    }

    #[test]
    fn selection_changes_the_goal_only_and_retains_direct_fixed_node_weights() {
        with_cx(|cx,_| {
            let mesh=mesh();
            let goal=prepare(cx,&mesh,&[9,10],9,&TEMPERATURE).unwrap();
            near(goal.value_k,303.0); assert_eq!(goal.elements,1);
            assert_eq!(goal.weights,[0.25,0.25,0.25,0.25,0.0,0.0,0.0]);
            // The load exists at a selected prescribed node as well as at free
            // nodes. The existing prescribed-temperature owner consumes it.
            let mut changed=TEMPERATURE; changed[0]+=0.01;
            let altered=prepare(cx,&mesh,&[9,10],9,&changed).unwrap();
            near((altered.value_k-goal.value_k)/0.01,goal.weights[0]);
            near(prepare(cx,&mesh,&[9,10],9,&[300.0;7]).unwrap().value_k,300.0);
        });
    }

    #[test]
    fn absent_regions_malformed_fields_and_cancellation_do_not_return_a_goal() {
        with_cx(|cx,gate| {
            let mesh=mesh();
            assert!(prepare(cx,&mesh,&[9,10],11,&TEMPERATURE).is_err());
            assert!(prepare(cx,&mesh,&[9],9,&TEMPERATURE).is_err());
            assert!(prepare(cx,&mesh,&[9,10],9,&TEMPERATURE[..6]).is_err());
            let mut changed=TEMPERATURE; changed[6]=f64::NAN;
            assert!(prepare(cx,&mesh,&[9,10],9,&changed).is_err());
            gate.request();
            assert_eq!(prepare(cx,&mesh,&[9,10],9,&TEMPERATURE).err().unwrap().code,
                "cli-solve-cancelled");
        });
    }
}
