use super::*;
use fs_couple::render::plate::impact::cavity::CavityCoupling;
use fs_couple::vibroacoustic::CavityModes;

fn piano(damping: bool) -> BridgeResponse {
    let course = super::super::geometry::demonstration_scale().unwrap()[48];
    let course = Course { unison:1, duplex_length_m:0.0, detune_cents:0.0, ..course };
    let mut board = vec![
        BoardMode { frequency_hz:170.0, damping_ratio:0.012, bridge:[0.0;88], volume:0.1 },
        BoardMode { frequency_hz:310.0, damping_ratio:0.018, bridge:[0.0;88], volume:-0.03 },
    ];
    board[0].bridge[48] = 0.08; board[1].bridge[48] = -0.04;
    BridgeResponse::new(&[course],&board,192_000,21_600.0,4,damping).unwrap()
}

struct Air {
    basis: CavityModes,
    overlaps: Vec<f64>,
    drag: Vec<f64>,
}
impl Air {
    fn standing(omega: f64, drag: f64) -> Self {
        Self { basis:CavityModes { omegas:vec![0.0,omega], lambdas:vec![0.3,0.15],
                interface:vec![vec![1.0];2], loss_factor:0.0, rho0:1.2, c0:343.0 },
            // Independent signed physical ports in the loaded board basis.
            overlaps:vec![0.06,0.04,-0.02,0.05], drag:vec![0.0,drag] }
    }
    fn model(&self) -> CavityCoupling {
        CavityCoupling::new(&self.basis,2,&self.overlaps,&self.drag).unwrap()
    }
    fn potential(&self, board: &[f64], acoustic: &[f64]) -> f64 {
        let count = self.basis.omegas.len(); let mut next = 0;
        self.basis.omegas.iter().enumerate().map(|(j,&omega)| {
            let a = self.basis.rho0*self.basis.c0*self.basis.c0/self.basis.lambdas[j];
            let mut volume = board.iter().enumerate().map(|(r,q)| self.overlaps[r*count+j]*q).sum::<f64>();
            if omega > 0.0 { volume += omega/a.sqrt()*acoustic[next]; next += 1; }
            0.5*a*volume*volume
        }).sum()
    }
}

// Independently polarize the original bank's physical potential plus the
// declared gas-volume energy. Do not use the new border, impedance or generic
// dynamic_stiffness assembly to construct this complete mechanical pencil.
fn full_pencil(m: &BridgeResponse, air: &Air, hz: f64, force: C64, exterior: Option<&[C64]>)
    -> Vec<C64> {
    let n = m.bank.modes.len(); let r = m.bank.board_count; let base = n+r;
    let positive = air.basis.omegas.iter().filter(|&&w| w > 0.0).count();
    let size = base+positive; let w = TAU*hz;
    let zero = vec![0.0;base]; let mut q = vec![0.0;size]; let mut energies = vec![0.0;size];
    let potential = |q: &[f64]| m.bank.energy_at(&q[..base],&zero)
        +air.potential(&q[n..base],&q[base..]);
    for i in 0..size { q[i] = 1.0; energies[i] = potential(&q); q[i] = 0.0; }
    let mut matrix = vec![C64::ZERO;size*size];
    for i in 0..size { for j in i..size {
        let stiffness = if i == j { 2.0*energies[i] } else {
            q[i] = 1.0; q[j] = 1.0;
            let k = potential(&q)-energies[i]-energies[j]; q[i] = 0.0; q[j] = 0.0; k
        };
        matrix[i*size+j] = C64::new(stiffness,0.0);
        matrix[j*size+i] = C64::new(stiffness,0.0);
    } }
    for i in 0..size { matrix[i*size+i].re -= w*w; }
    for i in 0..n { matrix[i*size+i].im -= w*m.string_c[i]; }
    for i in 0..r { for j in 0..r {
        let at = (n+i)*size+n+j;
        matrix[at].im -= w*m.board_c[i*r+j];
        if let Some(z) = exterior { matrix[at] = matrix[at]+C64::new(0.0,-w)*z[i*r+j]; }
    } }
    let mut next = base;
    for (j,&omega) in air.basis.omegas.iter().enumerate() {
        if omega > 0.0 { matrix[next*size+next].im -= w*air.drag[j]; next += 1; }
    }
    let mut rhs = vec![C64::ZERO;size];
    for (v,g) in rhs[n..base].iter_mut().zip(m.bridge_row(69).unwrap()) { *v = force.scale(*g); }
    lu_complex(&matrix,size).unwrap().solve(&mut rhs); rhs
}

fn compare(m: &BridgeResponse, air: &Air, hz: f64, exterior: Option<&[C64]>) -> Response {
    let force = C64::new(1.0,0.3);
    let expected = full_pencil(m,air,hz,force,exterior);
    let actual = m.solve(hz,69,force,exterior).unwrap();
    let scale = expected.iter().map(|q|q.abs()).fold(0.0_f64,f64::max);
    assert_eq!(actual.cavity_displacement.len(),air.basis.omegas.iter().filter(|&&w|w > 0.0).count());
    for (actual,expected) in actual.string_displacement.iter().chain(&actual.board_displacement)
        .chain(&actual.cavity_displacement).zip(&expected) {
        assert!((*actual-*expected).abs() < 1e-7*scale,
            "complete string/board/gas energy pencil differs at {hz} Hz: {actual:?} vs {expected:?}");
    }
    let w = TAU*hz; let mut positive = 0;
    let expected_loss = air.basis.omegas.iter().enumerate().map(|(j,&omega)| {
        if omega == 0.0 { 0.0 } else {
            let z = expected[m.bank.q.len()+positive]; positive += 1;
            0.5*w*w*air.drag[j]*z.abs().powi(2)
        }
    }).sum::<f64>();
    assert!((actual.cavity_loss_w-expected_loss).abs() < 1e-11+1e-7*expected_loss);
    assert!(actual.backward_error < 1e-11);
    assert!(actual.power_defect_w.abs() < 1e-10+1e-7*actual.input_w.abs());
    actual
}

#[test]
fn cavity_harmonics_match_the_complete_energy_pencil_and_keep_air_and_exterior_loss_separate() {
    let mut model = piano(true); let air = Air::standing(TAU*230.0,40.0);
    model.configure_cavity(&air.model()).unwrap();
    let exterior = [C64::new(3.0,-0.2),C64::new(0.4,0.1),
        C64::new(0.4,0.1),C64::new(2.0,-0.3)];
    for hz in [100.0,230.0,301.0,470.0] {
        let bare = compare(&model,&air,hz,None);
        assert_eq!(bare.radiation_w,0.0);
        assert!(bare.cavity_loss_w > 0.0 && bare.board_loss_w > 0.0 && bare.string_loss_w > 0.0);
        let loaded = compare(&model,&air,hz,Some(&exterior));
        assert!(loaded.cavity_loss_w > 0.0 && loaded.radiation_w > 0.0);
        assert!((loaded.input_w-loaded.board_loss_w-loaded.string_loss_w
            -loaded.cavity_loss_w-loaded.radiation_w).abs() < 1e-10);
    }
}

#[test]
fn lossless_air_pole_and_coincident_string_pole_keep_their_actual_inertias_in_the_same_solve() {
    for coincide in [false,true] {
        let mut model = piano(false);
        let omega = if coincide { model.bank.modes[0].omega } else { TAU*230.0 };
        let air = Air::standing(omega,0.0); let query = omega/TAU;
        model.configure_cavity(&air.model()).unwrap();
        assert!(air.model().impedance(omega).is_err(),"the eliminated fixed-wall pole really is singular");
        for perturbation in [0.0,-1e-9,1e-9] {
            let actual = compare(&model,&air,query*(1.0+perturbation),None);
            assert_eq!(actual.retained_string_poles,usize::from(coincide));
            assert_eq!(actual.cavity_loss_w,0.0);
            assert_eq!(actual.string_loss_w,0.0); assert_eq!(actual.board_loss_w,0.0);
            assert!(actual.input_w.abs() < 1e-10);
        }
        assert!(model.bank.q.iter().chain(&model.bank.v).all(|&x| x == 0.0));
        let zero = model.solve(query,69,C64::ZERO,None).unwrap();
        assert!(zero.board_displacement.iter().chain(&zero.string_displacement)
            .chain(&zero.cavity_displacement).all(|&x| x == C64::ZERO));
    }
}

#[test]
fn uniform_cavity_harmonics_add_only_the_original_compression_spring() {
    let mut model = piano(false); let mut air = Air::standing(TAU*230.0,0.0);
    air.basis.omegas.truncate(1); air.basis.lambdas.truncate(1); air.basis.interface.truncate(1);
    air.overlaps = vec![0.06,-0.02]; air.drag.truncate(1);
    model.configure_cavity(&air.model()).unwrap();
    for hz in [100.0,170.0,301.0] {
        let actual = compare(&model,&air,hz,None);
        assert!(actual.cavity_displacement.is_empty());
        assert_eq!(actual.cavity_loss_w,0.0);
    }
}
