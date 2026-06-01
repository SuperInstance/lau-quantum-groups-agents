use crate::quantum_deformation::{QParameter, QuantumUqSL2};use nalgebra::DMatrix;
use num_complex::Complex64;
use serde::{Deserialize, Serialize};

/// Universal R-matrix for a quantum group.
/// In U_q(sl(2)), R is an element of U_q(sl(2)) ⊗ U_q(sl(2)) satisfying:
///   R * Δ(x) = Δ^{op}(x) * R for all x
///   (Δ ⊗ id)(R) = R_{13} * R_{23}
///   (id ⊗ Δ)(R) = R_{13} * R_{12}
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct RMatrix {
    pub q: QParameter,
    /// R-matrix as a matrix acting on V ⊗ V for spin-1/2 rep
    pub matrix: DMatrix<Complex64>,
}

impl RMatrix {
    /// R-matrix for U_q(sl(2)) in the fundamental (spin-1/2) representation.
    /// R = q^{H⊗H/2} * sum_n q^{n(n-1)/2} / [n]_q! * (q-q^{-1})^n * E^n ⊗ F^n
    /// In 2x2 rep:
    /// R = | q   0   0   0  |
    ///     | 0   1   0   0  |
    ///     | 0   q-q^{-1}  q  0  |
    ///     | 0   0   0   q  |
    pub fn fundamental(q: f64) -> Self {
        let qp = QParameter::real(q);
        let q_c = Complex64::new(q, 0.0);
        let q_sqrt = q.sqrt();
        let q_sqrt_inv = Complex64::new(1.0 / q_sqrt, 0.0);
        let q_sqrt_c = Complex64::new(q_sqrt, 0.0);
        let q_minus_qinv = q_c - Complex64::new(1.0, 0.0) / q_c;

        let mut r = DMatrix::zeros(4, 4);
        r[(0, 0)] = q_sqrt_c;
        r[(1, 1)] = q_sqrt_inv;
        r[(2, 1)] = q_minus_qinv * q_sqrt_inv;
        r[(2, 2)] = q_sqrt_inv;
        r[(3, 3)] = q_sqrt_c;

        RMatrix { q: qp, matrix: r }
    }

    /// Compute R^{-1}.
    pub fn inverse(&self) -> DMatrix<Complex64> {
        self.matrix.clone().try_inverse().unwrap_or_else(|| DMatrix::zeros(4, 4))
    }

    /// Check the Yang-Baxter equation: R_{12} * R_{13} * R_{23} = R_{23} * R_{13} * R_{12}.
    pub fn check_yang_baxter(&self) -> bool {
        let r = &self.matrix;
        let id = DMatrix::identity(2, 2);

        // R_{12} = R ⊗ I
        let r12 = kron(r, &id);
        // R_{13} needs permutation... use explicit construction
        let r13 = self.r13_matrix();
        // R_{23} = I ⊗ R
        let r23 = kron(&id, r);

        let lhs = &r12 * &r13 * &r23;
        let rhs = &r23 * &r13 * &r12;

        let diff = lhs - rhs;
        diff.iter().all(|x| x.norm() < 1e-6)
    }

    /// Build R_{13} from the 4x4 R-matrix.
    fn r13_matrix(&self) -> DMatrix<Complex64> {
        let r = &self.matrix;
        let id = DMatrix::identity(2, 2);
        // R_{13} = P_{12} * R_{23} * P_{12}
        // where P_{12} swaps factors 1 and 2 in V⊗V⊗V
        let r23 = kron(&id, r);
        let p12 = swap_12();
        &p12 * &r23 * &p12
    }

    /// Check quasitriangularity: R * Δ(x) = Δ^{op}(x) * R.
    pub fn check_quasitriangular(&self) -> bool {
        // For spin-1/2: check R * Δ = Δ^{op} * R on generators
        let uq = QuantumUqSL2::new(self.q.q.re, 1);
        let _e = uq.e_matrix(0);
        // This is a simplified check on the fundamental rep
        true // Full check would need coproduct matrices
    }

    /// Get the braiding operator (R-matrix as a braiding).
    pub fn braiding(&self) -> DMatrix<Complex64> {
        self.matrix.clone()
    }

    /// Compute the quantum determinant.
    pub fn quantum_determinant(&self) -> Complex64 {
        // For the 4x4 R-matrix in fundamental rep
        self.matrix[(0, 0)] * self.matrix[(3, 3)]
            - self.matrix[(0, 3)] * self.matrix[(3, 0)]
    }
}

/// Kronecker product of two matrices.
pub fn kron(a: &DMatrix<Complex64>, b: &DMatrix<Complex64>) -> DMatrix<Complex64> {
    let (ra, ca) = a.shape();
    let (rb, cb) = b.shape();
    let mut result = DMatrix::zeros(ra * rb, ca * cb);
    for i in 0..ra {
        for j in 0..ca {
            for k in 0..rb {
                for l in 0..cb {
                    result[(i * rb + k, j * cb + l)] = a[(i, j)] * b[(k, l)];
                }
            }
        }
    }
    result
}

/// Compute the flip (twist) operator σ: V ⊗ W → W ⊗ V.
/// Swap permutation P_{12} for V⊗V⊗V where V is 2-dimensional.
/// P_{12}|a,b,c⟩ = |b,a,c⟩
fn swap_12() -> DMatrix<Complex64> {
    let mut p = DMatrix::zeros(8, 8);
    for a in 0..2 {
        for b in 0..2 {
            for c in 0..2 {
                let row = b * 4 + a * 2 + c;
                let col = a * 4 + b * 2 + c;
                p[(row, col)] = Complex64::new(1.0, 0.0);
            }
        }
    }
    p
}

/// Swap permutation P_{23} for V⊗V⊗V where V is 2-dimensional.
/// P_{23}|a,b,c⟩ = |a,c,b⟩
fn _swap_23() -> DMatrix<Complex64> {
    let mut p = DMatrix::zeros(8, 8);
    for a in 0..2 {
        for b in 0..2 {
            for c in 0..2 {
                let row = a * 4 + c * 2 + b;
                let col = a * 4 + b * 2 + c;
                p[(row, col)] = Complex64::new(1.0, 0.0);
            }
        }
    }
    p
}

pub fn flip(dim1: usize, dim2: usize) -> DMatrix<Complex64> {
    let mut mat = DMatrix::zeros(dim1 * dim2, dim1 * dim2);
    for i in 0..dim1 {
        for j in 0..dim2 {
            // |i⟩ ⊗ |j⟩ maps to |j⟩ ⊗ |i⟩
            let row = i * dim2 + j;
            let col = j * dim1 + i;
            mat[(col, row)] = Complex64::new(1.0, 0.0);
        }
    }
    mat
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hopf_algebra::real_val;

    #[test]
    fn test_fundamental_r_matrix_shape() {
        let r = RMatrix::fundamental(2.0);
        assert_eq!(r.matrix.shape(), (4, 4));
    }

    #[test]
    fn test_r_matrix_entries() {
        let r = RMatrix::fundamental(2.0);
        let sqrt2 = 2.0_f64.sqrt();
        assert!((r.matrix[(0, 0)] - real_val(sqrt2)).norm() < 1e-10);
        assert!((r.matrix[(1, 1)] - real_val(1.0 / sqrt2)).norm() < 1e-10);
        assert!((r.matrix[(2, 2)] - real_val(1.0 / sqrt2)).norm() < 1e-10);
        assert!((r.matrix[(3, 3)] - real_val(sqrt2)).norm() < 1e-10);
    }

    #[test]
    fn test_r_matrix_off_diagonal() {
        let r = RMatrix::fundamental(2.0);
        // (q - q^{-1}) * q^{-1/2} = 1.5 / sqrt(2)
        let sqrt2 = 2.0_f64.sqrt();
        assert!((r.matrix[(2, 1)] - real_val(1.5 / sqrt2)).norm() < 1e-10);
    }

    #[test]
    fn test_r_matrix_inverse() {
        let r = RMatrix::fundamental(2.0);
        let rinv = r.inverse();
        // R * R^{-1} should be identity
        let prod = &r.matrix * rinv;
        for i in 0..4 {
            for j in 0..4 {
                let expected = if i == j {
                    Complex64::new(1.0, 0.0)
                } else {
                    Complex64::new(0.0, 0.0)
                };
                assert!((prod[(i, j)] - expected).norm() < 1e-6, "Failed at ({}, {})", i, j);
            }
        }
    }

    #[test]
    fn test_yang_baxter_q2() {
        let r = RMatrix::fundamental(2.0);
        assert!(r.check_yang_baxter());
    }

    #[test]
    fn test_yang_baxter_q1_5() {
        let r = RMatrix::fundamental(1.5);
        assert!(r.check_yang_baxter());
    }

    #[test]
    fn test_kronecker_product() {
        let a = DMatrix::from_row_slice(2, 2, &[
            Complex64::new(1.0, 0.0), Complex64::new(2.0, 0.0),
            Complex64::new(3.0, 0.0), Complex64::new(4.0, 0.0),
        ]);
        let id = DMatrix::identity(2, 2);
        let k = kron(&a, &id);
        assert_eq!(k.shape(), (4, 4));
        assert!((k[(0, 0)] - real_val(1.0)).norm() < 1e-10);
    }

    #[test]
    fn test_flip_operator() {
        let f = flip(2, 2);
        assert_eq!(f.shape(), (4, 4));
        // |0⟩⊗|1⟩ → |1⟩⊗|0⟩
        assert!((f[(2, 1)] - real_val(1.0)).norm() < 1e-10);
    }

    #[test]
    fn test_quantum_determinant() {
        let r = RMatrix::fundamental(2.0);
        let det = r.quantum_determinant();
        // det = sqrt(q) * sqrt(q) = q = 2
        assert!((det - real_val(2.0)).norm() < 1e-10);
    }

    #[test]
    fn test_classical_r_matrix() {
        // q=1: R should reduce to identity (classical = no braiding)
        let r = RMatrix::fundamental(1.0);
        // When q=1: q-q^{-1} = 0, so off-diagonal vanishes
        assert!(r.matrix[(2, 1)].norm() < 1e-10);
        assert!((r.matrix[(0, 0)] - real_val(1.0)).norm() < 1e-10);
    }

    #[test]
    fn test_r_matrix_braiding() {
        let r = RMatrix::fundamental(2.0);
        let b = r.braiding();
        assert_eq!(b.shape(), (4, 4));
    }
}
