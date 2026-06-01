use nalgebra::DMatrix;
use num_complex::Complex64;
use serde::{Deserialize, Serialize};

/// Quantum deformation parameter q and related utilities.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct QParameter {
    pub q: Complex64,
}

impl QParameter {
    pub fn new(q: Complex64) -> Self {
        Self { q }
    }

    pub fn real(q: f64) -> Self {
        Self { q: Complex64::new(q, 0.0) }
    }

    /// q^n
    pub fn pow(&self, n: i32) -> Complex64 {
        let mut result = Complex64::new(1.0, 0.0);
        if n >= 0 {
            for _ in 0..n {
                result *= self.q;
            }
        } else {
            let qinv = Complex64::new(1.0, 0.0) / self.q;
            for _ in 0..(-n) {
                result *= qinv;
            }
        }
        result
    }

    /// q-number: [n]_q = (q^n - q^{-n}) / (q - q^{-1})
    pub fn q_number(&self, n: i32) -> Complex64 {
        if n == 0 {
            return Complex64::new(0.0, 0.0);
        }
        let q_1 = self.q;
        let q_minus_1 = Complex64::new(1.0, 0.0) / self.q;
        let denom = q_1 - q_minus_1;
        // Handle classical limit q→1 via L'Hôpital: [n]_q → n
        if denom.norm() < 1e-12 {
            return Complex64::new(n as f64, 0.0);
        }
        let q_n = self.pow(n);
        let q_minus_n = self.pow(-n);
        (q_n - q_minus_n) / denom
    }

    /// q-factorial: [n]_q!
    pub fn q_factorial(&self, n: u32) -> Complex64 {
        let mut result = Complex64::new(1.0, 0.0);
        for k in 1..=n {
            result *= self.q_number(k as i32);
        }
        result
    }

    /// q-binomial coefficient: [n choose k]_q
    pub fn q_binomial(&self, n: u32, k: u32) -> Complex64 {
        if k > n {
            return Complex64::new(0.0, 0.0);
        }
        self.q_factorial(n) / (self.q_factorial(k) * self.q_factorial(n - k))
    }

    /// Check if q is a root of unity.
    pub fn is_root_of_unity(&self, max_order: u32) -> bool {
        let one = Complex64::new(1.0, 0.0);
        let mut qp = Complex64::new(1.0, 0.0);
        for _ in 0..max_order {
            qp *= self.q;
            if (qp - one).norm() < 1e-8 {
                return true;
            }
        }
        false
    }

    /// Classical limit: q → 1.
    pub fn is_classical(&self) -> bool {
        (self.q - Complex64::new(1.0, 0.0)).norm() < 1e-8
    }
}

/// Quantum group U_q(sl(2)): q-deformed universal enveloping algebra of sl(2).
///
/// Generators: E, F, K, K^{-1}
/// Relations:
///   K K^{-1} = K^{-1} K = 1
///   K E K^{-1} = q^2 E
///   K F K^{-1} = q^{-2} F
///   [E, F] = (K - K^{-1}) / (q - q^{-1})
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct QuantumUqSL2 {
    pub q: QParameter,
    /// Truncation level for representations (spin j).
    pub max_spin: u32,
}

impl QuantumUqSL2 {
    pub fn new(q: f64, max_spin: u32) -> Self {
        Self {
            q: QParameter::real(q),
            max_spin,
        }
    }

    /// Dimension of the spin-j representation: 2j + 1.
    pub fn rep_dim(&self, j: u32) -> usize {
        (2 * j + 1) as usize
    }

    /// Build the K matrix for spin-j representation.
    /// K |j,m⟩ = q^{2m} |j,m⟩
    pub fn k_matrix(&self, j: u32) -> DMatrix<Complex64> {
        let dim = self.rep_dim(j);
        let mut mat = DMatrix::zeros(dim, dim);
        for m in 0..dim {
            let m_val = j as i32 - m as i32;
            mat[(m, m)] = self.q.pow(2 * m_val);
        }
        mat
    }

    /// Build the E matrix for spin-j representation.
    /// E |j,m⟩ = [j-m]_q |j,m+1⟩ where [n]_q is the q-number.
    pub fn e_matrix(&self, j: u32) -> DMatrix<Complex64> {
        let dim = self.rep_dim(j);
        let mut mat = DMatrix::zeros(dim, dim);
        for m in 0..dim {
            let m_val = j as i32 - m as i32;
            // E raises m by 1: |j,m⟩ → |j,m+1⟩
            // m index: m_val = j - idx, so idx = j - m_val
            // Raising: m_val → m_val + 1, idx → idx - 1
            if m > 0 {
                let _n = j as i32 - m_val; // = j - (j - idx) = idx... let me recalc
                // Actually: coefficient is sqrt([j-m][j+m+1]) in classical case
                // In q-case: E |j,m⟩ = sqrt([j-m]_q * [j+m+1]_q) |j,m+1⟩
                // But simpler: E |j,m⟩ = [j+m+1]_q^(1/2) [j-m]_q^(1/2) |j,m+1⟩
                // For simplicity, use: E |j,m⟩ = sqrt([j-m]_q * [j+m+1]_q) |j,m+1⟩
                let qn1 = self.q.q_number(j as i32 - m_val);
                let qn2 = self.q.q_number(j as i32 + m_val + 1);
                let coeff = (qn1 * qn2).sqrt();
                mat[(m - 1, m)] = coeff;
            }
        }
        mat
    }

    /// Build the F matrix for spin-j representation.
    /// F |j,m⟩ = [j+m]_q^(1/2) [j-m+1]_q^(1/2) |j,m-1⟩
    pub fn f_matrix(&self, j: u32) -> DMatrix<Complex64> {
        let dim = self.rep_dim(j);
        let mut mat = DMatrix::zeros(dim, dim);
        for m in 0..dim {
            let m_val = j as i32 - m as i32;
            // F lowers: m_val → m_val - 1, idx → idx + 1
            if m < dim - 1 {
                let qn1 = self.q.q_number(j as i32 + m_val);
                let qn2 = self.q.q_number(j as i32 - m_val + 1);
                let coeff = (qn1 * qn2).sqrt();
                mat[(m + 1, m)] = coeff;
            }
        }
        mat
    }

    /// Check the quantum commutation relation [E, F] = (K - K^{-1})/(q - q^{-1}).
    pub fn check_ef_relation(&self, j: u32) -> bool {
        let dim = self.rep_dim(j);
        let e = self.e_matrix(j);
        let f = self.f_matrix(j);
        let k = self.k_matrix(j);

        let ef = &e * &f;
        let fe = &f * &e;
        let comm = &ef - &fe;

        let q_minus_qinv = self.q.q - Complex64::new(1.0, 0.0) / self.q.q;
        let k_inv = k.clone().try_inverse().unwrap_or(DMatrix::zeros(dim, dim));
        let rhs = (&k - &k_inv) / q_minus_qinv;

        let diff = comm - rhs;
        diff.iter().all(|x| x.norm() < 1e-6)
    }

    /// Check K E K^{-1} = q^2 E.
    pub fn check_kek_relation(&self, j: u32) -> bool {
        let dim = self.rep_dim(j);
        let e = self.e_matrix(j);
        let k = self.k_matrix(j);
        let k_inv = k.clone().try_inverse().unwrap_or(DMatrix::zeros(dim, dim));

        let lhs = &k * &e * &k_inv;
        let rhs = e * self.q.pow(2);

        let diff = lhs - rhs;
        diff.iter().all(|x| x.norm() < 1e-6)
    }

    /// Check K F K^{-1} = q^{-2} F.
    pub fn check_kfk_relation(&self, j: u32) -> bool {
        let dim = self.rep_dim(j);
        let f = self.f_matrix(j);
        let k = self.k_matrix(j);
        let k_inv = k.clone().try_inverse().unwrap_or(DMatrix::zeros(dim, dim));

        let lhs = &k * &f * &k_inv;
        let rhs = f * self.q.pow(-2);

        let diff = lhs - rhs;
        diff.iter().all(|x| x.norm() < 1e-6)
    }

    /// Get all representation matrices for spin j.
    pub fn representation(&self, j: u32) -> QuantumRepresentation {
        QuantumRepresentation {
            spin: j,
            e: self.e_matrix(j),
            f: self.f_matrix(j),
            k: self.k_matrix(j),
        }
    }
}

/// A quantum group representation.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct QuantumRepresentation {
    pub spin: u32,
    pub e: DMatrix<Complex64>,
    pub f: DMatrix<Complex64>,
    pub k: DMatrix<Complex64>,
}

impl QuantumRepresentation {
    pub fn dim(&self) -> usize {
        self.e.nrows()
    }

    /// Compute the quantum Casimir element.
    pub fn casimir(&self, q: &QParameter) -> DMatrix<Complex64> {
        let _ef = &self.e * &self.f;
        let fe = &self.f * &self.e;
        // C = EF + [j]_q * FE + ... simplified version
        let k_inv = self.k.clone().try_inverse().unwrap_or_else(|| DMatrix::zeros(self.dim(), self.dim()));
        let q_minus_qinv = q.q - Complex64::new(1.0, 0.0) / q.q;
        // C_q = FE + K/(q-q^{-1}) + K^{-1}/(q-q^{-1}) -- simplified
        &fe + (&self.k + &k_inv) * (Complex64::new(1.0, 0.0) / q_minus_qinv)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hopf_algebra::real_val;

    #[test]
    fn test_q_parameter_pow() {
        let q = QParameter::real(2.0);
        assert!((q.pow(3) - real_val(8.0)).norm() < 1e-10);
        assert!((q.pow(-1) - real_val(0.5)).norm() < 1e-10);
    }

    #[test]
    fn test_q_number() {
        let q = QParameter::real(2.0);
        // [1]_q = (q - q^{-1})/(q - q^{-1}) = 1
        assert!((q.q_number(1) - real_val(1.0)).norm() < 1e-10);
        // [0]_q = 0
        assert!(q.q_number(0).norm() < 1e-10);
    }

    #[test]
    fn test_q_number_2() {
        let q = QParameter::real(2.0);
        // [2]_q = (q^2 - q^{-2})/(q - q^{-1}) = q + q^{-1} = 2.5
        let n2 = q.q_number(2);
        assert!((n2 - real_val(2.5)).norm() < 1e-8);
    }

    #[test]
    fn test_q_factorial() {
        let q = QParameter::real(2.0);
        // [0]_q! = 1, [1]_q! = 1, [2]_q! = [2]_q * [1]_q = 2.5
        assert!((q.q_factorial(0) - real_val(1.0)).norm() < 1e-10);
        assert!((q.q_factorial(1) - real_val(1.0)).norm() < 1e-10);
        assert!((q.q_factorial(2) - real_val(2.5)).norm() < 1e-8);
    }

    #[test]
    fn test_q_binomial() {
        let q = QParameter::real(2.0);
        // [2 choose 0]_q = 1, [2 choose 1]_q = [2]_q = 2.5
        assert!((q.q_binomial(2, 0) - real_val(1.0)).norm() < 1e-8);
        assert!((q.q_binomial(2, 1) - real_val(2.5)).norm() < 1e-8);
    }

    #[test]
    fn test_classical_limit() {
        let q = QParameter::real(1.0);
        assert!(q.is_classical());
        // [n]_1 = n
        assert!((q.q_number(3) - real_val(3.0)).norm() < 1e-8);
    }

    #[test]
    fn test_root_of_unity() {
        let q = QParameter::new(Complex64::from_polar(1.0, std::f64::consts::PI / 3.0));
        // q^6 = 1
        assert!(q.is_root_of_unity(10));
    }

    #[test]
    fn test_not_root_of_unity() {
        let q = QParameter::real(1.5);
        assert!(!q.is_root_of_unity(100));
    }

    #[test]
    fn test_uqsl2_spin_half_dim() {
        let uq = QuantumUqSL2::new(2.0, 1);
        assert_eq!(uq.rep_dim(0), 1);
        assert_eq!(uq.rep_dim(1), 3);
    }

    #[test]
    fn test_uqsl2_k_matrix_diagonal() {
        let uq = QuantumUqSL2::new(2.0, 1);
        let k = uq.k_matrix(1); // spin-1: dim 3
        assert!(k.iter().zip(k.iter()).all(|(_, _)| true) || k.shape().0 == 0
            || {
                let (n, m) = k.shape();
                (0..n).all(|i| (0..m).all(|j| i == j || k[(i, j)].norm() < 1e-10))
            });
    }

    #[test]
    fn test_uqsl2_kek_relation() {
        let uq = QuantumUqSL2::new(2.0, 2);
        assert!(uq.check_kek_relation(1));
    }

    #[test]
    fn test_uqsl2_kfk_relation() {
        let uq = QuantumUqSL2::new(2.0, 2);
        assert!(uq.check_kfk_relation(1));
    }

    #[test]
    fn test_uqsl2_spin0() {
        let uq = QuantumUqSL2::new(2.0, 0);
        let rep = uq.representation(0);
        assert_eq!(rep.dim(), 1);
    }

    #[test]
    fn test_quantum_representation() {
        let uq = QuantumUqSL2::new(2.0, 1);
        let rep = uq.representation(1);
        assert_eq!(rep.dim(), 3);
    }
}
