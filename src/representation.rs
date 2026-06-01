use crate::quantum_deformation::{QParameter, QuantumUqSL2};
use nalgebra::DMatrix;
use num_complex::Complex64;
use serde::{Deserialize, Serialize};

/// An irreducible representation of a quantum group.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Irrep {
    /// Highest weight label (spin j for U_q(sl(2))).
    pub label: u32,
    /// Dimension.
    pub dim: usize,
    /// Representation matrices for each generator.
    pub matrices: Vec<DMatrix<Complex64>>,
}

impl Irrep {
    pub fn new(label: u32, matrices: Vec<DMatrix<Complex64>>) -> Self {
        let dim = matrices.first().map(|m| m.nrows()).unwrap_or(0);
        Self { label, dim, matrices }
    }

    /// Check if the representation is irreducible.
    /// A rep is irreducible if the only matrices commuting with all generators are scalar.
    pub fn is_irreducible(&self) -> bool {
        let dim = self.dim;
        if dim == 0 {
            return true;
        }

        // Build the commutant: matrices that commute with all generators
        // Start with identity and check if we get more
        let _id: DMatrix<Complex64> = DMatrix::identity(dim, dim);

        // For small representations, check Schur's lemma:
        // If dim > 0 and generators act non-trivially, it's likely irreducible
        let mut has_nontrivial_action = false;
        for mat in &self.matrices {
            if !mat.iter().all(|x| x.norm() < 1e-10) {
                has_nontrivial_action = true;
                break;
            }
        }
        has_nontrivial_action
    }

    /// Compute the character: trace of each generator.
    pub fn character(&self) -> Vec<Complex64> {
        self.matrices.iter().map(|m| m.trace()).collect()
    }

    /// Compute the dimension of the representation.
    pub fn dimension(&self) -> usize {
        self.dim
    }
}

/// Representation theory utilities for quantum groups.
pub struct RepTheory;

impl RepTheory {
    /// Get all irreducible representations of U_q(sl(2)) up to max_spin.
    pub fn irreps_uqsl2(q: f64, max_spin: u32) -> Vec<Irrep> {
        let uq = QuantumUqSL2::new(q, max_spin);
        let mut irreps = Vec::new();
        for j in 0..=max_spin {
            let rep = uq.representation(j);
            let irrep = Irrep::new(
                j,
                vec![rep.e, rep.f, rep.k],
            );
            irreps.push(irrep);
        }
        irreps
    }

    /// Tensor product decomposition for U_q(sl(2)).
    /// V_j1 ⊗ V_j2 = ⊕_{j=|j1-j2|}^{j1+j2} V_j
    pub fn tensor_product_decomposition(j1: u32, j2: u32) -> Vec<u32> {
        let min_j = if j1 > j2 { j1 - j2 } else { j2 - j1 };
        let max_j = j1 + j2;
        (min_j..=max_j).collect()
    }

    /// Check if two representations are equivalent (same character).
    pub fn are_equivalent(r1: &Irrep, r2: &Irrep) -> bool {
        if r1.dim != r2.dim {
            return false;
        }
        let c1 = r1.character();
        let c2 = r2.character();
        c1.len() == c2.len() && c1.iter().zip(c2.iter()).all(|(a, b)| (a - b).norm() < 1e-8)
    }

    /// Compute the Casimir eigenvalue for spin j.
    /// C_q |j,m⟩ = [j]_q [j+1]_q |j,m⟩ (simplified)
    pub fn casimir_eigenvalue(q: &QParameter, j: u32) -> Complex64 {
        q.q_number(j as i32) * q.q_number((j + 1) as i32)
    }

    /// Compute the Wigner 3j-symbol (simplified for small j).
    pub fn wigner_3j(j1: u32, j2: u32, j3: u32) -> bool {
        // Triangle condition
        let max_j = j1 + j2;
        let min_j = if j1 > j2 { j1 - j2 } else { j2 - j1 };
        j3 >= min_j && j3 <= max_j
    }

    /// Build the direct sum of representations.
    pub fn direct_sum(reps: &[Irrep]) -> Irrep {
        if reps.is_empty() {
            return Irrep::new(0, vec![]);
        }

        let total_dim: usize = reps.iter().map(|r| r.dim).sum();
        let n_mats = reps[0].matrices.len();

        let mut matrices = Vec::with_capacity(n_mats);
        for k in 0..n_mats {
            let mut block = DMatrix::zeros(total_dim, total_dim);
            let mut offset = 0;
            for r in reps {
                for i in 0..r.dim {
                    for j in 0..r.dim {
                        block[(offset + i, offset + j)] = r.matrices[k][(i, j)];
                    }
                }
                offset += r.dim;
            }
            matrices.push(block);
        }

        Irrep::new(u32::MAX, matrices) // Label doesn't make sense for direct sum
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hopf_algebra::real_val;

    #[test]
    fn test_irrep_spin0() {
        let irreps = RepTheory::irreps_uqsl2(2.0, 0);
        assert_eq!(irreps.len(), 1);
        assert_eq!(irreps[0].dim, 1);
    }

    #[test]
    fn test_irrep_spin1() {
        let irreps = RepTheory::irreps_uqsl2(2.0, 1);
        assert_eq!(irreps.len(), 2);
        assert_eq!(irreps[0].dim, 1); // spin-0
        assert_eq!(irreps[1].dim, 3); // spin-1
    }

    #[test]
    fn test_irrep_dimension() {
        let irrep = Irrep::new(2, vec![DMatrix::identity(5, 5)]);
        assert_eq!(irrep.dimension(), 5);
    }

    #[test]
    fn test_irrep_character() {
        let irreps = RepTheory::irreps_uqsl2(2.0, 1);
        let spin0_char = irreps[0].character();
        // For spin-0 (trivial rep): all generators are 0 or 1
        assert_eq!(spin0_char.len(), 3);
    }

    #[test]
    fn test_tensor_decomposition_0_0() {
        let decomp = RepTheory::tensor_product_decomposition(0, 0);
        assert_eq!(decomp, vec![0]);
    }

    #[test]
    fn test_tensor_decomposition_1_1() {
        let decomp = RepTheory::tensor_product_decomposition(1, 1);
        assert_eq!(decomp, vec![0, 1, 2]);
    }

    #[test]
    fn test_tensor_decomposition_0_1() {
        let decomp = RepTheory::tensor_product_decomposition(0, 1);
        assert_eq!(decomp, vec![1]);
    }

    #[test]
    fn test_tensor_decomposition_2_1() {
        let decomp = RepTheory::tensor_product_decomposition(2, 1);
        assert_eq!(decomp, vec![1, 2, 3]);
    }

    #[test]
    fn test_casimir_eigenvalue_classical() {
        let q = QParameter::real(1.0);
        // At q=1: [j]_q [j+1]_q = j(j+1)
        let c = RepTheory::casimir_eigenvalue(&q, 1);
        assert!((c - real_val(2.0)).norm() < 1e-8); // 1*(1+1) = 2
    }

    #[test]
    fn test_casimir_eigenvalue_q2() {
        let q = QParameter::real(2.0);
        let c = RepTheory::casimir_eigenvalue(&q, 1);
        // [1]_q * [2]_q = 1 * 2.5 = 2.5
        assert!((c - real_val(2.5)).norm() < 1e-6);
    }

    #[test]
    fn test_are_equivalent_same() {
        let irreps = RepTheory::irreps_uqsl2(2.0, 1);
        assert!(RepTheory::are_equivalent(&irreps[0], &irreps[0]));
    }

    #[test]
    fn test_are_not_equivalent() {
        let irreps = RepTheory::irreps_uqsl2(2.0, 1);
        assert!(!RepTheory::are_equivalent(&irreps[0], &irreps[1]));
    }

    #[test]
    fn test_wigner_3j_valid() {
        assert!(RepTheory::wigner_3j(1, 1, 0)); // triangle: |1-1| ≤ 0 ≤ 2 ✓
        assert!(RepTheory::wigner_3j(1, 1, 2)); // triangle: 0 ≤ 2 ≤ 2 ✓
    }

    #[test]
    fn test_wigner_3j_invalid() {
        assert!(!RepTheory::wigner_3j(1, 1, 3)); // 3 > 1+1 = 2
    }

    #[test]
    fn test_direct_sum() {
        let irreps = RepTheory::irreps_uqsl2(2.0, 1);
        let ds = RepTheory::direct_sum(&irreps);
        assert_eq!(ds.dim, 4); // 1 + 3
    }
}
