use crate::hopf_algebra::{Algebra, Field};
use num_complex::Complex64;
use serde::{Deserialize, Serialize};

/// A Lie algebra element.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct LieAlgebraElement {
    pub coeffs: Vec<Field>,
}

impl LieAlgebraElement {
    pub fn zero(dim: usize) -> Self {
        Self { coeffs: vec![Complex64::new(0.0, 0.0); dim] }
    }

    pub fn basis(dim: usize, idx: usize) -> Self {
        let mut v = vec![Complex64::new(0.0, 0.0); dim];
        if idx < dim {
            v[idx] = Complex64::new(1.0, 0.0);
        }
        Self { coeffs: v }
    }
}

/// A Lie algebra with structure constants [e_i, e_j] = sum_k c_{ij}^k e_k.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct LieAlgebra {
    pub dim: usize,
    /// Structure constants: bracket[i][j][k] = c_{ij}^k
    pub bracket: Vec<Vec<Vec<Field>>>,
}

impl LieAlgebra {
    pub fn new(dim: usize, bracket: Vec<Vec<Vec<Field>>>) -> Self {
        Self { dim, bracket }
    }

    /// Compute the Lie bracket [x, y].
    pub fn bracket(&self, x: &LieAlgebraElement, y: &LieAlgebraElement) -> LieAlgebraElement {
        let mut result = vec![Complex64::new(0.0, 0.0); self.dim];
        for i in 0..self.dim {
            for j in 0..self.dim {
                let c = x.coeffs[i] * y.coeffs[j];
                for k in 0..self.dim {
                    result[k] = result[k] + c * self.bracket[i][j][k];
                }
            }
        }
        LieAlgebraElement { coeffs: result }
    }

    /// Check antisymmetry: [x, y] = -[y, x].
    pub fn is_antisymmetric(&self) -> bool {
        let eps = 1e-8;
        for i in 0..self.dim {
            for j in 0..self.dim {
                for k in 0..self.dim {
                    if (self.bracket[i][j][k] + self.bracket[j][i][k]).norm() > eps {
                        return false;
                    }
                }
            }
        }
        true
    }

    /// Check Jacobi identity: [x,[y,z]] + [y,[z,x]] + [z,[x,y]] = 0.
    pub fn satisfies_jacobi(&self) -> bool {
        let eps = 1e-8;
        for i in 0..self.dim {
            for j in 0..self.dim {
                for k in 0..self.dim {
                    // [e_i, [e_j, e_k]] + cyclic
                    for m in 0..self.dim {
                        let val: Field = (0..self.dim)
                            .map(|l| {
                                self.bracket[j][k][l] * self.bracket[i][l][m]
                                    + self.bracket[k][i][l] * self.bracket[j][l][m]
                                    + self.bracket[i][j][l] * self.bracket[k][l][m]
                            })
                            .sum();
                        if val.norm() > eps {
                            return false;
                        }
                    }
                }
            }
        }
        true
    }

    /// sl(2) Lie algebra with basis {e, f, h}.
    /// [h,e] = 2e, [h,f] = -2f, [e,f] = h
    pub fn sl2() -> Self {
        let dim = 3;
        let mut bracket = vec![vec![vec![Complex64::new(0.0, 0.0); dim]; dim]; dim];
        // Basis: e_0 = e, e_1 = f, e_2 = h
        // [e, f] = h → bracket[0][1][2] = 1
        bracket[0][1][2] = Complex64::new(1.0, 0.0);
        // [f, e] = -h → bracket[1][0][2] = -1
        bracket[1][0][2] = Complex64::new(-1.0, 0.0);
        // [h, e] = 2e → bracket[2][0][0] = 2
        bracket[2][0][0] = Complex64::new(2.0, 0.0);
        // [e, h] = -2e → bracket[0][2][0] = -2
        bracket[0][2][0] = Complex64::new(-2.0, 0.0);
        // [h, f] = -2f → bracket[2][1][1] = -2
        bracket[2][1][1] = Complex64::new(-2.0, 0.0);
        // [f, h] = 2f → bracket[1][2][1] = 2
        bracket[1][2][1] = Complex64::new(2.0, 0.0);

        Self { dim, bracket }
    }

    /// Build the universal enveloping algebra U(g).
    /// Uses PBW basis with ordered monomials.
    pub fn universal_enveloping(&self, max_degree: usize) -> UniversalEnvelopingAlgebra {
        UniversalEnvelopingAlgebra::from_lie_algebra(self, max_degree)
    }
}

/// Universal enveloping algebra U(g) of a Lie algebra g.
/// Represented as polynomials in PBW basis up to max_degree.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct UniversalEnvelopingAlgebra {
    pub lie_dim: usize,
    pub max_degree: usize,
    /// Monomials: list of exponent vectors [d_0, d_1, ..., d_{n-1}]
    pub monomials: Vec<Vec<usize>>,
    /// Multiplication table (modulo the Lie bracket relations).
    pub mult_table: Vec<Vec<Vec<Field>>>,
}

impl UniversalEnvelopingAlgebra {
    /// Construct U(g) from a Lie algebra.
    pub fn from_lie_algebra(lie: &LieAlgebra, max_degree: usize) -> Self {
        // Generate all monomials e_0^{d_0} e_1^{d_1} ... with sum(d_i) <= max_degree
        let monomials = generate_monomials(lie.dim, max_degree);
        let n = monomials.len();

        // For simplicity, build a free algebra multiplication table
        // (in a full implementation, we'd reduce using PBW relations)
        let mut mult_table = vec![vec![vec![Complex64::new(0.0, 0.0); n]; n]; n];

        for i in 0..n {
            for j in 0..n {
                // Multiply monomials: add exponent vectors
                let mut product = vec![0usize; lie.dim];
                for k in 0..lie.dim {
                    product[k] = monomials[i][k] + monomials[j][k];
                }
                let total: usize = product.iter().sum();
                if total <= max_degree {
                    if let Some(idx) = monomials.iter().position(|m| m == &product) {
                        mult_table[i][j][idx] = Complex64::new(1.0, 0.0);
                    }
                }
            }
        }

        Self { lie_dim: lie.dim, max_degree, monomials, mult_table }
    }

    /// Number of basis elements.
    pub fn dim(&self) -> usize {
        self.monomials.len()
    }

    /// Get the PBW algebra.
    pub fn to_algebra(&self) -> Algebra {
        let n = self.dim();
        let dim = n.max(1);
        let mut mult = vec![vec![vec![Complex64::new(0.0, 0.0); dim]; dim]; dim];
        for i in 0..n {
            for j in 0..n {
                for k in 0..n {
                    mult[i][j][k] = self.mult_table[i][j][k];
                }
            }
        }
        Algebra { dim, mult }
    }
}

fn generate_monomials(dim: usize, max_degree: usize) -> Vec<Vec<usize>> {
    let mut result = Vec::new();
    let mut current = vec![0usize; dim];
    loop {
        result.push(current.clone());
        if !increment_monomial(&mut current, dim, max_degree) {
            break;
        }
    }
    result
}

fn increment_monomial(m: &mut Vec<usize>, dim: usize, max_degree: usize) -> bool {
    let _total: usize = m.iter().sum();
    // Try to increment in lex order
    for i in (0..dim).rev() {
        let partial: usize = m[0..=i].iter().sum();
        if partial < max_degree {
            m[i] += 1;
            for j in i + 1..dim {
                m[j] = 0;
            }
            return true;
        }
    }
    false
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hopf_algebra::real_val;

    #[test]
    fn test_sl2_antisymmetry() {
        let sl2 = LieAlgebra::sl2();
        assert!(sl2.is_antisymmetric());
    }

    #[test]
    fn test_sl2_jacobi() {
        let sl2 = LieAlgebra::sl2();
        assert!(sl2.satisfies_jacobi());
    }

    #[test]
    fn test_sl2_bracket_ef() {
        let sl2 = LieAlgebra::sl2();
        let e = LieAlgebraElement::basis(3, 0);
        let f = LieAlgebraElement::basis(3, 1);
        let comm = sl2.bracket(&e, &f);
        // [e, f] = h
        assert!((comm.coeffs[2] - real_val(1.0)).norm() < 1e-10);
    }

    #[test]
    fn test_sl2_bracket_he() {
        let sl2 = LieAlgebra::sl2();
        let h = LieAlgebraElement::basis(3, 2);
        let e = LieAlgebraElement::basis(3, 0);
        let comm = sl2.bracket(&h, &e);
        // [h, e] = 2e
        assert!((comm.coeffs[0] - real_val(2.0)).norm() < 1e-10);
    }

    #[test]
    fn test_sl2_bracket_hf() {
        let sl2 = LieAlgebra::sl2();
        let h = LieAlgebraElement::basis(3, 2);
        let f = LieAlgebraElement::basis(3, 1);
        let comm = sl2.bracket(&h, &f);
        // [h, f] = -2f
        assert!((comm.coeffs[1] - real_val(-2.0)).norm() < 1e-10);
    }

    #[test]
    fn test_uea_sl2_construction() {
        let sl2 = LieAlgebra::sl2();
        let uea = sl2.universal_enveloping(2);
        assert!(uea.dim() > 0);
    }

    #[test]
    fn test_uea_monomials() {
        let sl2 = LieAlgebra::sl2();
        let uea = sl2.universal_enveloping(2);
        // With 3 generators and max_degree 2:
        // (0,0,0), (1,0,0), (2,0,0), (0,1,0), (0,2,0), (0,0,1), (0,0,2), (1,1,0), (1,0,1), (0,1,1)
        assert!(uea.dim() >= 10);
    }

    #[test]
    fn test_uea_to_algebra() {
        let sl2 = LieAlgebra::sl2();
        let uea = sl2.universal_enveloping(2);
        let alg = uea.to_algebra();
        assert_eq!(alg.dim, uea.dim());
    }

    #[test]
    fn test_lie_element_zero() {
        let z = LieAlgebraElement::zero(3);
        assert!(z.coeffs.iter().all(|c| c.norm() < 1e-10));
    }

    #[test]
    fn test_lie_element_basis() {
        let e = LieAlgebraElement::basis(3, 1);
        assert!((e.coeffs[1] - real_val(1.0)).norm() < 1e-10);
        assert!(e.coeffs[0].norm() < 1e-10);
    }
}
