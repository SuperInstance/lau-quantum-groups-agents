use num_complex::Complex64;
use serde::{Deserialize, Serialize};

/// Field element type used throughout the library.
pub type Field = Complex64;

/// Create a real field element.
pub fn real_val(x: f64) -> Field {
    Complex64::new(x, 0.0)
}

/// Create a complex field element.
pub fn complex_val(re: f64, im: f64) -> Field {
    Complex64::new(re, im)
}

/// Basis element index.
pub type BasisIndex = usize;

/// A finite-dimensional algebra element represented as a vector over a basis.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct AlgebraElement {
    /// Coefficients in the algebra basis.
    pub coeffs: Vec<Field>,
}

impl AlgebraElement {
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

    pub fn dim(&self) -> usize {
        self.coeffs.len()
    }

    pub fn add(&self, other: &Self) -> Self {
        let dim = self.dim().max(other.dim());
        let mut coeffs = Vec::with_capacity(dim);
        for i in 0..dim {
            let a = self.coeffs.get(i).copied().unwrap_or(Complex64::new(0.0, 0.0));
            let b = other.coeffs.get(i).copied().unwrap_or(Complex64::new(0.0, 0.0));
            coeffs.push(a + b);
        }
        Self { coeffs }
    }

    pub fn scale(&self, s: Field) -> Self {
        Self { coeffs: self.coeffs.iter().map(|c| s * c).collect() }
    }

    pub fn is_zero(&self) -> bool {
        self.coeffs.iter().all(|c| c.norm() < 1e-10)
    }
}

/// A finite-dimensional associative algebra with structure constants.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Algebra {
    /// Dimension of the algebra.
    pub dim: usize,
    /// Structure constants: mult[i][j][k] = coefficient of e_k in e_i * e_j.
    pub mult: Vec<Vec<Vec<Field>>>,
}

impl Algebra {
    /// Create a new algebra with given dimension and structure constants.
    pub fn new(dim: usize, mult: Vec<Vec<Vec<Field>>>) -> Self {
        assert_eq!(mult.len(), dim);
        for i in &mult {
            assert_eq!(i.len(), dim);
            for j in i {
                assert_eq!(j.len(), dim);
            }
        }
        Self { dim, mult }
    }

    /// Create the trivial 1-dimensional algebra.
    pub fn trivial() -> Self {
        let mut mult = vec![vec![vec![Complex64::new(0.0, 0.0); 1]; 1]; 1];
        mult[0][0][0] = Complex64::new(1.0, 0.0);
        Self { dim: 1, mult }
    }

    /// Multiply two algebra elements.
    pub fn multiply(&self, a: &AlgebraElement, b: &AlgebraElement) -> AlgebraElement {
        let mut result = vec![Complex64::new(0.0, 0.0); self.dim];
        for i in 0..self.dim {
            for j in 0..self.dim {
                let c = a.coeffs[i] * b.coeffs[j];
                for k in 0..self.dim {
                    result[k] += c * self.mult[i][j][k];
                }
            }
        }
        AlgebraElement { coeffs: result }
    }

    /// Get the unit element.
    pub fn unit(&self) -> AlgebraElement {
        // Find e_0 as unit (convention: first basis element is unit)
        let mut coeffs = vec![Complex64::new(0.0, 0.0); self.dim];
        coeffs[0] = Complex64::new(1.0, 0.0);
        AlgebraElement { coeffs }
    }

    /// Check associativity.
    pub fn is_associative(&self) -> bool {
        let eps = 1e-8;
        for i in 0..self.dim {
            for j in 0..self.dim {
                for k in 0..self.dim {
                    // (e_i * e_j) * e_k vs e_i * (e_j * e_k)
                    for m in 0..self.dim {
                        let lhs: Field = (0..self.dim)
                            .map(|l| self.mult[i][j][l] * self.mult[l][k][m])
                            .sum();
                        let rhs: Field = (0..self.dim)
                            .map(|l| self.mult[i][l][m] * self.mult[j][k][l])
                            .sum();
                        if (lhs - rhs).norm() > eps {
                            return false;
                        }
                    }
                }
            }
        }
        true
    }

    /// Matrix algebra M_n(k).
    pub fn matrix_algebra(n: usize) -> Self {
        let dim = n * n;
        let mut mult = vec![vec![vec![Complex64::new(0.0, 0.0); dim]; dim]; dim];

        for i1 in 0..n {
            for j1 in 0..n {
                for i2 in 0..n {
                    for j2 in 0..n {
                        let a = i1 * n + j1;
                        let b = i2 * n + j2;
                        if j1 == i2 {
                            let c = i1 * n + j2;
                            mult[a][b][c] = Complex64::new(1.0, 0.0);
                        }
                    }
                }
            }
        }

        Self { dim, mult }
    }
}

/// Multiplication in the tensor product algebra.
pub fn tensor_multiply(
    a: &AlgebraElement,
    b: &AlgebraElement,
    dim1: usize,
    dim2: usize,
) -> AlgebraElement {
    let dim = dim1 * dim2;
    let mut result = vec![Complex64::new(0.0, 0.0); dim];
    for i in 0..dim1 {
        for j in 0..dim2 {
            for k in 0..dim1 {
                for l in 0..dim2 {
                    let idx = i * dim2 + j;
                    let idx2 = k * dim2 + l;
                    let _out = i * dim2 + j + k * dim2 + l;
                    // Simple tensor product multiplication
                    let r = a.coeffs.get(idx).copied().unwrap_or(Complex64::new(0.0, 0.0))
                        * b.coeffs.get(idx2).copied().unwrap_or(Complex64::new(0.0, 0.0));
                    if r.norm() > 1e-15 {
                        result[idx] += r;
                    }
                }
            }
        }
    }
    AlgebraElement { coeffs: result }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_algebra_element_zero() {
        let z = AlgebraElement::zero(3);
        assert!(z.is_zero());
        assert_eq!(z.dim(), 3);
    }

    #[test]
    fn test_algebra_element_basis() {
        let e = AlgebraElement::basis(3, 1);
        assert!(!e.is_zero());
        assert_eq!(e.coeffs[1], Complex64::new(1.0, 0.0));
        assert_eq!(e.coeffs[0], Complex64::new(0.0, 0.0));
    }

    #[test]
    fn test_trivial_algebra() {
        let alg = Algebra::trivial();
        assert_eq!(alg.dim, 1);
        assert!(alg.is_associative());
    }

    #[test]
    fn test_trivial_algebra_multiply() {
        let alg = Algebra::trivial();
        let u = alg.unit();
        let prod = alg.multiply(&u, &u);
        assert!((prod.coeffs[0] - Complex64::new(1.0, 0.0)).norm() < 1e-10);
    }

    #[test]
    fn test_matrix_algebra_m2_associative() {
        let alg = Algebra::matrix_algebra(2);
        assert_eq!(alg.dim, 4);
        assert!(alg.is_associative());
    }

    #[test]
    fn test_matrix_algebra_m3() {
        let alg = Algebra::matrix_algebra(3);
        assert_eq!(alg.dim, 9);
        assert!(alg.is_associative());
    }

    #[test]
    fn test_algebra_add() {
        let a = AlgebraElement { coeffs: vec![real_val(1.0), real_val(2.0)] };
        let b = AlgebraElement { coeffs: vec![real_val(3.0), real_val(4.0)] };
        let c = a.add(&b);
        assert!((c.coeffs[0] - real_val(4.0)).norm() < 1e-10);
        assert!((c.coeffs[1] - real_val(6.0)).norm() < 1e-10);
    }

    #[test]
    fn test_algebra_scale() {
        let a = AlgebraElement { coeffs: vec![real_val(2.0), real_val(3.0)] };
        let b = a.scale(real_val(0.5));
        assert!((b.coeffs[0] - real_val(1.0)).norm() < 1e-10);
    }

    #[test]
    fn test_algebra_unit() {
        let alg = Algebra::trivial();
        let u = alg.unit();
        let e = AlgebraElement::basis(1, 0);
        let prod = alg.multiply(&u, &e);
        assert!((prod.coeffs[0] - real_val(1.0)).norm() < 1e-10);
    }
}
