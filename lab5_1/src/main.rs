use calculus_lib::Matrix;

fn qr_factorization(matrix: &mut Vec<Vec<f64>>) -> (Vec<Vec<f64>>, Vec<Vec<f64>>) {
    let n = matrix.len();
    let mut q = vec![vec![0.0; n]; n];
    let mut r = matrix.clone();

    // Initialize Q as the identity matrix
    for i in 0..n {
        q[i][i] = 1.0;
    }

    for k in 0..n - 1 {
        let mut sigma = 0.0;
        for i in k..n {
            sigma += r[i][k] * r[i][k];
        }
        sigma = sigma.sqrt();

        if sigma != 0.0 {
            if r[k][k] < 0.0 {
                sigma = -sigma;
            }

            let mut vk = vec![0.0; n];
            vk[k] = r[k][k] + sigma;

            for i in k + 1..n {
                vk[i] = r[i][k];
            }

            let beta = 0.5 * vk.iter().map(|&x| x * x).sum::<f64>();
            let mut h = vec![vec![0.0; n]; n];

            for i in 0..n {
                for j in 0..n {
                    h[i][j] = if i != j {
                        -vk[i] * vk[j] / beta
                    } else {
                        1.0 - vk[i] * vk[i] / beta
                    };
                }
            }

            r = _multiply_matrices(h.clone(), r);
            q = _multiply_matrices(q, h);
        }
    }

    (q, r)
}

fn _multiply_matrices(a: Vec<Vec<f64>>, b: Vec<Vec<f64>>) -> Vec<Vec<f64>> {
    let n = a.len();
    let m = b[0].len();
    let p = b.len();
    let mut c = vec![vec![0.0; m]; n];

    for i in 0..n {
        for j in 0..m {
            for k in 0..p {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }

    c
}

fn main() {
    let mut matrix: Vec<Vec<f64>> = vec![
        vec![1.0, 0.0, -1.0],
        vec![1.0, 1.0, 0.0],
        vec![0.0, -1.0, 1.0],
    ];

    let (q, r) = qr_factorization(&mut matrix);

    let Q = Matrix::new(q);
    let R = Matrix::new(r);

    println!("Q = {}", Q);
    println!("R = {}", R);
}
