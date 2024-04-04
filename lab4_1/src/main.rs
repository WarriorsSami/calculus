use calculus_lib::Matrix;

fn main() {
    let mut A = Matrix::new(vec![
        vec![2.0, 2.0, 3.0],
        vec![4.0, 5.0, 6.0],
        vec![1.0, 2.0, 2.0],
    ]);

    // LR factorization
    for k in 1..A.n() {
        if A.at(k, k) != 0.0 {
            for i in (k + 1)..=A.n() {
                for j in (k + 1)..=A.n() {
                    A.set(i, j, A.at(i, j) - A.at(i, k) * A.at(k, j) / A.at(k, k));
                }
                A.set(i, k, A.at(i, k) / A.at(k, k));
            }
        } else {
            println!(
                "Cannot continue factorization: zero pivot found at A[{}, {}]",
                k, k
            );
        }
    }

    println!("{}", A);
}
